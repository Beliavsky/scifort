! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Copyright (c) 2021-2022 Jacob Williams
!
! Inexact-Newton nonlinear solvers for the scipy.optimize translation.
!
! The iteration controller, Armijo line search, autoscaling, Broyden inverse
! updates/rank reduction, Anderson/diagonal/scalar mixing, and matrix-free
! Krylov behavior are adapted from SciPy 1.17.0 scipy/optimize/_nonlin.py under
! SciPy's BSD-3-Clause license. The dense Broyden first-Jacobian rank-one formula
! was also checked against Jacob Williams' BSD-3-Clause nlesolver-fortran 1.1.0.

module scifort_optimize_nonlin
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_iterative, only : krylov_iterative_solve
    use scifort_optimize_krylov_rank_ls, only : krylov_rank_least_squares
    use scifort_optimize_krylov_augmentation, only : retain_lgmres_direction, &
        retain_lgmres_image_direction, lgmres_pruning_workspace, &
        reset_lgmres_image_pruning_stream, retain_lgmres_image_direction_stream, &
        retain_lgmres_image_candidate_cached, &
        reset_lgmres_direction_pruning_stream, retain_lgmres_direction_stream, &
        retain_lgmres_direction_fifo
    use scifort_linalg, only : linalg_status_success, symmetric_eigen_jacobi
    use scifort_optimize_wolfe, only : scalar_search_wolfe1
    implicit none
    private

    abstract interface
        subroutine nonlin_function(x, fvec)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent variable vector
            real(dp), intent(out) :: fvec(:) !! square residual vector at x
        end subroutine nonlin_function

        function nonlin_norm_function(values) result(value)
            import dp
            implicit none
            real(dp), intent(in) :: values(:) !! vector to measure
            real(dp) :: value !! nonnegative vector norm
        end function nonlin_norm_function

        subroutine nonlin_linear_operator(vec, product, ok)
            import dp
            implicit none
            real(dp), intent(in) :: vec(:) !! vector to transform
            real(dp), intent(out) :: product(:) !! transformed vector
            logical, intent(out) :: ok !! false if operator application fails
        end subroutine nonlin_linear_operator

        subroutine nonlin_jacobian_product(x, vec, product, ok)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! current nonlinear point
            real(dp), intent(in) :: vec(:) !! Jacobian multiplication direction
            real(dp), intent(out) :: product(:) !! exact or user-computed J(x)*vec
            logical, intent(out) :: ok !! false when the action cannot be evaluated
        end subroutine nonlin_jacobian_product

        subroutine nonlin_base_diagonal(x, diagonal, ok)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! current nonlinear iterate
            real(dp), intent(out) :: diagonal(:) !! uncoupled diagonal D(x), not diag(J(x))
            logical, intent(out) :: ok !! false when the diagonal cannot be evaluated
        end subroutine nonlin_base_diagonal

        subroutine nonlin_preconditioner_update(x, fvec)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! accepted nonlinear iterate
            real(dp), intent(in) :: fvec(:) !! corresponding residual
        end subroutine nonlin_preconditioner_update

        subroutine nonlin_inner_solver(matvec, preconditioner, rhs, rtol, &
                maxiter, solution, ok)
            import dp, nonlin_linear_operator
            implicit none
            procedure(nonlin_linear_operator) :: matvec !! matrix-free Jacobian product
            procedure(nonlin_linear_operator) :: preconditioner !! inverse-Jacobian approximation or identity
            real(dp), intent(in) :: rhs(:) !! unnormalized right-hand side
            real(dp), intent(in) :: rtol !! requested relative linear tolerance
            integer, intent(in) :: maxiter !! inner iteration budget
            real(dp), intent(out) :: solution(:) !! approximate J**(-1) * rhs
            logical, intent(out) :: ok !! false if inner linear solve fails
        end subroutine nonlin_inner_solver
    end interface

    public :: nonlin_function
    public :: nonlin_norm_function
    public :: nonlin_linear_operator
    public :: nonlin_jacobian_product
    public :: nonlin_base_diagonal
    public :: nonlin_preconditioner_update
    public :: nonlin_inner_solver
    public :: solve_krylov_root
    public :: solve_nonlin_root

    integer, parameter :: jac_step_global = 0
    integer, parameter :: jac_step_componentwise = 1
    integer, parameter :: jac_step_adaptive = 2

    ! Keep the Arnoldi workspace for the lifetime of one nonlinear solve.
    ! It belongs to a local solver instance, not to module SAVE state, so
    ! independent or concurrent root solves cannot share mutable storage.
    type :: fgmres_workspace
        real(dp), allocatable :: v(:, :)
        real(dp), allocatable :: z(:, :)
        real(dp), allocatable :: h(:, :)
        real(dp), allocatable :: cs(:)
        real(dp), allocatable :: sn(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: raw_products(:, :)
        type(lgmres_pruning_workspace) :: pruning
    end type fgmres_workspace

    ! Current-Jacobian image filtering uses a separate workspace from
    ! Arnoldi.  This permits images to be passed to Arnoldi as read-only
    ! inputs without aliasing the Arnoldi scratch arrays.  These arrays live
    ! only for the duration of one nonlinear root solve.
    type :: lgmres_image_arrays
        real(dp), allocatable :: images(:, :)
        real(dp), allocatable :: image(:)
    end type lgmres_image_arrays

contains

    subroutine apply_nonlin_linear_operator(operator, vec, product, ok)
        ! Keep the invocation of an optional, host-associated procedure dummy
        ! outside the internal Krylov callback.  Some compilers diagnose that
        ! nested call as having an implicit interface even though the dummy
        ! has the explicit nonlin_linear_operator interface.
        procedure(nonlin_linear_operator) :: operator !! supplied linear operator callback
        real(dp), intent(in) :: vec(:) !! vector to transform
        real(dp), intent(out) :: product(:) !! operator action on the vector
        logical, intent(out) :: ok !! true if the operator application succeeded

        call operator(vec, product, ok)
    end subroutine apply_nonlin_linear_operator

    subroutine ensure_lgmres_image_arrays(arrays, n, ncols)
        type(lgmres_image_arrays), intent(inout) :: arrays !! reusable image-filter storage
        integer, intent(in) :: n !! nonlinear system dimension
        integer, intent(in) :: ncols !! required active image capacity

        if (allocated(arrays%images)) then
            if (size(arrays%images, 1) /= n .or. &
                    size(arrays%images, 2) < ncols) deallocate(arrays%images)
        end if
        if (.not. allocated(arrays%images)) allocate(arrays%images(n, ncols))
        if (allocated(arrays%image)) then
            if (size(arrays%image) /= n) deallocate(arrays%image)
        end if
        if (.not. allocated(arrays%image)) allocate(arrays%image(n))
    end subroutine ensure_lgmres_image_arrays

    subroutine solve_nonlin_root(fun, x0, method, x, fvec, nfev, nit, status, &
            success, message, tol, ftol, fatol, xtol, xatol, maxiter, iterations, &
            line_search, alpha, anderson_m, w0, alphamax, reduction_method, max_rank, &
            tol_norm)
        procedure(nonlin_function) :: fun !! square residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        character(len=*), intent(in) :: method !! inexact-Newton method name
        real(dp), intent(out) :: x(:) !! final solution estimate
        real(dp), intent(out) :: fvec(:) !! residual vector at x
        integer, intent(out) :: nfev !! residual evaluations
        integer, intent(out) :: nit !! termination-check iteration count
        integer, intent(out) :: status !! 1 success, 2 iteration limit, 0 invalid/failure
        logical, intent(out) :: success !! true when status is 1
        character(len=*), intent(out) :: message !! termination description
        real(dp), intent(in), optional :: tol !! generic root tolerance
        real(dp), intent(in), optional :: ftol !! relative residual tolerance
        real(dp), intent(in), optional :: fatol !! absolute max-norm residual tolerance
        real(dp), intent(in), optional :: xtol !! relative step tolerance
        real(dp), intent(in), optional :: xatol !! absolute max-norm step tolerance
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        integer, intent(in), optional :: iterations !! exact number of iterations to make
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial Jacobian scaling parameter
        integer, intent(in), optional :: anderson_m !! retained Anderson history, default 5
        real(dp), intent(in), optional :: w0 !! Anderson regularization, default 0.01
        real(dp), intent(in), optional :: alphamax !! exciting-mixing upper scale, default 1
        character(len=*), intent(in), optional :: reduction_method !! Broyden rank reduction
        integer, intent(in), optional :: max_rank !! maximum Broyden inverse rank
        procedure(nonlin_norm_function), optional :: tol_norm !! stopping-test norm, default max norm

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: broyden_c(:, :)
        real(dp), allocatable :: broyden_d(:, :)
        real(dp), allocatable :: c(:)
        real(dp), allocatable :: rank_d(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: df(:)
        real(dp), allocatable :: dx(:)
        real(dp), allocatable :: dx_actual(:)
        real(dp), allocatable :: fnew(:)
        real(dp), allocatable :: gamma(:)
        real(dp), allocatable :: h(:, :)
        real(dp), allocatable :: hist_df(:, :)
        real(dp), allocatable :: hist_dx(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: solve_vec(:)
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: xnew(:)
        real(dp) :: alpha_use
        real(dp) :: alphamax_use
        real(dp) :: denom
        real(dp) :: df_norm
        real(dp) :: dx_norm
        real(dp) :: f0_norm
        real(dp) :: f_abs_tol
        real(dp) :: f_rel_tol
        real(dp) :: fx_norm2
        real(dp) :: w0_use
        real(dp) :: x_abs_tol
        real(dp) :: x_rel_tol
        integer :: anderson_m_use
        integer :: broyden_count
        integer :: hist_count
        integer :: max_rank_use
        integer :: i
        integer :: iter_count
        integer :: j
        integer :: maxiter_use
        integer :: n
        integer :: step_status
        logical :: broyden_collapsed
        logical :: f_abs_active
        logical :: f_rel_active
        logical :: linear_ok
        logical :: rank_control
        logical :: use_armijo
        logical :: use_wolfe
        logical :: x_abs_active
        logical :: x_rel_active
        character(len=16) :: meth
        character(len=16) :: reduction
        character(len=16) :: search

        n = size(x0)
        x = 0.0_dp
        fvec = 0.0_dp
        nfev = 0
        nit = 0
        status = 0
        success = .false.
        message = 'improper input parameters'
        if (n <= 0 .or. size(x) /= n .or. size(fvec) /= n) return
        if (.not. all(ieee_is_finite(x0))) return

        meth = lowercase(trim(method))
        select case (trim(meth))
        case ('broyden1', 'broyden2', 'anderson', 'linearmixing', 'diagbroyden', &
                'excitingmixing')
        case default
            return
        end select

        reduction = 'restart'
        rank_control = .false.
        max_rank_use = 0
        broyden_count = 0
        broyden_collapsed = .false.
        if (trim(meth) == 'broyden1' .or. trim(meth) == 'broyden2') then
            if (present(reduction_method)) reduction = lowercase(trim(reduction_method))
            select case (trim(reduction))
            case ('restart', 'simple', 'svd')
            case default
                message = 'unsupported Broyden reduction; use restart, simple, or svd'
                return
            end select
            if (present(max_rank)) then
                if (max_rank < 2) then
                    message = 'Broyden max_rank must be at least 2'
                    return
                end if
                max_rank_use = max_rank
                rank_control = .true.
            end if
        end if

        search = 'armijo'
        if (present(line_search)) search = lowercase(trim(line_search))
        select case (trim(search))
        case ('armijo')
            use_armijo = .true.
            use_wolfe = .false.
        case ('wolfe')
            use_armijo = .false.
            use_wolfe = .true.
        case ('none', 'off', 'false')
            use_armijo = .false.
            use_wolfe = .false.
        case default
            message = 'unsupported line search; use armijo, wolfe, or none'
            return
        end select

        maxiter_use = 100 * (n + 1)
        if (present(iterations)) then
            if (iterations < 0) return
            maxiter_use = iterations + 1
        end if
        if (present(maxiter)) then
            if (maxiter < 0) return
            maxiter_use = maxiter
        end if
        if (maxiter_use < 1) then
            status = 2
            message = 'The maximum number of iterations allowed has been reached.'
            return
        end if

        call setup_tolerances(tol, ftol, fatol, xtol, xatol, f_abs_tol, &
            f_rel_tol, x_abs_tol, x_rel_tol, f_abs_active, f_rel_active, &
            x_abs_active, x_rel_active, linear_ok)
        if (.not. linear_ok) return

        x = x0
        call fun(x, fvec)
        nfev = 1
        if (.not. all(ieee_is_finite(fvec))) then
            message = 'nonfinite residual encountered'
            return
        end if

        fx_norm2 = norm2(fvec)
        alpha_use = 1.0_dp
        if (fx_norm2 > 0.0_dp) then
            alpha_use = 0.5_dp * max(norm2(x), 1.0_dp) / fx_norm2
        end if
        if (present(alpha)) alpha_use = alpha
        if (.not. ieee_is_finite(alpha_use) .or. same_real(alpha_use, 0.0_dp)) return

        anderson_m_use = 5
        if (present(anderson_m)) anderson_m_use = anderson_m
        if (anderson_m_use < 0) return
        w0_use = 0.01_dp
        if (present(w0)) w0_use = w0
        if (.not. ieee_is_finite(w0_use) .or. w0_use < 0.0_dp) return
        alphamax_use = 1.0_dp
        if (present(alphamax)) alphamax_use = alphamax
        if (.not. ieee_is_finite(alphamax_use) .or. alphamax_use < 0.0_dp) return

        allocate(dx(n), dx_actual(n), df(n), fnew(n), solve_vec(n), xnew(n))
        allocate(c(n), rank_d(n), v(n))
        dx = huge(1.0_dp)
        f0_norm = maxval(abs(fvec))
        if (present(tol_norm)) f0_norm = tol_norm(fvec)
        if (.not. ieee_is_finite(f0_norm)) then
            message = 'invalid tolerance norm value'
            return
        end if
        if (f0_norm < 0.0_dp) then
            message = 'invalid tolerance norm value'
            return
        end if
        iter_count = 0
        hist_count = 0

        select case (trim(meth))
        case ('broyden1', 'broyden2')
            allocate(h(n, n))
            call set_scaled_identity(h, -alpha_use)
            if (rank_control) then
                allocate(broyden_c(n, n + 1), broyden_d(n, n + 1))
                broyden_c = 0.0_dp
                broyden_d = 0.0_dp
            end if
        case ('diagbroyden')
            allocate(d(n))
            d = 1.0_dp / alpha_use
        case ('excitingmixing')
            allocate(beta(n))
            beta = alpha_use
        case ('anderson')
            if (anderson_m_use > 0) then
                allocate(hist_dx(n, anderson_m_use), hist_df(n, anderson_m_use))
                allocate(a(anderson_m_use, anderson_m_use))
                allocate(rhs(anderson_m_use), gamma(anderson_m_use))
                hist_dx = 0.0_dp
                hist_df = 0.0_dp
                a = 0.0_dp
            end if
        end select

        do i = 1, maxiter_use
            iter_count = iter_count + 1
            if (termination_met(fvec, x, dx, f0_norm, f_abs_tol, f_rel_tol, &
                    x_abs_tol, x_rel_tol, f_abs_active, f_rel_active, &
                    x_abs_active, x_rel_active, iterations, iter_count, status, &
                    tol_norm)) then
                exit
            end if

            select case (trim(meth))
            case ('broyden1', 'broyden2')
                solve_vec = matmul(h, fvec)
            case ('anderson')
                solve_vec = -alpha_use * fvec
                if (hist_count > 0) then
                    do j = 1, hist_count
                        rhs(j) = dot_product(hist_df(:, j), fvec)
                    end do
                    call solve_small_system(a(1:hist_count, 1:hist_count), &
                        rhs(1:hist_count), gamma(1:hist_count), linear_ok)
                    if (linear_ok) then
                        do j = 1, hist_count
                            solve_vec = solve_vec + gamma(j) * &
                                (hist_dx(:, j) + alpha_use * hist_df(:, j))
                        end do
                    else
                        hist_count = 0
                    end if
                end if
            case ('linearmixing')
                solve_vec = -alpha_use * fvec
            case ('diagbroyden')
                solve_vec = -fvec / d
            case ('excitingmixing')
                solve_vec = -fvec * beta
            end select

            if (.not. all(ieee_is_finite(solve_vec))) then
                if (trim(meth) == 'broyden1' .or. trim(meth) == 'broyden2') then
                    call set_scaled_identity(h, -alpha_use)
                    solve_vec = matmul(h, fvec)
                end if
            end if
            if (.not. all(ieee_is_finite(solve_vec))) then
                message = 'Jacobian inversion produced nonfinite values'
                return
            end if

            dx = -solve_vec
            if (same_real(norm2(dx), 0.0_dp)) then
                message = 'Jacobian inversion yielded zero vector'
                return
            end if

            if (use_armijo) then
                call armijo_step(fun, x, fvec, dx, xnew, fnew, nfev, step_status)
                if (step_status /= 0) then
                    message = 'nonfinite residual encountered'
                    return
                end if
            else if (use_wolfe) then
                call wolfe_step(fun, x, fvec, dx, xnew, fnew, nfev, step_status)
                if (step_status /= 0) then
                    message = 'nonfinite residual encountered'
                    return
                end if
            else
                xnew = x + dx
                call fun(xnew, fnew)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(fnew))) then
                    message = 'nonfinite residual encountered'
                    return
                end if
            end if

            dx_actual = xnew - x
            df = fnew - fvec
            dx_norm = norm2(dx_actual)
            df_norm = norm2(df)

            if (rank_control .and. .not. broyden_collapsed .and. &
                    (trim(meth) == 'broyden1' .or. trim(meth) == 'broyden2')) then
                call reduce_broyden_matrix(h, -alpha_use, broyden_c, broyden_d, &
                    broyden_count, max_rank_use, reduction, linear_ok)
                if (.not. linear_ok) then
                    message = 'Broyden rank reduction failed'
                    return
                end if
            end if

            select case (trim(meth))
            case ('broyden1')
                v = matmul(transpose(h), dx_actual)
                c = dx_actual - matmul(h, df)
                denom = dot_product(df, v)
                if (same_real(denom, 0.0_dp) .or. .not. ieee_is_finite(denom)) then
                    call set_scaled_identity(h, -alpha_use)
                    broyden_count = 0
                    broyden_collapsed = .false.
                else
                    rank_d = v / denom
                    do j = 1, n
                        h(:, j) = h(:, j) + c * rank_d(j)
                    end do
                    if (rank_control .and. .not. broyden_collapsed) then
                        call append_broyden_update(broyden_c, broyden_d, broyden_count, &
                            c, rank_d, broyden_collapsed)
                    end if
                end if
            case ('broyden2')
                c = dx_actual - matmul(h, df)
                denom = df_norm * df_norm
                if (same_real(denom, 0.0_dp) .or. .not. ieee_is_finite(denom)) then
                    call set_scaled_identity(h, -alpha_use)
                    broyden_count = 0
                    broyden_collapsed = .false.
                else
                    rank_d = df / denom
                    do j = 1, n
                        h(:, j) = h(:, j) + c * rank_d(j)
                    end do
                    if (rank_control .and. .not. broyden_collapsed) then
                        call append_broyden_update(broyden_c, broyden_d, broyden_count, &
                            c, rank_d, broyden_collapsed)
                    end if
                end if
            case ('anderson')
                if (anderson_m_use > 0) then
                    call append_history(hist_dx, hist_df, hist_count, dx_actual, df)
                    call build_anderson_matrix(hist_df, hist_count, w0_use, a)
                end if
            case ('diagbroyden')
                denom = dx_norm * dx_norm
                if (same_real(denom, 0.0_dp)) then
                    message = 'zero step in diagonal Broyden update'
                    return
                end if
                d = d - (df + d * dx_actual) * dx_actual / denom
            case ('excitingmixing')
                do j = 1, n
                    if (fnew(j) * fvec(j) > 0.0_dp) then
                        beta(j) = beta(j) + alpha_use
                    else
                        beta(j) = alpha_use
                    end if
                    beta(j) = min(alphamax_use, max(0.0_dp, beta(j)))
                end do
            case ('linearmixing')
            end select

            x = xnew
            fvec = fnew
            dx = dx_actual
        end do

        nit = iter_count
        if (status == 1) then
            success = .true.
            message = 'A solution was found at the specified tolerance.'
        else if (status == -1) then
            status = 0
            message = 'invalid tolerance norm value'
        else
            status = 2
            success = .false.
            message = 'The maximum number of iterations allowed has been reached.'
        end if
    end subroutine solve_nonlin_root

    subroutine solve_krylov_root(fun, x0, x, fvec, nfev, nit, status, success, &
            message, tol, ftol, fatol, xtol, xatol, maxiter, iterations, &
            line_search, rdiff, inner_method, inner_maxiter, outer_k, tol_norm, &
            inner_m, inner_m_update, inner_solver, inner_verify_residual, &
            inner_tfqmr_right_preconditioned, jacobian_diff_scheme, &
            jacobian_step_mode, jacobian_product, &
            inner_check_residual, inner_orthogonalization, outer_prune_tol, &
            outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, outer_image_basis_reuse, &
            outer_candidate_cache_max_mib, inner_restart, &
            inner_restart_adaptive, inner_restart_max, inner_m_diag, &
            inner_m_diag_auto, inner_m_diag_update_interval, inner_m_diag_floor, &
            inner_m_diag_probes, inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed, inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol, inner_m_diag_probe_max, &
            inner_m_diag_rank_one, inner_m_diag_rank_one_right, &
            inner_m_diag_rank_two, inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch, inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two_discover, inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base, jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive, inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth, inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost, inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost, inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock)
        procedure(nonlin_function) :: fun !! square residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        real(dp), intent(out) :: x(:) !! final solution estimate
        real(dp), intent(out) :: fvec(:) !! residual vector at x
        integer, intent(out) :: nfev !! residual evaluations
        integer, intent(out) :: nit !! termination-check iteration count
        integer, intent(out) :: status !! 1 success, 2 iteration limit, 0 invalid/failure
        logical, intent(out) :: success !! true when status is 1
        character(len=*), intent(out) :: message !! termination description
        real(dp), intent(in), optional :: tol !! generic root tolerance
        real(dp), intent(in), optional :: ftol !! relative residual tolerance
        real(dp), intent(in), optional :: fatol !! absolute max-norm residual tolerance
        real(dp), intent(in), optional :: xtol !! relative step tolerance
        real(dp), intent(in), optional :: xatol !! absolute max-norm step tolerance
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        integer, intent(in), optional :: iterations !! exact number of iterations to make
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: rdiff !! relative numerical-Jacobian step
        character(len=*), intent(in), optional :: inner_method !! lgmres, gmres, bicgstab, cgs, minres, or tfqmr
        integer, intent(in), optional :: inner_maxiter !! Krylov basis size
        integer, intent(in), optional :: outer_k !! retained LGMRES augmentation vectors
        procedure(nonlin_norm_function), optional :: tol_norm !! stopping-test norm, default max norm
        procedure(nonlin_linear_operator), optional :: inner_m !! approximate inverse-Jacobian action
        procedure(nonlin_preconditioner_update), optional :: inner_m_update !! update after accepted step
        procedure(nonlin_inner_solver), optional :: inner_solver !! custom matrix-free inner solve
        logical, intent(in), optional :: inner_verify_residual !! verify original TFQMR residual
        logical, intent(in), optional :: inner_tfqmr_right_preconditioned !! use consistent right TFQMR
        character(len=*), intent(in), optional :: jacobian_diff_scheme !! 2-point (default) or 3-point
        character(len=*), intent(in), optional :: jacobian_step_mode !! global (default), adaptive, or componentwise
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! analytic matrix-free J(x)*vec
        logical, intent(in), optional :: inner_check_residual !! certify the true inner residual before accepting a step
        character(len=*), intent(in), optional :: inner_orthogonalization !! once, selective, or twice for GMRES/LGMRES
        real(dp), intent(in), optional :: outer_prune_tol !! drop nearly dependent retained LGMRES vectors
        real(dp), intent(in), optional :: outer_image_prune_tol !! prune in current Jacobian-image space
        logical, intent(in), optional :: outer_image_prune_auto !! infer image-space pruning tolerance
        logical, intent(in), optional :: outer_image_reuse !! reuse current-Jacobian images in Arnoldi
        logical, intent(in), optional :: outer_candidate_image_reuse !! reconstruct candidate image from Arnoldi products
        logical, intent(in), optional :: outer_fifo_prune_fast !! enable full-history solution-space QR downdates
        logical, intent(in), optional :: outer_image_basis_reuse !! reuse image-filter basis for candidate insertion
        integer, intent(in), optional :: outer_candidate_cache_max_mib !! max raw Jacobian image cache MiB
        integer, intent(in), optional :: inner_restart !! GMRES/LGMRES restart width; inner_maxiter is total budget
        logical, intent(in), optional :: inner_restart_adaptive !! grow restarted LGMRES width on weak progress
        integer, intent(in), optional :: inner_restart_max !! maximum adaptive restart width; <= inner_maxiter
        real(dp), intent(in), optional :: inner_m_diag(:) !! diagonal of inverse-Jacobian preconditioner M
        logical, intent(in), optional :: inner_m_diag_auto !! estimate a diagonal inverse Jacobian (opt-in)
        integer, intent(in), optional :: inner_m_diag_update_interval !! nonlinear steps between updates; 0 initial only
        real(dp), intent(in), optional :: inner_m_diag_floor !! relative floor for small Jacobian diagonal entries
        integer, intent(in), optional :: inner_m_diag_probes !! 0 (default) exact, or number of Walsh probes
        character(len=*), intent(in), optional :: inner_m_diag_probe_scheme !! walsh (default), rademacher
        integer, intent(in), optional :: inner_m_diag_probe_seed !! reproducible Rademacher seed
        integer, intent(in), optional :: inner_m_diag_probe_max !! cap for progressive verified Rademacher probing
        integer, intent(in), optional :: inner_m_diag_verify_samples !! coordinate checks before use
        real(dp), intent(in), optional :: inner_m_diag_verify_rtol !! relative spot-check discrepancy limit
        logical, intent(in), optional :: inner_m_diag_rank_one !! checked mean-field rank-one correction
        real(dp), intent(in), optional :: inner_m_diag_rank_one_right(:) !! fixed right coupling weights
        logical, intent(in), optional :: inner_m_diag_rank_one_discover !! learn unknown right weights
        logical, intent(in), optional :: inner_m_diag_rank_two !! checked rank-two Woodbury correction
        real(dp), intent(in), optional :: inner_m_diag_rank_two_right(:,:) !! known n-by-2 right factors
        logical, intent(in), optional :: inner_m_diag_rank_two_sketch !! opt-in 3-probe fit before exact fallback
        logical, intent(in), optional :: inner_m_diag_rank_two_discover !! learn both unknown rank-two factors
        logical, intent(in), optional :: inner_m_diag_rank_two_stream !! stream unknown-factor discovery in O(n) space
        real(dp), intent(in), optional :: inner_m_diag_rank_two_symmetric_base(:) !! known D for symmetric rank-two L
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_base(:) !! known D for general rank-two L
        procedure(nonlin_jacobian_product), optional :: jacobian_transpose_product !! analytic transpose J(x)^T*v
        logical, intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_verify_transpose !! check two adjoint directions
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost !! J^T*v cost / Jv cost (> 0)
        procedure(nonlin_base_diagonal), optional :: inner_m_diag_rank_two_nonsymmetric_base_fun !! D(x) callback
        logical, intent(in), optional :: inner_m_diag_update_adaptive !! refresh early on weak residual reduction
        real(dp), intent(in), optional :: inner_m_diag_update_stall_ratio !! poor-step ratio threshold in (0,1]
        real(dp), intent(in), optional :: inner_m_diag_update_work_growth !! work-growth factor above one
        integer, intent(in), optional :: inner_m_diag_update_rebuild_cost !! estimated rebuild cost in inner actions
        logical, intent(in), optional :: inner_m_diag_update_measured_cost !! infer setup work from product counts
        real(dp), intent(in), optional :: inner_m_diag_update_transpose_cost !! cost of J^T*v relative to Jv
        logical, intent(in), optional :: inner_m_diag_update_timed_cost !! use CPU time for refresh decisions
        real(dp), intent(in), optional :: inner_m_diag_update_timed_min_seconds !! minimum reliable timing interval
        logical, intent(in), optional :: inner_m_diag_update_wall_clock !! wall instead of CPU time

        type(fgmres_workspace) :: gmres_work
        type(lgmres_image_arrays) :: image_arrays
        real(dp), allocatable :: dx(:)
        real(dp), allocatable :: dx_actual(:)
        real(dp), allocatable :: fnew(:)
        real(dp), allocatable :: outer_vec(:, :)
        real(dp), allocatable :: solve_vec(:)
        real(dp), allocatable :: linear_product(:)
        real(dp), allocatable :: diag_preconditioner(:)
        real(dp), allocatable :: rank_base(:), rank_correction(:), raw_jac_diag(:), rank_right(:)
        real(dp), allocatable :: rank_one_jac(:,:), rank_two_jac(:,:)
        real(dp), allocatable :: rank_two_base(:), rank_two_correction(:,:), rank_two_right(:,:)
        real(dp), allocatable :: rank_two_dynamic_base(:)
        real(dp), allocatable :: xnew(:)
        real(dp) :: eta
        real(dp) :: eta_a
        real(dp) :: inner_rtol_use
        real(dp) :: diag_floor_use
        real(dp) :: diag_verify_rtol_use
        real(dp) :: norm_ratio
        real(dp) :: outer_prune_tol_use
        real(dp) :: outer_image_prune_tol_use
        real(dp), parameter :: eta_max = 0.9999_dp
        real(dp), parameter :: eta_threshold = 0.1_dp
        real(dp) :: f0_norm
        real(dp) :: f_abs_tol
        real(dp) :: f_rel_tol
        real(dp) :: fx_norm
        real(dp) :: fx_norm_new
        real(dp), parameter :: gamma = 0.9_dp
        real(dp) :: inner_tol
        real(dp) :: omega
        real(dp) :: rdiff_use
        real(dp) :: x_abs_tol
        real(dp) :: x_rel_tol
        integer :: i
        integer :: jac_step_mode
        integer :: inner_maxiter_use
        integer :: inner_restart_use
        integer :: inner_restart_max_use
        integer :: diag_update_interval_use
        integer :: diag_last_refresh_iter
        integer :: diag_inner_actions, diag_work_baseline, diag_rebuild_cost_use
        integer :: diag_refresh_forward, diag_refresh_transpose
        real(dp) :: diag_measured_rebuild, diag_transpose_cost_use
        logical :: diag_measure_cost_use
        logical :: diag_timed_cost_use, diag_setup_clock_valid, diag_solve_clock_valid
        logical :: diag_wall_clock_use
        integer(int64) :: diag_wall_start, diag_wall_end, diag_wall_rate
        real(dp) :: diag_setup_seconds, diag_solve_seconds, diag_solve_baseline
        real(dp) :: diag_clock_start, diag_clock_end, diag_timed_min_seconds_use
        logical :: diag_work_use
        real(dp) :: diag_work_growth_use
        logical :: diag_adaptive_use, diag_refresh_pending, diag_refresh_now
        real(dp) :: diag_stall_ratio_use
        integer :: diag_probes_use
        integer :: diag_probe_seed_use
        integer :: diag_probe_max_use
        integer :: diag_verify_use
        character(len=16) :: diag_probe_scheme_use
        integer :: allocation_status
        integer :: candidate_cache_mib_use
        integer :: orth_mode
        integer :: iter_count
        integer :: maxiter_use
        integer :: n
        integer :: outer_count
        integer :: outer_k_use
        integer :: step_status
        logical :: f_abs_active
        logical :: f_rel_active
        logical :: inner_ok
        logical :: outer_image_prune_auto_use
        logical :: outer_image_reuse_use
        logical :: outer_candidate_image_reuse_use
        logical :: outer_fifo_prune_fast_use
        logical :: outer_image_basis_reuse_use
        logical :: check_inner
        logical :: inner_restart_adaptive_use
        logical :: diag_auto_use
        logical :: rank_one_use, rank_enabled, rank_one_discover
        logical :: rank_two_use, rank_two_enabled, rank_two_sketch, sketch_available
        logical :: rank_two_discover, rank_two_stream, rank_two_symmetric, rank_two_nonsymmetric
        logical :: rank_two_transpose_verify
        real(dp) :: rank_two_adaptive_cost
        logical :: use_central
        logical :: use_armijo
        logical :: use_wolfe
        logical :: x_abs_active
        logical :: x_rel_active
        character(len=16) :: inner
        character(len=16) :: orth_name
        character(len=16) :: search
        character(len=16) :: diff_scheme
        character(len=16) :: step_mode_name

        n = size(x0)
        x = 0.0_dp
        fvec = 0.0_dp
        nfev = 0
        nit = 0
        status = 0
        success = .false.
        message = 'improper input parameters'
        if (n <= 0 .or. size(x) /= n .or. size(fvec) /= n) return
        if (.not. all(ieee_is_finite(x0))) return

        ! Automatic diagonal estimation is deliberately opt-in.  Reject
        ! incompatible hooks and invalid controls before the first residual.
        diag_auto_use = .false.
        if (present(inner_m_diag_auto)) diag_auto_use = inner_m_diag_auto
        rank_one_use = .false.
        rank_enabled = .false.
        if (present(inner_m_diag_rank_one)) rank_one_use = inner_m_diag_rank_one
        rank_one_discover = .false.
        if (present(inner_m_diag_rank_one_discover)) then
            rank_one_discover = inner_m_diag_rank_one_discover
            if (rank_one_discover .and. .not. rank_one_use) then
                message = 'inner_m_diag_rank_one_discover requires inner_m_diag_rank_one'
                return
            end if
        end if
        rank_two_use = .false.
        rank_two_enabled = .false.
        if (present(inner_m_diag_rank_two)) rank_two_use = inner_m_diag_rank_two
        rank_two_discover = .false.
        if (present(inner_m_diag_rank_two_discover)) then
            rank_two_discover = inner_m_diag_rank_two_discover
            if (rank_two_discover .and. .not. rank_two_use) then
                message = 'inner_m_diag_rank_two_discover requires inner_m_diag_rank_two'
                return
            end if
        end if
        rank_two_stream = .false.
        if (present(inner_m_diag_rank_two_stream)) then
            rank_two_stream = inner_m_diag_rank_two_stream
            if (rank_two_stream .and. .not. rank_two_discover) then
                message = 'inner_m_diag_rank_two_stream requires rank-two discovery'
                return
            end if
        end if
        rank_two_transpose_verify = .false.
        if (present(inner_m_diag_rank_two_nonsymmetric_verify_transpose)) then
            rank_two_transpose_verify = inner_m_diag_rank_two_nonsymmetric_verify_transpose
            if (.not. (present(inner_m_diag_rank_two_nonsymmetric_base) .or. &
                    present(inner_m_diag_rank_two_nonsymmetric_base_fun))) then
                message = 'transpose verification requires nonsymmetric rank-two base'
                return
            end if
        end if
        rank_two_adaptive_cost = 0.0_dp
        if (present(inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost)) then
            if (.not. (present(inner_m_diag_rank_two_nonsymmetric_base) .or. &
                    present(inner_m_diag_rank_two_nonsymmetric_base_fun))) then
                message = 'adaptive transpose cost requires nonsymmetric rank-two base'
                return
            end if
            if (present(inner_m_diag_rank_two_nonsymmetric_verify_transpose)) then
                message = 'adaptive transpose cost conflicts with fixed transpose checks'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost)) then
                message = 'adaptive transpose cost must be positive and finite'
                return
            end if
            if (inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost <= 0.0_dp) then
                message = 'adaptive transpose cost must be positive and finite'
                return
            end if
            rank_two_adaptive_cost = inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost
        end if
        if (present(inner_m_diag_rank_two_nonsymmetric_base) .and. &
                present(inner_m_diag_rank_two_nonsymmetric_base_fun)) then
            message = 'nonsymmetric rank-two base array and callback are exclusive'
            return
        end if
        rank_two_nonsymmetric = present(inner_m_diag_rank_two_nonsymmetric_base) .or. &
            present(inner_m_diag_rank_two_nonsymmetric_base_fun)
        if (present(jacobian_transpose_product) .neqv. rank_two_nonsymmetric) then
            message = 'transpose product and nonsymmetric base must be supplied together'
            return
        end if
        if (rank_two_nonsymmetric) then
            if (.not. diag_auto_use .or. .not. rank_two_discover .or. rank_two_stream) then
                message = 'nonsymmetric rank-two base requires auto discovery without streaming'
                return
            end if
            if (present(inner_m_diag_rank_two_nonsymmetric_base)) then
                if (size(inner_m_diag_rank_two_nonsymmetric_base) /= n) then
                    message = 'nonsymmetric rank-two base must have size(x0) entries'
                    return
                end if
                if (.not. all(ieee_is_finite(inner_m_diag_rank_two_nonsymmetric_base))) then
                    message = 'nonsymmetric rank-two base must be finite and nonsingular'
                    return
                end if
                if (any(abs(inner_m_diag_rank_two_nonsymmetric_base) <= sqrt(tiny(1.0_dp)))) then
                    message = 'nonsymmetric rank-two base must be finite and nonsingular'
                    return
                end if
            end if
        end if
        rank_two_symmetric = present(inner_m_diag_rank_two_symmetric_base)
        if (rank_two_symmetric) then
            if (rank_two_nonsymmetric) then
                message = 'symmetric and nonsymmetric rank-two bases are exclusive'
                return
            end if
            if (.not. diag_auto_use .or. .not. rank_two_discover .or. rank_two_stream) then
                message = 'symmetric rank-two base requires auto discovery without streaming'
                return
            end if
            if (size(inner_m_diag_rank_two_symmetric_base) /= n) then
                message = 'symmetric rank-two base must have size(x0) entries'
                return
            end if
            if (.not. all(ieee_is_finite(inner_m_diag_rank_two_symmetric_base))) then
                message = 'symmetric rank-two base must be finite and nonsingular'
                return
            end if
            if (any(abs(inner_m_diag_rank_two_symmetric_base) <= sqrt(tiny(1.0_dp)))) then
                message = 'symmetric rank-two base must be finite and nonsingular'
                return
            end if
        end if
        rank_two_sketch = .false.
        if (present(inner_m_diag_rank_two_sketch)) then
            rank_two_sketch = inner_m_diag_rank_two_sketch
            if (rank_two_sketch .and. .not. rank_two_use) then
                message = 'inner_m_diag_rank_two_sketch requires inner_m_diag_rank_two'
                return
            end if
        end if
        diag_adaptive_use = .false.
        if (present(inner_m_diag_update_adaptive)) then
            diag_adaptive_use = inner_m_diag_update_adaptive
            if (diag_adaptive_use .and. .not. diag_auto_use) then
                message = 'inner_m_diag_update_adaptive requires inner_m_diag_auto'
                return
            end if
        end if
        diag_stall_ratio_use = 0.8_dp
        if (present(inner_m_diag_update_stall_ratio)) then
            if (.not. diag_adaptive_use) then
                message = 'inner_m_diag_update_stall_ratio requires adaptive refresh'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_update_stall_ratio)) then
                message = 'inner_m_diag_update_stall_ratio must be finite and in (0,1]'
                return
            end if
            if (inner_m_diag_update_stall_ratio <= 0.0_dp .or. &
                    inner_m_diag_update_stall_ratio > 1.0_dp) then
                message = 'inner_m_diag_update_stall_ratio must be finite and in (0,1]'
                return
            end if
            diag_stall_ratio_use = inner_m_diag_update_stall_ratio
        end if
        ! Optional work-sensitive scheduling compares the inner solver's
        ! preconditioner applications with the count at the last refresh.
        ! For GMRES, an application corresponds to an Arnoldi step; for
        ! other/custom solvers it is an inexpensive proxy for linear work.
        diag_work_use = present(inner_m_diag_update_work_growth)
        diag_work_growth_use = 2.0_dp
        diag_rebuild_cost_use = 0
        if (diag_work_use) then
            if (.not. diag_adaptive_use) then
                message = 'inner_m_diag_update_work_growth requires adaptive refresh'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_update_work_growth)) then
                message = 'inner_m_diag_update_work_growth must be finite and greater than one'
                return
            end if
            if (inner_m_diag_update_work_growth <= 1.0_dp) then
                message = 'inner_m_diag_update_work_growth must be finite and greater than one'
                return
            end if
            diag_work_growth_use = inner_m_diag_update_work_growth
        end if
        if (present(inner_m_diag_update_rebuild_cost)) then
            if (.not. diag_work_use .or. inner_m_diag_update_rebuild_cost < 0) then
                message = 'inner_m_diag_update_rebuild_cost requires work growth and must be nonnegative'
                return
            end if
            diag_rebuild_cost_use = inner_m_diag_update_rebuild_cost
        end if
        ! Count attempted setup Jacobian products in every automatic
        ! preconditioner, including rejected sketches and exact fallback.
        ! This is opt-in and cannot override a supplied cost.
        diag_measure_cost_use = .false.
        if (present(inner_m_diag_update_measured_cost)) then
            diag_measure_cost_use = inner_m_diag_update_measured_cost
            if (diag_measure_cost_use .and. (.not. diag_work_use .or. &
                    .not. diag_auto_use .or. &
                    present(inner_m_diag_update_rebuild_cost))) then
                message = 'measured refresh cost requires auto and work growth, without manual cost'
                return
            end if
        end if
        diag_transpose_cost_use = 1.0_dp
        if (present(inner_m_diag_update_transpose_cost)) then
            if (.not. diag_measure_cost_use) then
                message = 'transpose refresh cost requires measured refresh cost'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_update_transpose_cost)) then
                message = 'transpose refresh cost must be positive and finite'
                return
            end if
            if (inner_m_diag_update_transpose_cost <= 0.0_dp) then
                message = 'transpose refresh cost must be positive and finite'
                return
            end if
            diag_transpose_cost_use = inner_m_diag_update_transpose_cost
        end if
        ! Alternative to the product-count cost model: compare actual CPU
        ! seconds used by the inner solve against the last post-refresh solve
        ! and the preceding setup time. Never infer a refresh from a clock
        ! that is unavailable or too coarse to measure these operations.
        diag_timed_cost_use = .false.
        if (present(inner_m_diag_update_timed_cost)) then
            diag_timed_cost_use = inner_m_diag_update_timed_cost
            if (diag_timed_cost_use .and. (.not. diag_work_use .or. &
                    .not. diag_auto_use .or. diag_measure_cost_use .or. &
                    present(inner_m_diag_update_rebuild_cost) .or. &
                    present(inner_m_diag_update_transpose_cost))) then
                message = 'timed refresh cost requires auto and work growth, without another cost model'
                return
            end if
        end if
        diag_timed_min_seconds_use = 0.001_dp
        if (present(inner_m_diag_update_timed_min_seconds)) then
            if (.not. diag_timed_cost_use) then
                message = 'timed minimum seconds requires timed refresh cost'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_update_timed_min_seconds)) then
                message = 'timed minimum seconds must be finite and nonnegative'
                return
            end if
            if (inner_m_diag_update_timed_min_seconds < 0.0_dp) then
                message = 'timed minimum seconds must be finite and nonnegative'
                return
            end if
            diag_timed_min_seconds_use = inner_m_diag_update_timed_min_seconds
        end if
        diag_wall_clock_use = .false.
        if (present(inner_m_diag_update_wall_clock)) then
            diag_wall_clock_use = inner_m_diag_update_wall_clock
            if (diag_wall_clock_use .and. .not. diag_timed_cost_use) then
                message = 'wall-clock refresh requires timed refresh cost'
                return
            end if
        end if
        ! Adaptive mode defaults to a four-step maximum age.  Fixed mode
        ! continues to refresh on every iteration by default.
        diag_update_interval_use = 1
        if (diag_adaptive_use) diag_update_interval_use = 4
        diag_probes_use = 0
        diag_probe_seed_use = 1
        diag_probe_max_use = 0
        diag_verify_use = 0
        diag_verify_rtol_use = 0.25_dp
        diag_probe_scheme_use = 'walsh'
        diag_floor_use = sqrt(epsilon(1.0_dp))
        if (present(inner_m_diag_probes)) then
            if (.not. diag_auto_use .or. inner_m_diag_probes < 0) then
                message = 'inner_m_diag_probes requires auto and must be nonnegative'
                return
            end if
            diag_probes_use = inner_m_diag_probes
        end if
        if (present(inner_m_diag_probe_scheme)) then
            if (.not. diag_auto_use) then
                message = 'inner_m_diag_probe_scheme requires auto'
                return
            end if
            diag_probe_scheme_use = lowercase(trim(inner_m_diag_probe_scheme))
            if (trim(diag_probe_scheme_use) /= 'walsh' .and. &
                    trim(diag_probe_scheme_use) /= 'rademacher') then
                message = 'inner_m_diag_probe_scheme must be walsh or rademacher'
                return
            end if
        end if
        if (present(inner_m_diag_probe_seed)) then
            if (.not. diag_auto_use .or. trim(diag_probe_scheme_use) /= 'rademacher') then
                message = 'inner_m_diag_probe_seed requires rademacher auto probes'
                return
            end if
            diag_probe_seed_use = inner_m_diag_probe_seed
        end if
        if (present(inner_m_diag_verify_samples)) then
            if (.not. diag_auto_use .or. inner_m_diag_verify_samples < 0 .or. &
                    inner_m_diag_verify_samples > n) then
                message = 'inner_m_diag_verify_samples requires auto and must be in [0,n]'
                return
            end if
            diag_verify_use = inner_m_diag_verify_samples
        end if
        if (present(inner_m_diag_verify_rtol)) then
            if (.not. diag_auto_use .or. diag_verify_use <= 0) then
                message = 'inner_m_diag_verify_rtol requires positive verification samples'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_verify_rtol)) then
                message = 'inner_m_diag_verify_rtol must be finite and in (0, 1]'
                return
            end if
            if (inner_m_diag_verify_rtol <= 0.0_dp .or. inner_m_diag_verify_rtol > 1.0_dp) then
                message = 'inner_m_diag_verify_rtol must be finite and in (0, 1]'
                return
            end if
            diag_verify_rtol_use = inner_m_diag_verify_rtol
        end if
        if (present(inner_m_diag_probe_max)) then
            if (.not. diag_auto_use .or. trim(diag_probe_scheme_use) /= 'rademacher' .or. &
                    diag_verify_use <= 0 .or. diag_probes_use <= 0 .or. &
                    diag_probes_use >= n .or. inner_m_diag_probe_max < diag_probes_use .or. &
                    inner_m_diag_probe_max > n) then
                message = 'inner_m_diag_probe_max requires verified Rademacher probes in [1,n)'
                return
            end if
            diag_probe_max_use = inner_m_diag_probe_max
        end if
        if (present(inner_m_diag_update_interval)) then
            if (.not. diag_auto_use .or. inner_m_diag_update_interval < 0) then
                message = 'inner_m_diag_update_interval requires auto and must be nonnegative'
                return
            end if
            if (diag_adaptive_use .and. inner_m_diag_update_interval == 0) then
                message = 'adaptive diagonal refresh requires a positive update interval'
                return
            end if
            diag_update_interval_use = inner_m_diag_update_interval
        end if
        if (present(inner_m_diag_floor)) then
            if (.not. diag_auto_use) then
                message = 'inner_m_diag_floor requires inner_m_diag_auto'
                return
            end if
            if (.not. ieee_is_finite(inner_m_diag_floor)) then
                message = 'inner_m_diag_floor must be finite and in (0, 1]'
                return
            end if
            if (inner_m_diag_floor <= 0.0_dp .or. inner_m_diag_floor > 1.0_dp) then
                message = 'inner_m_diag_floor must be finite and in (0, 1]'
                return
            end if
            diag_floor_use = inner_m_diag_floor
        end if
        if (present(inner_m_diag_rank_one_right)) then
            if (rank_one_discover) then
                message = 'discovered and supplied rank-one right factors conflict'
                return
            end if
            if (.not. rank_one_use) then
                message = 'inner_m_diag_rank_one_right requires inner_m_diag_rank_one'
                return
            end if
            if (size(inner_m_diag_rank_one_right) /= n) then
                message = 'inner_m_diag_rank_one_right must have size(x0) entries'
                return
            end if
            if (.not. all(ieee_is_finite(inner_m_diag_rank_one_right))) then
                message = 'inner_m_diag_rank_one_right must be finite and nonzero'
                return
            end if
            if (all(inner_m_diag_rank_one_right == 0.0_dp)) then
                message = 'inner_m_diag_rank_one_right must be finite and nonzero'
                return
            end if
        end if
        if (present(inner_m_diag_rank_two_right)) then
            if (rank_two_discover) then
                message = 'discovered and supplied rank-two right factors conflict'
                return
            end if
            if (.not. rank_two_use) then
                message = 'inner_m_diag_rank_two_right requires inner_m_diag_rank_two'
                return
            end if
            if (size(inner_m_diag_rank_two_right,1) /= n .or. &
                    size(inner_m_diag_rank_two_right,2) /= 2) then
                message = 'inner_m_diag_rank_two_right must have shape (size(x0),2)'
                return
            end if
            if (.not. all(ieee_is_finite(inner_m_diag_rank_two_right))) then
                message = 'inner_m_diag_rank_two_right must have finite nonzero columns'
                return
            end if
            if (any(maxval(abs(inner_m_diag_rank_two_right),dim=1) == 0.0_dp)) then
                message = 'inner_m_diag_rank_two_right must have finite nonzero columns'
                return
            end if
        end if
        if (rank_two_discover .and. rank_two_sketch) then
            message = 'inner_m_diag_rank_two_discover conflicts with sketching'
            return
        end if
        if (rank_two_use) then
            if (rank_one_use) then
                message = 'rank-one and rank-two options cannot be combined'
                return
            end if
            if (.not. rank_two_discover .and. .not. present(inner_m_diag_rank_two_right)) then
                message = 'inner_m_diag_rank_two requires inner_m_diag_rank_two_right'
                return
            end if
            if (.not. diag_auto_use .or. diag_probes_use /= 0) then
                message = 'inner_m_diag_rank_two requires exact inner_m_diag_auto'
                return
            end if
            if (present(inner_method)) then
                if (lowercase(trim(inner_method)) == 'minres' .and. .not. present(inner_solver)) then
                    message = 'inner_m_diag_rank_two is unsupported with minres'
                    return
                end if
            end if
        end if
        if (rank_one_use) then
            if (.not. diag_auto_use .or. diag_probes_use /= 0) then
                message = 'inner_m_diag_rank_one requires exact inner_m_diag_auto'
                return
            end if
            if (present(inner_solver)) then
                ! Custom inner solvers can consume nonsymmetric corrections.
            else if (present(inner_method)) then
                if (lowercase(trim(inner_method)) == 'minres') then
                    message = 'inner_m_diag_rank_one is unsupported with minres'
                    return
                end if
            end if
        end if
        if (diag_auto_use) then
            if (present(inner_m_diag) .or. present(inner_m) .or. present(inner_m_update)) then
                message = 'inner_m_diag_auto conflicts with other preconditioners'
                return
            end if
        end if

        ! Static diagonal right preconditioning, an inexpensive alternative to
        ! the procedure callback.  Reject invalid entries before the first
        ! residual evaluation, so a failed setup has nfev=0.
        if (present(inner_m_diag)) then
            if (present(inner_m) .or. present(inner_m_update)) then
                message = 'inner_m_diag cannot be combined with inner_m or inner_m_update'
                return
            end if
            if (size(inner_m_diag) /= n) then
                message = 'inner_m_diag must have size(x0) entries'
                return
            end if
            if (.not. all(ieee_is_finite(inner_m_diag))) then
                message = 'inner_m_diag entries must be finite and nonzero'
                return
            end if
            if (any(inner_m_diag == 0.0_dp)) then
                message = 'inner_m_diag entries must be finite and nonzero'
                return
            end if
        end if

        search = 'armijo'
        if (present(line_search)) search = lowercase(trim(line_search))
        select case (trim(search))
        case ('armijo')
            use_armijo = .true.
            use_wolfe = .false.
        case ('wolfe')
            use_armijo = .false.
            use_wolfe = .true.
        case ('none', 'off', 'false')
            use_armijo = .false.
            use_wolfe = .false.
        case default
            message = 'unsupported line search; use armijo, wolfe, or none'
            return
        end select

        inner = 'lgmres'
        if (present(inner_method)) inner = lowercase(trim(inner_method))
        if (.not. present(inner_solver)) then
            select case (trim(inner))
            case ('lgmres', 'gmres', 'bicgstab', 'cgs', 'minres', 'tfqmr')
            case default
                message = 'unsupported Krylov method; use lgmres, gmres, bicgstab, cgs, minres, or tfqmr'
                return
            end select
        end if
        ! MINRES needs a positive-definite preconditioner. Other Krylov
        ! methods may accept a signed diagonal inverse-Jacobian estimate.
        if (present(inner_m_diag)) then
            if (trim(inner) == 'minres' .and. .not. present(inner_solver)) then
                if (any(inner_m_diag <= 0.0_dp)) then
                    message = 'minres requires a positive inner_m_diag'
                    return
                end if
            end if
        end if

        ! Default to one modified Gram-Schmidt pass, as in the existing
        ! SciPy-compatible implementation. Optional DGKS reorthogonalization
        ! improves the Arnoldi basis when a new vector nearly lies in its span.
        orth_mode = 0
        if (present(inner_orthogonalization)) then
            if (present(inner_solver) .or. &
                    (trim(inner) /= 'gmres' .and. trim(inner) /= 'lgmres')) then
                message = 'inner_orthogonalization requires built-in gmres or lgmres'
                return
            end if
            orth_name = lowercase(trim(inner_orthogonalization))
            select case (trim(orth_name))
            case ('once')
                orth_mode = 0
            case ('selective')
                orth_mode = 1
            case ('twice')
                orth_mode = 2
            case default
                message = 'inner_orthogonalization must be once, selective, or twice'
                return
            end select
        end if

        inner_maxiter_use = 20
        if (present(inner_maxiter)) inner_maxiter_use = inner_maxiter
        if (inner_maxiter_use < 1) return
        inner_restart_use = 0
        if (present(inner_restart)) then
            if (inner_restart < 1) then
                message = 'inner_restart must be positive'
                return
            end if
            if (present(inner_solver) .or. &
                    (trim(inner) /= 'gmres' .and. trim(inner) /= 'lgmres')) then
                message = 'inner_restart requires built-in gmres or lgmres'
                return
            end if
            if (inner_restart < inner_maxiter_use) inner_restart_use = inner_restart
        end if
        inner_restart_adaptive_use = .false.
        if (present(inner_restart_adaptive)) inner_restart_adaptive_use = inner_restart_adaptive
        inner_restart_max_use = inner_restart_use
        if (present(inner_restart_max) .and. .not. inner_restart_adaptive_use) then
            message = 'inner_restart_max requires inner_restart_adaptive'
            return
        end if
        if (inner_restart_adaptive_use) then
            if (present(inner_solver) .or. trim(inner) /= 'lgmres' .or. &
                    inner_restart_use == 0) then
                message = 'inner_restart_adaptive requires restarted built-in lgmres'
                return
            end if
            ! Limit default adaptive growth to four times the original width.
            ! Use int64 before multiplication to avoid integer overflow.
            inner_restart_max_use = int(min(int(inner_maxiter_use, int64), &
                4_int64 * int(inner_restart_use, int64)))
            if (present(inner_restart_max)) then
                if (inner_restart_max < inner_restart_use .or. &
                        inner_restart_max > inner_maxiter_use) then
                    message = 'inner_restart_max must be between inner_restart and inner_maxiter'
                    return
                end if
                inner_restart_max_use = inner_restart_max
            end if
        end if
        outer_k_use = 10
        if (present(outer_k)) outer_k_use = outer_k
        if (outer_k_use < 0) return
        if (trim(inner) /= 'lgmres') outer_k_use = 0
        outer_prune_tol_use = 0.0_dp
        if (present(outer_prune_tol)) then
            if (.not. ieee_is_finite(outer_prune_tol)) then
                message = 'outer_prune_tol must be finite and in [0, 1)'
                return
            end if
            if (outer_prune_tol < 0.0_dp .or. outer_prune_tol >= 1.0_dp) then
                message = 'outer_prune_tol must be finite and in [0, 1)'
                return
            end if
            if (outer_prune_tol > 0.0_dp) then
                if (trim(inner) /= 'lgmres' .or. present(inner_solver) .or. outer_k_use == 0) then
                    message = 'outer_prune_tol requires built-in lgmres with outer_k > 0'
                    return
                end if
            end if
            outer_prune_tol_use = outer_prune_tol
        end if
        outer_fifo_prune_fast_use = .false.
        if (present(outer_fifo_prune_fast)) then
            outer_fifo_prune_fast_use = outer_fifo_prune_fast
            if (outer_fifo_prune_fast_use) then
                if (trim(inner) /= 'lgmres' .or. present(inner_solver) .or. &
                        outer_k_use < 2 .or. outer_prune_tol_use <= 0.0_dp) then
                    message = 'outer_fifo_prune_fast requires built-in lgmres, outer_prune_tol > 0, outer_k >= 2'
                    return
                end if
            end if
        end if
        outer_image_prune_tol_use = 0.0_dp
        if (present(outer_image_prune_tol)) then
            if (.not. ieee_is_finite(outer_image_prune_tol)) then
                message = 'outer_image_prune_tol must be finite and in [0, 1)'
                return
            end if
            if (outer_image_prune_tol < 0.0_dp .or. outer_image_prune_tol >= 1.0_dp) then
                message = 'outer_image_prune_tol must be finite and in [0, 1)'
                return
            end if
            if (outer_image_prune_tol > 0.0_dp) then
                if (trim(inner) /= 'lgmres' .or. present(inner_solver) .or. outer_k_use == 0) then
                    message = 'outer_image_prune_tol requires built-in lgmres with outer_k > 0'
                    return
                end if
                if (outer_prune_tol_use > 0.0_dp) then
                    message = 'outer_image_prune_tol cannot be combined with outer_prune_tol'
                    return
                end if
            end if
            outer_image_prune_tol_use = outer_image_prune_tol
        end if
        outer_image_prune_auto_use = .false.
        if (present(outer_image_prune_auto)) then
            outer_image_prune_auto_use = outer_image_prune_auto
            if (outer_image_prune_auto_use) then
                if (trim(inner) /= 'lgmres' .or. present(inner_solver) .or. outer_k_use == 0) then
                    message = 'outer_image_prune_auto requires built-in lgmres with outer_k > 0'
                    return
                end if
                if (outer_image_prune_tol_use > 0.0_dp .or. outer_prune_tol_use > 0.0_dp) then
                    message = 'outer_image_prune_auto cannot be combined with fixed pruning thresholds'
                    return
                end if
            end if
        end if
        outer_image_basis_reuse_use = .false.
        if (present(outer_image_basis_reuse)) then
            outer_image_basis_reuse_use = outer_image_basis_reuse
            if (outer_image_basis_reuse_use) then
                if (.not. (outer_image_prune_auto_use .or. outer_image_prune_tol_use > 0.0_dp)) then
                    message = 'outer_image_basis_reuse requires Jacobian-image pruning'
                    return
                end if
            end if
        end if
        outer_image_reuse_use = .false.
        if (present(outer_image_reuse)) then
            outer_image_reuse_use = outer_image_reuse
            if (outer_image_reuse_use) then
                if (.not. (outer_image_prune_auto_use .or. outer_image_prune_tol_use > 0.0_dp)) then
                    message = 'outer_image_reuse requires Jacobian-image pruning'
                    return
                end if
            end if
        end if
        outer_candidate_image_reuse_use = .false.
        if (present(outer_candidate_image_reuse)) then
            outer_candidate_image_reuse_use = outer_candidate_image_reuse
            if (outer_candidate_image_reuse_use) then
                if (.not. (outer_image_prune_auto_use .or. outer_image_prune_tol_use > 0.0_dp)) then
                    message = 'outer_candidate_image_reuse requires Jacobian-image pruning'
                    return
                end if
                if (.not. present(jacobian_product)) then
                    message = 'outer_candidate_image_reuse requires analytic jacobian_product'
                    return
                end if
            end if
        end if
        candidate_cache_mib_use = 0
        if (present(outer_candidate_cache_max_mib)) then
            candidate_cache_mib_use = outer_candidate_cache_max_mib
            if (candidate_cache_mib_use < 0) then
                message = 'outer_candidate_cache_max_mib must be nonnegative'
                return
            end if
            if (candidate_cache_mib_use > 0 .and. .not. outer_candidate_image_reuse_use) then
                message = 'outer_candidate_cache_max_mib requires candidate-image reuse'
                return
            end if
        end if
        if (present(inner_verify_residual)) then
            if (inner_verify_residual .and. &
                    (trim(inner) /= 'tfqmr' .or. present(inner_solver))) then
                message = 'inner_verify_residual requires the built-in tfqmr method'
                return
            end if
        end if

        if (present(inner_tfqmr_right_preconditioned)) then
            if (inner_tfqmr_right_preconditioned .and. &
                    (trim(inner) /= 'tfqmr' .or. present(inner_solver))) then
                message = 'inner_tfqmr_right_preconditioned requires the built-in tfqmr method'
                return
            end if
        end if

        ! The default remains SciPy's forward difference. A centered stencil
        ! has O(h**2) truncation error and uses a larger optimal step.
        diff_scheme = '2-point'
        if (present(jacobian_diff_scheme)) diff_scheme = lowercase(trim(jacobian_diff_scheme))
        select case (trim(diff_scheme))
        case ('2-point')
            use_central = .false.
        case ('3-point')
            use_central = .true.
        case default
            message = 'jacobian_diff_scheme must be 2-point or 3-point'
            return
        end select

        ! The SciPy-compatible global step remains the default. Adaptive mode
        ! uses per-coordinate steps only if the global step is ill-scaled or
        ! unrepresentable for an active coordinate.
        step_mode_name = 'global'
        if (present(jacobian_step_mode)) step_mode_name = lowercase(trim(jacobian_step_mode))
        select case (trim(step_mode_name))
        case ('global')
            jac_step_mode = jac_step_global
        case ('adaptive')
            jac_step_mode = jac_step_adaptive
        case ('componentwise')
            jac_step_mode = jac_step_componentwise
        case default
            message = 'jacobian_step_mode must be global, adaptive, or componentwise'
            return
        end select

        if (present(jacobian_product)) then
            if (present(rdiff) .or. present(jacobian_diff_scheme) .or. &
                    present(jacobian_step_mode)) then
                message = 'jacobian_product cannot be combined with finite-difference options'
                return
            end if
        end if

        check_inner = .false.
        if (present(inner_check_residual)) check_inner = inner_check_residual

        rdiff_use = sqrt(epsilon(1.0_dp))
        if (use_central) rdiff_use = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
        if (present(rdiff)) rdiff_use = rdiff
        if (.not. ieee_is_finite(rdiff_use) .or. rdiff_use <= 0.0_dp) return

        maxiter_use = 100 * (n + 1)
        if (present(iterations)) then
            if (iterations < 0) return
            maxiter_use = iterations + 1
        end if
        if (present(maxiter)) then
            if (maxiter < 0) return
            maxiter_use = maxiter
        end if
        if (maxiter_use < 1) then
            status = 2
            message = 'The maximum number of iterations allowed has been reached.'
            return
        end if

        call setup_tolerances(tol, ftol, fatol, xtol, xatol, f_abs_tol, &
            f_rel_tol, x_abs_tol, x_rel_tol, f_abs_active, f_rel_active, &
            x_abs_active, x_rel_active, inner_ok)
        if (.not. inner_ok) return

        allocate(dx(n), dx_actual(n), fnew(n), solve_vec(n), xnew(n))
        if (check_inner) allocate(linear_product(n))
        if (diag_auto_use) allocate(diag_preconditioner(n))
        if (rank_one_use) then
            allocate(raw_jac_diag(n), rank_base(n), rank_correction(n), rank_right(n))
            if (rank_one_discover) then
                allocate(rank_one_jac(n,n), stat=allocation_status)
                if (allocation_status /= 0) then
                    message = 'unable to allocate rank-one discovery Jacobian'
                    return
                end if
            end if
            rank_right = 1.0_dp
            if (present(inner_m_diag_rank_one_right)) then
                ! The coupling u*v^T is unchanged by normalization of v.
                ! Avoid overflow from an arbitrary user-supplied weight scale.
                rank_right = inner_m_diag_rank_one_right / &
                    maxval(abs(inner_m_diag_rank_one_right))
            end if
        end if
        if (rank_two_nonsymmetric) then
            allocate(rank_two_dynamic_base(n))
            if (present(inner_m_diag_rank_two_nonsymmetric_base)) then
                rank_two_dynamic_base = inner_m_diag_rank_two_nonsymmetric_base
            end if
        end if
        if (rank_two_use) then
            allocate(raw_jac_diag(n), rank_two_base(n), &
                rank_two_correction(n,2), rank_two_right(n,2))
            if (rank_two_discover) then
                if (.not. rank_two_stream .and. .not. rank_two_symmetric .and. &
                        .not. rank_two_nonsymmetric) then
                    allocate(rank_two_jac(n,n), stat=allocation_status)
                    if (allocation_status /= 0) then
                        message = 'unable to allocate rank-two discovery Jacobian'
                        return
                    end if
                end if
            else
                rank_two_right(:,1) = inner_m_diag_rank_two_right(:,1) / &
                    maxval(abs(inner_m_diag_rank_two_right(:,1)))
                rank_two_right(:,2) = inner_m_diag_rank_two_right(:,2) / &
                    maxval(abs(inner_m_diag_rank_two_right(:,2)))
            end if
        end if
        allocate(outer_vec(n, outer_k_use))
        outer_vec = 0.0_dp
        outer_count = 0
        call reset_lgmres_direction_pruning_stream(gmres_work%pruning)

        x = x0
        call fun(x, fvec)
        nfev = 1
        if (.not. all(ieee_is_finite(fvec))) then
            message = 'nonfinite residual encountered'
            return
        end if

        dx = huge(1.0_dp)
        f0_norm = maxval(abs(fvec))
        if (present(tol_norm)) f0_norm = tol_norm(fvec)
        if (.not. ieee_is_finite(f0_norm)) then
            message = 'invalid tolerance norm value'
            return
        end if
        if (f0_norm < 0.0_dp) then
            message = 'invalid tolerance norm value'
            return
        end if
        fx_norm = safe_two_norm(fvec)
        eta = 1.0e-3_dp
        iter_count = 0
        diag_last_refresh_iter = 0
        diag_refresh_pending = .false.

        diag_inner_actions = 0
        diag_work_baseline = 0
        diag_measured_rebuild = 0.0_dp
        diag_setup_seconds = 0.0_dp
        diag_solve_seconds = 0.0_dp
        diag_solve_baseline = 0.0_dp
        diag_setup_clock_valid = .false.
        diag_solve_clock_valid = .false.
        do i = 1, maxiter_use
            iter_count = iter_count + 1
            if (termination_met(fvec, x, dx, f0_norm, f_abs_tol, f_rel_tol, &
                    x_abs_tol, x_rel_tol, f_abs_active, f_rel_active, &
                    x_abs_active, x_rel_active, iterations, iter_count, status, &
                    tol_norm)) then
                exit
            end if

            if (present(jacobian_product)) then
                ! No perturbations or finite-difference scaling are needed.
                omega = 0.0_dp
            else if (jac_step_mode == jac_step_componentwise) then
                ! The coordinatewise operator never needs a global scale.
                omega = rdiff_use
            else
                if (rdiff_use > huge(1.0_dp) / max(1.0_dp, maxval(abs(x)))) then
                    message = 'Krylov difference scale exceeds representable range'
                    return
                end if
                omega = rdiff_use * max(1.0_dp, maxval(abs(x))) / &
                    max(1.0_dp, maxval(abs(fvec)))
            end if
            if (diag_auto_use) then
                if (diag_adaptive_use) then
                    ! Never exceed the user-specified maximum refresh age.
                    ! A weak accepted step can cause an earlier refresh, but
                    ! a successful step leaves the current inverse cached.
                    diag_refresh_now = i == 1 .or. diag_refresh_pending .or. &
                        i - diag_last_refresh_iter >= diag_update_interval_use
                else
                    diag_refresh_now = i == 1 .or. (diag_update_interval_use > 0 .and. &
                        mod(i - 1, max(1, diag_update_interval_use)) == 0)
                end if
                if (diag_refresh_now) then
                    if (diag_timed_cost_use) call start_diag_timer()
                    diag_refresh_forward = 0
                    diag_refresh_transpose = 0
                    diag_last_refresh_iter = i
                    diag_refresh_pending = .false.
                    ! Sketch-based rank-two fitting can avoid the n exact
                    ! coordinate products. Verify the fit before using it; on
                    ! rejection, recover the original exact-diagonal path.
                    rank_two_enabled = .false.
                    if (rank_two_use .and. rank_two_sketch) then
                        call estimate_rank_two_diagonal_sketch(fun, x, fvec, omega, &
                            nfev, use_central, jac_step_mode, rdiff_use, &
                            rank_two_right, raw_jac_diag, sketch_available, inner_ok, &
                            jacobian_product, &
                            forward_calls=diag_refresh_forward)
                        if (.not. inner_ok) then
                            message = 'rank-two sketch Jacobian evaluation failed'
                            return
                        end if
                        if (sketch_available) then
                            call estimate_rank_two_inverse(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, raw_jac_diag, &
                                rank_two_right, rank_two_base, rank_two_correction, &
                                rank_two_enabled, inner_ok, jacobian_product, &
                                forward_calls=diag_refresh_forward)
                            if (.not. inner_ok) then
                                message = 'rank-two preconditioner verification failed'
                                return
                            end if
                        end if
                    end if
                    if (rank_two_nonsymmetric) then
                        if (present(inner_m_diag_rank_two_nonsymmetric_base_fun)) then
                            ! The supplied D can vary with x; refresh it at the
                            ! same point as J(x), never at a trial line-search point.
                            call inner_m_diag_rank_two_nonsymmetric_base_fun(x, &
                                rank_two_dynamic_base, inner_ok)
                            if (.not. inner_ok) then
                                message = 'nonsymmetric rank-two base callback failed'
                                return
                            end if
                            if (.not. all(ieee_is_finite(rank_two_dynamic_base))) then
                                message = 'nonsymmetric rank-two base callback returned invalid diagonal'
                                return
                            end if
                            if (any(abs(rank_two_dynamic_base) <= sqrt(tiny(1.0_dp)))) then
                                message = 'nonsymmetric rank-two base callback returned invalid diagonal'
                                return
                            end if
                        end if
                        call discover_rank_two_inverse_nonsymmetric(fun, x, fvec, omega, nfev, &
                            use_central, jac_step_mode, rdiff_use, &
                            rank_two_dynamic_base, rank_two_right, &
                            rank_two_base, rank_two_correction, rank_two_enabled, &
                            inner_ok, jacobian_transpose_product, jacobian_product, &
                            verify_transpose=rank_two_transpose_verify, &
                            adaptive_transpose_cost=rank_two_adaptive_cost, &
                            forward_calls=diag_refresh_forward, &
                            transpose_calls=diag_refresh_transpose)
                        if (.not. inner_ok) then
                            message = 'nonsymmetric rank-two Jacobian evaluation failed'
                            return
                        end if
                    end if
                    if (rank_two_symmetric) then
                        call discover_rank_two_inverse_symmetric(fun, x, fvec, omega, nfev, &
                            use_central, jac_step_mode, rdiff_use, &
                            inner_m_diag_rank_two_symmetric_base, rank_two_right, &
                            rank_two_base, rank_two_correction, rank_two_enabled, &
                            inner_ok, jacobian_product, &
                            forward_calls=diag_refresh_forward)
                        if (.not. inner_ok) then
                            message = 'symmetric rank-two Jacobian evaluation failed'
                            return
                        end if
                    end if
                    if (.not. rank_two_enabled) then
                        if (rank_two_stream) then
                            call discover_rank_two_inverse_stream(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, diag_floor_use, &
                                diag_preconditioner, raw_jac_diag, rank_two_right, &
                                rank_two_base, rank_two_correction, rank_two_enabled, &
                                inner_ok, jacobian_product, &
                                forward_calls=diag_refresh_forward)
                            if (.not. inner_ok) then
                                message = 'streamed rank-two Jacobian evaluation failed'
                                return
                            end if
                        else if (rank_two_symmetric .or. rank_two_nonsymmetric) then
                            ! A failed low-product fit uses the original exact
                            ! diagonal rather than trusting the supplied D.
                            call estimate_inverse_diagonal(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, diag_floor_use, &
                                .false., diag_preconditioner, inner_ok, jacobian_product, &
                                0, 'walsh', 1, 0, 0.25_dp, 0, &
                                forward_calls=diag_refresh_forward)
                            ! Count exact fallback probes in the same measured cost
                            ! accumulator as the rejected low-product fit.
                        else if (rank_one_discover) then
                            ! Reuse the exact coordinate products both for diagonal
                            ! scaling and the full off-diagonal structure fit.
                            call estimate_inverse_diagonal(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, diag_floor_use, &
                                .false., diag_preconditioner, inner_ok, jacobian_product, &
                                diag_probes_use, trim(diag_probe_scheme_use), &
                                diag_probe_seed_use, diag_verify_use, diag_verify_rtol_use, &
                                diag_probe_max_use, raw_diagonal=raw_jac_diag, &
                                jacobian_matrix=rank_one_jac, &
                                forward_calls=diag_refresh_forward)
                        else if (rank_two_discover .and. .not. rank_two_symmetric .and. &
                                .not. rank_two_nonsymmetric) then
                            call estimate_inverse_diagonal(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, diag_floor_use, &
                                .false., diag_preconditioner, inner_ok, jacobian_product, &
                                diag_probes_use, trim(diag_probe_scheme_use), &
                                diag_probe_seed_use, diag_verify_use, diag_verify_rtol_use, &
                                diag_probe_max_use, raw_diagonal=raw_jac_diag, &
                                jacobian_matrix=rank_two_jac, &
                                forward_calls=diag_refresh_forward)
                        else if (rank_one_use .or. rank_two_use) then
                            call estimate_inverse_diagonal(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, diag_floor_use, &
                                trim(inner) == 'minres' .and. .not. present(inner_solver), &
                                diag_preconditioner, inner_ok, jacobian_product, diag_probes_use, &
                                trim(diag_probe_scheme_use), diag_probe_seed_use, diag_verify_use, &
                                diag_verify_rtol_use, diag_probe_max_use, raw_diagonal=raw_jac_diag, &
                                forward_calls=diag_refresh_forward)
                        else
                            call estimate_inverse_diagonal(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, diag_floor_use, &
                                trim(inner) == 'minres' .and. .not. present(inner_solver), &
                                diag_preconditioner, inner_ok, jacobian_product, diag_probes_use, &
                                trim(diag_probe_scheme_use), diag_probe_seed_use, diag_verify_use, &
                                diag_verify_rtol_use, diag_probe_max_use, &
                                forward_calls=diag_refresh_forward)
                        end if
                        if (.not. inner_ok) then
                            message = 'automatic diagonal preconditioner estimation failed'
                            return
                        end if
                        if (rank_one_discover) then
                            call discover_rank_one_inverse(rank_one_jac, raw_jac_diag, &
                                rank_right, rank_base, rank_correction, rank_enabled)
                        else if (rank_one_use) then
                            call estimate_row_coupled_inverse(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, raw_jac_diag, &
                                rank_base, rank_correction, rank_enabled, inner_ok, &
                                jacobian_product, rank_right, forward_calls=diag_refresh_forward)
                            if (.not. inner_ok) then
                                message = 'rank-one preconditioner verification failed'
                                return
                            end if
                        else if (rank_two_discover .and. .not. rank_two_stream .and. &
                                .not. rank_two_symmetric .and. .not. rank_two_nonsymmetric) then
                            call discover_rank_two_inverse(rank_two_jac, raw_jac_diag, &
                                rank_two_right, rank_two_base, rank_two_correction, &
                                rank_two_enabled)
                        else if (rank_two_use .and. .not. rank_two_stream .and. &
                                .not. rank_two_symmetric .and. .not. rank_two_nonsymmetric) then
                            call estimate_rank_two_inverse(fun, x, fvec, omega, nfev, &
                                use_central, jac_step_mode, rdiff_use, raw_jac_diag, &
                                rank_two_right, rank_two_base, rank_two_correction, &
                                rank_two_enabled, inner_ok, jacobian_product, &
                                forward_calls=diag_refresh_forward)
                            if (.not. inner_ok) then
                                message = 'rank-two preconditioner verification failed'
                                return
                            end if
                        end if
                    end if
                    if (diag_timed_cost_use) then
                        call stop_diag_timer(diag_setup_seconds, diag_setup_clock_valid)
                    end if
                    if (diag_measure_cost_use) then
                        ! Saturate the cost model rather than overflow for an
                        ! extreme, but finite user-supplied transpose weight.
                        diag_measured_rebuild = min(huge(1.0_dp)/4.0_dp, &
                            real(diag_refresh_forward, dp) + &
                            real(diag_refresh_transpose, dp) * &
                            min(diag_transpose_cost_use, huge(1.0_dp)/ &
                                (4.0_dp * real(max(1,diag_refresh_transpose),dp))))
                    end if
                end if
            end if
            ! Count only operator work done by the inner linear solver, not
            ! the Jacobian products used to build the preconditioner.
            if (diag_timed_cost_use) call start_diag_timer()
            diag_inner_actions = 0
            inner_tol = min(eta, eta * fx_norm)
            inner_rtol_use = min(1.0_dp, max(0.0_dp, inner_tol))
            if (present(inner_solver)) then
                call inner_solver(inner_matvec, inner_precond, fvec, &
                    inner_rtol_use, inner_maxiter_use, &
                    solve_vec, inner_ok)
            else if (trim(inner) == 'bicgstab' .or. trim(inner) == 'cgs' .or. &
                    trim(inner) == 'minres' .or. trim(inner) == 'tfqmr') then
                call krylov_iterative_solve(trim(inner), inner_matvec, &
                    inner_precond, fvec, inner_rtol_use, &
                    inner_maxiter_use, solve_vec, inner_ok, inner_verify_residual, &
                    inner_tfqmr_right_preconditioned)
            else if (trim(inner) == 'lgmres') then
                if (present(inner_m_diag) .or. diag_auto_use) then
                    call lgmres_jacobian_solve(fun, x, fvec, fvec, omega, &
                        inner_maxiter_use, outer_vec, outer_count, outer_k_use, &
                        inner_tol, solve_vec, nfev, inner_ok, inner_precond, use_central, &
                        jac_step_mode, rdiff_use, jacobian_product, orth_mode, outer_prune_tol_use, &
                        outer_image_prune_tol_use, outer_image_prune_auto_use, &
                        outer_image_reuse_use, outer_candidate_image_reuse_use, &
                        outer_fifo_prune_fast_use, outer_image_basis_reuse_use, &
                        candidate_cache_mib_use, gmres_work, image_arrays, inner_restart_use, &
                        inner_restart_adaptive_use, inner_restart_max_use)
                else
                    call lgmres_jacobian_solve(fun, x, fvec, fvec, omega, &
                        inner_maxiter_use, outer_vec, outer_count, outer_k_use, &
                        inner_tol, solve_vec, nfev, inner_ok, inner_m, use_central, &
                        jac_step_mode, rdiff_use, jacobian_product, orth_mode, outer_prune_tol_use, &
                        outer_image_prune_tol_use, outer_image_prune_auto_use, &
                        outer_image_reuse_use, outer_candidate_image_reuse_use, &
                        outer_fifo_prune_fast_use, outer_image_basis_reuse_use, &
                        candidate_cache_mib_use, gmres_work, image_arrays, inner_restart_use, &
                        inner_restart_adaptive_use, inner_restart_max_use)
                end if
            else
                if (present(inner_m_diag) .or. diag_auto_use) then
                    call gmres_jacobian_solve(fun, x, fvec, fvec, omega, &
                        inner_maxiter_use, inner_tol, solve_vec, nfev, inner_ok, inner_precond, &
                        use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, gmres_work, &
                        inner_restart_use)
                else
                    call gmres_jacobian_solve(fun, x, fvec, fvec, omega, &
                        inner_maxiter_use, inner_tol, solve_vec, nfev, inner_ok, inner_m, &
                        use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, gmres_work, &
                        inner_restart_use)
                end if
            end if
            if (diag_timed_cost_use) then
                call stop_diag_timer(diag_solve_seconds, diag_solve_clock_valid)
                diag_solve_clock_valid = diag_solve_clock_valid .and. diag_solve_seconds > 0.0_dp
                if (diag_refresh_now) then
                    ! Keep a positive, measurable reference solve. A clock
                    ! resolution below the minimum disables early timing updates.
                    diag_solve_baseline = 0.0_dp
                    if (diag_solve_clock_valid .and. diag_setup_clock_valid .and. &
                            diag_solve_seconds >= diag_timed_min_seconds_use) then
                        diag_solve_baseline = diag_solve_seconds
                    end if
                end if
            end if
            if (diag_work_use .and. diag_refresh_now) then
                ! No preconditioner calls means the solver provided no
                ! usable linear-work signal; disable work-triggered refresh
                ! until a later refresh supplies a positive baseline.
                diag_work_baseline = diag_inner_actions
            end if
            ! A failed user solver need not define its intent(out) solution.
            ! Do not inspect that vector until the callback reports success.
            if (.not. inner_ok) then
                message = 'Krylov Jacobian inversion failed'
                return
            end if
            if (.not. all(ieee_is_finite(solve_vec))) then
                message = 'Krylov Jacobian inversion failed'
                return
            end if

            ! SciPy accepts an incomplete inner iterate. Opt-in strict checking
            ! evaluates the actual J*p-f residual before the line search, even
            ! when a custom inner solver reported ok=.true. Normalization takes
            ! place before subtraction, so two very large, finite operands do
            ! not overflow under IEEE floating-point traps.
            if (check_inner) then
                call inner_matvec(solve_vec, linear_product, inner_ok)
                if (.not. inner_ok) then
                    message = 'Krylov inner residual verification failed'
                    return
                end if
                if (.not. inner_residual_acceptable(fvec, linear_product, inner_rtol_use)) then
                    message = 'Krylov inner residual exceeds tolerance'
                    return
                end if
            end if

            dx = -solve_vec
            if (all(dx == 0.0_dp)) then
                message = 'Jacobian inversion yielded zero vector'
                return
            end if

            if (use_armijo) then
                call armijo_step(fun, x, fvec, dx, xnew, fnew, nfev, step_status)
                if (step_status /= 0) then
                    message = 'nonfinite residual encountered'
                    return
                end if
            else if (use_wolfe) then
                call wolfe_step(fun, x, fvec, dx, xnew, fnew, nfev, step_status)
                if (step_status /= 0) then
                    message = 'nonfinite residual encountered'
                    return
                end if
            else
                xnew = x + dx
                call fun(xnew, fnew)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(fnew))) then
                    message = 'nonfinite residual encountered'
                    return
                end if
            end if

            dx_actual = xnew - x
            x = xnew
            fvec = fnew
            dx = dx_actual
            if (present(inner_m_update)) call inner_m_update(x, fvec)

            fx_norm_new = safe_two_norm(fvec)
            if (diag_adaptive_use) then
                ! Compare without division, even when norms are extremely
                ! small.  Only accepted (not trial) steps affect scheduling.
                diag_refresh_pending = fx_norm_new >= diag_stall_ratio_use * fx_norm
                if (diag_timed_cost_use .and. .not. diag_refresh_now .and. &
                        diag_solve_baseline > 0.0_dp .and. diag_solve_clock_valid) then
                    ! The setup cost must be earned back by a slower inner
                    ! solve. Division avoids overflow for huge growth factors.
                    ! A minimum sample duration guards against timer noise.
                    if (diag_solve_seconds > diag_setup_seconds .and. &
                            diag_solve_seconds >= diag_timed_min_seconds_use) then
                        if ((diag_solve_seconds - diag_setup_seconds) / &
                                diag_solve_baseline >= diag_work_growth_use) then
                            diag_refresh_pending = .true.
                        end if
                    end if
                end if
                if (diag_work_use .and. .not. diag_timed_cost_use .and. &
                        .not. diag_refresh_now .and. diag_work_baseline > 0) then
                    ! Compare a quotient rather than multiplying by the user
                    ! growth factor, which may be close to HUGE(real). The
                    ! cost term represents extra inner applications that must
                    ! be saved to justify refreshing early.
                    if (diag_measure_cost_use) then
                        if (real(diag_inner_actions,dp) > diag_measured_rebuild) then
                            if ((real(diag_inner_actions,dp) - diag_measured_rebuild) / &
                                    real(diag_work_baseline,dp) >= diag_work_growth_use) then
                                diag_refresh_pending = .true.
                            end if
                        end if
                    else if (diag_inner_actions > diag_rebuild_cost_use) then
                        if (real(diag_inner_actions - diag_rebuild_cost_use, dp) / &
                                real(diag_work_baseline, dp) >= diag_work_growth_use) then
                            diag_refresh_pending = .true.
                        end if
                    end if
                end if
            end if
            if (same_real(fx_norm, 0.0_dp)) then
                eta_a = 0.0_dp
            else
                ! Squaring the norms first can overflow or underflow even when
                ! their ratio is representable (e.g. 1e180 or 1e-180).
                ! eta_a is capped by eta_max below, so do not form a ratio
                ! greater than one only to square it.
                if (fx_norm_new >= fx_norm) then
                    ! Inverting the ratio avoids overflow when a residual
                    ! jumps by hundreds of orders of magnitude. Preserve
                    ! gamma*(new/old)**2 exactly when it is below eta_max.
                    norm_ratio = fx_norm / fx_norm_new
                    if (norm_ratio < sqrt(gamma / eta_max)) then
                        eta_a = eta_max
                    else
                        eta_a = gamma / (norm_ratio * norm_ratio)
                    end if
                else
                    norm_ratio = fx_norm_new / fx_norm
                    eta_a = gamma * norm_ratio**2
                end if
            end if
            if (gamma * eta * eta < eta_threshold) then
                eta = min(eta_max, eta_a)
            else
                eta = min(eta_max, max(eta_a, gamma * eta * eta))
            end if
            fx_norm = fx_norm_new
        end do

        nit = iter_count
        if (status == 1) then
            success = .true.
            message = 'A solution was found at the specified tolerance.'
        else if (status == -1) then
            status = 0
            message = 'invalid tolerance norm value'
        else
            status = 2
            success = .false.
            message = 'The maximum number of iterations allowed has been reached.'
        end if
    contains
        ! A single timer implementation for both preconditioner setup and
        ! the following inner solve. A non-monotone or wrapped system_clock
        ! sample is invalid; fixed-interval refresh still remains active.
        subroutine start_diag_timer()
            if (diag_wall_clock_use) then
                call system_clock(count=diag_wall_start, count_rate=diag_wall_rate)
            else
                call cpu_time(diag_clock_start)
            end if
        end subroutine start_diag_timer

        subroutine stop_diag_timer(seconds, valid)
            real(dp), intent(out) :: seconds !! elapsed seconds, zero if invalid
            logical, intent(out) :: valid !! true if timing sample is usable

            seconds = 0.0_dp
            if (diag_wall_clock_use) then
                call system_clock(count=diag_wall_end)
                valid = diag_wall_rate > 0_int64 .and. diag_wall_start >= 0_int64 .and. &
                    diag_wall_end >= diag_wall_start
                if (valid) then
                    seconds = real(diag_wall_end - diag_wall_start, dp) / &
                        real(diag_wall_rate, dp)
                end if
            else
                call cpu_time(diag_clock_end)
                valid = ieee_is_finite(diag_clock_start) .and. &
                    ieee_is_finite(diag_clock_end) .and. &
                    diag_clock_start >= 0.0_dp .and. diag_clock_end >= diag_clock_start
                if (valid) seconds = diag_clock_end - diag_clock_start
            end if
            if (valid) valid = ieee_is_finite(seconds) .and. seconds >= 0.0_dp
            if (.not. valid) seconds = 0.0_dp
        end subroutine stop_diag_timer

        subroutine inner_matvec(vec, product, ok)
            real(dp), intent(in) :: vec(:) !! custom solver's vector
            real(dp), intent(out) :: product(:) !! numerical Jacobian product
            logical, intent(out) :: ok !! true if product is finite

            call jacobian_vector(fun, x, fvec, vec, omega, product, nfev, ok, &
                use_central, jac_step_mode, rdiff_use, jacobian_product)
        end subroutine inner_matvec

        subroutine inner_precond(vec, product, ok)
            real(dp), intent(in) :: vec(:) !! vector to precondition
            real(dp), intent(out) :: product(:) !! inverse-Jacobian approximation
            logical, intent(out) :: ok !! true if product is finite

            if (diag_work_use .and. diag_inner_actions < huge(0)) then
                diag_inner_actions = diag_inner_actions + 1
            end if
            if (present(inner_m_diag) .or. diag_auto_use) then
                call diagonal_precond(vec, product, ok)
            else if (present(inner_m)) then
                call apply_nonlin_linear_operator(inner_m, vec, product, ok)
                if (ok) ok = all(ieee_is_finite(product))
            else
                product = vec
                ok = .true.
            end if
        end subroutine inner_precond

        subroutine diagonal_precond(vec, product, ok)
            real(dp), intent(in) :: vec(:) !! vector to precondition
            real(dp), intent(out) :: product(:) !! diagonal inverse-Jacobian action
            logical, intent(out) :: ok !! false for nonfinite or unrepresentable products
            integer :: k
            real(dp) :: factor

            product = 0.0_dp
            ok = .false.
            if (.not. all(ieee_is_finite(vec))) return
            do k = 1, size(vec)
                ! Avoid raising IEEE overflow before a subsequent finiteness
                ! check.  A factor with magnitude <= 1 cannot overflow.
                if (rank_two_enabled) then
                    factor = rank_two_base(k)
                else if (rank_enabled) then
                    factor = rank_base(k)
                else if (present(inner_m_diag)) then
                    factor = inner_m_diag(k)
                else
                    factor = diag_preconditioner(k)
                end if
                if (abs(factor) > 1.0_dp) then
                    if (abs(vec(k)) > huge(1.0_dp) / abs(factor)) return
                end if
                product(k) = factor * vec(k)
            end do
            if (rank_enabled) then
                ! Woodbury inverse: B*z = D^-1*z - corr*(right^T D^-1*z).
                ! Guard both accumulation and elementwise corrections before
                ! they can overflow when floating point traps are enabled.
                block
                    real(dp) :: scale, normalized_sum, projection, term
                    scale = maxval(abs(product))
                    if (scale > 0.0_dp) then
                        normalized_sum = sum(rank_right * (product / scale))
                        if (abs(normalized_sum) > 1.0_dp) then
                            if (scale > huge(1.0_dp) / abs(normalized_sum)) return
                        end if
                        projection = scale * normalized_sum
                        do k = 1, size(vec)
                            factor = rank_correction(k)
                            if (abs(factor) > 1.0_dp) then
                                if (abs(projection) > huge(1.0_dp) / abs(factor)) return
                            end if
                            term = factor * projection
                            if (sign(1.0_dp, product(k)) /= sign(1.0_dp, term)) then
                                if (abs(product(k)) > huge(1.0_dp) - abs(term)) return
                            end if
                            product(k) = product(k) - term
                        end do
                    end if
                end block
            end if
            if (rank_two_enabled) then
                block
                    real(dp) :: scale, normalized_sum, projection(2), term(2), total
                    integer :: j
                    scale = maxval(abs(product))
                    if (scale > 0.0_dp) then
                        do j = 1, 2
                            normalized_sum = sum(rank_two_right(:,j) * (product / scale))
                            if (abs(normalized_sum) > 1.0_dp) then
                                if (scale > huge(1.0_dp)/abs(normalized_sum)) return
                            end if
                            projection(j) = scale * normalized_sum
                        end do
                        do k = 1, size(vec)
                            do j = 1, 2
                                factor = rank_two_correction(k,j)
                                if (abs(factor) > 1.0_dp) then
                                    if (abs(projection(j)) > huge(1.0_dp)/abs(factor)) return
                                end if
                                term(j) = factor * projection(j)
                            end do
                            if (sign(1.0_dp,term(1)) == sign(1.0_dp,term(2))) then
                                if (abs(term(1)) > huge(1.0_dp)-abs(term(2))) return
                            end if
                            total = term(1) + term(2)
                            if (sign(1.0_dp,product(k)) /= sign(1.0_dp,total)) then
                                if (abs(product(k)) > huge(1.0_dp)-abs(total)) return
                            end if
                            product(k) = product(k) - total
                        end do
                    end if
                end block
            end if
            ok = all(ieee_is_finite(product))
        end subroutine diagonal_precond

    end subroutine solve_krylov_root

    subroutine discover_rank_one_inverse(jac, diagonal, right, base, correction, enabled)
        ! Infer a full unknown rank-one right factor from the off-diagonal
        ! columns of an explicitly probed Jacobian. For rank-one off-diagonal
        ! structure J(i,j)=u(i)*v(j) (i/=j), a strong off-diagonal pivot
        ! determines all but one u and one v entry. Recover those entries
        ! from independent rows and columns, then verify EVERY off-diagonal
        ! element before using Sherman-Morrison. Rejection retains the
        ! ordinary, already estimated inverse diagonal.
        real(dp), intent(in) :: jac(:,:) !! complete Jacobian at current iterate
        real(dp), intent(in) :: diagonal(:) !! exact Jacobian diagonal
        real(dp), intent(out) :: right(:) !! fitted normalized right vector
        real(dp), intent(out) :: base(:) !! inverse fitted uncoupled diagonal
        real(dp), intent(out) :: correction(:) !! Sherman-Morrison factor
        logical, intent(out) :: enabled !! true only for a verified safe fit

        real(dp) :: u(size(diagonal)), v(size(diagonal)), a(size(diagonal))
        real(dp) :: work(size(diagonal)), pivot, maxoff, candidate, scale
        real(dp) :: threshold, denom, bound, value, vmax
        integer :: n, i, j, p, q, r, s

        enabled = .false.
        right = 1.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        n = size(diagonal)
        if (n < 4) return
        if (size(jac,1) /= n .or. size(jac,2) /= n) return
        if (.not. all(ieee_is_finite(jac))) return
        if (maxval(abs(jac)) > sqrt(huge(1.0_dp)/real(n,dp))) return

        p = 0
        q = 0
        maxoff = 0.0_dp
        do j = 1, n
            do i = 1, n
                if (i == j) cycle
                candidate = abs(jac(i,j))
                if (candidate > maxoff) then
                    maxoff = candidate
                    p = i
                    q = j
                end if
            end do
        end do
        if (maxoff <= sqrt(tiny(1.0_dp))) return
        pivot = jac(p,q)
        u = 0.0_dp
        v = 0.0_dp
        u(p) = 1.0_dp
        do i = 1, n
            if (i /= q) u(i) = jac(i,q)/pivot
            if (i /= p) v(i) = jac(p,i)
        end do

        ! Missing entries u(q) and v(p) cannot be read from the diagonal.
        ! Independent off-diagonal observations are needed to identify them.
        r = 0
        s = 0
        value = 0.0_dp
        do j = 1, n
            if (j == p .or. j == q) cycle
            if (abs(v(j)) > value) then
                r = j
                value = abs(v(j))
            end if
        end do
        threshold = 100.0_dp*epsilon(1.0_dp)*maxoff
        if (r == 0 .or. value <= threshold) return
        value = 0.0_dp
        do i = 1, n
            if (i == p .or. i == q) cycle
            if (abs(u(i)) > value) then
                s = i
                value = abs(u(i))
            end if
        end do
        if (s == 0 .or. value <= 100.0_dp*epsilon(1.0_dp)) return
        u(q) = jac(q,r)/v(r)
        v(p) = jac(s,p)/u(s)
        if (.not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) return

        ! A full off-diagonal check is possible because the coordinate Jv
        ! products have already been computed for exact diagonal scaling.
        do j = 1, n
            do i = 1, n
                if (i == j) cycle
                candidate = u(i)*v(j)
                scale = max(abs(candidate),abs(jac(i,j)),sqrt(tiny(1.0_dp)))
                if (abs(candidate/scale-jac(i,j)/scale) > 1.0e-5_dp) return
            end do
        end do
        vmax = maxval(abs(v))
        if (vmax == 0.0_dp) return
        right = v/vmax
        u = u*vmax
        a = diagonal-u*right
        if (.not. all(ieee_is_finite(a))) return
        if (any(abs(a) < sqrt(tiny(1.0_dp)))) return
        base = 1.0_dp/a

        ! Match the rank-one solver's guarded Woodbury evaluation. Preserve
        ! headroom for the subsequent weighted sum and correction products.
        do i = 1, n
            if (abs(base(i)) > 1.0_dp) then
                if (abs(u(i)) > (huge(1.0_dp)/(4.0_dp*real(n,dp)))/ &
                        abs(base(i))) return
            else
                if (abs(u(i)) > huge(1.0_dp)/(4.0_dp*real(n,dp))) return
            end if
        end do
        work = base*u
        bound = sum(abs(work*right))
        if (.not. ieee_is_finite(bound)) return
        denom = 1.0_dp+sum(work*right)
        threshold = 100.0_dp*epsilon(1.0_dp)*(1.0_dp+bound)
        if (.not. ieee_is_finite(denom)) return
        if (abs(denom) <= threshold) return
        correction = work/denom
        if (.not. all(ieee_is_finite(correction))) return
        enabled = .true.
    end subroutine discover_rank_one_inverse

    subroutine discover_rank_two_inverse(jac, diagonal, right, base, correction, enabled)
        ! Recover a diagonal-plus-rank-two Jacobian from full coordinate Jv
        ! products without assuming either coupling factor is known. The
        ! diagonal entries are masked while fitting the low-rank component.
        ! An off-diagonal 2-by-2 anchor and independent rows/columns complete
        ! both factors; check every off-diagonal entry before using Woodbury.
        ! Degenerate or numerically unsafe completions retain diagonal scaling.
        real(dp), intent(in) :: jac(:,:) !! full measured square Jacobian
        real(dp), intent(in) :: diagonal(:) !! measured unmodified diagonal
        real(dp), intent(out) :: right(:,:) !! normalized learned right factors (n,2)
        real(dp), intent(out) :: base(:) !! reciprocal uncoupled diagonal
        real(dp), intent(out) :: correction(:,:) !! Woodbury correction (n,2)
        logical, intent(out) :: enabled !! true for fully verified safe rank-two inverse

        integer :: n, i, j, k, l, r(2), c(2), pivot_row, pivot_col
        integer :: first, second
        real(dp) :: left(size(diagonal),2), weights(size(diagonal),2)
        real(dp) :: u(size(diagonal),2), a(size(diagonal)), work(size(diagonal),2)
        real(dp) :: pivot(2,2), small(2,2), normalized(2,2), cof(2,2)
        real(dp) :: maxoff, candidate, best, determinant, scale, val1, val2
        real(dp) :: vmax, max_allowed, threshold, predicted, measured

        n = size(diagonal)
        enabled = .false.
        right = 0.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        ! At least two independent off-diagonal observations are required
        ! for every unknown row and column, including the four anchors.
        if (n < 6) return
        if (size(jac,1) /= n .or. size(jac,2) /= n) return
        if (.not. all(ieee_is_finite(jac))) return
        max_allowed = sqrt(huge(1.0_dp))/64.0_dp
        if (maxval(abs(jac)) > max_allowed) return

        maxoff = 0.0_dp
        r = 0
        c = 0
        do j = 1, n
            do i = 1, n
                if (i == j) cycle
                if (abs(jac(i,j)) > maxoff) then
                    maxoff = abs(jac(i,j))
                    r(1) = i
                    c(1) = j
                end if
            end do
        end do
        if (maxoff <= sqrt(tiny(1.0_dp))) return
        ! Choose the most independent second row/column, with no anchor
        ! on the matrix diagonal. The divisions bound pivot arithmetic.
        best = 0.0_dp
        pivot(1,1) = jac(r(1),c(1))/maxoff
        do j = 1, n
            if (j == r(1) .or. j == c(1)) cycle
            do i = 1, n
                if (i == r(1) .or. i == c(1) .or. i == j) cycle
                candidate = abs(pivot(1,1)*(jac(i,j)/maxoff) - &
                    (jac(r(1),j)/maxoff)*(jac(i,c(1))/maxoff))
                if (candidate > best) then
                    best = candidate
                    r(2) = i
                    c(2) = j
                end if
            end do
        end do
        if (best <= 1.0e-8_dp) return
        pivot(1,2) = jac(r(1),c(2))/maxoff
        pivot(2,1) = jac(r(2),c(1))/maxoff
        pivot(2,2) = jac(r(2),c(2))/maxoff
        determinant = pivot(1,1)*pivot(2,2)-pivot(1,2)*pivot(2,1)
        if (abs(determinant) <= 1.0e-8_dp) return

        ! Fix the factor gauge with V(j,:) = J(r(:),j)/maxoff for all
        ! columns other than the two anchor rows. Recover U from two
        ! anchor columns wherever those observations are off-diagonal.
        left = 0.0_dp
        weights = 0.0_dp
        do j = 1, n
            if (j == r(1) .or. j == r(2)) cycle
            weights(j,1) = jac(r(1),j)/maxoff
            weights(j,2) = jac(r(2),j)/maxoff
        end do
        do i = 1, n
            if (i == c(1) .or. i == c(2)) cycle
            val1 = jac(i,c(1))/maxoff
            val2 = jac(i,c(2))/maxoff
            left(i,1) = (val1*pivot(2,2)-val2*pivot(2,1))/determinant
            left(i,2) = (val2*pivot(1,1)-val1*pivot(1,2))/determinant
        end do
        ! For either column-anchor row, its diagonal masked one fitting
        ! entry. Find a different pair of known independent right rows.
        do k = 1, 2
            i = c(k)
            best = 0.0_dp
            first = 0
            second = 0
            do j = 1, n-1
                if (j == i .or. any(j == r)) cycle
                do l = j+1, n
                    if (l == i .or. any(l == r)) cycle
                    candidate = abs(weights(j,1)*weights(l,2)- &
                        weights(j,2)*weights(l,1))
                    if (candidate > best) then
                        best = candidate
                        first = j
                        second = l
                    end if
                end do
            end do
            if (best <= 1.0e-8_dp) return
            determinant = weights(first,1)*weights(second,2) - &
                weights(first,2)*weights(second,1)
            val1 = jac(i,first)/maxoff
            val2 = jac(i,second)/maxoff
            left(i,1) = (val1*weights(second,2)-val2*weights(first,2))/determinant
            left(i,2) = (val2*weights(first,1)-val1*weights(second,1))/determinant
        end do
        if (.not. all(ieee_is_finite(left))) return
        if (maxval(abs(left)) > 1.0e8_dp) return

        ! Complete the right factors in the anchor-row columns using
        ! off-diagonal measurements from two other independent rows.
        do k = 1, 2
            j = r(k)
            best = 0.0_dp
            first = 0
            second = 0
            do i = 1, n-1
                if (i == j) cycle
                do l = i+1, n
                    if (l == j) cycle
                    candidate = abs(left(i,1)*left(l,2)-left(i,2)*left(l,1))
                    if (candidate > best) then
                        best = candidate
                        first = i
                        second = l
                    end if
                end do
            end do
            if (best <= 1.0e-8_dp) return
            determinant = left(first,1)*left(second,2) - &
                left(first,2)*left(second,1)
            val1 = jac(first,j)/maxoff
            val2 = jac(second,j)/maxoff
            weights(j,1) = (val1*left(second,2)-val2*left(first,2))/determinant
            weights(j,2) = (val2*left(first,1)-val1*left(second,1))/determinant
        end do
        if (.not. all(ieee_is_finite(weights))) return
        if (maxval(abs(weights)) > 1.0e8_dp) return

        ! All off-diagonal entries must agree. Diagonal entries are never
        ! included in the fit and will be absorbed into D afterwards.
        do j = 1, n
            do i = 1, n
                if (i == j) cycle
                predicted = left(i,1)*weights(j,1)+left(i,2)*weights(j,2)
                measured = jac(i,j)/maxoff
                scale = max(abs(predicted),abs(measured),sqrt(tiny(1.0_dp)))
                if (abs(predicted/scale-measured/scale) > 1.0e-5_dp) return
            end do
        end do

        u = left
        do k = 1, 2
            vmax = maxval(abs(weights(:,k)))
            if (vmax <= 1.0e-8_dp) return
            right(:,k) = weights(:,k)/vmax
            scale = maxoff*vmax
            if (maxval(abs(left(:,k))) > max_allowed/scale) return
            u(:,k) = left(:,k)*scale
        end do
        a = diagonal - u(:,1)*right(:,1) - u(:,2)*right(:,2)
        if (.not. all(ieee_is_finite(a))) return
        if (any(abs(a) <= sqrt(tiny(1.0_dp)))) return
        base = 1.0_dp/a
        do k = 1, 2
            do i = 1, n
                if (abs(base(i)) > 1.0_dp) then
                    if (abs(u(i,k)) > max_allowed/abs(base(i))) return
                else
                    if (abs(u(i,k)) > max_allowed) return
                end if
            end do
            work(:,k) = base*u(:,k)
        end do
        if (maxval(abs(work)) > max_allowed) return

        small = matmul(transpose(right),work)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        scale = max(1.0_dp,maxval(abs(small)))
        normalized = small/scale
        determinant = normalized(1,1)*normalized(2,2) - &
            normalized(1,2)*normalized(2,1)
        threshold = 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp,sum(abs(normalized(:,1)))*sum(abs(normalized(:,2))))
        if (abs(determinant) <= threshold) return
        cof(1,1) = normalized(2,2)
        cof(1,2) = -normalized(1,2)
        cof(2,1) = -normalized(2,1)
        cof(2,2) = normalized(1,1)
        correction = matmul(work,cof)/(scale*determinant)
        if (.not. all(ieee_is_finite(correction))) return
        enabled = .true.
    end subroutine discover_rank_two_inverse

    subroutine discover_rank_two_inverse_nonsymmetric(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, d, right, base, correction, &
            enabled, ok, jacobian_transpose_product, jacobian_product, &
            verify_transpose, adaptive_transpose_cost, forward_calls, transpose_calls)
        ! General rank-two cross approximation with a known diagonal D.
        ! L = J-D, C=L*Q, R=L^T*P, B=P^T*C give L=C*B^-1*R^T.
        ! Hold-out Jv products guard the heuristic before Woodbury inversion.
        procedure(nonlin_function) :: fun !! nonlinear residual callback
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(in) :: fvec(:) !! residual at current point
        real(dp), intent(in) :: omega !! finite difference scale
        integer, intent(inout) :: nfev !! residual evaluation count (not analytic products)
        logical, intent(in) :: use_central !! centered numerical Jacobian difference
        integer, intent(in) :: jac_step_mode !! global, adaptive, or componentwise step mode
        real(dp), intent(in) :: rdiff_use !! requested finite difference step
        real(dp), intent(in) :: d(:) !! known uncoupled diagonal of D
        real(dp), intent(out) :: right(:,:) !! fitted Woodbury right factor
        real(dp), intent(out) :: base(:) !! reciprocal uncoupled diagonal
        real(dp), intent(out) :: correction(:,:) !! Woodbury correction columns
        logical, intent(out) :: enabled !! true only for a verified, safe factorization
        logical, intent(out) :: ok !! false when a Jacobian action fails
        procedure(nonlin_jacobian_product) :: jacobian_transpose_product !! J(x)^T*v callback
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! analytic Jv callback
        logical, intent(in) :: verify_transpose !! verify two independent adjoint actions
        real(dp), intent(in) :: adaptive_transpose_cost !! J^T*v cost / Jv cost, 0 disables
        integer, intent(out), optional :: forward_calls !! attempted forward actions during setup
        integer, intent(out), optional :: transpose_calls !! attempted adjoint actions during setup

        integer :: n, i, j, k
        real(dp) :: q(size(x),2), p(size(x),2), c(size(x),2), r(size(x),2)
        real(dp) :: left(size(x),2), weights(size(x),2), basis(size(x),2)
        real(dp) :: transpose_basis(size(x),2)
        real(dp) :: trial(size(x)), observed(size(x)), predicted(size(x))
        real(dp) :: diagonal(size(x)), b(2,2), rcheck(2,2), invb(2,2)
        real(dp) :: bc, rc, determinant, scale, norm_value, error_scale, theta
        real(dp) :: forward_discrepancy, pivot_quality
        real(dp), parameter :: pivot_tol = 1.0e-8_dp
        real(dp), parameter :: verify_tol = 2.0e-6_dp
        real(dp), parameter :: max_allowed = sqrt(huge(1.0_dp))/64.0_dp
        logical :: valid
        integer :: transpose_checks

        enabled = .false.
        ok = .true.
        right = 0.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        n = size(x)
        if (present(forward_calls)) forward_calls = 0
        if (present(transpose_calls)) transpose_calls = 0
        if (n < 6) return
        if (maxval(abs(d)) > max_allowed) return

        ! Use different reproducible left and right projection subspaces.
        do i = 1, n
            theta = real(i,dp)
            q(i,1) = sin(0.67_dp*theta) + cos(1.17_dp*theta)
            q(i,2) = sin(1.41_dp*theta) - cos(0.43_dp*theta)
            p(i,1) = sin(0.53_dp*theta) + cos(0.93_dp*theta)
            p(i,2) = cos(1.31_dp*theta) - sin(0.79_dp*theta)
        end do
        do j = 1, 2
            do k = 1, j-1
                q(:,j) = q(:,j) - dot_product(q(:,j),q(:,k))*q(:,k)
                p(:,j) = p(:,j) - dot_product(p(:,j),p(:,k))*p(:,k)
            end do
            norm_value = sqrt(sum(q(:,j)**2))
            if (norm_value <= pivot_tol) return
            q(:,j) = q(:,j)/norm_value
            norm_value = sqrt(sum(p(:,j)**2))
            if (norm_value <= pivot_tol) return
            p(:,j) = p(:,j)/norm_value

            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, q(:,j), omega, observed, &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
            c(:,j) = observed - d*q(:,j)
            if (present(transpose_calls)) transpose_calls = transpose_calls + 1
            call jacobian_transpose_product(x, p(:,j), observed, valid)
            if (valid) valid = all(ieee_is_finite(observed))
            if (.not. valid) then
                ok = .false.
                return
            end if
            r(:,j) = observed - d*p(:,j)
        end do
        if (.not. all(ieee_is_finite(c)) .or. .not. all(ieee_is_finite(r))) return
        bc = maxval(abs(c))
        rc = maxval(abs(r))
        if (bc <= sqrt(tiny(1.0_dp)) .or. rc <= sqrt(tiny(1.0_dp))) return
        if (bc > max_allowed .or. rc > max_allowed) return
        left = c/bc
        r = r/rc
        b = matmul(transpose(p),left)
        rcheck = matmul(transpose(q),r)
        if (.not. all(ieee_is_finite(b)) .or. .not. all(ieee_is_finite(rcheck))) return
        ! This is an independent bilinear consistency check on Jv and J^Tv.
        scale = max(1.0_dp,maxval(abs(b*bc)),maxval(abs(rcheck*rc)))
        if (maxval(abs(b*bc-transpose(rcheck)*rc))/scale > verify_tol) return
        scale = maxval(abs(b))
        if (scale <= pivot_tol) return
        b = b/scale
        determinant = b(1,1)*b(2,2)-b(1,2)*b(2,1)
        if (abs(determinant) <= pivot_tol) return
        ! The normalized 2-by-2 determinant is a cheap fit-conditioning indicator.
        pivot_quality = abs(determinant)/ &
            max(1.0_dp, sum(abs(b(:,1)))*sum(abs(b(:,2))))
        invb(1,1) = b(2,2)/(determinant*scale)
        invb(1,2) = -b(1,2)/(determinant*scale)
        invb(2,1) = -b(2,1)/(determinant*scale)
        invb(2,2) = b(1,1)/(determinant*scale)
        if (.not. all(ieee_is_finite(invb))) return
        weights = matmul(r,transpose(invb))
        if (.not. all(ieee_is_finite(weights))) return
        if (maxval(abs(weights)) > 1.0e8_dp) return

        forward_discrepancy = 0.0_dp
        ! The two independent verification directions have no overlap with
        ! the fitting subspace; a higher-rank component may be detected here.
        do j = 1, 2
            do i = 1, n
                theta = real(i,dp)
                if (j == 1) then
                    trial(i) = sin(0.31_dp*theta) + cos(1.73_dp*theta)
                else
                    trial(i) = cos(0.89_dp*theta) - sin(1.03_dp*theta)
                end if
            end do
            do k = 1, 2
                trial = trial - dot_product(trial,q(:,k))*q(:,k)
            end do
            do k = 1, j-1
                trial = trial - dot_product(trial,basis(:,k))*basis(:,k)
            end do
            norm_value = sqrt(sum(trial**2))
            if (norm_value <= pivot_tol) return
            basis(:,j) = trial/norm_value
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, basis(:,j), omega, observed, &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
            trial = observed - d*basis(:,j)
            predicted = matmul(left,matmul(transpose(weights),basis(:,j)))*rc
            if (.not. all(ieee_is_finite(predicted))) return
            error_scale = max(1.0_dp,maxval(abs(trial)),maxval(abs(predicted)))
            forward_discrepancy = max(forward_discrepancy, &
                maxval(abs(trial/error_scale-predicted/error_scale)))
            if (forward_discrepancy > verify_tol) return
        end do
        transpose_checks = 0
        if (verify_transpose) then
            transpose_checks = 2
        else if (adaptive_transpose_cost > 0.0_dp) then
            ! Limit nominal verification cost to one quarter of n exact
            ! coordinate products, while always checking at least one adjoint.
            transpose_checks = 1
            if (adaptive_transpose_cost <= 0.125_dp*real(n,dp)) transpose_checks = 2
            ! Nearly dependent fitting directions or borderline forward
            ! consistency demand both checks regardless of callback cost.
            if (pivot_quality <= 1.0e-4_dp .or. &
                    forward_discrepancy >= 0.1_dp*verify_tol) transpose_checks = 2
        end if
        if (transpose_checks > 0) then
            ! Optional held-out adjoint checks catch J^T callbacks that agree
            ! on the fitting vectors but disagree on other directions.
            do j = 1, transpose_checks
                do i = 1, n
                    theta = real(i,dp)
                    if (j == 1) then
                        trial(i) = cos(0.27_dp*theta) + sin(1.91_dp*theta)
                    else
                        trial(i) = sin(0.97_dp*theta) - cos(0.81_dp*theta)
                    end if
                end do
                do k = 1, 2
                    trial = trial - dot_product(trial,p(:,k))*p(:,k)
                end do
                do k = 1, j-1
                    trial = trial - dot_product(trial,transpose_basis(:,k))*transpose_basis(:,k)
                end do
                norm_value = sqrt(sum(trial**2))
                if (norm_value <= pivot_tol) return
                transpose_basis(:,j) = trial/norm_value
                if (present(transpose_calls)) transpose_calls = transpose_calls + 1
                call jacobian_transpose_product(x, transpose_basis(:,j), observed, valid)
                if (valid) valid = all(ieee_is_finite(observed))
                if (.not. valid) then
                    ok = .false.
                    return
                end if
                trial = observed - d*transpose_basis(:,j)
                predicted = matmul(weights, &
                    matmul(transpose(left),transpose_basis(:,j)))*rc
                if (.not. all(ieee_is_finite(predicted))) return
                error_scale = max(1.0_dp,maxval(abs(trial)),maxval(abs(predicted)))
                if (maxval(abs(trial/error_scale-predicted/error_scale)) > verify_tol) return
            end do
        end if
        diagonal = d + rc*(left(:,1)*weights(:,1) + left(:,2)*weights(:,2))
        if (.not. all(ieee_is_finite(diagonal))) return
        ! Fitted factors carry sampling error: reject near-singular inverse
        ! systems more aggressively than the exact-factor discovery paths.
        call inverse_rank_two_factors(diagonal, left, weights, rc, &
            right, base, correction, enabled, min_det=1.0e-8_dp)
    end subroutine discover_rank_two_inverse_nonsymmetric

    subroutine discover_rank_two_inverse_symmetric(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, d, right, base, correction, &
            enabled, ok, jacobian_product, forward_calls)
        ! Low-product Nystrom fit for a symmetric, rank-two remainder J-diag(d).
        ! This is a heuristic: two held-out directions cannot certify the
        ! structure of an arbitrary n-by-n matrix. Reject unsafe fits.
        procedure(nonlin_function) :: fun !! nonlinear residual callback
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(in) :: fvec(:) !! residual at current point
        real(dp), intent(in) :: omega !! global finite-difference scaling
        integer, intent(inout) :: nfev !! accumulated residual evaluation count
        logical, intent(in) :: use_central !! centered finite differences
        integer, intent(in) :: jac_step_mode !! finite-difference step mode
        real(dp), intent(in) :: rdiff_use !! relative finite-difference step
        real(dp), intent(in) :: d(:) !! known uncoupled diagonal, not diag(J)
        real(dp), intent(out) :: right(:,:) !! learned right rank-two factor
        real(dp), intent(out) :: base(:) !! Woodbury reciprocal diagonal
        real(dp), intent(out) :: correction(:,:) !! Woodbury correction columns
        logical, intent(out) :: enabled !! verified and numerically safe fit
        logical, intent(out) :: ok !! false only on a failed Jacobian evaluation
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! optional analytic Jv
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products

        integer :: n, i, j, k
        real(dp) :: fit(size(x),2), left(size(x),2), weights(size(x),2)
        real(dp) :: basis(size(x),2), observed(size(x)), predicted(size(x))
        real(dp) :: residual(size(x)), diagonal(size(x)), b(2,2), invb(2,2)
        real(dp) :: determinant, scale, bscale, norm_value, error_scale
        real(dp) :: value, theta
        logical :: valid
        real(dp), parameter :: pivot_tol = 1.0e-8_dp
        real(dp), parameter :: verify_tol = 2.0e-6_dp

        enabled = .false.
        ok = .true.
        right = 0.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        n = size(x)
        if (n < 6) return
        if (maxval(abs(d)) > sqrt(huge(1.0_dp))/64.0_dp) return

        ! Distinct bounded trigonometric directions make a reproducible fit.
        ! Modified Gram-Schmidt prevents redundant fitting directions.
        do i = 1, n
            theta = real(i,dp)
            fit(i,1) = sin(0.67_dp*theta) + cos(1.17_dp*theta)
            fit(i,2) = sin(1.41_dp*theta) - cos(0.43_dp*theta)
        end do
        do j = 1, 2
            do k = 1, j-1
                fit(:,j) = fit(:,j) - dot_product(fit(:,j),fit(:,k))*fit(:,k)
            end do
            norm_value = sqrt(sum(fit(:,j)**2))
            if (norm_value <= pivot_tol) return
            fit(:,j) = fit(:,j)/norm_value
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, fit(:,j), omega, observed, &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
            left(:,j) = observed - d*fit(:,j)
        end do
        scale = maxval(abs(left))
        if (.not. ieee_is_finite(scale)) return
        if (scale <= sqrt(tiny(1.0_dp))) return
        if (scale > sqrt(huge(1.0_dp))/64.0_dp) return
        left = left / scale
        b = matmul(transpose(fit),left)
        bscale = maxval(abs(b))
        if (bscale <= sqrt(tiny(1.0_dp))) return
        if (abs(b(1,2)-b(2,1)) > verify_tol*bscale) return
        b = b/bscale
        determinant = b(1,1)*b(2,2) - b(1,2)*b(2,1)
        if (abs(determinant) <= pivot_tol) return
        invb(1,1) = b(2,2)/(determinant*bscale)
        invb(1,2) = -b(1,2)/(determinant*bscale)
        invb(2,1) = -b(2,1)/(determinant*bscale)
        invb(2,2) = b(1,1)/(determinant*bscale)
        if (.not. all(ieee_is_finite(invb))) return
        weights = matmul(left,transpose(invb))
        if (.not. all(ieee_is_finite(weights))) return
        if (maxval(abs(weights)) > 1.0e8_dp) return

        ! Holdout probes are orthogonal to both fitting directions and to
        ! each other; they can detect a rank > 2 or nonsymmetric remainder.
        do j = 1, 2
            do i = 1, n
                theta = real(i,dp)
                if (j == 1) then
                    residual(i) = sin(0.31_dp*theta) + cos(1.73_dp*theta)
                else
                    residual(i) = cos(0.89_dp*theta) - sin(1.03_dp*theta)
                end if
            end do
            do k = 1, 2
                residual = residual - dot_product(residual,fit(:,k))*fit(:,k)
            end do
            do k = 1, j-1
                residual = residual - dot_product(residual,basis(:,k))*basis(:,k)
            end do
            norm_value = sqrt(sum(residual**2))
            if (norm_value <= pivot_tol) return
            basis(:,j) = residual/norm_value
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, basis(:,j), omega, observed, &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
            residual = observed - d*basis(:,j)
            predicted = matmul(left,matmul(transpose(weights),basis(:,j))) * scale
            if (.not. all(ieee_is_finite(predicted))) return
            error_scale = max(1.0_dp,maxval(abs(residual)),maxval(abs(predicted)))
            if (maxval(abs(residual/error_scale-predicted/error_scale)) > verify_tol) return
        end do

        diagonal = d + scale*(left(:,1)*weights(:,1) + left(:,2)*weights(:,2))
        if (.not. all(ieee_is_finite(diagonal))) return
        call inverse_rank_two_factors(diagonal, left, weights, scale, &
            right, base, correction, enabled)
    end subroutine discover_rank_two_inverse_symmetric

    subroutine discover_rank_two_inverse_stream(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, floor_fraction, invdiag, &
            diagonal, right, base, correction, enabled, ok, jacobian_product, forward_calls)
        ! Exact streamed discovery for J = D + U*V^T with two unknown factors.
        ! Use at most eight spread-out pilot columns to find independent
        ! off-diagonal rows/columns. Reuse all n coordinate Jv products for
        ! both the diagonal and a complete off-diagonal verification. Store
        ! O(n) reals rather than a dense n-by-n Jacobian. If no reliable
        ! anchor exists, preserve the exact diagonal fallback.
        procedure(nonlin_function) :: fun !! square residual callback
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: fvec(:) !! residual at current iterate
        real(dp), intent(in) :: omega !! directional numerical-Jacobian scale
        integer, intent(inout) :: nfev !! number of residual evaluations
        logical, intent(in) :: use_central !! use centered differences
        integer, intent(in) :: jac_step_mode !! finite-difference scaling mode
        real(dp), intent(in) :: rdiff_use !! coordinatewise relative step
        real(dp), intent(in) :: floor_fraction !! relative floor for diagonal scaling
        real(dp), intent(out) :: invdiag(:) !! exact inverse-diagonal fallback
        real(dp), intent(out) :: diagonal(:) !! measured Jacobian diagonal
        real(dp), intent(out) :: right(:,:) !! fitted right factors (n,2)
        real(dp), intent(out) :: base(:) !! fitted reciprocal uncoupled diagonal
        real(dp), intent(out) :: correction(:,:) !! fitted Woodbury correction
        logical, intent(out) :: enabled !! true when the model passes exact checks
        logical, intent(out) :: ok !! false only on failed Jacobian evaluation
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! analytic Jv callback
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products

        real(dp), allocatable :: basis(:), product(:), pilot(:,:), left(:,:)
        real(dp), allocatable :: weights(:,:), saved_rows(:,:), special(:,:)
        real(dp) :: best, determinant, pivot(2,2), pairmax, global_scale
        real(dp) :: value, predicted, scale, maxabs, cutoff
        integer :: n, p, j, k, i, a, b, r(2), c(2), best_a, best_b
        integer :: rfirst, rsecond, pilot_index(8), trial, allocation_status
        logical :: valid, candidate, verified
        real(dp), parameter :: pivot_tol = 1.0e-8_dp
        real(dp), parameter :: limit = sqrt(huge(1.0_dp))/64.0_dp

        n = size(x)
        ok = .false.
        enabled = .false.
        invdiag = 1.0_dp
        diagonal = 0.0_dp
        right = 0.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        p = min(8,n)
        allocate(basis(n), product(n), pilot(n,p), left(n,2), &
            weights(n,2), saved_rows(n,2), special(n,2), stat=allocation_status)
        if (allocation_status /= 0) return
        pilot = 0.0_dp
        left = 0.0_dp
        weights = 0.0_dp
        saved_rows = 0.0_dp
        special = 0.0_dp
        basis = 0.0_dp
        pilot_index = 0
        if (p > 1) then
            do k = 1, p
                pilot_index(k) = 1 + (k-1)*(n-1)/(p-1)
                j = pilot_index(k)
                basis(j) = 1.0_dp
                if (present(forward_calls)) forward_calls = forward_calls + 1
                call jacobian_vector(fun, x, fvec, basis, omega, pilot(:,k), &
                    nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                basis(j) = 0.0_dp
                if (.not. valid) return
            end do
        else if (p == 1) then
            pilot_index(1) = 1
            basis(1) = 1.0_dp
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, basis, omega, pilot(:,1), &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            basis(1) = 0.0_dp
            if (.not. valid) return
        end if

        ! Find a well-conditioned off-diagonal 2-by-2 minor without retaining
        ! any columns other than the pilot set. All indices are distinct.
        candidate = .false.
        best = pivot_tol
        r = 0
        c = 0
        global_scale = 0.0_dp
        best_a = 0
        best_b = 0
        if (n >= 6) then
            do a = 1, p-1
                do b = a+1, p
                    pairmax = 0.0_dp
                    do i = 1, n
                        if (i == pilot_index(a) .or. i == pilot_index(b)) cycle
                        pairmax = max(pairmax,abs(pilot(i,a)),abs(pilot(i,b)))
                    end do
                    if (pairmax <= sqrt(tiny(1.0_dp))) cycle
                    rfirst = 0
                    value = 0.0_dp
                    do i = 1, n
                        if (i == pilot_index(a) .or. i == pilot_index(b)) cycle
                        scale = max(abs(pilot(i,a)),abs(pilot(i,b)))
                        if (scale > value) then
                            value = scale
                            rfirst = i
                        end if
                    end do
                    if (rfirst == 0) cycle
                    rsecond = 0
                    value = 0.0_dp
                    do i = 1, n
                        if (i == pilot_index(a) .or. i == pilot_index(b) .or. i == rfirst) cycle
                        determinant = abs((pilot(rfirst,a)/pairmax)*(pilot(i,b)/pairmax) - &
                            (pilot(rfirst,b)/pairmax)*(pilot(i,a)/pairmax))
                        if (determinant > value) then
                            value = determinant
                            rsecond = i
                        end if
                    end do
                    if (rsecond /= 0 .and. value > best) then
                        best = value
                        best_a = a
                        best_b = b
                        r = [rfirst,rsecond]
                        c = [pilot_index(a),pilot_index(b)]
                        global_scale = pairmax
                        candidate = .true.
                    end if
                end do
            end do
        end if

        if (candidate) then
            pivot(1,1) = pilot(r(1),best_a)/global_scale
            pivot(1,2) = pilot(r(1),best_b)/global_scale
            pivot(2,1) = pilot(r(2),best_a)/global_scale
            pivot(2,2) = pilot(r(2),best_b)/global_scale
            determinant = pivot(1,1)*pivot(2,2)-pivot(1,2)*pivot(2,1)
            if (abs(determinant) <= pivot_tol) then
                candidate = .false.
            else
                do i = 1, n
                    if (any(i == c)) cycle
                    value = pilot(i,best_a)/global_scale
                    scale = pilot(i,best_b)/global_scale
                    left(i,1) = (value*pivot(2,2)-scale*pivot(2,1))/determinant
                    left(i,2) = (scale*pivot(1,1)-value*pivot(1,2))/determinant
                end do
                if (.not. all(ieee_is_finite(left))) candidate = .false.
                if (maxval(abs(left)) > 1.0e8_dp) candidate = .false.
            end if
        end if

        ! Every column is evaluated only once. For non-anchor rows, validate
        ! its off-diagonal elements immediately. Cache only two unresolved
        ! rows and two unresolved columns for the final exact checks.
        do j = 1, n
            trial = 0
            do k = 1, p
                if (pilot_index(k) == j) trial = k
            end do
            if (trial > 0) then
                product = pilot(:,trial)
            else
                basis(j) = 1.0_dp
                if (present(forward_calls)) forward_calls = forward_calls + 1
                call jacobian_vector(fun, x, fvec, basis, omega, product, &
                    nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                basis(j) = 0.0_dp
                if (.not. valid) return
            end if
            diagonal(j) = product(j)
            if (.not. candidate) cycle
            if (global_scale < 1.0_dp) then
                if (maxval(abs(product)) > global_scale*limit) then
                    candidate = .false.
                    cycle
                end if
            end if
            saved_rows(j,:) = product(c)
            if (j == r(1)) then
                special(:,1) = product
                cycle
            else if (j == r(2)) then
                special(:,2) = product
                cycle
            end if
            weights(j,:) = product(r)/global_scale
            do i = 1, n
                if (i == j .or. any(i == c)) cycle
                predicted = dot_product(left(i,:),weights(j,:))
                value = product(i)/global_scale
                if (.not. stream_rank_two_match(value,predicted)) candidate = .false.
                if (.not. candidate) exit
            end do
        end do

        ! Preserve the original exact diagonal fallback on every rejection.
        maxabs = maxval(abs(diagonal))
        if (.not. ieee_is_finite(maxabs)) return
        if (maxabs /= 0.0_dp) then
            cutoff = max(maxabs*floor_fraction,tiny(1.0_dp))
            do j = 1, n
                value = abs(diagonal(j))
                if (value == 0.0_dp) cycle
                invdiag(j) = sign(1.0_dp/max(value,cutoff),diagonal(j))
            end do
        end if
        ok = all(ieee_is_finite(invdiag))
        if (.not. ok .or. .not. candidate) return

        call finish_rank_two_stream(left, weights, saved_rows, special, &
            diagonal, global_scale, r, c, right, base, correction, verified)
        enabled = verified
    end subroutine discover_rank_two_inverse_stream

    pure logical function stream_rank_two_match(measured, predicted) result(matches)
        ! Use ratios to avoid overflow in the comparison of finite elements.
        real(dp), intent(in) :: measured !! measured normalized off-diagonal entry
        real(dp), intent(in) :: predicted !! fitted normalized off-diagonal entry
        real(dp) :: scale
        matches = .false.
        if (.not. ieee_is_finite(measured) .or. .not. ieee_is_finite(predicted)) return
        scale = max(abs(measured),abs(predicted),sqrt(tiny(1.0_dp)))
        matches = abs(measured/scale-predicted/scale) <= 1.0e-5_dp
    end function stream_rank_two_match

    subroutine finish_rank_two_stream(left, weights, saved_rows, special, diagonal, &
            maxoff, r, c, right, base, correction, enabled)
        ! Complete the four missing factor rows from independent off-diagonal
        ! observations; then verify every entry omitted from the streaming
        ! column checks before performing the guarded Woodbury construction.
        real(dp), intent(inout) :: left(:,:) !! fitted left factor, c rows missing
        real(dp), intent(inout) :: weights(:,:) !! fitted right factor, r rows missing
        real(dp), intent(in) :: saved_rows(:,:) !! J(c(:),j) for every column j
        real(dp), intent(in) :: special(:,:) !! J(:,r(:)) for missing right rows
        real(dp), intent(in) :: diagonal(:) !! exact measured Jacobian diagonal
        real(dp), intent(in) :: maxoff !! scaling used for factor reconstruction
        integer, intent(in) :: r(2) !! two independent anchor rows
        integer, intent(in) :: c(2) !! two independent anchor columns
        real(dp), intent(out) :: right(:,:) !! normalized fitted right factors
        real(dp), intent(out) :: base(:) !! reciprocal diagonal in Woodbury inverse
        real(dp), intent(out) :: correction(:,:) !! Woodbury correction columns
        logical, intent(out) :: enabled !! verified and numerically safe fit
        integer :: i, j, k, n, first, second
        real(dp) :: best, value, determinant, y1, y2, predicted
        real(dp), parameter :: pivot_tol = 1.0e-8_dp

        enabled = .false.
        right = 0.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        n = size(diagonal)
        do k = 1, 2
            i = c(k)
            first = 0
            best = 0.0_dp
            do j = 1, n
                if (j == i .or. any(j == r)) cycle
                value = max(abs(weights(j,1)),abs(weights(j,2)))
                if (value > best) then
                    best = value
                    first = j
                end if
            end do
            if (first == 0 .or. best <= pivot_tol) return
            second = 0
            best = pivot_tol
            do j = 1, n
                if (j == i .or. j == first .or. any(j == r)) cycle
                value = abs(weights(first,1)*weights(j,2)- &
                    weights(first,2)*weights(j,1))
                if (value > best) then
                    best = value
                    second = j
                end if
            end do
            if (second == 0) return
            determinant = weights(first,1)*weights(second,2) - &
                weights(first,2)*weights(second,1)
            y1 = saved_rows(first,k)/maxoff
            y2 = saved_rows(second,k)/maxoff
            left(i,1) = (y1*weights(second,2)-y2*weights(first,2))/determinant
            left(i,2) = (y2*weights(first,1)-y1*weights(second,1))/determinant
        end do
        if (.not. all(ieee_is_finite(left))) return
        if (maxval(abs(left)) > 1.0e8_dp) return

        do k = 1, 2
            j = r(k)
            first = 0
            best = 0.0_dp
            do i = 1, n
                if (i == j) cycle
                value = max(abs(left(i,1)),abs(left(i,2)))
                if (value > best) then
                    best = value
                    first = i
                end if
            end do
            if (first == 0 .or. best <= pivot_tol) return
            second = 0
            best = pivot_tol
            do i = 1, n
                if (i == j .or. i == first) cycle
                value = abs(left(first,1)*left(i,2)-left(first,2)*left(i,1))
                if (value > best) then
                    best = value
                    second = i
                end if
            end do
            if (second == 0) return
            determinant = left(first,1)*left(second,2)-left(first,2)*left(second,1)
            y1 = special(first,k)/maxoff
            y2 = special(second,k)/maxoff
            weights(j,1) = (y1*left(second,2)-y2*left(first,2))/determinant
            weights(j,2) = (y2*left(first,1)-y1*left(second,1))/determinant
        end do
        if (.not. all(ieee_is_finite(weights))) return
        if (maxval(abs(weights)) > 1.0e8_dp) return

        do j = 1, n
            do k = 1, 2
                i = c(k)
                if (i == j) cycle
                predicted = dot_product(left(i,:),weights(j,:))
                if (.not. stream_rank_two_match(saved_rows(j,k)/maxoff,predicted)) return
            end do
        end do
        do k = 1, 2
            j = r(k)
            do i = 1, n
                if (i == j) cycle
                predicted = dot_product(left(i,:),weights(j,:))
                if (.not. stream_rank_two_match(special(i,k)/maxoff,predicted)) return
            end do
        end do
        call inverse_rank_two_factors(diagonal, left, weights, maxoff, &
            right, base, correction, enabled)
    end subroutine finish_rank_two_stream

    subroutine inverse_rank_two_factors(diagonal, left, weights, maxoff, &
            right, base, correction, enabled, min_det)
        ! Woodbury inverse from checked normalized off-diagonal factors.
        real(dp), intent(in) :: diagonal(:) !! measured full Jacobian diagonal
        real(dp), intent(in) :: left(:,:) !! normalized learned left factor
        real(dp), intent(in) :: weights(:,:) !! normalized learned right factor
        real(dp), intent(in) :: maxoff !! off-diagonal magnitude used for scaling
        real(dp), intent(out) :: right(:,:) !! normalized learned right factor
        real(dp), intent(out) :: base(:) !! reciprocal uncoupled diagonal
        real(dp), intent(out) :: correction(:,:) !! Woodbury correction columns
        logical, intent(out) :: enabled !! safe Woodbury factorization
        real(dp), intent(in), optional :: min_det !! extra guard for fitted 2-by-2 inverse

        integer :: n, i, k
        real(dp) :: u(size(diagonal),2), a(size(diagonal)), work(size(diagonal),2)
        real(dp) :: small(2,2), normalized(2,2), cof(2,2)
        real(dp) :: vmax, scale, determinant, threshold
        real(dp), parameter :: max_allowed = sqrt(huge(1.0_dp))/64.0_dp

        enabled = .false.
        right = 0.0_dp
        base = 0.0_dp
        correction = 0.0_dp
        n = size(diagonal)
        u = left
        do k = 1, 2
            vmax = maxval(abs(weights(:,k)))
            if (vmax <= 1.0e-8_dp) return
            right(:,k) = weights(:,k)/vmax
            scale = maxoff*vmax
            if (scale <= 0.0_dp .or. scale > max_allowed) return
            if (maxval(abs(left(:,k))) > max_allowed/scale) return
            u(:,k) = left(:,k)*scale
        end do
        a = diagonal - u(:,1)*right(:,1) - u(:,2)*right(:,2)
        if (.not. all(ieee_is_finite(a))) return
        if (any(abs(a) <= sqrt(tiny(1.0_dp)))) return
        base = 1.0_dp/a
        do k = 1, 2
            do i = 1, n
                if (abs(base(i)) > 1.0_dp) then
                    if (abs(u(i,k)) > max_allowed/abs(base(i))) return
                else
                    if (abs(u(i,k)) > max_allowed) return
                end if
            end do
            work(:,k) = base*u(:,k)
        end do
        if (maxval(abs(work)) > max_allowed) return
        small = matmul(transpose(right),work)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        scale = max(1.0_dp,maxval(abs(small)))
        normalized = small/scale
        determinant = normalized(1,1)*normalized(2,2) - &
            normalized(1,2)*normalized(2,1)
        threshold = 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp,sum(abs(normalized(:,1)))*sum(abs(normalized(:,2))))
        if (present(min_det)) threshold = max(threshold,min_det)
        if (abs(determinant) <= threshold) return
        cof(1,1) = normalized(2,2)
        cof(1,2) = -normalized(1,2)
        cof(2,1) = -normalized(2,1)
        cof(2,2) = normalized(1,1)
        correction = matmul(work,cof)/(scale*determinant)
        if (.not. all(ieee_is_finite(correction))) return
        enabled = .true.
    end subroutine inverse_rank_two_factors

    subroutine estimate_row_coupled_inverse(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, diagonal, base, correction, &
            enabled, ok, jacobian_product, right, forward_calls)
        ! Optional checked inverse for J = diag(a) + u*right^T.
        ! The bounded sign(right) probe identifies each row's u from the
        ! exact Jacobian diagonal. Two independent sign probes then test
        ! the model before enabling a Sherman-Morrison inverse.
        ! A rejected or ill-conditioned fit retains diagonal scaling.
        procedure(nonlin_function) :: fun !! nonlinear residual callback
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(in) :: fvec(:) !! residual at x
        real(dp), intent(in) :: omega !! directional finite-difference scale
        integer, intent(inout) :: nfev !! residual evaluation count
        logical, intent(in) :: use_central !! use centered differences
        integer, intent(in) :: jac_step_mode !! perturbation scaling rule
        real(dp), intent(in) :: rdiff_use !! coordinate relative step
        real(dp), intent(in) :: diagonal(:) !! exact unclipped Jacobian diagonal
        real(dp), intent(out) :: base(:) !! reciprocal uncoupled diagonal
        real(dp), intent(out) :: correction(:) !! Woodbury rank-one factor
        logical, intent(out) :: enabled !! true only if model is verified and nonsingular
        logical, intent(out) :: ok !! false if any requested Jacobian product fails
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! exact Jv action
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products
        real(dp), intent(in) :: right(:) !! normalized fixed rank-one right factor

        real(dp) :: direction(size(x)), observed(size(x)), u(size(x)), a(size(x))
        real(dp) :: fit_direction(size(x)), first_check(size(x))
        real(dp) :: model(size(x)), work(size(x)), bound, denom, threshold, scale
        real(dp) :: mass, row_mass, projection
        integer :: i, check_index, flip_index
        logical :: valid

        enabled = .false.
        ok = .false.
        base = 0.0_dp
        correction = 0.0_dp
        if (size(x) <= 1) then
            ok = .true.
            return
        end if
        ! The right weights are normalized by the caller to max(abs(v))=1.
        ! Using sign(v), rather than v, avoids tiny probe components.
        direction = 0.0_dp
        where (right > 0.0_dp) direction = 1.0_dp
        where (right < 0.0_dp) direction = -1.0_dp
        mass = sum(abs(right))
        fit_direction = direction
        if (present(forward_calls)) forward_calls = forward_calls + 1
        call jacobian_vector(fun, x, fvec, direction, omega, observed, nfev, &
            valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
        if (.not. valid) return
        ok = .true.
        bound = max(maxval(abs(observed)), maxval(abs(diagonal)))
        if (bound > sqrt(huge(1.0_dp)) .or. bound == 0.0_dp) return
        do i = 1, size(x)
            ! Jii = ai + ui*vi, while J*sign(v) contains
            ! ai*sign(vi) + ui*sum(abs(v)). Their difference isolates ui.
            row_mass = mass - abs(right(i))
            if (row_mass <= sqrt(epsilon(1.0_dp)) * mass) return
            u(i) = (observed(i) - diagonal(i)*direction(i)) / row_mass
        end do
        a = diagonal - u*right
        if (.not. all(ieee_is_finite(a))) return
        if (any(abs(a) < sqrt(tiny(1.0_dp)))) return
        base = 1.0_dp / a
        ! Reserve summation headroom, including when IEEE traps are enabled.
        do i = 1, size(x)
            if (abs(base(i)) > 1.0_dp) then
                if (abs(u(i)) > (huge(1.0_dp) / (4.0_dp * real(size(x), dp))) / &
                        abs(base(i))) return
            else
                if (abs(u(i)) > huge(1.0_dp) / (4.0_dp * real(size(x), dp))) return
            end if
        end do
        work = base * u
        bound = sum(abs(work * right))
        if (.not. ieee_is_finite(bound)) return
        denom = 1.0_dp + sum(work * right)
        threshold = 100.0_dp * epsilon(1.0_dp) * (1.0_dp + bound)
        if (abs(denom) <= threshold .or. .not. ieee_is_finite(denom)) return
        correction = work / denom
        if (.not. all(ieee_is_finite(correction))) return

        do check_index = 1, 2
            do i = 1, size(x)
                if (check_index == 1) then
                    if (modulo(i, 2) == 1) then
                        direction(i) = 1.0_dp
                    else
                        direction(i) = -1.0_dp
                    end if
                else
                    if (modulo(i - 1, 4) < 2) then
                        direction(i) = 1.0_dp
                    else
                        direction(i) = -1.0_dp
                    end if
                end if
            end do
            ! For signed right factors, a fixed sign check can be identical
            ! (or opposite) to the fitting direction or to another check.
            ! Rotate by deterministic coordinate flips to avoid redundant
            ! checks when there are at least three variables.
            if (size(x) > 2) then
                do flip_index = 1, size(x)
                    if (any(direction /= fit_direction) .and. &
                            any(direction /= -fit_direction)) then
                        if (check_index == 1) exit
                        if (any(direction /= first_check) .and. &
                                any(direction /= -first_check)) exit
                    end if
                    direction(flip_index) = -direction(flip_index)
                end do
            end if
            if (check_index == 1) first_check = direction
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, direction, omega, observed, nfev, &
                valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
            projection = sum(right * direction)
            model = a * direction + u * projection
            if (.not. all(ieee_is_finite(model))) return
            do i = 1, size(x)
                scale = max(abs(model(i)), abs(observed(i)), &
                    abs(diagonal(i)), sqrt(tiny(1.0_dp)))
                if (abs(model(i)/scale - observed(i)/scale) > 1.0e-5_dp) return
            end do
        end do
        enabled = .true.
    end subroutine estimate_row_coupled_inverse

    subroutine rank_two_sign_candidates(right, candidates)
        real(dp), intent(in) :: right(:,:) !! two normalized right coupling factors
        real(dp), intent(out) :: candidates(:,:) !! twelve deterministic bounded sign vectors
        integer :: i

        candidates(:,1) = 1.0_dp
        candidates(:,2) = merge(1.0_dp,-1.0_dp,right(:,1) >= 0.0_dp)
        candidates(:,3) = merge(1.0_dp,-1.0_dp,right(:,2) >= 0.0_dp)
        candidates(:,4) = merge(1.0_dp,-1.0_dp,right(:,1)+right(:,2) >= 0.0_dp)
        candidates(:,5) = merge(1.0_dp,-1.0_dp,right(:,1)-right(:,2) >= 0.0_dp)
        do i = 1, size(right,1)
            candidates(i,6) = merge(1.0_dp,-1.0_dp,modulo(i-1,2) == 0)
            candidates(i,7) = merge(1.0_dp,-1.0_dp,modulo(i-1,4) < 2)
            candidates(i,8) = merge(1.0_dp,-1.0_dp,modulo(i-1,6) < 3)
            candidates(i,9) = merge(1.0_dp,-1.0_dp,modulo(i-1,8) < 4)
            candidates(i,10) = merge(1.0_dp,-1.0_dp,modulo(3*i+1,7) < 3)
            candidates(i,11) = merge(1.0_dp,-1.0_dp,modulo(5*i+2,11) < 5)
            candidates(i,12) = merge(1.0_dp,-1.0_dp,modulo(7*i+3,13) < 6)
        end do
    end subroutine rank_two_sign_candidates

    pure function determinant_three(a) result(det)
        real(dp), intent(in) :: a(3,3) !! row-major coefficients of a small linear solve
        real(dp) :: det !! determinant of the 3-by-3 matrix

        det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - &
            a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + &
            a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    end function determinant_three

    subroutine estimate_rank_two_diagonal_sketch(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, right, diagonal, available, &
            ok, jacobian_product, forward_calls)
        ! Recover Jii in a D + U*V^T model using three bounded sign probes,
        ! rather than n exact coordinate products. Each row has three unknowns:
        ! Dii, Ui1, Ui2. The later guarded Woodbury estimator uses two more
        ! fitting and two verifying Jv products. If either stage fails, callers
        ! retain the exact-coordinate-diagonal fallback.
        procedure(nonlin_function) :: fun !! square residual callback
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(in) :: fvec(:) !! residual vector at the current point
        real(dp), intent(in) :: omega !! numerical Jacobian step scale
        integer, intent(inout) :: nfev !! number of residual callback evaluations
        logical, intent(in) :: use_central !! use a centered derivative
        integer, intent(in) :: jac_step_mode !! finite-difference step policy
        real(dp), intent(in) :: rdiff_use !! requested componentwise difference step
        real(dp), intent(in) :: right(:,:) !! normalized right factors, shape (n,2)
        real(dp), intent(out) :: diagonal(:) !! reconstructed exact-model Jii
        logical, intent(out) :: available !! true if all rows were identifiable
        logical, intent(out) :: ok !! false if a Jacobian callback fails
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! analytic Jv callback
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products

        integer, parameter :: nc = 12
        real(dp) :: candidates(size(x),nc), masses(nc,2)
        real(dp) :: observations(size(x),3), fit(size(x),3)
        real(dp) :: h(3,3), system(3,3), rhs(3)
        real(dp) :: best, worst, score, determinant, scale_m(2)
        real(dp) :: coefficients(3), left(3), max_allowed
        integer :: i, j, k, l, m, selected(3), indices(3)
        logical :: valid

        available = .false.
        ok = .true.
        diagonal = 0.0_dp
        if (size(x) < 5) return
        call rank_two_sign_candidates(right,candidates)
        masses = matmul(transpose(candidates),right)
        scale_m = maxval(abs(masses),dim=1)
        if (any(scale_m <= 1.0e-12_dp)) return
        masses(:,1) = masses(:,1)/scale_m(1)
        masses(:,2) = masses(:,2)/scale_m(2)

        ! Prefer the triple having the best worst-case row geometry. A
        ! relative determinant rejects triples where a row's three
        ! unknowns cannot be distinguished from bounded measurements.
        best = 0.0_dp
        selected = 0
        do j = 1, nc-2
            do k = j+1, nc-1
                do l = k+1, nc
                    indices = [j,k,l]
                    worst = 1.0_dp
                    do i = 1, size(x)
                        do m = 1, 3
                            h(m,:) = [candidates(i,indices(m)), masses(indices(m),:)]
                        end do
                        left(1) = sum(abs(h(1,:)))
                        left(2) = sum(abs(h(2,:)))
                        left(3) = sum(abs(h(3,:)))
                        determinant = determinant_three(h)
                        score = abs(determinant)/product(left)
                        worst = min(worst,score)
                        if (worst <= best) exit
                    end do
                    if (worst > best) then
                        best = worst
                        selected = indices
                    end if
                end do
            end do
        end do
        if (best < 1.0e-5_dp) return
        fit = candidates(:,selected)
        do j = 1, 3
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, fit(:,j), omega, observations(:,j), &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
        end do
        max_allowed = sqrt(huge(1.0_dp))/64.0_dp
        if (maxval(abs(observations)) > max_allowed) return
        do i = 1, size(x)
            do j = 1, 3
                h(j,:) = [fit(i,j), masses(selected(j),:)]
            end do
            determinant = determinant_three(h)
            if (abs(determinant) <= 1.0e-6_dp*product( &
                    [sum(abs(h(1,:))),sum(abs(h(2,:))),sum(abs(h(3,:)))]) ) return
            rhs = observations(i,:)
            do j = 1, 3
                system = h
                system(:,j) = rhs
                coefficients(j) = determinant_three(system)/determinant
            end do
            ! Moments were scaled in the fitted system; undo the scales
            ! before reconstructing the true Jacobian diagonal.
            if (maxval(abs(coefficients)) > max_allowed) return
            diagonal(i) = coefficients(1) + &
                coefficients(2)/scale_m(1)*right(i,1) + &
                coefficients(3)/scale_m(2)*right(i,2)
        end do
        if (.not. all(ieee_is_finite(diagonal))) return
        if (maxval(abs(diagonal)) > max_allowed) return
        available = .true.
    end subroutine estimate_rank_two_diagonal_sketch

    subroutine estimate_rank_two_inverse(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, diagonal, right, base, &
            correction, enabled, ok, jacobian_product, forward_calls)
        ! Fit J = diag(a) + U*transpose(V) with two known right factors.
        ! Exact Jii and two bounded probes identify each row of U. A pair of
        ! further independent probes checks the model. Bad fits, degenerate
        ! weight patterns and unsafe Woodbury inverses retain diagonal scaling.
        procedure(nonlin_function) :: fun !! square nonlinear residual callback
        real(dp), intent(in) :: x(:) !! point at which the Jacobian is evaluated
        real(dp), intent(in) :: fvec(:) !! residual at the current point
        real(dp), intent(in) :: omega !! directional differencing scale
        integer, intent(inout) :: nfev !! residual evaluation count
        logical, intent(in) :: use_central !! centered differences requested
        integer, intent(in) :: jac_step_mode !! finite-difference step rule
        real(dp), intent(in) :: rdiff_use !! coordinatewise difference step
        real(dp), intent(in) :: diagonal(:) !! exact, unclipped Jacobian diagonal
        real(dp), intent(in) :: right(:,:) !! normalized known right factors (n,2)
        real(dp), intent(out) :: base(:) !! inverse diagonal of fitted D
        real(dp), intent(out) :: correction(:,:) !! precomputed Woodbury correction (n,2)
        logical, intent(out) :: enabled !! true only after both checks pass
        logical, intent(out) :: ok !! false if a requested Jacobian action fails
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! exact Jacobian action
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products

        integer, parameter :: nc = 12
        real(dp) :: candidates(size(x),nc), fit(size(x),2)
        real(dp) :: observations(size(x),2), basis(size(x),4)
        real(dp) :: direction(size(x)), observed(size(x))
        real(dp) :: u(size(x),2), a(size(x)), work(size(x),2), model(size(x))
        real(dp) :: mass(nc,2), h(2,2), normalized(2,2), small(2,2)
        real(dp) :: projection(2), cof(2,2), residual(size(x))
        real(dp) :: best, score, worst, determinant, scale, norm1, norm2
        real(dp) :: tolerance, max_allowed, dot_value, norm_value
        integer :: i, j, k, left, right_index, selected(2), choice
        logical :: valid

        enabled = .false.
        ok = .true.
        base = 0.0_dp
        correction = 0.0_dp
        if (size(x) < 4) return

        ! These fixed signs have bounded magnitude and require no random state.
        ! Choose the fitting pair maximizing its worst per-row separation.
        call rank_two_sign_candidates(right,candidates)
        mass = matmul(transpose(candidates),right)
        best = 0.0_dp
        selected = 0
        do left = 1, nc-1
            do right_index = left+1, nc
                worst = 1.0_dp
                do i = 1, size(x)
                    h(1,:) = mass(left,:) - candidates(i,left)*right(i,:)
                    h(2,:) = mass(right_index,:) - candidates(i,right_index)*right(i,:)
                    norm1 = sum(abs(h(1,:)))
                    norm2 = sum(abs(h(2,:)))
                    if (min(norm1,norm2) <= tiny(1.0_dp)) then
                        worst = 0.0_dp
                        exit
                    end if
                    determinant = h(1,1)*h(2,2)-h(1,2)*h(2,1)
                    score = abs(determinant)/(norm1*norm2)
                    worst = min(worst,score)
                    if (worst <= best) exit
                end do
                if (worst > best) then
                    best = worst
                    selected = [left,right_index]
                end if
            end do
        end do
        ! Poorly identifiable rows cannot distinguish a from U using Jii.
        if (best <= 1.0e-8_dp) return
        fit(:,1) = candidates(:,selected(1))
        fit(:,2) = candidates(:,selected(2))
        do j = 1, 2
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, fit(:,j), omega, observations(:,j), &
                nfev, valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
        end do
        max_allowed = sqrt(huge(1.0_dp))/64.0_dp
        if (max(maxval(abs(diagonal)),maxval(abs(observations))) > max_allowed) return
        do i = 1, size(x)
            h(1,:) = mass(selected(1),:) - fit(i,1)*right(i,:)
            h(2,:) = mass(selected(2),:) - fit(i,2)*right(i,:)
            scale = maxval(abs(h))
            ! Tiny off-diagonal moment coefficients are not identifiable
            ! reliably and can overflow the subsequent rowwise solve.
            if (scale <= 1.0e-8_dp) return
            normalized = h/scale
            determinant = normalized(1,1)*normalized(2,2) - &
                normalized(1,2)*normalized(2,1)
            if (abs(determinant) <= 1.0e-8_dp) return
            ! Fit the two off-diagonal row contributions without subtracting
            ! the inferred rank-two diagonal a second time.
            projection = (observations(i,:) - diagonal(i)*fit(i,:))/scale
            if (maxval(abs(projection)) > 0.25_dp*max_allowed*abs(determinant)) return
            u(i,1) = (projection(1)*normalized(2,2) - &
                projection(2)*normalized(1,2))/determinant
            u(i,2) = (normalized(1,1)*projection(2) - &
                normalized(2,1)*projection(1))/determinant
        end do
        if (.not. all(ieee_is_finite(u))) return
        if (maxval(abs(u)) > max_allowed) return
        a = diagonal - u(:,1)*right(:,1) - u(:,2)*right(:,2)
        if (.not. all(ieee_is_finite(a))) return
        if (any(abs(a) <= sqrt(tiny(1.0_dp)))) return
        base = 1.0_dp/a
        do j = 1, 2
            do i = 1, size(x)
                if (abs(base(i)) > 1.0_dp) then
                    if (abs(u(i,j)) > max_allowed/abs(base(i))) return
                else
                    if (abs(u(i,j)) > max_allowed) return
                end if
            end do
            work(:,j) = base*u(:,j)
        end do
        if (maxval(abs(work)) > max_allowed) return

        ! M = I + V^T D^-1 U; correction = D^-1 U M^-1.
        ! Normalize the small matrix before inverting to guard its determinant.
        small = matmul(transpose(right),work)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        scale = max(1.0_dp,maxval(abs(small)))
        normalized = small/scale
        determinant = normalized(1,1)*normalized(2,2) - &
            normalized(1,2)*normalized(2,1)
        tolerance = 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp,sum(abs(normalized(:,1)))*sum(abs(normalized(:,2))))
        if (abs(determinant) <= tolerance) return
        cof(1,1) = normalized(2,2)
        cof(1,2) = -normalized(1,2)
        cof(2,1) = -normalized(2,1)
        cof(2,2) = normalized(1,1)
        correction = matmul(work,cof)/(scale*determinant)
        if (.not. all(ieee_is_finite(correction))) return

        ! Orthogonalize sign checks against the two fit directions and each
        ! other. A check inside their span cannot detect a fitting error.
        do j = 1, 4
            if (j <= 2) then
                direction = fit(:,j)
            else
                best = 0.0_dp
                choice = 0
                do k = 1, nc
                    residual = candidates(:,k)
                    do i = 1, j-1
                        dot_value = dot_product(basis(:,i),residual)
                        residual = residual - dot_value*basis(:,i)
                    end do
                    norm_value = sqrt(sum(residual**2))
                    if (norm_value > best) then
                        best = norm_value
                        choice = k
                    end if
                end do
                if (choice == 0) return
                direction = candidates(:,choice)
            end if
            residual = direction
            do i = 1, j-1
                residual = residual - dot_product(residual,basis(:,i))*basis(:,i)
            end do
            norm_value = sqrt(sum(residual**2))
            if (norm_value < 1.0e-5_dp*sqrt(real(size(x),dp))) return
            basis(:,j) = residual/norm_value
            if (j <= 2) cycle
            if (present(forward_calls)) forward_calls = forward_calls + 1
            call jacobian_vector(fun, x, fvec, direction, omega, observed, nfev, &
                valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. valid) then
                ok = .false.
                return
            end if
            projection = matmul(transpose(right),direction)
            model = a*direction + u(:,1)*projection(1) + u(:,2)*projection(2)
            if (.not. all(ieee_is_finite(model))) return
            do i = 1, size(x)
                scale = max(abs(model(i)),abs(observed(i)),abs(diagonal(i)), &
                    sqrt(tiny(1.0_dp)))
                if (abs(model(i)/scale - observed(i)/scale) > 1.0e-5_dp) return
            end do
        end do
        enabled = .true.
    end subroutine estimate_rank_two_inverse

    subroutine estimate_inverse_diagonal(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, floor_fraction, positive_only, &
            invdiag, ok, jacobian_product, nprobes, probe_scheme, probe_seed, &
            verify_samples, verify_rtol, probe_max, raw_diagonal, jacobian_matrix, forward_calls)
        ! Exact mode uses n coordinate Jv products. Walsh-sign probing can
        ! recover narrow-banded diagonals exactly; seeded pseudo-Rademacher
        ! signs avoid fixed Walsh alias patterns but remain approximate.
        ! Verification combines probe-batch consistency (for seeded signs)
        ! and coordinate spot checks. In either case a failed check causes
        ! exact fallback. These are diagnostics, not accuracy guarantees.
        ! A zero diagonal uses identity scaling; a tiny nonzero diagonal is
        ! clipped relative to the largest observed diagonal magnitude.
        procedure(nonlin_function) :: fun !! square residual callback
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(in) :: fvec(:) !! residual at x
        real(dp), intent(in) :: omega !! current directional finite-difference scale
        integer, intent(inout) :: nfev !! number of residual evaluations
        logical, intent(in) :: use_central !! use centered finite differences
        integer, intent(in) :: jac_step_mode !! finite-difference step selection
        real(dp), intent(in) :: rdiff_use !! coordinate relative perturbation
        real(dp), intent(in) :: floor_fraction !! nonzero-diagonal relative floor
        logical, intent(in) :: positive_only !! MINRES requires positive preconditioning
        real(dp), intent(out) :: invdiag(:) !! estimated inverse diagonal
        logical, intent(out) :: ok !! true for valid finite estimates
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! exact Jv action
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products
        integer, intent(in) :: nprobes !! zero exact; positive less than n uses sign probes
        character(len=*), intent(in) :: probe_scheme !! walsh or rademacher
        integer, intent(in) :: probe_seed !! reproducible seed, ignored for Walsh
        integer, intent(in) :: verify_samples !! evenly spaced coordinate spot checks
        real(dp), intent(in) :: verify_rtol !! relative discrepancy trigger
        integer, intent(in) :: probe_max !! zero preserves the fixed-budget estimator
        real(dp), intent(out), optional :: raw_diagonal(:) !! unclipped J diagonal estimate
        real(dp), intent(out), optional :: jacobian_matrix(:,:) !! exact J columns, when requested

        real(dp), allocatable :: basis(:), product(:), raw(:), first_group(:)
        real(dp), allocatable :: checked_values(:)
        logical, allocatable :: checked_mask(:)
        real(dp) :: maxabs, cutoff, value, weight, scale
        integer :: j, k, sample, index
        integer(int64) :: prng_state
        logical :: valid, fallback_exact


        ok = .false.
        invdiag = 1.0_dp
        allocate(basis(size(x)), product(size(x)), raw(size(x)))
        basis = 0.0_dp
        raw = 0.0_dp
        fallback_exact = (nprobes == 0 .or. nprobes >= size(x))
        if (.not. fallback_exact .and. probe_max > nprobes) then
            ! Progressive stages sample different coordinates, but all use
            ! the same Jacobian. Retain exact values for later verification
            ! and for fallback to the full diagonal.
            allocate(checked_mask(size(x)), checked_values(size(x)))
            call estimate_progressive_probes(fun, x, fvec, omega, nfev, use_central, &
                jac_step_mode, rdiff_use, raw, basis, product, valid, fallback_exact, &
                jacobian_product, nprobes, probe_max, probe_seed, verify_samples, verify_rtol, &
                checked_mask, checked_values, forward_calls=forward_calls)
            if (.not. valid) return
        else if (.not. fallback_exact) then
            ! For verified seeded probes, compare the first half of the
            ! signs against the full estimate. This adds only O(n) work and
            ! storage, with no extra Jacobian products. Substantial
            ! disagreement implies high estimator variance, so use exact
            ! coordinate probes rather than an unstable preconditioner.
            if (probe_scheme == 'rademacher' .and. verify_samples > 0 .and. &
                    nprobes > 1) then
                allocate(first_group(size(x)))
                first_group = 0.0_dp
            end if
            ! The seed and Park-Miller state are local to this call: no
            ! dependence on RANDOM_NUMBER or global mutable state. Arithmetic
            ! is performed in int64 so the recurrence cannot overflow.
            prng_state = 1_int64 + modulo(int(probe_seed, int64), 2147483646_int64)
            weight = 1.0_dp / real(nprobes, dp)
            do j = 1, nprobes
                do k = 1, size(x)
                    if (probe_scheme == 'walsh') then
                        if (btest(popcnt(iand(k - 1, j - 1)), 0)) then
                            basis(k) = -1.0_dp
                        else
                            basis(k) = 1.0_dp
                        end if
                    else
                        prng_state = modulo(48271_int64 * prng_state, 2147483647_int64)
                        if (prng_state < 1073741824_int64) then
                            basis(k) = -1.0_dp
                        else
                            basis(k) = 1.0_dp
                        end if
                    end if
                end do
                if (present(forward_calls)) forward_calls = forward_calls + 1
                call jacobian_vector(fun, x, fvec, basis, omega, product, nfev, &
                    valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                if (.not. valid) return
                ! Scale before summation to avoid overflow of finite products.
                raw = raw + basis * (product * weight)
                if (.not. all(ieee_is_finite(raw))) return
                if (allocated(first_group)) then
                    if (j <= nprobes / 2) then
                        first_group = first_group + basis * &
                            (product / real(nprobes / 2, dp))
                        if (.not. all(ieee_is_finite(first_group))) return
                    end if
                end if
            end do
            if (allocated(first_group)) then
                do k = 1, size(x)
                    scale = max(abs(raw(k)), abs(first_group(k)), tiny(1.0_dp))
                    if (abs(raw(k)/scale - first_group(k)/scale) > verify_rtol) then
                        fallback_exact = .true.
                        exit
                    end if
                end do
            end if

            ! Spot-check a uniform sample of exact coordinate diagonals.
            ! This detects some biased probe estimates, but not errors in
            ! unexamined coordinates. If any check disagrees by >25% of its
            ! local magnitude, rebuild the entire diagonal exactly.
            basis = 0.0_dp
            do sample = 1, verify_samples
                if (fallback_exact) exit
                index = 1 + int(real(sample - 1, dp) * real(size(x), dp) / &
                    real(verify_samples, dp))
                index = min(index, size(x))
                basis(index) = 1.0_dp
                if (present(forward_calls)) forward_calls = forward_calls + 1
                call jacobian_vector(fun, x, fvec, basis, omega, product, nfev, &
                    valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                basis(index) = 0.0_dp
                if (.not. valid) return
                scale = max(abs(product(index)), abs(raw(index)), tiny(1.0_dp))
                if (abs(product(index)/scale - raw(index)/scale) > verify_rtol) then
                    fallback_exact = .true.
                    exit
                end if
            end do
        end if
        if (fallback_exact) then
            ! The old fixed-budget path still evaluates every coordinate.
            ! In progressive mode reuse exact spot products already checked
            ! at this nonlinear iterate, with no additional Jv calls.
            basis = 0.0_dp
            do j = 1, size(x)
                if (allocated(checked_mask)) then
                    if (checked_mask(j)) then
                        raw(j) = checked_values(j)
                        cycle
                    end if
                end if
                basis(j) = 1.0_dp
                if (present(forward_calls)) forward_calls = forward_calls + 1
                call jacobian_vector(fun, x, fvec, basis, omega, product, nfev, &
                    valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                basis(j) = 0.0_dp
                if (.not. valid) return
                raw(j) = product(j)
                if (present(jacobian_matrix)) jacobian_matrix(:,j) = product
            end do
        end if
        maxabs = maxval(abs(raw))
        if (.not. ieee_is_finite(maxabs)) return
        if (present(raw_diagonal)) raw_diagonal = raw
        if (maxabs == 0.0_dp) then
            ok = .true.
            return
        end if
        cutoff = max(maxabs * floor_fraction, tiny(1.0_dp))
        do j = 1, size(x)
            value = abs(raw(j))
            if (value == 0.0_dp) cycle
            value = max(value, cutoff)
            if (positive_only) then
                invdiag(j) = 1.0_dp / value
            else
                invdiag(j) = sign(1.0_dp / value, raw(j))
            end if
        end do
        ok = all(ieee_is_finite(invdiag))
    end subroutine estimate_inverse_diagonal

    subroutine estimate_progressive_probes(fun, x, fvec, omega, nfev, &
            use_central, jac_step_mode, rdiff_use, raw, basis, product, ok, &
            fallback_exact, jacobian_product, nprobes, probe_max, probe_seed, &
            verify_samples, verify_rtol, checked_mask, checked_values, forward_calls)
        ! Grow verified seeded probes incrementally, reusing every sampled Jv.
        ! Compare earlier and expanded averages; failure at the budget limit
        ! triggers exact coordinate fallback in the calling estimator.
        procedure(nonlin_function) :: fun !! nonlinear residual
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: fvec(:) !! residual at x
        real(dp), intent(in) :: omega !! directional difference scale
        integer, intent(inout) :: nfev !! residual evaluation count
        logical, intent(in) :: use_central !! use centered differences
        integer, intent(in) :: jac_step_mode !! difference scaling mode
        real(dp), intent(in) :: rdiff_use !! relative difference scale
        real(dp), intent(out) :: raw(:) !! estimated Jacobian diagonal
        real(dp), intent(inout) :: basis(:) !! signed probe/coordinate workspace
        real(dp), intent(inout) :: product(:) !! Jacobian product workspace
        logical, intent(out) :: ok !! all products finite and valid
        logical, intent(out) :: fallback_exact !! estimate failed checks at max probes
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! analytic Jv
        integer, intent(inout), optional :: forward_calls !! attempted setup Jv products
        integer, intent(in) :: nprobes !! initial sign-probe budget
        integer, intent(in) :: probe_max !! maximum number of sign probes
        integer, intent(in) :: probe_seed !! local deterministic seed
        integer, intent(in) :: verify_samples !! coordinate checks per accepted batch
        real(dp), intent(in) :: verify_rtol !! normalized discrepancy cutoff
        logical, intent(out) :: checked_mask(:) !! true where an exact diagonal value is cached
        real(dp), intent(out) :: checked_values(:) !! exact diagonal at checked coordinates

        real(dp), allocatable :: first_group(:)
        integer, allocatable :: next_free(:)
        real(dp) :: denom, fraction, scale
        integer :: target, count, j, k, sample, index, stage, phase
        integer(int64) :: offset, numerator, denominator
        integer(int64) :: prng_state
        logical :: valid, disagreement

        allocate(first_group(size(x)), next_free(size(x) + 1))
        ! A local successor table finds an unchecked coordinate without a
        ! potentially quadratic scan when verification covers most of x.
        ! All entries initially point to themselves; after sampling an entry
        ! points to the next available coordinate (n+1 is the sentinel).
        do j = 1, size(next_free)
            next_free(j) = j
        end do
        raw = 0.0_dp
        first_group = 0.0_dp
        prng_state = 1_int64 + modulo(int(probe_seed, int64), 2147483646_int64)
        ok = .false.
        fallback_exact = .true.
        checked_mask = .false.
        target = nprobes
        count = 0
        stage = 0
        do
            do j = count + 1, target
                do k = 1, size(x)
                    prng_state = modulo(48271_int64 * prng_state, 2147483647_int64)
                    if (prng_state < 1073741824_int64) then
                        basis(k) = -1.0_dp
                    else
                        basis(k) = 1.0_dp
                    end if
                end do
                if (present(forward_calls)) forward_calls = forward_calls + 1
                call jacobian_vector(fun, x, fvec, basis, omega, product, nfev, &
                    valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                if (.not. valid) return
                ! A running average prevents large intermediate sums when the
                ! probe budget increases; never reconstruct earlier products.
                fraction = 1.0_dp / real(j, dp)
                raw = raw * (1.0_dp - fraction) + basis * (product * fraction)
                if (.not. all(ieee_is_finite(raw))) return
                if (count == 0 .and. j == max(1, nprobes / 2)) then
                    first_group = raw
                end if
            end do
            disagreement = .false.
            if (target > 1) then
                do k = 1, size(x)
                    scale = max(abs(raw(k)), abs(first_group(k)), tiny(1.0_dp))
                    if (abs(raw(k)/scale - first_group(k)/scale) > verify_rtol) then
                        disagreement = .true.
                        exit
                    end if
                end do
            end if
            if (.not. disagreement) then
                basis = 0.0_dp
                ! Recheck all previously visited coordinates against the
                ! expanded estimate without repeating their Jv products.
                ! A formerly bad coordinate must not be silently forgotten.
                do index = 1, size(x)
                    if (.not. checked_mask(index)) cycle
                    denom = max(abs(checked_values(index)), abs(raw(index)), tiny(1.0_dp))
                    if (abs(checked_values(index)/denom - raw(index)/denom) > verify_rtol) then
                        disagreement = .true.
                        exit
                    end if
                end do
                if (.not. disagreement) then
                    ! Van der Corput offsets (0, 1/2, 1/4, 3/4, ...)
                    ! spread successive samples through each stratum.
                    ! Use integer arithmetic for reproducibility and to
                    ! avoid overflow in large problem dimensions.
                    phase = stage
                    numerator = 0_int64
                    denominator = 1_int64
                    do while (phase > 0)
                        numerator = 2_int64 * numerator + int(modulo(phase, 2), int64)
                        denominator = 2_int64 * denominator
                        phase = phase / 2
                    end do
                    offset = (int(size(x) / verify_samples, int64) * numerator) / denominator
                    do sample = 1, verify_samples
                        index = 1 + int(real(sample - 1, dp) * real(size(x), dp) / &
                            real(verify_samples, dp))
                        index = min(index, size(x))
                        index = 1 + int(modulo(int(index - 1, int64) + offset, &
                            int(size(x), int64)))
                        ! The bit-reversed offset can round to zero when
                        ! verify_samples is large, causing entire stages to
                        ! revisit the same sites. Prefer the next unchecked
                        ! coordinate in cyclic order. Previously sampled
                        ! exact values were compared above without new Jv calls.
                        index = next_unchecked(index)
                        if (index > size(x)) index = next_unchecked(1)
                        if (index > size(x)) exit
                        basis(index) = 1.0_dp
                        if (present(forward_calls)) forward_calls = forward_calls + 1
                        call jacobian_vector(fun, x, fvec, basis, omega, product, nfev, &
                            valid, use_central, jac_step_mode, rdiff_use, jacobian_product)
                        basis(index) = 0.0_dp
                        if (.not. valid) return
                        checked_mask(index) = .true.
                        checked_values(index) = product(index)
                        next_free(index) = next_unchecked(index + 1)
                        denom = max(abs(checked_values(index)), abs(raw(index)), tiny(1.0_dp))
                        if (abs(checked_values(index)/denom - raw(index)/denom) > verify_rtol) then
                            disagreement = .true.
                            exit
                        end if
                    end do
                end if
            end if
            if (.not. disagreement) then
                fallback_exact = .false.
                ok = .true.
                return
            end if
            if (target >= probe_max) exit
            count = target
            stage = stage + 1
            first_group = raw
            ! Avoid integer overflow even for unusually large user budgets.
            if (target > probe_max / 2) then
                target = probe_max
            else
                target = 2 * target
            end if
        end do
        ok = .true.

    contains

        integer function next_unchecked(start) result(found)
            ! Path-compressed successor lookup; at most one Jv per coordinate.
            integer, intent(in) :: start !! first candidate in [1,n+1]
            integer :: node, successor

            node = start
            do while (next_free(node) /= node)
                node = next_free(node)
            end do
            found = node
            node = start
            do while (next_free(node) /= node)
                successor = next_free(node)
                next_free(node) = found
                node = successor
            end do
        end function next_unchecked

    end subroutine estimate_progressive_probes

    pure real(dp) function safe_two_norm(values) result(value)
        ! Avoid underflow/overflow in NORM2 for extreme finite magnitudes.
        ! Preserve intrinsic NORM2 rounding for normally scaled vectors.
        real(dp), intent(in) :: values(:) !! input vector

        real(dp) :: maxabs

        maxabs = maxval(abs(values))
        if (maxabs == 0.0_dp) then
            value = 0.0_dp
        else if (maxabs < sqrt(tiny(1.0_dp)) .or. &
                maxabs > sqrt(huge(1.0_dp)) / sqrt(real(size(values), dp))) then
            value = maxabs * sqrt(sum((values / maxabs)**2))
        else
            value = norm2(values)
        end if
    end function safe_two_norm

    pure logical function inner_residual_acceptable(rhs, product, rtol) result(acceptable)
        ! Compare scaled original-system residuals without overflowing during
        ! subtraction; tolerate a few ulps when the requested relative
        ! tolerance is smaller than working precision.
        real(dp), intent(in) :: rhs(:) !! right-hand side J*p = rhs
        real(dp), intent(in) :: product(:) !! computed Jacobian-vector product
        real(dp), intent(in) :: rtol !! requested inner relative tolerance

        real(dp) :: scale
        real(dp) :: normalized_rhs_norm
        real(dp) :: normalized_residual_norm

        acceptable = .false.
        if (size(rhs) /= size(product) .or. size(rhs) == 0) return
        if (.not. all(ieee_is_finite(rhs)) .or. &
                .not. all(ieee_is_finite(product))) return
        if (.not. ieee_is_finite(rtol) .or. rtol < 0.0_dp) return
        scale = max(maxval(abs(rhs)), maxval(abs(product)))
        if (scale == 0.0_dp) then
            acceptable = .true.
            return
        end if
        normalized_rhs_norm = norm2(rhs / scale)
        normalized_residual_norm = norm2(rhs / scale - product / scale)
        acceptable = normalized_residual_norm <= &
            max(rtol, 32.0_dp * epsilon(1.0_dp)) * normalized_rhs_norm
    end function inner_residual_acceptable

    subroutine lgmres_jacobian_solve(fun, xbase, fbase, rhs, omega, inner_m, &
            outer_vec, outer_count, outer_k, rtol, solution, nfev, ok, preconditioner, &
            use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, outer_prune_tol, &
            outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, outer_image_basis_reuse, &
            candidate_cache_mib, workspace, image_arrays, restart_size, &
            adaptive_restart, restart_max)
        procedure(nonlin_function) :: fun !! residual callback
        real(dp), intent(in) :: xbase(:) !! nonlinear point at which J is sampled
        real(dp), intent(in) :: fbase(:) !! residual at xbase
        real(dp), intent(in) :: rhs(:) !! linear-system right-hand side
        real(dp), intent(in) :: omega !! finite-difference scale
        integer, intent(in) :: inner_m !! GMRES Krylov dimension
        real(dp), intent(inout) :: outer_vec(:, :) !! retained augmentation vectors
        integer, intent(inout) :: outer_count !! active augmentation-vector count
        integer, intent(in) :: outer_k !! maximum retained augmentation vectors
        real(dp), intent(in) :: rtol !! relative linear tolerance
        real(dp), intent(out) :: solution(:) !! approximate solution of J*x = rhs
        integer, intent(inout) :: nfev !! running residual evaluation count
        logical, intent(out) :: ok !! true for a finite linear solve
        procedure(nonlin_linear_operator), optional :: preconditioner !! approximate J inverse
        logical, intent(in) :: use_central !! true for three-point directional difference
        integer, intent(in) :: jac_step_mode !! global, componentwise, or adaptive steps
        real(dp), intent(in) :: rdiff_use !! relative coordinatewise difference step
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! direct J*v callback
        integer, intent(in) :: orth_mode !! 0 one pass, 1 selective DGKS, 2 two passes
        real(dp), intent(in) :: outer_prune_tol !! 0 disables direction-space filtering
        real(dp), intent(in) :: outer_image_prune_tol !! 0 disables fixed Jacobian-image filtering
        logical, intent(in) :: outer_image_prune_auto !! choose a projection-dependent cutoff
        logical, intent(in) :: outer_image_reuse !! reuse accepted images in Arnoldi
        logical, intent(in) :: outer_candidate_image_reuse !! derive candidate Jv from Arnoldi products
        logical, intent(in) :: outer_fifo_prune_fast !! enable rank-revealing FIFO QR updates
        logical, intent(in) :: outer_image_basis_reuse !! reuse current-Jacobian image QR during insertion
        integer, intent(in) :: candidate_cache_mib !! upper limit on stored raw Jv columns, zero unlimited
        integer, intent(in) :: restart_size !! zero selects original single-cycle LGMRES
        logical, intent(in) :: adaptive_restart !! grow width if restarted LGMRES stalls
        integer, intent(in) :: restart_max !! maximum width allowed when growing

        type(fgmres_workspace), intent(inout), target :: workspace !! persistent Arnoldi workspace
        type(lgmres_image_arrays), intent(inout), target :: image_arrays !! persistent image storage
        real(dp), pointer :: images(:, :)
        real(dp), pointer :: candidate(:)
        real(dp), pointer :: image(:)
        real(dp) :: rhs_norm
        real(dp) :: sol_norm
        real(dp) :: tol_use
        real(dp), allocatable :: cycle_correction(:)
        real(dp), allocatable :: cycle_residual(:)
        real(dp) :: cycle_norm
        real(dp) :: cycle_target
        real(dp) :: previous_norm
        integer :: cycle_cols
        integer :: cycle_width
        integer :: used_cols
        integer :: remaining_cols
        integer :: j
        integer :: kept
        logical :: added
        logical :: basis_cached
        logical :: image_basis_in_place
        logical :: cache_valid
        logical :: cache_candidate
        logical :: cycle_matvec_ok
        logical :: widened_cycle
        integer :: potential_columns

        solution = 0.0_dp
        ok = .false.
        rhs_norm = norm2(rhs)
        if (same_real(rhs_norm, 0.0_dp)) then
            ok = .true.
            return
        end if
        tol_use = min(1.0_dp, max(0.0_dp, rtol))
        if (outer_image_prune_tol > 0.0_dp .or. outer_image_prune_auto) then
            ! Recompute old images at this Newton point. Cached J(x)*v from a
            ! previous step is invalid because the Jacobian may have changed.
            ! Compress old directions before Arnoldi, not just on insertion.
            ! Only old directions need cached images: the prospective
            ! Newton-step image has its own vector, and is never appended to
            ! images at this iterate. An empty history needs no image columns.
            ! Every live image is written by the screening loop below, so
            ! clearing the storage (and touching all its pages) is unnecessary.
            call ensure_lgmres_image_arrays(image_arrays, size(rhs), outer_count)
            images => image_arrays%images
            image => image_arrays%image
            ! Pre-Arnoldi screening and the Arnoldi solve never use their
            ! direction workspaces at the same time. Borrow the first z
            ! column instead of allocating another full-length candidate.
            ! fgmres_jacobian subsequently resets z before using it.
            if (.not. allocated(workspace%z)) then
                if (restart_size > 0) then
                    allocate(workspace%z(size(rhs), min(inner_m, restart_size)))
                else
                    allocate(workspace%z(size(rhs), inner_m))
                end if
            end if
            candidate => workspace%z(:, 1)
            ! Retain a separate basis only for modest systems; borrowed Arnoldi
            ! basis columns are overwritten by the subsequent linear solve.
            ! For large matrix-free systems, preserve the old low-memory path.
            basis_cached = outer_image_basis_reuse .and. outer_count > 0
            ! A separate basis is capped at 8 MiB. The in-place path below
            ! needs no extra n-by-history allocation and has no such cap.
            if (outer_image_reuse) basis_cached = basis_cached .and. &
                int(size(rhs), int64) * int(outer_count, int64) <= 1048576_int64
            ! If Arnoldi is not going to reuse the raw J*v columns, their
            ! existing storage can hold the current-Jacobian QR basis instead.
            ! This avoids retaining a second n-by-history array. With image
            ! reuse enabled, the original J*v values must remain untouched.
            image_basis_in_place = basis_cached .and. .not. outer_image_reuse
            ! Build the current-Jacobian image basis incrementally. Old
            ! filtering reconstructed this basis for every retained vector,
            ! making the compression pass cubic in the history length.
            call reset_lgmres_image_pruning_stream(workspace%pruning)
            kept = 0
            do j = 1, outer_count
                ! kept <= j-1: inserting an accepted column cannot overwrite
                ! an unprocessed direction in the existing retained history.
                candidate = outer_vec(:, j)
                call jacobian_vector(fun, xbase, fbase, candidate, omega, &
                    image, nfev, ok, use_central, jac_step_mode, rdiff_use, &
                    jacobian_product)
                if (.not. ok) return
                if (image_basis_in_place) then
                    call retain_lgmres_image_direction_stream(outer_vec, kept, &
                        candidate, image, outer_image_prune_tol, added, &
                        workspace%pruning, auto_threshold=outer_image_prune_auto, &
                        basis_work=images, stream_limit=outer_count)
                else if (basis_cached) then
                    call retain_lgmres_image_direction_stream(outer_vec, kept, &
                        candidate, image, outer_image_prune_tol, added, &
                        workspace%pruning, auto_threshold=outer_image_prune_auto, &
                        stream_limit=outer_count)
                else if (allocated(workspace%v)) then
                    if (size(workspace%v, 1) == size(rhs) .and. &
                            size(workspace%v, 2) >= outer_count) then
                        ! These Arnoldi columns are dead until the new solve.
                        call retain_lgmres_image_direction_stream(outer_vec, kept, &
                            candidate, image, outer_image_prune_tol, added, &
                            workspace%pruning, auto_threshold=outer_image_prune_auto, &
                            basis_work=workspace%v, stream_limit=outer_count)
                    else
                        call retain_lgmres_image_direction_stream(outer_vec, kept, &
                            candidate, image, outer_image_prune_tol, added, &
                            workspace%pruning, auto_threshold=outer_image_prune_auto, &
                            stream_limit=outer_count)
                    end if
                else
                    call retain_lgmres_image_direction_stream(outer_vec, kept, &
                        candidate, image, outer_image_prune_tol, added, &
                        workspace%pruning, auto_threshold=outer_image_prune_auto, &
                        stream_limit=outer_count)
                end if
                if (added) then
                    ! When images holds the QR basis, raw J*v would overwrite
                    ! its newly orthogonalized column. The raw image cache is
                    ! needed only if Arnoldi explicitly requests image reuse.
                    if (.not. image_basis_in_place) images(:, kept) = image
                    if (outer_image_reuse) then
                        ! The cached image is for the exact original direction.
                        ! Preserve that direction rather than re-normalizing it.
                        outer_vec(:, kept) = candidate
                    end if
                end if
            end do
            outer_count = kept
        end if
        ! Each restarted cycle checks J(x)*accumulated_step directly. The last
        ! check already supplies the candidate image, so caching all raw Arnoldi
        ! products is unnecessary in the restarted variant.
        cache_candidate = outer_candidate_image_reuse .and. restart_size == 0
        if (restart_size > 0 .and. allocated(workspace%raw_products)) &
            deallocate(workspace%raw_products)
        if (cache_candidate .and. candidate_cache_mib > 0) then
            ! The raw Arnoldi products have O(n * inner_maxiter) storage.
            ! An explicit budget may instead pay for one accurate, analytic
            ! J(x)*step call. No Hessenberg reconstruction is used: that
            ! approximation can lose precision for dependent Krylov bases.
            potential_columns = inner_m
            if (.not. outer_image_reuse) potential_columns = inner_m + outer_count
            if (8_int64 * int(size(rhs), int64) * int(potential_columns, int64) > &
                    1048576_int64 * int(candidate_cache_mib, int64)) then
                cache_candidate = .false.
                if (allocated(workspace%raw_products)) &
                    deallocate(workspace%raw_products)
            end if
        end if
        if (restart_size > 0) then
            ! Restart the Arnoldi portion, preserving the retained LGMRES
            ! augmentation vectors throughout the complete inner solve. Only
            ! after the final accumulated correction is formed will the new
            ! direction be screened and appended to the outer history.
            ! inner_m counts ordinary Krylov steps, not retained augmentation
            ! columns; every cycle uses at most restart_size new directions.
            allocate(cycle_correction(size(rhs)), cycle_residual(size(rhs)))
            cycle_residual = rhs
            cycle_norm = rhs_norm
            cycle_target = tol_use * rhs_norm
            remaining_cols = inner_m
            cycle_width = restart_size
            do while (remaining_cols > 0)
                cycle_cols = min(cycle_width, remaining_cols)
                previous_norm = cycle_norm
                if (outer_image_reuse) then
                    call fgmres_jacobian(fun, xbase, fbase, cycle_residual / cycle_norm, &
                        omega, cycle_cols, outer_vec, outer_count, &
                        min(1.0_dp, cycle_target / cycle_norm), cycle_correction, &
                        nfev, ok, preconditioner, use_central, jac_step_mode, &
                        rdiff_use, jacobian_product, orth_mode, workspace, &
                        image_cache=images(:, 1:outer_count), used_cols=used_cols)
                else
                    call fgmres_jacobian(fun, xbase, fbase, cycle_residual / cycle_norm, &
                        omega, cycle_cols, outer_vec, outer_count, &
                        min(1.0_dp, cycle_target / cycle_norm), cycle_correction, &
                        nfev, ok, preconditioner, use_central, jac_step_mode, &
                        rdiff_use, jacobian_product, orth_mode, workspace, &
                        used_cols=used_cols)
                end if
                if (.not. ok .or. used_cols < 1) then
                    ok = .false.
                    return
                end if
                solution = solution + cycle_norm * cycle_correction
                if (.not. all(ieee_is_finite(solution))) then
                    ok = .false.
                    return
                end if
                ! A cycle can stop while evaluating just the retained
                ! augmentation columns; still charge one inner iteration to
                ! guarantee forward progress in the total-budget accounting.
                remaining_cols = remaining_cols - max(1, used_cols - outer_count)
                call jacobian_vector(fun, xbase, fbase, solution, omega, &
                    workspace%v(:, 1), nfev, cycle_matvec_ok, use_central, &
                    jac_step_mode, rdiff_use, jacobian_product)
                if (.not. cycle_matvec_ok) then
                    ok = .false.
                    return
                end if
                cycle_residual = rhs - workspace%v(:, 1)
                if (.not. all(ieee_is_finite(cycle_residual))) then
                    ok = .false.
                    return
                end if
                cycle_norm = safe_two_norm(cycle_residual)
                if (.not. ieee_is_finite(cycle_norm)) then
                    ok = .false.
                    return
                end if
                if (cycle_norm <= cycle_target .or. cycle_norm == 0.0_dp) exit
                widened_cycle = .false.
                ! Opt-in adaptive restart: if the true linear residual fell
                ! by less than 20%, enlarge the next Arnoldi cycle. Only the
                ! current Newton step is accumulated; old augmentation vectors
                ! stay fixed until the complete inner solve returns. Never
                ! exceed the caller's restart width ceiling or iteration budget.
                if (adaptive_restart .and. cycle_norm >= 0.8_dp * previous_norm .and. &
                        cycle_width < restart_max) then
                    cycle_width = cycle_width + min(cycle_width, restart_max - cycle_width)
                    widened_cycle = .true.
                end if
                ! A truly stagnant cycle with no room to enlarge the basis
                ! cannot improve by repeating the identical cycle. Fixed
                ! restart follows exactly the original early-exit policy.
                if (cycle_norm >= previous_norm * (1.0_dp - 16.0_dp * epsilon(1.0_dp))) then
                    if (.not. widened_cycle) exit
                end if
            end do
            ! The verified J*step belongs to the final accumulated Newton
            ! direction and needs no further operator evaluation for pruning.
            if (outer_image_prune_tol > 0.0_dp .or. outer_image_prune_auto) &
                image = workspace%v(:, 1)
        else if (cache_candidate) then
            ! When J*v is supplied analytically, its linear action on the
            ! Arnoldi solution is the same combination of already-computed
            ! J*z columns. The returned image belongs to normalized RHS.
            if (outer_image_reuse) then
                call fgmres_jacobian(fun, xbase, fbase, rhs / rhs_norm, omega, inner_m, &
                    outer_vec, outer_count, tol_use, solution, nfev, ok, preconditioner, &
                    use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, &
                    workspace, image_cache=images(:, 1:outer_count), solution_image=image)
            else
                call fgmres_jacobian(fun, xbase, fbase, rhs / rhs_norm, omega, inner_m, &
                    outer_vec, outer_count, tol_use, solution, nfev, ok, preconditioner, &
                    use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, &
                    workspace, solution_image=image)
            end if
        else if (outer_image_reuse) then
            call fgmres_jacobian(fun, xbase, fbase, rhs / rhs_norm, omega, inner_m, &
                outer_vec, outer_count, tol_use, solution, nfev, ok, preconditioner, &
                use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, &
                workspace, image_cache=images(:, 1:outer_count))
        else
            call fgmres_jacobian(fun, xbase, fbase, rhs / rhs_norm, omega, inner_m, &
                outer_vec, outer_count, tol_use, solution, nfev, ok, preconditioner, &
                use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, workspace)
        end if
        if (.not. ok) return
        if (restart_size == 0) solution = rhs_norm * solution

        if (outer_image_prune_tol > 0.0_dp .or. outer_image_prune_auto) then
            if (.not. cache_candidate .and. restart_size == 0) then
                call jacobian_vector(fun, xbase, fbase, solution, omega, image, &
                    nfev, ok, use_central, jac_step_mode, rdiff_use, jacobian_product)
                if (.not. ok) return
            end if
            ! With candidate reuse, image is J*(solution/rhs_norm) rather
            ! than J*solution. Dependence filtering normalizes each image,
            ! so the positive RHS factor cancels. Avoiding that scaling also
            ! prevents overflow for very large (or subnormal) Newton steps.
            if (basis_cached) then
                if (image_basis_in_place) then
                    call retain_lgmres_image_candidate_cached(outer_vec, outer_count, &
                        solution, image, outer_image_prune_tol, added, &
                        workspace%pruning, cache_valid, auto_threshold=outer_image_prune_auto, &
                        basis_work=images)
                else
                    call retain_lgmres_image_candidate_cached(outer_vec, outer_count, &
                        solution, image, outer_image_prune_tol, added, &
                        workspace%pruning, cache_valid, auto_threshold=outer_image_prune_auto)
                end if
                if (cache_valid) return
            end if
            if (image_basis_in_place) then
                ! Exceptional invalidation of the incremental state: restore
                ! original J*v images before using the legacy reconstruction.
                ! This path is normally unreachable for a completed stream,
                ! but retains a correct fallback rather than treating Q as J*v.
                do j = 1, outer_count
                    call jacobian_vector(fun, xbase, fbase, outer_vec(:, j), omega, &
                        images(:, j), nfev, ok, use_central, jac_step_mode, &
                        rdiff_use, jacobian_product)
                    if (.not. ok) return
                end do
            end if
            if (allocated(workspace%v)) then
                call retain_lgmres_image_direction(outer_vec, outer_count, solution, &
                    images, image, outer_image_prune_tol, added, &
                    auto_threshold=outer_image_prune_auto, &
                    workspace=workspace%pruning, basis_work=workspace%v)
            else
                call retain_lgmres_image_direction(outer_vec, outer_count, solution, &
                    images, image, outer_image_prune_tol, added, &
                    auto_threshold=outer_image_prune_auto, workspace=workspace%pruning)
            end if
            return
        end if
        if (outer_prune_tol > 0.0_dp) then
            ! Before FIFO eviction begins, reuse the solution-space QR span
            ! across nonlinear iterations. Bound the additional basis memory
            ! to eight MiB; large matrix-free problems retain the old path.
            ! Once history is full, a FIFO deletion invalidates incremental QR
            ! and the original rank-revealing reconstruction remains in use.
            if (outer_count < outer_k .and. &
                    int(size(rhs), int64) * int(outer_k, int64) <= 1048576_int64) then
                call retain_lgmres_direction_stream(outer_vec, outer_count, &
                    solution, outer_prune_tol, added, workspace%pruning)
            else if (outer_fifo_prune_fast .and. outer_count == outer_k .and. &
                    int(size(rhs), int64) * int(outer_k, int64) <= 1048576_int64) then
                call retain_lgmres_direction_fifo(outer_vec, outer_count, &
                    solution, outer_prune_tol, added, workspace%pruning)
            else if (allocated(workspace%v)) then
                call retain_lgmres_direction(outer_vec, outer_count, solution, &
                    outer_prune_tol, added, workspace=workspace%pruning, &
                    basis_work=workspace%v)
            else
                call retain_lgmres_direction(outer_vec, outer_count, solution, &
                    outer_prune_tol, added, workspace=workspace%pruning)
            end if
            return
        end if
        sol_norm = norm2(solution)
        if (outer_k > 0 .and. sol_norm > 0.0_dp .and. all(ieee_is_finite(solution))) then
            if (outer_count < outer_k) then
                outer_count = outer_count + 1
            else
                if (outer_k > 1) outer_vec(:, 1:outer_k - 1) = outer_vec(:, 2:outer_k)
            end if
            outer_vec(:, outer_count) = solution / sol_norm
        end if
    end subroutine lgmres_jacobian_solve

    subroutine gmres_jacobian_solve(fun, xbase, fbase, rhs, omega, inner_m, &
            rtol, solution, nfev, ok, preconditioner, use_central, &
            jac_step_mode, rdiff_use, jacobian_product, orth_mode, workspace, restart_size)
        procedure(nonlin_function) :: fun !! residual callback
        real(dp), intent(in) :: xbase(:) !! nonlinear point at which J is sampled
        real(dp), intent(in) :: fbase(:) !! residual at xbase
        real(dp), intent(in) :: rhs(:) !! linear-system right-hand side
        real(dp), intent(in) :: omega !! finite-difference scale
        integer, intent(in) :: inner_m !! GMRES Krylov dimension
        real(dp), intent(in) :: rtol !! relative linear tolerance
        real(dp), intent(out) :: solution(:) !! approximate solution of J*x = rhs
        integer, intent(inout) :: nfev !! running residual evaluation count
        logical, intent(out) :: ok !! true for a finite linear solve
        procedure(nonlin_linear_operator), optional :: preconditioner !! approximate J inverse
        logical, intent(in) :: use_central !! true for three-point directional difference
        integer, intent(in) :: jac_step_mode !! global, componentwise, or adaptive steps
        real(dp), intent(in) :: rdiff_use !! relative coordinatewise difference step
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! direct J*v callback
        integer, intent(in) :: orth_mode !! 0 one pass, 1 selective DGKS, 2 two passes

        type(fgmres_workspace), intent(inout) :: workspace !! persistent Arnoldi workspace
        integer, intent(in) :: restart_size !! zero for original single-cycle GMRES
        real(dp), allocatable :: correction(:)
        real(dp), allocatable :: residual_vector(:)
        real(dp), allocatable :: no_outer(:, :)
        real(dp) :: rhs_norm
        real(dp) :: remaining_norm
        real(dp) :: target_norm
        real(dp) :: tol_use
        integer :: cols_used
        integer :: remaining
        integer :: cycle_cols
        integer :: no_count
        logical :: check_ok

        solution = 0.0_dp
        ok = .false.
        rhs_norm = norm2(rhs)
        if (same_real(rhs_norm, 0.0_dp)) then
            ok = .true.
            return
        end if
        allocate(no_outer(size(rhs), 0))
        no_count = 0
        tol_use = min(1.0_dp, max(0.0_dp, rtol))
        if (restart_size == 0) then
            ! Keep the original SciPy-compatible single cycle (including its
            ! Jacobian evaluation count) when no restart was requested.
            call fgmres_jacobian(fun, xbase, fbase, rhs / rhs_norm, omega, inner_m, &
                no_outer, no_count, tol_use, solution, nfev, ok, preconditioner, &
                use_central, jac_step_mode, rdiff_use, jacobian_product, orth_mode, workspace)
            if (.not. ok) return
            solution = rhs_norm * solution
            call jacobian_vector(fun, xbase, fbase, solution, omega, workspace%v(:, 1), &
                nfev, check_ok, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. check_ok) ok = .false.
            return
        end if

        ! Restarted right-preconditioned GMRES. A cycle keeps at most
        ! restart_size Arnoldi vectors, independent of the total inner_maxiter
        ! budget. Check the unpreconditioned linear residual after every cycle;
        ! this prevents accepting optimistic Hessenberg residual estimates.
        ! Workspaces persist across cycles, so memory is O(n*restart_size).
        allocate(correction(size(rhs)), residual_vector(size(rhs)))
        residual_vector = rhs
        remaining_norm = rhs_norm
        target_norm = tol_use * rhs_norm
        remaining = inner_m
        do while (remaining > 0)
            cycle_cols = min(restart_size, remaining)
            call fgmres_jacobian(fun, xbase, fbase, residual_vector / remaining_norm, &
                omega, cycle_cols, no_outer, no_count, &
                min(1.0_dp, target_norm / remaining_norm), &
                correction, nfev, ok, preconditioner, use_central, &
                jac_step_mode, rdiff_use, jacobian_product, orth_mode, workspace, &
                used_cols=cols_used)
            if (.not. ok .or. cols_used < 1) then
                ok = .false.
                return
            end if
            solution = solution + remaining_norm * correction
            if (.not. all(ieee_is_finite(solution))) then
                ok = .false.
                return
            end if
            remaining = remaining - cols_used
            call jacobian_vector(fun, xbase, fbase, solution, omega, workspace%v(:, 1), &
                nfev, check_ok, use_central, jac_step_mode, rdiff_use, jacobian_product)
            if (.not. check_ok) then
                ok = .false.
                return
            end if
            residual_vector = rhs - workspace%v(:, 1)
            if (.not. all(ieee_is_finite(residual_vector))) then
                ok = .false.
                return
            end if
            remaining_norm = safe_two_norm(residual_vector)
            if (.not. ieee_is_finite(remaining_norm)) then
                ok = .false.
                return
            end if
            if (remaining_norm <= target_norm .or. remaining_norm == 0.0_dp) exit
        end do
    end subroutine gmres_jacobian_solve

    subroutine fgmres_jacobian(fun, xbase, fbase, v0, omega, inner_m, outer_vec, &
            outer_count, atol, solution, nfev, ok, preconditioner, use_central, &
            jac_step_mode, rdiff_use, jacobian_product, orth_mode, workspace, &
            image_cache, solution_image, used_cols)
        procedure(nonlin_function) :: fun !! residual callback
        real(dp), intent(in) :: xbase(:) !! nonlinear point at which J is sampled
        real(dp), intent(in) :: fbase(:) !! residual at xbase
        real(dp), intent(in) :: v0(:) !! normalized right-hand side
        real(dp), intent(in) :: omega !! finite-difference scale
        integer, intent(in) :: inner_m !! ordinary Krylov-vector count
        real(dp), intent(in) :: outer_vec(:, :) !! augmentation vectors
        integer, intent(in) :: outer_count !! active augmentation-vector count
        real(dp), intent(in) :: atol !! normalized inner residual tolerance
        real(dp), intent(out) :: solution(:) !! normalized linear-system solution
        integer, intent(inout) :: nfev !! running residual evaluation count
        logical, intent(out) :: ok !! true for a finite solve
        procedure(nonlin_linear_operator), optional :: preconditioner !! right inverse-Jacobian action
        logical, intent(in) :: use_central !! true for three-point directional difference
        integer, intent(in) :: jac_step_mode !! global, componentwise, or adaptive steps
        real(dp), intent(in) :: rdiff_use !! relative coordinatewise difference step
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! direct J*v callback
        integer, intent(in) :: orth_mode !! 0 one pass, 1 selective DGKS, 2 two passes
        real(dp), intent(in), optional :: image_cache(:, :) !! J(xbase)*old outer vectors, same directions
        real(dp), intent(out), optional :: solution_image(:) !! J(xbase)*normalized Newton step
        integer, intent(out), optional :: used_cols !! actual Arnoldi steps (for restarted GMRES)

        type(fgmres_workspace), intent(inout) :: workspace !! reusable Arnoldi storage
        real(dp) :: alpha
        real(dp) :: breakdown_limit
        real(dp) :: hnorm
        real(dp) :: residual
        real(dp) :: rho
        real(dp) :: tmp
        real(dp) :: w_norm
        integer :: i
        integer :: j
        integer :: k
        integer :: max_cols
        integer :: new_workspace_cols
        integer :: extra_workspace_cols
        integer :: numerical_rank
        integer :: n
        logical :: backsolve_ok
        logical :: breakdown
        logical :: matvec_ok
        logical :: precond_ok

        n = size(v0)
        solution = 0.0_dp
        ok = .false.
        if (present(used_cols)) used_cols = 0
        if (present(solution_image)) then
            solution_image = 0.0_dp
            if (size(solution_image) /= n) return
        end if
        max_cols = inner_m + outer_count
        if (n <= 0 .or. max_cols <= 0) return
        if (present(image_cache)) then
            if (size(image_cache, 1) /= n .or. size(image_cache, 2) /= outer_count) return
            if (.not. all(ieee_is_finite(image_cache))) return
        end if
        ! Keep the separately stored ordinary Krylov directions, but do not
        ! reserve all outer_k extra Arnoldi columns before they are needed.
        ! An oversized outer_k can otherwise consume O(n*outer_k) memory even
        ! when the nonlinear solve uses only one or two retained directions.
        ! Grow geometrically in the number of outer columns, capped at outer_k,
        ! to avoid reallocating at every nonlinear iteration. The old Arnoldi
        ! basis is no longer needed when a new inner solve begins.
        new_workspace_cols = 0
        if (allocated(workspace%v)) then
            if (size(workspace%v, 2) < max_cols + 1) then
                extra_workspace_cols = max(2, 2 * max(0, &
                    size(workspace%v, 2) - 1 - inner_m))
                new_workspace_cols = inner_m + min(size(outer_vec, 2), &
                    max(outer_count, extra_workspace_cols))
                deallocate(workspace%v, workspace%h, workspace%cs, &
                    workspace%sn, workspace%g, workspace%y)
            end if
        else
            new_workspace_cols = inner_m + min(size(outer_vec, 2), &
                max(outer_count, 2))
        end if
        if (new_workspace_cols > 0) then
            allocate(workspace%v(n, new_workspace_cols + 1))
            allocate(workspace%h(new_workspace_cols + 1, new_workspace_cols))
            allocate(workspace%cs(new_workspace_cols), workspace%sn(new_workspace_cols))
            allocate(workspace%g(new_workspace_cols + 1), workspace%y(new_workspace_cols))
        end if
        ! Adaptive restarts can enlarge the Arnoldi subspace between cycles.
        ! A pre-existing z buffer from screening or a smaller cycle must grow
        ! before its new live columns are associated below.
        if (allocated(workspace%z)) then
            if (size(workspace%z, 1) /= n .or. &
                    size(workspace%z, 2) < inner_m) deallocate(workspace%z)
        end if
        if (.not. allocated(workspace%z)) allocate(workspace%z(n, inner_m))
        if (present(solution_image)) then
            if (allocated(workspace%raw_products)) then
                if (present(image_cache)) then
                    if (size(workspace%raw_products, 2) < inner_m) &
                        deallocate(workspace%raw_products)
                else
                    if (size(workspace%raw_products, 2) < max_cols) &
                        deallocate(workspace%raw_products)
                end if
            end if
            if (.not. allocated(workspace%raw_products)) then
                if (present(image_cache)) then
                    ! The retained images are stored by the caller.
                    allocate(workspace%raw_products(n, inner_m))
                else
                    allocate(workspace%raw_products(n, max_cols))
                end if
            end if
        end if
        ! Work only on live columns. The next unused Arnoldi basis column
        ! doubles as the matvec and orthogonalization scratch vector, so the
        ! previous n-element w allocation is no longer needed.
        associate (v => workspace%v(:, 1:max_cols + 1), &
                z => workspace%z(:, 1:inner_m), &
                h => workspace%h(1:max_cols + 1, 1:max_cols), &
                cs => workspace%cs(1:max_cols), &
                sn => workspace%sn(1:max_cols), &
                g => workspace%g(1:max_cols + 1), &
                y => workspace%y(1:max_cols))
            v = 0.0_dp
            z = 0.0_dp
            h = 0.0_dp
            cs = 0.0_dp
            sn = 0.0_dp
            g = 0.0_dp
            y = 0.0_dp
            v(:, 1) = v0
            g(1) = 1.0_dp
            k = 0

            do j = 1, max_cols
                associate (w => v(:, j + 1))
                    if (j <= outer_count) then
                        ! Retained solution directions already live in outer_vec.
                        ! Do not duplicate them in the ordinary Krylov basis z.
                        if (present(image_cache)) then
                            w = image_cache(:, j)
                        else
                            call jacobian_vector(fun, xbase, fbase, outer_vec(:, j), &
                                omega, w, nfev, matvec_ok, use_central, &
                                jac_step_mode, rdiff_use, jacobian_product)
                            if (.not. matvec_ok) return
                        end if
                    else
                        if (present(preconditioner)) then
                            ! Write M*v directly into its persistent Krylov
                            ! column. Reading from v (not z) avoids copying a
                            ! full n-vector before and after preconditioning.
                            if (j == outer_count + 1) then
                                call preconditioner(v(:, 1), z(:, j - outer_count), precond_ok)
                            else
                                call preconditioner(v(:, j), z(:, j - outer_count), precond_ok)
                            end if
                            if (.not. precond_ok) return
                            if (.not. all(ieee_is_finite(z(:, j - outer_count)))) return
                        else if (j == outer_count + 1) then
                            z(:, j - outer_count) = v(:, 1)
                        else
                            z(:, j - outer_count) = v(:, j)
                        end if
                        call jacobian_vector(fun, xbase, fbase, z(:, j - outer_count), &
                            omega, w, nfev, matvec_ok, use_central, &
                            jac_step_mode, rdiff_use, jacobian_product)
                        if (.not. matvec_ok) return
                    end if
                    ! Retain unmodified operator products before Gram-Schmidt changes
                    ! w. This is only allocated by the explicit candidate-reuse mode.
                    if (present(solution_image)) then
                        if (present(image_cache)) then
                            if (j > outer_count) workspace%raw_products(:, j - outer_count) = w
                        else
                            workspace%raw_products(:, j) = w
                        end if
                    end if
                    w_norm = safe_two_norm(w)
                    if (.not. ieee_is_finite(w_norm)) return
                    do i = 1, j
                        h(i, j) = dot_product(v(:, i), w)
                        w = w - h(i, j) * v(:, i)
                    end do
                    h(j + 1, j) = safe_two_norm(w)
                    ! DGKS criterion: a large norm reduction signals loss of
                    ! orthogonality in a one-pass modified Gram-Schmidt projection.
                    ! Accumulate second-pass coefficients into the same Hessenberg
                    ! column; no additional matrix-vector products are required.
                    if (orth_mode == 2 .or. &
                            (orth_mode == 1 .and. h(j + 1, j) <= 0.5_dp * w_norm)) then
                        do i = 1, j
                            tmp = dot_product(v(:, i), w)
                            h(i, j) = h(i, j) + tmp
                            w = w - tmp * v(:, i)
                        end do
                        h(j + 1, j) = safe_two_norm(w)
                    end if
                    if (.not. ieee_is_finite(h(j + 1, j))) return
                    if (.not. all(ieee_is_finite(h(1:j, j)))) return
                    breakdown_limit = epsilon(1.0_dp) * w_norm
                    breakdown = .not. (h(j + 1, j) > breakdown_limit)
                    if (.not. breakdown) then
                        v(:, j + 1) = w / h(j + 1, j)
                    else
                        v(:, j + 1) = 0.0_dp
                    end if

                    do i = 1, j - 1
                        tmp = cs(i) * h(i, j) + sn(i) * h(i + 1, j)
                        h(i + 1, j) = -sn(i) * h(i, j) + cs(i) * h(i + 1, j)
                        h(i, j) = tmp
                    end do
                    rho = hypot(h(j, j), h(j + 1, j))
                    if (same_real(rho, 0.0_dp) .or. .not. ieee_is_finite(rho)) then
                        cs(j) = 1.0_dp
                        sn(j) = 0.0_dp
                    else
                        cs(j) = h(j, j) / rho
                        sn(j) = h(j + 1, j) / rho
                    end if
                    h(j, j) = cs(j) * h(j, j) + sn(j) * h(j + 1, j)
                    h(j + 1, j) = 0.0_dp
                    tmp = cs(j) * g(j) + sn(j) * g(j + 1)
                    g(j + 1) = -sn(j) * g(j) + cs(j) * g(j + 1)
                    g(j) = tmp
                    residual = abs(g(j + 1))
                    k = j
                end associate
                if (residual < atol .or. breakdown) exit
            end do

            if (k <= 0) return
            backsolve_ok = .true.
            do i = k, 1, -1
                tmp = g(i)
                if (i < k) tmp = tmp - dot_product(h(i, i + 1:k), y(i + 1:k))
                hnorm = abs(h(i, i))
                if (.not. ieee_is_finite(hnorm) .or. &
                        hnorm <= sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(h)))) then
                    backsolve_ok = .false.
                    exit
                end if
                y(i) = tmp / h(i, i)
            end do

            if (.not. backsolve_ok) then
                ! The triangular Arnoldi system can be rank deficient when an
                ! LGMRES augmentation vector is nearly dependent on earlier
                ! directions. Solve the small least-squares problem by a direct
                ! one-sided Jacobi SVD, not normal equations (which square the
                ! condition number and lose small singular directions).
                call krylov_rank_least_squares(h(1:k, 1:k), g(1:k), &
                    y(1:k), numerical_rank, backsolve_ok)
                if (.not. backsolve_ok) return
            end if

            solution = 0.0_dp
            do i = 1, k
                alpha = y(i)
                if (i <= outer_count) then
                    solution = solution + alpha * outer_vec(:, i)
                else
                    solution = solution + alpha * z(:, i - outer_count)
                end if
            end do
            ok = all(ieee_is_finite(solution))
            if (present(used_cols)) used_cols = k
            if (present(solution_image)) then
                if (present(image_cache)) then
                    ! The retained images are already held by the caller. Form
                    ! the same candidate combination without a second copy of
                    ! those columns. Never reuse products across Newton steps.
                    solution_image = 0.0_dp
                    do i = 1, k
                        if (i <= outer_count) then
                            solution_image = solution_image + y(i) * image_cache(:, i)
                        else
                            solution_image = solution_image + y(i) * workspace%raw_products(:, i - outer_count)
                        end if
                    end do
                else
                    solution_image = matmul(workspace%raw_products(:, 1:k), y(1:k))
                end if
                ok = ok .and. all(ieee_is_finite(solution_image))
            end if
        end associate
    end subroutine fgmres_jacobian

    subroutine jacobian_vector(fun, xbase, fbase, vec, omega, product, nfev, ok, &
            use_central, jac_step_mode, rdiff_use, jacobian_product)
        procedure(nonlin_function) :: fun !! residual callback
        real(dp), intent(in) :: xbase(:) !! nonlinear point at which J is sampled
        real(dp), intent(in) :: fbase(:) !! residual at xbase
        real(dp), intent(in) :: vec(:) !! vector multiplied by the numerical Jacobian
        real(dp), intent(in) :: omega !! finite-difference scale
        real(dp), intent(out) :: product(:) !! numerical J*vec product
        integer, intent(inout) :: nfev !! running residual evaluation count
        logical, intent(out) :: ok !! true for a finite product
        logical, intent(in) :: use_central !! true for three-point directional difference
        integer, intent(in) :: jac_step_mode !! global, componentwise, or adaptive steps
        real(dp), intent(in) :: rdiff_use !! relative coordinatewise difference step
        procedure(nonlin_jacobian_product), optional :: jacobian_product !! direct J*v callback

        real(dp), allocatable :: fminus(:)
        real(dp), allocatable :: ftrial(:)
        real(dp), allocatable :: xtrial(:)
        real(dp) :: scale
        real(dp) :: vec_norm
        logical :: use_coordinate_steps
        integer :: j

        product = 0.0_dp
        ok = .false.
        if (present(jacobian_product)) then
            if (.not. all(ieee_is_finite(vec))) return
            if (all(vec == 0.0_dp)) then
                ok = .true.
                return
            end if
            call jacobian_product(xbase, vec, product, ok)
            if (ok) ok = all(ieee_is_finite(product))
            ! An analytic J*v call is not a residual evaluation.
            return
        end if
        ! Forced coordinatewise mode does not depend on the global omega.
        if (jac_step_mode == jac_step_componentwise) then
            call coordinate_jacobian_vector(fun, xbase, fbase, vec, rdiff_use, &
                product, nfev, ok, use_central)
            return
        end if
        vec_norm = norm2(vec)
        if (same_real(vec_norm, 0.0_dp)) then
            ok = .true.
            return
        end if
        scale = omega / vec_norm
        if (.not. ieee_is_finite(scale) .or. same_real(scale, 0.0_dp)) then
            if (jac_step_mode == jac_step_adaptive) then
                call coordinate_jacobian_vector(fun, xbase, fbase, vec, rdiff_use, &
                    product, nfev, ok, use_central)
            end if
            return
        end if
        use_coordinate_steps = .false.
        if (jac_step_mode == jac_step_adaptive) then
            ! In an ill-scaled state the global x norm can dictate an
            ! inappropriate step for an unrelated (small) active coordinate.
            ! Also switch when a requested directional perturbation vanishes
            ! on addition due to floating-point rounding.
            do j = 1, size(vec)
                if (vec(j) == 0.0_dp) cycle
                if (max(1.0_dp, abs(xbase(j))) < &
                        1.0e-5_dp * max(1.0_dp, maxval(abs(xbase)))) then
                    use_coordinate_steps = .true.
                    exit
                end if
                ! A derivative probe may be representable in x but far
                ! too small compared with the local rdiff scale because
                ! the global step divides by a very large residual norm.
                if (abs(scale) > huge(1.0_dp) / max(1.0_dp, abs(vec(j)))) then
                    use_coordinate_steps = .true.
                    exit
                end if
                if (abs(scale * vec(j)) / max(1.0_dp, abs(xbase(j))) < &
                        1.0e-3_dp * rdiff_use) then
                    use_coordinate_steps = .true.
                    exit
                end if
                if (xbase(j) + scale * vec(j) == xbase(j)) then
                    use_coordinate_steps = .true.
                    exit
                end if
            end do
        end if
        if (use_coordinate_steps) then
            call coordinate_jacobian_vector(fun, xbase, fbase, vec, rdiff_use, &
                product, nfev, ok, use_central)
            return
        end if
        allocate(xtrial(size(xbase)), ftrial(size(fbase)))
        xtrial = xbase + scale * vec
        if (use_central) then
            ! If the requested step rounds away, a zero directional
            ! difference is not evidence that the Jacobian action is zero.
            if (.not. all(ieee_is_finite(xtrial))) return
            if (all(xtrial == xbase)) return
        end if
        call fun(xtrial, ftrial)
        nfev = nfev + 1
        if (.not. all(ieee_is_finite(ftrial))) return
        if (use_central) then
            allocate(fminus(size(fbase)))
            xtrial = xbase - scale * vec
            if (.not. all(ieee_is_finite(xtrial))) return
            if (all(xtrial == xbase)) return
            call fun(xtrial, fminus)
            nfev = nfev + 1
            if (.not. all(ieee_is_finite(fminus))) return
            ! Halve before subtracting to avoid overflow when finite
            ! opposite-signed residual samples have a large difference.
            product = (0.5_dp * ftrial - 0.5_dp * fminus) / scale
        else
            product = (ftrial - fbase) / scale
        end if
        ok = all(ieee_is_finite(product))
    end subroutine jacobian_vector

    subroutine coordinate_jacobian_vector(fun, xbase, fbase, vec, rdiff, &
            product, nfev, ok, use_central)
        ! Matrix-free but coordinate-separated J*v: each nonzero component of
        ! v uses a perturbation sized for that coordinate, rather than for
        ! the largest unrelated state component. This costs O(nnz(v)) calls.
        procedure(nonlin_function) :: fun !! square residual callback
        real(dp), intent(in) :: xbase(:) !! point at which the Jacobian is sampled
        real(dp), intent(in) :: fbase(:) !! already evaluated residual
        real(dp), intent(in) :: vec(:) !! direction to multiply by the Jacobian
        real(dp), intent(in) :: rdiff !! base relative finite-difference step
        real(dp), intent(out) :: product(:) !! approximate J*vec
        integer, intent(inout) :: nfev !! residual function evaluation count
        logical, intent(out) :: ok !! true if every component was evaluated
        logical, intent(in) :: use_central !! centered rather than forward differences

        real(dp) :: h
        real(dp) :: coordinate_scale
        real(dp), allocatable :: fplus(:)
        real(dp), allocatable :: fminus(:)
        real(dp), allocatable :: xtrial(:)
        integer :: j

        product = 0.0_dp
        ok = .false.
        allocate(xtrial(size(xbase)), fplus(size(fbase)))
        if (use_central) allocate(fminus(size(fbase)))
        do j = 1, size(vec)
            if (vec(j) == 0.0_dp) cycle
            coordinate_scale = max(1.0_dp, abs(xbase(j)))
            ! Avoid overflow even when an excessive user rdiff is passed.
            if (rdiff > (huge(1.0_dp) / coordinate_scale)) return
            ! A very small user step must not silently return a zero
            ! derivative when it rounds to the current coordinate.
            h = max(rdiff, 8.0_dp * epsilon(1.0_dp)) * coordinate_scale
            if (h > huge(1.0_dp) - abs(xbase(j))) return
            xtrial = xbase
            xtrial(j) = xbase(j) + h
            if (xtrial(j) == xbase(j) .or. .not. ieee_is_finite(xtrial(j))) return
            call fun(xtrial, fplus)
            nfev = nfev + 1
            if (.not. all(ieee_is_finite(fplus))) return
            if (use_central) then
                xtrial(j) = xbase(j) - h
                if (xtrial(j) == xbase(j) .or. .not. ieee_is_finite(xtrial(j))) return
                call fun(xtrial, fminus)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(fminus))) return
                product = product + vec(j) * ((0.5_dp * fplus - 0.5_dp * fminus) / h)
            else
                product = product + vec(j) * ((fplus - fbase) / h)
            end if
            if (.not. all(ieee_is_finite(product))) return
        end do
        ok = .true.
    end subroutine coordinate_jacobian_vector


    subroutine setup_tolerances(tol, ftol, fatol, xtol, xatol, f_abs_tol, &
            f_rel_tol, x_abs_tol, x_rel_tol, f_abs_active, f_rel_active, &
            x_abs_active, x_rel_active, ok)
        real(dp), intent(in), optional :: tol !! generic root tolerance
        real(dp), intent(in), optional :: ftol !! relative residual tolerance
        real(dp), intent(in), optional :: fatol !! absolute residual tolerance
        real(dp), intent(in), optional :: xtol !! relative step tolerance
        real(dp), intent(in), optional :: xatol !! absolute step tolerance
        real(dp), intent(out) :: f_abs_tol !! resolved absolute residual tolerance
        real(dp), intent(out) :: f_rel_tol !! resolved relative residual tolerance
        real(dp), intent(out) :: x_abs_tol !! resolved absolute step tolerance
        real(dp), intent(out) :: x_rel_tol !! resolved relative step tolerance
        logical, intent(out) :: f_abs_active !! whether absolute residual test is active
        logical, intent(out) :: f_rel_active !! whether relative residual test is active
        logical, intent(out) :: x_abs_active !! whether absolute step test is active
        logical, intent(out) :: x_rel_active !! whether relative step test is active
        logical, intent(out) :: ok !! true for valid tolerance inputs

        f_abs_tol = epsilon(1.0_dp) ** (1.0_dp / 3.0_dp)
        f_rel_tol = 0.0_dp
        x_abs_tol = 0.0_dp
        x_rel_tol = 0.0_dp
        f_abs_active = .true.
        f_rel_active = .false.
        x_abs_active = .false.
        x_rel_active = .false.
        ok = .true.

        if (present(tol)) then
            if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) then
                ok = .false.
                return
            end if
            f_abs_active = .false.
            f_rel_active = .false.
            x_abs_active = .false.
            x_rel_active = .true.
            x_rel_tol = tol
        end if

        if (present(fatol)) then
            if (.not. ieee_is_finite(fatol) .or. fatol < 0.0_dp) then
                ok = .false.
                return
            end if
            f_abs_tol = fatol
            f_abs_active = .true.
        end if
        if (present(ftol)) then
            if (.not. ieee_is_finite(ftol) .or. ftol < 0.0_dp) then
                ok = .false.
                return
            end if
            f_rel_tol = ftol
            f_rel_active = .true.
        end if
        if (present(xatol)) then
            if (.not. ieee_is_finite(xatol) .or. xatol < 0.0_dp) then
                ok = .false.
                return
            end if
            x_abs_tol = xatol
            x_abs_active = .true.
        end if
        if (present(xtol)) then
            if (.not. ieee_is_finite(xtol) .or. xtol < 0.0_dp) then
                ok = .false.
                return
            end if
            x_rel_tol = xtol
            x_rel_active = .true.
        end if
    end subroutine setup_tolerances

    logical function termination_met(fvec, x, dx, f0_norm, f_abs_tol, f_rel_tol, &
            x_abs_tol, x_rel_tol, f_abs_active, f_rel_active, x_abs_active, &
            x_rel_active, iterations, iter_count, status, tol_norm) result(done)
        real(dp), intent(in) :: fvec(:) !! current residual vector
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: dx(:) !! previous accepted step
        real(dp), intent(in) :: f0_norm !! initial residual norm under selected metric
        real(dp), intent(in) :: f_abs_tol !! absolute residual tolerance
        real(dp), intent(in) :: f_rel_tol !! relative residual tolerance
        real(dp), intent(in) :: x_abs_tol !! absolute step tolerance
        real(dp), intent(in) :: x_rel_tol !! relative step tolerance
        logical, intent(in) :: f_abs_active !! enable absolute residual test
        logical, intent(in) :: f_rel_active !! enable relative residual test
        logical, intent(in) :: x_abs_active !! enable absolute step test
        logical, intent(in) :: x_rel_active !! enable relative step test
        integer, intent(in), optional :: iterations !! exact iteration count request
        integer, intent(in) :: iter_count !! current one-based check count
        integer, intent(out) :: status !! status when termination occurs
        procedure(nonlin_norm_function), optional :: tol_norm !! custom residual/step norm

        real(dp) :: dx_norm
        real(dp) :: f_norm
        real(dp) :: x_norm
        logical :: f_ok
        logical :: x_ok

        done = .false.
        status = 0
        f_norm = maxval(abs(fvec))
        x_norm = maxval(abs(x))
        dx_norm = maxval(abs(dx))
        if (present(tol_norm)) then
            f_norm = tol_norm(fvec)
            x_norm = tol_norm(x)
            ! No actual step exists before the first iteration. Do not pass
            ! the huge sentinel to an arbitrary user norm (which could
            ! overflow under IEEE traps), because it cannot affect stopping.
            if (iter_count > 1) dx_norm = tol_norm(dx)
            if (.not. ieee_is_finite(f_norm) .or. &
                    .not. ieee_is_finite(x_norm) .or. &
                    .not. ieee_is_finite(dx_norm)) then
                status = -1
                done = .true.
                return
            end if
            if (f_norm < 0.0_dp .or. x_norm < 0.0_dp .or. dx_norm < 0.0_dp) then
                status = -1
                done = .true.
                return
            end if
        end if

        if (same_real(f_norm, 0.0_dp)) then
            status = 1
            done = .true.
            return
        end if

        if (present(iterations)) then
            if (iter_count > iterations) then
                status = 2
                done = .true.
            end if
            return
        end if

        ! SciPy tests the initial point using an infinite sentinel step. Unless
        ! the residual is exactly zero, this initial check cannot converge.
        ! Explicitly skip tolerance tests here to avoid Inf/Inf and overflow
        ! with user-defined norms, while keeping iteration accounting intact.
        if (iter_count == 1) return

        f_ok = .true.
        if (f_abs_active) f_ok = f_ok .and. f_norm <= f_abs_tol
        if (f_rel_active) f_ok = f_ok .and. f_norm <= f_rel_tol * f0_norm
        x_ok = .true.
        if (x_abs_active) x_ok = x_ok .and. dx_norm <= x_abs_tol
        if (x_rel_active) x_ok = x_ok .and. dx_norm <= x_rel_tol * x_norm
        if (f_ok .and. x_ok) then
            status = 1
            done = .true.
        end if
    end function termination_met

    subroutine armijo_step(fun, x, fvec, dx, xnew, fnew, nfev, status)
        procedure(nonlin_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: fvec(:) !! residual at x
        real(dp), intent(in) :: dx(:) !! proposed full step
        real(dp), intent(out) :: xnew(:) !! accepted iterate
        real(dp), intent(out) :: fnew(:) !! residual at xnew
        integer, intent(inout) :: nfev !! running residual evaluation count
        integer, intent(out) :: status !! zero on finite accepted residual

        real(dp), allocatable :: cache_f(:)
        real(dp) :: a
        real(dp) :: alpha0
        real(dp) :: alpha1
        real(dp) :: alpha2
        real(dp) :: b
        real(dp) :: cache_phi
        real(dp) :: cache_s
        real(dp), parameter :: c1 = 1.0e-4_dp
        real(dp) :: derphi0
        real(dp) :: discr
        real(dp) :: factor
        real(dp), parameter :: amin = 1.0e-2_dp
        real(dp) :: phi0
        real(dp) :: phi_a0
        real(dp) :: phi_a1
        real(dp) :: phi_a2
        real(dp) :: s
        logical :: accepted

        allocate(cache_f(size(fvec)))
        cache_s = 0.0_dp
        cache_f = fvec
        phi0 = safe_squared_norm(fvec)
        cache_phi = phi0
        derphi0 = -phi0
        alpha0 = 1.0_dp
        accepted = .false.

        call eval_phi(alpha0, phi_a0)
        if (phi_a0 <= phi0 + c1 * alpha0 * derphi0) then
            s = alpha0
            accepted = .true.
        else
            denom: block
                real(dp) :: qden
                qden = 2.0_dp * (phi_a0 - phi0 - derphi0 * alpha0)
                if (same_real(qden, 0.0_dp)) then
                    alpha1 = 0.5_dp * alpha0
                else
                    alpha1 = -derphi0 * alpha0 * alpha0 / qden
                end if
            end block denom
            call eval_phi(alpha1, phi_a1)
            if (phi_a1 <= phi0 + c1 * alpha1 * derphi0) then
                s = alpha1
                accepted = .true.
            else
                do while (alpha1 > amin)
                    factor = alpha0 * alpha0 * alpha1 * alpha1 * (alpha1 - alpha0)
                    if (same_real(factor, 0.0_dp)) exit
                    a = alpha0 * alpha0 * (phi_a1 - phi0 - derphi0 * alpha1) - &
                        alpha1 * alpha1 * (phi_a0 - phi0 - derphi0 * alpha0)
                    a = a / factor
                    b = -alpha0**3 * (phi_a1 - phi0 - derphi0 * alpha1) + &
                        alpha1**3 * (phi_a0 - phi0 - derphi0 * alpha0)
                    b = b / factor
                    if (same_real(a, 0.0_dp)) then
                        alpha2 = 0.5_dp * alpha1
                    else
                        discr = abs(b * b - 3.0_dp * a * derphi0)
                        alpha2 = (-b + sqrt(discr)) / (3.0_dp * a)
                    end if
                    call eval_phi(alpha2, phi_a2)
                    if (phi_a2 <= phi0 + c1 * alpha2 * derphi0) then
                        s = alpha2
                        accepted = .true.
                        exit
                    end if
                    if ((alpha1 - alpha2) > alpha1 / 2.0_dp .or. &
                            (1.0_dp - alpha2 / alpha1) < 0.96_dp) then
                        alpha2 = alpha1 / 2.0_dp
                    end if
                    alpha0 = alpha1
                    alpha1 = alpha2
                    phi_a0 = phi_a1
                    phi_a1 = phi_a2
                end do
            end if
        end if

        if (.not. accepted) s = 1.0_dp
        xnew = x + s * dx
        if (same_real(s, cache_s)) then
            fnew = cache_f
        else
            call fun(xnew, fnew)
            nfev = nfev + 1
        end if
        if (all(ieee_is_finite(fnew))) then
            status = 0
        else
            status = 1
        end if

    contains

        subroutine eval_phi(step, value)
            real(dp), intent(in) :: step !! trial line-search step length
            real(dp), intent(out) :: value !! squared residual 2-norm at trial point

            real(dp), allocatable :: trial_f(:)
            real(dp), allocatable :: trial_x(:)

            if (same_real(step, cache_s)) then
                value = cache_phi
                return
            end if
            allocate(trial_x(size(x)), trial_f(size(fvec)))
            trial_x = x + step * dx
            call fun(trial_x, trial_f)
            nfev = nfev + 1
            value = safe_squared_norm(trial_f)
            cache_s = step
            cache_phi = value
            cache_f = trial_f
        end subroutine eval_phi

    end subroutine armijo_step

    subroutine wolfe_step(fun, x, fvec, dx, xnew, fnew, nfev, status)
        procedure(nonlin_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: fvec(:) !! residual at x
        real(dp), intent(in) :: dx(:) !! proposed full step
        real(dp), intent(out) :: xnew(:) !! accepted iterate
        real(dp), intent(out) :: fnew(:) !! residual at xnew
        integer, intent(inout) :: nfev !! running residual evaluation count
        integer, intent(out) :: status !! zero on finite accepted residual

        real(dp), allocatable :: cache_f(:)
        real(dp) :: cache_phi
        real(dp) :: cache_s
        real(dp) :: phi0
        real(dp) :: phi1
        real(dp) :: s
        real(dp) :: s_norm
        logical :: found

        allocate(cache_f(size(fvec)))
        cache_s = 0.0_dp
        cache_f = fvec
        phi0 = safe_squared_norm(fvec)
        cache_phi = phi0
        s_norm = norm2(x) / norm2(dx)

        call scalar_search_wolfe1(phi, derphi, phi0, s, phi1, found, &
            amin=1.0e-2_dp, xtol=1.0e-2_dp)
        if (.not. found) s = 1.0_dp

        xnew = x + s * dx
        if (same_real(s, cache_s)) then
            fnew = cache_f
        else
            call fun(xnew, fnew)
            nfev = nfev + 1
        end if
        if (all(ieee_is_finite(fnew))) then
            status = 0
        else
            status = 1
        end if

    contains

        function phi(step) result(value)
            real(dp), intent(in) :: step !! trial line-search step
            real(dp) :: value !! squared residual norm

            call eval_phi(step, .true., value)
        end function phi

        function derphi(step) result(value)
            real(dp), intent(in) :: step !! trial line-search step
            real(dp) :: value !! finite-difference directional derivative

            real(dp) :: ds
            real(dp) :: phi_base
            real(dp) :: phi_plus

            ds = (abs(step) + s_norm + 1.0_dp) * 1.0e-8_dp
            call eval_phi(step + ds, .false., phi_plus)
            call eval_phi(step, .true., phi_base)
            value = (phi_plus - phi_base) / ds
        end function derphi

        subroutine eval_phi(step, store, value)
            real(dp), intent(in) :: step !! trial line-search step
            logical, intent(in) :: store !! whether to update the one-entry cache
            real(dp), intent(out) :: value !! squared residual norm

            real(dp), allocatable :: trial_f(:)
            real(dp), allocatable :: trial_x(:)

            if (same_real(step, cache_s)) then
                value = cache_phi
                return
            end if
            allocate(trial_x(size(x)), trial_f(size(fvec)))
            trial_x = x + step * dx
            call fun(trial_x, trial_f)
            nfev = nfev + 1
            value = safe_squared_norm(trial_f)
            if (store) then
                cache_s = step
                cache_phi = value
                cache_f = trial_f
            end if
        end subroutine eval_phi

    end subroutine wolfe_step

    pure function safe_squared_norm(v) result(value)
        real(dp), intent(in) :: v(:) !! vector whose squared Euclidean norm is needed
        real(dp) :: value

        if (.not. all(ieee_is_finite(v))) then
            value = huge(1.0_dp)
        else
            value = dot_product(v, v)
            if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
        end if
    end function safe_squared_norm

    subroutine append_broyden_update(cs, ds, count, c, d, collapsed)
        real(dp), intent(inout) :: cs(:, :) !! low-rank left factors
        real(dp), intent(inout) :: ds(:, :) !! low-rank right factors
        integer, intent(inout) :: count !! retained rank-one updates
        real(dp), intent(in) :: c(:) !! new left factor
        real(dp), intent(in) :: d(:) !! new right factor
        logical, intent(inout) :: collapsed !! true after rank exceeds dimension

        count = count + 1
        if (count <= size(cs, 2)) then
            cs(:, count) = c
            ds(:, count) = d
        end if
        if (count > size(c)) collapsed = .true.
    end subroutine append_broyden_update

    subroutine reduce_broyden_matrix(h, base_scale, cs, ds, count, max_rank, &
            reduction, ok)
        real(dp), intent(inout) :: h(:, :) !! dense inverse-Jacobian approximation
        real(dp), intent(in) :: base_scale !! scalar identity base
        real(dp), intent(inout) :: cs(:, :) !! low-rank left factors
        real(dp), intent(inout) :: ds(:, :) !! low-rank right factors
        integer, intent(inout) :: count !! retained rank-one update count
        integer, intent(in) :: max_rank !! SciPy-compatible maximum rank setting
        character(len=*), intent(in) :: reduction !! restart, simple, or svd
        logical, intent(out) :: ok !! true when reduction succeeds

        real(dp), allocatable :: ata(:, :)
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: update(:, :)
        real(dp), allocatable :: new_c(:, :)
        real(dp), allocatable :: new_d(:, :)
        integer :: i
        integer :: j
        integer :: keep
        integer :: linalg_status
        integer :: n
        integer :: p
        integer :: target

        ok = .true.
        n = size(h, 1)
        target = max_rank - 1
        if (target < 1 .or. count <= 0) return

        select case (trim(reduction))
        case ('restart')
            if (count > target) then
                count = 0
                call set_scaled_identity(h, base_scale)
            end if
        case ('simple')
            if (count > target) then
                keep = target
                cs(:, 1:keep) = cs(:, count - keep + 1:count)
                ds(:, 1:keep) = ds(:, count - keep + 1:count)
                count = keep
                call rebuild_broyden_matrix(h, base_scale, cs, ds, count)
            end if
        case ('svd')
            p = min(target, n)
            keep = max(0, min(target - 2, p - 1))
            if (count < p) return
            if (keep == 0) then
                count = 0
                call set_scaled_identity(h, base_scale)
                return
            end if

            allocate(update(n, n), ata(n, n), eig(n), vec(n, n))
            update = h
            do i = 1, n
                update(i, i) = update(i, i) - base_scale
            end do
            ata = matmul(transpose(update), update)
            call symmetric_eigen_jacobi(ata, eig, vec, linalg_status)
            if (linalg_status /= linalg_status_success) then
                ok = .false.
                return
            end if
            allocate(new_c(n, keep), new_d(n, keep))
            do j = 1, keep
                i = n - keep + j
                new_d(:, j) = vec(:, i)
                new_c(:, j) = matmul(update, new_d(:, j))
            end do
            cs(:, 1:keep) = new_c
            ds(:, 1:keep) = new_d
            count = keep
            call rebuild_broyden_matrix(h, base_scale, cs, ds, count)
        case default
            ok = .false.
        end select
    end subroutine reduce_broyden_matrix

    subroutine rebuild_broyden_matrix(h, base_scale, cs, ds, count)
        real(dp), intent(out) :: h(:, :) !! reconstructed dense inverse Jacobian
        real(dp), intent(in) :: base_scale !! scalar identity base
        real(dp), intent(in) :: cs(:, :) !! low-rank left factors
        real(dp), intent(in) :: ds(:, :) !! low-rank right factors
        integer, intent(in) :: count !! number of active factors

        integer :: i
        integer :: j

        call set_scaled_identity(h, base_scale)
        do i = 1, count
            do j = 1, size(h, 2)
                h(:, j) = h(:, j) + cs(:, i) * ds(j, i)
            end do
        end do
    end subroutine rebuild_broyden_matrix

    subroutine append_history(hist_dx, hist_df, count, dx, df)
        real(dp), intent(inout) :: hist_dx(:, :) !! retained accepted steps by column
        real(dp), intent(inout) :: hist_df(:, :) !! retained residual differences by column
        integer, intent(inout) :: count !! current retained history count
        real(dp), intent(in) :: dx(:) !! newest accepted step
        real(dp), intent(in) :: df(:) !! newest residual difference

        integer :: m

        m = size(hist_dx, 2)
        if (m <= 0) return
        if (count < m) then
            count = count + 1
        else
            hist_dx(:, 1:m - 1) = hist_dx(:, 2:m)
            hist_df(:, 1:m - 1) = hist_df(:, 2:m)
        end if
        hist_dx(:, count) = dx
        hist_df(:, count) = df
    end subroutine append_history

    subroutine build_anderson_matrix(hist_df, count, w0, a)
        real(dp), intent(in) :: hist_df(:, :) !! retained residual differences
        integer, intent(in) :: count !! active history columns
        real(dp), intent(in) :: w0 !! diagonal regularization parameter
        real(dp), intent(inout) :: a(:, :) !! Anderson normal matrix workspace

        integer :: i
        integer :: j
        real(dp) :: value

        if (count <= 0) return
        a(1:count, 1:count) = 0.0_dp
        do i = 1, count
            do j = i, count
                value = dot_product(hist_df(:, i), hist_df(:, j))
                if (i == j) value = value * (1.0_dp + w0 * w0)
                a(i, j) = value
                a(j, i) = value
            end do
        end do
    end subroutine build_anderson_matrix

    subroutine solve_small_system(a, b, x, ok)
        real(dp), intent(in) :: a(:, :) !! square coefficient matrix
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! solution when ok is true
        logical, intent(out) :: ok !! false for singular or invalid input

        real(dp), allocatable :: aa(:, :)
        real(dp), allocatable :: bb(:)
        real(dp) :: factor
        real(dp) :: pivot
        real(dp) :: tmp
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: p

        n = size(b)
        x = 0.0_dp
        ok = .false.
        if (n <= 0 .or. size(x) /= n .or. size(a, 1) /= n .or. size(a, 2) /= n) return
        allocate(aa(n, n), bb(n))
        aa = a
        bb = b

        do k = 1, n - 1
            p = k - 1 + maxloc(abs(aa(k:n, k)), dim=1)
            pivot = aa(p, k)
            if (.not. ieee_is_finite(pivot) .or. abs(pivot) <= tiny(1.0_dp)) return
            if (p /= k) then
                do j = k, n
                    tmp = aa(k, j)
                    aa(k, j) = aa(p, j)
                    aa(p, j) = tmp
                end do
                tmp = bb(k)
                bb(k) = bb(p)
                bb(p) = tmp
            end if
            do i = k + 1, n
                factor = aa(i, k) / aa(k, k)
                aa(i, k) = 0.0_dp
                aa(i, k + 1:n) = aa(i, k + 1:n) - factor * aa(k, k + 1:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do
        if (.not. ieee_is_finite(aa(n, n)) .or. abs(aa(n, n)) <= tiny(1.0_dp)) return

        do i = n, 1, -1
            tmp = bb(i)
            if (i < n) tmp = tmp - dot_product(aa(i, i + 1:n), x(i + 1:n))
            if (.not. ieee_is_finite(aa(i, i)) .or. abs(aa(i, i)) <= tiny(1.0_dp)) return
            x(i) = tmp / aa(i, i)
        end do
        ok = all(ieee_is_finite(x))
    end subroutine solve_small_system

    pure subroutine set_scaled_identity(a, scale)
        real(dp), intent(out) :: a(:, :) !! square matrix overwritten by scale times I
        real(dp), intent(in) :: scale !! diagonal value

        integer :: i

        a = 0.0_dp
        do i = 1, min(size(a, 1), size(a, 2))
            a(i, i) = scale
        end do
    end subroutine set_scaled_identity

    pure logical function same_real(a, b) result(equal)
        real(dp), intent(in) :: a !! first value
        real(dp), intent(in) :: b !! second value

        equal = a <= b .and. a >= b
    end function same_real

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text !! source text
        character(len=len(text)) :: lower
        integer :: code
        integer :: i

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module scifort_optimize_nonlin
