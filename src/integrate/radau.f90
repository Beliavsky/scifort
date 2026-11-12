! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Real-valued fifth-order Radau IIA stepper adapted from SciPy 1.17.0
! scipy/integrate/_ivp/radau.py and its dense finite-difference Jacobian helper.
! The surrounding solve_ivp orchestration remains in scifort_integrate_ivp.
! See CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.

module scifort_integrate_radau
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_integrate_jacobian_fd, only : group_fd_columns
    use scifort_integrate_sparse_jacobian, only : ivp_csc_from_pattern, ivp_sparse_matrix
    use scifort_integrate_sparse_lu, only : ivp_sparse_lu_real, &
        ivp_sparse_lu_complex, sparse_lu_factor
    implicit none
    private

    integer, parameter, public :: radau_status_success = 0
    integer, parameter, public :: radau_status_too_small_step = 1
    integer, parameter, public :: radau_status_singular = 2
    integer, parameter, public :: radau_status_nonfinite = 3
    integer, parameter, public :: radau_status_max_attempts = 4
    integer, parameter, public :: radau_status_invalid_jacobian = 5

    integer, parameter :: newton_maxiter = 6
    real(dp), parameter :: min_factor = 0.2_dp
    real(dp), parameter :: max_factor = 10.0_dp
    real(dp), parameter :: num_jac_diff_reject = 2.0097183471152322e-14_dp
    real(dp), parameter :: num_jac_diff_small = 1.8189894035458565e-12_dp
    real(dp), parameter :: num_jac_diff_big = 1.220703125e-4_dp
    real(dp), parameter :: num_jac_factor_increase = 10.0_dp
    real(dp), parameter :: num_jac_factor_decrease = 0.1_dp
    real(dp), parameter :: num_jac_min_factor = epsilon(1.0_dp) * 1000.0_dp

    abstract interface
        subroutine radau_rhs(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t
        end subroutine radau_rhs

        subroutine radau_rhs_vectorized(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable shared by all state columns
            real(dp), intent(in) :: y(:, :) !! state columns, shape (n,k)
            real(dp), intent(out) :: dydt(:, :) !! derivatives for each state column
        end subroutine radau_rhs_vectorized

        subroutine radau_jacobian(t, y, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: jac(:, :) !! Jacobian d fun / d y, shape (n,n)
        end subroutine radau_jacobian

        subroutine radau_sparse_jacobian(t, y, jac)
            import :: dp, ivp_sparse_matrix
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            type(ivp_sparse_matrix), intent(out) :: jac !! one-based CSC Jacobian d fun / d y
        end subroutine radau_sparse_jacobian
    end interface

    type, public :: radau_state
        real(dp), allocatable :: jac(:, :)
        real(dp), allocatable :: jac_factor(:)
        real(dp), allocatable :: lu_real(:, :)
        integer, allocatable :: piv_real(:)
        complex(dp), allocatable :: lu_complex(:, :)
        type(ivp_sparse_matrix) :: sparse_jac
        type(ivp_sparse_lu_real) :: sparse_lu_real
        type(ivp_sparse_lu_complex) :: sparse_lu_complex
        integer, allocatable :: piv_complex(:)
        real(dp), allocatable :: dense_q(:, :)
        real(dp), allocatable :: dense_y_old(:)
        real(dp) :: dense_t_old = 0.0_dp
        real(dp) :: dense_h = 0.0_dp
        real(dp) :: h_abs = 0.0_dp
        real(dp) :: h_abs_old = 0.0_dp
        real(dp) :: error_norm_old = 0.0_dp
        real(dp) :: newton_tol = 0.0_dp
        logical :: have_error_history = .false.
        logical :: have_dense = .false.
        logical :: lu_valid = .false.
        logical :: current_jac = .true.
        logical :: jacobian_constant = .false.
        logical :: jacobian_sparse = .false.
    end type radau_state

    public :: radau_advance
    public :: radau_dense_state
    public :: radau_initialize
    public :: radau_jacobian
    public :: radau_sparse_jacobian

contains

    subroutine radau_initialize(state, fun, t, y, f, h_abs, rtol, atol, njev, jac, &
            atol_vector, jac_matrix, jac_sparsity, fun_vectorized, jac_sparse_matrix, &
            jac_sparse, rtol_vector)
        type(radau_state), intent(out) :: state !! initialized Radau Jacobian and iteration state
        procedure(radau_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! initial independent-variable value
        real(dp), intent(in) :: y(:) !! initial state
        real(dp), intent(in) :: f(:) !! derivative at the initial state
        real(dp), intent(in) :: h_abs !! selected initial step magnitude
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol !! scalar absolute tolerance
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        procedure(radau_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        real(dp), intent(in), optional :: jac_matrix(:, :) !! constant dense Jacobian, shape (n,n)
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(radau_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for finite-difference Jacobians
        type(ivp_sparse_matrix), intent(in), optional :: &
            jac_sparse_matrix !! constant CSC Jacobian
        procedure(radau_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent CSC Jacobian callback
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        integer :: n
        real(dp) :: atol_use(size(y))
        real(dp) :: rtol_use(size(y))
        logical :: use_sparse

        n = size(y)
        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector
        allocate(state%jac_factor(n), state%dense_q(n, 3), state%dense_y_old(n))
        use_sparse = present(jac_sparse_matrix) .or. present(jac_sparse) .or. &
            (present(jac_sparsity) .and. .not. present(jac) .and. .not. present(jac_matrix))
        if (use_sparse) then
            allocate(state%lu_real(0, 0), state%piv_real(0), &
                state%lu_complex(0, 0), state%piv_complex(0))
        else
            allocate(state%jac(n, n), state%lu_real(n, n), state%piv_real(n), &
                state%lu_complex(n, n), state%piv_complex(n))
        end if
        state%h_abs = h_abs
        state%newton_tol = max(10.0_dp * epsilon(1.0_dp) / maxval(rtol_use), &
            min(0.03_dp, sqrt(maxval(rtol_use))))
        state%jac_factor = sqrt(epsilon(1.0_dp))
        state%dense_q = 0.0_dp
        state%dense_y_old = y
        state%have_error_history = .false.
        state%have_dense = .false.
        state%lu_valid = .false.
        state%current_jac = .true.
        state%jacobian_sparse = use_sparse
        state%jacobian_constant = present(jac_matrix) .or. present(jac_sparse_matrix)
        if (present(jac_sparse_matrix)) then
            state%sparse_jac = jac_sparse_matrix
        else if (present(jac_matrix)) then
            state%jac = jac_matrix
        else if (present(jac_sparse)) then
            call jac_sparse(t, y, state%sparse_jac)
            njev = njev + 1
        else if (present(jac)) then
            call jac(t, y, state%jac)
            njev = njev + 1
        else if (state%jacobian_sparse) then
            call ivp_csc_from_pattern(jac_sparsity, state%sparse_jac)
            call numerical_sparse_jacobian(fun, t, y, f, atol_use, state%jac_factor, &
                state%sparse_jac, jac_sparsity, fun_vectorized)
            njev = njev + 1
        else
            call numerical_jacobian(fun, t, y, f, atol_use, state%jac_factor, state%jac, &
                jac_sparsity, fun_vectorized)
            njev = njev + 1
        end if
    end subroutine radau_initialize

    subroutine radau_advance(state, fun, t, y, f, t_bound, direction, rtol, atol, &
            max_step, y_new, f_new, t_new, dense_q, dense_h, nfev, njev, nlu, &
            attempt_limit, attempts, rejected, status, jac, atol_vector, jac_sparsity, &
            fun_vectorized, jac_sparse, rtol_vector)
        type(radau_state), intent(inout) :: state !! mutable Radau Jacobian and iteration state
        procedure(radau_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: y(:) !! current state
        real(dp), intent(in) :: f(:) !! derivative at the current state
        real(dp), intent(in) :: t_bound !! final integration bound
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol !! scalar absolute tolerance
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        real(dp), intent(out) :: y_new(:) !! accepted state
        real(dp), intent(out) :: f_new(:) !! derivative at the accepted state
        real(dp), intent(out) :: t_new !! accepted independent-variable value
        real(dp), intent(out) :: dense_q(:, :) !! cubic dense-output coefficients, first 3 columns
        real(dp), intent(out) :: dense_h !! signed accepted step used by dense output
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        integer, intent(inout) :: nlu !! accumulated LU factorization count
        integer, intent(in) :: attempt_limit !! maximum trial steps available to this call
        integer, intent(out) :: attempts !! trial steps consumed by this accepted step
        integer, intent(out) :: rejected !! rejected trial-step count
        integer, intent(out) :: status !! Radau status code
        procedure(radau_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(radau_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for finite-difference Jacobians
        procedure(radau_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent CSC Jacobian callback
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        real(dp), parameter :: s6 = sqrt(6.0_dp)
        real(dp), parameter :: c(3) = [(4.0_dp - s6) / 10.0_dp, &
            (4.0_dp + s6) / 10.0_dp, 1.0_dp]
        real(dp), parameter :: e(3) = [(-13.0_dp - 7.0_dp * s6) / 3.0_dp, &
            (-13.0_dp + 7.0_dp * s6) / 3.0_dp, -1.0_dp / 3.0_dp]
        real(dp), parameter :: p(3, 3) = reshape([ &
            13.0_dp / 3.0_dp + 7.0_dp * s6 / 3.0_dp, &
            13.0_dp / 3.0_dp - 7.0_dp * s6 / 3.0_dp, 1.0_dp / 3.0_dp, &
            -23.0_dp / 3.0_dp - 22.0_dp * s6 / 3.0_dp, &
            -23.0_dp / 3.0_dp + 22.0_dp * s6 / 3.0_dp, -8.0_dp / 3.0_dp, &
            10.0_dp / 3.0_dp + 5.0_dp * s6, &
            10.0_dp / 3.0_dp - 5.0_dp * s6, 10.0_dp / 3.0_dp], [3, 3])
        real(dp), parameter :: mu_real = 3.0_dp + 3.0_dp**(2.0_dp / 3.0_dp) - &
            3.0_dp**(1.0_dp / 3.0_dp)
        complex(dp), parameter :: mu_complex = cmplx( &
            3.0_dp + 0.5_dp * (3.0_dp**(1.0_dp / 3.0_dp) - &
            3.0_dp**(2.0_dp / 3.0_dp)), &
            -0.5_dp * (3.0_dp**(5.0_dp / 6.0_dp) + 3.0_dp**(7.0_dp / 6.0_dp)), &
            kind=dp)

        real(dp) :: z(size(y), 3)
        real(dp) :: z0(size(y), 3)
        real(dp) :: scale(size(y))
        real(dp) :: ze(size(y))
        real(dp) :: error(size(y))
        real(dp) :: atol_use(size(y))
        real(dp) :: rtol_use(size(y))
        real(dp) :: f_error(size(y))
        real(dp) :: h_abs
        real(dp) :: h_abs_entry
        real(dp) :: h_abs_old
        real(dp) :: error_norm_old
        real(dp) :: h
        real(dp) :: min_step
        real(dp) :: error_norm
        real(dp) :: safety_newton
        real(dp) :: factor
        real(dp) :: rate
        integer :: n_iter
        integer :: info
        integer :: i
        logical :: converged
        logical :: step_accepted
        logical :: recompute_jac
        logical :: have_history
        logical :: current_jac

        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector

        status = radau_status_success
        attempts = 0
        rejected = 0
        dense_q = 0.0_dp
        dense_h = 0.0_dp
        y_new = y
        f_new = f
        t_new = t

        h_abs_entry = state%h_abs
        h_abs = state%h_abs
        have_history = state%have_error_history
        h_abs_old = state%h_abs_old
        error_norm_old = state%error_norm_old
        min_step = 10.0_dp * abs(nearest(t, direction) - t)
        if (h_abs > max_step) then
            h_abs = max_step
            have_history = .false.
        else if (h_abs < min_step) then
            h_abs = min_step
            have_history = .false.
        end if

        current_jac = state%current_jac .or. state%jacobian_constant
        step_accepted = .false.
        do while (.not. step_accepted)
            if (attempts >= attempt_limit) then
                status = radau_status_max_attempts
                return
            end if
            attempts = attempts + 1
            if (h_abs < min_step) then
                status = radau_status_too_small_step
                return
            end if

            h = h_abs * direction
            t_new = t + h
            if (direction * (t_new - t_bound) > 0.0_dp) t_new = t_bound
            h = t_new - t
            h_abs = abs(h)

            if (state%have_dense) then
                do i = 1, 3
                    call radau_dense_state(state%dense_t_old, state%dense_y_old, &
                        state%dense_h, state%dense_q, t + h * c(i), z0(:, i))
                    z0(:, i) = z0(:, i) - y
                end do
            else
                z0 = 0.0_dp
            end if
            scale = max(atol_use + abs(y) * rtol_use, tiny(1.0_dp))

            converged = .false.
            do while (.not. converged)
                if (.not. state%lu_valid) then
                    if (state%jacobian_sparse) then
                        call sparse_lu_factor(state%sparse_jac, mu_real / h, -1.0_dp, &
                            state%sparse_lu_real, info)
                    else
                        call factor_real_matrix(state%jac, mu_real / h, state%lu_real, &
                            state%piv_real, info)
                    end if
                    nlu = nlu + 1
                    if (info /= 0) then
                        status = radau_status_singular
                        return
                    end if
                    if (state%jacobian_sparse) then
                        call sparse_lu_factor(state%sparse_jac, mu_complex / h, &
                            cmplx(-1.0_dp, 0.0_dp, kind=dp), &
                            state%sparse_lu_complex, info)
                    else
                        call factor_complex_matrix(state%jac, mu_complex / h, &
                            state%lu_complex, state%piv_complex, info)
                    end if
                    nlu = nlu + 1
                    if (info /= 0) then
                        status = radau_status_singular
                        return
                    end if
                    state%lu_valid = .true.
                end if

                call solve_collocation_system(fun, t, y, h, z0, scale, &
                    state%newton_tol, state%lu_real, state%piv_real, &
                    state%lu_complex, state%piv_complex, z, converged, n_iter, &
                    rate, nfev, status, state%sparse_lu_real, state%sparse_lu_complex)
                if (status /= radau_status_success) return
                if (converged) exit
                if (current_jac) exit

                if (present(jac_sparse)) then
                    call jac_sparse(t, y, state%sparse_jac)
                    if (.not. sparse_jacobian_is_valid(state%sparse_jac, size(y))) then
                        status = radau_status_invalid_jacobian
                        return
                    end if
                else if (present(jac)) then
                    call jac(t, y, state%jac)
                else if (state%jacobian_sparse) then
                    call numerical_sparse_jacobian(fun, t, y, f, atol_use, &
                        state%jac_factor, state%sparse_jac, jac_sparsity, fun_vectorized)
                else
                    call numerical_jacobian(fun, t, y, f, atol_use, state%jac_factor, &
                        state%jac, jac_sparsity, fun_vectorized)
                end if
                njev = njev + 1
                if (.not. state%jacobian_sparse) then
                    if (.not. all(ieee_is_finite(state%jac))) then
                        status = radau_status_nonfinite
                        return
                    end if
                end if
                current_jac = .true.
                state%lu_valid = .false.
            end do

            if (.not. converged) then
                h_abs = 0.5_dp * h_abs
                state%lu_valid = .false.
                rejected = rejected + 1
                cycle
            end if

            y_new = y + z(:, 3)
            ze = matmul(z, e) / h
            if (state%jacobian_sparse) then
                call state%sparse_lu_real%solve(f + ze, error, info)
            else
                call lu_solve_real(state%lu_real, state%piv_real, f + ze, error, info)
            end if
            if (info /= 0) then
                status = radau_status_singular
                return
            end if
            scale = max(atol_use + max(abs(y), abs(y_new)) * rtol_use, tiny(1.0_dp))
            error_norm = rms_norm(error / scale)
            safety_newton = 0.9_dp * real(2 * newton_maxiter + 1, dp) / &
                real(2 * newton_maxiter + n_iter, dp)

            if (rejected > 0 .and. error_norm > 1.0_dp) then
                call fun(t, y + error, f_error)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(f_error))) then
                    error_norm = huge(1.0_dp)
                else
                    if (state%jacobian_sparse) then
                        call state%sparse_lu_real%solve(f_error + ze, error, info)
                    else
                        call lu_solve_real(state%lu_real, state%piv_real, f_error + ze, &
                            error, info)
                    end if
                    if (info /= 0) then
                        status = radau_status_singular
                        return
                    end if
                    error_norm = rms_norm(error / scale)
                end if
            end if

            if (.not. ieee_is_finite(error_norm)) then
                factor = min_factor
                h_abs = h_abs * factor
                state%lu_valid = .false.
                rejected = rejected + 1
            else if (error_norm > 1.0_dp) then
                factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old, &
                    have_history)
                h_abs = h_abs * max(min_factor, safety_newton * factor)
                state%lu_valid = .false.
                rejected = rejected + 1
            else
                step_accepted = .true.
            end if
        end do

        recompute_jac = .not. state%jacobian_constant .and. &
            n_iter > 2 .and. rate > 1.0e-3_dp
        factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old, have_history)
        factor = min(max_factor, safety_newton * factor)
        if (.not. recompute_jac .and. factor < 1.2_dp) then
            factor = 1.0_dp
        else
            state%lu_valid = .false.
        end if

        call fun(t_new, y_new, f_new)
        nfev = nfev + 1
        if (.not. all(ieee_is_finite(f_new))) then
            status = radau_status_nonfinite
            return
        end if

        if (recompute_jac) then
            if (present(jac_sparse)) then
                call jac_sparse(t_new, y_new, state%sparse_jac)
                if (.not. sparse_jacobian_is_valid(state%sparse_jac, size(y_new))) then
                    status = radau_status_invalid_jacobian
                    return
                end if
            else if (present(jac)) then
                call jac(t_new, y_new, state%jac)
            else if (state%jacobian_sparse) then
                call numerical_sparse_jacobian(fun, t_new, y_new, f_new, atol_use, &
                    state%jac_factor, state%sparse_jac, jac_sparsity, fun_vectorized)
            else
                call numerical_jacobian(fun, t_new, y_new, f_new, atol_use, &
                    state%jac_factor, state%jac, jac_sparsity, fun_vectorized)
            end if
            njev = njev + 1
            if (.not. state%jacobian_sparse) then
                if (.not. all(ieee_is_finite(state%jac))) then
                    status = radau_status_nonfinite
                    return
                end if
            end if
            state%current_jac = .true.
        else
            state%current_jac = state%jacobian_constant
        end if

        state%h_abs_old = h_abs_entry
        state%error_norm_old = error_norm
        state%have_error_history = .true.
        state%h_abs = h_abs * factor

        dense_q(:, 1:3) = matmul(z, p)
        dense_h = h
        state%dense_t_old = t
        state%dense_y_old = y
        state%dense_h = h
        state%dense_q = dense_q(:, 1:3)
        state%have_dense = .true.
    end subroutine radau_advance

    subroutine solve_collocation_system(fun, t, y, h, z0, scale, tol, lu_real, &
            piv_real, lu_complex, piv_complex, z, converged, n_iter, rate, nfev, status, &
            sparse_lu_real, sparse_lu_complex)
        procedure(radau_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! step start time
        real(dp), intent(in) :: y(:) !! state at the step start
        real(dp), intent(in) :: h !! signed trial step
        real(dp), intent(in) :: z0(:, :) !! initial stage increments, shape (n,3)
        real(dp), intent(in) :: scale(:) !! component Newton scales
        real(dp), intent(in) :: tol !! Newton convergence tolerance
        real(dp), intent(in) :: lu_real(:, :) !! LU factors of the real transformed system
        integer, intent(in) :: piv_real(:) !! pivots for the real transformed system
        complex(dp), intent(in) :: lu_complex(:, :) !! LU factors of complex transformed system
        integer, intent(in) :: piv_complex(:) !! pivots for the complex transformed system
        real(dp), intent(out) :: z(:, :) !! converged stage increments, shape (n,3)
        logical, intent(out) :: converged !! true when Newton convergence is achieved
        integer, intent(out) :: n_iter !! Newton iterations performed
        real(dp), intent(out) :: rate !! final estimated Newton convergence rate
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(out) :: status !! Radau status code
        type(ivp_sparse_lu_real), intent(in) :: &
            sparse_lu_real !! sparse real factors when a constant CSC Jacobian is active
        type(ivp_sparse_lu_complex), intent(in) :: &
            sparse_lu_complex !! sparse complex factors for the transformed Newton system

        real(dp), parameter :: s6 = sqrt(6.0_dp)
        real(dp), parameter :: c(3) = [(4.0_dp - s6) / 10.0_dp, &
            (4.0_dp + s6) / 10.0_dp, 1.0_dp]
        real(dp), parameter :: tmat(3, 3) = reshape([ &
            0.09443876248897524_dp, 0.25021312296533332_dp, 1.0_dp, &
            -0.14125529502095421_dp, 0.20412935229379994_dp, 1.0_dp, &
            0.03002919410514742_dp, -0.38294211275726192_dp, 0.0_dp], [3, 3])
        real(dp), parameter :: ti_real(3) = [4.17871859155190428_dp, &
            0.32768282076106237_dp, 0.52337644549944951_dp]
        complex(dp), parameter :: ti_complex(3) = [ &
            cmplx(-4.17871859155190428_dp, 0.50287263494578682_dp, kind=dp), &
            cmplx(-0.32768282076106237_dp, -2.57192694985560522_dp, kind=dp), &
            cmplx(0.47662355450055044_dp, 0.59603920482822492_dp, kind=dp)]
        real(dp), parameter :: mu_real = 3.0_dp + 3.0_dp**(2.0_dp / 3.0_dp) - &
            3.0_dp**(1.0_dp / 3.0_dp)
        complex(dp), parameter :: mu_complex = cmplx( &
            3.0_dp + 0.5_dp * (3.0_dp**(1.0_dp / 3.0_dp) - &
            3.0_dp**(2.0_dp / 3.0_dp)), &
            -0.5_dp * (3.0_dp**(5.0_dp / 6.0_dp) + 3.0_dp**(7.0_dp / 6.0_dp)), &
            kind=dp)

        real(dp) :: w(size(y), 3)
        real(dp) :: fstage(size(y), 3)
        real(dp) :: f_real(size(y))
        real(dp) :: dw_real(size(y))
        real(dp) :: dw(size(y), 3)
        complex(dp) :: f_complex(size(y))
        complex(dp) :: dw_complex(size(y))
        real(dp) :: dwnorm
        real(dp) :: dwnorm_old
        integer :: k
        integer :: i
        integer :: info
        logical :: have_old_norm

        call transform_z_to_w(z0, w)
        z = z0
        converged = .false.
        have_old_norm = .false.
        dwnorm_old = 0.0_dp
        rate = 0.0_dp
        status = radau_status_success
        n_iter = 0

        do k = 1, newton_maxiter
            n_iter = k
            do i = 1, 3
                call fun(t + h * c(i), y + z(:, i), fstage(:, i))
                nfev = nfev + 1
            end do
            if (.not. all(ieee_is_finite(fstage))) return

            f_real = matmul(fstage, ti_real) - (mu_real / h) * w(:, 1)
            f_complex = matmul(cmplx(fstage, 0.0_dp, kind=dp), ti_complex) - &
                (mu_complex / h) * cmplx(w(:, 2), w(:, 3), kind=dp)
            if (sparse_lu_real%n > 0) then
                call sparse_lu_real%solve(f_real, dw_real, info)
            else
                call lu_solve_real(lu_real, piv_real, f_real, dw_real, info)
            end if
            if (info /= 0) then
                status = radau_status_singular
                return
            end if
            if (sparse_lu_complex%n > 0) then
                call sparse_lu_complex%solve(f_complex, dw_complex, info)
            else
                call lu_solve_complex(lu_complex, piv_complex, f_complex, dw_complex, info)
            end if
            if (info /= 0) then
                status = radau_status_singular
                return
            end if
            dw(:, 1) = dw_real
            dw(:, 2) = real(dw_complex, kind=dp)
            dw(:, 3) = aimag(dw_complex)
            if (.not. all(ieee_is_finite(dw))) return

            dwnorm = rms_norm_2d_scaled(dw, scale)
            if (have_old_norm) then
                if (dwnorm_old == 0.0_dp) then
                    rate = 0.0_dp
                else
                    rate = dwnorm / dwnorm_old
                end if
                if (rate >= 1.0_dp) return
                if (rate**real(newton_maxiter - k, dp) / (1.0_dp - rate) * &
                        dwnorm > tol) return
            else
                rate = 0.0_dp
            end if

            w = w + dw
            call transform_w_to_z(w, z)
            if (dwnorm == 0.0_dp) then
                converged = .true.
                return
            end if
            if (have_old_norm) then
                if (rate / (1.0_dp - rate) * dwnorm < tol) then
                    converged = .true.
                    return
                end if
            end if
            dwnorm_old = dwnorm
            have_old_norm = .true.
        end do

    contains

        subroutine transform_z_to_w(zin, wout)
            real(dp), intent(in) :: zin(:, :) !! stage increments in collocation coordinates
            real(dp), intent(out) :: wout(:, :) !! transformed stage increments
            integer :: j

            do j = 1, 3
                wout(:, j) = ti_row_dot(j, zin)
            end do
        end subroutine transform_z_to_w

        function ti_row_dot(row, zin) result(value)
            integer, intent(in) :: row !! row of SciPy's inverse transformation matrix
            real(dp), intent(in) :: zin(:, :) !! stage increments
            real(dp) :: value(size(zin, 1))
            real(dp), parameter :: ti(3, 3) = reshape([ &
                4.17871859155190428_dp, -4.17871859155190428_dp, &
                0.50287263494578682_dp, &
                0.32768282076106237_dp, -0.32768282076106237_dp, &
                -2.57192694985560522_dp, &
                0.52337644549944951_dp, 0.47662355450055044_dp, &
                0.59603920482822492_dp], [3, 3])
            integer :: j

            value = 0.0_dp
            do j = 1, 3
                value = value + ti(row, j) * zin(:, j)
            end do
        end function ti_row_dot

        subroutine transform_w_to_z(win, zout)
            real(dp), intent(in) :: win(:, :) !! transformed stage increments
            real(dp), intent(out) :: zout(:, :) !! stage increments in collocation coordinates
            integer :: ii
            integer :: jj

            zout = 0.0_dp
            do ii = 1, 3
                do jj = 1, 3
                    zout(:, ii) = zout(:, ii) + tmat(ii, jj) * win(:, jj)
                end do
            end do
        end subroutine transform_w_to_z

    end subroutine solve_collocation_system

    subroutine radau_dense_state(t_old, y_old, h, q, tq, yq)
        real(dp), intent(in) :: t_old !! accepted-step start time
        real(dp), intent(in) :: y_old(:) !! state at the accepted-step start
        real(dp), intent(in) :: h !! signed accepted step size
        real(dp), intent(in) :: q(:, :) !! cubic dense-output coefficient matrix, shape (n,3)
        real(dp), intent(in) :: tq !! interpolation point
        real(dp), intent(out) :: yq(:) !! interpolated state

        real(dp) :: x
        real(dp) :: p1
        real(dp) :: p2
        real(dp) :: p3

        if (h == 0.0_dp) then
            yq = y_old
            return
        end if
        x = (tq - t_old) / h
        p1 = x
        p2 = p1 * x
        p3 = p2 * x
        yq = y_old + q(:, 1) * p1 + q(:, 2) * p2 + q(:, 3) * p3
    end subroutine radau_dense_state

    pure real(dp) function predict_factor(h_abs, h_abs_old, error_norm, &
            error_norm_old, have_history) result(value)
        real(dp), intent(in) :: h_abs !! current accepted step magnitude
        real(dp), intent(in) :: h_abs_old !! previous step magnitude used by predictor
        real(dp), intent(in) :: error_norm !! current normalized local error
        real(dp), intent(in) :: error_norm_old !! previous normalized local error
        logical, intent(in) :: have_history !! whether both previous quantities are valid
        real(dp) :: multiplier

        if (.not. have_history .or. error_norm == 0.0_dp) then
            multiplier = 1.0_dp
        else
            multiplier = h_abs / h_abs_old * (error_norm_old / error_norm)**0.25_dp
        end if
        if (error_norm == 0.0_dp) then
            value = huge(1.0_dp)
        else if (.not. ieee_is_finite(error_norm) .or. error_norm < 0.0_dp) then
            value = 0.0_dp
        else
            value = min(1.0_dp, multiplier) * error_norm**(-0.25_dp)
        end if
    end function predict_factor

    subroutine numerical_jacobian(fun, t, y, f, atol, factor, jac, jac_sparsity, &
            fun_vectorized)
        procedure(radau_rhs) :: fun !! scalar right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        real(dp), intent(out) :: jac(:, :) !! finite-difference Jacobian
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(radau_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS accepting state columns for Jacobian differences

        real(dp) :: f_batch(size(y), size(y))
        real(dp) :: f_new(size(y), size(y))
        real(dp) :: f_retry(size(y))
        real(dp) :: h(size(y))
        real(dp) :: h_retry
        real(dp) :: max_diff
        real(dp) :: max_diff_new
        real(dp) :: new_factor
        real(dp) :: scale_new
        real(dp) :: scale_ref
        real(dp) :: y_batch(size(y), size(y))
        real(dp) :: y_scale(size(y))
        real(dp) :: y_trial(size(y))
        integer :: g
        integer :: groups(size(y))
        integer :: i
        integer :: imax
        integer :: j
        integer :: n_groups
        logical :: use_retry
        logical :: use_sparsity

        use_sparsity = present(jac_sparsity)
        do j = 1, size(y)
            if (f(j) >= 0.0_dp) then
                y_scale(j) = max(atol(j), abs(y(j)))
            else
                y_scale(j) = -max(atol(j), abs(y(j)))
            end if
            h(j) = (y(j) + factor(j) * y_scale(j)) - y(j)
            do while (h(j) == 0.0_dp)
                factor(j) = factor(j) * 10.0_dp
                h(j) = (y(j) + factor(j) * y_scale(j)) - y(j)
            end do
        end do

        if (use_sparsity) then
            call group_fd_columns(jac_sparsity, groups, n_groups)
        else
            n_groups = size(y)
            groups = [(j, j = 1, size(y))]
        end if
        f_new = spread(f, dim=2, ncopies=size(y))
        if (n_groups > 0) then
            if (present(fun_vectorized)) then
                do g = 1, n_groups
                    y_batch(:, g) = y
                    do j = 1, size(y)
                        if (groups(j) == g) y_batch(j, g) = y(j) + h(j)
                    end do
                end do
                call fun_vectorized(t, y_batch(:, 1:n_groups), f_batch(:, 1:n_groups))
                do j = 1, size(y)
                    if (groups(j) > 0) f_new(:, j) = f_batch(:, groups(j))
                end do
            else if (use_sparsity) then
                do g = 1, n_groups
                    y_trial = y
                    do j = 1, size(y)
                        if (groups(j) == g) y_trial(j) = y(j) + h(j)
                    end do
                    call fun(t, y_trial, f_retry)
                    do j = 1, size(y)
                        if (groups(j) == g) f_new(:, j) = f_retry
                    end do
                end do
            else
                do j = 1, size(y)
                    y_trial = y
                    y_trial(j) = y(j) + h(j)
                    call fun(t, y_trial, f_new(:, j))
                end do
            end if
        end if

        jac = 0.0_dp
        do j = 1, size(y)
            if (use_sparsity) then
                if (groups(j) == 0) cycle
                imax = maxloc(abs(f_new(:, j) - f), dim=1, mask=jac_sparsity(:, j))
            else
                imax = maxloc(abs(f_new(:, j) - f), dim=1)
            end if
            max_diff = abs(f_new(imax, j) - f(imax))
            scale_ref = max(abs(f(imax)), abs(f_new(imax, j)))
            use_retry = max_diff < num_jac_diff_reject * scale_ref
            if (use_retry) then
                new_factor = num_jac_factor_increase * factor(j)
                h_retry = (y(j) + new_factor * y_scale(j)) - y(j)
                y_trial = y
                y_trial(j) = y(j) + h_retry
                if (present(fun_vectorized)) then
                    y_batch(:, 1) = y_trial
                    call fun_vectorized(t, y_batch(:, 1:1), f_batch(:, 1:1))
                    f_retry = f_batch(:, 1)
                else
                    call fun(t, y_trial, f_retry)
                end if
                if (use_sparsity) then
                    i = maxloc(abs(f_retry - f), dim=1, mask=jac_sparsity(:, j))
                else
                    i = maxloc(abs(f_retry - f), dim=1)
                end if
                max_diff_new = abs(f_retry(i) - f(i))
                scale_new = max(abs(f(i)), abs(f_retry(i)))
                if (max_diff * scale_new < max_diff_new * scale_ref) then
                    factor(j) = new_factor
                    h(j) = h_retry
                    f_new(:, j) = f_retry
                    max_diff = max_diff_new
                    scale_ref = scale_new
                end if
            end if
            if (use_sparsity) then
                where (jac_sparsity(:, j))
                    jac(:, j) = (f_new(:, j) - f) / h(j)
                elsewhere
                    jac(:, j) = 0.0_dp
                end where
            else
                jac(:, j) = (f_new(:, j) - f) / h(j)
            end if
            if (max_diff < num_jac_diff_small * scale_ref) then
                factor(j) = factor(j) * num_jac_factor_increase
            else if (max_diff > num_jac_diff_big * scale_ref) then
                factor(j) = factor(j) * num_jac_factor_decrease
            end if
            factor(j) = max(factor(j), num_jac_min_factor)
        end do
    end subroutine numerical_jacobian


    subroutine numerical_sparse_jacobian(fun, t, y, f, atol, factor, jac, sparsity, &
            fun_vectorized)
        procedure(radau_rhs) :: fun !! scalar right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        type(ivp_sparse_matrix), intent(inout) :: jac !! CSC Jacobian with fixed sparsity structure
        logical, intent(in) :: sparsity(:, :) !! finite-difference nonzero pattern
        procedure(radau_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS accepting state columns for Jacobian differences

        real(dp), allocatable :: f_group(:, :)
        real(dp), allocatable :: y_batch(:, :)
        real(dp) :: f_col(size(y))
        real(dp) :: f_retry(size(y))
        real(dp) :: h(size(y))
        real(dp) :: h_retry
        real(dp) :: max_diff
        real(dp) :: max_diff_new
        real(dp) :: new_factor
        real(dp) :: scale_new
        real(dp) :: scale_ref
        real(dp) :: y_scale(size(y))
        real(dp) :: y_trial(size(y))
        integer :: g
        integer :: groups(size(y))
        integer :: i
        integer :: imax
        integer :: j
        integer :: k
        integer :: n_groups
        logical :: use_retry

        do j = 1, size(y)
            if (f(j) >= 0.0_dp) then
                y_scale(j) = max(atol(j), abs(y(j)))
            else
                y_scale(j) = -max(atol(j), abs(y(j)))
            end if
            h(j) = (y(j) + factor(j) * y_scale(j)) - y(j)
            do while (h(j) == 0.0_dp)
                factor(j) = factor(j) * 10.0_dp
                h(j) = (y(j) + factor(j) * y_scale(j)) - y(j)
            end do
        end do

        call group_fd_columns(sparsity, groups, n_groups)
        allocate(f_group(size(y), max(1, n_groups)), y_batch(size(y), max(1, n_groups)))
        f_group = spread(f, dim=2, ncopies=max(1, n_groups))
        if (n_groups > 0) then
            if (present(fun_vectorized)) then
                do g = 1, n_groups
                    y_batch(:, g) = y
                    do j = 1, size(y)
                        if (groups(j) == g) y_batch(j, g) = y(j) + h(j)
                    end do
                end do
                call fun_vectorized(t, y_batch(:, 1:n_groups), f_group(:, 1:n_groups))
            else
                do g = 1, n_groups
                    y_trial = y
                    do j = 1, size(y)
                        if (groups(j) == g) y_trial(j) = y(j) + h(j)
                    end do
                    call fun(t, y_trial, f_group(:, g))
                end do
            end if
        end if

        jac%data = 0.0_dp
        do j = 1, size(y)
            g = groups(j)
            if (g == 0) cycle
            f_col = f_group(:, g)
            imax = maxloc(abs(f_col - f), dim=1, mask=sparsity(:, j))
            max_diff = abs(f_col(imax) - f(imax))
            scale_ref = max(abs(f(imax)), abs(f_col(imax)))
            use_retry = max_diff < num_jac_diff_reject * scale_ref
            if (use_retry) then
                new_factor = num_jac_factor_increase * factor(j)
                h_retry = (y(j) + new_factor * y_scale(j)) - y(j)
                y_trial = y
                y_trial(j) = y(j) + h_retry
                if (present(fun_vectorized)) then
                    y_batch(:, 1) = y_trial
                    call fun_vectorized(t, y_batch(:, 1:1), f_group(:, 1:1))
                    f_retry = f_group(:, 1)
                else
                    call fun(t, y_trial, f_retry)
                end if
                i = maxloc(abs(f_retry - f), dim=1, mask=sparsity(:, j))
                max_diff_new = abs(f_retry(i) - f(i))
                scale_new = max(abs(f(i)), abs(f_retry(i)))
                if (max_diff * scale_new < max_diff_new * scale_ref) then
                    factor(j) = new_factor
                    h(j) = h_retry
                    f_col = f_retry
                    max_diff = max_diff_new
                    scale_ref = scale_new
                end if
            end if
            do k = jac%indptr(j), jac%indptr(j + 1) - 1
                i = jac%indices(k)
                jac%data(k) = (f_col(i) - f(i)) / h(j)
            end do
            if (max_diff < num_jac_diff_small * scale_ref) then
                factor(j) = factor(j) * num_jac_factor_increase
            else if (max_diff > num_jac_diff_big * scale_ref) then
                factor(j) = factor(j) * num_jac_factor_decrease
            end if
            factor(j) = max(factor(j), num_jac_min_factor)
        end do
    end subroutine numerical_sparse_jacobian

    logical function sparse_jacobian_is_valid(jac, n) result(valid)
        type(ivp_sparse_matrix), intent(in) :: jac !! sparse Jacobian returned by a callback
        integer, intent(in) :: n !! required square matrix order

        valid = jac%is_valid()
        if (.not. valid) return
        valid = jac%nrows == n .and. jac%ncols == n
    end function sparse_jacobian_is_valid



    subroutine factor_real_matrix(jac, diagonal, lu, piv, info)
        real(dp), intent(in) :: jac(:, :) !! Jacobian matrix
        real(dp), intent(in) :: diagonal !! scalar diagonal term mu_real/h
        real(dp), intent(out) :: lu(:, :) !! LU factors of diagonal*I-J
        integer, intent(out) :: piv(:) !! row pivots for LU solve
        integer, intent(out) :: info !! zero on success; nonzero if numerically singular
        integer :: i

        lu = -jac
        do i = 1, size(jac, 1)
            lu(i, i) = lu(i, i) + diagonal
        end do
        call lu_factor_real(lu, piv, info)
    end subroutine factor_real_matrix

    subroutine factor_complex_matrix(jac, diagonal, lu, piv, info)
        real(dp), intent(in) :: jac(:, :) !! Jacobian matrix
        complex(dp), intent(in) :: diagonal !! scalar diagonal term mu_complex/h
        complex(dp), intent(out) :: lu(:, :) !! LU factors of diagonal*I-J
        integer, intent(out) :: piv(:) !! row pivots for LU solve
        integer, intent(out) :: info !! zero on success; nonzero if numerically singular
        integer :: i

        lu = -cmplx(jac, 0.0_dp, kind=dp)
        do i = 1, size(jac, 1)
            lu(i, i) = lu(i, i) + diagonal
        end do
        call lu_factor_complex(lu, piv, info)
    end subroutine factor_complex_matrix

    subroutine lu_factor_real(a, piv, info)
        real(dp), intent(inout) :: a(:, :) !! square matrix, overwritten by compact LU factors
        integer, intent(out) :: piv(:) !! row pivot selected at each elimination column
        integer, intent(out) :: info !! zero on success; nonzero for numerical singularity
        real(dp) :: row_tmp(size(a, 2))
        real(dp) :: pivot_scale
        integer :: i
        integer :: k
        integer :: n
        integer :: pivot

        n = size(a, 1)
        info = 0
        piv = [(i, i = 1, n)]
        do k = 1, n
            pivot = k - 1 + maxloc(abs(a(k:n, k)), dim=1)
            pivot_scale = max(1.0_dp, maxval(abs(a(k:n, k))))
            if (abs(a(pivot, k)) <= epsilon(1.0_dp) * pivot_scale) then
                info = k
                return
            end if
            piv(k) = pivot
            if (pivot /= k) then
                row_tmp = a(k, :)
                a(k, :) = a(pivot, :)
                a(pivot, :) = row_tmp
            end if
            if (k < n) then
                a(k + 1:n, k) = a(k + 1:n, k) / a(k, k)
                do i = k + 1, n
                    a(i, k + 1:n) = a(i, k + 1:n) - a(i, k) * a(k, k + 1:n)
                end do
            end if
        end do
    end subroutine lu_factor_real

    subroutine lu_factor_complex(a, piv, info)
        complex(dp), intent(inout) :: a(:, :) !! square matrix, overwritten by compact LU factors
        integer, intent(out) :: piv(:) !! row pivot selected at each elimination column
        integer, intent(out) :: info !! zero on success; nonzero for numerical singularity
        complex(dp) :: row_tmp(size(a, 2))
        real(dp) :: pivot_scale
        integer :: i
        integer :: k
        integer :: n
        integer :: pivot

        n = size(a, 1)
        info = 0
        piv = [(i, i = 1, n)]
        do k = 1, n
            pivot = k - 1 + maxloc(abs(a(k:n, k)), dim=1)
            pivot_scale = max(1.0_dp, maxval(abs(a(k:n, k))))
            if (abs(a(pivot, k)) <= epsilon(1.0_dp) * pivot_scale) then
                info = k
                return
            end if
            piv(k) = pivot
            if (pivot /= k) then
                row_tmp = a(k, :)
                a(k, :) = a(pivot, :)
                a(pivot, :) = row_tmp
            end if
            if (k < n) then
                a(k + 1:n, k) = a(k + 1:n, k) / a(k, k)
                do i = k + 1, n
                    a(i, k + 1:n) = a(i, k + 1:n) - a(i, k) * a(k, k + 1:n)
                end do
            end if
        end do
    end subroutine lu_factor_complex

    subroutine lu_solve_real(lu, piv, b, x, info)
        real(dp), intent(in) :: lu(:, :) !! compact LU factors
        integer, intent(in) :: piv(:) !! row pivots from factorization
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! solution of the factored linear system
        integer, intent(out) :: info !! zero on success; nonzero for singular diagonal
        real(dp) :: tmp
        integer :: i
        integer :: k
        integer :: n

        n = size(b)
        x = b
        info = 0
        do k = 1, n
            if (piv(k) /= k) then
                tmp = x(k)
                x(k) = x(piv(k))
                x(piv(k)) = tmp
            end if
        end do
        do i = 2, n
            x(i) = x(i) - sum(lu(i, 1:i - 1) * x(1:i - 1))
        end do
        do i = n, 1, -1
            if (lu(i, i) == 0.0_dp) then
                info = i
                return
            end if
            if (i < n) x(i) = x(i) - sum(lu(i, i + 1:n) * x(i + 1:n))
            x(i) = x(i) / lu(i, i)
        end do
    end subroutine lu_solve_real

    subroutine lu_solve_complex(lu, piv, b, x, info)
        complex(dp), intent(in) :: lu(:, :) !! compact LU factors
        integer, intent(in) :: piv(:) !! row pivots from factorization
        complex(dp), intent(in) :: b(:) !! right-hand side
        complex(dp), intent(out) :: x(:) !! solution of the factored linear system
        integer, intent(out) :: info !! zero on success; nonzero for singular diagonal
        complex(dp) :: tmp
        integer :: i
        integer :: k
        integer :: n

        n = size(b)
        x = b
        info = 0
        do k = 1, n
            if (piv(k) /= k) then
                tmp = x(k)
                x(k) = x(piv(k))
                x(piv(k)) = tmp
            end if
        end do
        do i = 2, n
            x(i) = x(i) - sum(lu(i, 1:i - 1) * x(1:i - 1))
        end do
        do i = n, 1, -1
            if (lu(i, i) == cmplx(0.0_dp, 0.0_dp, kind=dp)) then
                info = i
                return
            end if
            if (i < n) x(i) = x(i) - sum(lu(i, i + 1:n) * x(i + 1:n))
            x(i) = x(i) / lu(i, i)
        end do
    end subroutine lu_solve_complex

    pure real(dp) function rms_norm(x) result(value)
        real(dp), intent(in) :: x(:) !! scaled vector whose RMS norm is requested
        value = sqrt(sum(x * x) / real(size(x), dp))
    end function rms_norm

    pure real(dp) function rms_norm_2d_scaled(x, scale) result(value)
        real(dp), intent(in) :: x(:, :) !! Newton correction matrix, shape (n,3)
        real(dp), intent(in) :: scale(:) !! component scales, shape (n)
        integer :: j

        value = 0.0_dp
        do j = 1, size(x, 2)
            value = value + sum((x(:, j) / scale)**2)
        end do
        value = sqrt(value / real(size(x), dp))
    end function rms_norm_2d_scaled

end module scifort_integrate_radau
