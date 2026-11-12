! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Complex-valued variable-order BDF/NDF stepper adapted from SciPy 1.17.0
! scipy/integrate/_ivp/bdf.py and its dense finite-difference Jacobian helper.
! The surrounding solve_ivp orchestration remains in scifort_integrate_ivp_complex.
! See CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.

module scifort_integrate_bdf_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_integrate_jacobian_fd, only : group_fd_columns
    use scifort_integrate_sparse_jacobian, only : ivp_complex_sparse_matrix, &
        ivp_csc_from_pattern
    use scifort_integrate_sparse_lu, only : ivp_sparse_lu_complex, sparse_lu_factor
    implicit none
    private

    integer, parameter, public :: bdf_complex_status_success = 0
    integer, parameter, public :: bdf_complex_status_too_small_step = 1
    integer, parameter, public :: bdf_complex_status_singular = 2
    integer, parameter, public :: bdf_complex_status_nonfinite = 3
    integer, parameter, public :: bdf_complex_status_max_attempts = 4
    integer, parameter, public :: bdf_complex_status_invalid_jacobian = 5
    integer, parameter, public :: bdf_complex_max_order = 5

    integer, parameter :: newton_maxiter = 4
    real(dp), parameter :: min_factor = 0.2_dp
    real(dp), parameter :: max_factor = 10.0_dp
    real(dp), parameter :: num_jac_diff_reject = 2.0097183471152322e-14_dp
    real(dp), parameter :: num_jac_diff_small = 1.8189894035458565e-12_dp
    real(dp), parameter :: num_jac_diff_big = 1.220703125e-4_dp
    real(dp), parameter :: num_jac_factor_increase = 10.0_dp
    real(dp), parameter :: num_jac_factor_decrease = 0.1_dp
    real(dp), parameter :: num_jac_min_factor = epsilon(1.0_dp) * 1000.0_dp

    abstract interface
        subroutine bdf_complex_rhs(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            complex(dp), intent(in) :: y(:) !! complex state at t
            complex(dp), intent(out) :: dydt(:) !! complex derivative dy/dt at t
        end subroutine bdf_complex_rhs

        subroutine bdf_complex_rhs_vectorized(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable shared by all state columns
            complex(dp), intent(in) :: y(:, :) !! complex state columns, shape (n,k)
            complex(dp), intent(out) :: dydt(:, :) !! derivatives for each state column
        end subroutine bdf_complex_rhs_vectorized

        subroutine bdf_complex_jacobian(t, y, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            complex(dp), intent(in) :: y(:) !! complex state at t
            complex(dp), intent(out) :: jac(:, :) !! complex Jacobian d fun / d y, shape (n,n)
        end subroutine bdf_complex_jacobian

        subroutine bdf_complex_sparse_jacobian(t, y, jac)
            import :: dp, ivp_complex_sparse_matrix
            implicit none
            real(dp), intent(in) :: t !! independent variable
            complex(dp), intent(in) :: y(:) !! complex state at t
            type(ivp_complex_sparse_matrix), intent(out) :: &
                jac !! one-based CSC complex Jacobian d fun / d y
        end subroutine bdf_complex_sparse_jacobian
    end interface

    type, public :: bdf_complex_state
        complex(dp), allocatable :: d(:, :)
        complex(dp), allocatable :: jac(:, :)
        real(dp), allocatable :: jac_factor(:)
        complex(dp), allocatable :: lu(:, :)
        type(ivp_complex_sparse_matrix) :: sparse_jac
        type(ivp_sparse_lu_complex) :: sparse_lu
        integer, allocatable :: piv(:)
        real(dp) :: h_abs = 0.0_dp
        real(dp) :: newton_tol = 0.0_dp
        integer :: order = 1
        integer :: n_equal_steps = 0
        logical :: lu_valid = .false.
        logical :: jacobian_constant = .false.
        logical :: jacobian_sparse = .false.
    end type bdf_complex_state

    public :: bdf_complex_advance
    public :: bdf_complex_dense_state
    public :: bdf_complex_initialize
    public :: bdf_complex_jacobian
    public :: bdf_complex_sparse_jacobian

contains

    subroutine bdf_complex_initialize(state, fun, t, y, f, direction, h_abs, rtol, atol, &
            njev, jac, atol_vector, jac_matrix, jac_sparsity, fun_vectorized, &
            jac_sparse_matrix, jac_sparse, rtol_vector)
        type(bdf_complex_state), intent(out) :: state !! initialized BDF recurrence and Jacobian state
        procedure(bdf_complex_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! initial independent-variable value
        complex(dp), intent(in) :: y(:) !! initial complex state
        complex(dp), intent(in) :: f(:) !! complex derivative at the initial state
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: h_abs !! selected initial step magnitude
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol !! scalar absolute tolerance
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        procedure(bdf_complex_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        complex(dp), intent(in), optional :: jac_matrix(:, :) !! constant dense Jacobian
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_complex_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for finite-difference Jacobians
        type(ivp_complex_sparse_matrix), intent(in), optional :: &
            jac_sparse_matrix !! constant complex CSC Jacobian
        procedure(bdf_complex_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent complex CSC Jacobian callback
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        real(dp) :: atol_use(size(y))
        real(dp) :: rtol_use(size(y))
        integer :: n
        logical :: use_sparse

        n = size(y)
        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector
        allocate(state%d(n, bdf_complex_max_order + 3), state%jac_factor(n))
        use_sparse = present(jac_sparse_matrix) .or. present(jac_sparse) .or. &
            (present(jac_sparsity) .and. .not. present(jac) .and. .not. present(jac_matrix))
        if (use_sparse) then
            allocate(state%lu(0, 0), state%piv(0))
        else
            allocate(state%jac(n, n), state%lu(n, n), state%piv(n))
        end if
        state%d = 0.0_dp
        state%d(:, 1) = y
        state%d(:, 2) = f * h_abs * direction
        state%h_abs = h_abs
        state%newton_tol = max(10.0_dp * epsilon(1.0_dp) / maxval(rtol_use), &
            min(0.03_dp, sqrt(maxval(rtol_use))))
        state%order = 1
        state%n_equal_steps = 0
        state%lu_valid = .false.
        state%jacobian_sparse = use_sparse
        state%jacobian_constant = present(jac_matrix) .or. present(jac_sparse_matrix)
        state%jac_factor = sqrt(epsilon(1.0_dp))
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
    end subroutine bdf_complex_initialize

    subroutine bdf_complex_advance(state, fun, t, t_bound, direction, rtol, atol, max_step, &
            y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, nfev, njev, &
            nlu, attempt_limit, attempts, rejected, status, jac, atol_vector, jac_sparsity, &
            fun_vectorized, jac_sparse, rtol_vector)
        type(bdf_complex_state), intent(inout) :: state !! mutable BDF recurrence/Jacobian state
        procedure(bdf_complex_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: t_bound !! final integration bound
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol !! scalar absolute tolerance
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        complex(dp), intent(out) :: y_new(:) !! accepted complex state
        real(dp), intent(out) :: t_new !! accepted independent-variable value
        complex(dp), intent(out) :: dense_d(:, :) !! complex BDF dense difference coefficients
        real(dp), intent(out) :: dense_h !! signed BDF dense-output step spacing
        integer, intent(out) :: dense_order !! polynomial order for dense output
        real(dp), intent(out) :: dense_anchor !! full accepted endpoint used by dense output
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        integer, intent(inout) :: nlu !! accumulated LU factorization count
        integer, intent(in) :: attempt_limit !! maximum trial steps available to this call
        integer, intent(out) :: attempts !! trial steps consumed by this accepted step
        integer, intent(out) :: rejected !! rejected trial-step count
        integer, intent(out) :: status !! BDF status code
        procedure(bdf_complex_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_complex_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for finite-difference Jacobians
        procedure(bdf_complex_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent complex CSC Jacobian callback
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        real(dp), parameter :: kappa(0:bdf_complex_max_order) = [0.0_dp, -0.1850_dp, &
            -1.0_dp / 9.0_dp, -0.0823_dp, -0.0415_dp, 0.0_dp]
        real(dp) :: gamma(0:bdf_complex_max_order)
        real(dp) :: alpha(0:bdf_complex_max_order)
        real(dp) :: error_const(0:bdf_complex_max_order)
        real(dp) :: h_abs
        real(dp) :: h
        real(dp) :: min_step
        complex(dp) :: y_predict(size(y_new))
        real(dp) :: scale(size(y_new))
        complex(dp) :: psi(size(y_new))
        complex(dp) :: d_corr(size(y_new))
        complex(dp) :: error(size(y_new))
        complex(dp) :: error_m(size(y_new))
        complex(dp) :: error_p(size(y_new))
        real(dp) :: atol_use(size(y_new))
        real(dp) :: rtol_use(size(y_new))
        real(dp) :: factor
        real(dp) :: safety_newton
        real(dp) :: error_norm
        real(dp) :: error_m_norm
        real(dp) :: error_p_norm
        real(dp) :: factors(3)
        real(dp) :: c
        integer :: order
        integer :: i
        integer :: delta_order
        integer :: n_iter
        integer :: info
        logical :: converged
        logical :: accepted
        logical :: current_jac

        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector

        status = bdf_complex_status_success
        attempts = 0
        rejected = 0
        dense_d = 0.0_dp
        dense_h = 0.0_dp
        dense_order = 0
        dense_anchor = t

        gamma = 0.0_dp
        do i = 1, bdf_complex_max_order
            gamma(i) = gamma(i - 1) + 1.0_dp / real(i, dp)
        end do
        alpha = (1.0_dp - kappa) * gamma
        do i = 0, bdf_complex_max_order
            error_const(i) = kappa(i) * gamma(i) + 1.0_dp / real(i + 1, dp)
        end do

        h_abs = state%h_abs
        order = state%order
        min_step = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(t))
        if (h_abs > max_step) then
            factor = max_step / h_abs
            h_abs = max_step
            call change_d(state%d, order, factor)
            state%n_equal_steps = 0
        else if (h_abs < min_step) then
            factor = min_step / h_abs
            h_abs = min_step
            call change_d(state%d, order, factor)
            state%n_equal_steps = 0
        end if

        current_jac = state%jacobian_constant
        accepted = .false.
        do while (.not. accepted)
            if (attempts >= attempt_limit) then
                status = bdf_complex_status_max_attempts
                return
            end if
            attempts = attempts + 1
            if (h_abs < min_step) then
                status = bdf_complex_status_too_small_step
                return
            end if

            h = h_abs * direction
            t_new = t + h
            if (direction * (t_new - t_bound) > 0.0_dp) then
                t_new = t_bound
                factor = abs(t_new - t) / h_abs
                call change_d(state%d, order, factor)
                state%n_equal_steps = 0
                state%lu_valid = .false.
            end if
            h = t_new - t
            h_abs = abs(h)

            y_predict = sum(state%d(:, 1:order + 1), dim=2)
            scale = max(atol_use + rtol_use * abs(y_predict), tiny(1.0_dp))
            psi = 0.0_dp
            do i = 1, order
                psi = psi + gamma(i) * state%d(:, i + 1)
            end do
            psi = psi / alpha(order)
            c = h / alpha(order)

            converged = .false.
            do while (.not. converged)
                if (.not. state%lu_valid) then
                    if (state%jacobian_sparse) then
                        call sparse_lu_factor(state%sparse_jac, &
                            cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(-c, 0.0_dp, kind=dp), &
                            state%sparse_lu, info)
                    else
                        call factor_newton_matrix(state%jac, c, state%lu, state%piv, info)
                    end if
                    nlu = nlu + 1
                    if (info /= 0) then
                        status = bdf_complex_status_singular
                        return
                    end if
                    state%lu_valid = .true.
                end if

                call solve_bdf_complex_system(fun, t_new, y_predict, c, psi, state%lu, &
                    state%piv, scale, state%newton_tol, y_new, d_corr, converged, &
                    n_iter, nfev, status, sparse_lu=state%sparse_lu)
                if (status /= bdf_complex_status_success) return
                if (converged) exit
                if (current_jac) exit

                if (present(jac_sparse)) then
                    call jac_sparse(t_new, y_predict, state%sparse_jac)
                    if (.not. sparse_jacobian_is_valid(state%sparse_jac, size(y_predict))) then
                        status = bdf_complex_status_invalid_jacobian
                        return
                    end if
                else if (present(jac)) then
                    call jac(t_new, y_predict, state%jac)
                else if (state%jacobian_sparse) then
                    call numerical_sparse_jacobian_at_point(fun, t_new, y_predict, &
                        atol_use, state%jac_factor, state%sparse_jac, jac_sparsity, &
                        fun_vectorized)
                else
                    call numerical_jacobian_at_point(fun, t_new, y_predict, atol_use, &
                        state%jac_factor, state%jac, jac_sparsity, fun_vectorized)
                end if
                njev = njev + 1
                if (.not. state%jacobian_sparse) then
                    if (.not. all(complex_is_finite(state%jac))) then
                        status = bdf_complex_status_nonfinite
                        return
                    end if
                end if
                state%lu_valid = .false.
                current_jac = .true.
            end do

            if (.not. converged) then
                factor = 0.5_dp
                h_abs = h_abs * factor
                call change_d(state%d, order, factor)
                state%n_equal_steps = 0
                state%lu_valid = .false.
                rejected = rejected + 1
                cycle
            end if

            safety_newton = 0.9_dp * real(2 * newton_maxiter + 1, dp) / &
                real(2 * newton_maxiter + n_iter, dp)
            scale = max(atol_use + rtol_use * abs(y_new), tiny(1.0_dp))
            error = error_const(order) * d_corr
            error_norm = rms_norm(error / scale)
            if (error_norm > 1.0_dp) then
                factor = max(min_factor, &
                    safety_newton * error_norm**(-1.0_dp / real(order + 1, dp)))
                h_abs = h_abs * factor
                call change_d(state%d, order, factor)
                state%n_equal_steps = 0
                rejected = rejected + 1
            else
                accepted = .true.
            end if
        end do

        state%n_equal_steps = state%n_equal_steps + 1
        state%h_abs = h_abs

        state%d(:, order + 3) = d_corr - state%d(:, order + 2)
        state%d(:, order + 2) = d_corr
        do i = order + 1, 1, -1
            state%d(:, i) = state%d(:, i) + state%d(:, i + 1)
        end do

        if (state%n_equal_steps >= order + 1) then
            if (order > 1) then
                error_m = error_const(order - 1) * state%d(:, order + 1)
                error_m_norm = rms_norm(error_m / scale)
            else
                error_m_norm = huge(1.0_dp)
            end if
            if (order < bdf_complex_max_order) then
                error_p = error_const(order + 1) * state%d(:, order + 3)
                error_p_norm = rms_norm(error_p / scale)
            else
                error_p_norm = huge(1.0_dp)
            end if

            factors(1) = error_factor(error_m_norm, order)
            factors(2) = error_factor(error_norm, order + 1)
            factors(3) = error_factor(error_p_norm, order + 2)
            delta_order = maxloc(factors, dim=1) - 2
            order = order + delta_order
            state%order = order
            factor = min(max_factor, safety_newton * maxval(factors))
            state%h_abs = state%h_abs * factor
            call change_d(state%d, order, factor)
            state%n_equal_steps = 0
            state%lu_valid = .false.
        else
            state%order = order
        end if

        dense_order = state%order
        dense_h = state%h_abs * direction
        dense_anchor = t_new
        dense_d = 0.0_dp
        dense_d(:, 1:dense_order + 1) = state%d(:, 1:dense_order + 1)
    end subroutine bdf_complex_advance

    subroutine solve_bdf_complex_system(fun, t_new, y_predict, c, psi, lu, piv, scale, tol, &
            y, d, converged, n_iter, nfev, status, sparse_lu)
        procedure(bdf_complex_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t_new !! trial-step endpoint
        complex(dp), intent(in) :: y_predict(:) !! complex extrapolated predictor at t_new
        real(dp), intent(in) :: c !! BDF nonlinear-system scale h/alpha(order)
        complex(dp), intent(in) :: psi(:) !! complex BDF accumulated backward-difference term
        complex(dp), intent(in) :: lu(:, :) !! complex LU factors of I-cJ
        integer, intent(in) :: piv(:) !! row pivots paired with lu
        real(dp), intent(in) :: scale(:) !! component error scales
        real(dp), intent(in) :: tol !! Newton convergence tolerance
        complex(dp), intent(out) :: y(:) !! final complex Newton iterate
        complex(dp), intent(out) :: d(:) !! accumulated complex Newton correction
        logical, intent(out) :: converged !! true when Newton stopping criterion is satisfied
        integer, intent(out) :: n_iter !! Newton iterations performed
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(out) :: status !! BDF status code
        type(ivp_sparse_lu_complex), intent(in), optional :: &
            sparse_lu !! sparse LU factors used when a constant CSC Jacobian is active

        complex(dp) :: f(size(y))
        complex(dp) :: rhs(size(y))
        complex(dp) :: dy(size(y))
        real(dp) :: dy_norm
        real(dp) :: dy_norm_old
        real(dp) :: rate
        integer :: k
        integer :: info
        logical :: have_old_norm

        y = y_predict
        d = 0.0_dp
        converged = .false.
        have_old_norm = .false.
        dy_norm_old = 0.0_dp
        status = bdf_complex_status_success
        n_iter = 0

        do k = 1, newton_maxiter
            n_iter = k
            call fun(t_new, y, f)
            nfev = nfev + 1
            if (.not. all(complex_is_finite(f))) then
                status = bdf_complex_status_nonfinite
                return
            end if
            rhs = c * f - psi - d
            if (present(sparse_lu) .and. sparse_lu%n > 0) then
                call sparse_lu%solve(rhs, dy, info)
            else
                call lu_solve(lu, piv, rhs, dy, info)
            end if
            if (info /= 0) then
                status = bdf_complex_status_singular
                return
            end if
            dy_norm = rms_norm(dy / scale)

            if (have_old_norm) then
                if (dy_norm_old == 0.0_dp) then
                    rate = 0.0_dp
                else
                    rate = dy_norm / dy_norm_old
                end if
                if (rate >= 1.0_dp) exit
                if (rate**real(newton_maxiter - k, dp) / (1.0_dp - rate) * &
                        dy_norm > tol) exit
            else
                rate = 0.0_dp
            end if

            y = y + dy
            d = d + dy
            if (dy_norm == 0.0_dp) then
                converged = .true.
                exit
            end if
            if (have_old_norm) then
                if (rate / (1.0_dp - rate) * dy_norm < tol) then
                    converged = .true.
                    exit
                end if
            end if
            dy_norm_old = dy_norm
            have_old_norm = .true.
        end do
    end subroutine solve_bdf_complex_system

    subroutine bdf_complex_dense_state(t_anchor, h, order, d, tq, yq)
        real(dp), intent(in) :: t_anchor !! full accepted BDF endpoint anchoring the polynomial
        real(dp), intent(in) :: h !! signed BDF difference spacing used by the dense polynomial
        integer, intent(in) :: order !! BDF dense polynomial order, 1 through 5
        complex(dp), intent(in) :: d(:, :) !! complex backward-difference coefficients
        real(dp), intent(in) :: tq !! interpolation point in the accepted interval
        complex(dp), intent(out) :: yq(:) !! interpolated complex state

        real(dp) :: p
        real(dp) :: x
        integer :: j

        yq = d(:, 1)
        if (order <= 0 .or. h == 0.0_dp) return
        p = 1.0_dp
        do j = 1, order
            x = (tq - (t_anchor - h * real(j - 1, dp))) / (h * real(j, dp))
            p = p * x
            yq = yq + d(:, j + 1) * p
        end do
    end subroutine bdf_complex_dense_state

    subroutine numerical_jacobian(fun, t, y, f, atol, factor, jac, jac_sparsity, &
            fun_vectorized)
        procedure(bdf_complex_rhs) :: fun !! scalar complex RHS used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        complex(dp), intent(in) :: y(:) !! complex Jacobian evaluation state
        complex(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        complex(dp), intent(out) :: jac(:, :) !! finite-difference complex Jacobian
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_complex_rhs_vectorized), optional :: &
            fun_vectorized !! batched complex RHS for Jacobian differences

        complex(dp) :: f_batch(size(y), size(y))
        complex(dp) :: f_new(size(y), size(y))
        complex(dp) :: f_retry(size(y))
        complex(dp) :: h(size(y))
        complex(dp) :: h_retry
        complex(dp) :: y_batch(size(y), size(y))
        complex(dp) :: y_trial(size(y))
        real(dp) :: max_diff
        real(dp) :: max_diff_new
        real(dp) :: new_factor
        real(dp) :: scale_new
        real(dp) :: scale_ref
        real(dp) :: y_scale(size(y))
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
            if (real(f(j), dp) >= 0.0_dp) then
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

        jac = cmplx(0.0_dp, 0.0_dp, dp)
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
                    jac(:, j) = cmplx(0.0_dp, 0.0_dp, dp)
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

    subroutine numerical_jacobian_at_point(fun, t, y, atol, factor, jac, jac_sparsity, &
            fun_vectorized)
        procedure(bdf_complex_rhs) :: fun !! scalar complex RHS used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        complex(dp), intent(in) :: y(:) !! complex Jacobian evaluation state
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        complex(dp), intent(out) :: jac(:, :) !! finite-difference complex Jacobian
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_complex_rhs_vectorized), optional :: &
            fun_vectorized !! batched complex RHS for Jacobian differences

        complex(dp) :: f(size(y))
        complex(dp) :: f_batch(size(y), 1)
        complex(dp) :: y_batch(size(y), 1)

        if (present(fun_vectorized)) then
            y_batch(:, 1) = y
            call fun_vectorized(t, y_batch, f_batch)
            f = f_batch(:, 1)
        else
            call fun(t, y, f)
        end if
        call numerical_jacobian(fun, t, y, f, atol, factor, jac, jac_sparsity, &
            fun_vectorized)
    end subroutine numerical_jacobian_at_point


    subroutine numerical_sparse_jacobian(fun, t, y, f, atol, factor, jac, sparsity, &
            fun_vectorized)
        procedure(bdf_complex_rhs) :: fun !! scalar complex RHS used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        complex(dp), intent(in) :: y(:) !! complex Jacobian evaluation state
        complex(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        type(ivp_complex_sparse_matrix), intent(inout) :: &
            jac !! CSC complex Jacobian with fixed sparsity structure
        logical, intent(in) :: sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_complex_rhs_vectorized), optional :: &
            fun_vectorized !! batched complex RHS for Jacobian differences

        complex(dp), allocatable :: f_group(:, :)
        complex(dp), allocatable :: y_batch(:, :)
        complex(dp) :: f_col(size(y))
        complex(dp) :: f_retry(size(y))
        complex(dp) :: h(size(y))
        complex(dp) :: h_retry
        complex(dp) :: y_trial(size(y))
        real(dp) :: max_diff
        real(dp) :: max_diff_new
        real(dp) :: new_factor
        real(dp) :: scale_new
        real(dp) :: scale_ref
        real(dp) :: y_scale(size(y))
        integer :: g
        integer :: groups(size(y))
        integer :: i
        integer :: imax
        integer :: j
        integer :: k
        integer :: n_groups
        logical :: use_retry

        do j = 1, size(y)
            if (real(f(j), dp) >= 0.0_dp) then
                y_scale(j) = max(atol(j), abs(y(j)))
            else
                y_scale(j) = -max(atol(j), abs(y(j)))
            end if
            h(j) = (y(j) + factor(j) * y_scale(j)) - y(j)
            do while (h(j) == cmplx(0.0_dp, 0.0_dp, kind=dp))
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

        jac%data = cmplx(0.0_dp, 0.0_dp, kind=dp)
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

    subroutine numerical_sparse_jacobian_at_point(fun, t, y, atol, factor, jac, sparsity, &
            fun_vectorized)
        procedure(bdf_complex_rhs) :: fun !! scalar complex RHS used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        complex(dp), intent(in) :: y(:) !! complex Jacobian evaluation state
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        type(ivp_complex_sparse_matrix), intent(inout) :: &
            jac !! CSC complex Jacobian with fixed sparsity structure
        logical, intent(in) :: sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_complex_rhs_vectorized), optional :: &
            fun_vectorized !! batched complex RHS for Jacobian differences

        complex(dp) :: f(size(y))
        complex(dp) :: f_batch(size(y), 1)
        complex(dp) :: y_batch(size(y), 1)

        if (present(fun_vectorized)) then
            y_batch(:, 1) = y
            call fun_vectorized(t, y_batch, f_batch)
            f = f_batch(:, 1)
        else
            call fun(t, y, f)
        end if
        call numerical_sparse_jacobian(fun, t, y, f, atol, factor, jac, sparsity, &
            fun_vectorized)
    end subroutine numerical_sparse_jacobian_at_point

    logical function sparse_jacobian_is_valid(jac, n) result(valid)
        type(ivp_complex_sparse_matrix), intent(in) :: jac !! sparse Jacobian returned by callback
        integer, intent(in) :: n !! required square matrix order

        valid = jac%is_valid()
        if (.not. valid) return
        valid = jac%nrows == n .and. jac%ncols == n
    end function sparse_jacobian_is_valid


    subroutine factor_newton_matrix(jac, c, lu, piv, info)
        complex(dp), intent(in) :: jac(:, :) !! complex Jacobian matrix
        real(dp), intent(in) :: c !! BDF nonlinear-system scale h/alpha(order)
        complex(dp), intent(out) :: lu(:, :) !! complex LU factors of I-cJ
        integer, intent(out) :: piv(:) !! row pivots for LU solve
        integer, intent(out) :: info !! zero on success; nonzero if numerically singular

        integer :: i

        lu = -c * jac
        do i = 1, size(jac, 1)
            lu(i, i) = lu(i, i) + 1.0_dp
        end do
        call lu_factor(lu, piv, info)
    end subroutine factor_newton_matrix

    subroutine lu_factor(a, piv, info)
        complex(dp), intent(inout) :: a(:, :) !! complex square matrix, overwritten by compact LU factors
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
    end subroutine lu_factor

    subroutine lu_solve(lu, piv, b, x, info)
        complex(dp), intent(in) :: lu(:, :) !! compact complex LU factors
        integer, intent(in) :: piv(:) !! row pivots from lu_factor
        complex(dp), intent(in) :: b(:) !! complex right-hand side
        complex(dp), intent(out) :: x(:) !! complex linear-system solution
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
            if (lu(i, i) == 0.0_dp) then
                info = i
                return
            end if
            if (i < n) x(i) = x(i) - sum(lu(i, i + 1:n) * x(i + 1:n))
            x(i) = x(i) / lu(i, i)
        end do
    end subroutine lu_solve

    subroutine change_d(d, order, factor)
        complex(dp), intent(inout) :: d(:, :) !! complex backward-difference table to transform
        integer, intent(in) :: order !! current BDF order
        real(dp), intent(in) :: factor !! new step size divided by old step size

        real(dp) :: r(0:bdf_complex_max_order, 0:bdf_complex_max_order)
        real(dp) :: u(0:bdf_complex_max_order, 0:bdf_complex_max_order)
        real(dp) :: ru(0:bdf_complex_max_order, 0:bdf_complex_max_order)
        real(dp) :: m
        integer :: i
        integer :: j

        r = 0.0_dp
        u = 0.0_dp
        r(0, 0:order) = 1.0_dp
        u(0, 0:order) = 1.0_dp
        do j = 1, order
            do i = 1, order
                m = (real(i - 1, dp) - factor * real(j, dp)) / real(i, dp)
                r(i, j) = r(i - 1, j) * m
                m = (real(i - 1, dp) - real(j, dp)) / real(i, dp)
                u(i, j) = u(i - 1, j) * m
            end do
        end do
        ru(0:order, 0:order) = matmul(r(0:order, 0:order), u(0:order, 0:order))
        d(:, 1:order + 1) = matmul(d(:, 1:order + 1), &
            cmplx(ru(0:order, 0:order), 0.0_dp, kind=dp))
    end subroutine change_d

    pure real(dp) function rms_norm(x) result(value)
        complex(dp), intent(in) :: x(:) !! scaled complex vector whose RMS norm is requested

        value = sqrt(sum(abs(x)**2) / real(size(x), dp))
    end function rms_norm

    elemental logical function complex_is_finite(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value to test for finite components

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure real(dp) function error_factor(error_norm, exponent_denominator) result(value)
        real(dp), intent(in) :: error_norm !! normalized local-error estimate
        integer, intent(in) :: exponent_denominator !! positive denominator for adaptation exponent

        if (error_norm == 0.0_dp) then
            value = huge(1.0_dp)
        else if (.not. ieee_is_finite(error_norm)) then
            value = 0.0_dp
        else
            value = error_norm**(-1.0_dp / real(exponent_denominator, dp))
        end if
    end function error_factor

end module scifort_integrate_bdf_complex
