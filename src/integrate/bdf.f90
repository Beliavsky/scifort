! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Real-valued variable-order BDF/NDF stepper adapted from SciPy 1.17.0
! scipy/integrate/_ivp/bdf.py and its dense finite-difference Jacobian helper.
! The surrounding solve_ivp orchestration remains in scifort_integrate_ivp.
! An LSODA-only state flag enables ODEPACK-style pure-BDF Nordsieck stepping,
! chord convergence/order control, accepted-state weighted maximum norms, Jacobian
! refresh, and numerical-Jacobian accounting; standalone BDF leaves it disabled.
! See CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.

module scifort_integrate_bdf
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_integrate_jacobian_fd, only : group_fd_columns
    use scifort_integrate_sparse_jacobian, only : ivp_csc_from_pattern, ivp_sparse_matrix
    use scifort_integrate_sparse_lu, only : ivp_sparse_lu_real, sparse_lu_factor
    implicit none
    private

    integer, parameter, public :: bdf_status_success = 0
    integer, parameter, public :: bdf_status_too_small_step = 1
    integer, parameter, public :: bdf_status_singular = 2
    integer, parameter, public :: bdf_status_nonfinite = 3
    integer, parameter, public :: bdf_status_max_attempts = 4
    integer, parameter, public :: bdf_status_invalid_jacobian = 5
    integer, parameter, public :: bdf_max_order = 5

    integer, parameter :: newton_maxiter = 4
    integer, parameter :: odepack_jacobian_refresh_steps = 20
    real(dp), parameter :: min_factor = 0.2_dp
    real(dp), parameter :: max_factor = 10.0_dp
    real(dp), parameter :: odepack_jacobian_refresh_ratio = 0.3_dp
    real(dp), parameter :: num_jac_diff_reject = 2.0097183471152322e-14_dp
    real(dp), parameter :: num_jac_diff_small = 1.8189894035458565e-12_dp
    real(dp), parameter :: num_jac_diff_big = 1.220703125e-4_dp
    real(dp), parameter :: num_jac_factor_increase = 10.0_dp
    real(dp), parameter :: num_jac_factor_decrease = 0.1_dp
    real(dp), parameter :: num_jac_min_factor = epsilon(1.0_dp) * 1000.0_dp

    abstract interface
        subroutine bdf_rhs(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t
        end subroutine bdf_rhs

        subroutine bdf_rhs_vectorized(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable shared by all state columns
            real(dp), intent(in) :: y(:, :) !! state columns, shape (n,k)
            real(dp), intent(out) :: dydt(:, :) !! derivatives for each state column
        end subroutine bdf_rhs_vectorized

        subroutine bdf_jacobian(t, y, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: jac(:, :) !! Jacobian d fun / d y, shape (n,n)
        end subroutine bdf_jacobian

        subroutine bdf_sparse_jacobian(t, y, jac)
            import :: dp, ivp_sparse_matrix
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            type(ivp_sparse_matrix), intent(out) :: jac !! one-based CSC Jacobian d fun / d y
        end subroutine bdf_sparse_jacobian
    end interface

    type, public :: bdf_state
        real(dp), allocatable :: d(:, :)
        real(dp), allocatable :: jac(:, :)
        real(dp), allocatable :: jac_factor(:)
        real(dp), allocatable :: lu(:, :)
        type(ivp_sparse_matrix) :: sparse_jac
        type(ivp_sparse_lu_real) :: sparse_lu
        real(dp), allocatable :: odepack_yh(:, :)
        real(dp), allocatable :: odepack_saved_acor(:)
        integer, allocatable :: piv(:)
        real(dp) :: h_abs = 0.0_dp
        real(dp) :: newton_tol = 0.0_dp
        real(dp) :: last_error_norm = 0.0_dp
        real(dp) :: last_step_h_abs = 0.0_dp
        real(dp) :: jacobian_c = 0.0_dp
        real(dp) :: odepack_crate = 0.7_dp
        real(dp) :: odepack_rmax = 10.0_dp
        integer :: order = 1
        integer :: max_order = bdf_max_order
        integer :: odepack_ialth = 2
        integer :: last_step_order = 1
        integer :: n_equal_steps = 0
        integer :: jacobian_age = 0
        logical :: lu_valid = .false.
        logical :: jacobian_c_valid = .false.
        logical :: jacobian_initialized = .false.
        logical :: jacobian_constant = .false.
        logical :: jacobian_sparse = .false.
        logical :: odepack_jacobian_refresh = .false.
        logical :: odepack_saved_acor_valid = .false.
        logical :: odepack_yh_valid = .false.
    end type bdf_state

    public :: bdf_advance
    public :: bdf_dense_state
    public :: bdf_initialize
    public :: bdf_jacobian
    public :: bdf_sparse_jacobian

contains

    subroutine bdf_initialize(state, fun, t, y, f, direction, h_abs, rtol, atol, &
            njev, jac, lband, uband, defer_jacobian, atol_vector, jac_matrix, &
            jac_sparsity, fun_vectorized, jac_sparse_matrix, jac_sparse, rtol_vector, &
            max_order)
        type(bdf_state), intent(out) :: state !! initialized BDF recurrence and Jacobian state
        procedure(bdf_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! initial independent-variable value
        real(dp), intent(in) :: y(:) !! initial state
        real(dp), intent(in) :: f(:) !! derivative at the initial state
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: h_abs !! selected initial step magnitude
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol !! scalar absolute tolerance
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        procedure(bdf_jacobian), optional :: jac !! optional analytic Jacobian callback
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense
        logical, intent(in), optional :: defer_jacobian !! postpone initial Jacobian evaluation
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        real(dp), intent(in), optional :: jac_matrix(:, :) !! constant dense Jacobian, shape (n,n)
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for finite-difference Jacobians
        type(ivp_sparse_matrix), intent(in), optional :: &
            jac_sparse_matrix !! constant CSC Jacobian for ordinary BDF
        procedure(bdf_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent CSC Jacobian callback for ordinary BDF
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances
        integer, intent(in), optional :: max_order !! maximum BDF order, 1 through 5; default 5

        integer :: n
        real(dp) :: atol_use(size(y))
        real(dp) :: rtol_use(size(y))
        logical :: defer_jac
        logical :: use_sparse

        n = size(y)
        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector
        allocate(state%d(n, bdf_max_order + 3), state%jac_factor(n), &
            state%odepack_yh(n, bdf_max_order + 1), state%odepack_saved_acor(n))
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
        state%odepack_yh = 0.0_dp
        state%odepack_yh(:, 1) = y
        state%odepack_yh(:, 2) = f * h_abs * direction
        state%odepack_saved_acor = 0.0_dp
        state%h_abs = h_abs
        state%newton_tol = max(10.0_dp * epsilon(1.0_dp) / maxval(rtol_use), &
            min(0.03_dp, sqrt(maxval(rtol_use))))
        state%last_error_norm = 0.0_dp
        state%last_step_h_abs = 0.0_dp
        state%jacobian_c = 0.0_dp
        state%odepack_crate = 0.7_dp
        state%odepack_rmax = 10.0_dp
        state%order = 1
        state%max_order = bdf_max_order
        if (present(max_order)) state%max_order = min(max(1, max_order), bdf_max_order)
        state%odepack_ialth = 2
        state%last_step_order = 1
        state%n_equal_steps = 0
        state%jacobian_age = 0
        state%lu_valid = .false.
        state%jacobian_c_valid = .false.
        state%jacobian_initialized = .false.
        state%jacobian_sparse = use_sparse
        state%jacobian_constant = present(jac_matrix) .or. present(jac_sparse_matrix)
        state%odepack_jacobian_refresh = .false.
        state%odepack_saved_acor_valid = .false.
        state%odepack_yh_valid = .true.
        state%jac_factor = sqrt(epsilon(1.0_dp))
        defer_jac = .false.
        if (present(defer_jacobian)) defer_jac = defer_jacobian
        if (present(jac_sparse_matrix)) then
            state%sparse_jac = jac_sparse_matrix
            state%jacobian_initialized = .true.
        else if (present(jac_matrix)) then
            state%jac = jac_matrix
            state%jacobian_initialized = .true.
        else if (defer_jac) then
            state%jac = 0.0_dp
        else
            if (present(jac_sparse)) then
                call jac_sparse(t, y, state%sparse_jac)
            else if (present(jac)) then
                call evaluate_user_jacobian(jac, t, y, state%jac, lband, uband)
            else if (state%jacobian_sparse) then
                call ivp_csc_from_pattern(jac_sparsity, state%sparse_jac)
                call numerical_sparse_jacobian(fun, t, y, f, atol_use, state%jac_factor, &
                    state%sparse_jac, jac_sparsity, fun_vectorized)
            else
                call numerical_jacobian(fun, t, y, f, atol_use, state%jac_factor, state%jac, &
                    lband, uband, jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
            end if
            state%jacobian_initialized = .true.
            njev = njev + 1
        end if
    end subroutine bdf_initialize

    subroutine bdf_advance(state, fun, t, t_bound, direction, rtol, atol, max_step, &
            y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, nfev, njev, &
            nlu, attempt_limit, attempts, rejected, status, jac, min_step_user, lband, uband, &
            f_new_out, atol_vector, jac_sparsity, fun_vectorized, jac_sparse, rtol_vector)
        type(bdf_state), intent(inout) :: state !! mutable BDF recurrence/Jacobian state
        procedure(bdf_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: t_bound !! final integration bound
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol !! scalar absolute tolerance
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        real(dp), intent(out) :: y_new(:) !! accepted state
        real(dp), intent(out) :: t_new !! accepted independent-variable value
        real(dp), intent(out) :: dense_d(:, :) !! BDF dense difference coefficients
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
        procedure(bdf_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: min_step_user !! minimum retry step magnitude; default 0
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense
        real(dp), intent(out), optional :: f_new_out(:) !! derivative implied by the converged
            !! BDF equation
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for finite-difference Jacobians
        procedure(bdf_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent CSC Jacobian callback
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        real(dp), parameter :: kappa(0:bdf_max_order) = [0.0_dp, -0.1850_dp, &
            -1.0_dp / 9.0_dp, -0.0823_dp, -0.0415_dp, 0.0_dp]
        real(dp) :: gamma(0:bdf_max_order)
        real(dp) :: alpha(0:bdf_max_order)
        real(dp) :: error_const(0:bdf_max_order)
        real(dp) :: h_abs
        real(dp) :: h
        real(dp) :: min_step
        real(dp) :: retry_min_step
        real(dp) :: y_predict(size(y_new))
        real(dp) :: scale(size(y_new))
        real(dp) :: psi(size(y_new))
        real(dp) :: d_corr(size(y_new))
        real(dp) :: error(size(y_new))
        real(dp) :: error_m(size(y_new))
        real(dp) :: error_p(size(y_new))
        real(dp) :: factor
        real(dp) :: safety_newton
        real(dp) :: error_norm
        real(dp) :: error_m_norm
        real(dp) :: error_p_norm
        real(dp) :: factors(3)
        real(dp) :: c
        real(dp) :: jac_f(size(y_new))
        real(dp) :: atol_use(size(y_new))
        real(dp) :: rtol_use(size(y_new))
        integer :: order
        integer :: i
        integer :: delta_order
        integer :: n_iter
        integer :: info
        integer :: jac_nfev
        logical :: converged
        logical :: accepted
        logical :: current_jac
        logical :: refresh_jacobian
        logical :: have_jac_f

        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector

        if (state%odepack_jacobian_refresh) then
            call bdf_advance_odepack(state, fun, t, t_bound, direction, rtol, atol_use, &
                max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                nfev, njev, nlu, attempt_limit, attempts, rejected, status, jac, &
                min_step_user, lband, uband, f_new_out, rtol_vector=rtol_use)
            return
        end if

        status = bdf_status_success
        attempts = 0
        rejected = 0
        dense_d = 0.0_dp
        dense_h = 0.0_dp
        dense_order = 0
        dense_anchor = t

        gamma = 0.0_dp
        do i = 1, bdf_max_order
            gamma(i) = gamma(i - 1) + 1.0_dp / real(i, dp)
        end do
        if (state%odepack_jacobian_refresh) then
            alpha = gamma
            do i = 0, bdf_max_order
                error_const(i) = 1.0_dp / real(i + 1, dp)
            end do
        else
            alpha = (1.0_dp - kappa) * gamma
            do i = 0, bdf_max_order
                error_const(i) = kappa(i) * gamma(i) + 1.0_dp / real(i + 1, dp)
            end do
        end if

        h_abs = state%h_abs
        order = state%order
        min_step = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(t))
        retry_min_step = min_step
        if (present(min_step_user)) retry_min_step = max(retry_min_step, min_step_user)
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
                status = bdf_status_max_attempts
                return
            end if
            attempts = attempts + 1
            if (h_abs < min_step) then
                status = bdf_status_too_small_step
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
            if (state%odepack_jacobian_refresh) then
                scale = max(atol_use + rtol_use * abs(state%d(:, 1)), tiny(1.0_dp))
            else
                scale = max(atol_use + rtol_use * abs(y_predict), tiny(1.0_dp))
            end if
            psi = 0.0_dp
            do i = 1, order
                psi = psi + gamma(i) * state%d(:, i + 1)
            end do
            psi = psi / alpha(order)
            c = h / alpha(order)

            converged = .false.
            have_jac_f = .false.
            do while (.not. converged)
                refresh_jacobian = .not. state%jacobian_initialized
                if (state%odepack_jacobian_refresh .and. .not. current_jac .and. &
                        state%jacobian_initialized) then
                    refresh_jacobian = state%jacobian_age >= &
                        odepack_jacobian_refresh_steps
                    if (state%jacobian_c_valid) then
                        refresh_jacobian = refresh_jacobian .or. &
                            abs(c / state%jacobian_c - 1.0_dp) > &
                            odepack_jacobian_refresh_ratio
                    end if
                end if
                if (refresh_jacobian) then
                    if (present(jac_sparse)) then
                        call jac_sparse(t_new, y_predict, state%sparse_jac)
                        if (.not. sparse_jacobian_is_valid(state%sparse_jac, size(y_predict))) then
                            status = bdf_status_invalid_jacobian
                            return
                        end if
                    else if (present(jac)) then
                        call evaluate_user_jacobian(jac, t_new, y_predict, state%jac, &
                            lband, uband)
                    else if (state%odepack_jacobian_refresh) then
                        call fun(t_new, y_predict, jac_f)
                        nfev = nfev + 1
                        if (.not. all(ieee_is_finite(jac_f))) then
                            status = bdf_status_nonfinite
                            return
                        end if
                        jac_nfev = 0
                        call odepack_numerical_jacobian(fun, t_new, y_predict, &
                            state%d(:, 1), jac_f, h, rtol_use, atol_use, state%jac, lband, &
                            uband, jac_nfev)
                        nfev = nfev + jac_nfev
                        have_jac_f = .true.
                    else if (state%jacobian_sparse) then
                        call numerical_sparse_jacobian_at_point(fun, t_new, y_predict, &
                            atol_use, state%jac_factor, state%sparse_jac, jac_sparsity, &
                            fun_vectorized)
                    else
                        call numerical_jacobian_at_point(fun, t_new, y_predict, atol_use, &
                            state%jac_factor, state%jac, lband, uband, &
                            jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
                    end if
                    state%jacobian_initialized = .true.
                    njev = njev + 1
                    if (.not. state%jacobian_sparse) then
                        if (.not. all(ieee_is_finite(state%jac))) then
                            status = bdf_status_nonfinite
                            return
                        end if
                    end if
                    state%lu_valid = .false.
                    state%jacobian_age = 0
                    state%jacobian_c = c
                    state%jacobian_c_valid = .true.
                    current_jac = .true.
                end if

                if (.not. state%lu_valid) then
                    if (state%jacobian_sparse) then
                        call sparse_lu_factor(state%sparse_jac, 1.0_dp, -c, &
                            state%sparse_lu, info)
                    else
                        call factor_newton_matrix(state%jac, c, state%lu, state%piv, info)
                    end if
                    nlu = nlu + 1
                    if (info /= 0) then
                        status = bdf_status_singular
                        return
                    end if
                    state%lu_valid = .true.
                    if (state%odepack_jacobian_refresh .and. &
                            .not. state%jacobian_c_valid) then
                        state%jacobian_c = c
                        state%jacobian_c_valid = .true.
                    end if
                end if

                if (have_jac_f) then
                    call solve_bdf_system(fun, t_new, y_predict, c, psi, state%lu, &
                        state%piv, scale, state%newton_tol, y_new, d_corr, converged, &
                        n_iter, nfev, status, jac_f, state%odepack_jacobian_refresh, &
                        sparse_lu=state%sparse_lu)
                    have_jac_f = .false.
                else
                    call solve_bdf_system(fun, t_new, y_predict, c, psi, state%lu, &
                        state%piv, scale, state%newton_tol, y_new, d_corr, converged, &
                        n_iter, nfev, status, use_max_norm=state%odepack_jacobian_refresh, &
                        sparse_lu=state%sparse_lu)
                end if
                if (status /= bdf_status_success) return
                if (converged) exit
                if (current_jac) exit

                if (present(jac_sparse)) then
                    call jac_sparse(t_new, y_predict, state%sparse_jac)
                    if (.not. sparse_jacobian_is_valid(state%sparse_jac, size(y_predict))) then
                        status = bdf_status_invalid_jacobian
                        return
                    end if
                else if (present(jac)) then
                    call evaluate_user_jacobian(jac, t_new, y_predict, state%jac, &
                        lband, uband)
                else if (state%odepack_jacobian_refresh) then
                    call fun(t_new, y_predict, jac_f)
                    nfev = nfev + 1
                    if (.not. all(ieee_is_finite(jac_f))) then
                        status = bdf_status_nonfinite
                        return
                    end if
                    jac_nfev = 0
                    call odepack_numerical_jacobian(fun, t_new, y_predict, &
                        state%d(:, 1), jac_f, h, rtol_use, atol_use, state%jac, lband, &
                        uband, jac_nfev)
                    nfev = nfev + jac_nfev
                    have_jac_f = .true.
                else if (state%jacobian_sparse) then
                    call numerical_sparse_jacobian_at_point(fun, t_new, y_predict, &
                        atol_use, state%jac_factor, state%sparse_jac, jac_sparsity, &
                        fun_vectorized)
                else
                    call numerical_jacobian_at_point(fun, t_new, y_predict, atol_use, &
                        state%jac_factor, state%jac, lband, uband, &
                        jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
                end if
                state%jacobian_initialized = .true.
                njev = njev + 1
                if (.not. state%jacobian_sparse) then
                    if (.not. all(ieee_is_finite(state%jac))) then
                        status = bdf_status_nonfinite
                        return
                    end if
                end if
                state%lu_valid = .false.
                state%jacobian_age = 0
                state%jacobian_c = c
                state%jacobian_c_valid = .true.
                current_jac = .true.
            end do

            if (.not. converged) then
                factor = 0.5_dp
                if (h_abs * factor < retry_min_step) then
                    status = bdf_status_too_small_step
                    return
                end if
                h_abs = h_abs * factor
                call change_d(state%d, order, factor)
                state%n_equal_steps = 0
                state%lu_valid = .false.
                if (state%odepack_jacobian_refresh) current_jac = .false.
                rejected = rejected + 1
                cycle
            end if

            safety_newton = 0.9_dp * real(2 * newton_maxiter + 1, dp) / &
                real(2 * newton_maxiter + n_iter, dp)
            if (.not. state%odepack_jacobian_refresh) then
                scale = max(atol_use + rtol_use * abs(y_new), tiny(1.0_dp))
            end if
            error = error_const(order) * d_corr
            if (state%odepack_jacobian_refresh) then
                error_norm = maxval(abs(error / scale))
            else
                error_norm = rms_norm(error / scale)
            end if
            if (error_norm > 1.0_dp) then
                factor = max(min_factor, &
                    safety_newton * error_norm**(-1.0_dp / real(order + 1, dp)))
                if (h_abs * factor < retry_min_step) then
                    status = bdf_status_too_small_step
                    return
                end if
                h_abs = h_abs * factor
                call change_d(state%d, order, factor)
                state%n_equal_steps = 0
                if (state%odepack_jacobian_refresh) current_jac = .false.
                rejected = rejected + 1
            else
                accepted = .true.
            end if
        end do

        state%n_equal_steps = state%n_equal_steps + 1
        state%h_abs = h_abs
        state%last_error_norm = error_norm
        state%last_step_h_abs = h_abs
        state%last_step_order = order
        if (state%odepack_jacobian_refresh) then
            state%jacobian_age = state%jacobian_age + 1
        end if

        state%d(:, order + 3) = d_corr - state%d(:, order + 2)
        state%d(:, order + 2) = d_corr
        do i = order + 1, 1, -1
            state%d(:, i) = state%d(:, i) + state%d(:, i + 1)
        end do

        if (state%n_equal_steps >= order + 1) then
            if (order > 1) then
                error_m = error_const(order - 1) * state%d(:, order + 1)
                if (state%odepack_jacobian_refresh) then
                    error_m_norm = maxval(abs(error_m / scale))
                else
                    error_m_norm = rms_norm(error_m / scale)
                end if
            else
                error_m_norm = huge(1.0_dp)
            end if
            if (order < state%max_order) then
                error_p = error_const(order + 1) * state%d(:, order + 3)
                if (state%odepack_jacobian_refresh) then
                    error_p_norm = maxval(abs(error_p / scale))
                else
                    error_p_norm = rms_norm(error_p / scale)
                end if
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
            if (.not. state%odepack_jacobian_refresh) state%lu_valid = .false.
        else
            state%order = order
        end if

        dense_order = state%order
        dense_h = state%h_abs * direction
        dense_anchor = t_new
        dense_d = 0.0_dp
        dense_d(:, 1:dense_order + 1) = state%d(:, 1:dense_order + 1)
        if (present(f_new_out)) f_new_out = (d_corr + psi) / c
    end subroutine bdf_advance

    subroutine bdf_advance_odepack(state, fun, t, t_bound, direction, rtol, atol, &
            max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, nfev, &
            njev, nlu, attempt_limit, attempts, rejected, status, jac, min_step_user, &
            lband, uband, f_new_out, rtol_vector)
        type(bdf_state), intent(inout) :: state !! mutable LSODA BDF Nordsieck state
        procedure(bdf_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: t_bound !! final integration bound
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: rtol !! scalar relative tolerance
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        real(dp), intent(out) :: y_new(:) !! accepted state
        real(dp), intent(out) :: t_new !! accepted independent-variable value
        real(dp), intent(out) :: dense_d(:, :) !! accepted-step backward differences
        real(dp), intent(out) :: dense_h !! signed accepted-step spacing
        integer, intent(out) :: dense_order !! accepted BDF order
        real(dp), intent(out) :: dense_anchor !! accepted endpoint anchoring dense output
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        integer, intent(inout) :: nlu !! accumulated LU factorization count
        integer, intent(in) :: attempt_limit !! maximum trial steps available to this call
        integer, intent(out) :: attempts !! trial steps consumed by this accepted step
        integer, intent(out) :: rejected !! rejected trial-step count
        integer, intent(out) :: status !! BDF status code
        procedure(bdf_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: min_step_user !! minimum retry step magnitude
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth
        real(dp), intent(out), optional :: f_new_out(:) !! accepted polynomial derivative
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        real(dp) :: acor(size(y_new))
        real(dp) :: current_f(size(y_new))
        real(dp) :: el(bdf_max_order + 1)
        real(dp) :: scale(size(y_new))
        real(dp) :: rtol_use(size(y_new))
        real(dp) :: tesco(3)
        real(dp) :: yh_start(size(y_new), bdf_max_order + 1)
        real(dp) :: yh_trial(size(y_new), bdf_max_order + 1)
        real(dp) :: y_predict(size(y_new))
        real(dp) :: c
        real(dp) :: ddn
        real(dp) :: dsm
        real(dp) :: dup
        real(dp) :: factor
        real(dp) :: h
        real(dp) :: h_abs
        real(dp) :: min_step
        real(dp) :: retry_min_step
        real(dp) :: rhdn
        real(dp) :: rhsm
        real(dp) :: rhup
        integer :: accepted_order
        integer :: error_failures
        integer :: info
        integer :: jac_nfev
        integer :: n_iter
        integer :: ncf
        integer :: new_order
        integer :: order
        logical :: converged
        logical :: current_jac
        logical :: refresh_jacobian

        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector

        status = bdf_status_success
        attempts = 0
        rejected = 0
        dense_d = 0.0_dp
        dense_h = 0.0_dp
        dense_order = 0
        dense_anchor = t
        if (.not. state%odepack_yh_valid) then
            status = bdf_status_nonfinite
            return
        end if

        order = state%order
        h_abs = state%h_abs
        min_step = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(t))
        retry_min_step = min_step
        if (present(min_step_user)) retry_min_step = max(retry_min_step, min_step_user)
        yh_start = state%odepack_yh

        if (h_abs > max_step) then
            factor = max_step / h_abs
            call rescale_nordsieck(yh_start, order, factor)
            h_abs = max_step
            state%odepack_ialth = order + 1
            state%odepack_saved_acor_valid = .false.
        else if (h_abs < min_step) then
            factor = min_step / h_abs
            call rescale_nordsieck(yh_start, order, factor)
            h_abs = min_step
            state%odepack_ialth = order + 1
            state%odepack_saved_acor_valid = .false.
        end if

        error_failures = 0
        ncf = 0
        do
            if (attempts >= attempt_limit) then
                status = bdf_status_max_attempts
                return
            end if
            attempts = attempts + 1
            if (h_abs < min_step) then
                status = bdf_status_too_small_step
                return
            end if

            h = direction * h_abs
            t_new = t + h
            if (direction * (t_new - t_bound) > 0.0_dp) then
                t_new = t_bound
                factor = abs(t_new - t) / h_abs
                call rescale_nordsieck(yh_start, order, factor)
                h_abs = abs(t_new - t)
                h = t_new - t
                state%odepack_ialth = order + 1
                state%odepack_saved_acor_valid = .false.
            end if

            scale = max(atol + rtol_use * abs(yh_start(:, 1)), tiny(1.0_dp))
            call odepack_bdf_coefficients(order, el, tesco)
            yh_trial = yh_start
            call predict_nordsieck(yh_trial, order)
            y_predict = yh_trial(:, 1)
            c = h * el(1)

            call fun(t_new, y_predict, current_f)
            nfev = nfev + 1
            if (.not. all(ieee_is_finite(current_f))) then
                status = bdf_status_nonfinite
                return
            end if

            current_jac = .false.
            do
                refresh_jacobian = .not. state%jacobian_initialized
                if (.not. current_jac .and. state%jacobian_initialized) then
                    refresh_jacobian = state%jacobian_age >= odepack_jacobian_refresh_steps
                    if (state%jacobian_c_valid) then
                        refresh_jacobian = refresh_jacobian .or. &
                            abs(c / state%jacobian_c - 1.0_dp) > &
                            odepack_jacobian_refresh_ratio
                    end if
                end if
                if (refresh_jacobian) then
                    if (present(jac)) then
                        call evaluate_user_jacobian(jac, t_new, y_predict, state%jac, &
                            lband, uband)
                    else
                        jac_nfev = 0
                        call odepack_numerical_jacobian(fun, t_new, y_predict, &
                            yh_start(:, 1), current_f, h, rtol_use, atol, state%jac, lband, &
                            uband, jac_nfev)
                        nfev = nfev + jac_nfev
                    end if
                    state%jacobian_initialized = .true.
                    njev = njev + 1
                    if (.not. all(ieee_is_finite(state%jac))) then
                        status = bdf_status_nonfinite
                        return
                    end if
                    call factor_newton_matrix(state%jac, c, state%lu, state%piv, info)
                    nlu = nlu + 1
                    if (info /= 0) then
                        status = bdf_status_singular
                        return
                    end if
                    state%lu_valid = .true.
                    state%jacobian_age = 0
                    state%jacobian_c = c
                    state%jacobian_c_valid = .true.
                    state%odepack_crate = 0.7_dp
                    current_jac = .true.
                else if (.not. state%lu_valid) then
                    call factor_newton_matrix(state%jac, c, state%lu, state%piv, info)
                    nlu = nlu + 1
                    if (info /= 0) then
                        status = bdf_status_singular
                        return
                    end if
                    state%lu_valid = .true.
                end if

                call solve_odepack_chord(fun, t_new, yh_trial, h, el(1), tesco(2), &
                    scale, order, state%lu, state%piv, state%odepack_crate, current_f, y_new, &
                    acor, converged, n_iter, nfev, status)
                if (status /= bdf_status_success) return
                if (converged) exit
                if (current_jac) exit

                call fun(t_new, y_predict, current_f)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(current_f))) then
                    status = bdf_status_nonfinite
                    return
                end if
                if (present(jac)) then
                    call evaluate_user_jacobian(jac, t_new, y_predict, state%jac, &
                        lband, uband)
                else
                    jac_nfev = 0
                    call odepack_numerical_jacobian(fun, t_new, y_predict, &
                        yh_start(:, 1), current_f, h, rtol_use, atol, state%jac, lband, &
                        uband, jac_nfev)
                    nfev = nfev + jac_nfev
                end if
                state%jacobian_initialized = .true.
                njev = njev + 1
                call factor_newton_matrix(state%jac, c, state%lu, state%piv, info)
                nlu = nlu + 1
                if (info /= 0) then
                    status = bdf_status_singular
                    return
                end if
                state%lu_valid = .true.
                state%jacobian_age = 0
                state%jacobian_c = c
                state%jacobian_c_valid = .true.
                state%odepack_crate = 0.7_dp
                current_jac = .true.
            end do

            if (.not. converged) then
                ncf = ncf + 1
                rejected = rejected + 1
                state%odepack_rmax = 2.0_dp
                if (ncf >= 10) then
                    status = bdf_status_max_attempts
                    return
                end if
                factor = 0.25_dp
                if (h_abs * factor < retry_min_step) then
                    status = bdf_status_too_small_step
                    return
                end if
                h_abs = h_abs * factor
                call rescale_nordsieck(yh_start, order, factor)
                state%odepack_ialth = order + 1
                state%odepack_saved_acor_valid = .false.
                state%lu_valid = .false.
                cycle
            end if

            dsm = maxval(abs(acor / scale)) / tesco(2)
            if (dsm <= 1.0_dp) exit

            error_failures = error_failures + 1
            rejected = rejected + 1
            state%odepack_rmax = 2.0_dp
            if (h_abs <= retry_min_step * (1.0_dp + 1.0e-12_dp)) then
                status = bdf_status_too_small_step
                return
            end if
            if (error_failures >= 3) then
                if (error_failures >= 10) then
                    status = bdf_status_max_attempts
                    return
                end if
                factor = max(0.1_dp, retry_min_step / h_abs)
                h_abs = h_abs * factor
                call rescale_nordsieck(yh_start, order, factor)
                call fun(t, yh_start(:, 1), current_f)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(current_f))) then
                    status = bdf_status_nonfinite
                    return
                end if
                order = 1
                yh_start(:, 2) = direction * h_abs * current_f
                yh_start(:, 3:) = 0.0_dp
                state%odepack_ialth = 5
                state%odepack_saved_acor_valid = .false.
                state%lu_valid = .false.
                cycle
            end if

            rhsm = odepack_step_factor(dsm, order + 1, 1.2_dp)
            new_order = order
            factor = rhsm
            if (order > 1) then
                ddn = maxval(abs(yh_start(:, order + 1) / scale)) / tesco(1)
                rhdn = odepack_step_factor(ddn, order, 1.3_dp)
                if (rhdn > rhsm) then
                    new_order = order - 1
                    factor = min(1.0_dp, rhdn)
                end if
            end if
            if (error_failures >= 2) factor = min(factor, 0.2_dp)
            factor = max(factor, retry_min_step / h_abs)
            order = new_order
            h_abs = h_abs * factor
            call rescale_nordsieck(yh_start, order, factor)
            state%odepack_ialth = order + 1
            state%odepack_saved_acor_valid = .false.
        end do

        accepted_order = order
        do info = 1, order + 1
            yh_trial(:, info) = yh_trial(:, info) + el(info) * acor
        end do
        y_new = yh_trial(:, 1)
        dense_order = order
        dense_h = h
        dense_anchor = t_new
        call nordsieck_to_backward_differences(yh_trial, order, dense_d)

        state%last_error_norm = dsm
        state%last_step_h_abs = h_abs
        state%last_step_order = accepted_order
        state%jacobian_age = state%jacobian_age + 1
        state%odepack_yh = yh_trial
        state%h_abs = h_abs
        state%order = order

        state%odepack_ialth = state%odepack_ialth - 1
        if (state%odepack_ialth <= 0) then
            rhup = 0.0_dp
            if (order < state%max_order .and. state%odepack_saved_acor_valid) then
                dup = maxval(abs((acor - state%odepack_saved_acor) / scale)) / tesco(3)
                rhup = odepack_step_factor(dup, order + 2, 1.4_dp)
            end if
            rhsm = odepack_step_factor(dsm, order + 1, 1.2_dp)
            rhdn = 0.0_dp
            if (order > 1) then
                ddn = maxval(abs(state%odepack_yh(:, order + 1) / scale)) / tesco(1)
                rhdn = odepack_step_factor(ddn, order, 1.3_dp)
            end if

            new_order = order
            factor = rhsm
            if (rhup > factor .and. rhup > rhdn) then
                new_order = order + 1
                factor = rhup
            else if (rhdn > factor) then
                new_order = order - 1
                factor = rhdn
            end if

            if (factor < 1.1_dp) then
                state%odepack_ialth = 3
            else
                if (new_order == order + 1) then
                    state%odepack_yh(:, new_order + 1) = &
                        acor * el(order + 1) / real(order + 1, dp)
                end if
                order = new_order
                factor = min(factor, state%odepack_rmax)
                if (h_abs > max_step / factor) factor = max_step / h_abs
                call rescale_nordsieck(state%odepack_yh, order, factor)
                state%h_abs = state%h_abs * factor
                state%order = order
                state%odepack_ialth = order + 1
                state%odepack_rmax = max_factor
                state%odepack_saved_acor_valid = .false.
            end if
        else if (state%odepack_ialth <= 1 .and. order < state%max_order) then
            state%odepack_saved_acor = acor
            state%odepack_saved_acor_valid = .true.
        end if

        state%d = 0.0_dp
        call nordsieck_to_backward_differences(state%odepack_yh, state%order, state%d)
        state%n_equal_steps = max(0, state%order + 1 - state%odepack_ialth)
        if (present(f_new_out)) then
            f_new_out = state%odepack_yh(:, 2) / (direction * state%h_abs)
        end if
    end subroutine bdf_advance_odepack

    subroutine solve_odepack_chord(fun, t_new, yh, h, el0, tesco2, scale, order, lu, piv, &
            crate, f_initial, y, acor, converged, n_iter, nfev, status)
        procedure(bdf_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t_new !! trial endpoint
        real(dp), intent(in) :: yh(:, :) !! predicted BDF Nordsieck array
        real(dp), intent(in) :: h !! signed trial-step size
        real(dp), intent(in) :: el0 !! leading BDF correction coefficient
        real(dp), intent(in) :: tesco2 !! local-error/convergence test constant
        real(dp), intent(in) :: scale(:) !! accepted-state error scales
        integer, intent(in) :: order !! active BDF method order
        real(dp), intent(in) :: lu(:, :) !! factored chord matrix
        integer, intent(in) :: piv(:) !! LU row pivots
        real(dp), intent(inout) :: crate !! chord convergence-rate estimate
        real(dp), intent(in) :: f_initial(:) !! RHS at the predictor
        real(dp), intent(out) :: y(:) !! final chord iterate
        real(dp), intent(out) :: acor(:) !! accumulated raw ODEPACK correction
        logical, intent(out) :: converged !! true when the chord test passes
        integer, intent(out) :: n_iter !! number of corrections performed
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(out) :: status !! BDF status code

        real(dp) :: dy(size(y))
        real(dp) :: f(size(y))
        real(dp) :: rhs(size(y))
        real(dp) :: dcon
        real(dp) :: del
        real(dp) :: delp
        real(dp) :: pnorm
        real(dp) :: rm
        real(dp) :: conit
        integer :: info
        integer :: m

        y = yh(:, 1)
        acor = 0.0_dp
        f = f_initial
        converged = .false.
        status = bdf_status_success
        n_iter = 0
        delp = 0.0_dp
        pnorm = maxval(abs(yh(:, 1) / scale))
        conit = 0.5_dp / real(order + 2, dp)

        do m = 0, 2
            n_iter = m + 1
            if (m > 0) then
                call fun(t_new, y, f)
                nfev = nfev + 1
                if (.not. all(ieee_is_finite(f))) then
                    status = bdf_status_nonfinite
                    return
                end if
            end if
            rhs = h * f - (yh(:, 2) + acor)
            call lu_solve(lu, piv, rhs, dy, info)
            if (info /= 0) then
                status = bdf_status_singular
                return
            end if
            del = maxval(abs(dy / scale))
            acor = acor + dy
            y = yh(:, 1) + el0 * acor

            if (del <= 100.0_dp * pnorm * epsilon(1.0_dp)) then
                converged = .true.
                exit
            end if
            if (m > 0) then
                rm = 1024.0_dp
                if (del <= 1024.0_dp * delp .and. delp > 0.0_dp) rm = del / delp
                crate = max(0.2_dp * crate, rm)
            end if
            dcon = del * min(1.0_dp, 1.5_dp * crate) / (tesco2 * conit)
            if (dcon <= 1.0_dp) then
                converged = .true.
                exit
            end if
            if (m == 2) exit
            if (m + 1 >= 2 .and. del > 2.0_dp * delp) exit
            delp = del
        end do
    end subroutine solve_odepack_chord


    subroutine odepack_bdf_coefficients(order, el, tesco)
        integer, intent(in) :: order !! BDF order, 1 through 5
        real(dp), intent(out) :: el(:) !! correction polynomial coefficients
        real(dp), intent(out) :: tesco(3) !! lower/current/upper error-test constants

        real(dp) :: pc(bdf_max_order + 1)
        real(dp) :: fnq
        real(dp) :: rq1fac
        integer :: i
        integer :: ib
        integer :: nq

        el = 0.0_dp
        tesco = 0.0_dp
        pc = 0.0_dp
        pc(1) = 1.0_dp
        rq1fac = 1.0_dp
        do nq = 1, order
            fnq = real(nq, dp)
            pc(nq + 1) = 0.0_dp
            do ib = 1, nq
                i = nq + 2 - ib
                pc(i) = pc(i - 1) + fnq * pc(i)
            end do
            pc(1) = fnq * pc(1)
            if (nq == order) then
                el(1:nq + 1) = pc(1:nq + 1) / pc(2)
                el(2) = 1.0_dp
                tesco(1) = rq1fac
                tesco(2) = real(nq + 1, dp) / el(1)
                tesco(3) = real(nq + 2, dp) / el(1)
            end if
            rq1fac = rq1fac / fnq
        end do
    end subroutine odepack_bdf_coefficients

    subroutine predict_nordsieck(yh, order)
        real(dp), intent(inout) :: yh(:, :) !! Nordsieck array to advance one unit step
        integer, intent(in) :: order !! active polynomial order

        integer :: jb
        integer :: j

        do jb = 1, order
            do j = order - jb + 1, order
                yh(:, j) = yh(:, j) + yh(:, j + 1)
            end do
        end do
    end subroutine predict_nordsieck

    subroutine rescale_nordsieck(yh, order, factor)
        real(dp), intent(inout) :: yh(:, :) !! Nordsieck array to rescale in place
        integer, intent(in) :: order !! active polynomial order
        real(dp), intent(in) :: factor !! new step divided by old step

        real(dp) :: power
        integer :: j

        power = 1.0_dp
        do j = 1, order
            power = power * factor
            yh(:, j + 1) = yh(:, j + 1) * power
        end do
    end subroutine rescale_nordsieck

    subroutine nordsieck_to_backward_differences(yh, order, d)
        real(dp), intent(in) :: yh(:, :) !! endpoint Nordsieck polynomial
        integer, intent(in) :: order !! polynomial order
        real(dp), intent(out) :: d(:, :) !! backward-difference representation

        real(dp) :: differences(size(yh, 1), 0:bdf_max_order)
        real(dp) :: samples(size(yh, 1), 0:bdf_max_order)
        real(dp) :: x
        integer :: j
        integer :: k

        d = 0.0_dp
        samples = 0.0_dp
        do j = 0, order
            x = -real(j, dp)
            samples(:, j) = yh(:, order + 1)
            do k = order - 1, 0, -1
                samples(:, j) = samples(:, j) * x + yh(:, k + 1)
            end do
        end do
        differences = samples
        d(:, 1) = differences(:, 0)
        do k = 1, order
            do j = 0, order - k
                differences(:, j) = differences(:, j) - differences(:, j + 1)
            end do
            d(:, k + 1) = differences(:, 0)
        end do
    end subroutine nordsieck_to_backward_differences

    pure real(dp) function odepack_step_factor(error_norm, exponent_denominator, &
            coefficient) result(value)
        real(dp), intent(in) :: error_norm !! normalized local-error estimate
        integer, intent(in) :: exponent_denominator !! positive adaptation exponent denominator
        real(dp), intent(in) :: coefficient !! DSTODA denominator coefficient

        if (error_norm == 0.0_dp) then
            value = huge(1.0_dp)
        else if (.not. ieee_is_finite(error_norm)) then
            value = 0.0_dp
        else
            value = 1.0_dp / (coefficient * error_norm**(1.0_dp / &
                real(exponent_denominator, dp)) + coefficient * 1.0e-6_dp)
        end if
    end function odepack_step_factor

    subroutine solve_bdf_system(fun, t_new, y_predict, c, psi, lu, piv, scale, tol, &
            y, d, converged, n_iter, nfev, status, f_initial, use_max_norm, sparse_lu)
        procedure(bdf_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t_new !! trial-step endpoint
        real(dp), intent(in) :: y_predict(:) !! extrapolated predictor at t_new
        real(dp), intent(in) :: c !! BDF nonlinear-system scale h/alpha(order)
        real(dp), intent(in) :: psi(:) !! BDF accumulated backward-difference term
        real(dp), intent(in) :: lu(:, :) !! LU factors of I-cJ
        integer, intent(in) :: piv(:) !! row pivots paired with lu
        real(dp), intent(in) :: scale(:) !! component error scales
        real(dp), intent(in) :: tol !! Newton convergence tolerance
        real(dp), intent(out) :: y(:) !! final Newton iterate
        real(dp), intent(out) :: d(:) !! accumulated Newton correction
        logical, intent(out) :: converged !! true when Newton stopping criterion is satisfied
        integer, intent(out) :: n_iter !! Newton iterations performed
        integer, intent(inout) :: nfev !! accumulated RHS evaluation count
        integer, intent(out) :: status !! BDF status code
        real(dp), intent(in), optional :: f_initial(:) !! reusable RHS at the predictor
        logical, intent(in), optional :: use_max_norm !! use ODEPACK weighted maximum norm
        type(ivp_sparse_lu_real), intent(in), optional :: &
            sparse_lu !! sparse LU factors used when a constant CSC Jacobian is active

        real(dp) :: f(size(y))
        real(dp) :: rhs(size(y))
        real(dp) :: dy(size(y))
        real(dp) :: dy_norm
        real(dp) :: dy_norm_old
        real(dp) :: rate
        integer :: k
        integer :: info
        logical :: have_old_norm
        logical :: use_initial_f
        logical :: use_max

        y = y_predict
        d = 0.0_dp
        converged = .false.
        have_old_norm = .false.
        dy_norm_old = 0.0_dp
        status = bdf_status_success
        n_iter = 0
        use_initial_f = present(f_initial)
        use_max = .false.
        if (present(use_max_norm)) use_max = use_max_norm

        do k = 1, newton_maxiter
            n_iter = k
            if (k == 1 .and. use_initial_f) then
                f = f_initial
            else
                call fun(t_new, y, f)
                nfev = nfev + 1
            end if
            if (.not. all(ieee_is_finite(f))) then
                status = bdf_status_nonfinite
                return
            end if
            rhs = c * f - psi - d
            if (present(sparse_lu) .and. sparse_lu%n > 0) then
                call sparse_lu%solve(rhs, dy, info)
            else
                call lu_solve(lu, piv, rhs, dy, info)
            end if
            if (info /= 0) then
                status = bdf_status_singular
                return
            end if
            if (use_max) then
                dy_norm = maxval(abs(dy / scale))
            else
                dy_norm = rms_norm(dy / scale)
            end if

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
    end subroutine solve_bdf_system

    subroutine bdf_dense_state(t_anchor, h, order, d, tq, yq)
        real(dp), intent(in) :: t_anchor !! full accepted BDF endpoint anchoring the polynomial
        real(dp), intent(in) :: h !! signed BDF difference spacing used by the dense polynomial
        integer, intent(in) :: order !! BDF dense polynomial order, 1 through 5
        real(dp), intent(in) :: d(:, :) !! backward-difference coefficients, columns 1:order+1
        real(dp), intent(in) :: tq !! interpolation point in the accepted interval
        real(dp), intent(out) :: yq(:) !! interpolated state

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
    end subroutine bdf_dense_state

    subroutine evaluate_user_jacobian(jac_fun, t, y, jac, lband, uband)
        procedure(bdf_jacobian) :: jac_fun !! user Jacobian callback
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(out) :: jac(:, :) !! unpacked dense Jacobian
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense

        real(dp), allocatable :: packed(:, :)
        integer :: lower
        integer :: upper

        if (.not. present(lband) .and. .not. present(uband)) then
            call jac_fun(t, y, jac)
            return
        end if

        lower = 0
        upper = 0
        if (present(lband)) lower = lband
        if (present(uband)) upper = uband
        allocate(packed(lower + upper + 1, size(y)))
        call jac_fun(t, y, packed)
        call unpack_banded_jacobian(packed, lower, upper, jac)
    end subroutine evaluate_user_jacobian

    subroutine unpack_banded_jacobian(packed, lband, uband, jac)
        real(dp), intent(in) :: packed(:, :) !! packed band matrix, diagonals by row
        integer, intent(in) :: lband !! lower Jacobian bandwidth
        integer, intent(in) :: uband !! upper Jacobian bandwidth
        real(dp), intent(out) :: jac(:, :) !! unpacked dense Jacobian

        integer :: i
        integer :: j

        jac = 0.0_dp
        do j = 1, size(jac, 2)
            do i = max(1, j - uband), min(size(jac, 1), j + lband)
                jac(i, j) = packed(uband + i - j + 1, j)
            end do
        end do
    end subroutine unpack_banded_jacobian


    subroutine odepack_numerical_jacobian(fun, t, y, y_weight, f, h_step, rtol, atol, &
            jac, lband, uband, nfev_jac)
        procedure(bdf_rhs) :: fun !! right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! predicted state where the Jacobian is evaluated
        real(dp), intent(in) :: y_weight(:) !! accepted state used to form LSODA error weights
        real(dp), intent(in) :: f(:) !! fun(t,y), already counted by the caller
        real(dp), intent(in) :: h_step !! signed trial-step size used by the stiff method
        real(dp), intent(in) :: rtol(:) !! componentwise relative tolerances
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        real(dp), intent(out) :: jac(:, :) !! finite-difference Jacobian in dense storage
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense
        integer, intent(inout) :: nfev_jac !! RHS calls made for Jacobian perturbations

        real(dp) :: f_trial(size(y))
        real(dp) :: y_trial(size(y))
        real(dp) :: scale(size(y))
        real(dp) :: perturb(size(y))
        real(dp) :: f_norm
        real(dp) :: r0
        real(dp) :: srur
        integer :: lower
        integer :: upper
        integer :: width
        integer :: group
        integer :: i
        integer :: j
        integer :: ilo
        integer :: ihi

        scale = max(atol + rtol * abs(y_weight), tiny(1.0_dp))
        f_norm = maxval(abs(f) / scale)
        r0 = 1000.0_dp * abs(h_step) * epsilon(1.0_dp) * real(size(y), dp) * f_norm
        if (r0 == 0.0_dp) r0 = 1.0_dp
        srur = sqrt(epsilon(1.0_dp))
        perturb = max(srur * abs(y), r0 * scale)
        do j = 1, size(y)
            if (perturb(j) == 0.0_dp) perturb(j) = srur
        end do

        lower = 0
        upper = 0
        if (present(lband)) lower = lband
        if (present(uband)) upper = uband
        if (.not. present(lband) .and. .not. present(uband)) then
            do j = 1, size(y)
                y_trial = y
                y_trial(j) = y_trial(j) + perturb(j)
                call fun(t, y_trial, f_trial)
                nfev_jac = nfev_jac + 1
                jac(:, j) = (f_trial - f) / perturb(j)
            end do
            return
        end if

        width = lower + upper + 1
        jac = 0.0_dp
        do group = 1, min(width, size(y))
            y_trial = y
            do j = group, size(y), width
                y_trial(j) = y_trial(j) + perturb(j)
            end do
            call fun(t, y_trial, f_trial)
            nfev_jac = nfev_jac + 1
            do j = group, size(y), width
                ilo = max(1, j - upper)
                ihi = min(size(y), j + lower)
                do i = ilo, ihi
                    jac(i, j) = (f_trial(i) - f(i)) / perturb(j)
                end do
            end do
        end do
    end subroutine odepack_numerical_jacobian

    subroutine numerical_banded_jacobian(fun, t, y, f, atol, factor, jac, lband, uband, &
            nfev_jac)
        procedure(bdf_rhs) :: fun !! right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        real(dp), intent(out) :: jac(:, :) !! banded finite-difference Jacobian in dense storage
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means zero
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means zero
        integer, intent(inout), optional :: nfev_jac !! RHS calls made for Jacobian differences

        real(dp) :: f_group(size(y))
        real(dp) :: f_retry(size(y))
        real(dp) :: y_trial(size(y))
        real(dp) :: h(size(y))
        real(dp) :: y_scale(size(y))
        real(dp) :: max_diff
        real(dp) :: scale_ref
        real(dp) :: max_diff_new
        real(dp) :: scale_new
        real(dp) :: new_factor
        real(dp) :: h_retry
        integer :: lower
        integer :: upper
        integer :: width
        integer :: group
        integer :: i
        integer :: j
        integer :: ilo
        integer :: ihi
        integer :: imax
        logical :: use_retry

        lower = 0
        upper = 0
        if (present(lband)) lower = lband
        if (present(uband)) upper = uband
        width = lower + upper + 1
        jac = 0.0_dp

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

        do group = 1, width
            y_trial = y
            do j = group, size(y), width
                y_trial(j) = y(j) + h(j)
            end do
            call fun(t, y_trial, f_group)
            if (present(nfev_jac)) nfev_jac = nfev_jac + 1

            do j = group, size(y), width
                ilo = max(1, j - upper)
                ihi = min(size(y), j + lower)
                jac(ilo:ihi, j) = (f_group(ilo:ihi) - f(ilo:ihi)) / h(j)
                imax = maxloc(abs(f_group(ilo:ihi) - f(ilo:ihi)), dim=1) + ilo - 1
                max_diff = abs(f_group(imax) - f(imax))
                scale_ref = max(abs(f(imax)), abs(f_group(imax)))
                use_retry = max_diff < num_jac_diff_reject * scale_ref
                if (use_retry) then
                    new_factor = num_jac_factor_increase * factor(j)
                    h_retry = (y(j) + new_factor * y_scale(j)) - y(j)
                    y_trial = y
                    y_trial(j) = y(j) + h_retry
                    call fun(t, y_trial, f_retry)
                    if (present(nfev_jac)) nfev_jac = nfev_jac + 1
                    i = maxloc(abs(f_retry(ilo:ihi) - f(ilo:ihi)), dim=1) + ilo - 1
                    max_diff_new = abs(f_retry(i) - f(i))
                    scale_new = max(abs(f(i)), abs(f_retry(i)))
                    if (max_diff * scale_new < max_diff_new * scale_ref) then
                        factor(j) = new_factor
                        h(j) = h_retry
                        jac(ilo:ihi, j) = (f_retry(ilo:ihi) - f(ilo:ihi)) / h(j)
                        max_diff = max_diff_new
                        scale_ref = scale_new
                    end if
                end if
                if (max_diff < num_jac_diff_small * scale_ref) then
                    factor(j) = factor(j) * num_jac_factor_increase
                else if (max_diff > num_jac_diff_big * scale_ref) then
                    factor(j) = factor(j) * num_jac_factor_decrease
                end if
                factor(j) = max(factor(j), num_jac_min_factor)
            end do
        end do
    end subroutine numerical_banded_jacobian

    subroutine numerical_jacobian(fun, t, y, f, atol, factor, jac, lband, uband, &
            nfev_jac, jac_sparsity, fun_vectorized)
        procedure(bdf_rhs) :: fun !! scalar right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        real(dp), intent(out) :: jac(:, :) !! finite-difference Jacobian
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense
        integer, intent(inout), optional :: nfev_jac !! RHS calls made for Jacobian differences
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_rhs_vectorized), optional :: &
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

        if (present(lband) .or. present(uband)) then
            call numerical_banded_jacobian(fun, t, y, f, atol, factor, jac, lband, &
                uband, nfev_jac)
            return
        end if

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
                if (present(nfev_jac)) nfev_jac = nfev_jac + 1
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
                    if (present(nfev_jac)) nfev_jac = nfev_jac + 1
                    do j = 1, size(y)
                        if (groups(j) == g) f_new(:, j) = f_retry
                    end do
                end do
            else
                do j = 1, size(y)
                    y_trial = y
                    y_trial(j) = y(j) + h(j)
                    call fun(t, y_trial, f_new(:, j))
                    if (present(nfev_jac)) nfev_jac = nfev_jac + 1
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
                if (present(nfev_jac)) nfev_jac = nfev_jac + 1
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

    subroutine numerical_jacobian_at_point(fun, t, y, atol, factor, jac, lband, uband, &
            nfev_jac, jac_sparsity, fun_vectorized)
        procedure(bdf_rhs) :: fun !! scalar right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        real(dp), intent(out) :: jac(:, :) !! finite-difference Jacobian
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense
        integer, intent(inout), optional :: nfev_jac !! RHS calls including the base evaluation
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS accepting state columns for Jacobian differences

        real(dp) :: f(size(y))
        real(dp) :: f_batch(size(y), 1)
        real(dp) :: y_batch(size(y), 1)

        if (present(fun_vectorized)) then
            y_batch(:, 1) = y
            call fun_vectorized(t, y_batch, f_batch)
            f = f_batch(:, 1)
        else
            call fun(t, y, f)
        end if
        if (present(nfev_jac)) nfev_jac = nfev_jac + 1
        call numerical_jacobian(fun, t, y, f, atol, factor, jac, lband, uband, nfev_jac, &
            jac_sparsity, fun_vectorized)
    end subroutine numerical_jacobian_at_point


    subroutine numerical_sparse_jacobian(fun, t, y, f, atol, factor, jac, sparsity, &
            fun_vectorized)
        procedure(bdf_rhs) :: fun !! scalar right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: f(:) !! fun(t,y), supplied without incrementing solver nfev
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        type(ivp_sparse_matrix), intent(inout) :: jac !! CSC Jacobian with fixed sparsity structure
        logical, intent(in) :: sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_rhs_vectorized), optional :: &
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

    subroutine numerical_sparse_jacobian_at_point(fun, t, y, atol, factor, jac, sparsity, &
            fun_vectorized)
        procedure(bdf_rhs) :: fun !! scalar right-hand side callback used for finite differences
        real(dp), intent(in) :: t !! Jacobian evaluation time
        real(dp), intent(in) :: y(:) !! Jacobian evaluation state
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances for perturbations
        real(dp), intent(inout) :: factor(:) !! persistent finite-difference step factors
        type(ivp_sparse_matrix), intent(inout) :: jac !! CSC Jacobian with fixed sparsity structure
        logical, intent(in) :: sparsity(:, :) !! finite-difference nonzero pattern
        procedure(bdf_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS accepting state columns for Jacobian differences

        real(dp) :: f(size(y))
        real(dp) :: f_batch(size(y), 1)
        real(dp) :: y_batch(size(y), 1)

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
        type(ivp_sparse_matrix), intent(in) :: jac !! sparse Jacobian returned by a callback
        integer, intent(in) :: n !! required square matrix order

        valid = jac%is_valid()
        if (.not. valid) return
        valid = jac%nrows == n .and. jac%ncols == n
    end function sparse_jacobian_is_valid


    subroutine factor_newton_matrix(jac, c, lu, piv, info)
        real(dp), intent(in) :: jac(:, :) !! Jacobian matrix
        real(dp), intent(in) :: c !! BDF nonlinear-system scale h/alpha(order)
        real(dp), intent(out) :: lu(:, :) !! LU factors of I-cJ
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
    end subroutine lu_factor

    subroutine lu_solve(lu, piv, b, x, info)
        real(dp), intent(in) :: lu(:, :) !! compact LU factors
        integer, intent(in) :: piv(:) !! row pivots from lu_factor
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
            x(i) = x(i) - dot_product(lu(i, 1:i - 1), x(1:i - 1))
        end do
        do i = n, 1, -1
            if (lu(i, i) == 0.0_dp) then
                info = i
                return
            end if
            if (i < n) x(i) = x(i) - dot_product(lu(i, i + 1:n), x(i + 1:n))
            x(i) = x(i) / lu(i, i)
        end do
    end subroutine lu_solve

    subroutine change_d(d, order, factor)
        real(dp), intent(inout) :: d(:, :) !! backward-difference table to transform in place
        integer, intent(in) :: order !! current BDF order
        real(dp), intent(in) :: factor !! new step size divided by old step size

        real(dp) :: r(0:bdf_max_order, 0:bdf_max_order)
        real(dp) :: u(0:bdf_max_order, 0:bdf_max_order)
        real(dp) :: ru(0:bdf_max_order, 0:bdf_max_order)
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
        d(:, 1:order + 1) = matmul(d(:, 1:order + 1), ru(0:order, 0:order))
    end subroutine change_d

    pure real(dp) function rms_norm(x) result(value)
        real(dp), intent(in) :: x(:) !! scaled vector whose RMS norm is requested

        value = sqrt(sum(x * x) / real(size(x), dp))
    end function rms_norm

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

end module scifort_integrate_bdf
