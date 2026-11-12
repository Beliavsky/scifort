! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real-valued automatic Adams/BDF initial-value stepper for the SciFort
! solve_ivp LSODA compatibility path. The persistent Adams Nordsieck update,
! DCFODE correction/error constants, step rescaling, Pascal prediction,
! IALTH/TESCO step/order selection, Adams stability limits, switch hysteresis,
! roundoff/failure recovery, method-state transfer, Jacobian refresh rules,
! and ideal-step comparison constants are adapted from public-domain ODEPACK
! DSTODA/DCFODE. SciFort retains independent derivative-history reconstruction
! and reuses its BDF backward-difference machinery with LSODA-only pure-BDF rules.
! See CODE_PROVENANCE.md.

module scifort_integrate_lsoda
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_integrate_bdf, only : bdf_advance, bdf_initialize, bdf_max_order, &
        bdf_state, bdf_status_max_attempts, bdf_status_nonfinite, &
        bdf_status_singular, bdf_status_success, bdf_status_too_small_step
    implicit none
    private

    integer, parameter, public :: lsoda_status_success = 0
    integer, parameter, public :: lsoda_status_too_small_step = 1
    integer, parameter, public :: lsoda_status_singular = 2
    integer, parameter, public :: lsoda_status_nonfinite = 3
    integer, parameter, public :: lsoda_status_max_attempts = 4

    integer, parameter :: lsoda_status_switch_bdf = 5

    integer, parameter :: mode_adams = 1
    integer, parameter :: mode_bdf = 2
    integer, parameter :: adams_history_min = 4
    integer, parameter :: adams_method_min = 1
    integer, parameter :: adams_history_max = 11
    integer, parameter :: adams_method_max = adams_history_max + 1
    real(dp), parameter :: safety = 0.9_dp
    real(dp), parameter :: min_factor = 0.2_dp
    real(dp), parameter :: adams_initial_rmax = 10000.0_dp
    real(dp), parameter :: adams_regular_rmax = 10.0_dp
    real(dp), parameter :: adams_failure_rmax = 2.0_dp
    real(dp), parameter :: stiff_switch_threshold = 0.5_dp
    real(dp), parameter :: nonstiff_switch_threshold = 0.15_dp
    integer, parameter :: switch_wait_steps = 20
    integer, parameter :: adams_failure_restart_count = 4
    integer, parameter :: adams_switch_method_order_max = 5
    real(dp), parameter :: method_switch_ratio = 5.0_dp
    real(dp), parameter :: odepack_adams_error_constant(5) = [ &
        2.0_dp, 6.0_dp, 4.0_dp, 30.0_dp / 19.0_dp, 4.0_dp / 9.0_dp]
    real(dp), parameter :: odepack_bdf_error_constant(5) = [ &
        2.0_dp, 3.0_dp / 2.0_dp, 2.0_dp / 3.0_dp, 5.0_dp / 24.0_dp, 1.0_dp / 20.0_dp]
    real(dp), parameter :: adams_stability_limit(12) = [ &
        0.5_dp, 0.575_dp, 0.55_dp, 0.45_dp, 0.35_dp, 0.25_dp, &
        0.20_dp, 0.15_dp, 0.10_dp, 0.075_dp, 0.050_dp, 0.025_dp]

    abstract interface
        subroutine lsoda_rhs(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t
        end subroutine lsoda_rhs

        subroutine lsoda_jacobian(t, y, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: jac(:, :) !! Jacobian d fun / d y, shape (n,n)
        end subroutine lsoda_jacobian
    end interface

    type, public :: lsoda_state
        real(dp), allocatable :: f_hist(:, :)
        real(dp), allocatable :: history_offset(:)
        real(dp), allocatable :: adams_yh(:, :)
        real(dp), allocatable :: adams_saved_acor(:)
        type(bdf_state) :: bdf
        real(dp) :: h_abs = 0.0_dp
        real(dp) :: adams_yh_h = 0.0_dp
        integer :: mode = mode_adams
        integer :: history_count = 0
        integer :: adams_order = adams_history_min
        integer :: adams_method_order = adams_method_min
        integer :: adams_max_order = adams_method_max
        integer :: bdf_max_order_allowed = bdf_max_order
        integer :: adams_yh_order = 0
        integer :: adams_ialth = 2
        integer :: steps_at_order = 0
        integer :: low_stiff_count = 0
        integer :: high_stiff_count = 0
        integer :: switch_countdown = switch_wait_steps
        integer :: lband = -1
        integer :: uband = -1
        real(dp) :: min_step = 0.0_dp
        real(dp) :: adams_pdlast = 0.0_dp
        real(dp) :: adams_rmax = adams_initial_rmax
        real(dp) :: last_step_h = 0.0_dp
        real(dp) :: last_switch_time = 0.0_dp
        integer :: last_step_order = 1
        integer :: last_method_used = 1
        logical :: force_history_build = .false.
        logical :: adams_yh_valid = .false.
        logical :: adams_saved_acor_valid = .false.
    end type lsoda_state

    public :: lsoda_adams_dense_state
    public :: lsoda_advance
    public :: lsoda_initialize
    public :: lsoda_jacobian

contains

    subroutine lsoda_initialize(state, y, f, h_abs, min_step, lband, uband, &
            max_order_adams, max_order_bdf, t0)
        type(lsoda_state), intent(out) :: state !! initialized automatic-switching solver state
        real(dp), intent(in) :: y(:) !! initial state
        real(dp), intent(in) :: f(:) !! derivative at the initial state
        real(dp), intent(in) :: h_abs !! selected initial step magnitude
        real(dp), intent(in), optional :: min_step !! minimum retry step magnitude; default 0
        integer, intent(in), optional :: lband !! lower Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: uband !! upper Jacobian bandwidth; omitted means dense
        integer, intent(in), optional :: max_order_adams !! maximum Adams order, 1 through 12
        integer, intent(in), optional :: max_order_bdf !! maximum BDF order, 1 through 5
        real(dp), intent(in), optional :: t0 !! initial independent-variable value for diagnostics

        state%h_abs = h_abs
        allocate(state%f_hist(size(y), adams_history_max))
        allocate(state%history_offset(adams_history_max))
        allocate(state%adams_yh(size(y), adams_method_max + 1))
        allocate(state%adams_saved_acor(size(y)))
        state%f_hist = 0.0_dp
        state%history_offset = 0.0_dp
        state%adams_yh = 0.0_dp
        state%adams_saved_acor = 0.0_dp
        state%f_hist(:, 1) = f
        state%adams_yh(:, 1) = y
        state%mode = mode_adams
        state%history_count = 1
        state%adams_order = adams_history_min
        state%adams_method_order = adams_method_min
        state%adams_max_order = adams_method_max
        state%bdf_max_order_allowed = bdf_max_order
        if (present(max_order_adams)) then
            state%adams_max_order = min(max(adams_method_min, max_order_adams), &
                adams_method_max)
        end if
        if (present(max_order_bdf)) then
            state%bdf_max_order_allowed = min(max(1, max_order_bdf), bdf_max_order)
        end if
        state%adams_yh_order = 0
        state%adams_ialth = 2
        state%adams_yh_h = 0.0_dp
        state%steps_at_order = 0
        state%low_stiff_count = 0
        state%high_stiff_count = 0
        state%switch_countdown = switch_wait_steps
        state%lband = -1
        state%uband = -1
        state%min_step = 0.0_dp
        state%adams_pdlast = 0.0_dp
        state%adams_rmax = adams_initial_rmax
        state%last_step_h = 0.0_dp
        state%last_switch_time = 0.0_dp
        if (present(t0)) state%last_switch_time = t0
        state%last_step_order = adams_method_min
        state%last_method_used = mode_adams
        if (present(min_step)) state%min_step = min_step
        if (present(lband) .or. present(uband)) then
            state%lband = 0
            state%uband = 0
            if (present(lband)) state%lband = lband
            if (present(uband)) state%uband = uband
        end if
        state%force_history_build = .false.
        state%adams_yh_valid = .false.
        state%adams_saved_acor_valid = .false.
    end subroutine lsoda_initialize

    subroutine lsoda_advance(state, fun, t, y, f, t_bound, direction, rtol, atol, &
            max_step, y_new, f_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
            dense_valid, nfev, njev, nlu, attempt_limit, attempts, rejected, status, jac, &
            atol_vector, rtol_vector)
        type(lsoda_state), intent(inout) :: state !! mutable Adams/BDF switching state
        procedure(lsoda_rhs) :: fun !! right-hand side callback
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
        real(dp), intent(out) :: dense_d(:, :) !! BDF dense coefficients when stiff branch is used
        real(dp), intent(out) :: dense_h !! signed BDF dense spacing; zero for Adams branch
        integer, intent(out) :: dense_order !! BDF dense order; zero for Adams branch
        real(dp), intent(out) :: dense_anchor !! BDF dense endpoint anchor
        logical, intent(out) :: dense_valid !! true when dense_d describes this accepted interval
        integer, intent(inout) :: nfev !! accumulated right-hand-side evaluation count
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        integer, intent(inout) :: nlu !! accumulated LU factorization count
        integer, intent(in) :: attempt_limit !! maximum trial steps available to this call
        integer, intent(out) :: attempts !! trial steps consumed by this accepted step
        integer, intent(out) :: rejected !! rejected trial-step count
        integer, intent(out) :: status !! LSODA compatibility status code
        procedure(lsoda_jacobian), optional :: jac !! optional analytic Jacobian callback
        real(dp), intent(in), optional :: atol_vector(:) !! per-component absolute tolerances
        real(dp), intent(in), optional :: rtol_vector(:) !! per-component relative tolerances

        real(dp) :: atol_use(size(y))
        real(dp) :: rtol_use(size(y))
        integer :: bdf_attempts
        integer :: bdf_rejected
        integer :: bdf_status

        atol_use = atol
        if (present(atol_vector)) atol_use = atol_vector
        rtol_use = rtol
        if (present(rtol_vector)) rtol_use = rtol_vector

        attempts = 0
        rejected = 0
        status = lsoda_status_success
        dense_d = 0.0_dp
        dense_h = 0.0_dp
        dense_order = 0
        dense_anchor = t
        dense_valid = .false.

        if (state%mode == mode_bdf) then
            if (present(jac)) then
                if (state%lband >= 0) then
                    call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, &
                        max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                        nfev, njev, nlu, attempt_limit, bdf_attempts, bdf_rejected, &
                        bdf_status, jac, min_step_user=state%min_step, &
                        lband=state%lband, uband=state%uband, f_new_out=f_new, atol_vector=atol_use, rtol_vector=rtol_use)
                else
                    call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, &
                        max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                        nfev, njev, nlu, attempt_limit, bdf_attempts, bdf_rejected, &
                        bdf_status, jac, min_step_user=state%min_step, f_new_out=f_new, &
                        atol_vector=atol_use, rtol_vector=rtol_use)
                end if
            else if (state%lband >= 0) then
                call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, &
                    max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                    nfev, njev, nlu, attempt_limit, bdf_attempts, bdf_rejected, &
                    bdf_status, min_step_user=state%min_step, lband=state%lband, &
                    uband=state%uband, f_new_out=f_new, atol_vector=atol_use, rtol_vector=rtol_use)
            else
                call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, &
                    max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                    nfev, njev, nlu, attempt_limit, bdf_attempts, bdf_rejected, &
                    bdf_status, min_step_user=state%min_step, f_new_out=f_new, atol_vector=atol_use, rtol_vector=rtol_use)
            end if
            attempts = bdf_attempts
            rejected = bdf_rejected
            if (bdf_status /= bdf_status_success) then
                status = map_bdf_status(bdf_status)
                return
            end if
            if (.not. all(ieee_is_finite(f_new))) then
                status = lsoda_status_nonfinite
                return
            end if
            dense_valid = .true.
            state%last_step_h = t_new - t
            state%last_step_order = state%bdf%last_step_order
            state%last_method_used = mode_bdf
            state%h_abs = state%bdf%h_abs
            call consider_switch_to_adams(state, f_new, direction, t_new)
            return
        end if

        if (present(jac)) then
            call adams_advance(state, fun, t, y, f, t_bound, direction, rtol_use, atol_use, &
                max_step, y_new, f_new, t_new, dense_d, dense_h, dense_order, &
                dense_anchor, dense_valid, nfev, njev, attempt_limit, attempts, &
                rejected, status, jac)
        else
            call adams_advance(state, fun, t, y, f, t_bound, direction, rtol_use, atol_use, &
                max_step, y_new, f_new, t_new, dense_d, dense_h, dense_order, &
                dense_anchor, dense_valid, nfev, njev, attempt_limit, attempts, &
                rejected, status)
        end if
        if (status == lsoda_status_success) then
            state%last_step_h = t_new - t
            if (dense_order < 0) then
                state%last_step_order = -dense_order
            else
                state%last_step_order = max(adams_method_min, state%adams_method_order)
            end if
            state%last_method_used = mode_adams
        end if
        if (status /= lsoda_status_switch_bdf) return

        if (attempts >= attempt_limit) then
            status = lsoda_status_max_attempts
            return
        end if
        if (present(jac)) then
            if (state%lband >= 0) then
                call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, &
                    max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                    nfev, njev, nlu, attempt_limit - attempts, bdf_attempts, bdf_rejected, &
                    bdf_status, jac, min_step_user=state%min_step, &
                    lband=state%lband, uband=state%uband, f_new_out=f_new, atol_vector=atol_use, rtol_vector=rtol_use)
            else
                call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, &
                    max_step, y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
                    nfev, njev, nlu, attempt_limit - attempts, bdf_attempts, bdf_rejected, &
                    bdf_status, jac, min_step_user=state%min_step, f_new_out=f_new, &
                    atol_vector=atol_use, rtol_vector=rtol_use)
            end if
        else if (state%lband >= 0) then
            call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, max_step, &
                y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, nfev, njev, &
                nlu, attempt_limit - attempts, bdf_attempts, bdf_rejected, bdf_status, &
                min_step_user=state%min_step, lband=state%lband, uband=state%uband, &
                f_new_out=f_new, atol_vector=atol_use, rtol_vector=rtol_use)
        else
            call bdf_advance(state%bdf, fun, t, t_bound, direction, rtol, atol, max_step, &
                y_new, t_new, dense_d, dense_h, dense_order, dense_anchor, nfev, njev, &
                nlu, attempt_limit - attempts, bdf_attempts, bdf_rejected, bdf_status, &
                min_step_user=state%min_step, f_new_out=f_new, atol_vector=atol_use, rtol_vector=rtol_use)
        end if
        attempts = attempts + bdf_attempts
        rejected = rejected + bdf_rejected
        status = map_bdf_status(bdf_status)
        if (status /= lsoda_status_success) return
        if (.not. all(ieee_is_finite(f_new))) then
            status = lsoda_status_nonfinite
            return
        end if
        dense_valid = .true.
        state%last_step_h = t_new - t
        state%last_step_order = state%bdf%last_step_order
        state%last_method_used = mode_bdf
        state%h_abs = state%bdf%h_abs
        call consider_switch_to_adams(state, f_new, direction, t_new)
    end subroutine lsoda_advance

    subroutine adams_advance(state, fun, t, y, f, t_bound, direction, rtol, atol, &
            max_step, y_new, f_new, t_new, dense_d, dense_h, dense_order, dense_anchor, &
            dense_valid, nfev, njev, attempt_limit, attempts, rejected, status, jac)
        type(lsoda_state), intent(inout) :: state !! mutable Adams predictor-corrector state
        procedure(lsoda_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: y(:) !! current state
        real(dp), intent(in) :: f(:) !! derivative at the current state
        real(dp), intent(in) :: t_bound !! final integration bound
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: rtol(:) !! componentwise relative tolerances
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        real(dp), intent(out) :: y_new(:) !! accepted state
        real(dp), intent(out) :: f_new(:) !! derivative retained at the accepted state
        real(dp), intent(out) :: t_new !! accepted independent-variable value
        real(dp), intent(out) :: dense_d(:, :) !! accepted Adams Nordsieck dense coefficients
        real(dp), intent(out) :: dense_h !! signed scaling step for Adams dense output
        integer, intent(out) :: dense_order !! negative Adams order when dense output is valid
        real(dp), intent(out) :: dense_anchor !! accepted Adams endpoint for dense output
        logical, intent(out) :: dense_valid !! true when dense_d contains Adams Nordsieck data
        integer, intent(inout) :: nfev !! accumulated right-hand-side evaluation count
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        integer, intent(in) :: attempt_limit !! maximum trial steps available to this call
        integer, intent(out) :: attempts !! trial steps consumed by this accepted step
        integer, intent(out) :: rejected !! rejected trial-step count
        integer, intent(out) :: status !! LSODA compatibility status code
        procedure(lsoda_jacobian), optional :: jac !! optional analytic Jacobian callback

        real(dp) :: acor(size(y))
        real(dp) :: err(size(y))
        real(dp) :: scale(size(y))
        real(dp) :: yh_accepted(size(y), adams_method_max + 1)
        real(dp) :: h_abs
        real(dp) :: h
        real(dp) :: min_step
        real(dp) :: retry_min_step
        real(dp) :: err_norm
        real(dp) :: factor
        real(dp) :: proposed_h_abs
        real(dp) :: stiffness_index
        real(dp) :: stiffness_rate
        real(dp) :: local_lipschitz
        real(dp) :: stability_h_abs
        real(dp) :: stability_rate
        real(dp) :: weighted_solution_norm
        real(dp) :: bdf_switch_factor
        integer :: order
        logical :: accepted
        logical :: corrector_converged
        logical :: roundoff_polluted
        logical :: stability_limited

        attempts = 0
        rejected = 0
        status = lsoda_status_success
        dense_d = 0.0_dp
        dense_h = 0.0_dp
        dense_order = 0
        dense_anchor = t
        dense_valid = .false.
        h_abs = min(state%h_abs, max_step)
        min_step = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(t))
        retry_min_step = max(min_step, state%min_step)

        do
            if (attempts >= attempt_limit) then
                status = lsoda_status_max_attempts
                return
            end if
            if (h_abs < min_step) then
                status = lsoda_status_too_small_step
                return
            end if

            h_abs = min(h_abs, max_step, abs(t_bound - t))
            h = direction * h_abs
            t_new = t + h
            attempts = attempts + 1

            order = min(max(state%adams_method_order, adams_method_min), state%adams_max_order)
            if (.not. state%adams_yh_valid .or. state%adams_yh_order /= order) then
                call initialize_adams_nordsieck(state, y, f, h, order)
                order = state%adams_method_order
            else
                call rescale_adams_nordsieck(state, h)
            end if
            if (.not. state%adams_yh_valid) then
                status = lsoda_status_nonfinite
                return
            end if

            call adams_nordsieck_functional_step(fun, t, y, h, state%adams_yh, order, &
                rtol, atol, y_new, f_new, acor, err, yh_accepted, local_lipschitz, &
                corrector_converged, nfev)

            if (.not. corrector_converged) then
                rejected = rejected + 1
                factor = 0.25_dp
                if (h_abs * factor < retry_min_step) then
                    status = lsoda_status_too_small_step
                    return
                end if
                h_abs = h_abs * factor
                state%h_abs = h_abs
                state%adams_rmax = adams_failure_rmax
                state%adams_ialth = order + 1
                cycle
            end if

            scale = max(atol + rtol * abs(y), tiny(1.0_dp))
            err_norm = odepack_max_norm(err / scale)
            if (.not. all(ieee_is_finite(y_new)) .or. &
                    .not. all(ieee_is_finite(f_new))) then
                err_norm = huge(1.0_dp)
            end if

            accepted = err_norm <= 1.0_dp
            if (accepted) then
                if (err_norm == 0.0_dp) then
                    factor = state%adams_rmax
                else
                    factor = safety * err_norm**(-1.0_dp / real(order + 1, dp))
                end if
                if (rejected > 0) factor = min(1.0_dp, factor)
                factor = min(state%adams_rmax, max(0.5_dp, factor))
                proposed_h_abs = min(max_step, h_abs * factor)

                if (local_lipschitz > 0.0_dp) then
                    stiffness_rate = local_lipschitz
                    state%adams_pdlast = local_lipschitz
                    stiffness_index = h_abs * stiffness_rate
                else
                    stiffness_index = secant_stiffness_index(y, f, y_new, f_new, h_abs)
                    stiffness_rate = stiffness_index / max(h_abs, tiny(1.0_dp))
                end if
                ! Adams stability restrictions use PDLAST, the most recent
                ! nonzero corrector-rate estimate, not a secant fallback.
                stability_rate = state%adams_pdlast
                stability_limited = .false.
                if (stability_rate > 0.0_dp) then
                    stability_h_abs = adams_stability_limit(order) / stability_rate
                    if (proposed_h_abs > stability_h_abs) then
                        proposed_h_abs = max(h_abs * min_factor, stability_h_abs)
                        stability_limited = .true.
                    end if
                end if

                if (stiffness_index > stiff_switch_threshold .or. stability_limited) then
                    state%high_stiff_count = state%high_stiff_count + 1
                else
                    state%high_stiff_count = 0
                end if
                state%switch_countdown = state%switch_countdown - 1

                call update_history(state, f_new, h)
                state%adams_yh = yh_accepted
                state%adams_yh_h = h
                state%adams_yh_order = order
                state%adams_yh_valid = .true.
                if (size(dense_d, 2) >= order + 1) then
                    dense_d = 0.0_dp
                    dense_d(:, 1:order + 1) = yh_accepted(:, 1:order + 1)
                    dense_h = h
                    dense_order = -order
                    dense_anchor = t_new
                    dense_valid = .true.
                end if

                if (state%switch_countdown < 0 .and. &
                        order <= adams_switch_method_order_max) then
                    weighted_solution_norm = odepack_max_norm(y_new / scale)
                    roundoff_polluted = err_norm <= 100.0_dp * epsilon(1.0_dp) * &
                        max(weighted_solution_norm, 1.0_dp) .or. state%adams_pdlast <= 0.0_dp
                    if (roundoff_polluted) then
                        if (stability_limited) then
                            state%h_abs = min(max_step, 2.0_dp * h_abs)
                            if (present(jac)) then
                                call switch_to_bdf(state, fun, t_new, y_new, f_new, &
                                    direction, rtol, atol, njev, 2.0_dp, order, jac)
                            else
                                call switch_to_bdf(state, fun, t_new, y_new, f_new, &
                                    direction, rtol, atol, njev, 2.0_dp, order)
                            end if
                            return
                        end if
                    else if (adams_to_bdf_advantage(err_norm, h_abs, stability_rate, &
                            order, bdf_switch_factor)) then
                        if (present(jac)) then
                            call switch_to_bdf(state, fun, t_new, y_new, f_new, direction, &
                                rtol, atol, njev, bdf_switch_factor, order, jac)
                        else
                            call switch_to_bdf(state, fun, t_new, y_new, f_new, direction, &
                                rtol, atol, njev, bdf_switch_factor, order)
                        end if
                        return
                    end if
                end if

                call select_adams_step_order(state, acor, err_norm, scale, h_abs, &
                    max_step, stability_rate, rejected)
                return
            end if

            rejected = rejected + 1
            if (ieee_is_finite(local_lipschitz) .and. local_lipschitz > 0.0_dp) then
                state%adams_pdlast = local_lipschitz
            end if
            if (rejected >= adams_failure_restart_count) then
                if (0.1_dp * h_abs < retry_min_step) then
                    status = lsoda_status_too_small_step
                    return
                end if
                h_abs = 0.1_dp * h_abs
                state%h_abs = h_abs
                call reset_history(state, f)
                call initialize_adams_nordsieck(state, y, f, direction * h_abs, &
                    adams_method_min)
                cycle
            end if

            stability_rate = state%adams_pdlast
            if (.not. ieee_is_finite(stability_rate) .or. stability_rate < 0.0_dp) then
                stability_rate = 0.0_dp
            end if
            state%adams_rmax = adams_failure_rmax
            call select_rejected_adams_step_order(state, err_norm, scale, h_abs, &
                stability_rate, rejected, factor)
            if (h_abs * factor < retry_min_step) then
                status = lsoda_status_too_small_step
                return
            end if
            h_abs = h_abs * factor
            state%h_abs = h_abs
        end do
    end subroutine adams_advance



    subroutine initialize_adams_nordsieck(state, y, f, h, requested_order)
        type(lsoda_state), intent(inout) :: state !! Adams state to initialize at current t
        real(dp), intent(in) :: y(:) !! current solution value
        real(dp), intent(in) :: f(:) !! current derivative value
        real(dp), intent(in) :: h !! signed step scaling for the Nordsieck columns
        integer, intent(in) :: requested_order !! requested Adams method order, 1 through 12

        integer :: order
        logical :: ok

        order = min(max(requested_order, adams_method_min), state%adams_max_order)
        state%adams_yh = 0.0_dp
        state%adams_yh(:, 1) = y
        if (order == 1) then
            state%adams_yh(:, 2) = h * f
            ok = all(ieee_is_finite(state%adams_yh(:, 1:2)))
        else if (order <= state%history_count) then
            call construct_adams_nordsieck(y, state%f_hist, state%history_offset, &
                order, h, state%adams_yh, ok)
        else
            order = min(state%history_count, min(adams_history_max, state%adams_max_order))
            order = max(order, adams_method_min)
            call construct_adams_nordsieck(y, state%f_hist, state%history_offset, &
                order, h, state%adams_yh, ok)
        end if
        if (.not. ok) then
            state%adams_yh = 0.0_dp
            state%adams_yh(:, 1) = y
            state%adams_yh(:, 2) = h * f
            order = adams_method_min
            ok = all(ieee_is_finite(state%adams_yh(:, 1:2)))
        end if
        state%adams_method_order = order
        state%adams_order = max(adams_history_min, order - 1)
        state%adams_yh_order = order
        state%adams_yh_h = h
        state%adams_yh_valid = ok
        state%adams_ialth = max(2, order + 1)
        state%adams_saved_acor_valid = .false.
    end subroutine initialize_adams_nordsieck

    subroutine adams_nordsieck_functional_step(fun, t, y, h, yh, order, rtol, atol, &
            y_new, f_new, acor, err, yh_new, local_lipschitz, converged, nfev)
        procedure(lsoda_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: y(:) !! state at current time
        real(dp), intent(in) :: h !! signed trial step size
        real(dp), intent(in) :: yh(:, :) !! scaled derivatives at current time
        integer, intent(in) :: order !! Adams method order, 1 through 12
        real(dp), intent(in) :: rtol(:) !! componentwise relative tolerances
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        real(dp), intent(out) :: y_new(:) !! corrected state at t+h
        real(dp), intent(out) :: f_new(:) !! derivative from the converged iteration
        real(dp), intent(out) :: acor(:) !! accumulated Nordsieck correction
        real(dp), intent(out) :: err(:) !! one-step local-error estimate
        real(dp), intent(out) :: yh_new(:, :) !! corrected Nordsieck state at t+h
        real(dp), intent(out) :: local_lipschitz !! directional Lipschitz estimate
        logical, intent(out) :: converged !! true when functional iteration converged
        integer, intent(inout) :: nfev !! accumulated right-hand-side evaluation count

        real(dp) :: acor_old(size(y))
        real(dp) :: acor_trial(size(y))
        real(dp) :: conit
        real(dp) :: crate
        real(dp) :: dcon
        real(dp) :: del
        real(dp) :: delp
        real(dp) :: el(adams_method_max + 1)
        real(dp) :: pnorm
        real(dp) :: rate
        real(dp) :: scale(size(y))
        real(dp) :: tesco(3)
        real(dp) :: y_iter(size(y))
        real(dp) :: yh_predict(size(y), adams_method_max + 1)
        integer :: iteration
        integer :: j

        call predict_adams_nordsieck(yh, order, yh_predict)
        call odepack_adams_coefficients(order, el, tesco)
        scale = max(atol + rtol * abs(y), tiny(1.0_dp))
        pnorm = odepack_max_norm(yh_predict(:, 1) / scale)
        conit = 0.5_dp / real(order + 2, dp)
        crate = 0.7_dp
        rate = 0.0_dp
        delp = 0.0_dp
        acor = 0.0_dp
        acor_old = 0.0_dp
        y_iter = yh_predict(:, 1)
        y_new = y_iter
        f_new = 0.0_dp
        converged = .false.

        do iteration = 0, 2
            call fun(t + h, y_iter, f_new)
            nfev = nfev + 1
            acor_trial = h * f_new - yh_predict(:, 2)
            del = odepack_max_norm((acor_trial - acor_old) / scale)
            y_new = yh_predict(:, 1) + el(1) * acor_trial

            if (del <= 100.0_dp * pnorm * epsilon(1.0_dp)) then
                acor = acor_trial
                converged = .true.
                exit
            end if

            if (iteration > 0) then
                if (delp > 0.0_dp) then
                    rate = min(1024.0_dp, del / delp)
                else
                    rate = 1024.0_dp
                end if
                crate = max(0.2_dp * crate, rate)
                dcon = del * min(1.0_dp, 1.5_dp * crate) / (tesco(2) * conit)
                if (dcon <= 1.0_dp) then
                    acor = acor_trial
                    converged = .true.
                    exit
                end if
            end if

            if (iteration == 2) exit
            if (iteration >= 1 .and. del > 2.0_dp * delp) exit
            delp = del
            acor_old = acor_trial
            y_iter = y_new
        end do

        yh_new = yh_predict
        err = acor / max(tesco(2), tiny(1.0_dp))
        if (converged) then
            do j = 1, order + 1
                yh_new(:, j) = yh_new(:, j) + el(j) * acor
            end do
        end if
        if (rate > 0.0_dp .and. abs(h * el(1)) > tiny(1.0_dp)) then
            local_lipschitz = rate / abs(h * el(1))
        else
            ! DSTODA updates PDEST/PDLAST only from the measured corrector
            ! convergence rate.  In particular, a roundoff-level first
            ! correction leaves the Lipschitz estimate unchanged.
            local_lipschitz = 0.0_dp
        end if
    end subroutine adams_nordsieck_functional_step

    subroutine select_adams_step_order(state, acor, err_norm, scale, h_abs, &
            max_step, stiffness_rate, rejected)
        type(lsoda_state), intent(inout) :: state !! accepted Adams Nordsieck state
        real(dp), intent(in) :: acor(:) !! accumulated correction on the accepted step
        real(dp), intent(in) :: err_norm !! accepted ODEPACK-scaled local-error norm
        real(dp), intent(in) :: scale(:) !! componentwise error scale at the endpoint
        real(dp), intent(in) :: h_abs !! magnitude of the accepted step
        real(dp), intent(in) :: max_step !! maximum permitted next-step magnitude
        real(dp), intent(in) :: stiffness_rate !! current Adams stability-rate estimate
        integer, intent(in) :: rejected !! rejected attempts preceding this accepted step

        real(dp) :: ddn
        real(dp) :: dup
        real(dp) :: el(adams_method_max + 1)
        real(dp) :: pdh
        real(dp) :: rh
        real(dp) :: rhdn
        real(dp) :: rhsm
        real(dp) :: rhup
        real(dp) :: tesco(3)
        integer :: new_order
        integer :: order
        logical :: stability_forced

        order = state%adams_method_order
        state%adams_ialth = state%adams_ialth - 1
        if (state%adams_ialth > 0) then
            state%h_abs = min(max_step, h_abs)
            state%steps_at_order = state%steps_at_order + 1
            if (state%adams_ialth == 1 .and. order < state%adams_max_order) then
                state%adams_saved_acor = acor
                state%adams_saved_acor_valid = .true.
            end if
            return
        end if

        call odepack_adams_coefficients(order, el, tesco)
        rhsm = ideal_method_step_factor(max(err_norm, tiny(1.0_dp)), order)
        rhdn = 0.0_dp
        if (order > adams_method_min .and. tesco(1) > 0.0_dp) then
            ddn = odepack_max_norm(state%adams_yh(:, order + 1) / scale) / tesco(1)
            rhdn = 1.0_dp / &
                (1.3_dp * max(ddn, tiny(1.0_dp))**(1.0_dp / real(order, dp)) + 1.3e-6_dp)
        end if
        rhup = 0.0_dp
        if (order < state%adams_max_order .and. state%adams_saved_acor_valid .and. &
                tesco(3) > 0.0_dp) then
            dup = odepack_max_norm((acor - state%adams_saved_acor) / scale) / tesco(3)
            rhup = 1.0_dp / (1.4_dp * max(dup, tiny(1.0_dp))** &
                (1.0_dp / real(order + 2, dp)) + 1.4e-6_dp)
        end if

        pdh = max(h_abs * stiffness_rate, 1.0e-6_dp)
        if (order < state%adams_max_order) then
            rhup = min(rhup, adams_stability_limit(order + 1) / pdh)
        end if
        rhsm = min(rhsm, adams_stability_limit(order) / pdh)
        if (order > adams_method_min) then
            rhdn = min(rhdn, adams_stability_limit(order - 1) / pdh)
        end if

        new_order = order
        rh = rhsm
        if (rhup > rh .and. rhup > rhdn) then
            new_order = order + 1
            rh = rhup
        else if (rhdn > rh) then
            new_order = order - 1
            rh = rhdn
        end if
        if (rejected > 0) rh = min(1.0_dp, rh)
        rh = min(state%adams_rmax, max(min_factor, rh))
        stability_forced = rh * pdh * 1.00001_dp >= adams_stability_limit(new_order)

        if (rh < 1.1_dp .and. .not. stability_forced) then
            state%h_abs = min(max_step, h_abs)
            state%adams_ialth = 3
            state%adams_saved_acor_valid = .false.
            state%steps_at_order = state%steps_at_order + 1
            return
        end if

        if (new_order == order + 1) then
            state%adams_yh(:, new_order + 1) = &
                acor * el(order + 1) / real(order + 1, dp)
        else if (new_order < order) then
            state%adams_yh(:, new_order + 2:) = 0.0_dp
        end if
        state%adams_method_order = new_order
        state%adams_order = max(adams_history_min, new_order - 1)
        state%adams_yh_order = new_order
        state%adams_ialth = new_order + 1
        state%adams_saved_acor_valid = .false.
        if (new_order /= order) then
            state%steps_at_order = 0
        else
            state%steps_at_order = state%steps_at_order + 1
        end if
        state%h_abs = min(max_step, h_abs * rh)
        state%adams_rmax = adams_regular_rmax
    end subroutine select_adams_step_order

    subroutine select_rejected_adams_step_order(state, err_norm, scale, h_abs, &
            stiffness_rate, failure_count, factor)
        type(lsoda_state), intent(inout) :: state !! rejected Adams Nordsieck state
        real(dp), intent(in) :: err_norm !! failed ODEPACK-scaled local-error norm
        real(dp), intent(in) :: scale(:) !! componentwise error scale for the failed trial
        real(dp), intent(in) :: h_abs !! magnitude of the failed trial step
        real(dp), intent(in) :: stiffness_rate !! latest Adams stability-rate estimate
        integer, intent(in) :: failure_count !! consecutive error-test failures in this call
        real(dp), intent(out) :: factor !! multiplier for the next retry step magnitude

        real(dp) :: ddn
        real(dp) :: pdh
        real(dp) :: rhdn
        real(dp) :: rhsm
        real(dp) :: tesco(3)
        real(dp) :: unused_el(adams_method_max + 1)
        integer :: new_order
        integer :: order

        order = state%adams_method_order
        call odepack_adams_coefficients(order, unused_el, tesco)

        if (ieee_is_finite(err_norm) .and. err_norm > 0.0_dp) then
            rhsm = ideal_method_step_factor(err_norm, order)
        else
            rhsm = 0.0_dp
        end if
        rhdn = 0.0_dp
        if (order > adams_method_min .and. tesco(1) > 0.0_dp) then
            ddn = odepack_max_norm(state%adams_yh(:, order + 1) / scale) / tesco(1)
            if (ieee_is_finite(ddn)) then
                rhdn = 1.0_dp / (1.3_dp * max(ddn, tiny(1.0_dp))** &
                    (1.0_dp / real(order, dp)) + 1.3e-6_dp)
            end if
        end if

        pdh = max(h_abs * stiffness_rate, 1.0e-6_dp)
        rhsm = min(rhsm, adams_stability_limit(order) / pdh)
        if (order > adams_method_min) then
            rhdn = min(rhdn, adams_stability_limit(order - 1) / pdh)
        end if

        new_order = order
        factor = rhsm
        if (rhdn > rhsm) then
            new_order = order - 1
            factor = min(1.0_dp, rhdn)
        end if
        if (failure_count >= 2) factor = min(factor, 0.2_dp)
        factor = max(0.0_dp, factor)

        if (new_order < order) then
            state%adams_method_order = new_order
            state%adams_order = max(adams_history_min, new_order - 1)
            state%adams_yh_order = new_order
            state%adams_yh(:, new_order + 2:) = 0.0_dp
            state%steps_at_order = 0
        end if
        state%adams_ialth = new_order + 1
        state%adams_saved_acor_valid = .false.
    end subroutine select_rejected_adams_step_order


    subroutine predict_adams_nordsieck(yh, order, yh_predict)
        real(dp), intent(in) :: yh(:, :) !! scaled derivatives at the current Adams endpoint
        integer, intent(in) :: order !! highest scaled derivative retained in yh
        real(dp), intent(out) :: yh_predict(:, :) !! Pascal-shifted state at one step ahead

        integer :: first_column
        integer :: j
        integer :: pass

        yh_predict = yh
        do pass = 1, order
            first_column = order - pass + 1
            do j = first_column, order
                yh_predict(:, j) = yh_predict(:, j) + yh_predict(:, j + 1)
            end do
        end do
    end subroutine predict_adams_nordsieck

    subroutine rescale_adams_nordsieck(state, h_new)
        type(lsoda_state), intent(inout) :: state !! Adams Nordsieck state to rescale in place
        real(dp), intent(in) :: h_new !! new signed step represented by the scaled derivatives

        real(dp) :: ratio
        real(dp) :: ratio_power
        integer :: j

        if (.not. state%adams_yh_valid) return
        if (state%adams_yh_h == 0.0_dp .or. h_new == 0.0_dp) then
            state%adams_yh_valid = .false.
            return
        end if
        ratio = h_new / state%adams_yh_h
        if (.not. ieee_is_finite(ratio)) then
            state%adams_yh_valid = .false.
            return
        end if
        ratio_power = 1.0_dp
        do j = 1, state%adams_yh_order
            ratio_power = ratio_power * ratio
            state%adams_yh(:, j + 1) = ratio_power * state%adams_yh(:, j + 1)
        end do
        state%adams_yh_h = h_new
        if (.not. all(ieee_is_finite( &
                state%adams_yh(:, 1:state%adams_yh_order + 1)))) then
            state%adams_yh_valid = .false.
        end if
    end subroutine rescale_adams_nordsieck


    subroutine construct_adams_nordsieck(y_current, f_hist, history_offset, order, &
            h_reference, yh, ok)
        real(dp), intent(in) :: y_current(:) !! state at the derivative-history origin
        real(dp), intent(in) :: f_hist(:, :) !! newest-to-oldest derivative samples
        real(dp), intent(in) :: history_offset(:) !! sample times relative to current time
        integer, intent(in) :: order !! number of derivative samples used, 1 through 11
        real(dp), intent(in) :: h_reference !! signed step used to scale derivative columns
        real(dp), intent(out) :: yh(:, :) !! reconstructed scaled-derivative state
        logical, intent(out) :: ok !! true when interpolation produced finite columns

        real(dp) :: basis(0:adams_history_max - 1)
        real(dp) :: next_basis(0:adams_history_max - 1)
        real(dp) :: nodes(adams_history_max)
        real(dp) :: denominator
        integer :: degree
        integer :: i
        integer :: j
        integer :: k

        ok = .false.
        yh = 0.0_dp
        yh(:, 1) = y_current
        if (order < 1 .or. order > min(size(f_hist, 2), adams_history_max)) return
        if (h_reference == 0.0_dp .or. .not. ieee_is_finite(h_reference)) return

        nodes = 0.0_dp
        nodes(1:order) = history_offset(1:order) / h_reference
        do j = 1, order
            basis = 0.0_dp
            basis(0) = 1.0_dp
            degree = 0
            do i = 1, order
                if (i == j) cycle
                denominator = nodes(j) - nodes(i)
                if (denominator == 0.0_dp .or. .not. ieee_is_finite(denominator)) return
                next_basis = 0.0_dp
                do k = 0, degree
                    next_basis(k) = next_basis(k) - nodes(i) * basis(k) / denominator
                    next_basis(k + 1) = next_basis(k + 1) + basis(k) / denominator
                end do
                basis = next_basis
                degree = degree + 1
            end do
            do k = 0, order - 1
                yh(:, k + 2) = yh(:, k + 2) + &
                    h_reference * basis(k) * f_hist(:, j) / real(k + 1, dp)
            end do
        end do
        if (.not. all(ieee_is_finite(yh(:, 1:order + 1)))) return
        ok = .true.
    end subroutine construct_adams_nordsieck

    subroutine odepack_adams_coefficients(order, el, tesco_current)
        integer, intent(in) :: order !! Adams method order, 1 through 12
        real(dp), intent(out) :: el(:) !! ODEPACK correction polynomial coefficients
        real(dp), intent(out) :: tesco_current(3) !! ODEPACK error-test constants at order

        real(dp) :: agamq
        real(dp) :: elco(adams_method_max + 1, adams_method_max)
        real(dp) :: fnqm1
        real(dp) :: pc(adams_method_max)
        real(dp) :: pint
        real(dp) :: ragq
        real(dp) :: rqfac
        real(dp) :: rq1fac
        real(dp) :: tesco(3, adams_method_max)
        real(dp) :: tsign
        real(dp) :: xpin
        integer :: i
        integer :: ib
        integer :: nq
        integer :: nqm1
        integer :: nqp1

        el = 0.0_dp
        tesco_current = 0.0_dp
        if (order < 1 .or. order > adams_method_max) return

        elco = 0.0_dp
        tesco = 0.0_dp
        pc = 0.0_dp
        elco(1, 1) = 1.0_dp
        elco(2, 1) = 1.0_dp
        tesco(1, 1) = 0.0_dp
        tesco(2, 1) = 2.0_dp
        tesco(1, 2) = 1.0_dp
        tesco(3, adams_method_max) = 0.0_dp
        pc(1) = 1.0_dp
        rqfac = 1.0_dp
        do nq = 2, adams_method_max
            rq1fac = rqfac
            rqfac = rqfac / real(nq, dp)
            nqm1 = nq - 1
            fnqm1 = real(nqm1, dp)
            nqp1 = nq + 1
            pc(nq) = 0.0_dp
            do ib = 1, nqm1
                i = nqp1 - ib
                pc(i) = pc(i - 1) + fnqm1 * pc(i)
            end do
            pc(1) = fnqm1 * pc(1)
            pint = pc(1)
            xpin = pc(1) / 2.0_dp
            tsign = 1.0_dp
            do i = 2, nq
                tsign = -tsign
                pint = pint + tsign * pc(i) / real(i, dp)
                xpin = xpin + tsign * pc(i) / real(i + 1, dp)
            end do
            elco(1, nq) = pint * rq1fac
            elco(2, nq) = 1.0_dp
            do i = 2, nq
                elco(i + 1, nq) = rq1fac * pc(i) / real(i, dp)
            end do
            agamq = rqfac * xpin
            ragq = 1.0_dp / agamq
            tesco(2, nq) = ragq
            if (nq < adams_method_max) then
                tesco(1, nqp1) = ragq * rqfac / real(nqp1, dp)
            end if
            tesco(3, nqm1) = ragq
        end do

        el(1:order + 1) = elco(1:order + 1, order)
        tesco_current = tesco(:, order)
    end subroutine odepack_adams_coefficients



    subroutine switch_to_bdf(state, fun, t, y, f, direction, rtol, atol, njev, &
            step_factor, requested_order, jac)
        type(lsoda_state), intent(inout) :: state !! state to switch from Adams to BDF
        procedure(lsoda_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! current independent-variable value
        real(dp), intent(in) :: y(:) !! current state
        real(dp), intent(in) :: f(:) !! derivative at the current state
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: rtol(:) !! componentwise relative tolerances
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        integer, intent(inout) :: njev !! accumulated Jacobian evaluation count
        real(dp), intent(in) :: step_factor !! requested BDF step divided by accepted Adams step
        integer, intent(in) :: requested_order !! BDF order selected by the method-switch test
        procedure(lsoda_jacobian), optional :: jac !! optional analytic Jacobian callback

        real(dp) :: target_h_abs
        logical :: transfer_ok

        target_h_abs = abs(state%adams_yh_h) * &
            min(adams_regular_rmax, max(min_factor, step_factor))
        if (.not. ieee_is_finite(target_h_abs) .or. target_h_abs <= 0.0_dp) then
            target_h_abs = state%h_abs
        end if
        if (present(jac)) then
            if (state%lband >= 0) then
                call bdf_initialize(state%bdf, fun, t, y, f, direction, target_h_abs, &
                    maxval(rtol), atol(1), njev, jac, lband=state%lband, uband=state%uband, &
                    defer_jacobian=.true., atol_vector=atol, rtol_vector=rtol, &
                    max_order=state%bdf_max_order_allowed)
            else
                call bdf_initialize(state%bdf, fun, t, y, f, direction, target_h_abs, &
                    maxval(rtol), atol(1), njev, jac, defer_jacobian=.true., &
                    atol_vector=atol, rtol_vector=rtol, max_order=state%bdf_max_order_allowed)
            end if
        else if (state%lband >= 0) then
            call bdf_initialize(state%bdf, fun, t, y, f, direction, target_h_abs, maxval(rtol), &
                atol(1), njev, lband=state%lband, uband=state%uband, &
                defer_jacobian=.true., atol_vector=atol, rtol_vector=rtol, &
                max_order=state%bdf_max_order_allowed)
        else
            call bdf_initialize(state%bdf, fun, t, y, f, direction, target_h_abs, maxval(rtol), &
                atol(1), njev, defer_jacobian=.true., atol_vector=atol, &
                rtol_vector=rtol, max_order=state%bdf_max_order_allowed)
        end if
        call transfer_adams_state_to_bdf(state, direction, requested_order, &
            target_h_abs, transfer_ok)
        state%mode = mode_bdf
        state%last_switch_time = t
        state%h_abs = state%bdf%h_abs
        state%bdf%odepack_jacobian_refresh = .true.
        state%adams_yh_valid = .false.
        state%adams_yh_order = 0
        state%adams_yh_h = 0.0_dp
        state%adams_method_order = adams_method_min
        state%adams_order = adams_history_min
        state%adams_ialth = 2
        state%adams_saved_acor_valid = .false.
        state%low_stiff_count = 0
        state%high_stiff_count = 0
        state%adams_pdlast = 0.0_dp
        state%switch_countdown = switch_wait_steps
        state%steps_at_order = 0
    end subroutine switch_to_bdf

    subroutine transfer_adams_state_to_bdf(state, direction, requested_order, &
            target_h_abs, ok)
        type(lsoda_state), intent(inout) :: state !! solver state containing accepted Adams data
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        integer, intent(in) :: requested_order !! BDF order selected by the switch test
        real(dp), intent(in) :: target_h_abs !! desired BDF step magnitude after switching
        logical, intent(out) :: ok !! true when a finite BDF difference table was formed

        real(dp) :: samples(size(state%adams_yh, 1), 0:adams_switch_method_order_max)
        real(dp) :: differences(size(state%adams_yh, 1), 0:adams_switch_method_order_max)
        real(dp) :: ratio
        real(dp) :: x
        integer :: j
        integer :: k
        integer :: order

        ok = .false.
        if (.not. state%adams_yh_valid) return
        if (state%adams_yh_h == 0.0_dp) return
        if (.not. ieee_is_finite(state%adams_yh_h)) return
        if (target_h_abs <= 0.0_dp .or. .not. ieee_is_finite(target_h_abs)) return

        order = min(max(requested_order, 1), &
            min(state%bdf_max_order_allowed, state%adams_yh_order))
        ratio = direction * target_h_abs / state%adams_yh_h
        samples = 0.0_dp
        do j = 0, order
            x = -real(j, dp) * ratio
            samples(:, j) = state%adams_yh(:, order + 1)
            do k = order - 1, 0, -1
                samples(:, j) = samples(:, j) * x + state%adams_yh(:, k + 1)
            end do
        end do
        if (.not. all(ieee_is_finite(samples(:, 0:order)))) return

        state%bdf%odepack_yh = 0.0_dp
        state%bdf%odepack_yh(:, 1) = state%adams_yh(:, 1)
        x = 1.0_dp
        do k = 1, order
            x = x * ratio
            state%bdf%odepack_yh(:, k + 1) = x * state%adams_yh(:, k + 1)
        end do
        if (.not. all(ieee_is_finite(state%bdf%odepack_yh(:, 1:order + 1)))) return

        state%bdf%d = 0.0_dp
        differences = samples
        state%bdf%d(:, 1) = differences(:, 0)
        do k = 1, order
            do j = 0, order - k
                differences(:, j) = differences(:, j) - differences(:, j + 1)
            end do
            state%bdf%d(:, k + 1) = differences(:, 0)
        end do
        if (.not. all(ieee_is_finite(state%bdf%d(:, 1:order + 1)))) return

        state%bdf%h_abs = target_h_abs
        state%bdf%order = order
        state%bdf%last_step_order = order
        state%bdf%n_equal_steps = 0
        state%bdf%odepack_ialth = order + 1
        state%bdf%odepack_rmax = adams_regular_rmax
        state%bdf%odepack_crate = 0.7_dp
        state%bdf%odepack_saved_acor_valid = .false.
        state%bdf%odepack_yh_valid = .true.
        state%bdf%lu_valid = .false.
        ok = .true.
    end subroutine transfer_adams_state_to_bdf

    subroutine consider_switch_to_adams(state, f_new, direction, t_new)
        type(lsoda_state), intent(inout) :: state !! BDF state after an accepted step
        real(dp), intent(in) :: f_new(:) !! accepted derivative used to seed Adams history
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: t_new !! accepted time at which a method switch can occur

        real(dp) :: dm_adams
        real(dp) :: dsm
        real(dp) :: jacobian_index
        real(dp) :: rh_adams
        real(dp) :: rh_bdf
        real(dp) :: rh_stability
        integer :: order
        logical :: transfer_ok

        state%switch_countdown = state%switch_countdown - 1
        if (state%switch_countdown >= 0) return

        order = min(max(state%bdf%last_step_order, 1), &
            min(adams_switch_method_order_max, state%adams_max_order))
        dsm = max(state%bdf%last_error_norm, tiny(1.0_dp))
        rh_bdf = ideal_method_step_factor(dsm, order)
        dm_adams = dsm * odepack_bdf_error_constant(order) / &
            odepack_adams_error_constant(order)
        rh_adams = ideal_method_step_factor(dm_adams, order)
        jacobian_index = state%bdf%last_step_h_abs * matrix_inf_norm(state%bdf%jac)
        rh_stability = 2.0_dp * rh_adams
        if (jacobian_index * rh_adams > 1.0e-5_dp) then
            rh_stability = adams_stability_limit(order) / &
                max(jacobian_index, tiny(1.0_dp))
        end if
        rh_adams = min(rh_adams, rh_stability)

        if (rh_adams * method_switch_ratio >= 5.0_dp * rh_bdf) then
            state%mode = mode_adams
            state%last_switch_time = t_new
            state%h_abs = state%bdf%last_step_h_abs * rh_adams
            state%low_stiff_count = 0
            state%high_stiff_count = 0
            state%adams_pdlast = 0.0_dp
            state%switch_countdown = switch_wait_steps
            call reset_history(state, f_new)
            call transfer_bdf_state_to_adams(state, f_new, direction, order, transfer_ok)
            state%force_history_build = .not. transfer_ok
        end if
    end subroutine consider_switch_to_adams

    subroutine transfer_bdf_state_to_adams(state, f_current, direction, requested_order, ok)
        type(lsoda_state), intent(inout) :: state !! solver state containing current BDF data
        real(dp), intent(in) :: f_current(:) !! derivative at the BDF endpoint
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        integer, intent(in) :: requested_order !! Adams order selected by the switch test
        logical, intent(out) :: ok !! true when a finite Nordsieck transfer was formed

        real(dp) :: basis(0:adams_switch_method_order_max)
        real(dp) :: next_basis(0:adams_switch_method_order_max)
        real(dp) :: h_reference
        real(dp) :: scale
        integer :: degree
        integer :: j
        integer :: k
        integer :: order

        ok = .false.
        order = min(max(requested_order, adams_method_min), &
            min(state%adams_max_order, min(adams_switch_method_order_max, state%bdf%order)))
        if (order < adams_method_min) return
        if (.not. allocated(state%bdf%d)) return
        if (size(state%bdf%d, 2) < order + 1) return

        h_reference = direction * state%bdf%h_abs
        if (h_reference == 0.0_dp .or. .not. ieee_is_finite(h_reference)) return

        state%adams_yh = 0.0_dp
        state%adams_yh(:, 1) = state%bdf%d(:, 1)
        basis = 0.0_dp
        basis(0) = 1.0_dp
        do degree = 1, order
            next_basis = 0.0_dp
            scale = 1.0_dp / real(degree, dp)
            do k = 0, degree - 1
                next_basis(k) = next_basis(k) + &
                    real(degree - 1, dp) * scale * basis(k)
                next_basis(k + 1) = next_basis(k + 1) + scale * basis(k)
            end do
            basis = next_basis
            do j = 1, degree
                state%adams_yh(:, j + 1) = state%adams_yh(:, j + 1) + &
                    basis(j) * state%bdf%d(:, degree + 1)
            end do
        end do

        ! SciFort's stiff branch is the SciPy BDF/NDF recurrence rather than
        ! ODEPACK's pure BDF Nordsieck form.  Preserve the transferred higher
        ! derivatives, but make the Adams first-derivative column exactly h*f
        ! as required by its functional corrector.
        state%adams_yh(:, 2) = h_reference * f_current
        ok = all(ieee_is_finite(state%adams_yh(:, 1:order + 1)))
        if (.not. ok) then
            state%adams_yh = 0.0_dp
            return
        end if

        state%adams_method_order = order
        state%adams_order = max(adams_history_min, order - 1)
        state%adams_yh_order = order
        state%adams_yh_h = h_reference
        state%adams_yh_valid = .true.
        state%adams_ialth = order + 1
        state%adams_saved_acor_valid = .false.
        state%steps_at_order = 0
    end subroutine transfer_bdf_state_to_adams

    logical function adams_to_bdf_advantage(error_norm, h_abs, stiffness_rate, &
            method_order, bdf_factor) result(switch_method)
        real(dp), intent(in) :: error_norm !! accepted Adams local-error norm
        real(dp), intent(in) :: h_abs !! magnitude of the accepted Adams step
        real(dp), intent(in) :: stiffness_rate !! local Lipschitz or secant rate estimate
        integer, intent(in) :: method_order !! Adams-Moulton method order, 1 through 5
        real(dp), intent(out) :: bdf_factor !! ideal BDF step-size multiplier if switching

        real(dp) :: dm_bdf
        real(dp) :: pdh
        real(dp) :: rh_adams
        real(dp) :: rh_bdf
        real(dp) :: rh_stability

        switch_method = .false.
        bdf_factor = 1.0_dp
        if (method_order < 1 .or. method_order > adams_switch_method_order_max) return

        rh_adams = ideal_method_step_factor(max(error_norm, tiny(1.0_dp)), method_order)
        pdh = max(h_abs * stiffness_rate, 1.0e-6_dp)
        rh_stability = 2.0_dp * rh_adams
        if (pdh * rh_adams > 1.0e-5_dp) then
            rh_stability = adams_stability_limit(method_order) / pdh
        end if
        rh_adams = min(rh_adams, rh_stability)

        dm_bdf = max(error_norm, tiny(1.0_dp)) * &
            odepack_adams_error_constant(method_order) / &
            odepack_bdf_error_constant(method_order)
        rh_bdf = ideal_method_step_factor(dm_bdf, method_order)
        bdf_factor = min(adams_regular_rmax, max(min_factor, rh_bdf))
        switch_method = rh_bdf >= method_switch_ratio * rh_adams
    end function adams_to_bdf_advantage

    pure real(dp) function ideal_method_step_factor(error_norm, method_order) result(value)
        real(dp), intent(in) :: error_norm !! positive local-error norm for the method
        integer, intent(in) :: method_order !! local-error method order

        real(dp) :: exponent

        exponent = 1.0_dp / real(method_order + 1, dp)
        value = 1.0_dp / (1.2_dp * max(error_norm, tiny(1.0_dp))**exponent + 1.2e-6_dp)
    end function ideal_method_step_factor

    subroutine update_history(state, f_new, h)
        type(lsoda_state), intent(inout) :: state !! Adams derivative/time history to update
        real(dp), intent(in) :: f_new(:) !! derivative at the newly accepted endpoint
        real(dp), intent(in) :: h !! signed accepted step from the old to the new endpoint

        integer :: j
        integer :: new_count

        new_count = min(adams_history_max, state%history_count + 1)
        do j = new_count, 2, -1
            state%f_hist(:, j) = state%f_hist(:, j - 1)
            state%history_offset(j) = state%history_offset(j - 1) - h
        end do
        state%f_hist(:, 1) = f_new
        state%history_offset(1) = 0.0_dp
        state%history_count = new_count
    end subroutine update_history

    subroutine reset_history(state, f_current)
        type(lsoda_state), intent(inout) :: state !! Adams derivative/time history to reset
        real(dp), intent(in) :: f_current(:) !! derivative at the new history origin

        state%f_hist = 0.0_dp
        state%history_offset = 0.0_dp
        state%f_hist(:, 1) = f_current
        state%history_count = 1
        state%adams_order = adams_history_min
        state%adams_method_order = adams_method_min
        state%adams_ialth = 2
        state%adams_saved_acor_valid = .false.
        state%steps_at_order = 0
        state%adams_yh_valid = .false.
        state%adams_yh_order = 0
        state%adams_yh_h = 0.0_dp
    end subroutine reset_history


    pure real(dp) function secant_stiffness_index(y0, f0, y1, f1, h_abs) result(value)
        real(dp), intent(in) :: y0(:) !! state at the start of the secant interval
        real(dp), intent(in) :: f0(:) !! derivative at the start of the secant interval
        real(dp), intent(in) :: y1(:) !! finite state at the secant endpoint
        real(dp), intent(in) :: f1(:) !! derivative at the finite secant endpoint
        real(dp), intent(in) :: h_abs !! magnitude of the secant time interval

        real(dp) :: dy_norm
        real(dp) :: df_norm

        if (.not. all(ieee_is_finite(y1)) .or. .not. all(ieee_is_finite(f1))) then
            value = 0.0_dp
            return
        end if
        dy_norm = rms_norm(y1 - y0)
        df_norm = rms_norm(f1 - f0)
        if (dy_norm <= sqrt(tiny(1.0_dp)) * max(1.0_dp, rms_norm(y0))) then
            value = 0.0_dp
        else
            value = h_abs * df_norm / dy_norm
        end if
    end function secant_stiffness_index


    subroutine lsoda_adams_dense_state(t_anchor, h, order, yh, tq, yq)
        real(dp), intent(in) :: t_anchor !! accepted Adams endpoint represented by yh
        real(dp), intent(in) :: h !! signed scaling step used in yh
        integer, intent(in) :: order !! Adams Nordsieck polynomial order
        real(dp), intent(in) :: yh(:, :) !! endpoint scaled derivatives, columns 0:order
        real(dp), intent(in) :: tq !! interpolation point inside the accepted step
        real(dp), intent(out) :: yq(:) !! Nordsieck polynomial value at tq

        real(dp) :: theta
        integer :: j

        if (h == 0.0_dp .or. order < 0 .or. size(yh, 2) < order + 1) then
            yq = yh(:, 1)
            return
        end if
        theta = (tq - t_anchor) / h
        yq = yh(:, order + 1)
        do j = order, 1, -1
            yq = yh(:, j) + theta * yq
        end do
    end subroutine lsoda_adams_dense_state

    pure real(dp) function matrix_inf_norm(a) result(value)
        real(dp), intent(in) :: a(:, :) !! dense matrix whose infinity norm is requested

        integer :: i

        value = 0.0_dp
        do i = 1, size(a, 1)
            value = max(value, sum(abs(a(i, :))))
        end do
    end function matrix_inf_norm

    pure real(dp) function rms_norm(x) result(value)
        real(dp), intent(in) :: x(:) !! scaled vector whose RMS norm is requested

        value = sqrt(sum(x * x) / real(size(x), dp))
    end function rms_norm

    pure real(dp) function odepack_max_norm(x) result(value)
        real(dp), intent(in) :: x(:) !! tolerance-scaled vector for ODEPACK DMNORM semantics

        value = maxval(abs(x))
    end function odepack_max_norm

    pure integer function map_bdf_status(status) result(mapped)
        integer, intent(in) :: status !! BDF status code to translate

        select case (status)
        case (bdf_status_success)
            mapped = lsoda_status_success
        case (bdf_status_too_small_step)
            mapped = lsoda_status_too_small_step
        case (bdf_status_singular)
            mapped = lsoda_status_singular
        case (bdf_status_nonfinite)
            mapped = lsoda_status_nonfinite
        case (bdf_status_max_attempts)
            mapped = lsoda_status_max_attempts
        case default
            mapped = lsoda_status_nonfinite
        end select
    end function map_bdf_status

end module scifort_integrate_lsoda
