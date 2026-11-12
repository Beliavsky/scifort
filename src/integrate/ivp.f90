! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2022-2023, Jacob Williams
!
! Adaptive initial-value solver orchestration. The Bogacki-Shampine
! 3(2) and Dormand-Prince 5(4) tableaux are adapted from Jacob Williams'
! BSD-3-Clause rklib implementation. RK45 quartic dense-output coefficients are
! adapted from SciPy 1.17.0 under BSD-3-Clause. DOP853 uses the Hairer/Wanner
! tableau and seventh-order dense output adapted from Jacob Williams' permissively
! licensed modern Fortran DOP853. The solve_ivp driver, result/storage types, and
! event orchestration are SciFort code; event occurrence-limit semantics and
! per-event result grouping follow SciPy 1.17.0. Radau and BDF are adapted
! separately in radau.f90 and bdf.f90; LSODA compatibility is in lsoda.f90.
! See CODE_PROVENANCE.md.

module scifort_integrate_ivp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : error_unit, output_unit
    use scifort_kinds, only : dp
    use scifort_integrate_bdf, only : bdf_advance, bdf_dense_state, bdf_initialize, &
        bdf_jacobian, bdf_sparse_jacobian, bdf_state, bdf_status_invalid_jacobian, &
        bdf_status_max_attempts, bdf_status_nonfinite, bdf_status_singular, &
        bdf_status_success, bdf_status_too_small_step
    use scifort_integrate_dop853, only : dop853_dense_coefficients, &
        dop853_dense_state, dop853_step
    use scifort_integrate_lsoda, only : lsoda_adams_dense_state, lsoda_advance, &
        lsoda_initialize, lsoda_state, &
        lsoda_status_max_attempts, lsoda_status_nonfinite, lsoda_status_singular, &
        lsoda_status_success, lsoda_status_too_small_step
    use scifort_integrate_radau, only : radau_advance, radau_dense_state, &
        radau_initialize, radau_jacobian, radau_sparse_jacobian, radau_state, &
        radau_status_invalid_jacobian, radau_status_max_attempts, radau_status_nonfinite, &
        radau_status_singular, radau_status_success, radau_status_too_small_step
    use scifort_integrate_sparse_jacobian, only : ivp_sparse_matrix
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: ivp_status_success = 0
    integer, parameter, public :: ivp_status_event = 1
    integer, parameter, public :: ivp_status_failure = -1
    integer, parameter, public :: ivp_status_invalid_argument = -2
    integer, parameter, public :: ivp_status_max_steps = -3

    integer, parameter :: initial_capacity = 256
    integer, parameter :: default_max_steps = 100000
    real(dp), parameter :: safety = 0.9_dp
    real(dp), parameter :: min_factor = 0.2_dp
    real(dp), parameter :: max_factor = 10.0_dp

    abstract interface
        subroutine ode_rhs(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t
        end subroutine ode_rhs

        subroutine ode_rhs_vectorized(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable shared by all state columns
            real(dp), intent(in) :: y(:, :) !! state columns, shape (n,k)
            real(dp), intent(out) :: dydt(:, :) !! derivatives for each state column
        end subroutine ode_rhs_vectorized

        subroutine ode_jacobian(t, y, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: jac(:, :) !! Jacobian d fun / d y, shape (n,n)
        end subroutine ode_jacobian

        subroutine ode_sparse_jacobian(t, y, jac)
            import :: dp, ivp_sparse_matrix
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            type(ivp_sparse_matrix), intent(out) :: jac !! one-based CSC Jacobian d fun / d y
        end subroutine ode_sparse_jacobian

        function ode_event(t, y) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp) :: value
        end function ode_event

        subroutine ode_events_vector(t, y, values)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: values(:) !! event-function values at t
        end subroutine ode_events_vector
    end interface

    type, public :: ivp_event_set
        integer :: n_events = 0
        procedure(ode_events_vector), pointer, nopass :: fun => null()
        logical, allocatable :: terminal(:)
        real(dp), allocatable :: direction(:)
        integer, allocatable :: max_occurrences(:)
    end type ivp_event_set

    type, public :: ivp_dense_solution
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: y(:, :)
        real(dp), allocatable :: f(:, :)
        integer :: method_id = 0
        real(dp), allocatable :: rk45_q(:, :, :)
        real(dp), allocatable :: rk45_h(:)
        logical, allocatable :: rk45_valid(:)
        real(dp), allocatable :: dop853_f(:, :, :)
        real(dp), allocatable :: dop853_h(:)
        logical, allocatable :: dop853_valid(:)
        real(dp), allocatable :: bdf_d(:, :, :)
        real(dp), allocatable :: bdf_h(:)
        real(dp), allocatable :: bdf_anchor(:)
        integer, allocatable :: bdf_order(:)
        logical, allocatable :: bdf_valid(:)
        real(dp), allocatable :: radau_q(:, :, :)
        real(dp), allocatable :: radau_h(:)
        logical, allocatable :: radau_valid(:)
    contains
        procedure :: evaluate => ivp_dense_evaluate
    end type ivp_dense_solution

    type, public :: solve_ivp_result
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: y(:, :)
        real(dp), allocatable :: t_events(:)
        real(dp), allocatable :: y_events(:, :)
        integer, allocatable :: event_indices(:)
        integer, allocatable :: event_counts(:)
        integer :: n_event_functions = 0
        integer :: status = ivp_status_failure
        logical :: success = .false.
        integer :: nfev = 0
        integer :: njev = 0
        integer :: nlu = 0
        integer :: nsteps = 0
        integer :: naccepted = 0
        integer :: nrejected = 0
        real(dp), allocatable :: lsoda_hu(:)
        real(dp), allocatable :: lsoda_tcur(:)
        real(dp), allocatable :: lsoda_tolsf(:)
        real(dp), allocatable :: lsoda_tsw(:)
        integer, allocatable :: lsoda_nst(:)
        integer, allocatable :: lsoda_nfe(:)
        integer, allocatable :: lsoda_nje(:)
        integer, allocatable :: lsoda_nqu(:)
        integer, allocatable :: lsoda_mused(:)
        character(len=:), allocatable :: message
        type(ivp_dense_solution) :: sol
    contains
        procedure :: event_times_for => solve_ivp_event_times_for
        procedure :: event_states_for => solve_ivp_event_states_for
    end type solve_ivp_result

    public :: solve_ivp

contains

    function solve_ivp(fun, t_span, y0, method, t_eval, dense_output, events, &
            event_terminal, event_direction, rtol, atol, first_step, max_step, &
            max_steps, event_set, jac, event_max_occurrences, min_step, lband, uband, &
            jac_matrix, jac_sparsity, fun_vectorized, jac_sparse_matrix, jac_sparse, &
            lsoda_tcrit, lsoda_max_order_adams, lsoda_max_order_bdf, &
            lsoda_output_max_steps, lsoda_collect_output_info, lsoda_print_switches, &
            lsoda_max_hnil_messages) result(res)
        procedure(ode_rhs) :: fun !! right-hand side callback dy/dt = fun(t,y)
        real(dp), intent(in) :: t_span(2) !! initial and final independent-variable values
        real(dp), intent(in) :: y0(:) !! initial state at t_span(1)
        character(len=*), intent(in), optional :: &
            method !! RK23, RK45, DOP853, Radau, BDF, or LSODA; default RK45
        real(dp), intent(in), optional :: t_eval(:) !! sorted output points inside t_span
        logical, intent(in), optional :: dense_output !! retain continuous solution; default false
        procedure(ode_event), optional :: events !! one scalar zero-crossing event callback
        logical, intent(in), optional :: event_terminal !! stop at scalar event; default false
        real(dp), intent(in), optional :: event_direction !! scalar event direction; default 0
        real(dp), intent(in), optional :: &
            rtol(..) !! scalar or length-n relative tolerance; default 1e-3
        real(dp), intent(in), optional :: &
            atol(..) !! scalar or length-n absolute tolerance; default 1e-6
        real(dp), intent(in), optional :: first_step !! positive initial step magnitude
        real(dp), intent(in), optional :: max_step !! positive maximum step magnitude
        integer, intent(in), optional :: max_steps !! accepted/rejected step limit; default 100000
        type(ivp_event_set), intent(in), optional :: event_set !! vector event callback and controls
        procedure(ode_jacobian), optional :: &
            jac !! analytic Jacobian; LSODA uses packed storage when bands are supplied
        integer, intent(in), optional :: event_max_occurrences !! scalar event limit; 0 means unlimited
        real(dp), intent(in), optional :: min_step !! LSODA minimum retry step magnitude; default 0
        integer, intent(in), optional :: lband !! LSODA lower Jacobian bandwidth; default dense
        integer, intent(in), optional :: uband !! LSODA upper Jacobian bandwidth; default dense
        real(dp), intent(in), optional :: jac_matrix(:, :) !! constant Jacobian for Radau/BDF
        logical, intent(in), optional :: jac_sparsity(:, :) !! finite-difference Jacobian pattern
        procedure(ode_rhs_vectorized), optional :: &
            fun_vectorized !! batched RHS for BDF/Radau finite-difference Jacobians
        type(ivp_sparse_matrix), intent(in), optional :: &
            jac_sparse_matrix !! constant CSC Jacobian for BDF/Radau
        procedure(ode_sparse_jacobian), optional :: &
            jac_sparse !! state-dependent CSC Jacobian callback for BDF/Radau
        real(dp), intent(in), optional :: &
            lsoda_tcrit(:) !! LSODA hard critical-time sequence used by odeint compatibility
        integer, intent(in), optional :: &
            lsoda_max_order_adams !! LSODA Adams order cap, 1 through 12
        integer, intent(in), optional :: &
            lsoda_max_order_bdf !! LSODA BDF order cap, 1 through 5
        integer, intent(in), optional :: &
            lsoda_output_max_steps !! accepted LSODA-step cap before each t_eval output
        logical, intent(in), optional :: &
            lsoda_collect_output_info !! collect odeint-style LSODA diagnostics at t_eval outputs
        logical, intent(in), optional :: &
            lsoda_print_switches !! print LSODA Adams/BDF method switches for odeint compatibility
        integer, intent(in), optional :: &
            lsoda_max_hnil_messages !! maximum LSODA too-small-step warning messages; zero suppresses
        type(solve_ivp_result) :: res

        real(dp), allocatable :: t_store(:)
        real(dp), allocatable :: y_store(:, :)
        real(dp), allocatable :: f_store(:, :)
        real(dp), allocatable :: dense_store(:, :, :)
        real(dp), allocatable :: dense_h_store(:)
        real(dp), allocatable :: dense_anchor_store(:)
        integer, allocatable :: dense_order_store(:)
        logical, allocatable :: dense_valid_store(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: y_new(:)
        real(dp), allocatable :: f(:)
        real(dp), allocatable :: f_new(:)
        real(dp), allocatable :: err(:)
        real(dp), allocatable :: err_aux(:)
        real(dp), allocatable :: scale(:)
        real(dp), allocatable :: y_event(:)
        real(dp), allocatable :: event_times(:)
        real(dp), allocatable :: event_states(:, :)
        integer, allocatable :: event_indices_store(:)
        real(dp), allocatable :: g_old_vec(:)
        real(dp), allocatable :: g_new_vec(:)
        real(dp), allocatable :: event_dirs(:)
        integer, allocatable :: event_limits(:)
        integer, allocatable :: event_occurrences(:)
        real(dp), allocatable :: root_times(:)
        real(dp), allocatable :: root_states(:, :)
        integer, allocatable :: root_indices(:)
        real(dp), allocatable :: dop_k(:, :)
        real(dp), allocatable :: dense_step_coeff(:, :)
        type(bdf_state) :: bdf
        type(lsoda_state) :: lsoda
        type(radau_state) :: radau
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: t_bound
        real(dp) :: step_bound
        real(dp) :: direction
        real(dp) :: h_abs
        real(dp) :: h
        real(dp) :: hmax
        real(dp) :: min_step_use
        real(dp), allocatable :: rtol_use(:)
        real(dp) :: rtol_scalar_use
        real(dp), allocatable :: atol_use(:)
        real(dp) :: atol_scalar_use
        real(dp) :: err_norm
        real(dp) :: exponent
        real(dp) :: factor
        real(dp) :: machine_min_step
        real(dp) :: t_event
        real(dp) :: dense_h_step
        real(dp) :: dense_anchor_step
        integer :: dense_order_step
        integer :: method_id
        integer :: lband_use
        integer :: uband_use
        integer :: order_for_initial
        integer :: n
        integer :: capacity
        integer :: used
        integer :: event_capacity
        integer :: nevents
        integer :: n_event_functions
        integer :: steps_limit
        integer :: attempt
        integer :: iev
        integer :: ir
        integer :: nroots
        integer :: bdf_status
        integer :: bdf_attempts
        integer :: bdf_rejected
        integer :: lsoda_status
        integer :: lsoda_attempts
        integer :: lsoda_rejected
        integer :: radau_status
        integer :: radau_attempts
        integer :: radau_rejected
        integer :: current_event_index
        integer :: terminal_root_pos
        integer :: lsoda_tcrit_index
        integer :: lsoda_eval_index
        integer :: lsoda_budget_eval_index
        integer :: lsoda_call_steps
        integer :: lsoda_max_order_adams_use
        integer :: lsoda_max_order_bdf_use
        integer :: lsoda_output_max_steps_use
        integer :: lsoda_info_eval_index
        integer :: lsoda_mode_before
        integer :: lsoda_hnil_messages
        integer :: lsoda_max_hnil_messages_use
        logical :: keep_dense
        logical :: accepted
        logical :: have_scalar_event
        logical :: have_vector_events
        logical :: have_events
        logical :: event_hit
        logical :: event_failure
        logical :: root_ok
        logical :: last_rejected
        logical :: step_dense_valid
        logical :: need_stored_dense
        logical :: lsoda_dense_valid
        logical :: any_crossing
        logical :: lsoda_collect_output_info_use
        logical :: lsoda_print_switches_use

        call initialize_result(res)

        n = size(y0)
        if (n <= 0) then
            call set_failure(res, ivp_status_invalid_argument, 'y0 must be nonempty')
            return
        end if
        if (.not. all(ieee_is_finite(y0))) then
            call set_failure(res, ivp_status_invalid_argument, 'y0 must contain finite values')
            return
        end if
        if (.not. all(ieee_is_finite(t_span))) then
            call set_failure(res, ivp_status_invalid_argument, 't_span must be finite')
            return
        end if

        method_id = parse_method(method)
        if (method_id == 0) then
            call set_failure(res, ivp_status_invalid_argument, &
                'method must be RK23, RK45, DOP853, Radau, BDF, or LSODA')
            return
        end if

        if (present(fun_vectorized)) then
            if (method_id /= 15 .and. method_id /= 50) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'fun_vectorized is supported only for BDF and Radau')
                return
            end if
        end if

        if (merge(1, 0, present(jac)) + merge(1, 0, present(jac_matrix)) + &
                merge(1, 0, present(jac_sparse_matrix)) + merge(1, 0, present(jac_sparse)) > 1) then
            call set_failure(res, ivp_status_invalid_argument, &
                'use only one of jac, jac_matrix, jac_sparse_matrix, or jac_sparse')
            return
        end if
        if (present(jac_matrix)) then
            if (method_id /= 15 .and. method_id /= 50) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_matrix is supported only for BDF and Radau')
                return
            end if
            if (size(jac_matrix, 1) /= n .or. size(jac_matrix, 2) /= n) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_matrix must have shape (size(y0),size(y0))')
                return
            end if
            if (.not. all(ieee_is_finite(jac_matrix))) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_matrix must contain finite values')
                return
            end if
        end if
        if (present(jac_sparse_matrix)) then
            if (method_id /= 15 .and. method_id /= 50) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_sparse_matrix is supported only for BDF and Radau')
                return
            end if
            if (.not. jac_sparse_matrix%is_valid()) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_sparse_matrix must be a valid finite one-based CSC matrix')
                return
            end if
            if (jac_sparse_matrix%nrows /= n .or. jac_sparse_matrix%ncols /= n) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_sparse_matrix must have shape (size(y0),size(y0))')
                return
            end if
        end if
        if (present(jac_sparse)) then
            if (method_id /= 15 .and. method_id /= 50) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_sparse is supported only for BDF and Radau')
                return
            end if
        end if
        if (present(jac_sparsity)) then
            if (method_id /= 15 .and. method_id /= 50) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_sparsity is supported only for BDF and Radau')
                return
            end if
            if (size(jac_sparsity, 1) /= n .or. size(jac_sparsity, 2) /= n) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'jac_sparsity must have shape (size(y0),size(y0))')
                return
            end if
        end if

        allocate(rtol_use(n))
        rtol_use = 1.0e-3_dp
        if (present(rtol)) then
            select rank (rtol)
            rank (0)
                rtol_use = rtol
            rank (1)
                if (size(rtol) /= n) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'vector rtol must have size(y0) elements')
                    return
                end if
                rtol_use = rtol
            rank default
                call set_failure(res, ivp_status_invalid_argument, &
                    'rtol must be a scalar or rank-one array')
                return
            end select
        end if
        allocate(atol_use(n))
        atol_use = 1.0e-6_dp
        if (present(atol)) then
            select rank (atol)
            rank (0)
                atol_use = atol
            rank (1)
                if (size(atol) /= n) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'vector atol must have size(y0) elements')
                    return
                end if
                atol_use = atol
            rank default
                call set_failure(res, ivp_status_invalid_argument, &
                    'atol must be a scalar or rank-one array')
                return
            end select
        end if
        if (.not. all(ieee_is_finite(rtol_use)) .or. any(rtol_use <= 0.0_dp) .or. &
                .not. all(ieee_is_finite(atol_use)) .or. any(atol_use < 0.0_dp)) then
            call set_failure(res, ivp_status_invalid_argument, &
                'rtol must be positive and atol nonnegative')
            return
        end if
        rtol_use = max(rtol_use, 100.0_dp * epsilon(1.0_dp))
        rtol_scalar_use = maxval(rtol_use)
        atol_scalar_use = atol_use(1)

        hmax = huge(1.0_dp)
        if (present(max_step)) hmax = max_step
        if (ieee_is_nan(hmax) .or. hmax <= 0.0_dp) then
            call set_failure(res, ivp_status_invalid_argument, 'max_step must be positive')
            return
        end if

        min_step_use = 0.0_dp
        lband_use = -1
        uband_use = -1
        if (method_id == 12) then
            if (present(min_step)) min_step_use = min_step
            if (.not. ieee_is_finite(min_step_use) .or. min_step_use < 0.0_dp) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA min_step must be nonnegative and finite')
                return
            end if
            if (present(lband) .or. present(uband)) then
                lband_use = 0
                uband_use = 0
                if (present(lband)) lband_use = lband
                if (present(uband)) uband_use = uband
                if (lband_use < 0 .or. lband_use >= n .or. &
                        uband_use < 0 .or. uband_use >= n) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'LSODA lband and uband must be between 0 and size(y0)-1')
                    return
                end if
            end if
        end if

        lsoda_max_order_adams_use = 12
        if (present(lsoda_max_order_adams)) lsoda_max_order_adams_use = lsoda_max_order_adams
        lsoda_max_order_bdf_use = 5
        if (present(lsoda_max_order_bdf)) lsoda_max_order_bdf_use = lsoda_max_order_bdf
        lsoda_collect_output_info_use = .false.
        if (present(lsoda_collect_output_info)) then
            lsoda_collect_output_info_use = lsoda_collect_output_info
        end if
        lsoda_print_switches_use = .false.
        if (present(lsoda_print_switches)) lsoda_print_switches_use = lsoda_print_switches
        lsoda_max_hnil_messages_use = 0
        if (present(lsoda_max_hnil_messages)) then
            if (lsoda_max_hnil_messages < 0) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA maximum hnil-message count must be nonnegative')
                return
            end if
            lsoda_max_hnil_messages_use = lsoda_max_hnil_messages
        end if
        if (present(lsoda_max_order_adams) .or. present(lsoda_max_order_bdf) .or. &
                present(lsoda_tcrit) .or. present(lsoda_output_max_steps) .or. &
                lsoda_collect_output_info_use .or. present(lsoda_print_switches) .or. &
                present(lsoda_max_hnil_messages)) then
            if (method_id /= 12) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA-specific controls require method=LSODA')
                return
            end if
        end if
        if (lsoda_max_order_adams_use < 1 .or. lsoda_max_order_adams_use > 12) then
            call set_failure(res, ivp_status_invalid_argument, &
                'LSODA Adams maximum order must be between 1 and 12')
            return
        end if
        if (lsoda_max_order_bdf_use < 1 .or. lsoda_max_order_bdf_use > 5) then
            call set_failure(res, ivp_status_invalid_argument, &
                'LSODA BDF maximum order must be between 1 and 5')
            return
        end if
        lsoda_output_max_steps_use = 0
        if (present(lsoda_output_max_steps)) then
            if (.not. present(t_eval)) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA per-output step limits require t_eval')
                return
            end if
            if (lsoda_output_max_steps <= 0) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA per-output step limit must be positive')
                return
            end if
            lsoda_output_max_steps_use = lsoda_output_max_steps
        end if
        if (lsoda_collect_output_info_use .and. .not. present(t_eval)) then
            call set_failure(res, ivp_status_invalid_argument, &
                'LSODA output diagnostics require t_eval')
            return
        end if

        steps_limit = default_max_steps
        if (present(max_steps)) steps_limit = max_steps
        if (steps_limit <= 0) then
            call set_failure(res, ivp_status_invalid_argument, 'max_steps must be positive')
            return
        end if

        keep_dense = .false.
        if (present(dense_output)) keep_dense = dense_output

        have_scalar_event = present(events)
        have_vector_events = present(event_set)
        if (have_scalar_event .and. have_vector_events) then
            call set_failure(res, ivp_status_invalid_argument, &
                'use either events or event_set, not both')
            return
        end if
        have_events = have_scalar_event .or. have_vector_events

        if (have_scalar_event) then
            n_event_functions = 1
            allocate(g_old_vec(1), g_new_vec(1), event_dirs(1), &
                event_limits(1), event_occurrences(1))
            event_limits = 0
            if (present(event_terminal)) then
                if (event_terminal) event_limits(1) = 1
            end if
            if (present(event_max_occurrences)) then
                if (event_max_occurrences < 0) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'event_max_occurrences must be nonnegative')
                    return
                end if
                event_limits(1) = event_max_occurrences
            end if
            event_occurrences = 0
            event_dirs = 0.0_dp
            if (present(event_direction)) event_dirs(1) = event_direction
        else if (have_vector_events) then
            if (event_set%n_events <= 0 .or. .not. associated(event_set%fun)) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'event_set requires n_events > 0 and an associated callback')
                return
            end if
            n_event_functions = event_set%n_events
            allocate(g_old_vec(n_event_functions), g_new_vec(n_event_functions), &
                event_dirs(n_event_functions), &
                event_limits(n_event_functions), event_occurrences(n_event_functions))
            event_dirs = 0.0_dp
            event_limits = 0
            event_occurrences = 0
            if (allocated(event_set%direction)) then
                if (size(event_set%direction) /= n_event_functions) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'event_set direction must have n_events elements')
                    return
                end if
                event_dirs = event_set%direction
            end if
            if (allocated(event_set%terminal)) then
                if (size(event_set%terminal) /= n_event_functions) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'event_set terminal must have n_events elements')
                    return
                end if
                where (event_set%terminal)
                    event_limits = 1
                end where
            end if
            if (allocated(event_set%max_occurrences)) then
                if (size(event_set%max_occurrences) /= n_event_functions) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'event_set max_occurrences must have n_events elements')
                    return
                end if
                if (any(event_set%max_occurrences < 0)) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'event_set max_occurrences must be nonnegative')
                    return
                end if
                event_limits = event_set%max_occurrences
            end if
        else
            n_event_functions = 0
            allocate(g_old_vec(0), g_new_vec(0), event_dirs(0), &
                event_limits(0), event_occurrences(0))
        end if
        if (n_event_functions > 0) then
            if (.not. all(ieee_is_finite(event_dirs))) then
                call set_failure(res, ivp_status_invalid_argument, 'event directions must be finite')
                return
            end if
            where (event_dirs > 0.0_dp)
                event_dirs = 1.0_dp
            elsewhere (event_dirs < 0.0_dp)
                event_dirs = -1.0_dp
            end where
        end if

        t = t_span(1)
        t_bound = t_span(2)
        if (t_bound > t) then
            direction = 1.0_dp
        else if (t_bound < t) then
            direction = -1.0_dp
        else
            direction = 0.0_dp
        end if

        if (present(t_eval)) then
            if (.not. valid_t_eval(t_eval, t, t_bound, direction)) then
                call set_failure(res, ivp_status_invalid_argument, &
                    't_eval must be finite, sorted, and contained in t_span')
                return
            end if
        end if
        if (present(lsoda_tcrit)) then
            if (.not. present(t_eval)) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA critical times require t_eval')
                return
            end if
            if (.not. valid_ordered_points(lsoda_tcrit, direction)) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'LSODA critical times must be finite and monotone in the integration direction')
                return
            end if
        end if

        allocate(y(n), y_new(n), f(n), f_new(n), err(n), err_aux(n), scale(n), y_event(n))
        allocate(dop_k(n, 13))
        if (method_id == 12) then
            allocate(dense_step_coeff(n, 13))
        else
            allocate(dense_step_coeff(n, 7))
        end if
        dop_k = 0.0_dp
        dense_step_coeff = 0.0_dp
        y = y0
        call fun(t, y, f)
        res%nfev = 1
        if (.not. all(ieee_is_finite(f))) then
            call set_failure(res, ivp_status_failure, 'initial derivative is not finite')
            return
        end if

        if (lsoda_collect_output_info_use) then
            call initialize_lsoda_output_info(t)
        end if

        capacity = initial_capacity
        need_stored_dense = (method_id == 45 .or. method_id == 853 .or. method_id == 15 .or. &
            method_id == 50 .or. method_id == 12) .and. &
            (keep_dense .or. present(t_eval))
        allocate(t_store(capacity), y_store(n, capacity), f_store(n, capacity))
        if (need_stored_dense) then
            if (method_id == 12) then
                allocate(dense_store(n, 13, capacity))
            else
                allocate(dense_store(n, 7, capacity))
            end if
            allocate(dense_h_store(capacity), dense_anchor_store(capacity), &
                dense_order_store(capacity), dense_valid_store(capacity))
            dense_store = 0.0_dp
            dense_h_store = 0.0_dp
            dense_anchor_store = 0.0_dp
            dense_order_store = 0
            dense_valid_store = .false.
        else
            allocate(dense_store(0, 0, 0), dense_h_store(0), dense_anchor_store(0), &
                dense_order_store(0), dense_valid_store(0))
        end if
        used = 1
        t_store(1) = t
        y_store(:, 1) = y
        f_store(:, 1) = f

        event_capacity = 8
        allocate(event_times(event_capacity), event_states(n, event_capacity), &
            event_indices_store(event_capacity))
        nevents = 0
        if (have_events) then
            if (have_scalar_event) then
                g_old_vec(1) = events(t, y)
            else
                call event_set%fun(t, y, g_old_vec)
            end if
            if (.not. all(ieee_is_finite(g_old_vec))) then
                call set_failure(res, ivp_status_failure, 'initial event value is not finite')
                return
            end if
        end if

        if (direction == 0.0_dp) then
            if (have_events) then
                do iev = 1, n_event_functions
                    if (g_old_vec(iev) /= 0.0_dp) cycle
                    event_occurrences(iev) = event_occurrences(iev) + 1
                    call append_event(event_times, event_states, event_indices_store, &
                        event_capacity, nevents, t, y, iev)
                    if (event_limits(iev) > 0 .and. &
                            event_occurrences(iev) >= event_limits(iev)) then
                        res%status = ivp_status_event
                        res%success = .true.
                        res%message = 'A termination event occurred.'
                        exit
                    end if
                end do
            end if
            call finish_result(res, t_store, y_store, f_store, dense_store, &
                dense_h_store, dense_anchor_store, dense_order_store, dense_valid_store, method_id, used, t_eval, keep_dense, &
                event_times, event_states, event_indices_store, nevents, &
                n_event_functions, direction)
            if (.not. allocated(res%message)) then
                res%status = ivp_status_success
                res%success = .true.
                res%message = 'The solver successfully reached the end of the integration interval.'
            end if
            return
        end if

        lsoda_tcrit_index = 1
        lsoda_eval_index = 1
        lsoda_budget_eval_index = 1
        lsoda_info_eval_index = 2
        lsoda_call_steps = 0
        lsoda_hnil_messages = 0
        step_bound = t_bound
        if (method_id == 12 .and. present(lsoda_tcrit)) then
            call update_lsoda_task_bound(t, step_bound, lsoda_eval_index, lsoda_tcrit_index)
        end if
        if (method_id == 12 .and. present(lsoda_output_max_steps)) then
            call update_lsoda_output_budget(t, lsoda_budget_eval_index, lsoda_call_steps)
        end if

        if (method_id == 23) then
            exponent = -1.0_dp / 3.0_dp
            order_for_initial = 2
        else if (method_id == 45) then
            exponent = -1.0_dp / 5.0_dp
            order_for_initial = 4
        else if (method_id == 853) then
            exponent = -1.0_dp / 8.0_dp
            order_for_initial = 7
        else if (method_id == 12) then
            exponent = 0.0_dp
            order_for_initial = 1
        else if (method_id == 50) then
            exponent = 0.0_dp
            order_for_initial = 3
        else
            exponent = 0.0_dp
            order_for_initial = 1
        end if

        if (present(first_step)) then
            h_abs = first_step
            if (.not. ieee_is_finite(h_abs) .or. h_abs <= 0.0_dp) then
                call set_failure(res, ivp_status_invalid_argument, 'first_step must be positive and finite')
                return
            end if
            if (h_abs > abs(t_bound - t)) then
                call set_failure(res, ivp_status_invalid_argument, &
                    'first_step must not exceed the integration interval')
                return
            end if
            h_abs = min(h_abs, hmax)
        else if (method_id == 12) then
            call select_lsoda_initial_step(t, y, f, step_bound, rtol_use, atol_use, hmax, h_abs)
        else
            call select_initial_step(fun, t, y, f, direction, t_bound, order_for_initial, &
                rtol_use, atol_use, hmax, h_abs, res%nfev)
        end if

        if (method_id == 12) then
            if (lband_use >= 0) then
                call lsoda_initialize(lsoda, y, f, h_abs, min_step_use, &
                    lband=lband_use, uband=uband_use, &
                    max_order_adams=lsoda_max_order_adams_use, &
                    max_order_bdf=lsoda_max_order_bdf_use, t0=t)
            else
                call lsoda_initialize(lsoda, y, f, h_abs, min_step_use, &
                    max_order_adams=lsoda_max_order_adams_use, &
                    max_order_bdf=lsoda_max_order_bdf_use, t0=t)
            end if
            if (lsoda_collect_output_info_use) call update_lsoda_output_info(t)
        end if

        if (method_id == 15) then
            if (present(jac_sparse)) then
                call bdf_initialize(bdf, fun, t, y, f, direction, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized, &
                    jac_sparse=jac_sparse)
            else if (present(jac)) then
                call bdf_initialize(bdf, fun, t, y, f, direction, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, jac, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
            else if (present(jac_matrix)) then
                call bdf_initialize(bdf, fun, t, y, f, direction, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_matrix=jac_matrix, jac_sparsity=jac_sparsity, &
                    fun_vectorized=fun_vectorized)
            else if (present(jac_sparse_matrix)) then
                call bdf_initialize(bdf, fun, t, y, f, direction, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparse_matrix=jac_sparse_matrix, jac_sparsity=jac_sparsity, &
                    fun_vectorized=fun_vectorized)
            else
                call bdf_initialize(bdf, fun, t, y, f, direction, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
            end if
            if (bdf%jacobian_sparse) then
                if (.not. bdf%sparse_jac%is_valid() .or. bdf%sparse_jac%nrows /= n .or. &
                        bdf%sparse_jac%ncols /= n) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'initial BDF sparse Jacobian must be a valid finite n-by-n CSC matrix')
                    return
                end if
            else if (allocated(bdf%jac)) then
                if (.not. all(ieee_is_finite(bdf%jac))) then
                    call set_failure(res, ivp_status_failure, &
                        'initial BDF Jacobian is not finite')
                    return
                end if
            end if
        end if

        if (method_id == 50) then
            if (present(jac_sparse)) then
                call radau_initialize(radau, fun, t, y, f, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized, &
                    jac_sparse=jac_sparse)
            else if (present(jac)) then
                call radau_initialize(radau, fun, t, y, f, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, jac, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
            else if (present(jac_matrix)) then
                call radau_initialize(radau, fun, t, y, f, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_matrix=jac_matrix, jac_sparsity=jac_sparsity, &
                    fun_vectorized=fun_vectorized)
            else if (present(jac_sparse_matrix)) then
                call radau_initialize(radau, fun, t, y, f, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparse_matrix=jac_sparse_matrix, jac_sparsity=jac_sparsity, &
                    fun_vectorized=fun_vectorized)
            else
                call radau_initialize(radau, fun, t, y, f, h_abs, rtol_scalar_use, &
                    atol_scalar_use, res%njev, atol_vector=atol_use, rtol_vector=rtol_use, &
                    jac_sparsity=jac_sparsity, fun_vectorized=fun_vectorized)
            end if
            if (radau%jacobian_sparse) then
                if (.not. radau%sparse_jac%is_valid() .or. radau%sparse_jac%nrows /= n .or. &
                        radau%sparse_jac%ncols /= n) then
                    call set_failure(res, ivp_status_invalid_argument, &
                        'initial Radau sparse Jacobian must be a valid finite n-by-n CSC matrix')
                    return
                end if
            else if (allocated(radau%jac)) then
                if (.not. all(ieee_is_finite(radau%jac))) then
                    call set_failure(res, ivp_status_failure, &
                        'initial Radau Jacobian is not finite')
                    return
                end if
            end if
        end if

        allocate(root_times(max(1, n_event_functions)), &
            root_states(n, max(1, n_event_functions)), &
            root_indices(max(1, n_event_functions)))
        current_event_index = 1
        attempt = 0
        last_rejected = .false.
        do while (direction * (t_bound - t) > 0.0_dp)
            dense_step_coeff = 0.0_dp
            dense_h_step = 0.0_dp
            dense_anchor_step = t
            dense_order_step = 0
            step_dense_valid = .false.

            if (method_id == 12) then
                if (attempt >= steps_limit) then
                    call set_failure(res, ivp_status_max_steps, &
                        'maximum number of step attempts exceeded')
                    exit
                end if
                if (present(lsoda_output_max_steps)) then
                    if (lsoda_budget_eval_index <= size(t_eval) .and. &
                            lsoda_call_steps >= lsoda_output_max_steps_use) then
                        call set_failure(res, ivp_status_max_steps, &
                            'maximum number of LSODA steps exceeded before next output time')
                        exit
                    end if
                end if
                lsoda_mode_before = lsoda%mode
                if (present(jac)) then
                    call lsoda_advance(lsoda, fun, t, y, f, step_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, f_new, t_new, dense_step_coeff, &
                        dense_h_step, dense_order_step, dense_anchor_step, lsoda_dense_valid, &
                        res%nfev, res%njev, res%nlu, steps_limit - attempt, lsoda_attempts, &
                        lsoda_rejected, lsoda_status, jac, atol_vector=atol_use, rtol_vector=rtol_use)
                else
                    call lsoda_advance(lsoda, fun, t, y, f, step_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, f_new, t_new, dense_step_coeff, &
                        dense_h_step, dense_order_step, dense_anchor_step, lsoda_dense_valid, &
                        res%nfev, res%njev, res%nlu, steps_limit - attempt, lsoda_attempts, &
                        lsoda_rejected, lsoda_status, atol_vector=atol_use, rtol_vector=rtol_use)
                end if
                attempt = attempt + lsoda_attempts
                res%nrejected = res%nrejected + lsoda_rejected
                res%nsteps = res%naccepted + res%nrejected
                if (lsoda_status == lsoda_status_success .and. lsoda_print_switches_use) then
                    if (lsoda%mode /= lsoda_mode_before) then
                        call print_lsoda_method_switch(lsoda%last_switch_time, lsoda%mode)
                    end if
                end if
                if (lsoda_status /= lsoda_status_success) then
                    if (lsoda_status == lsoda_status_too_small_step .and. &
                            lsoda_hnil_messages < lsoda_max_hnil_messages_use) then
                        lsoda_hnil_messages = lsoda_hnil_messages + 1
                        call print_lsoda_hnil_warning(t, direction * lsoda%h_abs, &
                            lsoda_hnil_messages, lsoda_max_hnil_messages_use)
                    end if
                    select case (lsoda_status)
                    case (lsoda_status_max_attempts)
                        call set_failure(res, ivp_status_max_steps, &
                            'maximum number of step attempts exceeded')
                    case (lsoda_status_too_small_step)
                        call set_failure(res, ivp_status_failure, &
                            'required LSODA step size is too small')
                    case (lsoda_status_singular)
                        call set_failure(res, ivp_status_failure, &
                            'LSODA BDF Newton matrix is numerically singular')
                    case (lsoda_status_nonfinite)
                        call set_failure(res, ivp_status_failure, &
                            'LSODA right-hand side or Jacobian is not finite')
                    case default
                        call set_failure(res, ivp_status_failure, 'LSODA step failed')
                    end select
                    exit
                end if
                accepted = .true.
                h = t_new - t
                factor = 1.0_dp
                last_rejected = lsoda_rejected > 0
                step_dense_valid = lsoda_dense_valid
                if (present(lsoda_output_max_steps)) lsoda_call_steps = lsoda_call_steps + 1
            else if (method_id == 50) then
                if (attempt >= steps_limit) then
                    call set_failure(res, ivp_status_max_steps, &
                        'maximum number of step attempts exceeded')
                    exit
                end if
                if (present(jac_sparse)) then
                    call radau_advance(radau, fun, t, y, f, t_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, f_new, t_new, dense_step_coeff, &
                        dense_h_step, res%nfev, res%njev, res%nlu, steps_limit - attempt, &
                        radau_attempts, radau_rejected, radau_status, &
                        atol_vector=atol_use, rtol_vector=rtol_use, jac_sparsity=jac_sparsity, &
                        fun_vectorized=fun_vectorized, jac_sparse=jac_sparse)
                else if (present(jac)) then
                    call radau_advance(radau, fun, t, y, f, t_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, f_new, t_new, dense_step_coeff, &
                        dense_h_step, res%nfev, res%njev, res%nlu, steps_limit - attempt, &
                        radau_attempts, radau_rejected, radau_status, jac, &
                        atol_vector=atol_use, rtol_vector=rtol_use, jac_sparsity=jac_sparsity, &
                        fun_vectorized=fun_vectorized)
                else
                    call radau_advance(radau, fun, t, y, f, t_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, f_new, t_new, dense_step_coeff, &
                        dense_h_step, res%nfev, res%njev, res%nlu, steps_limit - attempt, &
                        radau_attempts, radau_rejected, radau_status, &
                        atol_vector=atol_use, rtol_vector=rtol_use, jac_sparsity=jac_sparsity, &
                        fun_vectorized=fun_vectorized)
                end if
                attempt = attempt + radau_attempts
                res%nrejected = res%nrejected + radau_rejected
                res%nsteps = res%naccepted + res%nrejected
                if (radau_status /= radau_status_success) then
                    select case (radau_status)
                    case (radau_status_max_attempts)
                        call set_failure(res, ivp_status_max_steps, &
                            'maximum number of step attempts exceeded')
                    case (radau_status_too_small_step)
                        call set_failure(res, ivp_status_failure, &
                            'required Radau step size is too small')
                    case (radau_status_singular)
                        call set_failure(res, ivp_status_failure, &
                            'Radau Newton matrix is numerically singular')
                    case (radau_status_nonfinite)
                        call set_failure(res, ivp_status_failure, &
                            'Radau right-hand side or Jacobian is not finite')
                    case (radau_status_invalid_jacobian)
                        call set_failure(res, ivp_status_failure, &
                            'Radau sparse Jacobian callback returned an invalid CSC matrix')
                    case default
                        call set_failure(res, ivp_status_failure, 'Radau step failed')
                    end select
                    exit
                end if
                accepted = .true.
                h = t_new - t
                factor = 1.0_dp
                last_rejected = radau_rejected > 0
                dense_anchor_step = t_new
                dense_order_step = 3
                step_dense_valid = .true.
            else if (method_id == 15) then
                if (attempt >= steps_limit) then
                    call set_failure(res, ivp_status_max_steps, &
                        'maximum number of step attempts exceeded')
                    exit
                end if
                if (present(jac_sparse)) then
                    call bdf_advance(bdf, fun, t, t_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, t_new, dense_step_coeff, dense_h_step, &
                        dense_order_step, dense_anchor_step, res%nfev, res%njev, res%nlu, &
                        steps_limit - attempt, bdf_attempts, bdf_rejected, bdf_status, &
                        atol_vector=atol_use, rtol_vector=rtol_use, jac_sparsity=jac_sparsity, &
                        fun_vectorized=fun_vectorized, jac_sparse=jac_sparse)
                else if (present(jac)) then
                    call bdf_advance(bdf, fun, t, t_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, t_new, dense_step_coeff, dense_h_step, &
                        dense_order_step, dense_anchor_step, res%nfev, res%njev, res%nlu, &
                        steps_limit - attempt, bdf_attempts, bdf_rejected, bdf_status, jac, &
                        atol_vector=atol_use, rtol_vector=rtol_use, jac_sparsity=jac_sparsity, &
                        fun_vectorized=fun_vectorized)
                else
                    call bdf_advance(bdf, fun, t, t_bound, direction, rtol_scalar_use, &
                        atol_scalar_use, hmax, y_new, t_new, dense_step_coeff, dense_h_step, &
                        dense_order_step, dense_anchor_step, res%nfev, res%njev, res%nlu, &
                        steps_limit - attempt, bdf_attempts, bdf_rejected, bdf_status, &
                        atol_vector=atol_use, rtol_vector=rtol_use, jac_sparsity=jac_sparsity, &
                        fun_vectorized=fun_vectorized)
                end if
                attempt = attempt + bdf_attempts
                res%nrejected = res%nrejected + bdf_rejected
                res%nsteps = res%naccepted + res%nrejected
                if (bdf_status /= bdf_status_success) then
                    select case (bdf_status)
                    case (bdf_status_max_attempts)
                        call set_failure(res, ivp_status_max_steps, &
                            'maximum number of step attempts exceeded')
                    case (bdf_status_too_small_step)
                        call set_failure(res, ivp_status_failure, &
                            'required BDF step size is too small')
                    case (bdf_status_singular)
                        call set_failure(res, ivp_status_failure, &
                            'BDF Newton matrix is numerically singular')
                    case (bdf_status_nonfinite)
                        call set_failure(res, ivp_status_failure, &
                            'BDF right-hand side or Jacobian is not finite')
                    case (bdf_status_invalid_jacobian)
                        call set_failure(res, ivp_status_failure, &
                            'BDF sparse Jacobian callback returned an invalid CSC matrix')
                    case default
                        call set_failure(res, ivp_status_failure, 'BDF step failed')
                    end select
                    exit
                end if
                accepted = .true.
                h = t_new - t
                f_new = 0.0_dp
                factor = 1.0_dp
                last_rejected = bdf_rejected > 0
                step_dense_valid = .true.
            else
                attempt = attempt + 1
                if (attempt > steps_limit) then
                    call set_failure(res, ivp_status_max_steps, &
                        'maximum number of step attempts exceeded')
                    exit
                end if

                machine_min_step = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(t))
                h_abs = min(h_abs, hmax, abs(t_bound - t))
                if (abs(t_bound - t) <= machine_min_step) then
                    h_abs = abs(t_bound - t)
                else if (h_abs < machine_min_step) then
                    call set_failure(res, ivp_status_failure, 'required step size is too small')
                    exit
                end if
                h = direction * h_abs

                if (method_id == 23) then
                    call rk23_step(fun, t, y, f, h, y_new, f_new, err, res%nfev)
                else if (method_id == 45) then
                    call rk45_step(fun, t, y, f, h, y_new, f_new, err, dense_step_coeff, &
                        res%nfev)
                else
                    call dop853_step(fun, t, y, f, h, y_new, f_new, err, err_aux, &
                        res%nfev, dop_k)
                end if

                if (.not. all(ieee_is_finite(y_new)) .or. &
                        .not. all(ieee_is_finite(f_new))) then
                    err_norm = huge(1.0_dp)
                else
                    scale = max(atol_use + rtol_use * max(abs(y), abs(y_new)), tiny(1.0_dp))
                    if (method_id == 853) then
                        err_norm = dop853_error_norm(err, err_aux, scale, h)
                    else
                        err_norm = rms_norm(err / scale)
                    end if
                end if

                accepted = err_norm <= 1.0_dp
                if (accepted) then
                    if (err_norm == 0.0_dp) then
                        factor = max_factor
                    else
                        factor = min(max_factor, safety * err_norm**exponent)
                    end if
                    if (last_rejected) factor = min(1.0_dp, factor)

                    t_new = t + h
                    dense_h_step = h
                    dense_anchor_step = t_new
                    step_dense_valid = method_id == 45
                    if (method_id == 853 .and. need_stored_dense) then
                        call dop853_dense_coefficients(fun, t, y, h, y_new, f_new, dop_k, &
                            dense_step_coeff, res%nfev)
                        if (.not. all(ieee_is_finite(dense_step_coeff))) then
                            call set_failure(res, ivp_status_failure, &
                                'DOP853 dense-output derivative is not finite')
                            exit
                        end if
                        step_dense_valid = .true.
                    end if
                end if
            end if

            if (accepted) then
                event_hit = .false.
                event_failure = .false.
                nroots = 0
                terminal_root_pos = 0
                if (have_events) then
                    if (have_scalar_event) then
                        g_new_vec(1) = events(t_new, y_new)
                    else
                        call event_set%fun(t_new, y_new, g_new_vec)
                    end if
                    if (.not. all(ieee_is_finite(g_new_vec))) then
                        call set_failure(res, ivp_status_failure, 'event value is not finite')
                        exit
                    end if
                    any_crossing = .false.
                    do iev = 1, n_event_functions
                        any_crossing = any_crossing .or. &
                            crossing_detected(g_old_vec(iev), g_new_vec(iev), event_dirs(iev))
                    end do
                    if (method_id == 853 .and. any_crossing .and. .not. step_dense_valid) then
                        call dop853_dense_coefficients(fun, t, y, h, y_new, f_new, dop_k, &
                            dense_step_coeff, res%nfev)
                        if (.not. all(ieee_is_finite(dense_step_coeff))) then
                            call set_failure(res, ivp_status_failure, &
                                'DOP853 dense-output derivative is not finite')
                            exit
                        end if
                        step_dense_valid = .true.
                    end if

                    do iev = 1, n_event_functions
                        if (.not. crossing_detected(g_old_vec(iev), g_new_vec(iev), &
                                event_dirs(iev))) cycle
                        current_event_index = iev
                        call locate_event(selected_event, method_id, t, y, f, t_new, y_new, &
                            f_new, dense_step_coeff, dense_h_step, dense_anchor_step, &
                            dense_order_step, step_dense_valid, g_old_vec(iev), &
                            g_new_vec(iev), t_event, y_event, root_ok)
                        if (.not. root_ok) then
                            call set_failure(res, ivp_status_failure, &
                                'event value is not finite during root location')
                            event_failure = .true.
                            exit
                        end if
                        nroots = nroots + 1
                        root_times(nroots) = t_event
                        root_states(:, nroots) = y_event
                        root_indices(nroots) = iev
                        event_occurrences(iev) = event_occurrences(iev) + 1
                    end do
                    if (event_failure) exit

                    call sort_event_roots(root_times, root_states, root_indices, nroots, direction)
                    do ir = 1, nroots
                        call append_event(event_times, event_states, event_indices_store, &
                            event_capacity, nevents, root_times(ir), root_states(:, ir), &
                            root_indices(ir))
                        if (event_limits(root_indices(ir)) > 0 .and. &
                                event_occurrences(root_indices(ir)) >= &
                                event_limits(root_indices(ir))) then
                            terminal_root_pos = ir
                            exit
                        end if
                    end do
                    event_hit = terminal_root_pos > 0
                    if (event_hit) then
                        t_new = root_times(terminal_root_pos)
                        y_new = root_states(:, terminal_root_pos)
                        if (method_id /= 15 .and. method_id /= 50) then
                            call fun(t_new, y_new, f_new)
                            res%nfev = res%nfev + 1
                        else
                            f_new = 0.0_dp
                        end if
                    end if
                end if

                call append_state(t_store, y_store, f_store, dense_store, dense_h_store, &
                    dense_anchor_store, dense_order_store, dense_valid_store, capacity, used, &
                    t_new, y_new, f_new, dense_step_coeff, dense_h_step, dense_anchor_step, &
                    dense_order_step, step_dense_valid, need_stored_dense)
                res%naccepted = res%naccepted + 1
                last_rejected = .false.
                res%nsteps = res%naccepted + res%nrejected
                t = t_new
                y = y_new
                f = f_new
                if (method_id == 12 .and. lsoda_collect_output_info_use) then
                    call update_lsoda_output_info(t)
                end if
                if (method_id == 12 .and. present(lsoda_tcrit)) then
                    call update_lsoda_task_bound(t, step_bound, lsoda_eval_index, &
                        lsoda_tcrit_index)
                end if
                if (method_id == 12 .and. present(lsoda_output_max_steps)) then
                    call update_lsoda_output_budget(t, lsoda_budget_eval_index, &
                        lsoda_call_steps)
                end if
                if (have_events .and. .not. event_hit) g_old_vec = g_new_vec
                if (method_id /= 15 .and. method_id /= 50) h_abs = h_abs * factor

                if (event_hit) then
                    res%status = ivp_status_event
                    res%success = .true.
                    res%message = 'A termination event occurred.'
                    exit
                end if
            else
                res%nrejected = res%nrejected + 1
                last_rejected = .true.
                res%nsteps = res%naccepted + res%nrejected
                if (.not. ieee_is_finite(err_norm)) then
                    factor = min_factor
                else
                    factor = max(min_factor, safety * err_norm**exponent)
                end if
                h_abs = h_abs * factor
            end if
        end do

        if (.not. allocated(res%message)) then
            res%status = ivp_status_success
            res%success = .true.
            res%message = 'The solver successfully reached the end of the integration interval.'
        end if

        call finish_result(res, t_store, y_store, f_store, dense_store, dense_h_store, &
            dense_anchor_store, dense_order_store, dense_valid_store, method_id, used, &
            t_eval, keep_dense, event_times, &
            event_states, event_indices_store, nevents, n_event_functions, direction)

    contains

        subroutine print_lsoda_method_switch(t_switch, new_mode)
            real(dp), intent(in) :: t_switch !! time at which LSODA selected a new method
            integer, intent(in) :: new_mode !! LSODA mode selected for the following step

            if (new_mode == 2) then
                write(output_unit, '(a,es24.16)') 'lsoda: switch to BDF at t = ', t_switch
            else
                write(output_unit, '(a,es24.16)') 'lsoda: switch to Adams at t = ', t_switch
            end if
        end subroutine print_lsoda_method_switch

        subroutine print_lsoda_hnil_warning(t_now, h_now, message_index, message_limit)
            real(dp), intent(in) :: t_now !! current LSODA time at the too-small-step failure
            real(dp), intent(in) :: h_now !! signed LSODA step that cannot make reliable progress
            integer, intent(in) :: message_index !! one-based warning number
            integer, intent(in) :: message_limit !! configured maximum number of warnings

            write(error_unit, '(a,es24.16,a,es24.16)') &
                'lsoda warning: step is too small to advance reliably at t = ', t_now, &
                ', h = ', h_now
            if (message_index == message_limit .and. message_limit > 1) then
                write(error_unit, '(a,i0,a)') 'lsoda warning: maximum of ', message_limit, &
                    ' hnil messages reached; further messages are suppressed'
            end if
        end subroutine print_lsoda_hnil_warning

        subroutine initialize_lsoda_output_info(t_initial)
            real(dp), intent(in) :: t_initial !! initial LSODA time for deterministic diagnostics

            integer :: ninfo

            ninfo = max(0, size(t_eval) - 1)
            if (allocated(res%lsoda_hu)) deallocate(res%lsoda_hu)
            if (allocated(res%lsoda_tcur)) deallocate(res%lsoda_tcur)
            if (allocated(res%lsoda_tolsf)) deallocate(res%lsoda_tolsf)
            if (allocated(res%lsoda_tsw)) deallocate(res%lsoda_tsw)
            if (allocated(res%lsoda_nst)) deallocate(res%lsoda_nst)
            if (allocated(res%lsoda_nfe)) deallocate(res%lsoda_nfe)
            if (allocated(res%lsoda_nje)) deallocate(res%lsoda_nje)
            if (allocated(res%lsoda_nqu)) deallocate(res%lsoda_nqu)
            if (allocated(res%lsoda_mused)) deallocate(res%lsoda_mused)
            allocate(res%lsoda_hu(ninfo), res%lsoda_tcur(ninfo), &
                res%lsoda_tolsf(ninfo), res%lsoda_tsw(ninfo), res%lsoda_nst(ninfo), &
                res%lsoda_nfe(ninfo), res%lsoda_nje(ninfo), res%lsoda_nqu(ninfo), &
                res%lsoda_mused(ninfo))
            res%lsoda_hu = 0.0_dp
            res%lsoda_tcur = t_initial
            res%lsoda_tolsf = 0.0_dp
            res%lsoda_tsw = t_initial
            res%lsoda_nst = 0
            res%lsoda_nfe = res%nfev
            res%lsoda_nje = 0
            res%lsoda_nqu = 1
            res%lsoda_mused = 1
        end subroutine initialize_lsoda_output_info

        subroutine update_lsoda_output_info(t_now)
            real(dp), intent(in) :: t_now !! current LSODA internal time after an accepted step

            integer :: info_index

            if (.not. lsoda_collect_output_info_use) return
            do while (lsoda_info_eval_index <= size(t_eval))
                if (direction * (t_eval(lsoda_info_eval_index) - t_now) > 0.0_dp) exit
                info_index = lsoda_info_eval_index - 1
                res%lsoda_hu(info_index) = lsoda%last_step_h
                res%lsoda_tcur(info_index) = t_now
                res%lsoda_tolsf(info_index) = 0.0_dp
                res%lsoda_tsw(info_index) = lsoda%last_switch_time
                res%lsoda_nst(info_index) = res%naccepted
                res%lsoda_nfe(info_index) = res%nfev
                res%lsoda_nje(info_index) = res%njev
                res%lsoda_nqu(info_index) = lsoda%last_step_order
                res%lsoda_mused(info_index) = lsoda%last_method_used
                lsoda_info_eval_index = lsoda_info_eval_index + 1
            end do
        end subroutine update_lsoda_output_info

        subroutine update_lsoda_output_budget(t_now, eval_index, call_steps)
            real(dp), intent(in) :: t_now !! current LSODA internal time
            integer, intent(inout) :: eval_index !! next odeint-style output index
            integer, intent(inout) :: call_steps !! accepted steps since the previous output call

            logical :: output_reached

            if (.not. present(t_eval)) return
            output_reached = .false.
            do while (eval_index <= size(t_eval))
                if (direction * (t_eval(eval_index) - t_now) > 0.0_dp) exit
                eval_index = eval_index + 1
                output_reached = .true.
            end do
            if (output_reached) call_steps = 0
        end subroutine update_lsoda_output_budget

        subroutine update_lsoda_task_bound(t_now, task_bound, eval_index, crit_index)
            real(dp), intent(in) :: t_now !! current LSODA time
            real(dp), intent(out) :: task_bound !! hard bound for the next LSODA internal step
            integer, intent(inout) :: eval_index !! next unevaluated t_eval index
            integer, intent(inout) :: crit_index !! current critical-time index

            task_bound = t_bound
            if (.not. present(t_eval) .or. .not. present(lsoda_tcrit)) return
            if (direction == 0.0_dp) return

            do while (eval_index <= size(t_eval))
                if (direction * (t_eval(eval_index) - t_now) > 0.0_dp) exit
                eval_index = eval_index + 1
            end do
            if (eval_index > size(t_eval)) return

            do while (crit_index <= size(lsoda_tcrit))
                if (direction * (t_eval(eval_index) - lsoda_tcrit(crit_index)) <= 0.0_dp) exit
                crit_index = crit_index + 1
            end do
            if (crit_index > size(lsoda_tcrit)) return

            if (direction * (lsoda_tcrit(crit_index) - t_now) > 0.0_dp .and. &
                    direction * (t_bound - lsoda_tcrit(crit_index)) >= 0.0_dp) then
                task_bound = lsoda_tcrit(crit_index)
            end if
        end subroutine update_lsoda_task_bound

        function selected_event(tq, yq) result(value)
            real(dp), intent(in) :: tq !! event-root trial time
            real(dp), intent(in) :: yq(:) !! interpolated trial state
            real(dp) :: value
            real(dp) :: values(max(1, n_event_functions))

            if (have_scalar_event) then
                value = events(tq, yq)
            else
                call event_set%fun(tq, yq, values(1:n_event_functions))
                value = values(current_event_index)
            end if
        end function selected_event

    end function solve_ivp

    subroutine rk23_step(fun, t, y, f0, h, ynew, fnew, err, nfev)
        procedure(ode_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! start time
        real(dp), intent(in) :: y(:) !! state at start time
        real(dp), intent(in) :: f0(:) !! derivative at start time
        real(dp), intent(in) :: h !! signed step size
        real(dp), intent(out) :: ynew(:) !! third-order state estimate
        real(dp), intent(out) :: fnew(:) !! derivative at the new state
        real(dp), intent(out) :: err(:) !! third-minus-second-order estimate
        integer, intent(inout) :: nfev !! accumulated right-hand-side evaluations

        real(dp) :: k2(size(y))
        real(dp) :: k3(size(y))
        real(dp) :: y2(size(y))

        call fun(t + 0.5_dp * h, y + h * (0.5_dp * f0), k2)
        call fun(t + 0.75_dp * h, y + h * (0.75_dp * k2), k3)
        ynew = y + h * ((2.0_dp / 9.0_dp) * f0 + (1.0_dp / 3.0_dp) * k2 + &
            (4.0_dp / 9.0_dp) * k3)
        call fun(t + h, ynew, fnew)
        y2 = y + h * ((7.0_dp / 24.0_dp) * f0 + 0.25_dp * k2 + &
            (1.0_dp / 3.0_dp) * k3 + 0.125_dp * fnew)
        err = ynew - y2
        nfev = nfev + 3
    end subroutine rk23_step

    subroutine rk45_step(fun, t, y, f0, h, ynew, fnew, err, dense_coeff, nfev)
        procedure(ode_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t !! start time
        real(dp), intent(in) :: y(:) !! state at start time
        real(dp), intent(in) :: f0(:) !! derivative at start time
        real(dp), intent(in) :: h !! signed step size
        real(dp), intent(out) :: ynew(:) !! fifth-order state estimate
        real(dp), intent(out) :: fnew(:) !! derivative at the new state
        real(dp), intent(out) :: err(:) !! fifth-minus-fourth-order estimate
        real(dp), intent(out) :: dense_coeff(:, :) !! quartic coefficients; shape (n,>=4)
        integer, intent(inout) :: nfev !! accumulated right-hand-side evaluations

        real(dp) :: k2(size(y))
        real(dp) :: k3(size(y))
        real(dp) :: k4(size(y))
        real(dp) :: k5(size(y))
        real(dp) :: k6(size(y))
        real(dp) :: y4(size(y))

        call fun(t + h / 5.0_dp, y + h * (f0 / 5.0_dp), k2)
        call fun(t + 3.0_dp * h / 10.0_dp, &
            y + h * ((3.0_dp / 40.0_dp) * f0 + (9.0_dp / 40.0_dp) * k2), k3)
        call fun(t + 4.0_dp * h / 5.0_dp, &
            y + h * ((44.0_dp / 45.0_dp) * f0 - (56.0_dp / 15.0_dp) * k2 + &
            (32.0_dp / 9.0_dp) * k3), k4)
        call fun(t + 8.0_dp * h / 9.0_dp, &
            y + h * ((19372.0_dp / 6561.0_dp) * f0 - (25360.0_dp / 2187.0_dp) * k2 + &
            (64448.0_dp / 6561.0_dp) * k3 - (212.0_dp / 729.0_dp) * k4), k5)
        call fun(t + h, y + h * ((9017.0_dp / 3168.0_dp) * f0 - &
            (355.0_dp / 33.0_dp) * k2 + (46732.0_dp / 5247.0_dp) * k3 + &
            (49.0_dp / 176.0_dp) * k4 - (5103.0_dp / 18656.0_dp) * k5), k6)
        ynew = y + h * ((35.0_dp / 384.0_dp) * f0 + (500.0_dp / 1113.0_dp) * k3 + &
            (125.0_dp / 192.0_dp) * k4 - (2187.0_dp / 6784.0_dp) * k5 + &
            (11.0_dp / 84.0_dp) * k6)
        call fun(t + h, ynew, fnew)
        y4 = y + h * ((5179.0_dp / 57600.0_dp) * f0 + &
            (7571.0_dp / 16695.0_dp) * k3 + (393.0_dp / 640.0_dp) * k4 - &
            (92097.0_dp / 339200.0_dp) * k5 + (187.0_dp / 2100.0_dp) * k6 + &
            (1.0_dp / 40.0_dp) * fnew)
        err = ynew - y4

        dense_coeff = 0.0_dp
        dense_coeff(:, 1) = f0
        dense_coeff(:, 2) = &
            (-8048581381.0_dp / 2820520608.0_dp) * f0 + &
            (131558114200.0_dp / 32700410799.0_dp) * k3 - &
            (1754552775.0_dp / 470086768.0_dp) * k4 + &
            (127303824393.0_dp / 49829197408.0_dp) * k5 - &
            (282668133.0_dp / 205662961.0_dp) * k6 + &
            (40617522.0_dp / 29380423.0_dp) * fnew
        dense_coeff(:, 3) = &
            (8663915743.0_dp / 2820520608.0_dp) * f0 - &
            (68118460800.0_dp / 10900136933.0_dp) * k3 + &
            (14199869525.0_dp / 1410260304.0_dp) * k4 - &
            (318862633887.0_dp / 49829197408.0_dp) * k5 + &
            (2019193451.0_dp / 616988883.0_dp) * k6 - &
            (110615467.0_dp / 29380423.0_dp) * fnew
        dense_coeff(:, 4) = &
            (-12715105075.0_dp / 11282082432.0_dp) * f0 + &
            (87487479700.0_dp / 32700410799.0_dp) * k3 - &
            (10690763975.0_dp / 1880347072.0_dp) * k4 + &
            (701980252875.0_dp / 199316789632.0_dp) * k5 - &
            (1453857185.0_dp / 822651844.0_dp) * k6 + &
            (69997945.0_dp / 29380423.0_dp) * fnew
        nfev = nfev + 6
    end subroutine rk45_step

    subroutine select_lsoda_initial_step(t0, y0, f0, t_bound, rtol, atol, &
            max_step, h_abs)
        real(dp), intent(in) :: t0 !! initial independent-variable value
        real(dp), intent(in) :: y0(:) !! initial state
        real(dp), intent(in) :: f0(:) !! derivative at the initial state
        real(dp), intent(in) :: t_bound !! final independent-variable value
        real(dp), intent(in) :: rtol(:) !! componentwise relative tolerances
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        real(dp), intent(out) :: h_abs !! ODEPACK-style initial LSODA step magnitude

        real(dp) :: scale(size(y0))
        real(dp) :: f_norm
        real(dp) :: interval
        real(dp) :: sum_term
        real(dp) :: tol
        real(dp) :: w0

        interval = abs(t_bound - t0)
        scale = max(atol + rtol * abs(y0), tiny(1.0_dp))
        tol = min(1.0e-3_dp, max(maxval(rtol), 100.0_dp * epsilon(1.0_dp)))
        f_norm = maxval(abs(f0) / scale)
        w0 = max(abs(t0), abs(t_bound))
        if (w0 <= 0.0_dp) then
            h_abs = min(interval, max_step)
            return
        end if
        sum_term = 1.0_dp / (tol * w0 * w0) + tol * f_norm * f_norm
        h_abs = min(1.0_dp / sqrt(sum_term), interval, max_step)
    end subroutine select_lsoda_initial_step

    subroutine select_initial_step(fun, t0, y0, f0, direction, t_bound, order, &
            rtol, atol, max_step, h_abs, nfev)
        procedure(ode_rhs) :: fun !! right-hand side callback
        real(dp), intent(in) :: t0 !! initial independent-variable value
        real(dp), intent(in) :: y0(:) !! initial state
        real(dp), intent(in) :: f0(:) !! derivative at the initial state
        real(dp), intent(in) :: direction !! integration direction, +1 or -1
        real(dp), intent(in) :: t_bound !! final independent-variable value
        integer, intent(in) :: order !! embedded error-estimator order
        real(dp), intent(in) :: rtol(:) !! componentwise relative tolerances
        real(dp), intent(in) :: atol(:) !! componentwise absolute tolerances
        real(dp), intent(in) :: max_step !! maximum permitted step magnitude
        real(dp), intent(out) :: h_abs !! selected initial step magnitude
        integer, intent(inout) :: nfev !! accumulated right-hand-side evaluations

        real(dp) :: scale(size(y0))
        real(dp) :: y1(size(y0))
        real(dp) :: f1(size(y0))
        real(dp) :: d0
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: h0
        real(dp) :: h1
        real(dp) :: interval

        interval = abs(t_bound - t0)
        scale = max(atol + rtol * abs(y0), tiny(1.0_dp))
        d0 = rms_norm(y0 / scale)
        d1 = rms_norm(f0 / scale)
        if (d0 < 1.0e-5_dp .or. d1 < 1.0e-5_dp) then
            h0 = 1.0e-6_dp
        else
            h0 = 0.01_dp * d0 / d1
        end if
        h0 = min(h0, interval, max_step)
        if (h0 <= 0.0_dp) then
            h_abs = min(interval, max_step)
            return
        end if

        y1 = y0 + direction * h0 * f0
        call fun(t0 + direction * h0, y1, f1)
        nfev = nfev + 1
        d2 = rms_norm((f1 - f0) / scale) / h0
        if (max(d1, d2) <= 1.0e-15_dp) then
            h1 = max(1.0e-6_dp, 1.0e-3_dp * h0)
        else
            h1 = (0.01_dp / max(d1, d2))**(1.0_dp / real(order + 1, dp))
        end if
        h_abs = min(100.0_dp * h0, h1, interval, max_step)
    end subroutine select_initial_step

    pure real(dp) function rms_norm(x) result(value)
        real(dp), intent(in) :: x(:) !! scaled vector whose RMS norm is requested
        value = sqrt(sum(x * x) / real(size(x), dp))
    end function rms_norm

    pure real(dp) function dop853_error_norm(err5, err3, scale, h) result(value)
        real(dp), intent(in) :: err5(:) !! fifth-order embedded error derivative
        real(dp), intent(in) :: err3(:) !! third-order embedded error derivative
        real(dp), intent(in) :: scale(:) !! component error scales
        real(dp), intent(in) :: h !! signed step size

        real(dp) :: err5_norm2
        real(dp) :: err3_norm2
        real(dp) :: denom

        err5_norm2 = sum((err5 / scale)**2)
        err3_norm2 = sum((err3 / scale)**2)
        if (err5_norm2 == 0.0_dp .and. err3_norm2 == 0.0_dp) then
            value = 0.0_dp
            return
        end if
        denom = err5_norm2 + 0.01_dp * err3_norm2
        value = abs(h) * err5_norm2 / sqrt(denom * real(size(scale), dp))
    end function dop853_error_norm

    pure integer function parse_method(method) result(method_id)
        character(len=*), intent(in), optional :: method !! requested solver name
        character(len=:), allocatable :: name

        if (.not. present(method)) then
            method_id = 45
            return
        end if
        name = lowercase(trim(adjustl(method)))
        select case (name)
        case ('rk23')
            method_id = 23
        case ('rk45')
            method_id = 45
        case ('dop853')
            method_id = 853
        case ('radau')
            method_id = 50
        case ('bdf')
            method_id = 15
        case ('lsoda')
            method_id = 12
        case default
            method_id = 0
        end select
    end function parse_method

    pure function lowercase(text) result(out)
        character(len=*), intent(in) :: text !! input ASCII text
        character(len=len(text)) :: out
        integer :: i
        integer :: code

        out = text
        do i = 1, len(text)
            code = iachar(out(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

    pure logical function valid_t_eval(values, t0, t1, direction) result(valid)
        real(dp), intent(in) :: values(:) !! requested output points
        real(dp), intent(in) :: t0 !! interval start
        real(dp), intent(in) :: t1 !! interval end
        real(dp), intent(in) :: direction !! integration direction
        integer :: i

        valid = all(ieee_is_finite(values))
        if (.not. valid) return
        if (size(values) == 0) return
        if (direction > 0.0_dp) then
            valid = values(1) >= t0 .and. values(size(values)) <= t1
            do i = 2, size(values)
                valid = valid .and. values(i) >= values(i - 1)
            end do
        else if (direction < 0.0_dp) then
            valid = values(1) <= t0 .and. values(size(values)) >= t1
            do i = 2, size(values)
                valid = valid .and. values(i) <= values(i - 1)
            end do
        else
            valid = all(values == t0)
        end if
    end function valid_t_eval

    pure logical function valid_ordered_points(values, direction) result(valid)
        real(dp), intent(in) :: values(:) !! points to validate
        real(dp), intent(in) :: direction !! integration direction, +1, -1, or 0

        integer :: i

        valid = all(ieee_is_finite(values))
        if (.not. valid .or. size(values) <= 1 .or. direction == 0.0_dp) return
        do i = 2, size(values)
            if (direction * (values(i) - values(i - 1)) < 0.0_dp) then
                valid = .false.
                return
            end if
        end do
    end function valid_ordered_points


    pure logical function crossing_detected(g0, g1, direction) result(hit)
        real(dp), intent(in) :: g0 !! event value at step start
        real(dp), intent(in) :: g1 !! event value at step end
        real(dp), intent(in) :: direction !! requested crossing direction

        if (direction > 0.0_dp) then
            hit = g0 <= 0.0_dp .and. g1 >= 0.0_dp
        else if (direction < 0.0_dp) then
            hit = g0 >= 0.0_dp .and. g1 <= 0.0_dp
        else
            hit = (g0 <= 0.0_dp .and. g1 >= 0.0_dp) .or. &
                (g0 >= 0.0_dp .and. g1 <= 0.0_dp)
        end if
    end function crossing_detected

    subroutine locate_event(event, method_id, ta, ya, fa, tb, yb, fb, &
            dense_coeff, dense_h, dense_anchor, dense_order, dense_valid, ga, gb, &
            root_t, root_y, ok)
        procedure(ode_event) :: event !! event callback whose zero is sought
        integer, intent(in) :: method_id !! active solver identifier
        real(dp), intent(in) :: ta !! left step endpoint
        real(dp), intent(in) :: ya(:) !! state at left endpoint
        real(dp), intent(in) :: fa(:) !! derivative at left endpoint
        real(dp), intent(in) :: tb !! right step endpoint
        real(dp), intent(in) :: yb(:) !! state at right endpoint
        real(dp), intent(in) :: fb(:) !! derivative at right endpoint
        real(dp), intent(in) :: dense_coeff(:, :) !! method-native dense coefficients
        real(dp), intent(in) :: dense_h !! full step size or BDF difference spacing
        real(dp), intent(in) :: dense_anchor !! full accepted endpoint for BDF dense output
        integer, intent(in) :: dense_order !! BDF dense polynomial order; zero otherwise
        logical, intent(in) :: dense_valid !! whether dense coefficients are valid
        real(dp), intent(in) :: ga !! event value at left endpoint
        real(dp), intent(in) :: gb !! event value at right endpoint
        real(dp), intent(out) :: root_t !! located event time
        real(dp), intent(out) :: root_y(:) !! interpolated state at event time
        logical, intent(out) :: ok !! false if an event evaluation becomes nonfinite

        real(dp) :: left
        real(dp) :: right
        real(dp) :: gleft
        real(dp) :: gright
        real(dp) :: mid
        real(dp) :: gmid
        real(dp) :: ymid(size(ya))
        integer :: iter

        ok = .true.
        left = ta
        right = tb
        gleft = ga
        gright = gb
        if (gleft == 0.0_dp) then
            root_t = left
            root_y = ya
            return
        end if
        if (gright == 0.0_dp) then
            root_t = right
            root_y = yb
            return
        end if

        mid = 0.5_dp * (left + right)
        do iter = 1, 64
            mid = 0.5_dp * (left + right)
            call step_state(method_id, ta, ya, fa, tb, yb, fb, dense_coeff, dense_h, &
                dense_anchor, dense_order, dense_valid, mid, ymid)
            gmid = event(mid, ymid)
            if (.not. ieee_is_finite(gmid)) then
                ok = .false.
                root_t = mid
                root_y = ymid
                return
            end if
            if (gmid == 0.0_dp .or. abs(right - left) <= &
                    8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
            if ((gleft < 0.0_dp .and. gmid > 0.0_dp) .or. &
                    (gleft > 0.0_dp .and. gmid < 0.0_dp)) then
                right = mid
                gright = gmid
            else
                left = mid
                gleft = gmid
            end if
        end do
        root_t = mid
        root_y = ymid
    end subroutine locate_event

    subroutine step_state(method_id, t0, y0, f0, t1, y1, f1, dense_coeff, dense_h, &
            dense_anchor, dense_order, dense_valid, tq, yq)
        integer, intent(in) :: method_id !! active solver identifier
        real(dp), intent(in) :: t0 !! step start
        real(dp), intent(in) :: y0(:) !! state at step start
        real(dp), intent(in) :: f0(:) !! derivative at step start
        real(dp), intent(in) :: t1 !! stored step end
        real(dp), intent(in) :: y1(:) !! state at stored step end
        real(dp), intent(in) :: f1(:) !! derivative at stored step end
        real(dp), intent(in) :: dense_coeff(:, :) !! method-native dense coefficients
        real(dp), intent(in) :: dense_h !! full step size or BDF difference spacing
        real(dp), intent(in) :: dense_anchor !! full accepted endpoint for BDF dense output
        integer, intent(in) :: dense_order !! BDF dense polynomial order; zero otherwise
        logical, intent(in) :: dense_valid !! whether dense coefficients are valid
        real(dp), intent(in) :: tq !! interpolation point
        real(dp), intent(out) :: yq(:) !! interpolated state

        if (method_id == 45 .and. dense_valid) then
            call rk45_dense_state(t0, y0, dense_h, dense_coeff(:, 1:4), tq, yq)
        else if (method_id == 853 .and. dense_valid) then
            call dop853_dense_state(t0, y0, dense_h, dense_coeff, tq, yq)
        else if (method_id == 50 .and. dense_valid) then
            call radau_dense_state(t0, y0, dense_h, dense_coeff(:, 1:3), tq, yq)
        else if (method_id == 12 .and. dense_valid .and. dense_order < 0) then
            call lsoda_adams_dense_state(dense_anchor, dense_h, -dense_order, &
                dense_coeff, tq, yq)
        else if ((method_id == 15 .or. method_id == 12) .and. dense_valid .and. &
                dense_order > 0) then
            call bdf_dense_state(dense_anchor, dense_h, dense_order, dense_coeff, tq, yq)
        else
            call hermite_state(t0, y0, f0, t1, y1, f1, tq, yq)
        end if
    end subroutine step_state

    subroutine rk45_dense_state(t0, y0, h, coeff, tq, yq)
        real(dp), intent(in) :: t0 !! full step start
        real(dp), intent(in) :: y0(:) !! state at full step start
        real(dp), intent(in) :: h !! signed full RK45 step size
        real(dp), intent(in) :: coeff(:, :) !! quartic derivative coefficients, shape (n,4)
        real(dp), intent(in) :: tq !! interpolation point inside the accepted step
        real(dp), intent(out) :: yq(:) !! quartic-interpolated state

        real(dp) :: x

        if (h == 0.0_dp) then
            yq = y0
            return
        end if
        x = (tq - t0) / h
        yq = y0 + h * x * (coeff(:, 1) + x * (coeff(:, 2) + &
            x * (coeff(:, 3) + x * coeff(:, 4))))
    end subroutine rk45_dense_state

    subroutine hermite_state(t0, y0, f0, t1, y1, f1, tq, yq)
        real(dp), intent(in) :: t0 !! step start
        real(dp), intent(in) :: y0(:) !! state at step start
        real(dp), intent(in) :: f0(:) !! derivative at step start
        real(dp), intent(in) :: t1 !! step end
        real(dp), intent(in) :: y1(:) !! state at step end
        real(dp), intent(in) :: f1(:) !! derivative at step end
        real(dp), intent(in) :: tq !! interpolation point inside the step
        real(dp), intent(out) :: yq(:) !! cubic-Hermite interpolated state

        real(dp) :: h
        real(dp) :: s
        real(dp) :: h00
        real(dp) :: h10
        real(dp) :: h01
        real(dp) :: h11

        h = t1 - t0
        if (h == 0.0_dp) then
            yq = y0
            return
        end if
        s = (tq - t0) / h
        h00 = (1.0_dp + 2.0_dp * s) * (1.0_dp - s)**2
        h10 = s * (1.0_dp - s)**2
        h01 = s * s * (3.0_dp - 2.0_dp * s)
        h11 = s * s * (s - 1.0_dp)
        yq = h00 * y0 + h10 * h * f0 + h01 * y1 + h11 * h * f1
    end subroutine hermite_state

    subroutine ivp_dense_evaluate(self, t_eval, y_eval, status)
        class(ivp_dense_solution), intent(in) :: self
        real(dp), intent(in) :: t_eval(:) !! points at which to evaluate the stored solution
        real(dp), allocatable, intent(out) :: y_eval(:, :) !! interpolated states, shape (n,nt)
        integer, intent(out), optional :: status !! 0 on success, invalid-argument otherwise

        integer :: i
        integer :: j
        integer :: n
        integer :: nt
        logical :: increasing

        if (present(status)) status = ivp_status_success
        if (.not. allocated(self%t) .or. .not. allocated(self%y) .or. &
                .not. allocated(self%f)) then
            allocate(y_eval(0, 0))
            if (present(status)) status = ivp_status_invalid_argument
            return
        end if
        n = size(self%y, 1)
        nt = size(t_eval)
        allocate(y_eval(n, nt))
        if (size(self%t) == 1) then
            do i = 1, nt
                if (t_eval(i) /= self%t(1)) then
                    y_eval(:, i) = quiet_nan(0.0_dp)
                    if (present(status)) status = ivp_status_invalid_argument
                else
                    y_eval(:, i) = self%y(:, 1)
                end if
            end do
            return
        end if

        increasing = self%t(size(self%t)) > self%t(1)
        do i = 1, nt
            if (.not. point_in_interval(t_eval(i), self%t(1), self%t(size(self%t)))) then
                y_eval(:, i) = quiet_nan(0.0_dp)
                if (present(status)) status = ivp_status_invalid_argument
                cycle
            end if
            if (t_eval(i) == self%t(size(self%t))) then
                y_eval(:, i) = self%y(:, size(self%t))
                cycle
            end if
            j = locate_interval(self%t, t_eval(i), increasing)
            if (self%method_id == 45 .and. allocated(self%rk45_valid)) then
                if (self%rk45_valid(j)) then
                    call rk45_dense_state(self%t(j), self%y(:, j), self%rk45_h(j), &
                        self%rk45_q(:, :, j), t_eval(i), y_eval(:, i))
                else
                    call hermite_state(self%t(j), self%y(:, j), self%f(:, j), &
                        self%t(j + 1), self%y(:, j + 1), self%f(:, j + 1), &
                        t_eval(i), y_eval(:, i))
                end if
            else if (self%method_id == 853 .and. allocated(self%dop853_valid)) then
                if (self%dop853_valid(j)) then
                    call dop853_dense_state(self%t(j), self%y(:, j), self%dop853_h(j), &
                        self%dop853_f(:, :, j), t_eval(i), y_eval(:, i))
                else
                    call hermite_state(self%t(j), self%y(:, j), self%f(:, j), &
                        self%t(j + 1), self%y(:, j + 1), self%f(:, j + 1), &
                        t_eval(i), y_eval(:, i))
                end if
            else if (self%method_id == 50 .and. allocated(self%radau_valid)) then
                if (self%radau_valid(j)) then
                    call radau_dense_state(self%t(j), self%y(:, j), self%radau_h(j), &
                        self%radau_q(:, :, j), t_eval(i), y_eval(:, i))
                else
                    call hermite_state(self%t(j), self%y(:, j), self%f(:, j), &
                        self%t(j + 1), self%y(:, j + 1), self%f(:, j + 1), &
                        t_eval(i), y_eval(:, i))
                end if
            else if ((self%method_id == 15 .or. self%method_id == 12) .and. &
                    allocated(self%bdf_valid)) then
                if (self%bdf_valid(j) .and. self%method_id == 12 .and. &
                        self%bdf_order(j) < 0) then
                    call lsoda_adams_dense_state(self%bdf_anchor(j), self%bdf_h(j), &
                        -self%bdf_order(j), self%bdf_d(:, :, j), t_eval(i), y_eval(:, i))
                else if (self%bdf_valid(j)) then
                    call bdf_dense_state(self%bdf_anchor(j), self%bdf_h(j), &
                        self%bdf_order(j), self%bdf_d(:, :, j), t_eval(i), y_eval(:, i))
                else
                    call hermite_state(self%t(j), self%y(:, j), self%f(:, j), &
                        self%t(j + 1), self%y(:, j + 1), self%f(:, j + 1), &
                        t_eval(i), y_eval(:, i))
                end if
            else
                call hermite_state(self%t(j), self%y(:, j), self%f(:, j), &
                    self%t(j + 1), self%y(:, j + 1), self%f(:, j + 1), &
                    t_eval(i), y_eval(:, i))
            end if
        end do
    end subroutine ivp_dense_evaluate

    pure logical function point_in_interval(t, a, b) result(inside)
        real(dp), intent(in) :: t !! point to test
        real(dp), intent(in) :: a !! first endpoint
        real(dp), intent(in) :: b !! second endpoint
        inside = t >= min(a, b) .and. t <= max(a, b)
    end function point_in_interval

    pure integer function locate_interval(t, value, increasing) result(index)
        real(dp), intent(in) :: t(:) !! monotone solution nodes
        real(dp), intent(in) :: value !! interpolation point
        logical, intent(in) :: increasing !! true for increasing node order
        integer :: lo
        integer :: hi
        integer :: mid

        lo = 1
        hi = size(t) - 1
        do while (lo < hi)
            mid = (lo + hi + 1) / 2
            if (increasing) then
                if (t(mid) <= value) then
                    lo = mid
                else
                    hi = mid - 1
                end if
            else
                if (t(mid) >= value) then
                    lo = mid
                else
                    hi = mid - 1
                end if
            end if
        end do
        index = lo
    end function locate_interval

    subroutine append_state(t_store, y_store, f_store, dense_store, dense_h_store, &
            dense_anchor_store, dense_order_store, dense_valid_store, capacity, used, &
            t, y, f, dense_coeff, dense_h, dense_anchor, dense_order, dense_valid, &
            store_dense)
        real(dp), allocatable, intent(inout) :: t_store(:) !! stored solution times
        real(dp), allocatable, intent(inout) :: y_store(:, :) !! stored solution states
        real(dp), allocatable, intent(inout) :: f_store(:, :) !! stored derivatives
        real(dp), allocatable, intent(inout) :: dense_store(:, :, :) !! method dense coefficients
        real(dp), allocatable, intent(inout) :: dense_h_store(:) !! dense-output step spacings
        real(dp), allocatable, intent(inout) :: dense_anchor_store(:) !! BDF dense anchors
        integer, allocatable, intent(inout) :: dense_order_store(:) !! BDF dense orders
        logical, allocatable, intent(inout) :: dense_valid_store(:) !! native-dense validity flags
        integer, intent(inout) :: capacity !! current allocated point capacity
        integer, intent(inout) :: used !! number of occupied point slots
        real(dp), intent(in) :: t !! time to append
        real(dp), intent(in) :: y(:) !! state to append
        real(dp), intent(in) :: f(:) !! derivative to append
        real(dp), intent(in) :: dense_coeff(:, :) !! dense coefficients for preceding interval
        real(dp), intent(in) :: dense_h !! dense step spacing for preceding interval
        real(dp), intent(in) :: dense_anchor !! BDF full-step endpoint anchor
        integer, intent(in) :: dense_order !! BDF polynomial order; zero otherwise
        logical, intent(in) :: dense_valid !! method-native coefficients are valid
        logical, intent(in) :: store_dense !! whether native coefficients must be retained

        real(dp), allocatable :: t_new(:)
        real(dp), allocatable :: y_new(:, :)
        real(dp), allocatable :: f_new(:, :)
        real(dp), allocatable :: dense_new(:, :, :)
        real(dp), allocatable :: dense_h_new(:)
        real(dp), allocatable :: dense_anchor_new(:)
        integer, allocatable :: dense_order_new(:)
        logical, allocatable :: dense_valid_new(:)
        integer :: new_capacity
        integer :: interval_index

        if (used == capacity) then
            new_capacity = 2 * capacity
            allocate(t_new(new_capacity), y_new(size(y, 1), new_capacity), &
                f_new(size(f, 1), new_capacity))
            t_new(1:used) = t_store(1:used)
            y_new(:, 1:used) = y_store(:, 1:used)
            f_new(:, 1:used) = f_store(:, 1:used)
            call move_alloc(t_new, t_store)
            call move_alloc(y_new, y_store)
            call move_alloc(f_new, f_store)
            if (store_dense) then
                allocate(dense_new(size(dense_coeff, 1), size(dense_coeff, 2), new_capacity), &
                    dense_h_new(new_capacity), dense_anchor_new(new_capacity), &
                    dense_order_new(new_capacity), dense_valid_new(new_capacity))
                dense_new = 0.0_dp
                dense_h_new = 0.0_dp
                dense_anchor_new = 0.0_dp
                dense_order_new = 0
                dense_valid_new = .false.
                if (used > 1) then
                    dense_new(:, :, 1:used - 1) = dense_store(:, :, 1:used - 1)
                    dense_h_new(1:used - 1) = dense_h_store(1:used - 1)
                    dense_anchor_new(1:used - 1) = dense_anchor_store(1:used - 1)
                    dense_order_new(1:used - 1) = dense_order_store(1:used - 1)
                    dense_valid_new(1:used - 1) = dense_valid_store(1:used - 1)
                end if
                call move_alloc(dense_new, dense_store)
                call move_alloc(dense_h_new, dense_h_store)
                call move_alloc(dense_anchor_new, dense_anchor_store)
                call move_alloc(dense_order_new, dense_order_store)
                call move_alloc(dense_valid_new, dense_valid_store)
            end if
            capacity = new_capacity
        end if
        interval_index = used
        if (store_dense) then
            dense_valid_store(interval_index) = dense_valid
            dense_h_store(interval_index) = dense_h
            dense_anchor_store(interval_index) = dense_anchor
            dense_order_store(interval_index) = dense_order
            if (dense_valid) dense_store(:, :, interval_index) = dense_coeff
        end if
        used = used + 1
        t_store(used) = t
        y_store(:, used) = y
        f_store(:, used) = f
    end subroutine append_state

    subroutine append_event(times, states, indices, capacity, used, t, y, event_index)
        real(dp), allocatable, intent(inout) :: times(:) !! stored event times
        real(dp), allocatable, intent(inout) :: states(:, :) !! stored event states
        integer, allocatable, intent(inout) :: indices(:) !! event-function indices
        integer, intent(inout) :: capacity !! event storage capacity
        integer, intent(inout) :: used !! number of stored events
        real(dp), intent(in) :: t !! event time to append
        real(dp), intent(in) :: y(:) !! event state to append
        integer, intent(in) :: event_index !! one-based event-function index

        real(dp), allocatable :: times_new(:)
        real(dp), allocatable :: states_new(:, :)
        integer, allocatable :: indices_new(:)
        integer :: new_capacity

        if (used == capacity) then
            new_capacity = 2 * capacity
            allocate(times_new(new_capacity), states_new(size(y), new_capacity), &
                indices_new(new_capacity))
            times_new(1:used) = times(1:used)
            states_new(:, 1:used) = states(:, 1:used)
            indices_new(1:used) = indices(1:used)
            call move_alloc(times_new, times)
            call move_alloc(states_new, states)
            call move_alloc(indices_new, indices)
            capacity = new_capacity
        end if
        used = used + 1
        times(used) = t
        states(:, used) = y
        indices(used) = event_index
    end subroutine append_event

    subroutine sort_event_roots(times, states, indices, nroots, direction)
        real(dp), intent(inout) :: times(:) !! root times to sort
        real(dp), intent(inout) :: states(:, :) !! root states paired with times
        integer, intent(inout) :: indices(:) !! event indices paired with times
        integer, intent(in) :: nroots !! number of valid roots in arrays
        real(dp), intent(in) :: direction !! integration direction, +1 or -1

        real(dp) :: tmp_t
        real(dp) :: tmp_y(size(states, 1))
        integer :: tmp_index
        integer :: i
        integer :: j
        logical :: swap_needed

        do i = 2, nroots
            j = i
            do while (j > 1)
                if (direction >= 0.0_dp) then
                    swap_needed = times(j) < times(j - 1)
                else
                    swap_needed = times(j) > times(j - 1)
                end if
                if (.not. swap_needed) exit
                tmp_t = times(j - 1)
                times(j - 1) = times(j)
                times(j) = tmp_t
                tmp_y = states(:, j - 1)
                states(:, j - 1) = states(:, j)
                states(:, j) = tmp_y
                tmp_index = indices(j - 1)
                indices(j - 1) = indices(j)
                indices(j) = tmp_index
                j = j - 1
            end do
        end do
    end subroutine sort_event_roots

    subroutine finish_result(res, t_store, y_store, f_store, dense_store, &
            dense_h_store, dense_anchor_store, dense_order_store, dense_valid_store, &
            method_id, used, t_eval, keep_dense, &
            event_times, event_states, event_indices_store, nevents, &
            n_event_functions, direction)
        type(solve_ivp_result), intent(inout) :: res
        real(dp), intent(in) :: t_store(:) !! accepted solution times
        real(dp), intent(in) :: y_store(:, :) !! accepted solution states
        real(dp), intent(in) :: f_store(:, :) !! accepted derivatives
        real(dp), intent(in) :: dense_store(:, :, :) !! stored DOP853 dense coefficients
        real(dp), intent(in) :: dense_h_store(:) !! stored dense-output step spacings
        real(dp), intent(in) :: dense_anchor_store(:) !! stored BDF full-step endpoints
        integer, intent(in) :: dense_order_store(:) !! stored BDF dense polynomial orders
        logical, intent(in) :: dense_valid_store(:) !! dense-coefficient validity flags
        integer, intent(in) :: method_id !! active solver identifier
        integer, intent(in) :: used !! number of accepted stored points
        real(dp), intent(in), optional :: t_eval(:) !! requested output points
        logical, intent(in) :: keep_dense !! retain dense solution object
        real(dp), intent(in) :: event_times(:) !! located event times
        real(dp), intent(in) :: event_states(:, :) !! located event states
        integer, intent(in) :: event_indices_store(:) !! event-function indices
        integer, intent(in) :: nevents !! number of located events
        integer, intent(in) :: n_event_functions !! number of event callbacks in the call
        real(dp), intent(in) :: direction !! integration direction

        type(ivp_dense_solution) :: local_sol
        real(dp), allocatable :: selected_t(:)
        real(dp), allocatable :: selected_y(:, :)
        integer :: count
        integer :: i
        integer :: interp_status

        allocate(local_sol%t(used), local_sol%y(size(y_store, 1), used), &
            local_sol%f(size(f_store, 1), used))
        local_sol%t = t_store(1:used)
        local_sol%y = y_store(:, 1:used)
        local_sol%f = f_store(:, 1:used)
        local_sol%method_id = method_id
        if (method_id == 45 .and. used > 1) then
            if (size(dense_valid_store) >= used - 1) then
                allocate(local_sol%rk45_q(size(y_store, 1), 4, used - 1), &
                    local_sol%rk45_h(used - 1), local_sol%rk45_valid(used - 1))
                local_sol%rk45_q = dense_store(:, 1:4, 1:used - 1)
                local_sol%rk45_h = dense_h_store(1:used - 1)
                local_sol%rk45_valid = dense_valid_store(1:used - 1)
            end if
        else if (method_id == 853 .and. used > 1) then
            if (size(dense_valid_store) >= used - 1) then
                allocate(local_sol%dop853_f(size(y_store, 1), 7, used - 1), &
                    local_sol%dop853_h(used - 1), local_sol%dop853_valid(used - 1))
                local_sol%dop853_f = dense_store(:, :, 1:used - 1)
                local_sol%dop853_h = dense_h_store(1:used - 1)
                local_sol%dop853_valid = dense_valid_store(1:used - 1)
            end if
        else if (method_id == 50 .and. used > 1) then
            if (size(dense_valid_store) >= used - 1) then
                allocate(local_sol%radau_q(size(y_store, 1), 3, used - 1), &
                    local_sol%radau_h(used - 1), local_sol%radau_valid(used - 1))
                local_sol%radau_q = dense_store(:, 1:3, 1:used - 1)
                local_sol%radau_h = dense_h_store(1:used - 1)
                local_sol%radau_valid = dense_valid_store(1:used - 1)
            end if
        else if ((method_id == 15 .or. method_id == 12) .and. used > 1) then
            if (size(dense_valid_store) >= used - 1) then
                allocate(local_sol%bdf_d(size(y_store, 1), 7, used - 1), &
                    local_sol%bdf_h(used - 1), local_sol%bdf_anchor(used - 1), &
                    local_sol%bdf_order(used - 1), local_sol%bdf_valid(used - 1))
                local_sol%bdf_d = dense_store(:, :, 1:used - 1)
                local_sol%bdf_h = dense_h_store(1:used - 1)
                local_sol%bdf_anchor = dense_anchor_store(1:used - 1)
                local_sol%bdf_order = dense_order_store(1:used - 1)
                local_sol%bdf_valid = dense_valid_store(1:used - 1)
            end if
        end if

        if (present(t_eval)) then
            count = 0
            do i = 1, size(t_eval)
                if (direction >= 0.0_dp) then
                    if (t_eval(i) <= local_sol%t(used)) count = count + 1
                else
                    if (t_eval(i) >= local_sol%t(used)) count = count + 1
                end if
            end do
            allocate(selected_t(count))
            if (count > 0) selected_t = t_eval(1:count)
            call local_sol%evaluate(selected_t, selected_y, interp_status)
            call move_alloc(selected_t, res%t)
            call move_alloc(selected_y, res%y)
        else
            allocate(res%t(used), res%y(size(y_store, 1), used))
            res%t = local_sol%t
            res%y = local_sol%y
        end if

        allocate(res%t_events(nevents), res%y_events(size(y_store, 1), nevents), &
            res%event_indices(nevents))
        allocate(res%event_counts(n_event_functions))
        res%n_event_functions = n_event_functions
        res%event_counts = 0
        if (nevents > 0) then
            res%t_events = event_times(1:nevents)
            res%y_events = event_states(:, 1:nevents)
            res%event_indices = event_indices_store(1:nevents)
            do i = 1, nevents
                res%event_counts(res%event_indices(i)) = &
                    res%event_counts(res%event_indices(i)) + 1
            end do
        end if

        if (keep_dense) then
            res%sol = local_sol
        end if
    end subroutine finish_result

    function solve_ivp_event_times_for(self, event_index) result(times)
        class(solve_ivp_result), intent(in) :: self !! completed integration result
        integer, intent(in) :: event_index !! one-based event-function index
        real(dp), allocatable :: times(:)

        integer :: count
        integer :: i
        integer :: j

        if (event_index < 1 .or. event_index > self%n_event_functions .or. &
                .not. allocated(self%event_indices)) then
            allocate(times(0))
            return
        end if

        count = 0
        if (allocated(self%event_counts)) count = self%event_counts(event_index)
        allocate(times(count))
        j = 0
        do i = 1, size(self%event_indices)
            if (self%event_indices(i) /= event_index) cycle
            j = j + 1
            times(j) = self%t_events(i)
        end do
    end function solve_ivp_event_times_for

    function solve_ivp_event_states_for(self, event_index) result(states)
        class(solve_ivp_result), intent(in) :: self !! completed integration result
        integer, intent(in) :: event_index !! one-based event-function index
        real(dp), allocatable :: states(:, :)

        integer :: count
        integer :: i
        integer :: j
        integer :: nstate

        nstate = 0
        if (allocated(self%y_events)) nstate = size(self%y_events, 1)
        if (event_index < 1 .or. event_index > self%n_event_functions .or. &
                .not. allocated(self%event_indices)) then
            allocate(states(nstate, 0))
            return
        end if

        count = 0
        if (allocated(self%event_counts)) count = self%event_counts(event_index)
        allocate(states(nstate, count))
        j = 0
        do i = 1, size(self%event_indices)
            if (self%event_indices(i) /= event_index) cycle
            j = j + 1
            states(:, j) = self%y_events(:, i)
        end do
    end function solve_ivp_event_states_for

    subroutine initialize_result(res)
        type(solve_ivp_result), intent(out) :: res !! result object to initialize
        res%status = ivp_status_failure
        res%success = .false.
        res%nfev = 0
        res%njev = 0
        res%nlu = 0
        res%nsteps = 0
        res%naccepted = 0
        res%nrejected = 0
        allocate(res%lsoda_hu(0), res%lsoda_tcur(0), res%lsoda_tolsf(0), &
            res%lsoda_tsw(0), res%lsoda_nst(0), res%lsoda_nfe(0), res%lsoda_nje(0), &
            res%lsoda_nqu(0), res%lsoda_mused(0))
    end subroutine initialize_result

    subroutine set_failure(res, status, message)
        type(solve_ivp_result), intent(inout) :: res !! result receiving failure state
        integer, intent(in) :: status !! failure status code
        character(len=*), intent(in) :: message !! human-readable failure description
        res%status = status
        res%success = .false.
        res%message = message
    end subroutine set_failure

end module scifort_integrate_ivp
