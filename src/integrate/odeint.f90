! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Legacy odeint-style facade over SciFort's LSODA solve_ivp implementation.
! The public calling conventions follow scipy.integrate.odeint where they map
! naturally to a type-safe Fortran API. See CODE_PROVENANCE.md.

module scifort_integrate_odeint
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : error_unit
    use scifort_kinds, only : dp
    use scifort_integrate_ivp, only : ivp_status_invalid_argument, solve_ivp, &
        solve_ivp_result
    implicit none
    private

    abstract interface
        subroutine odeint_rhs_yt(y, t, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t
        end subroutine odeint_rhs_yt

        subroutine odeint_rhs_ty(t, y, dydt)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t
        end subroutine odeint_rhs_ty

        subroutine odeint_jacobian_yt(y, t, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(out) :: jac(:, :) !! dense or packed Jacobian
        end subroutine odeint_jacobian_yt

        subroutine odeint_jacobian_ty(t, y, jac)
            import :: dp
            implicit none
            real(dp), intent(in) :: t !! independent variable
            real(dp), intent(in) :: y(:) !! state at t
            real(dp), intent(out) :: jac(:, :) !! dense or packed Jacobian
        end subroutine odeint_jacobian_ty
    end interface

    type, public :: odeint_result
        real(dp), allocatable :: y(:, :)
        integer :: status = ivp_status_invalid_argument
        logical :: success = .false.
        integer :: nfev = 0
        integer :: njev = 0
        integer :: nlu = 0
        integer :: nsteps = 0
        integer :: naccepted = 0
        integer :: nrejected = 0
        real(dp), allocatable :: hu(:)
        real(dp), allocatable :: tcur(:)
        real(dp), allocatable :: tolsf(:)
        real(dp), allocatable :: tsw(:)
        integer, allocatable :: nst(:)
        integer, allocatable :: nfe(:)
        integer, allocatable :: nje(:)
        integer, allocatable :: nqu(:)
        integer, allocatable :: mused(:)
        integer :: imxer = -1
        integer :: lenrw = 0
        integer :: leniw = 0
        character(len=:), allocatable :: message
    end type odeint_result

    public :: odeint
    public :: odeint_tfirst

contains

    function odeint(func, y0, t, rtol, atol, dfun, col_deriv, ml, mu, h0, hmax, &
            hmin, mxstep, tcrit, mxordn, mxords, full_output, ixpr, mxhnil, &
            printmessg) result(res)
        procedure(odeint_rhs_yt) :: func !! legacy right-hand side callback func(y,t)
        real(dp), intent(in) :: y0(:) !! initial state corresponding to t(1)
        real(dp), intent(in) :: t(:) !! monotone requested output times; repeats are allowed
        real(dp), intent(in), optional :: rtol(..) !! scalar or length-n relative tolerance; default sqrt(epsilon)
        real(dp), intent(in), optional :: atol(..) !! scalar or length-n absolute tolerance
        procedure(odeint_jacobian_yt), optional :: dfun !! optional dense or packed Jacobian
        logical, intent(in), optional :: col_deriv !! transpose Dfun output before use; default false
        integer, intent(in), optional :: ml !! lower Jacobian bandwidth; default dense
        integer, intent(in), optional :: mu !! upper Jacobian bandwidth; default dense
        real(dp), intent(in), optional :: h0 !! initial step magnitude; zero means automatic
        real(dp), intent(in), optional :: hmax !! maximum step magnitude; zero means unlimited
        real(dp), intent(in), optional :: hmin !! LSODA minimum retry step; default zero
        integer, intent(in), optional :: mxstep !! accepted-step cap per output call; 0 means 500
        real(dp), intent(in), optional :: tcrit(:) !! critical times that LSODA must not overshoot
        integer, intent(in), optional :: mxordn !! maximum nonstiff Adams order; 0 means 12
        integer, intent(in), optional :: mxords !! maximum stiff BDF order; 0 means 5
        logical, intent(in), optional :: full_output !! collect SciPy-style LSODA diagnostics
        logical, intent(in), optional :: ixpr !! print Adams/BDF method-switch messages; default false
        integer, intent(in), optional :: mxhnil !! maximum too-small-step messages; 0 means 10
        logical, intent(in), optional :: printmessg !! print the final integration message; default false
        type(odeint_result) :: res

        logical :: collect_output
        logical :: transpose_jac
        procedure(odeint_jacobian_yt), pointer :: dfun_cb

        collect_output = .false.
        if (present(full_output)) collect_output = full_output
        transpose_jac = .false.
        if (present(col_deriv)) transpose_jac = col_deriv

        if (present(dfun)) then
            dfun_cb => dfun
            res = odeint_core(rhs_adapter, y0, t, rtol, atol, jac_adapter, ml, mu, &
                h0, hmax, hmin, mxstep, tcrit, mxordn, mxords, collect_output, &
                ixpr, mxhnil, printmessg)
        else
            res = odeint_core(rhs_adapter, y0, t, rtol, atol, ml=ml, mu=mu, h0=h0, &
                hmax=hmax, hmin=hmin, mxstep=mxstep, tcrit=tcrit, mxordn=mxordn, &
                mxords=mxords, full_output=collect_output, ixpr=ixpr, mxhnil=mxhnil, &
                printmessg=printmessg)
        end if

    contains

        subroutine rhs_adapter(t_now, y_now, dydt)
            real(dp), intent(in) :: t_now !! independent variable
            real(dp), intent(in) :: y_now(:) !! state at t_now
            real(dp), intent(out) :: dydt(:) !! derivative dy/dt at t_now

            call func(y_now, t_now, dydt)
        end subroutine rhs_adapter

        subroutine jac_adapter(t_now, y_now, jac_out)
            real(dp), intent(in) :: t_now !! independent variable
            real(dp), intent(in) :: y_now(:) !! state at t_now
            real(dp), intent(out) :: jac_out(:, :) !! dense or packed Jacobian for LSODA

            real(dp), allocatable :: transposed(:, :)

            if (transpose_jac) then
                allocate(transposed(size(jac_out, 2), size(jac_out, 1)))
                call dfun_cb(y_now, t_now, transposed)
                jac_out = transpose(transposed)
            else
                call dfun_cb(y_now, t_now, jac_out)
            end if
        end subroutine jac_adapter

    end function odeint

    function odeint_tfirst(func, y0, t, rtol, atol, dfun, col_deriv, ml, mu, h0, &
            hmax, hmin, mxstep, tcrit, mxordn, mxords, full_output, ixpr, mxhnil, &
            printmessg) result(res)
        procedure(odeint_rhs_ty) :: func !! right-hand side callback func(t,y)
        real(dp), intent(in) :: y0(:) !! initial state corresponding to t(1)
        real(dp), intent(in) :: t(:) !! monotone requested output times; repeats are allowed
        real(dp), intent(in), optional :: rtol(..) !! scalar or length-n relative tolerance; default sqrt(epsilon)
        real(dp), intent(in), optional :: atol(..) !! scalar or length-n absolute tolerance
        procedure(odeint_jacobian_ty), optional :: dfun !! optional dense or packed Jacobian
        logical, intent(in), optional :: col_deriv !! transpose Dfun output before use; default false
        integer, intent(in), optional :: ml !! lower Jacobian bandwidth; default dense
        integer, intent(in), optional :: mu !! upper Jacobian bandwidth; default dense
        real(dp), intent(in), optional :: h0 !! initial step magnitude; zero means automatic
        real(dp), intent(in), optional :: hmax !! maximum step magnitude; zero means unlimited
        real(dp), intent(in), optional :: hmin !! LSODA minimum retry step; default zero
        integer, intent(in), optional :: mxstep !! accepted-step cap per output call; 0 means 500
        real(dp), intent(in), optional :: tcrit(:) !! critical times that LSODA must not overshoot
        integer, intent(in), optional :: mxordn !! maximum nonstiff Adams order; 0 means 12
        integer, intent(in), optional :: mxords !! maximum stiff BDF order; 0 means 5
        logical, intent(in), optional :: full_output !! collect SciPy-style LSODA diagnostics
        logical, intent(in), optional :: ixpr !! print Adams/BDF method-switch messages; default false
        integer, intent(in), optional :: mxhnil !! maximum too-small-step messages; 0 means 10
        logical, intent(in), optional :: printmessg !! print the final integration message; default false
        type(odeint_result) :: res

        logical :: collect_output
        logical :: transpose_jac
        procedure(odeint_jacobian_ty), pointer :: dfun_cb

        collect_output = .false.
        if (present(full_output)) collect_output = full_output
        transpose_jac = .false.
        if (present(col_deriv)) transpose_jac = col_deriv

        if (present(dfun)) then
            dfun_cb => dfun
            res = odeint_core(func, y0, t, rtol, atol, jac_adapter, ml, mu, h0, hmax, &
                hmin, mxstep, tcrit, mxordn, mxords, collect_output, ixpr, mxhnil, &
                printmessg)
        else
            res = odeint_core(func, y0, t, rtol, atol, ml=ml, mu=mu, h0=h0, hmax=hmax, &
                hmin=hmin, mxstep=mxstep, tcrit=tcrit, mxordn=mxordn, mxords=mxords, &
                full_output=collect_output, ixpr=ixpr, mxhnil=mxhnil, &
                printmessg=printmessg)
        end if

    contains

        subroutine jac_adapter(t_now, y_now, jac_out)
            real(dp), intent(in) :: t_now !! independent variable
            real(dp), intent(in) :: y_now(:) !! state at t_now
            real(dp), intent(out) :: jac_out(:, :) !! dense or packed Jacobian for LSODA

            real(dp), allocatable :: transposed(:, :)

            if (transpose_jac) then
                allocate(transposed(size(jac_out, 2), size(jac_out, 1)))
                call dfun_cb(t_now, y_now, transposed)
                jac_out = transpose(transposed)
            else
                call dfun_cb(t_now, y_now, jac_out)
            end if
        end subroutine jac_adapter

    end function odeint_tfirst

    function odeint_core(func, y0, t, rtol, atol, dfun, ml, mu, h0, hmax, hmin, &
            mxstep, tcrit, mxordn, mxords, full_output, ixpr, mxhnil, printmessg) result(res)
        procedure(odeint_rhs_ty) :: func !! t-first right-hand side callback
        real(dp), intent(in) :: y0(:) !! initial state corresponding to t(1)
        real(dp), intent(in) :: t(:) !! requested output times
        real(dp), intent(in), optional :: rtol(..) !! scalar or length-n relative tolerance
        real(dp), intent(in), optional :: atol(..) !! scalar or length-n absolute tolerance
        procedure(odeint_jacobian_ty), optional :: dfun !! optional dense or packed Jacobian
        integer, intent(in), optional :: ml !! lower Jacobian bandwidth
        integer, intent(in), optional :: mu !! upper Jacobian bandwidth
        real(dp), intent(in), optional :: h0 !! initial step magnitude; zero means automatic
        real(dp), intent(in), optional :: hmax !! maximum step magnitude; zero means unlimited
        real(dp), intent(in), optional :: hmin !! LSODA minimum retry step
        integer, intent(in), optional :: mxstep !! accepted-step cap per output call; 0 means 500
        real(dp), intent(in), optional :: tcrit(:) !! critical times passed to LSODA
        integer, intent(in), optional :: mxordn !! maximum nonstiff Adams order; 0 means default
        integer, intent(in), optional :: mxords !! maximum stiff BDF order; 0 means default
        logical, intent(in), optional :: full_output !! collect per-output LSODA diagnostics
        logical, intent(in), optional :: ixpr !! print Adams/BDF method-switch messages
        integer, intent(in), optional :: mxhnil !! maximum too-small-step messages; 0 means 10
        logical, intent(in), optional :: printmessg !! print the final integration message
        type(odeint_result) :: res

        type(solve_ivp_result) :: ivp
        real(dp), allocatable :: atol_use(:)
        real(dp), allocatable :: rtol_use(:)
        real(dp) :: h0_use
        real(dp) :: hmax_use
        real(dp) :: hmin_use
        integer :: max_steps_use
        integer :: mxstep_use
        integer :: lband_use
        integer :: uband_use
        integer :: mxordn_use
        integer :: mxords_use
        integer :: mxhnil_use
        integer :: n
        logical :: banded
        logical :: full_output_use
        logical :: ixpr_use
        logical :: printmessg_use
        logical :: use_first_step

        call initialize_odeint_result(res)
        full_output_use = .false.
        if (present(full_output)) full_output_use = full_output
        ixpr_use = .false.
        if (present(ixpr)) ixpr_use = ixpr
        printmessg_use = .false.
        if (present(printmessg)) printmessg_use = printmessg
        mxhnil_use = 10
        if (present(mxhnil)) then
            if (mxhnil < 0) then
                call set_odeint_failure(res, 'mxhnil must be nonnegative')
                return
            end if
            if (mxhnil > 0) mxhnil_use = mxhnil
        end if

        n = size(y0)
        if (n <= 0) then
            call set_odeint_failure(res, 'y0 must be nonempty')
            return
        end if
        if (size(t) <= 0) then
            call set_odeint_failure(res, 't must be nonempty')
            return
        end if
        if (.not. all(ieee_is_finite(t))) then
            call set_odeint_failure(res, 't must contain finite values')
            return
        end if

        allocate(rtol_use(n))
        rtol_use = sqrt(epsilon(1.0_dp))
        if (present(rtol)) then
            select rank (rtol)
            rank (0)
                rtol_use = rtol
            rank (1)
                if (size(rtol) /= n) then
                    call set_odeint_failure(res, 'vector rtol must have size(y0) elements')
                    return
                end if
                rtol_use = rtol
            rank default
                call set_odeint_failure(res, 'rtol must be a scalar or rank-one array')
                return
            end select
        end if

        allocate(atol_use(n))
        atol_use = sqrt(epsilon(1.0_dp))
        if (present(atol)) then
            select rank (atol)
            rank (0)
                atol_use = atol
            rank (1)
                if (size(atol) /= n) then
                    call set_odeint_failure(res, 'vector atol must have size(y0) elements')
                    return
                end if
                atol_use = atol
            rank default
                call set_odeint_failure(res, 'atol must be a scalar or rank-one array')
                return
            end select
        end if

        h0_use = 0.0_dp
        if (present(h0)) h0_use = h0
        if (.not. ieee_is_finite(h0_use) .or. h0_use < 0.0_dp) then
            call set_odeint_failure(res, 'h0 must be nonnegative and finite')
            return
        end if
        use_first_step = h0_use > 0.0_dp

        hmax_use = huge(1.0_dp)
        if (present(hmax)) then
            if (.not. ieee_is_finite(hmax) .or. hmax < 0.0_dp) then
                call set_odeint_failure(res, 'hmax must be nonnegative and finite')
                return
            end if
            if (hmax > 0.0_dp) hmax_use = hmax
        end if

        hmin_use = 0.0_dp
        if (present(hmin)) hmin_use = hmin
        if (.not. ieee_is_finite(hmin_use) .or. hmin_use < 0.0_dp) then
            call set_odeint_failure(res, 'hmin must be nonnegative and finite')
            return
        end if

        mxstep_use = 500
        if (present(mxstep)) then
            if (mxstep < 0) then
                call set_odeint_failure(res, 'mxstep must be nonnegative')
                return
            end if
            if (mxstep > 0) mxstep_use = mxstep
        end if
        max_steps_use = huge(max_steps_use)

        mxordn_use = 12
        if (present(mxordn)) then
            if (mxordn < 0) then
                call set_odeint_failure(res, 'mxordn must be nonnegative')
                return
            end if
            if (mxordn > 0) mxordn_use = min(mxordn, 12)
        end if
        mxords_use = 5
        if (present(mxords)) then
            if (mxords < 0) then
                call set_odeint_failure(res, 'mxords must be nonnegative')
                return
            end if
            if (mxords > 0) mxords_use = min(mxords, 5)
        end if

        banded = present(ml) .or. present(mu)
        lband_use = 0
        uband_use = 0
        if (present(ml)) lband_use = ml
        if (present(mu)) uband_use = mu
        if (banded) then
            if (lband_use < 0 .or. lband_use >= n .or. uband_use < 0 .or. uband_use >= n) then
                call set_odeint_failure(res, 'ml and mu must be between 0 and size(y0)-1')
                return
            end if
        end if

        if (banded) then
            if (use_first_step) then
                ivp = solve_ivp(func, [t(1), t(size(t))], y0, method='LSODA', t_eval=t, &
                    rtol=rtol_use, atol=atol_use, first_step=h0_use, max_step=hmax_use, &
                    max_steps=max_steps_use, jac=dfun, min_step=hmin_use, &
                    lband=lband_use, uband=uband_use, lsoda_tcrit=tcrit, &
                    lsoda_max_order_adams=mxordn_use, lsoda_max_order_bdf=mxords_use, &
                    lsoda_output_max_steps=mxstep_use, &
                    lsoda_collect_output_info=full_output_use, &
                    lsoda_print_switches=ixpr_use, lsoda_max_hnil_messages=mxhnil_use)
            else
                ivp = solve_ivp(func, [t(1), t(size(t))], y0, method='LSODA', t_eval=t, &
                    rtol=rtol_use, atol=atol_use, max_step=hmax_use, max_steps=max_steps_use, &
                    jac=dfun, min_step=hmin_use, lband=lband_use, uband=uband_use, &
                    lsoda_tcrit=tcrit, lsoda_max_order_adams=mxordn_use, &
                    lsoda_max_order_bdf=mxords_use, lsoda_output_max_steps=mxstep_use, &
                    lsoda_collect_output_info=full_output_use, &
                    lsoda_print_switches=ixpr_use, lsoda_max_hnil_messages=mxhnil_use)
            end if
        else
            if (use_first_step) then
                ivp = solve_ivp(func, [t(1), t(size(t))], y0, method='LSODA', t_eval=t, &
                    rtol=rtol_use, atol=atol_use, first_step=h0_use, max_step=hmax_use, &
                    max_steps=max_steps_use, jac=dfun, min_step=hmin_use, &
                    lsoda_tcrit=tcrit, lsoda_max_order_adams=mxordn_use, &
                    lsoda_max_order_bdf=mxords_use, lsoda_output_max_steps=mxstep_use, &
                    lsoda_collect_output_info=full_output_use, &
                    lsoda_print_switches=ixpr_use, lsoda_max_hnil_messages=mxhnil_use)
            else
                ivp = solve_ivp(func, [t(1), t(size(t))], y0, method='LSODA', t_eval=t, &
                    rtol=rtol_use, atol=atol_use, max_step=hmax_use, max_steps=max_steps_use, &
                    jac=dfun, min_step=hmin_use, lsoda_tcrit=tcrit, &
                    lsoda_max_order_adams=mxordn_use, lsoda_max_order_bdf=mxords_use, &
                    lsoda_output_max_steps=mxstep_use, &
                    lsoda_collect_output_info=full_output_use, &
                    lsoda_print_switches=ixpr_use, lsoda_max_hnil_messages=mxhnil_use)
            end if
        end if

        call copy_ivp_result(ivp, res)
        if (full_output_use) then
            call set_odeint_work_lengths(res, n, banded, lband_use, uband_use, &
                mxordn_use, mxords_use)
        end if
        if (printmessg_use) call print_odeint_message(res, t)
    end function odeint_core

    subroutine print_odeint_message(res, t)
        type(odeint_result), intent(in) :: res !! completed odeint result to report
        real(dp), intent(in) :: t(:) !! requested output times associated with the result

        if (res%success) then
            if (maxval(t) <= minval(t)) then
                write(error_unit, '(a)') 'Nothing was done; the integration time was 0.'
            else
                write(error_unit, '(a)') 'Integration successful.'
            end if
        else if (allocated(res%message)) then
            write(error_unit, '(a)') trim(res%message)
        else
            write(error_unit, '(a)') 'Integration failed.'
        end if
    end subroutine print_odeint_message

    subroutine initialize_odeint_result(res)
        type(odeint_result), intent(out) :: res !! result object reset to an empty failure state

        allocate(res%y(0, 0))
        res%status = ivp_status_invalid_argument
        res%success = .false.
        res%nfev = 0
        res%njev = 0
        res%nlu = 0
        res%nsteps = 0
        res%naccepted = 0
        res%nrejected = 0
        allocate(res%hu(0), res%tcur(0), res%tolsf(0), res%tsw(0), res%nst(0), &
            res%nfe(0), res%nje(0), res%nqu(0), res%mused(0))
        res%imxer = -1
        res%lenrw = 0
        res%leniw = 0
        res%message = ''
    end subroutine initialize_odeint_result

    subroutine set_odeint_failure(res, message)
        type(odeint_result), intent(inout) :: res !! result object to mark invalid
        character(len=*), intent(in) :: message !! human-readable validation failure

        res%status = ivp_status_invalid_argument
        res%success = .false.
        res%message = message
    end subroutine set_odeint_failure

    subroutine copy_ivp_result(ivp, res)
        type(solve_ivp_result), intent(in) :: ivp !! underlying LSODA solve_ivp result
        type(odeint_result), intent(inout) :: res !! odeint-style result to populate

        if (allocated(res%y)) deallocate(res%y)
        if (allocated(ivp%y)) then
            allocate(res%y(size(ivp%y, 2), size(ivp%y, 1)))
            res%y = transpose(ivp%y)
        else
            allocate(res%y(0, 0))
        end if
        res%status = ivp%status
        res%success = ivp%success
        res%nfev = ivp%nfev
        res%njev = ivp%njev
        res%nlu = ivp%nlu
        res%nsteps = ivp%nsteps
        res%naccepted = ivp%naccepted
        res%nrejected = ivp%nrejected
        call copy_real_vector(ivp%lsoda_hu, res%hu)
        call copy_real_vector(ivp%lsoda_tcur, res%tcur)
        call copy_real_vector(ivp%lsoda_tolsf, res%tolsf)
        call copy_real_vector(ivp%lsoda_tsw, res%tsw)
        call copy_integer_vector(ivp%lsoda_nst, res%nst)
        call copy_integer_vector(ivp%lsoda_nfe, res%nfe)
        call copy_integer_vector(ivp%lsoda_nje, res%nje)
        call copy_integer_vector(ivp%lsoda_nqu, res%nqu)
        call copy_integer_vector(ivp%lsoda_mused, res%mused)
        res%imxer = -1
        if (allocated(ivp%message)) then
            res%message = ivp%message
        else
            res%message = ''
        end if
    end subroutine copy_ivp_result

    subroutine copy_real_vector(source, destination)
        real(dp), allocatable, intent(in) :: source(:) !! source diagnostic vector
        real(dp), allocatable, intent(inout) :: destination(:) !! copied diagnostic vector

        if (allocated(destination)) deallocate(destination)
        if (allocated(source)) then
            allocate(destination(size(source)))
            destination = source
        else
            allocate(destination(0))
        end if
    end subroutine copy_real_vector

    subroutine copy_integer_vector(source, destination)
        integer, allocatable, intent(in) :: source(:) !! source diagnostic vector
        integer, allocatable, intent(inout) :: destination(:) !! copied diagnostic vector

        if (allocated(destination)) deallocate(destination)
        if (allocated(source)) then
            allocate(destination(size(source)))
            destination = source
        else
            allocate(destination(0))
        end if
    end subroutine copy_integer_vector

    subroutine set_odeint_work_lengths(res, n, banded, ml, mu, mxordn, mxords)
        type(odeint_result), intent(inout) :: res !! result receiving ODEPACK work lengths
        integer, intent(in) :: n !! number of equations
        logical, intent(in) :: banded !! true for packed banded Jacobian storage
        integer, intent(in) :: ml !! lower Jacobian bandwidth
        integer, intent(in) :: mu !! upper Jacobian bandwidth
        integer, intent(in) :: mxordn !! maximum Adams order
        integer, intent(in) :: mxords !! maximum BDF order

        integer :: lmat
        integer :: lrn
        integer :: lrs

        lrn = 20 + n * (mxordn + 1) + 3 * n
        if (banded) then
            lmat = (2 * ml + mu + 1) * n + 2
        else
            lmat = n * n + 2
        end if
        lrs = 20 + n * (mxords + 1) + 3 * n + lmat
        res%lenrw = max(lrn, lrs)
        res%leniw = 20 + n
    end subroutine set_odeint_work_lengths

end module scifort_integrate_odeint
