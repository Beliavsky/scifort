! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Continuous-time LTI simulation and step/impulse responses.
! Input/output samples are rows, channels are columns; multiple independent
! input-specific response experiments occupy the third array dimension.
module scifort_signal_continuous_simulation
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_state_space, only : tf2ss, zpk2ss, state_space_result, state_space_success
    use scifort_signal_cont2discrete, only : matrix_exp_pade
    use scifort_signal_analog, only : polynomial_roots
    implicit none
    private

    integer, parameter, public :: continuous_sim_success = 0
    integer, parameter, public :: continuous_sim_invalid_input = 1
    integer, parameter, public :: continuous_sim_numerical_failure = 2

    type, public :: continuous_simulation_result
        real(dp), allocatable :: t(:) !! simulation sample times
        complex(dp), allocatable :: y(:,:) !! outputs (nt,nout)
        complex(dp), allocatable :: x(:,:) !! states (nt,nstates)
        integer :: status = continuous_sim_invalid_input !! status code
    end type continuous_simulation_result

    type, public :: continuous_response_result
        real(dp), allocatable :: t(:) !! requested response times
        complex(dp), allocatable :: y(:,:,:) !! response (nt,nout,nin)
        integer :: status = continuous_sim_invalid_input !! status code
    end type continuous_response_result

    interface lsim
        module procedure lsim_ss_real, lsim_ss_complex
        module procedure lsim_ss_real_complex_input, lsim_ss_complex_real_input
        module procedure lsim_tf_real, lsim_tf_complex
        module procedure lsim_tf_real_complex_input, lsim_tf_complex_real_input
    end interface lsim
    interface step
        module procedure step_ss_real, step_ss_complex
        module procedure step_tf_real, step_tf_complex
    end interface step
    interface impulse
        module procedure impulse_ss_real, impulse_ss_complex
        module procedure impulse_tf_real, impulse_tf_complex
    end interface impulse
    interface lsim_zpk
        module procedure lsim_zpk_real, lsim_zpk_complex
        module procedure lsim_zpk_real_complex_input, lsim_zpk_complex_real_input
        module procedure lsim_zpk_real_complex_gain, lsim_zpk_real_complex_gain_real_input
        module procedure lsim_zpk_complex_real_gain, lsim_zpk_complex_real_gain_real_input
    end interface lsim_zpk
    interface step_zpk
        module procedure step_zpk_real, step_zpk_complex
        module procedure step_zpk_real_complex_gain, step_zpk_complex_real_gain
    end interface step_zpk
    interface impulse_zpk
        module procedure impulse_zpk_real, impulse_zpk_complex
        module procedure impulse_zpk_real_complex_gain, impulse_zpk_complex_real_gain
    end interface impulse_zpk
    public :: lsim, step, impulse, lsim_zpk, step_zpk, impulse_zpk

contains

    pure logical function finite_matrix(v) result(ok)
        complex(dp), intent(in) :: v(:,:) !! matrix values
        ok = all(ieee_is_finite(real(v,dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_matrix

    pure logical function finite_vector(v) result(ok)
        complex(dp), intent(in) :: v(:) !! vector values
        ok = all(ieee_is_finite(real(v,dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_vector

    pure logical function valid_time(t,uniform) result(ok)
        real(dp), intent(in) :: t(:) !! nonnegative strictly increasing times
        logical, intent(in) :: uniform !! require uniform intervals
        real(dp) :: dt
        integer :: nt
        ok = .false.
        nt = size(t)
        if (nt < 1) return
        if (.not. all(ieee_is_finite(t))) return
        if (t(1) < 0.0_dp) return
        if (nt > 1) then
            if (any(t(2:) <= t(:nt-1))) return
            if (uniform .and. nt > 2) then
                dt = (t(nt)-t(1))/real(nt-1,dp)
                if (any(abs((t(2:)-t(:nt-1))-dt) > 1.0e-8_dp*dt)) return
            end if
        end if
        ok = .true.
    end function valid_time

    pure logical function valid_system(a,b,c,d) result(ok)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough
        integer :: n, ni, no
        n = size(a,1)
        ni = size(b,2)
        no = size(c,1)
        ok = size(a,2)==n .and. size(b,1)==n .and. size(c,2)==n .and. &
            size(d,1)==no .and. size(d,2)==ni .and. ni>0 .and. no>0
        if (.not. ok) return
        ok = finite_matrix(a) .and. finite_matrix(b) .and. finite_matrix(c) .and. finite_matrix(d)
    end function valid_system

    ! Compute x(t+dt) = Ad*x(t) + Bleft*u(t) + Bright*u(t+dt).
    ! The block exponential gives exact propagation for piecewise constant
    ! inputs or inputs interpolated linearly between consecutive samples.
    subroutine transition(a,b,dt,interp,ad,bleft,bright,info)
        complex(dp), intent(in) :: a(:,:) !! continuous state matrix
        complex(dp), intent(in) :: b(:,:) !! continuous input matrix
        real(dp), intent(in) :: dt !! interval length, nonnegative
        logical, intent(in) :: interp !! linear if true, zero-order hold otherwise
        complex(dp), allocatable, intent(out) :: ad(:,:) !! state transition
        complex(dp), allocatable, intent(out) :: bleft(:,:) !! left input weight
        complex(dp), allocatable, intent(out) :: bright(:,:) !! right input weight
        integer, intent(out) :: info !! exponential status
        complex(dp), allocatable :: augmented(:,:), ex(:,:)
        integer :: n, ni, j
        n = size(a,1)
        ni = size(b,2)
        allocate(bleft(n,ni),bright(n,ni))
        bleft = (0.0_dp,0.0_dp)
        bright = (0.0_dp,0.0_dp)
        if (interp) then
            allocate(augmented(n+2*ni,n+2*ni))
            augmented = (0.0_dp,0.0_dp)
            augmented(1:n,1:n) = dt*a
            augmented(1:n,n+1:n+ni) = dt*b
            do j = 1, ni
                augmented(n+j,n+ni+j) = (1.0_dp,0.0_dp)
            end do
            call matrix_exp_pade(augmented,ex,info)
            if (info /= 0) return
            ad = ex(1:n,1:n)
            bright = ex(1:n,n+ni+1:)
            bleft = ex(1:n,n+1:n+ni)-bright
        else
            allocate(augmented(n+ni,n+ni))
            augmented = (0.0_dp,0.0_dp)
            augmented(1:n,1:n) = dt*a
            augmented(1:n,n+1:) = dt*b
            call matrix_exp_pade(augmented,ex,info)
            if (info /= 0) return
            ad = ex(1:n,1:n)
            bleft = ex(1:n,n+1:)
        end if
    end subroutine transition

    ! Select SciPy-style default response times when t is absent. The slowest
    ! exponential mode sets a seven-time-constant range; purely imaginary
    ! and zero eigenvalues use a unit time constant. For n=1 return t=0.
    subroutine response_time_grid(a,t,n,grid,info)
        complex(dp), intent(in) :: a(:,:) !! continuous state matrix
        real(dp), intent(in), optional :: t(:) !! explicit response times
        integer, intent(in), optional :: n !! number of automatic samples (default 100)
        real(dp), allocatable, intent(out) :: grid(:) !! selected response times
        integer, intent(out) :: info !! continuous_sim_* status
        complex(dp), allocatable :: coeff(:), factors(:,:), prod(:,:), poles(:)
        complex(dp) :: ck
        real(dp) :: rate, endpoint, pole_real, pole_scale
        integer :: ns, nt, i, j, k
        logical :: ok, upper, lower

        info = continuous_sim_invalid_input
        if (present(t)) then
            if (.not. valid_time(t,.false.)) return
            grid = t
            info = continuous_sim_success
            return
        end if
        nt = 100
        if (present(n)) nt = n
        if (nt < 1) return
        ns = size(a,1)
        if (size(a,2) /= ns) return
        rate = 1.0_dp
        if (ns > 0) then
            allocate(poles(ns))
            upper = .true.
            lower = .true.
            do j = 1, ns
                do i = 1, ns
                    if (i > j .and. abs(a(i,j)) > 0.0_dp) upper = .false.
                    if (i < j .and. abs(a(i,j)) > 0.0_dp) lower = .false.
                end do
            end do
            if (upper .or. lower) then
                do i = 1, ns
                    poles(i) = a(i,i)
                end do
            else
                allocate(coeff(ns+1),factors(ns,ns),prod(ns,ns))
                coeff(1) = (1.0_dp,0.0_dp)
                factors = (0.0_dp,0.0_dp)
                do i = 1, ns
                    factors(i,i) = (1.0_dp,0.0_dp)
                end do
                do k = 1, ns
                    prod = matmul(a,factors)
                    ck = (0.0_dp,0.0_dp)
                    do i = 1, ns
                        ck = ck-prod(i,i)
                    end do
                    ck = ck/real(k,dp)
                    coeff(k+1) = ck
                    factors = prod
                    do i = 1, ns
                        factors(i,i) = factors(i,i)+ck
                    end do
                end do
                if (.not. finite_vector(coeff)) then
                    info = continuous_sim_numerical_failure
                    return
                end if
                deallocate(poles)
                call polynomial_roots(coeff,poles,ok)
                if (.not. ok) then
                    info = continuous_sim_numerical_failure
                    return
                end if
                if (.not. allocated(poles)) then
                    info = continuous_sim_numerical_failure
                    return
                end if
                if (size(poles) /= ns) then
                    info = continuous_sim_numerical_failure
                    return
                end if
            end if
            if (.not. finite_vector(poles)) then
                info = continuous_sim_numerical_failure
                return
            end if
            rate = huge(1.0_dp)
            do i = 1, ns
                pole_real = abs(real(poles(i),dp))
                pole_scale = max(1.0_dp,abs(poles(i)))
                ! Suppress spurious real parts introduced by finding roots
                ! of a polynomial with purely imaginary eigenvalues.
                if (pole_real < 1.0e-12_dp*pole_scale) pole_real = 0.0_dp
                rate = min(rate,pole_real)
            end do
            if (rate == 0.0_dp) rate = 1.0_dp
        end if
        endpoint = 7.0_dp/rate
        if (.not. ieee_is_finite(endpoint) .or. endpoint <= 0.0_dp) then
            info = continuous_sim_numerical_failure
            return
        end if
        allocate(grid(nt))
        if (nt == 1) then
            grid(1) = 0.0_dp
        else
            do i = 1, nt
                grid(i) = endpoint*real(i-1,dp)/real(nt-1,dp)
            end do
        end if
        info = continuous_sim_success
    end subroutine response_time_grid

    function response_auto_core(a,b,c,d,t,x0,is_impulse,n) result(out)
        complex(dp), intent(in) :: a(:,:) !! continuous state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough
        real(dp), intent(in), optional :: t(:) !! response time grid, automatic if absent
        complex(dp), intent(in), optional :: x0(:) !! initial state
        logical, intent(in) :: is_impulse !! true for an impulse, false for step
        integer, intent(in), optional :: n !! automatic time grid length
        type(continuous_response_result) :: out
        real(dp), allocatable :: selected(:)
        integer :: info
        if (.not. valid_system(a,b,c,d)) return
        call response_time_grid(a,t,n,selected,info)
        if (info /= continuous_sim_success) then
            out%status = info
            return
        end if
        out = response_core(a,b,c,d,selected,x0,is_impulse)
    end function response_auto_core

    function lsim_core(a,b,c,d,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix (n,n)
        complex(dp), intent(in) :: b(:,:) !! input matrix (n,nin)
        complex(dp), intent(in) :: c(:,:) !! output matrix (nout,n)
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough (nout,nin)
        complex(dp), intent(in) :: u(:,:) !! sample-major input (nt,nin)
        real(dp), intent(in) :: t(:) !! equally spaced times starting at or after zero
        complex(dp), intent(in), optional :: x0(:) !! state at t=0
        logical, intent(in), optional :: interp !! true: FOH (default); false: ZOH
        type(continuous_simulation_result) :: out
        complex(dp), allocatable :: ad(:,:), bl(:,:), br(:,:), initial(:)
        integer :: n, ni, no, nt, i, info
        logical :: linear
        if (.not. valid_system(a,b,c,d)) return
        if (.not. valid_time(t,.true.)) return
        n = size(a,1)
        ni = size(b,2)
        no = size(c,1)
        nt = size(t)
        if (size(u,1) /= nt .or. size(u,2) /= ni) return
        if (.not. finite_matrix(u)) return
        if (present(x0)) then
            if (size(x0) /= n) return
            if (.not. finite_vector(x0)) return
        end if
        allocate(out%t(nt),out%y(nt,no),out%x(nt,n),initial(n))
        out%t = t
        initial = (0.0_dp,0.0_dp)
        if (present(x0)) initial = x0
        if (t(1) > 0.0_dp) then
            call matrix_exp_pade(a*t(1),ad,info)
            if (info /= 0) then
                out%status = continuous_sim_numerical_failure
                return
            end if
            initial = matmul(ad,initial)
        end if
        out%x(1,:) = initial
        linear = .true.
        if (present(interp)) linear = interp
        if (nt > 1) then
            call transition(a,b,(t(nt)-t(1))/real(nt-1,dp),linear,ad,bl,br,info)
            if (info /= 0) then
                out%status = continuous_sim_numerical_failure
                return
            end if
        end if
        do i = 1, nt
            out%y(i,:) = matmul(c,out%x(i,:)) + matmul(d,u(i,:))
            if (i < nt) out%x(i+1,:) = matmul(ad,out%x(i,:)) + &
                matmul(bl,u(i,:)) + matmul(br,u(i+1,:))
        end do
        if (.not. finite_matrix(out%y) .or. .not. finite_matrix(out%x)) then
            out%status = continuous_sim_numerical_failure
            return
        end if
        out%status = continuous_sim_success
    end function lsim_core

    function response_core(a,b,c,d,t,x0,is_impulse) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: t(:) !! strictly increasing, nonnegative times
        complex(dp), intent(in), optional :: x0(:) !! initial state at time zero
        logical, intent(in) :: is_impulse !! true: Dirac impulse; false: unit step
        type(continuous_response_result) :: out
        complex(dp), allocatable :: states(:,:), ad(:,:), bl(:,:), br(:,:)
        integer :: n, ni, no, nt, i, j, info
        real(dp) :: prev_t, dt
        if (.not. valid_system(a,b,c,d)) return
        if (.not. valid_time(t,.false.)) return
        n = size(a,1)
        ni = size(b,2)
        no = size(c,1)
        nt = size(t)
        if (present(x0)) then
            if (size(x0) /= n) return
            if (.not. finite_vector(x0)) return
        end if
        allocate(states(n,ni),out%t(nt),out%y(nt,no,ni))
        states = (0.0_dp,0.0_dp)
        if (present(x0)) then
            do j = 1, ni
                states(:,j) = x0
            end do
        end if
        if (is_impulse) then
            states = states+b
        end if
        out%t = t
        prev_t = 0.0_dp
        do i = 1, nt
            dt = t(i)-prev_t
            if (dt > 0.0_dp) then
                if (is_impulse) then
                    call matrix_exp_pade(a*dt,ad,info)
                else
                    call transition(a,b,dt,.false.,ad,bl,br,info)
                end if
                if (info /= 0) then
                    out%status = continuous_sim_numerical_failure
                    return
                end if
                states = matmul(ad,states)
                if (.not. is_impulse) states = states+bl
            end if
            out%y(i,:,:) = matmul(c,states)
            if (.not. is_impulse) then
                do j = 1, ni
                    out%y(i,:,j) = out%y(i,:,j)+d(:,j)
                end do
            end if
            prev_t = t(i)
        end do
        if (.not. all(ieee_is_finite(real(out%y,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%y)))) then
            out%status = continuous_sim_numerical_failure
            return
        end if
        out%status = continuous_sim_success
    end function response_core

    function lsim_tf_core(num,den,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: num(:) !! numerator descending powers of s
        complex(dp), intent(in) :: den(:) !! denominator descending powers of s
        complex(dp), intent(in) :: u(:) !! scalar input
        real(dp), intent(in) :: t(:) !! sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolation policy
        type(continuous_simulation_result) :: out
        type(state_space_result) :: model
        model = tf2ss(num,den)
        if (model%status /= state_space_success) return
        out = lsim_core(model%a,model%b,model%c,model%d,reshape(u,[size(u),1]),t,x0,interp)
    end function lsim_tf_core

    function response_tf_core(num,den,t,x0,is_impulse,n) result(out)
        complex(dp), intent(in) :: num(:) !! numerator coefficients
        complex(dp), intent(in) :: den(:) !! denominator coefficients
        real(dp), intent(in), optional :: t(:) !! response times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in) :: is_impulse !! use impulse rather than step
        integer, intent(in), optional :: n !! automatic grid size
        type(continuous_response_result) :: out
        type(state_space_result) :: model
        model = tf2ss(num,den)
        if (model%status /= state_space_success) return
        out = response_auto_core(model%a,model%b,model%c,model%d,t,x0,is_impulse,n)
    end function response_tf_core

    function lsim_zpk_core(z,p,g,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        complex(dp), intent(in) :: u(:) !! scalar input
        real(dp), intent(in) :: t(:) !! sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! input interpolation policy
        type(continuous_simulation_result) :: out
        type(state_space_result) :: model
        model = zpk2ss(z,p,g)
        if (model%status /= state_space_success) return
        out = lsim_core(model%a,model%b,model%c,model%d,reshape(u,[size(u),1]),t,x0,interp)
    end function lsim_zpk_core

    function response_zpk_core(z,p,g,t,x0,is_impulse,n) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! response times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in) :: is_impulse !! select impulse
        integer, intent(in), optional :: n !! automatic grid size
        type(continuous_response_result) :: out
        type(state_space_result) :: model
        model = zpk2ss(z,p,g)
        if (model%status /= state_space_success) return
        out = response_auto_core(model%a,model%b,model%c,model%d,t,x0,is_impulse,n)
    end function response_zpk_core

    function lsim_ss_real(a,b,c,d,u,t,x0,interp) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: u(:,:) !! input samples, time first
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        real(dp), intent(in), optional :: x0(:) !! initial state
        logical, intent(in), optional :: interp !! linear interpolation by default
        type(continuous_simulation_result) :: out
        if (present(x0)) then
            out = lsim_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),cmplx(u,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),interp)
        else
            out = lsim_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),cmplx(u,0.0_dp,dp),t,interp=interp)
        end if
    end function lsim_ss_real

    function lsim_ss_complex(a,b,c,d,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        complex(dp), intent(in) :: u(:,:) !! input samples, time first
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        logical, intent(in), optional :: interp !! linear interpolation by default
        type(continuous_simulation_result) :: out
        out = lsim_core(a,b,c,d,u,t,x0,interp)
    end function lsim_ss_complex

    function lsim_ss_real_complex_input(a,b,c,d,u,t,x0,interp) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        complex(dp), intent(in) :: u(:,:) !! input samples, time first
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        logical, intent(in), optional :: interp !! linear interpolation by default
        type(continuous_simulation_result) :: out
        out = lsim_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),u,t,x0,interp)
    end function lsim_ss_real_complex_input

    function lsim_ss_complex_real_input(a,b,c,d,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: u(:,:) !! input samples, time first
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        logical, intent(in), optional :: interp !! linear interpolation by default
        type(continuous_simulation_result) :: out
        out = lsim_core(a,b,c,d,cmplx(u,0.0_dp,dp),t,x0,interp)
    end function lsim_ss_complex_real_input

    function lsim_tf_real(b,a,u,t,x0,interp) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients, descending powers
        real(dp), intent(in) :: a(:) !! denominator coefficients, descending powers
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        real(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolation flag
        type(continuous_simulation_result) :: out
        if (present(x0)) then
            out = lsim_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),cmplx(u,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),interp)
        else
            out = lsim_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),cmplx(u,0.0_dp,dp),t,interp=interp)
        end if
    end function lsim_tf_real

    function lsim_tf_complex(b,a,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, descending powers
        complex(dp), intent(in) :: a(:) !! denominator coefficients, descending powers
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolation flag
        type(continuous_simulation_result) :: out
        out = lsim_tf_core(b,a,u,t,x0,interp)
    end function lsim_tf_complex

    function lsim_tf_real_complex_input(b,a,u,t,x0,interp) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients, descending powers
        real(dp), intent(in) :: a(:) !! denominator coefficients, descending powers
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolation flag
        type(continuous_simulation_result) :: out
        out = lsim_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),u,t,x0,interp)
    end function lsim_tf_real_complex_input

    function lsim_tf_complex_real_input(b,a,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, descending powers
        complex(dp), intent(in) :: a(:) !! denominator coefficients, descending powers
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolation flag
        type(continuous_simulation_result) :: out
        out = lsim_tf_core(b,a,cmplx(u,0.0_dp,dp),t,x0,interp)
    end function lsim_tf_complex_real_input

    function step_ss_real(a,b,c,d,t,x0,n) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        real(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        if (present(x0)) then
            out = response_auto_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),.false.,n)
        else
            out = response_auto_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),t,is_impulse=.false.,n=n)
        end if
    end function step_ss_real

    function impulse_ss_real(a,b,c,d,t,x0,n) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        real(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        if (present(x0)) then
            out = response_auto_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),.true.,n)
        else
            out = response_auto_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),cmplx(c,0.0_dp,dp), &
                cmplx(d,0.0_dp,dp),t,is_impulse=.true.,n=n)
        end if
    end function impulse_ss_real

    function step_ss_complex(a,b,c,d,t,x0,n) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        out = response_auto_core(a,b,c,d,t,x0,.false.,n)
    end function step_ss_complex

    function impulse_ss_complex(a,b,c,d,t,x0,n) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        out = response_auto_core(a,b,c,d,t,x0,.true.,n)
    end function impulse_ss_complex

    function step_tf_real(b,a,t,x0,n) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        real(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        if (present(x0)) then
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),.false.,n)
        else
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),t,is_impulse=.false.,n=n)
        end if
    end function step_tf_real

    function impulse_tf_real(b,a,t,x0,n) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        real(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        if (present(x0)) then
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),.true.,n)
        else
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),t,is_impulse=.true.,n=n)
        end if
    end function impulse_tf_real

    function step_tf_complex(b,a,t,x0,n) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        out = response_tf_core(b,a,t,x0,.false.,n)
    end function step_tf_complex

    function impulse_tf_complex(b,a,t,x0,n) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        out = response_tf_core(b,a,t,x0,.true.,n)
    end function impulse_tf_complex

    function step_zpk_real(z,p,g,t,x0,n) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        real(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        if (present(x0)) then
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),.false.,n)
        else
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp),t,is_impulse=.false.,n=n)
        end if
    end function step_zpk_real

    function impulse_zpk_real(z,p,g,t,x0,n) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        real(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        if (present(x0)) then
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),.true.,n)
        else
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp),t,is_impulse=.true.,n=n)
        end if
    end function impulse_zpk_real

    function step_zpk_complex(z,p,g,t,x0,n) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        out = response_zpk_core(z,p,g,t,x0,.false.,n)
    end function step_zpk_complex

    function impulse_zpk_complex(z,p,g,t,x0,n) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! strictly increasing response times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        integer, intent(in), optional :: n !! automatic grid size (default 100)
        type(continuous_response_result) :: out
        out = response_zpk_core(z,p,g,t,x0,.true.,n)
    end function impulse_zpk_complex

    function lsim_zpk_real(z,p,g,u,t,x0,interp) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        real(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! input interpolation
        type(continuous_simulation_result) :: out
        if (present(x0)) then
            out = lsim_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp), &
                cmplx(u,0.0_dp,dp),t, &
                cmplx(x0,0.0_dp,dp),interp)
        else
            out = lsim_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp), &
                cmplx(u,0.0_dp,dp),t,interp=interp)
        end if
    end function lsim_zpk_real

    function lsim_zpk_complex(z,p,g,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! input interpolation
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(z,p,g,u,t,x0,interp)
    end function lsim_zpk_complex

    function step_zpk_real_complex_gain(z,p,g,t,x0,n) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! requested times, automatic if absent
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        integer, intent(in), optional :: n !! number of automatic time samples
        type(continuous_response_result) :: out
        out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),g,t,x0,.false.,n)
    end function step_zpk_real_complex_gain

    function step_zpk_complex_real_gain(z,p,g,t,x0,n) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! requested times, automatic if absent
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        integer, intent(in), optional :: n !! number of automatic time samples
        type(continuous_response_result) :: out
        out = response_zpk_core(z,p,cmplx(g,0.0_dp,dp),t,x0,.false.,n)
    end function step_zpk_complex_real_gain

    function impulse_zpk_real_complex_gain(z,p,g,t,x0,n) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! requested times, automatic if absent
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        integer, intent(in), optional :: n !! number of automatic time samples
        type(continuous_response_result) :: out
        out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),g,t,x0,.true.,n)
    end function impulse_zpk_real_complex_gain

    function impulse_zpk_complex_real_gain(z,p,g,t,x0,n) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        real(dp), intent(in), optional :: t(:) !! requested times, automatic if absent
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        integer, intent(in), optional :: n !! number of automatic time samples
        type(continuous_response_result) :: out
        out = response_zpk_core(z,p,cmplx(g,0.0_dp,dp),t,x0,.true.,n)
    end function impulse_zpk_complex_real_gain

    function lsim_zpk_real_complex_input(z,p,g,u,t,x0,interp) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolate inputs linearly
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(g,0.0_dp,dp),u,t,x0,interp)
    end function lsim_zpk_real_complex_input

    function lsim_zpk_complex_real_input(z,p,g,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolate inputs linearly
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(z,p,g,cmplx(u,0.0_dp,dp),t,x0,interp)
    end function lsim_zpk_complex_real_input

    function lsim_zpk_real_complex_gain(z,p,g,u,t,x0,interp) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolate inputs linearly
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),g,u,t,x0,interp)
    end function lsim_zpk_real_complex_gain

    function lsim_zpk_real_complex_gain_real_input(z,p,g,u,t,x0,interp) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: g !! gain
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolate inputs linearly
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),g,cmplx(u,0.0_dp,dp),t,x0,interp)
    end function lsim_zpk_real_complex_gain_real_input

    function lsim_zpk_complex_real_gain(z,p,g,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolate inputs linearly
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(z,p,cmplx(g,0.0_dp,dp),u,t,x0,interp)
    end function lsim_zpk_complex_real_gain

    function lsim_zpk_complex_real_gain_real_input(z,p,g,u,t,x0,interp) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: g !! gain
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! equally spaced sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in), optional :: interp !! interpolate inputs linearly
        type(continuous_simulation_result) :: out
        out = lsim_zpk_core(z,p,cmplx(g,0.0_dp,dp),cmplx(u,0.0_dp,dp),t,x0,interp)
    end function lsim_zpk_complex_real_gain_real_input

end module scifort_signal_continuous_simulation
