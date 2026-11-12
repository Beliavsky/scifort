! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Discrete-time linear time-invariant systems with SciPy-style conversions,
! simulation, and frequency analysis. The source representation is retained;
! polynomial conversion is only performed on explicit request.
!
! Polynomials are in descending powers of z. A time step of 1 represents
! an unspecified unit sampling interval. MIMO systems use state space.
module scifort_signal_dlti
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_state_space, only : tf2ss, ss2tf, zpk2ss, ss2zpk, &
        state_space_result, state_space_tf_result, state_space_success
    use scifort_signal_zpk_conversion, only : tf2zpk, zpk2tf, &
        filter_tf_result, filter_zpk_result, zpk_conversion_success
    use scifort_signal_discrete_simulation, only : dlsim, dstep, dimpulse, &
        dlsim_zpk, dstep_zpk, dimpulse_zpk, &
        discrete_simulation_result, discrete_response_result
    use scifort_signal_design, only : frequency_response
    use scifort_signal_lti_frequency, only : dfreqresp, dbode, dfreqresp_zpk, dbode_zpk, &
        dfreqresp_mimo, dbode_mimo, bode_result, mimo_frequency_response, mimo_bode_result
    implicit none
    private

    integer, parameter, public :: dlti_success = 0
    integer, parameter, public :: dlti_invalid_input = 1
    integer, parameter, public :: dlti_numerical_failure = 2
    integer, parameter, public :: dlti_transfer_function = 1
    integer, parameter, public :: dlti_zeros_poles_gain = 2
    integer, parameter, public :: dlti_state_space = 3

    ! Solutions of the discrete Lyapunov equations for a stable realization.
    ! wc = A wc A^H + B B^H and wo = A^H wo A + C^H C.
    type, public :: dlti_gramian_result
        complex(dp), allocatable :: wc(:,:) !! controllability Gramian
        complex(dp), allocatable :: wo(:,:) !! observability Gramian
        real(dp) :: h2_norm = 0.0_dp !! impulse-response energy norm, including D
        integer :: iterations = 0 !! maximum number of doubling iterations used
        integer :: status = dlti_invalid_input !! dlti_success or error status
    end type dlti_gramian_result

    type, public :: dlti_system
        ! Treat the stored coefficients as read-only after construction; use
        ! to_tf/to_zpk/to_ss or a constructor to change representations.
        real(dp) :: dt = 1.0_dp !! sampling interval in seconds (unit by default)
        integer :: status = dlti_invalid_input !! construction status
        integer :: representation = 0 !! one of dlti_* representation constants
        complex(dp), allocatable :: numerator(:) !! descending z powers
        complex(dp), allocatable :: denominator(:) !! descending z powers
        complex(dp), allocatable :: zeros(:) !! finite transfer-function zeros
        complex(dp), allocatable :: poles(:) !! finite transfer-function poles
        complex(dp) :: gain = (0.0_dp, 0.0_dp) !! ZPK gain
        complex(dp), allocatable :: a(:,:) !! state matrix (n,n)
        complex(dp), allocatable :: b(:,:) !! input matrix (n,nin)
        complex(dp), allocatable :: c(:,:) !! output matrix (nout,n)
        complex(dp), allocatable :: d(:,:) !! feedthrough (nout,nin)
    contains
        procedure :: to_tf => dlti_to_tf
        procedure :: to_zpk => dlti_to_zpk
        procedure :: to_ss => dlti_to_ss
        procedure :: freqresp => dlti_freqresp
        procedure :: bode => dlti_bode
        procedure :: freqresp_mimo => dlti_freqresp_mimo
        procedure :: bode_mimo => dlti_bode_mimo
        procedure :: output_real_vector => dlti_output_real_vector
        procedure :: output_complex_vector => dlti_output_complex_vector
        procedure :: output_real_matrix => dlti_output_real_matrix
        procedure :: output_complex_matrix => dlti_output_complex_matrix
        generic :: output => output_real_vector, output_complex_vector, &
            output_real_matrix, output_complex_matrix
        procedure :: step => dlti_step
        procedure :: impulse => dlti_impulse
        procedure :: ninputs => dlti_ninputs
        procedure :: noutputs => dlti_noutputs
        procedure :: nstates => dlti_nstates
        procedure :: stable => dlti_stable
        procedure :: series => dlti_series
        procedure :: parallel => dlti_parallel
        procedure :: append => dlti_append
        procedure :: select_io => dlti_select_io
        procedure :: feedback => dlti_feedback
        procedure :: controllability => dlti_controllability
        procedure :: observability => dlti_observability
        procedure :: gramians => dlti_gramians
        procedure :: minreal => dlti_minreal
        procedure :: balanced_truncate => dlti_balanced_truncate
    end type dlti_system

    ! Stable model-order reduction by square-root balanced truncation.
    type, public :: dlti_balanced_result
        type(dlti_system) :: model !! reduced state-space system
        real(dp), allocatable :: hsv(:) !! descending Hankel singular values
        real(dp) :: error_bound = 0.0_dp !! 2*sum of neglected values
        complex(dp), allocatable :: transform(:,:) !! n by retained order, x = transform*xr
        complex(dp), allocatable :: inverse(:,:) !! retained order by n, xr = inverse*x
        integer :: status = dlti_invalid_input !! success/invalid/numerical status
    end type dlti_balanced_result

    interface dlti_tf
        module procedure dlti_tf_real, dlti_tf_complex
    end interface
    interface dlti_zpk
        module procedure dlti_zpk_real, dlti_zpk_complex, &
            dlti_zpk_real_complex_gain, dlti_zpk_complex_real_gain
    end interface
    interface dlti_ss
        module procedure dlti_ss_real, dlti_ss_complex
    end interface
    ! Fortran cannot disambiguate dlti(num,den,dt) from dlti(z,p,k)
    ! when dt and k are both real. Use the named dlti_zpk constructor.
    interface dlti
        module procedure dlti_tf_real, dlti_tf_complex
        module procedure dlti_ss_real, dlti_ss_complex
    end interface
    public :: dlti, dlti_tf, dlti_zpk, dlti_ss

contains

    pure logical function valid_dt(dt) result(ok)
        real(dp), intent(in) :: dt !! sample interval
        ok = ieee_is_finite(dt) .and. dt > 0.0_dp
    end function valid_dt

    pure logical function finite_vector(values) result(ok)
        complex(dp), intent(in) :: values(:) !! complex entries
        ok = all(ieee_is_finite(real(values,dp))) .and. all(ieee_is_finite(aimag(values)))
    end function finite_vector

    pure logical function finite_matrix(values) result(ok)
        complex(dp), intent(in) :: values(:,:) !! complex entries
        ok = all(ieee_is_finite(real(values,dp))) .and. all(ieee_is_finite(aimag(values)))
    end function finite_matrix

    function dlti_tf_real(num,den,dt) result(sys)
        real(dp), intent(in) :: num(:) !! numerator coefficients
        real(dp), intent(in) :: den(:) !! denominator coefficients
        real(dp), optional, intent(in) :: dt !! sample interval, default 1
        type(dlti_system) :: sys
        sys = dlti_tf_complex(cmplx(num,0.0_dp,dp),cmplx(den,0.0_dp,dp),dt)
    end function dlti_tf_real

    function dlti_tf_complex(num,den,dt) result(sys)
        complex(dp), intent(in) :: num(:) !! numerator in descending z powers
        complex(dp), intent(in) :: den(:) !! denominator in descending z powers
        real(dp), optional, intent(in) :: dt !! positive sample interval
        type(dlti_system) :: sys
        real(dp) :: timestep
        integer :: lead

        timestep = 1.0_dp
        if (present(dt)) timestep = dt
        if (.not. valid_dt(timestep)) return
        if (size(num) == 0 .or. size(den) == 0) return
        if (.not. finite_vector(num) .or. .not. finite_vector(den)) return
        lead = 1
        do while (lead <= size(den))
            if (abs(den(lead)) /= 0.0_dp) exit
            lead = lead + 1
        end do
        if (lead > size(den)) return
        sys%numerator = num / den(lead)
        sys%denominator = den(lead:) / den(lead)
        if (.not. finite_vector(sys%numerator) .or. &
            .not. finite_vector(sys%denominator)) then
            sys = dlti_system()
            sys%status = dlti_numerical_failure
            return
        end if
        sys%dt = timestep
        sys%representation = dlti_transfer_function
        sys%status = dlti_success
    end function dlti_tf_complex

    function dlti_zpk_real(z,p,k,dt) result(sys)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        real(dp), intent(in) :: k !! real gain
        real(dp), optional, intent(in) :: dt !! sample interval
        type(dlti_system) :: sys
        sys = dlti_zpk_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp),dt)
    end function dlti_zpk_real

    function dlti_zpk_real_complex_gain(z,p,k,dt) result(sys)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        complex(dp), intent(in) :: k !! complex gain
        real(dp), optional, intent(in) :: dt !! sample interval
        type(dlti_system) :: sys
        sys = dlti_zpk_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),k,dt)
    end function dlti_zpk_real_complex_gain

    function dlti_zpk_complex_real_gain(z,p,k,dt) result(sys)
        complex(dp), intent(in) :: z(:) !! complex zeros
        complex(dp), intent(in) :: p(:) !! complex poles
        real(dp), intent(in) :: k !! real gain
        real(dp), optional, intent(in) :: dt !! sample interval
        type(dlti_system) :: sys
        sys = dlti_zpk_complex(z,p,cmplx(k,0.0_dp,dp),dt)
    end function dlti_zpk_complex_real_gain

    function dlti_zpk_complex(z,p,k,dt) result(sys)
        complex(dp), intent(in) :: z(:) !! complex zeros
        complex(dp), intent(in) :: p(:) !! complex poles
        complex(dp), intent(in) :: k !! complex gain
        real(dp), optional, intent(in) :: dt !! sample interval
        type(dlti_system) :: sys
        real(dp) :: timestep

        timestep = 1.0_dp
        if (present(dt)) timestep = dt
        if (.not. valid_dt(timestep)) return
        if (.not. finite_vector(z) .or. .not. finite_vector(p)) return
        if (.not. ieee_is_finite(real(k,dp)) .or. .not. ieee_is_finite(aimag(k))) return
        sys%zeros = z
        sys%poles = p
        sys%gain = k
        sys%dt = timestep
        sys%representation = dlti_zeros_poles_gain
        sys%status = dlti_success
    end function dlti_zpk_complex

    function dlti_ss_real(a,b,c,d,dt) result(sys)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), optional, intent(in) :: dt !! sample interval
        type(dlti_system) :: sys
        sys = dlti_ss_complex(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),dt)
    end function dlti_ss_real

    function dlti_ss_complex(a,b,c,d,dt) result(sys)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), optional, intent(in) :: dt !! sample interval
        type(dlti_system) :: sys
        real(dp) :: timestep
        integer :: n

        timestep = 1.0_dp
        if (present(dt)) timestep = dt
        if (.not. valid_dt(timestep)) return
        n = size(a,1)
        if (size(a,2) /= n .or. size(b,1) /= n .or. size(c,2) /= n) return
        if (size(b,2) < 1 .or. size(c,1) < 1) return
        if (size(d,1) /= size(c,1) .or. size(d,2) /= size(b,2)) return
        if (.not. finite_matrix(a) .or. .not. finite_matrix(b) .or. &
            .not. finite_matrix(c) .or. .not. finite_matrix(d)) return
        sys%a = a
        sys%b = b
        sys%c = c
        sys%d = d
        sys%dt = timestep
        sys%representation = dlti_state_space
        sys%status = dlti_success
    end function dlti_ss_complex

    integer function dlti_ninputs(self) result(count)
        class(dlti_system), intent(in) :: self !! system
        count = 0
        if (self%status /= dlti_success) return
        if (self%representation == dlti_state_space) then
            count = size(self%b,2)
        else
            count = 1
        end if
    end function dlti_ninputs

    integer function dlti_noutputs(self) result(count)
        class(dlti_system), intent(in) :: self !! system
        count = 0
        if (self%status /= dlti_success) return
        if (self%representation == dlti_state_space) then
            count = size(self%c,1)
        else
            count = 1
        end if
    end function dlti_noutputs

    integer function dlti_nstates(self) result(count)
        class(dlti_system), intent(in) :: self !! system
        type(dlti_system) :: ss
        count = 0
        if (self%status /= dlti_success) return
        if (self%representation == dlti_state_space) then
            count = size(self%a,1)
        else
            ss = self%to_ss()
            if (ss%status == dlti_success) count = size(ss%a,1)
        end if
    end function dlti_nstates

    function dlti_to_ss(self) result(out)
        class(dlti_system), intent(in) :: self !! original system
        type(dlti_system) :: out
        type(state_space_result) :: converted
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_state_space)
            out = self
            return
        case (dlti_transfer_function)
            converted = tf2ss(self%numerator,self%denominator)
        case (dlti_zeros_poles_gain)
            converted = zpk2ss(self%zeros,self%poles,self%gain)
        case default
            return
        end select
        if (converted%status /= state_space_success) then
            out%status = dlti_numerical_failure
            return
        end if
        out = dlti_ss(converted%a,converted%b,converted%c,converted%d,dt=self%dt)
    end function dlti_to_ss

    function dlti_to_tf(self) result(out)
        class(dlti_system), intent(in) :: self !! original system
        type(dlti_system) :: out
        type(filter_tf_result) :: converted
        type(state_space_tf_result) :: converted_ss
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_transfer_function)
            out = self
            return
        case (dlti_zeros_poles_gain)
            converted = zpk2tf(self%zeros,self%poles,self%gain)
            if (converted%status /= zpk_conversion_success) then
                out%status = dlti_numerical_failure
                return
            end if
            out = dlti_tf(converted%b,converted%a,dt=self%dt)
        case (dlti_state_space)
            if (self%ninputs() /= 1 .or. self%noutputs() /= 1) return
            converted_ss = ss2tf(self%a,self%b,self%c,self%d)
            if (converted_ss%status /= state_space_success) then
                out%status = dlti_numerical_failure
                return
            end if
            out = dlti_tf(converted_ss%b(1,:),converted_ss%a,dt=self%dt)
        end select
    end function dlti_to_tf

    function dlti_to_zpk(self) result(out)
        class(dlti_system), intent(in) :: self !! original system
        type(dlti_system) :: out
        type(filter_zpk_result) :: converted
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_zeros_poles_gain)
            out = self
            return
        case (dlti_transfer_function)
            converted = tf2zpk(self%numerator,self%denominator)
        case (dlti_state_space)
            if (self%ninputs() /= 1 .or. self%noutputs() /= 1) return
            converted = ss2zpk(self%a,self%b,self%c,self%d)
        case default
            return
        end select
        if (converted%status /= zpk_conversion_success) then
            out%status = dlti_numerical_failure
            return
        end if
        out = dlti_zpk(converted%zeros,converted%poles,converted%gain,dt=self%dt)
    end function dlti_to_zpk

    function dlti_freqresp(self,w,n,whole) result(out)
        class(dlti_system), intent(in) :: self !! system
        real(dp), optional, intent(in) :: w(:) !! angular frequencies in rad/sample
        integer, optional, intent(in) :: n !! automatic grid size
        logical, optional, intent(in) :: whole !! whole unit circle
        type(frequency_response) :: out
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_transfer_function)
            out = dfreqresp(self%numerator,self%denominator,w,n,whole,dt=self%dt)
        case (dlti_zeros_poles_gain)
            out = dfreqresp_zpk(self%zeros,self%poles,self%gain,w,n,whole,dt=self%dt)
        case (dlti_state_space)
            if (self%ninputs() /= 1 .or. self%noutputs() /= 1) return
            out = dfreqresp(self%a,self%b,self%c,self%d,w,n,whole,dt=self%dt)
        end select
    end function dlti_freqresp

    function dlti_bode(self,w,n) result(out)
        class(dlti_system), intent(in) :: self !! system
        real(dp), optional, intent(in) :: w(:) !! rad/sample grid
        integer, optional, intent(in) :: n !! automatic grid size
        type(bode_result) :: out
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_transfer_function)
            out = dbode(self%numerator,self%denominator,w,n,dt=self%dt)
        case (dlti_zeros_poles_gain)
            out = dbode_zpk(self%zeros,self%poles,self%gain,w,n,dt=self%dt)
        case (dlti_state_space)
            if (self%ninputs() /= 1 .or. self%noutputs() /= 1) return
            out = dbode(self%a,self%b,self%c,self%d,w,n,dt=self%dt)
        end select
    end function dlti_bode

    function dlti_freqresp_mimo(self,w,n,whole,solver,inputs,outputs) result(out)
        class(dlti_system), intent(in) :: self !! system
        real(dp), optional, intent(in) :: w(:) !! radians per sample
        integer, optional, intent(in) :: n !! automatic grid size
        logical, optional, intent(in) :: whole !! use whole unit circle
        character(*), optional, intent(in) :: solver !! direct, hessenberg, schur, auto
        integer, optional, intent(in) :: inputs(:) !! selected input channels, one-based
        integer, optional, intent(in) :: outputs(:) !! selected output channels, one-based
        type(mimo_frequency_response) :: out
        type(dlti_system) :: ss
        if (self%status /= dlti_success) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        out = dfreqresp_mimo(ss%a,ss%b,ss%c,ss%d,w,n,whole,dt=self%dt, &
            solver=solver,inputs=inputs,outputs=outputs)
    end function dlti_freqresp_mimo

    function dlti_bode_mimo(self,w,n,solver,inputs,outputs) result(out)
        class(dlti_system), intent(in) :: self !! system
        real(dp), optional, intent(in) :: w(:) !! radians per sample
        integer, optional, intent(in) :: n !! automatic grid size
        character(*), optional, intent(in) :: solver !! direct, hessenberg, schur, auto
        integer, optional, intent(in) :: inputs(:) !! selected input channels, one-based
        integer, optional, intent(in) :: outputs(:) !! selected output channels, one-based
        type(mimo_bode_result) :: out
        type(dlti_system) :: ss
        if (self%status /= dlti_success) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        out = dbode_mimo(ss%a,ss%b,ss%c,ss%d,w,n,dt=self%dt, &
            solver=solver,inputs=inputs,outputs=outputs)
    end function dlti_bode_mimo

    function dlti_output_complex_matrix(self,u,t,x0) result(out)
        class(dlti_system), intent(in) :: self !! system
        complex(dp), intent(in) :: u(:,:) !! inputs as (samples,channels)
        real(dp), optional, intent(in) :: t(:) !! input timestamps
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        if (self%status /= dlti_success) return
        if (size(u,2) /= self%ninputs()) return
        select case (self%representation)
        case (dlti_state_space)
            out = dlsim(self%a,self%b,self%c,self%d,u,self%dt,t,x0)
        case (dlti_transfer_function)
            out = dlsim(self%numerator,self%denominator,u(:,1),self%dt,t,x0)
        case (dlti_zeros_poles_gain)
            out = dlsim_zpk(self%zeros,self%poles,self%gain,u(:,1),self%dt,t,x0)
        end select
    end function dlti_output_complex_matrix

    function dlti_output_real_matrix(self,u,t,x0) result(out)
        class(dlti_system), intent(in) :: self !! system
        real(dp), intent(in) :: u(:,:) !! real inputs as (samples,channels)
        real(dp), optional, intent(in) :: t(:) !! input timestamps
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        out = self%output_complex_matrix(cmplx(u,0.0_dp,dp),t,x0)
    end function dlti_output_real_matrix

    function dlti_output_complex_vector(self,u,t,x0) result(out)
        class(dlti_system), intent(in) :: self !! system
        complex(dp), intent(in) :: u(:) !! SISO input samples
        real(dp), optional, intent(in) :: t(:) !! input timestamps
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        if (self%ninputs() /= 1) return
        out = self%output_complex_matrix(reshape(u,[size(u),1]),t,x0)
    end function dlti_output_complex_vector

    function dlti_output_real_vector(self,u,t,x0) result(out)
        class(dlti_system), intent(in) :: self !! system
        real(dp), intent(in) :: u(:) !! real SISO input samples
        real(dp), optional, intent(in) :: t(:) !! input timestamps
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        out = self%output_complex_vector(cmplx(u,0.0_dp,dp),t,x0)
    end function dlti_output_real_vector

    function dlti_step(self,n,t,x0) result(out)
        class(dlti_system), intent(in) :: self !! system
        integer, optional, intent(in) :: n !! response sample count (default 100)
        real(dp), optional, intent(in) :: t(:) !! response times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        type(discrete_response_result) :: out
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_state_space)
            out = dstep(self%a,self%b,self%c,self%d,self%dt,n,t,x0)
        case (dlti_transfer_function)
            out = dstep(self%numerator,self%denominator,self%dt,n,t,x0)
        case (dlti_zeros_poles_gain)
            out = dstep_zpk(self%zeros,self%poles,self%gain,self%dt,n,t,x0)
        end select
    end function dlti_step

    function dlti_impulse(self,n,t,x0) result(out)
        class(dlti_system), intent(in) :: self !! system
        integer, optional, intent(in) :: n !! response sample count (default 100)
        real(dp), optional, intent(in) :: t(:) !! response times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        type(discrete_response_result) :: out
        if (self%status /= dlti_success) return
        select case (self%representation)
        case (dlti_state_space)
            out = dimpulse(self%a,self%b,self%c,self%d,self%dt,n,t,x0)
        case (dlti_transfer_function)
            out = dimpulse(self%numerator,self%denominator,self%dt,n,t,x0)
        case (dlti_zeros_poles_gain)
            out = dimpulse_zpk(self%zeros,self%poles,self%gain,self%dt,n,t,x0)
        end select
    end function dlti_impulse

    function dlti_series(self,other) result(out)
        ! Cascade: the signal passes through OTHER first, then SELF.
        ! This matches scipy.signal.StateSpace multiplication self*other.
        class(dlti_system), intent(in) :: self !! downstream discrete system
        class(dlti_system), intent(in) :: other !! upstream discrete system
        type(dlti_system) :: out
        type(dlti_system) :: first, second
        complex(dp), allocatable :: aa(:,:), bb(:,:), cc(:,:), dd(:,:)
        integer :: n1, n2, nu, ny

        if (self%status /= dlti_success .or. other%status /= dlti_success) return
        if (self%dt /= other%dt) return
        if (self%ninputs() /= other%noutputs()) return
        first = self%to_ss()
        second = other%to_ss()
        if (first%status /= dlti_success .or. second%status /= dlti_success) return
        n1 = size(first%a,1)
        n2 = size(second%a,1)
        nu = other%ninputs()
        ny = self%noutputs()
        allocate(aa(n1+n2,n1+n2),bb(n1+n2,nu),cc(ny,n1+n2),dd(ny,nu))
        aa = (0.0_dp,0.0_dp)
        aa(:n1,:n1) = first%a
        aa(:n1,n1+1:) = matmul(first%b,second%c)
        aa(n1+1:,n1+1:) = second%a
        bb(:n1,:) = matmul(first%b,second%d)
        bb(n1+1:,:) = second%b
        cc(:,:n1) = first%c
        cc(:,n1+1:) = matmul(first%d,second%c)
        dd = matmul(first%d,second%d)
        out = dlti_ss(aa,bb,cc,dd,dt=self%dt)
    end function dlti_series

    function dlti_parallel(self,other) result(out)
        ! Sum two systems driven by the same input channels.
        class(dlti_system), intent(in) :: self !! first parallel system
        class(dlti_system), intent(in) :: other !! second parallel system
        type(dlti_system) :: out
        type(dlti_system) :: first, second
        complex(dp), allocatable :: aa(:,:), bb(:,:), cc(:,:), dd(:,:)
        integer :: n1, n2, nu, ny

        if (self%status /= dlti_success .or. other%status /= dlti_success) return
        if (self%dt /= other%dt) return
        if (self%ninputs() /= other%ninputs()) return
        if (self%noutputs() /= other%noutputs()) return
        first = self%to_ss()
        second = other%to_ss()
        if (first%status /= dlti_success .or. second%status /= dlti_success) return
        n1 = size(first%a,1)
        n2 = size(second%a,1)
        nu = self%ninputs()
        ny = self%noutputs()
        allocate(aa(n1+n2,n1+n2),bb(n1+n2,nu),cc(ny,n1+n2),dd(ny,nu))
        aa = (0.0_dp,0.0_dp)
        aa(:n1,:n1) = first%a
        aa(n1+1:,n1+1:) = second%a
        bb(:n1,:) = first%b
        bb(n1+1:,:) = second%b
        cc(:,:n1) = first%c
        cc(:,n1+1:) = second%c
        dd = first%d + second%d
        out = dlti_ss(aa,bb,cc,dd,dt=self%dt)
    end function dlti_parallel


    function dlti_append(self,other) result(out)
        ! Independent systems in block diagonal form, retaining all I/O channels.
        class(dlti_system), intent(in) :: self !! first system
        class(dlti_system), intent(in) :: other !! second system
        type(dlti_system) :: out
        type(dlti_system) :: first, second
        complex(dp), allocatable :: aa(:,:), bb(:,:), cc(:,:), dd(:,:)
        integer :: n1,n2,nu1,nu2,ny1,ny2

        if (self%status /= dlti_success .or. other%status /= dlti_success) return
        if (self%dt /= other%dt) return
        first = self%to_ss()
        second = other%to_ss()
        if (first%status /= dlti_success .or. second%status /= dlti_success) return
        n1 = size(first%a,1)
        n2 = size(second%a,1)
        nu1 = size(first%b,2)
        nu2 = size(second%b,2)
        ny1 = size(first%c,1)
        ny2 = size(second%c,1)
        allocate(aa(n1+n2,n1+n2),bb(n1+n2,nu1+nu2), &
            cc(ny1+ny2,n1+n2),dd(ny1+ny2,nu1+nu2))
        aa = (0.0_dp,0.0_dp)
        bb = (0.0_dp,0.0_dp)
        cc = (0.0_dp,0.0_dp)
        dd = (0.0_dp,0.0_dp)
        aa(:n1,:n1) = first%a
        aa(n1+1:,n1+1:) = second%a
        bb(:n1,:nu1) = first%b
        bb(n1+1:,nu1+1:) = second%b
        cc(:ny1,:n1) = first%c
        cc(ny1+1:,n1+1:) = second%c
        dd(:ny1,:nu1) = first%d
        dd(ny1+1:,nu1+1:) = second%d
        out = dlti_ss(aa,bb,cc,dd,self%dt)
    end function dlti_append

    function dlti_select_io(self,inputs,outputs) result(out)
        ! Select and/or reorder I/O channels; indices are one-based.
        class(dlti_system), intent(in) :: self !! original system
        integer, intent(in) :: inputs(:) !! requested input channels
        integer, intent(in) :: outputs(:) !! requested output channels
        type(dlti_system) :: out
        type(dlti_system) :: ss

        if (self%status /= dlti_success) return
        if (size(inputs) < 1 .or. size(outputs) < 1) return
        if (any(inputs < 1) .or. any(inputs > self%ninputs())) return
        if (any(outputs < 1) .or. any(outputs > self%noutputs())) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        out = dlti_ss(ss%a,ss%b(:,inputs),ss%c(outputs,:), &
            ss%d(outputs,inputs),dt=self%dt)
    end function dlti_select_io

    subroutine invert_complex_matrix(a,inverse,ok)
        ! Pivoted Gauss-Jordan inversion. This is used only for the small
        ! algebraic-loop matrix I - sign*Dk*Dp, not for state propagation.
        complex(dp), intent(in) :: a(:,:) !! square matrix
        complex(dp), intent(out) :: inverse(:,:) !! inverse, if successful
        logical, intent(out) :: ok !! false if singular or nonfinite
        complex(dp), allocatable :: work(:,:), temp_row(:), temp_inv(:)
        complex(dp) :: pivot, multiplier
        real(dp) :: scale, threshold
        integer :: n,i,j,k

        ok = .false.
        n = size(a,1)
        if (size(a,2) /= n .or. any(shape(inverse) /= [n,n])) return
        if (.not. finite_matrix(a)) return
        allocate(work(n,n),temp_row(n),temp_inv(n))
        work = a
        inverse = (0.0_dp,0.0_dp)
        do i=1,n
            inverse(i,i) = (1.0_dp,0.0_dp)
        end do
        if (n == 0) then
            ok = .true.
            return
        end if
        scale = maxval(abs(a))
        threshold = real(n,dp)*epsilon(1.0_dp)*max(1.0_dp,scale)
        do k=1,n
            j = k-1+maxloc(abs(work(k:n,k)),dim=1)
            if (abs(work(j,k)) <= threshold) return
            if (j /= k) then
                temp_row = work(k,:)
                work(k,:) = work(j,:)
                work(j,:) = temp_row
                temp_inv = inverse(k,:)
                inverse(k,:) = inverse(j,:)
                inverse(j,:) = temp_inv
            end if
            pivot = work(k,k)
            work(k,:) = work(k,:)/pivot
            inverse(k,:) = inverse(k,:)/pivot
            do i=1,n
                if (i == k) cycle
                multiplier = work(i,k)
                work(i,:) = work(i,:)-multiplier*work(k,:)
                inverse(i,:) = inverse(i,:)-multiplier*inverse(k,:)
            end do
        end do
        ok = finite_matrix(inverse)
    end subroutine invert_complex_matrix

    function dlti_feedback(self,controller,sign) result(out)
        ! Closed loop from the plant's external input r to plant output y:
        !     u = r + sign * K*y,       y = P*u.
        ! sign=-1 (default) is negative feedback. P is (ny,nu), K is (nu,ny).
        ! Nonzero Dp and Dk are handled by solving the algebraic loop.
        class(dlti_system), intent(in) :: self !! plant P
        class(dlti_system), intent(in) :: controller !! feedback system K
        integer, optional, intent(in) :: sign !! -1 (negative) or +1 (positive)
        type(dlti_system) :: out
        type(dlti_system) :: plant, ctrl
        complex(dp), allocatable :: loop(:,:), inverse(:,:), up(:,:), uc(:,:), &
            aa(:,:), bb(:,:), cc(:,:), dd(:,:)
        integer :: n1,n2,nu,ny,i,s
        logical :: ok

        if (self%status /= dlti_success .or. controller%status /= dlti_success) return
        if (self%dt /= controller%dt) return
        if (self%ninputs() /= controller%noutputs()) return
        if (self%noutputs() /= controller%ninputs()) return
        s = -1
        if (present(sign)) s = sign
        if (abs(s) /= 1) return
        plant = self%to_ss()
        ctrl = controller%to_ss()
        if (plant%status /= dlti_success .or. ctrl%status /= dlti_success) return
        n1 = size(plant%a,1)
        n2 = size(ctrl%a,1)
        nu = self%ninputs()
        ny = self%noutputs()
        allocate(loop(nu,nu),inverse(nu,nu),up(nu,n1),uc(nu,n2))
        loop = -real(s,dp)*matmul(ctrl%d,plant%d)
        do i=1,nu
            loop(i,i) = loop(i,i)+(1.0_dp,0.0_dp)
        end do
        call invert_complex_matrix(loop,inverse,ok)
        if (.not. ok) then
            out%status = dlti_numerical_failure
            return
        end if
        up = real(s,dp)*matmul(matmul(inverse,ctrl%d),plant%c)
        uc = real(s,dp)*matmul(inverse,ctrl%c)
        allocate(aa(n1+n2,n1+n2),bb(n1+n2,nu), &
            cc(ny,n1+n2),dd(ny,nu))
        aa(:n1,:n1) = plant%a+matmul(plant%b,up)
        aa(:n1,n1+1:) = matmul(plant%b,uc)
        aa(n1+1:,:n1) = matmul(ctrl%b,plant%c+matmul(plant%d,up))
        aa(n1+1:,n1+1:) = ctrl%a+matmul(matmul(ctrl%b,plant%d),uc)
        bb(:n1,:) = matmul(plant%b,inverse)
        bb(n1+1:,:) = matmul(matmul(ctrl%b,plant%d),inverse)
        cc(:,:n1) = plant%c+matmul(plant%d,up)
        cc(:,n1+1:) = matmul(plant%d,uc)
        dd = matmul(plant%d,inverse)
        out = dlti_ss(aa,bb,cc,dd,dt=self%dt)
        if (out%status /= dlti_success) out%status = dlti_numerical_failure
    end function dlti_feedback

    function dlti_controllability(self) result(matrix)
        ! [B, AB, A^2 B, ..., A^(n-1) B].
        class(dlti_system), intent(in) :: self !! discrete system
        complex(dp), allocatable :: matrix(:,:) !! n by n*nu
        type(dlti_system) :: ss
        complex(dp), allocatable :: power_b(:,:)
        integer :: n,nu,k

        if (self%status /= dlti_success) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        n = size(ss%a,1)
        nu = size(ss%b,2)
        allocate(matrix(n,n*nu))
        power_b = ss%b
        do k=0,n-1
            matrix(:,k*nu+1:(k+1)*nu) = power_b
            power_b = matmul(ss%a,power_b)
        end do
    end function dlti_controllability

    function dlti_observability(self) result(matrix)
        ! [C; CA; CA^2; ...; CA^(n-1)].
        class(dlti_system), intent(in) :: self !! discrete system
        complex(dp), allocatable :: matrix(:,:) !! n*ny by n
        type(dlti_system) :: ss
        complex(dp), allocatable :: power_c(:,:)
        integer :: n,ny,k

        if (self%status /= dlti_success) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        n = size(ss%a,1)
        ny = size(ss%c,1)
        allocate(matrix(n*ny,n))
        power_c = ss%c
        do k=0,n-1
            matrix(k*ny+1:(k+1)*ny,:) = power_c
            power_c = matmul(power_c,ss%a)
        end do
    end function dlti_observability

    subroutine doubling_gramian(a,q,tol,max_iter,w,iterations,ok)
        ! Compute W = A W A^H + Q by Smith's doubling recurrence.
        ! Repeated squaring requires O(log horizon) matrix products.
        complex(dp), intent(in) :: a(:,:) !! stable transition matrix
        complex(dp), intent(in) :: q(:,:) !! positive semidefinite right-hand side
        real(dp), intent(in) :: tol !! scaled error tolerance
        integer, intent(in) :: max_iter !! maximum doubling passes
        complex(dp), intent(out) :: w(:,:) !! Gramian on success
        integer, intent(out) :: iterations !! number of passes
        logical, intent(out) :: ok !! convergence status
        complex(dp), allocatable :: power(:,:), increment(:,:)
        real(dp) :: max_update, max_power
        integer :: k,n

        n = size(a,1)
        ok = .false.
        iterations = 0
        w = q
        if (n == 0) then
            ok = .true.
            return
        end if
        power = a
        do k=1,max_iter
            increment = matmul(matmul(power,w),transpose(conjg(power)))
            w = w + increment
            power = matmul(power,power)
            iterations = k
            if (.not. finite_matrix(w) .or. .not. finite_matrix(power)) return
            max_update = maxval(abs(increment))
            max_power = maxval(abs(power))
            if (max_update <= tol*max(1.0_dp,maxval(abs(w))) .and. &
                max_power <= sqrt(tol)) then
                w = 0.5_dp*(w+transpose(conjg(w)))
                ok = .true.
                return
            end if
        end do
    end subroutine doubling_gramian

    function dlti_gramians(self,tol,max_iter) result(out)
        ! Infinite-horizon discrete Gramians, and H2 norm from an impulse.
        ! Requires an asymptotically stable realization (including hidden modes).
        class(dlti_system), intent(in) :: self !! system to analyze
        real(dp), optional, intent(in) :: tol !! default 1e-12
        integer, optional, intent(in) :: max_iter !! default 60
        type(dlti_gramian_result) :: out
        type(dlti_system) :: ss
        complex(dp), allocatable :: q(:,:), ah(:,:), power(:,:)
        real(dp) :: tolerance, energy
        logical :: ok
        integer :: n,iter,iter2,limit,i

        if (self%status /= dlti_success) return
        tolerance = 1.0e-12_dp
        if (present(tol)) tolerance = tol
        limit = 60
        if (present(max_iter)) limit = max_iter
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance < 4.0_dp*epsilon(1.0_dp) .or. tolerance > 0.1_dp .or. limit < 1) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        n = size(ss%a,1)
        allocate(out%wc(n,n),out%wo(n,n))
        ah = transpose(conjg(ss%a))
        q = matmul(ss%b,transpose(conjg(ss%b)))
        call doubling_gramian(ss%a,q,tolerance,limit,out%wc,iter,ok)
        if (.not. ok) then
            out%status = dlti_numerical_failure
            return
        end if
        q = matmul(transpose(conjg(ss%c)),ss%c)
        call doubling_gramian(ah,q,tolerance,limit,out%wo,iter2,ok)
        if (.not. ok) then
            out%status = dlti_numerical_failure
            return
        end if
        power = matmul(matmul(ss%c,out%wc),transpose(conjg(ss%c)))+ &
            matmul(ss%d,transpose(conjg(ss%d)))
        energy = 0.0_dp
        do i=1,size(power,1)
            energy = energy + real(power(i,i),dp)
        end do
        if (.not. ieee_is_finite(energy) .or. energy < -100.0_dp*epsilon(energy)) then
            out%status = dlti_numerical_failure
            return
        end if
        out%h2_norm = sqrt(max(0.0_dp,energy))
        out%iterations = max(iter,iter2)
        out%status = dlti_success
    end function dlti_gramians


    function orthonormal_range(x,tol) result(q)
        ! Column-pivoted, twice-reorthogonalized modified Gram-Schmidt.
        ! Only numerical ranks well separated from the tolerance are reliable.
        complex(dp), intent(in) :: x(:,:) !! input columns
        real(dp), intent(in) :: tol !! relative rank cutoff
        complex(dp), allocatable :: q(:,:) !! orthonormal basis for range(X)
        complex(dp), allocatable :: work(:,:), basis(:,:), col(:)
        real(dp) :: initial_norm, greatest, projection_norm
        integer :: n,m,r,k,j,imax,pass

        n = size(x,1)
        m = size(x,2)
        allocate(q(n,0))
        if (n == 0 .or. m == 0) return
        if (.not. finite_matrix(x)) return
        initial_norm = maxval(maxval(abs(x),dim=1))
        if (initial_norm == 0.0_dp) return
        work = x
        allocate(basis(n,min(n,m)))
        r = 0
        do k=1,min(n,m)
            greatest = -1.0_dp
            imax = 1
            do j=1,m
                projection_norm = sqrt(sum(abs(work(:,j))**2))
                if (projection_norm > greatest) then
                    greatest = projection_norm
                    imax = j
                end if
            end do
            if (greatest <= tol*initial_norm) exit
            col = work(:,imax)/greatest
            r = r+1
            basis(:,r) = col
            do pass=1,2
                do j=1,m
                    work(:,j) = work(:,j)-col*dot_product(col,work(:,j))
                end do
            end do
        end do
        q = basis(:,:r)
    end function orthonormal_range

    function dlti_minreal(self,tol) result(out)
        ! Numerically minimal controllable/observable realization (Kalman
        ! decomposition). Transfers are preserved to the specified rank cutoff.
        ! This is not approximate balanced model-order reduction.
        class(dlti_system), intent(in) :: self !! original system
        real(dp), optional, intent(in) :: tol !! relative rank tolerance, default 1e-10
        type(dlti_system) :: out
        type(dlti_system) :: ss, reduced
        complex(dp), allocatable :: ctr(:,:), obs(:,:), qr(:,:), qo(:,:), &
            a(:,:), b(:,:), c(:,:)
        real(dp) :: threshold
        integer :: n

        if (self%status /= dlti_success) return
        threshold = 1.0e-10_dp
        if (present(tol)) threshold = tol
        if (.not. ieee_is_finite(threshold)) return
        if (threshold < 4.0_dp*epsilon(1.0_dp) .or. threshold >= 0.1_dp) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        n = size(ss%a,1)
        if (n == 0) then
            out = ss
            return
        end if
        ctr = ss%controllability()
        if (.not. allocated(ctr)) return
        if (.not. finite_matrix(ctr)) then
            out%status = dlti_numerical_failure
            return
        end if
        qr = orthonormal_range(ctr,threshold)
        a = matmul(matmul(transpose(conjg(qr)),ss%a),qr)
        b = matmul(transpose(conjg(qr)),ss%b)
        c = matmul(ss%c,qr)
        reduced = dlti_ss(a,b,c,ss%d,dt=self%dt)
        if (reduced%status /= dlti_success) then
            out%status = dlti_numerical_failure
            return
        end if
        obs = reduced%observability()
        if (.not. allocated(obs)) return
        if (.not. finite_matrix(obs)) then
            out%status = dlti_numerical_failure
            return
        end if
        qo = orthonormal_range(transpose(conjg(obs)),threshold)
        a = matmul(matmul(transpose(conjg(qo)),reduced%a),qo)
        b = matmul(transpose(conjg(qo)),reduced%b)
        c = matmul(reduced%c,qo)
        out = dlti_ss(a,b,c,ss%d,dt=self%dt)
        if (out%status /= dlti_success) out%status = dlti_numerical_failure
    end function dlti_minreal

    subroutine hermitian_jacobi(a,values,vectors,ok)
        ! Complex Hermitian Jacobi diagonalization, with no external LAPACK.
        ! Each plane rotation first makes the off-diagonal entry real.
        complex(dp), intent(in) :: a(:,:) !! Hermitian input matrix
        real(dp), intent(out) :: values(:) !! eigenvalues, decreasing order
        complex(dp), intent(out) :: vectors(:,:) !! orthonormal eigenvectors
        logical, intent(out) :: ok !! convergence and finite-input flag
        complex(dp), allocatable :: h(:,:)
        complex(dp) :: phase, hp, hq, vp, vq
        real(dp) :: scale, threshold, rho, tau, t, cs, sn, ap, aq, swapval
        integer :: n, i, j, k, p, q, pass, pick

        ok = .false.
        n = size(a,1)
        if (size(a,2) /= n .or. size(values) /= n) return
        if (any(shape(vectors) /= [n,n])) return
        if (.not. finite_matrix(a)) return
        vectors = (0.0_dp,0.0_dp)
        values = 0.0_dp
        do i=1,n
            vectors(i,i) = (1.0_dp,0.0_dp)
        end do
        if (n == 0) then
            ok = .true.
            return
        end if
        h = 0.5_dp*(a+transpose(conjg(a)))
        scale = maxval(abs(h))
        threshold = 16.0_dp*epsilon(1.0_dp)*scale
        if (scale == 0.0_dp) then
            ok = .true.
            return
        end if
        do pass=1,100
            rho = 0.0_dp
            do p=1,n-1
                do q=p+1,n
                    rho = max(rho,abs(h(p,q)))
                    if (abs(h(p,q)) <= threshold) cycle
                    phase = h(p,q)/abs(h(p,q))
                    tau = (real(h(q,q),dp)-real(h(p,p),dp))/(2.0_dp*abs(h(p,q)))
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
                    else
                        t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
                    end if
                    cs = 1.0_dp/sqrt(1.0_dp+t*t)
                    sn = cs*t
                    ap = real(h(p,p),dp)
                    aq = real(h(q,q),dp)
                    h(p,p) = cmplx(ap-t*abs(h(p,q)),0.0_dp,dp)
                    h(q,q) = cmplx(aq+t*abs(h(p,q)),0.0_dp,dp)
                    h(p,q) = (0.0_dp,0.0_dp)
                    h(q,p) = (0.0_dp,0.0_dp)
                    do k=1,n
                        if (k == p .or. k == q) cycle
                        hp = h(k,p)
                        hq = h(k,q)
                        h(k,p) = cs*hp-sn*conjg(phase)*hq
                        h(k,q) = sn*hp+cs*conjg(phase)*hq
                        h(p,k) = conjg(h(k,p))
                        h(q,k) = conjg(h(k,q))
                    end do
                    do k=1,n
                        vp = vectors(k,p)
                        vq = vectors(k,q)
                        vectors(k,p) = cs*vp-sn*conjg(phase)*vq
                        vectors(k,q) = sn*vp+cs*conjg(phase)*vq
                    end do
                end do
            end do
            if (rho <= threshold) exit
        end do
        if (pass > 100) return
        do i=1,n
            values(i) = real(h(i,i),dp)
        end do
        do i=1,n-1
            pick = i
            do j=i+1,n
                if (values(j) > values(pick)) pick = j
            end do
            if (pick == i) cycle
            swapval = values(i)
            values(i) = values(pick)
            values(pick) = swapval
            do k=1,n
                vp = vectors(k,i)
                vectors(k,i) = vectors(k,pick)
                vectors(k,pick) = vp
            end do
        end do
        ok = all(ieee_is_finite(values)) .and. finite_matrix(vectors)
    end subroutine hermitian_jacobi

    function dlti_balanced_truncate(self,order,tol) result(out)
        ! Stable square-root balanced truncation. For a minimal realization,
        ! ||G-Gr||_infinity <= 2*sum(hsv(order+1:)).
        ! This algorithm forms the SVD through an eigensystem of H^H*H;
        ! tiny Hankel singular values may therefore be lost. Reject
        ! borderline ranks rather than claim unreliable reductions.
        class(dlti_system), intent(in) :: self !! original discrete system
        integer, intent(in) :: order !! retained state count, 0..nstates
        real(dp), optional, intent(in) :: tol !! relative singular-value cutoff
        type(dlti_balanced_result) :: out
        type(dlti_system) :: ss
        type(dlti_gramian_result) :: gr
        complex(dp), allocatable :: uc(:,:), uo(:,:), lc(:,:), lo(:,:), &
            overlap(:,:), v(:,:), u(:,:), transform(:,:), inverse(:,:), &
            ah(:,:), ab(:,:), ac(:,:)
        real(dp), allocatable :: ec(:),eo(:),ev(:),singular(:)
        real(dp) :: cutoff, scale, small
        integer :: n,j,rank_c,rank_o,effective
        logical :: ok

        if (self%status /= dlti_success) return
        cutoff = 1.0e-10_dp
        if (present(tol)) cutoff = tol
        if (.not. ieee_is_finite(cutoff)) return
        if (cutoff < 32.0_dp*epsilon(1.0_dp) .or. cutoff >= 0.1_dp) return
        ss = self%to_ss()
        if (ss%status /= dlti_success) return
        n = size(ss%a,1)
        if (order < 0 .or. order > n) return
        if (n == 0) then
            out%hsv = [real(dp) ::]
            out%model = ss
            allocate(out%transform(0,0),out%inverse(0,0))
            out%status = dlti_success
            return
        end if
        gr = ss%gramians(tol=1.0e-13_dp)
        if (gr%status /= dlti_success) then
            out%status = dlti_numerical_failure
            return
        end if
        allocate(ec(n),eo(n),uc(n,n),uo(n,n))
        call hermitian_jacobi(gr%wc,ec,uc,ok)
        if (.not. ok) then
            out%status = dlti_numerical_failure
            return
        end if
        call hermitian_jacobi(gr%wo,eo,uo,ok)
        if (.not. ok) then
            out%status = dlti_numerical_failure
            return
        end if
        scale = max(1.0e-300_dp,max(maxval(abs(ec)),maxval(abs(eo))))
        small = 256.0_dp*epsilon(1.0_dp)*scale
        if (minval(ec) < -small .or. minval(eo) < -small) then
            out%status = dlti_numerical_failure
            return
        end if
        rank_c = count(ec > cutoff*maxval(abs(ec)))
        rank_o = count(eo > cutoff*maxval(abs(eo)))
        if (rank_c == 0 .or. rank_o == 0) then
            out%hsv = spread(0.0_dp,1,n)
            if (order /= 0) then
                out%status = dlti_numerical_failure
                return
            end if
            out%model = dlti_ss(ss%a(:0,:0),ss%b(:0,:),ss%c(:,:0),ss%d,ss%dt)
            allocate(out%transform(n,0),out%inverse(0,n))
            out%status = out%model%status
            return
        end if
        lc = uc(:,:rank_c)
        lo = uo(:,:rank_o)
        do j=1,rank_c
            lc(:,j) = lc(:,j)*sqrt(max(0.0_dp,ec(j)))
        end do
        do j=1,rank_o
            lo(:,j) = lo(:,j)*sqrt(max(0.0_dp,eo(j)))
        end do
        overlap = matmul(transpose(conjg(lo)),lc)
        ah = matmul(transpose(conjg(overlap)),overlap)
        allocate(ev(rank_c),v(rank_c,rank_c),singular(rank_c))
        call hermitian_jacobi(ah,ev,v,ok)
        if (.not. ok) then
            out%status = dlti_numerical_failure
            return
        end if
        if (minval(ev) < -256.0_dp*epsilon(1.0_dp)*maxval(abs(ev))) then
            out%status = dlti_numerical_failure
            return
        end if
        singular = sqrt(max(0.0_dp,ev))
        allocate(out%hsv(n))
        out%hsv = 0.0_dp
        out%hsv(:rank_c) = singular
        out%error_bound = 2.0_dp*sum(out%hsv(order+1:))
        effective = count(singular > cutoff*maxval(singular))
        if (order > effective) then
            out%status = dlti_numerical_failure
            return
        end if
        if (order == 0) then
            out%model = dlti_ss(ss%a(:0,:0),ss%b(:0,:),ss%c(:,:0),ss%d,ss%dt)
            allocate(out%transform(n,0),out%inverse(0,n))
            out%status = out%model%status
            return
        end if
        allocate(transform(n,order),inverse(order,n),u(rank_o,order))
        do j=1,order
            transform(:,j) = matmul(lc,v(:,j))/sqrt(singular(j))
            u(:,j) = matmul(overlap,v(:,j))/singular(j)
            inverse(j,:) = matmul(conjg(u(:,j)),transpose(conjg(lo)))/sqrt(singular(j))
        end do
        if (.not. finite_matrix(transform) .or. .not. finite_matrix(inverse)) then
            out%status = dlti_numerical_failure
            return
        end if
        ah = matmul(matmul(inverse,ss%a),transform)
        ab = matmul(inverse,ss%b)
        ac = matmul(ss%c,transform)
        out%model = dlti_ss(ah,ab,ac,ss%d,ss%dt)
        out%status = out%model%status
        if (out%status == dlti_success) then
            out%transform = transform
            out%inverse = inverse
        end if
    end function dlti_balanced_truncate

    logical function dlti_stable(self,tol) result(ok)
        class(dlti_system), intent(in) :: self !! system
        real(dp), optional, intent(in) :: tol !! subtract from unit-circle threshold
        type(dlti_system) :: zpk
        real(dp) :: margin
        ok = .false.
        if (self%status /= dlti_success) return
        margin = 0.0_dp
        if (present(tol)) margin = tol
        if (.not. ieee_is_finite(margin) .or. margin < 0.0_dp .or. margin >= 1.0_dp) return
        if (self%representation == dlti_state_space) then
            ! Internal stability includes uncontrollable and unobservable modes.
            ! ss2tf returns det(z*I-A), also for MIMO realizations.
            block
                type(state_space_tf_result) :: tf
                type(filter_zpk_result) :: roots
                complex(dp) :: unit_num(1)
                tf = ss2tf(self%a,self%b,self%c,self%d)
                if (tf%status /= state_space_success) return
                unit_num(1) = (1.0_dp,0.0_dp)
                roots = tf2zpk(unit_num,tf%a)
                if (roots%status /= zpk_conversion_success) return
                ok = all(abs(roots%poles) < 1.0_dp-margin)
            end block
        else
            zpk = self%to_zpk()
            if (zpk%status /= dlti_success) return
            ok = all(abs(zpk%poles) < 1.0_dp-margin)
        end if
    end function dlti_stable

end module scifort_signal_dlti
