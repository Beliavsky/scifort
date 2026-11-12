! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Continuous-time LTI objects. Algebraic state-space conversion and
! interconnection reuse the discrete LTI type; this internal representation
! has dt=1 only as a placeholder and never implies sampled dynamics.
! Simulation, stability, frequency response, and discretization use analog
! routines exclusively.
module scifort_signal_clti
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system, dlti_tf, dlti_zpk, dlti_ss, &
        dlti_success, dlti_gramian_result, dlti_balanced_result, &
        dlti_transfer_function, dlti_zeros_poles_gain, dlti_state_space
    use scifort_signal_design, only : frequency_response
    use scifort_signal_lti_frequency, only : freqresp, bode, freqresp_zpk, bode_zpk, &
        freqresp_mimo, bode_mimo, bode_result, mimo_frequency_response, mimo_bode_result
    use scifort_signal_continuous_simulation, only : lsim, step, impulse, &
        continuous_simulation_result, continuous_response_result
    use scifort_signal_cont2discrete, only : cont2discrete, discrete_state_space_result, &
        cont2discrete_success
    use scifort_signal_state_space, only : ss2tf, state_space_tf_result, state_space_success
    use scifort_signal_zpk_conversion, only : tf2zpk, filter_zpk_result, zpk_conversion_success
    implicit none
    private

    integer, parameter, public :: clti_success = 0
    integer, parameter, public :: clti_invalid_input = 1
    integer, parameter, public :: clti_numerical_failure = 2

    ! Continuous-time infinite-horizon solutions:
    ! A Wc + Wc A^H + B B^H = 0
    ! A^H Wo + Wo A + C^H C = 0
    type, public :: clti_gramian_result
        complex(dp), allocatable :: wc(:,:) !! controllability Gramian
        complex(dp), allocatable :: wo(:,:) !! observability Gramian
        real(dp) :: h2_norm = 0.0_dp !! infinity for nonzero direct feedthrough
        integer :: iterations = 0 !! discrete doubling iteration count
        integer :: status = clti_invalid_input !! success/invalid/numerical failure
    end type clti_gramian_result

    type, public :: clti_system
        private
        type(dlti_system) :: core !! shared polynomial/state-space representation
    contains
        procedure :: status => clti_status
        procedure :: nstates => clti_nstates
        procedure :: ninputs => clti_ninputs
        procedure :: noutputs => clti_noutputs
        procedure :: to_ss => clti_to_ss
        procedure :: to_tf => clti_to_tf
        procedure :: to_zpk => clti_to_zpk
        procedure :: freqresp => clti_freqresp
        procedure :: bode => clti_bode
        procedure :: freqresp_mimo => clti_freqresp_mimo
        procedure :: bode_mimo => clti_bode_mimo
        procedure :: output_real_vector => clti_output_real_vector
        procedure :: output_complex_vector => clti_output_complex_vector
        procedure :: output_real_matrix => clti_output_real_matrix
        procedure :: output_complex_matrix => clti_output_complex_matrix
        generic :: output => output_real_vector,output_complex_vector, &
            output_real_matrix,output_complex_matrix
        procedure :: step => clti_step
        procedure :: impulse => clti_impulse
        procedure :: to_discrete => clti_to_discrete
        procedure :: stable => clti_stable
        procedure :: series => clti_series
        procedure :: parallel => clti_parallel
        procedure :: append => clti_append
        procedure :: feedback => clti_feedback
        procedure :: select_io => clti_select_io
        procedure :: minreal => clti_minreal
        procedure :: controllability => clti_controllability
        procedure :: observability => clti_observability
        procedure :: gramians => clti_gramians
        procedure :: balanced_truncate => clti_balanced_truncate
    end type clti_system

    type, public :: clti_balanced_result
        type(clti_system) :: model !! reduced continuous-time model
        real(dp), allocatable :: hsv(:) !! descending Hankel singular values
        real(dp) :: error_bound = 0.0_dp !! H-infinity error bound
        integer :: status = clti_invalid_input !! success/invalid/numerical failure
    end type clti_balanced_result

    interface clti_tf
        module procedure clti_tf_real,clti_tf_complex
    end interface
    interface clti_zpk
        module procedure clti_zpk_real,clti_zpk_complex,clti_zpk_real_complex_gain, &
            clti_zpk_complex_real_gain
    end interface
    interface clti_ss
        module procedure clti_ss_real,clti_ss_complex
    end interface
    interface clti
        module procedure clti_tf_real,clti_tf_complex,clti_ss_real,clti_ss_complex
    end interface
    public :: clti,clti_tf,clti_zpk,clti_ss

contains

    function clti_tf_real(num,den) result(out)
        real(dp), intent(in) :: num(:) !! analog numerator, descending powers of s
        real(dp), intent(in) :: den(:) !! analog denominator, descending powers of s
        type(clti_system) :: out
        out%core = dlti_tf(num,den)
    end function clti_tf_real

    function clti_tf_complex(num,den) result(out)
        complex(dp), intent(in) :: num(:) !! analog numerator, descending powers of s
        complex(dp), intent(in) :: den(:) !! analog denominator, descending powers of s
        type(clti_system) :: out
        out%core = dlti_tf(num,den)
    end function clti_tf_complex

    function clti_zpk_real(z,p,k) result(out)
        real(dp), intent(in) :: z(:) !! real Laplace-plane zeros
        real(dp), intent(in) :: p(:) !! real Laplace-plane poles
        real(dp), intent(in) :: k !! real gain
        type(clti_system) :: out
        out%core = dlti_zpk(z,p,k)
    end function clti_zpk_real

    function clti_zpk_complex(z,p,k) result(out)
        complex(dp), intent(in) :: z(:) !! complex Laplace-plane zeros
        complex(dp), intent(in) :: p(:) !! complex Laplace-plane poles
        complex(dp), intent(in) :: k !! complex gain
        type(clti_system) :: out
        out%core = dlti_zpk(z,p,k)
    end function clti_zpk_complex

    function clti_zpk_real_complex_gain(z,p,k) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        complex(dp), intent(in) :: k !! complex gain
        type(clti_system) :: out
        out%core = dlti_zpk(z,p,k)
    end function clti_zpk_real_complex_gain

    function clti_zpk_complex_real_gain(z,p,k) result(out)
        complex(dp), intent(in) :: z(:) !! complex zeros
        complex(dp), intent(in) :: p(:) !! complex poles
        real(dp), intent(in) :: k !! real gain
        type(clti_system) :: out
        out%core = dlti_zpk(z,p,k)
    end function clti_zpk_complex_real_gain

    function clti_ss_real(a,b,c,d) result(out)
        real(dp), intent(in) :: a(:,:) !! state generator (n,n)
        real(dp), intent(in) :: b(:,:) !! input map (n,nu)
        real(dp), intent(in) :: c(:,:) !! output map (ny,n)
        real(dp), intent(in) :: d(:,:) !! direct feedthrough (ny,nu)
        type(clti_system) :: out
        out%core = dlti_ss(a,b,c,d)
    end function clti_ss_real

    function clti_ss_complex(a,b,c,d) result(out)
        complex(dp), intent(in) :: a(:,:) !! state generator (n,n)
        complex(dp), intent(in) :: b(:,:) !! input map (n,nu)
        complex(dp), intent(in) :: c(:,:) !! output map (ny,n)
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough (ny,nu)
        type(clti_system) :: out
        out%core = dlti_ss(a,b,c,d)
    end function clti_ss_complex

    integer function clti_status(self) result(status)
        class(clti_system), intent(in) :: self !! continuous system
        status = self%core%status
    end function clti_status

    integer function clti_ninputs(self) result(n)
        class(clti_system), intent(in) :: self !! continuous system
        n = self%core%ninputs()
    end function clti_ninputs

    integer function clti_noutputs(self) result(n)
        class(clti_system), intent(in) :: self !! continuous system
        n = self%core%noutputs()
    end function clti_noutputs

    integer function clti_nstates(self) result(n)
        class(clti_system), intent(in) :: self !! continuous system
        n = self%core%nstates()
    end function clti_nstates

    function clti_to_ss(self) result(out)
        class(clti_system), intent(in) :: self !! source system
        type(clti_system) :: out
        out%core = self%core%to_ss()
    end function clti_to_ss

    function clti_to_tf(self) result(out)
        class(clti_system), intent(in) :: self !! source SISO system
        type(clti_system) :: out
        out%core = self%core%to_tf()
    end function clti_to_tf

    function clti_to_zpk(self) result(out)
        class(clti_system), intent(in) :: self !! source SISO system
        type(clti_system) :: out
        out%core = self%core%to_zpk()
    end function clti_to_zpk

    function clti_freqresp(self,w,n) result(out)
        class(clti_system), intent(in) :: self !! analog SISO system
        real(dp), optional, intent(in) :: w(:) !! frequencies in rad/s
        integer, optional, intent(in) :: n !! automatic log-frequency grid size
        type(frequency_response) :: out
        if (self%status() /= clti_success) return
        select case(self%core%representation)
        case(dlti_transfer_function)
            out = freqresp(self%core%numerator,self%core%denominator,w,n)
        case(dlti_zeros_poles_gain)
            out = freqresp_zpk(self%core%zeros,self%core%poles,self%core%gain,w,n)
        case(dlti_state_space)
            if (self%ninputs() /= 1 .or. self%noutputs() /= 1) return
            out = freqresp(self%core%a,self%core%b,self%core%c,self%core%d,w,n)
        end select
    end function clti_freqresp

    function clti_bode(self,w,n) result(out)
        class(clti_system), intent(in) :: self !! analog SISO system
        real(dp), optional, intent(in) :: w(:) !! angular frequencies, rad/s
        integer, optional, intent(in) :: n !! logarithmic-grid size
        type(bode_result) :: out
        if (self%status() /= clti_success) return
        select case(self%core%representation)
        case(dlti_transfer_function)
            out = bode(self%core%numerator,self%core%denominator,w,n)
        case(dlti_zeros_poles_gain)
            out = bode_zpk(self%core%zeros,self%core%poles,self%core%gain,w,n)
        case(dlti_state_space)
            if (self%ninputs() /= 1 .or. self%noutputs() /= 1) return
            out = bode(self%core%a,self%core%b,self%core%c,self%core%d,w,n)
        end select
    end function clti_bode

    function clti_freqresp_mimo(self,w,n,solver,inputs,outputs) result(out)
        class(clti_system), intent(in) :: self !! analog system
        real(dp), optional, intent(in) :: w(:) !! angular frequencies, rad/s
        integer, optional, intent(in) :: n !! automatic grid size
        character(*), optional, intent(in) :: solver !! direct/hessenberg/schur/auto
        integer, optional, intent(in) :: inputs(:) !! one-based input selection
        integer, optional, intent(in) :: outputs(:) !! one-based output selection
        type(mimo_frequency_response) :: out
        type(clti_system) :: ss
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        out = freqresp_mimo(ss%core%a,ss%core%b,ss%core%c,ss%core%d,w,n, &
            solver=solver,inputs=inputs,outputs=outputs)
    end function clti_freqresp_mimo

    function clti_bode_mimo(self,w,n,solver,inputs,outputs) result(out)
        class(clti_system), intent(in) :: self !! analog system
        real(dp), optional, intent(in) :: w(:) !! angular frequencies, rad/s
        integer, optional, intent(in) :: n !! automatic grid size
        character(*), optional, intent(in) :: solver !! direct/hessenberg/schur/auto
        integer, optional, intent(in) :: inputs(:) !! one-based input selection
        integer, optional, intent(in) :: outputs(:) !! one-based output selection
        type(mimo_bode_result) :: out
        type(clti_system) :: ss
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        out = bode_mimo(ss%core%a,ss%core%b,ss%core%c,ss%core%d,w,n, &
            solver=solver,inputs=inputs,outputs=outputs)
    end function clti_bode_mimo

    function clti_output_complex_matrix(self,u,t,x0,interp) result(out)
        class(clti_system), intent(in) :: self !! analog system
        complex(dp), intent(in) :: u(:,:) !! samples by input channel
        real(dp), intent(in) :: t(:) !! nonnegative, equally spaced time grid
        complex(dp), optional, intent(in) :: x0(:) !! initial continuous state
        logical, optional, intent(in) :: interp !! linear (default) or zero-order hold
        type(continuous_simulation_result) :: out
        type(clti_system) :: ss
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        out = lsim(ss%core%a,ss%core%b,ss%core%c,ss%core%d,u,t,x0,interp)
    end function clti_output_complex_matrix

    function clti_output_real_matrix(self,u,t,x0,interp) result(out)
        class(clti_system), intent(in) :: self !! analog system
        real(dp), intent(in) :: u(:,:) !! sampled inputs, time first
        real(dp), intent(in) :: t(:) !! nonnegative equally spaced times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        logical, optional, intent(in) :: interp !! linear or stepwise input
        type(continuous_simulation_result) :: out
        out = self%output_complex_matrix(cmplx(u,0.0_dp,dp),t,x0,interp)
    end function clti_output_real_matrix

    function clti_output_complex_vector(self,u,t,x0,interp) result(out)
        class(clti_system), intent(in) :: self !! SISO analog system
        complex(dp), intent(in) :: u(:) !! sampled scalar input
        real(dp), intent(in) :: t(:) !! uniformly spaced times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        logical, optional, intent(in) :: interp !! interpolation policy
        type(continuous_simulation_result) :: out
        if (self%ninputs() /= 1) return
        out = self%output_complex_matrix(reshape(u,[size(u),1]),t,x0,interp)
    end function clti_output_complex_vector

    function clti_output_real_vector(self,u,t,x0,interp) result(out)
        class(clti_system), intent(in) :: self !! SISO analog system
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: t(:) !! uniformly spaced times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        logical, optional, intent(in) :: interp !! interpolation policy
        type(continuous_simulation_result) :: out
        if (self%ninputs() /= 1) return
        out = self%output_complex_matrix(reshape(cmplx(u,0.0_dp,dp),[size(u),1]),t,x0,interp)
    end function clti_output_real_vector

    function clti_step(self,t,x0,n) result(out)
        class(clti_system), intent(in) :: self !! analog system
        real(dp), optional, intent(in) :: t(:) !! response sampling times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        integer, optional, intent(in) :: n !! automatic grid size
        type(continuous_response_result) :: out
        type(clti_system) :: ss
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        out = step(ss%core%a,ss%core%b,ss%core%c,ss%core%d,t,x0,n)
    end function clti_step

    function clti_impulse(self,t,x0,n) result(out)
        class(clti_system), intent(in) :: self !! analog system
        real(dp), optional, intent(in) :: t(:) !! response sampling times
        complex(dp), optional, intent(in) :: x0(:) !! initial state
        integer, optional, intent(in) :: n !! automatic grid size
        type(continuous_response_result) :: out
        type(clti_system) :: ss
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        out = impulse(ss%core%a,ss%core%b,ss%core%c,ss%core%d,t,x0,n)
    end function clti_impulse

    function clti_to_discrete(self,dt,method,alpha) result(out)
        class(clti_system), intent(in) :: self !! analog system
        real(dp), intent(in) :: dt !! positive sampling interval in seconds
        character(*), optional, intent(in) :: method !! zoh, foh, impulse, gbt and variants
        real(dp), optional, intent(in) :: alpha !! generalized bilinear parameter
        type(dlti_system) :: out
        type(clti_system) :: ss
        type(discrete_state_space_result) :: disc
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        disc = cont2discrete(ss%core%a,ss%core%b,ss%core%c,ss%core%d,dt,method,alpha)
        if (disc%status /= cont2discrete_success) return
        out = dlti_ss(disc%a,disc%b,disc%c,disc%d,dt)
    end function clti_to_discrete

    logical function clti_stable(self,tol) result(ok)
        ! Continuous stability requires strictly negative real pole parts.
        class(clti_system), intent(in) :: self !! analog system
        real(dp), optional, intent(in) :: tol !! nonnegative real-part margin
        type(clti_system) :: ss
        type(state_space_tf_result) :: poly
        type(filter_zpk_result) :: zp
        real(dp) :: margin
        complex(dp) :: one(1)
        ok = .false.
        margin = 0.0_dp
        if (present(tol)) margin = tol
        if (.not. ieee_is_finite(margin) .or. margin < 0.0_dp) return
        if (self%status() /= clti_success) return
        if (self%core%representation == dlti_state_space) then
            ss = self%to_ss()
            if (ss%status() /= clti_success) return
            if (size(ss%core%a,1) == 0) then
                ok = .true.
                return
            end if
            poly = ss2tf(ss%core%a,ss%core%b,ss%core%c,ss%core%d)
            if (poly%status /= state_space_success) return
            one = (1.0_dp,0.0_dp)
            zp = tf2zpk(one,poly%a)
        else
            ss = self%to_zpk()
            if (ss%status() /= clti_success) return
            ok = all(real(ss%core%poles,dp) < -margin)
            return
        end if
        if (zp%status /= zpk_conversion_success) return
        ok = all(real(zp%poles,dp) < -margin)
    end function clti_stable

    function clti_series(self,other) result(out)
        class(clti_system), intent(in) :: self !! downstream system
        class(clti_system), intent(in) :: other !! upstream system
        type(clti_system) :: out
        out%core = self%core%series(other%core)
    end function clti_series

    function clti_parallel(self,other) result(out)
        class(clti_system), intent(in) :: self !! first parallel system
        class(clti_system), intent(in) :: other !! second parallel system
        type(clti_system) :: out
        out%core = self%core%parallel(other%core)
    end function clti_parallel

    function clti_append(self,other) result(out)
        class(clti_system), intent(in) :: self !! first independent system
        class(clti_system), intent(in) :: other !! second independent system
        type(clti_system) :: out
        out%core = self%core%append(other%core)
    end function clti_append

    function clti_feedback(self,controller,sign) result(out)
        class(clti_system), intent(in) :: self !! plant
        class(clti_system), intent(in) :: controller !! feedback controller
        integer, optional, intent(in) :: sign !! negative (-1) or positive (+1)
        type(clti_system) :: out
        out%core = self%core%feedback(controller%core,sign)
    end function clti_feedback

    function clti_select_io(self,inputs,outputs) result(out)
        class(clti_system), intent(in) :: self !! MIMO analog system
        integer, intent(in) :: inputs(:) !! selected one-based inputs
        integer, intent(in) :: outputs(:) !! selected one-based outputs
        type(clti_system) :: out
        out%core = self%core%select_io(inputs,outputs)
    end function clti_select_io

    function clti_minreal(self,tol) result(out)
        class(clti_system), intent(in) :: self !! analog state-space system
        real(dp), optional, intent(in) :: tol !! rank tolerance
        type(clti_system) :: out
        out%core = self%core%minreal(tol)
    end function clti_minreal

    function clti_controllability(self) result(matrix)
        class(clti_system), intent(in) :: self !! continuous system
        complex(dp), allocatable :: matrix(:,:) !! Kalman controllability matrix
        matrix = self%core%controllability()
    end function clti_controllability

    function clti_observability(self) result(matrix)
        class(clti_system), intent(in) :: self !! continuous system
        complex(dp), allocatable :: matrix(:,:) !! Kalman observability matrix
        matrix = self%core%observability()
    end function clti_observability

    subroutine clti_cayley_proxy(self,proxy)
        ! The generalized bilinear map converts a Hurwitz A to a Schur A_d.
        ! With M = I - h*A/2, the exact Gramian-preserving realization is
        ! A_d = M**(-1)*(I + h*A/2), B_d = sqrt(h)*M**(-1)*B,
        ! C_d = sqrt(h)*C*M**(-1). The two Gramians are unchanged.
        ! h adapts to the size of A to avoid an unnecessarily slow
        ! discrete doubling recurrence. The value of D is immaterial here.
        class(clti_system), intent(in) :: self !! analog model
        type(dlti_system), intent(out) :: proxy !! Gramian-equivalent discrete model
        type(clti_system) :: ss
        type(discrete_state_space_result) :: mapped
        real(dp) :: step, scale
        integer :: n

        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        n = size(ss%core%a,1)
        if (n == 0) then
            proxy = dlti_ss(ss%core%a,ss%core%b,ss%core%c,ss%core%d)
            return
        end if
        scale = maxval(abs(ss%core%a))
        if (.not. ieee_is_finite(scale)) return
        if (scale > huge(scale)/4.0_dp) return
        step = 1.0_dp/max(1.0e-100_dp,scale)
        mapped = cont2discrete(ss%core%a,ss%core%b,ss%core%c,ss%core%d, &
            step,method='bilinear')
        if (mapped%status /= cont2discrete_success) return
        proxy = dlti_ss(mapped%a,mapped%b/sqrt(step), &
            mapped%c*sqrt(step),ss%core%d)
    end subroutine clti_cayley_proxy

    function clti_gramians(self,tol,max_iter) result(out)
        ! Solve the continuous Lyapunov equations using exact Cayley
        ! equivalence and the discrete-time doubling solver.
        ! The H2 norm is finite only for strictly proper systems (D=0).
        class(clti_system), intent(in) :: self !! continuous system
        real(dp), optional, intent(in) :: tol !! relative doubling tolerance
        integer, optional, intent(in) :: max_iter !! doubling iteration limit
        type(clti_gramian_result) :: out
        type(clti_system) :: ss
        type(dlti_system) :: proxy
        type(dlti_gramian_result) :: result_d
        complex(dp), allocatable :: energy_mat(:,:)
        real(dp) :: energy, threshold
        integer :: i

        if (self%status() /= clti_success) return
        call clti_cayley_proxy(self,proxy)
        if (proxy%status /= dlti_success) then
            out%status = clti_numerical_failure
            return
        end if
        result_d = proxy%gramians(tol=tol,max_iter=max_iter)
        if (result_d%status /= dlti_success) then
            out%status = result_d%status
            return
        end if
        ss = self%to_ss()
        out%wc = result_d%wc
        out%wo = result_d%wo
        out%iterations = result_d%iterations
        out%h2_norm = ieee_value(0.0_dp,ieee_positive_inf)
        if (all(abs(ss%core%d) == 0.0_dp)) then
            energy_mat = matmul(matmul(ss%core%c,out%wc),transpose(conjg(ss%core%c)))
            energy = 0.0_dp
            do i=1,size(energy_mat,1)
                energy = energy + real(energy_mat(i,i),dp)
            end do
            threshold = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp,sum(abs(energy_mat)))
            if (.not. ieee_is_finite(energy) .or. energy < -threshold) then
                out%status = clti_numerical_failure
                return
            end if
            out%h2_norm = sqrt(max(0.0_dp,energy))
        end if
        out%status = clti_success
    end function clti_gramians

    function clti_balanced_truncate(self,order,tol) result(out)
        ! Use Cayley-equivalent Gramians for a square-root balancing basis.
        ! IMPORTANT: apply the projections to the *continuous* A, B and C,
        ! rather than mapping a truncated discrete A back with Tustin.
        ! This is true continuous-time balanced truncation and satisfies
        ! ||G-Gr||_infinity <= 2*sum(hsv(order+1:)) for stable minimal LTI.
        class(clti_system), intent(in) :: self !! continuous system
        integer, intent(in) :: order !! number of states to retain
        real(dp), optional, intent(in) :: tol !! numerical rank tolerance
        type(clti_balanced_result) :: out
        type(clti_system) :: ss
        type(dlti_system) :: proxy
        type(dlti_balanced_result) :: reduction
        complex(dp), allocatable :: a(:,:),b(:,:),c(:,:)

        if (self%status() /= clti_success) return
        ss = self%to_ss()
        if (ss%status() /= clti_success) return
        if (order < 0 .or. order > size(ss%core%a,1)) return
        call clti_cayley_proxy(self,proxy)
        if (proxy%status /= dlti_success) then
            out%status = clti_numerical_failure
            return
        end if
        reduction = proxy%balanced_truncate(order,tol)
        if (allocated(reduction%hsv)) out%hsv = reduction%hsv
        out%error_bound = reduction%error_bound
        if (reduction%status /= dlti_success) then
            out%status = reduction%status
            return
        end if
        if (.not. allocated(reduction%transform) .or. &
            .not. allocated(reduction%inverse)) then
            out%status = clti_numerical_failure
            return
        end if
        a = matmul(matmul(reduction%inverse,ss%core%a),reduction%transform)
        b = matmul(reduction%inverse,ss%core%b)
        c = matmul(ss%core%c,reduction%transform)
        out%model = clti_ss(a,b,c,ss%core%d)
        out%status = out%model%status()
    end function clti_balanced_truncate
end module scifort_signal_clti
