! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Discrete-time LTI simulation and responses in the spirit of scipy.signal.
! Input samples occupy rows, channels occupy columns. The step and impulse
! response stores one channel-specific experiment per third-array index.
module scifort_signal_discrete_simulation
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_sos, only : sosfilt, sos_status_success
    use scifort_signal_state_space, only : tf2ss, zpk2ss, state_space_result, state_space_success
    implicit none
    private

    integer, parameter, public :: discrete_sim_success = 0
    integer, parameter, public :: discrete_sim_invalid_input = 1
    integer, parameter, public :: discrete_sim_numerical_failure = 2

    type, public :: discrete_simulation_result
        real(dp), allocatable :: t(:) !! output times, starting at zero
        complex(dp), allocatable :: y(:,:) !! output samples (nt,nout)
        complex(dp), allocatable :: x(:,:) !! state history (nt,nstates)
        complex(dp), allocatable :: zf(:,:) !! final SOS states (nsections,2), SOS only
        integer :: status = discrete_sim_invalid_input !! result status
    end type discrete_simulation_result

    type, public :: discrete_sos_bank_result
        real(dp), allocatable :: t(:) !! output times (nt), starting at zero
        complex(dp), allocatable :: y(:,:) !! samples (nt,nchannels)
        complex(dp), allocatable :: x(:,:,:) !! pre-sample states (nt,2*nsections,nchannels)
        complex(dp), allocatable :: zf(:,:,:) !! final SOS states (nsections,2,nchannels)
        integer :: status = discrete_sim_invalid_input !! result status
    end type discrete_sos_bank_result

    type, public :: discrete_response_result
        real(dp), allocatable :: t(:) !! output times
        complex(dp), allocatable :: y(:,:,:) !! output samples (nt,nout,nin)
        integer :: status = discrete_sim_invalid_input !! result status
    end type discrete_response_result

    interface dlsim
        module procedure dlsim_ss_real, dlsim_ss_complex
        module procedure dlsim_ss_real_complex_input, dlsim_ss_complex_real_input
        module procedure dlsim_tf_real, dlsim_tf_complex
        module procedure dlsim_tf_real_complex_input, dlsim_tf_complex_real_input
        module procedure dlsim_sos_real, dlsim_sos_complex
        module procedure dlsim_sos_real_complex_input, dlsim_sos_complex_real_input
        module procedure dlsim_sos_bank_shared_rr, dlsim_sos_bank_shared_rc
        module procedure dlsim_sos_bank_shared_cr, dlsim_sos_bank_shared_cc
        module procedure dlsim_sos_bank_independent_rr, dlsim_sos_bank_independent_rc
        module procedure dlsim_sos_bank_independent_cr, dlsim_sos_bank_independent_cc
    end interface dlsim
    interface dstep
        module procedure dstep_ss_real, dstep_ss_complex
        module procedure dstep_tf_real, dstep_tf_complex
        module procedure dstep_sos_real, dstep_sos_complex
    end interface dstep
    interface dimpulse
        module procedure dimpulse_ss_real, dimpulse_ss_complex
        module procedure dimpulse_tf_real, dimpulse_tf_complex
        module procedure dimpulse_sos_real, dimpulse_sos_complex
    end interface dimpulse
    interface dlsim_zpk
        module procedure dlsim_zpk_real, dlsim_zpk_complex
    end interface dlsim_zpk
    interface dstep_zpk
        module procedure dstep_zpk_real, dstep_zpk_complex
    end interface dstep_zpk
    interface dimpulse_zpk
        module procedure dimpulse_zpk_real, dimpulse_zpk_complex
    end interface dimpulse_zpk
    interface dlsim_sos
        module procedure dlsim_sos_real, dlsim_sos_complex
        module procedure dlsim_sos_real_complex_input, dlsim_sos_complex_real_input
        module procedure dlsim_sos_bank_shared_rr, dlsim_sos_bank_shared_rc
        module procedure dlsim_sos_bank_shared_cr, dlsim_sos_bank_shared_cc
        module procedure dlsim_sos_bank_independent_rr, dlsim_sos_bank_independent_rc
        module procedure dlsim_sos_bank_independent_cr, dlsim_sos_bank_independent_cc
    end interface dlsim_sos
    interface dlsim_sos_bank
        module procedure dlsim_sos_bank_shared_rr, dlsim_sos_bank_shared_rc
        module procedure dlsim_sos_bank_shared_cr, dlsim_sos_bank_shared_cc
        module procedure dlsim_sos_bank_independent_rr, dlsim_sos_bank_independent_rc
        module procedure dlsim_sos_bank_independent_cr, dlsim_sos_bank_independent_cc
    end interface dlsim_sos_bank
    interface dstep_sos
        module procedure dstep_sos_real, dstep_sos_complex
    end interface dstep_sos
    interface dimpulse_sos
        module procedure dimpulse_sos_real, dimpulse_sos_complex
    end interface dimpulse_sos
    public :: dlsim, dstep, dimpulse, dlsim_zpk, dstep_zpk, dimpulse_zpk
    public :: dlsim_sos, dlsim_sos_bank, dstep_sos, dimpulse_sos

contains

    pure logical function finite_matrix(v) result(ok)
        complex(dp), intent(in) :: v(:,:) !! values to validate
        ok = all(ieee_is_finite(real(v, dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_matrix

    pure logical function finite_vector(v) result(ok)
        complex(dp), intent(in) :: v(:) !! values to validate
        ok = all(ieee_is_finite(real(v, dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_vector

    function simulate_core(a,b,c,d,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix (n,n)
        complex(dp), intent(in) :: b(:,:) !! input matrix (n,nin)
        complex(dp), intent(in) :: c(:,:) !! output matrix (nout,n)
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix (nout,nin)
        complex(dp), intent(in) :: u(:,:) !! input samples (nin samples in columns)
        real(dp), intent(in) :: dt !! positive sampling period
        real(dp), optional, intent(in) :: t(:) !! input sampling times, strictly increasing
        complex(dp), optional, intent(in) :: x0(:) !! initial state vector
        type(discrete_simulation_result) :: out
        complex(dp), allocatable :: inputs(:,:)
        integer :: ns, ni, no, nin, nt, i, j, pos
        real(dp) :: stop_time, frac, interval

        ns = size(a,1)
        ni = size(b,2)
        no = size(c,1)
        nin = size(u,1)
        if (size(a,2) /= ns .or. size(b,1) /= ns .or. size(c,2) /= ns) return
        if (size(d,1) /= no .or. size(d,2) /= ni .or. size(u,2) /= ni) return
        if (ni < 1 .or. no < 1 .or. nin < 1) return
        if (.not. ieee_is_finite(dt)) return
        if (dt <= 0.0_dp) return
        if (.not. finite_matrix(a) .or. .not. finite_matrix(b) .or. &
            .not. finite_matrix(c) .or. .not. finite_matrix(d) .or. &
            .not. finite_matrix(u)) return
        if (present(x0)) then
            if (size(x0) /= ns) return
            if (.not. finite_vector(x0)) return
        end if

        nt = nin
        if (present(t)) then
            if (size(t) /= nin .or. nin < 2) return
            if (.not. all(ieee_is_finite(t))) return
            if (any(t(2:) <= t(:nin-1))) return
            stop_time = t(nin)
            if (stop_time < 0.0_dp) return
            if (stop_time/dt >= real(huge(nt)-1,dp)) return
            nt = int(floor(stop_time/dt)) + 1
            if (nt < 1) return
        end if
        allocate(out%t(nt),out%y(nt,no),out%x(nt,ns),inputs(nt,ni))
        if (present(t)) then
            ! Match scipy.signal.dlsim: times are linspace(0,t(end),nt),
            ! not rounded integer multiples of dt if the endpoint is off-grid.
            if (nt == 1) then
                out%t(1) = 0.0_dp
            else
                do i = 1, nt
                    out%t(i) = stop_time * (real(i-1,dp)/real(nt-1,dp))
                end do
            end if
            pos = 1
            do i = 1, nt
                do while (pos < nin-1)
                    if (t(pos+1) >= out%t(i)) exit
                    pos = pos+1
                end do
                interval = t(pos+1)-t(pos)
                frac = (out%t(i)-t(pos))/interval
                inputs(i,:) = (1.0_dp-frac)*u(pos,:) + frac*u(pos+1,:)
            end do
        else
            do i = 1, nt
                out%t(i) = real(i-1,dp)*dt
            end do
            inputs = u
        end if
        out%x = (0.0_dp,0.0_dp)
        if (present(x0)) out%x(1,:) = x0
        do i = 1, nt
            out%y(i,:) = matmul(c,out%x(i,:)) + matmul(d,inputs(i,:))
            if (i < nt) out%x(i+1,:) = matmul(a,out%x(i,:)) + matmul(b,inputs(i,:))
        end do
        if (.not. finite_matrix(out%x) .or. .not. finite_matrix(out%y)) then
            out%status = discrete_sim_numerical_failure
            return
        end if
        out%status = discrete_sim_success
    end function simulate_core

    function dlsim_ss_real(a,b,c,d,u,dt,t,x0) result(out)
        real(dp), intent(in) :: a(:,:) !! real state matrix
        real(dp), intent(in) :: b(:,:) !! real input matrix
        real(dp), intent(in) :: c(:,:) !! real output matrix
        real(dp), intent(in) :: d(:,:) !! real feedthrough matrix
        real(dp), intent(in) :: u(:,:) !! real input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sampling times
        real(dp), intent(in), optional :: x0(:) !! real initial state
        type(discrete_simulation_result) :: out
        if (present(x0)) then
            out = simulate_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
                cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),cmplx(u,0.0_dp,dp),dt,t,cmplx(x0,0.0_dp,dp))
        else
            out = simulate_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
                cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),cmplx(u,0.0_dp,dp),dt,t)
        end if
    end function dlsim_ss_real

    function dlsim_ss_complex(a,b,c,d,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        complex(dp), intent(in) :: u(:,:) !! input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sampling times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        out = simulate_core(a,b,c,d,u,dt,t,x0)
    end function dlsim_ss_complex

    function dlsim_ss_real_complex_input(a,b,c,d,u,dt,t,x0) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        complex(dp), intent(in) :: u(:,:) !! complex input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sampling times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        out = simulate_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),u,dt,t,x0)
    end function dlsim_ss_real_complex_input

    function dlsim_ss_complex_real_input(a,b,c,d,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: u(:,:) !! real input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sampling times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_simulation_result) :: out
        out = simulate_core(a,b,c,d,cmplx(u,0.0_dp,dp),dt,t,x0)
    end function dlsim_ss_complex_real_input

    function simulate_tf_core(b,a,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: b(:) !! numerator, descending powers of z
        complex(dp), intent(in) :: a(:) !! denominator, descending powers of z
        complex(dp), intent(in) :: u(:) !! single-channel input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        type(discrete_simulation_result) :: out
        type(state_space_result) :: ss
        ss = tf2ss(b,a)
        if (ss%status /= state_space_success) return
        out = simulate_core(ss%a,ss%b,ss%c,ss%d,reshape(u,[size(u),1]),dt,t,x0)
    end function simulate_tf_core

    function dlsim_tf_real(b,a,u,dt,t,x0) result(out)
        real(dp), intent(in) :: b(:) !! numerator polynomial
        real(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        real(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_simulation_result) :: out
        if (present(x0)) then
            out = simulate_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
                cmplx(u,0.0_dp,dp),dt,t,cmplx(x0,0.0_dp,dp))
        else
            out = simulate_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),cmplx(u,0.0_dp,dp),dt,t)
        end if
    end function dlsim_tf_real

    function dlsim_tf_complex(b,a,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial
        complex(dp), intent(in) :: a(:) !! denominator polynomial
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        complex(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_simulation_result) :: out
        out = simulate_tf_core(b,a,u,dt,t,x0)
    end function dlsim_tf_complex

    function dlsim_tf_real_complex_input(b,a,u,dt,t,x0) result(out)
        real(dp), intent(in) :: b(:) !! real numerator polynomial
        real(dp), intent(in) :: a(:) !! real denominator polynomial
        complex(dp), intent(in) :: u(:) !! complex scalar input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        type(discrete_simulation_result) :: out
        out = simulate_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),u,dt,t,x0)
    end function dlsim_tf_real_complex_input

    function dlsim_tf_complex_real_input(b,a,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial
        complex(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: u(:) !! real scalar input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        type(discrete_simulation_result) :: out
        out = simulate_tf_core(b,a,cmplx(u,0.0_dp,dp),dt,t,x0)
    end function dlsim_tf_complex_real_input

    function dlsim_zpk_real(z,p,k,u,dt,t,x0) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: k !! gain
        real(dp), intent(in) :: u(:) !! input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        real(dp), intent(in), optional :: x0(:) !! canonical initial state
        type(discrete_simulation_result) :: out
        type(state_space_result) :: ss
        ss = zpk2ss(z,p,k)
        if (ss%status /= state_space_success) return
        if (present(x0)) then
            out = simulate_core(ss%a,ss%b,ss%c,ss%d,reshape(cmplx(u,0.0_dp,dp),[size(u),1]), &
                dt,t,cmplx(x0,0.0_dp,dp))
        else
            out = simulate_core(ss%a,ss%b,ss%c,ss%d,reshape(cmplx(u,0.0_dp,dp),[size(u),1]),dt,t)
        end if
    end function dlsim_zpk_real

    function dlsim_zpk_complex(z,p,k,u,dt,t,x0) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: k !! gain
        complex(dp), intent(in) :: u(:) !! input samples
        real(dp), intent(in) :: dt !! sampling period
        real(dp), intent(in), optional :: t(:) !! input sample times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        type(discrete_simulation_result) :: out
        type(state_space_result) :: ss
        ss = zpk2ss(z,p,k)
        if (ss%status /= state_space_success) return
        out = simulate_core(ss%a,ss%b,ss%c,ss%d,reshape(u,[size(u),1]),dt,t,x0)
    end function dlsim_zpk_complex

    function response_core(a,b,c,d,dt,n,t,x0,impulse) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! sample count, default 100
        real(dp), intent(in), optional :: t(:) !! input times overriding n
        complex(dp), intent(in), optional :: x0(:) !! initial state
        logical, intent(in) :: impulse !! true for impulse, false for step
        type(discrete_response_result) :: out
        type(discrete_simulation_result) :: sim
        complex(dp), allocatable :: u(:,:)
        integer :: nsamples, ni, idx
        nsamples = 100
        if (present(n)) nsamples = n
        if (present(t)) nsamples = size(t)
        if (nsamples < 1) return
        ni = size(b,2)
        if (ni < 1) return
        allocate(u(nsamples,ni))
        do idx = 1, ni
            u = (0.0_dp,0.0_dp)
            if (impulse) then
                u(1,idx) = (1.0_dp,0.0_dp)
            else
                u(:,idx) = (1.0_dp,0.0_dp)
            end if
            sim = simulate_core(a,b,c,d,u,dt,t,x0)
            if (sim%status /= discrete_sim_success) then
                out%status = sim%status
                return
            end if
            if (idx == 1) then
                out%t = sim%t
                allocate(out%y(size(sim%t),size(c,1),ni))
            end if
            out%y(:,:,idx) = sim%y
        end do
        out%status = discrete_sim_success
    end function response_core

    function dstep_ss_real(a,b,c,d,dt,n,t,x0) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        real(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_response_result) :: out
        if (present(x0)) then
            out = response_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
                cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),dt,n,t,cmplx(x0,0.0_dp,dp),.false.)
        else
            out = response_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
                cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),dt,n,t,impulse=.false.)
        end if
    end function dstep_ss_real

    function dstep_ss_complex(a,b,c,d,dt,n,t,x0) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_response_result) :: out
        out = response_core(a,b,c,d,dt,n,t,x0,.false.)
    end function dstep_ss_complex

    function dimpulse_ss_real(a,b,c,d,dt,n,t,x0) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        real(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_response_result) :: out
        if (present(x0)) then
            out = response_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
                cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),dt,n,t,cmplx(x0,0.0_dp,dp),.true.)
        else
            out = response_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
                cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),dt,n,t,impulse=.true.)
        end if
    end function dimpulse_ss_real

    function dimpulse_ss_complex(a,b,c,d,dt,n,t,x0) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! initial state
        type(discrete_response_result) :: out
        out = response_core(a,b,c,d,dt,n,t,x0,.true.)
    end function dimpulse_ss_complex

    ! Direct-form II transposed SOS simulation.  Use biquad states rather than
    ! reconstructing transfer-function coefficients, which is ill-conditioned
    ! for high-order filters.  x(i,:) contains the state *before* sample i:
    ! [section1 delay1, section1 delay2, section2 delay1, ...].
    ! The final states are also returned in zf(nsections,2), suitable as zi
    ! for the next call (including a streaming continuation).
    function sos_simulate_core(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! normalized sections (nsections,6), a0=1
        complex(dp), intent(in) :: u(:) !! scalar input samples
        real(dp), intent(in) :: dt !! positive simulation step
        real(dp), optional, intent(in) :: t(:) !! strictly increasing input timestamps
        complex(dp), optional, intent(in) :: zi(:,:) !! initial delays (nsections,2)
        type(discrete_simulation_result) :: out
        complex(dp), allocatable :: inputs(:,:), state(:,:)
        complex(dp) :: sample, next_sample
        real(dp) :: stop_time, fraction, interval
        integer :: nin, nt, nsec, i, k, pos

        nin = size(u)
        nsec = size(sos,1)
        if (nin < 1 .or. nsec < 1 .or. size(sos,2) /= 6) return
        if (.not. ieee_is_finite(dt)) return
        if (dt <= 0.0_dp) return
        if (.not. finite_matrix(sos) .or. .not. finite_vector(u)) return
        if (any(sos(:,4) /= cmplx(1.0_dp,0.0_dp,dp))) return
        if (present(zi)) then
            if (any(shape(zi) /= [nsec,2])) return
            if (.not. finite_matrix(zi)) return
        end if
        nt = nin
        if (present(t)) then
            if (nin < 2 .or. size(t) /= nin) return
            if (.not. all(ieee_is_finite(t))) return
            if (any(t(2:) <= t(:nin-1))) return
            stop_time = t(nin)
            if (stop_time < 0.0_dp) return
            if (stop_time/dt >= real(huge(nt)-1,dp)) return
            nt = int(floor(stop_time/dt)) + 1
            if (nt < 1) return
        end if
        allocate(out%t(nt),out%y(nt,1),out%x(nt,2*nsec), &
            out%zf(nsec,2),inputs(nt,1),state(nsec,2))
        if (present(t)) then
            if (nt == 1) then
                out%t(1) = 0.0_dp
            else
                do i = 1, nt
                    out%t(i) = stop_time*(real(i-1,dp)/real(nt-1,dp))
                end do
            end if
            pos = 1
            do i = 1, nt
                do while (pos < nin-1)
                    if (t(pos+1) >= out%t(i)) exit
                    pos = pos+1
                end do
                interval = t(pos+1)-t(pos)
                fraction = (out%t(i)-t(pos))/interval
                inputs(i,1) = (1.0_dp-fraction)*u(pos)+fraction*u(pos+1)
            end do
        else
            do i = 1, nt
                out%t(i) = real(i-1,dp)*dt
            end do
            inputs(:,1) = u
        end if
        if (.not. all(ieee_is_finite(out%t))) then
            out%status = discrete_sim_numerical_failure
            return
        end if
        state = (0.0_dp,0.0_dp)
        if (present(zi)) state = zi
        do i = 1, nt
            sample = inputs(i,1)
            do k = 1, nsec
                out%x(i,2*k-1) = state(k,1)
                out%x(i,2*k) = state(k,2)
                next_sample = sos(k,1)*sample + state(k,1)
                state(k,1) = sos(k,2)*sample - sos(k,5)*next_sample + state(k,2)
                state(k,2) = sos(k,3)*sample - sos(k,6)*next_sample
                sample = next_sample
            end do
            out%y(i,1) = sample
        end do
        out%zf = state
        if (.not. finite_matrix(out%y) .or. .not. finite_matrix(out%x) .or. &
            .not. finite_matrix(out%zf)) then
            out%status = discrete_sim_numerical_failure
            return
        end if
        out%status = discrete_sim_success
    end function sos_simulate_core

    function dlsim_sos_real(sos, u, dt, t, zi) result(out)
        real(dp), intent(in) :: sos(:,:) !! normalized digital SOS (nsections,6)
        real(dp), intent(in) :: u(:) !! real samples at t or at multiples of dt
        real(dp), intent(in) :: dt !! positive sampling period
        real(dp), optional, intent(in) :: t(:) !! optional input sample times
        real(dp), optional, intent(in) :: zi(:,:) !! initial SOS delays (nsections,2)
        type(discrete_simulation_result) :: out
        if (present(zi)) then
            out = sos_simulate_core(cmplx(sos,0.0_dp,dp),cmplx(u,0.0_dp,dp), &
                dt,t,cmplx(zi,0.0_dp,dp))
        else
            out = sos_simulate_core(cmplx(sos,0.0_dp,dp),cmplx(u,0.0_dp,dp),dt,t)
        end if
    end function dlsim_sos_real

    function dlsim_sos_complex(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! normalized complex SOS (nsections,6)
        complex(dp), intent(in) :: u(:) !! complex input samples
        real(dp), intent(in) :: dt !! positive sampling period
        real(dp), optional, intent(in) :: t(:) !! optional input sample times
        complex(dp), optional, intent(in) :: zi(:,:) !! initial SOS delays (nsections,2)
        type(discrete_simulation_result) :: out
        out = sos_simulate_core(sos,u,dt,t,zi)
    end function dlsim_sos_complex

    function dlsim_sos_real_complex_input(sos, u, dt, t, zi) result(out)
        real(dp), intent(in) :: sos(:,:) !! normalized real SOS (nsections,6)
        complex(dp), intent(in) :: u(:) !! complex input samples
        real(dp), intent(in) :: dt !! positive sampling period
        real(dp), optional, intent(in) :: t(:) !! optional input sample times
        complex(dp), optional, intent(in) :: zi(:,:) !! initial SOS delays (nsections,2)
        type(discrete_simulation_result) :: out
        out = sos_simulate_core(cmplx(sos,0.0_dp,dp),u,dt,t,zi)
    end function dlsim_sos_real_complex_input

    function dlsim_sos_complex_real_input(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! normalized complex SOS (nsections,6)
        real(dp), intent(in) :: u(:) !! real input samples
        real(dp), intent(in) :: dt !! positive sampling period
        real(dp), optional, intent(in) :: t(:) !! optional input sample times
        complex(dp), optional, intent(in) :: zi(:,:) !! initial SOS delays (nsections,2)
        type(discrete_simulation_result) :: out
        out = sos_simulate_core(sos,cmplx(u,0.0_dp,dp),dt,t,zi)
    end function dlsim_sos_complex_real_input

    ! Simulate independent input channels without expanding their SOS
    ! cascades into transfer-function polynomials. All channels have the
    ! same number of SOS sections, though their coefficients can differ.
    ! The pre-sample state x(i,:,channel) and the final state zf(:,:,channel)
    ! preserve the same ordering as the SISO dlsim_sos routine.
    function sos_bank_simulate_core(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:,:) !! (nsections,6,nchannels), unit a0
        complex(dp), intent(in) :: u(:,:) !! input samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample interval
        real(dp), optional, intent(in) :: t(:) !! input times, strictly increasing
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        complex(dp), allocatable :: state(:,:,:), samples(:)
        complex(dp) :: value, next_value
        real(dp) :: stop_time, fraction, interval
        integer :: nin, nt, nsec, nchannels, i, channel, section, pos

        nin = size(u,1)
        nchannels = size(u,2)
        nsec = size(sos,1)
        if (nin < 1 .or. nchannels < 1 .or. nsec < 1) return
        if (size(sos,2) /= 6 .or. size(sos,3) /= nchannels) return
        if (.not. ieee_is_finite(dt)) return
        if (dt <= 0.0_dp) return
        if (.not. finite_matrix(u)) return
        if (.not. all(ieee_is_finite(real(sos,dp))) .or. &
            .not. all(ieee_is_finite(aimag(sos)))) return
        if (any(sos(:,4,:) /= cmplx(1.0_dp,0.0_dp,dp))) return
        if (present(zi)) then
            if (any(shape(zi) /= [nsec,2,nchannels])) return
            if (.not. all(ieee_is_finite(real(zi,dp))) .or. &
                .not. all(ieee_is_finite(aimag(zi)))) return
        end if

        nt = nin
        if (present(t)) then
            if (nin < 2 .or. size(t) /= nin) return
            if (.not. all(ieee_is_finite(t))) return
            if (any(t(2:) <= t(:nin-1))) return
            stop_time = t(nin)
            if (stop_time < 0.0_dp) return
            if (stop_time/dt >= real(huge(nt)-1,dp)) return
            nt = int(floor(stop_time/dt)) + 1
            if (nt < 1) return
        end if

        allocate(out%t(nt), out%y(nt,nchannels), &
            out%x(nt,2*nsec,nchannels), out%zf(nsec,2,nchannels), &
            state(nsec,2,nchannels), samples(nchannels))
        if (present(t)) then
            if (nt == 1) then
                out%t(1) = 0.0_dp
            else
                do i = 1, nt
                    out%t(i) = stop_time*(real(i-1,dp)/real(nt-1,dp))
                end do
            end if
        else
            do i = 1, nt
                out%t(i) = real(i-1,dp)*dt
            end do
        end if
        if (.not. all(ieee_is_finite(out%t))) then
            out%status = discrete_sim_numerical_failure
            return
        end if

        state = (0.0_dp,0.0_dp)
        if (present(zi)) state = zi
        pos = 1
        do i = 1, nt
            if (present(t)) then
                do while (pos < nin-1)
                    if (t(pos+1) >= out%t(i)) exit
                    pos = pos+1
                end do
                interval = t(pos+1)-t(pos)
                fraction = (out%t(i)-t(pos))/interval
                samples = (1.0_dp-fraction)*u(pos,:) + fraction*u(pos+1,:)
            else
                samples = u(i,:)
            end if
            do channel = 1, nchannels
                value = samples(channel)
                do section = 1, nsec
                    out%x(i,2*section-1,channel) = state(section,1,channel)
                    out%x(i,2*section,channel) = state(section,2,channel)
                    next_value = sos(section,1,channel)*value + state(section,1,channel)
                    state(section,1,channel) = sos(section,2,channel)*value &
                        - sos(section,5,channel)*next_value + state(section,2,channel)
                    state(section,2,channel) = sos(section,3,channel)*value &
                        - sos(section,6,channel)*next_value
                    value = next_value
                end do
                out%y(i,channel) = value
            end do
        end do
        out%zf = state
        if (.not. finite_matrix(out%y) .or. &
            .not. all(ieee_is_finite(real(out%x,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%x))) .or. &
            .not. all(ieee_is_finite(real(out%zf,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%zf)))) then
            out%status = discrete_sim_numerical_failure
            return
        end if
        out%status = discrete_sim_success
    end function sos_bank_simulate_core

    function dlsim_sos_bank_shared_rr(sos, u, dt, t, zi) result(out)
        real(dp), intent(in) :: sos(:,:) !! normalized SOS coefficients
        real(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        real(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        if (present(zi)) then
            out = sos_bank_simulate_core(spread(cmplx(sos,0.0_dp,dp),dim=3, &
                ncopies=size(u,2)),cmplx(u,0.0_dp,dp),dt,t,cmplx(zi,0.0_dp,dp))
        else
            out = sos_bank_simulate_core(spread(cmplx(sos,0.0_dp,dp),dim=3,ncopies=size(u,2)),cmplx(u,0.0_dp,dp),dt,t)
        end if
    end function dlsim_sos_bank_shared_rr

    function dlsim_sos_bank_shared_rc(sos, u, dt, t, zi) result(out)
        real(dp), intent(in) :: sos(:,:) !! normalized SOS coefficients
        complex(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        out = sos_bank_simulate_core(spread(cmplx(sos,0.0_dp,dp),dim=3,ncopies=size(u,2)),u,dt,t,zi)
    end function dlsim_sos_bank_shared_rc

    function dlsim_sos_bank_shared_cr(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! normalized SOS coefficients
        real(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        out = sos_bank_simulate_core(spread(sos,dim=3,ncopies=size(u,2)),cmplx(u,0.0_dp,dp),dt,t,zi)
    end function dlsim_sos_bank_shared_cr

    function dlsim_sos_bank_shared_cc(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! normalized SOS coefficients
        complex(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        out = sos_bank_simulate_core(spread(sos,dim=3,ncopies=size(u,2)),u,dt,t,zi)
    end function dlsim_sos_bank_shared_cc

    function dlsim_sos_bank_independent_rr(sos, u, dt, t, zi) result(out)
        real(dp), intent(in) :: sos(:,:,:) !! normalized SOS coefficients
        real(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        real(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        if (present(zi)) then
            out = sos_bank_simulate_core(cmplx(sos,0.0_dp,dp),cmplx(u,0.0_dp,dp),dt,t,cmplx(zi,0.0_dp,dp))
        else
            out = sos_bank_simulate_core(cmplx(sos,0.0_dp,dp),cmplx(u,0.0_dp,dp),dt,t)
        end if
    end function dlsim_sos_bank_independent_rr

    function dlsim_sos_bank_independent_rc(sos, u, dt, t, zi) result(out)
        real(dp), intent(in) :: sos(:,:,:) !! normalized SOS coefficients
        complex(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        out = sos_bank_simulate_core(cmplx(sos,0.0_dp,dp),u,dt,t,zi)
    end function dlsim_sos_bank_independent_rc

    function dlsim_sos_bank_independent_cr(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:,:) !! normalized SOS coefficients
        real(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        out = sos_bank_simulate_core(sos,cmplx(u,0.0_dp,dp),dt,t,zi)
    end function dlsim_sos_bank_independent_cr

    function dlsim_sos_bank_independent_cc(sos, u, dt, t, zi) result(out)
        complex(dp), intent(in) :: sos(:,:,:) !! normalized SOS coefficients
        complex(dp), intent(in) :: u(:,:) !! samples (nsamples,nchannels)
        real(dp), intent(in) :: dt !! strictly positive sample spacing
        real(dp), optional, intent(in) :: t(:) !! optional input timestamps
        complex(dp), optional, intent(in) :: zi(:,:,:) !! initial states (nsections,2,nchannels)
        type(discrete_sos_bank_result) :: out
        out = sos_bank_simulate_core(sos,u,dt,t,zi)
    end function dlsim_sos_bank_independent_cc

    ! SOS impulse/step responses directly traverse the individual biquad
    ! state recurrences.  Expanding a long cascade to TF coefficients is
    ! unnecessarily ill-conditioned and can overflow even for stable systems.
    function sos_response_real(sos, dt, n, zi, impulse) result(out)
        real(dp), intent(in) :: sos(:,:) !! digital SOS (sections,6), with a0=1
        real(dp), intent(in) :: dt !! positive time interval between samples
        integer, intent(in), optional :: n !! number of samples; default 100
        real(dp), intent(in), optional :: zi(:,:) !! initial DF-II states (sections,2)
        logical, intent(in) :: impulse !! true for impulse, false for unit step
        type(discrete_response_result) :: out
        real(dp), allocatable :: excitation(:), values(:)
        integer :: nsamples, state, i

        nsamples = 100
        if (present(n)) nsamples = n
        if (nsamples < 1 .or. .not. ieee_is_finite(dt)) return
        if (dt <= 0.0_dp) return
        allocate(excitation(nsamples))
        excitation = 1.0_dp
        if (impulse) then
            excitation = 0.0_dp
            excitation(1) = 1.0_dp
        end if
        values = sosfilt(sos, excitation, zi=zi, status=state)
        if (state /= sos_status_success) return
        if (.not. all(ieee_is_finite(values))) then
            out%status = discrete_sim_numerical_failure
            return
        end if
        allocate(out%t(nsamples), out%y(nsamples,1,1))
        do i = 1, nsamples
            out%t(i) = dt*real(i-1,dp)
        end do
        if (.not. all(ieee_is_finite(out%t))) then
            deallocate(out%t,out%y)
            return
        end if
        out%y(:,1,1) = cmplx(values,0.0_dp,dp)
        out%status = discrete_sim_success
    end function sos_response_real

    function sos_response_complex(sos, dt, n, zi, impulse) result(out)
        complex(dp), intent(in) :: sos(:,:) !! complex digital SOS, with a0=1
        real(dp), intent(in) :: dt !! positive time interval between samples
        integer, intent(in), optional :: n !! number of samples; default 100
        complex(dp), intent(in), optional :: zi(:,:) !! initial DF-II states (sections,2)
        logical, intent(in) :: impulse !! true for impulse, false for unit step
        type(discrete_response_result) :: out
        complex(dp), allocatable :: excitation(:), values(:)
        integer :: nsamples, state, i

        nsamples = 100
        if (present(n)) nsamples = n
        if (nsamples < 1 .or. .not. ieee_is_finite(dt)) return
        if (dt <= 0.0_dp) return
        allocate(excitation(nsamples))
        excitation = (1.0_dp,0.0_dp)
        if (impulse) then
            excitation = (0.0_dp,0.0_dp)
            excitation(1) = (1.0_dp,0.0_dp)
        end if
        values = sosfilt(sos, excitation, zi=zi, status=state)
        if (state /= sos_status_success) return
        if (.not. all(ieee_is_finite(real(values,dp))) .or. &
            .not. all(ieee_is_finite(aimag(values)))) then
            out%status = discrete_sim_numerical_failure
            return
        end if
        allocate(out%t(nsamples),out%y(nsamples,1,1))
        do i = 1, nsamples
            out%t(i) = dt*real(i-1,dp)
        end do
        if (.not. all(ieee_is_finite(out%t))) then
            deallocate(out%t,out%y)
            return
        end if
        out%y(:,1,1) = values
        out%status = discrete_sim_success
    end function sos_response_complex

    function dstep_sos_real(sos, dt, n, zi) result(out)
        real(dp), intent(in) :: sos(:,:) !! digital SOS coefficients (sections,6)
        real(dp), intent(in) :: dt !! positive sample spacing
        integer, intent(in), optional :: n !! number of response samples
        real(dp), intent(in), optional :: zi(:,:) !! initial DF-II states
        type(discrete_response_result) :: out
        out = sos_response_real(sos,dt,n,zi,.false.)
    end function dstep_sos_real

    function dstep_sos_complex(sos, dt, n, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! digital SOS coefficients (sections,6)
        real(dp), intent(in) :: dt !! positive sample spacing
        integer, intent(in), optional :: n !! number of response samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial DF-II states
        type(discrete_response_result) :: out
        out = sos_response_complex(sos,dt,n,zi,.false.)
    end function dstep_sos_complex

    function dimpulse_sos_real(sos, dt, n, zi) result(out)
        real(dp), intent(in) :: sos(:,:) !! digital SOS coefficients (sections,6)
        real(dp), intent(in) :: dt !! positive sample spacing
        integer, intent(in), optional :: n !! number of response samples
        real(dp), intent(in), optional :: zi(:,:) !! initial DF-II states
        type(discrete_response_result) :: out
        out = sos_response_real(sos,dt,n,zi,.true.)
    end function dimpulse_sos_real

    function dimpulse_sos_complex(sos, dt, n, zi) result(out)
        complex(dp), intent(in) :: sos(:,:) !! digital SOS coefficients (sections,6)
        real(dp), intent(in) :: dt !! positive sample spacing
        integer, intent(in), optional :: n !! number of response samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial DF-II states
        type(discrete_response_result) :: out
        out = sos_response_complex(sos,dt,n,zi,.true.)
    end function dimpulse_sos_complex

    function response_tf_core(b,a,dt,n,t,x0,impulse) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial
        complex(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in) :: impulse !! impulse or step response selector
        type(discrete_response_result) :: out
        type(state_space_result) :: ss
        ss = tf2ss(b,a)
        if (ss%status /= state_space_success) return
        out = response_core(ss%a,ss%b,ss%c,ss%d,dt,n,t,x0,impulse)
    end function response_tf_core

    function dstep_tf_real(b,a,dt,n,t,x0) result(out)
        real(dp), intent(in) :: b(:) !! numerator polynomial
        real(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        real(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        if (present(x0)) then
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),dt,n,t, &
                cmplx(x0,0.0_dp,dp),.false.)
        else
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),dt,n,t,impulse=.false.)
        end if
    end function dstep_tf_real

    function dstep_tf_complex(b,a,dt,n,t,x0) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial
        complex(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        out = response_tf_core(b,a,dt,n,t,x0,.false.)
    end function dstep_tf_complex

    function dimpulse_tf_real(b,a,dt,n,t,x0) result(out)
        real(dp), intent(in) :: b(:) !! numerator polynomial
        real(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        real(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        if (present(x0)) then
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),dt,n,t, &
                cmplx(x0,0.0_dp,dp),.true.)
        else
            out = response_tf_core(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),dt,n,t,impulse=.true.)
        end if
    end function dimpulse_tf_real

    function dimpulse_tf_complex(b,a,dt,n,t,x0) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial
        complex(dp), intent(in) :: a(:) !! denominator polynomial
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        out = response_tf_core(b,a,dt,n,t,x0,.true.)
    end function dimpulse_tf_complex

    function response_zpk_core(z,p,k,dt,n,t,x0,impulse) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: k !! gain
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! canonical initial state
        logical, intent(in) :: impulse !! impulse or step response selector
        type(discrete_response_result) :: out
        type(state_space_result) :: ss
        ss = zpk2ss(z,p,k)
        if (ss%status /= state_space_success) return
        out = response_core(ss%a,ss%b,ss%c,ss%d,dt,n,t,x0,impulse)
    end function response_zpk_core

    function dstep_zpk_real(z,p,k,dt,n,t,x0) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: k !! gain
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        real(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        if (present(x0)) then
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp), &
                dt,n,t,cmplx(x0,0.0_dp,dp),.false.)
        else
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp), &
                dt,n,t,impulse=.false.)
        end if
    end function dstep_zpk_real

    function dstep_zpk_complex(z,p,k,dt,n,t,x0) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: k !! gain
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        out = response_zpk_core(z,p,k,dt,n,t,x0,.false.)
    end function dstep_zpk_complex

    function dimpulse_zpk_real(z,p,k,dt,n,t,x0) result(out)
        real(dp), intent(in) :: z(:) !! zeros
        real(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: k !! gain
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        real(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        if (present(x0)) then
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp), &
                dt,n,t,cmplx(x0,0.0_dp,dp),.true.)
        else
            out = response_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp), &
                dt,n,t,impulse=.true.)
        end if
    end function dimpulse_zpk_real

    function dimpulse_zpk_complex(z,p,k,dt,n,t,x0) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: k !! gain
        real(dp), intent(in) :: dt !! sampling period
        integer, intent(in), optional :: n !! number of samples
        real(dp), intent(in), optional :: t(:) !! input times
        complex(dp), intent(in), optional :: x0(:) !! initial canonical state
        type(discrete_response_result) :: out
        out = response_zpk_core(z,p,k,dt,n,t,x0,.true.)
    end function dimpulse_zpk_complex

end module scifort_signal_discrete_simulation
