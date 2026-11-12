! SPDX-License-Identifier: MIT
! Reusable weighted Lomb-Scargle plans: ranks, axes, atomic writes and recovery.
program test_signal_irregular_plan
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, irregular_status_success, &
        irregular_status_invalid_input, irregular_status_singular, &
        lombscargle_plan, lombscargle_into, lombscargle_amplitude_into, &
        lombscargle_plan_into, lombscargle_plan_amplitude_into
    implicit none
    integer, parameter :: nt=27, nf=13
    real(dp) :: t(nt), freqs(nf), wt(nt), sig(2,nt,3), p(2,nf,3), pref(2,nf,3)
    real(dp) :: high(2,1,1,1,1,nt,1,1,1,1,1,1,1,1,2)
    real(dp) :: highp(2,1,1,1,1,nf,1,1,1,1,1,1,1,1,2)
    real(dp) :: long_data(2*nt,4), long_power(2*nf,4), good(2,nf,3)
    complex(dp) :: amp(2,nf,3), aref(2,nf,3)
    type(lombscargle_plan) :: plan
    real(dp) :: nan_value
    integer :: st, i, j, k, mode

    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    do i=1,nt
        t(i) = 0.23_dp*i+0.007_dp*i*i
        wt(i) = 0.3_dp+0.1_dp*mod(3*i,11)
    end do
    do i=1,nf
        freqs(i) = 0.12_dp+0.23_dp*i
    end do
    do k=1,3
        do j=1,2
            do i=1,nt
                sig(j,i,k) = sin((0.34_dp+0.08_dp*j)*t(i)+0.2_dp*k) &
                    +0.13_dp*cos(1.3_dp*t(i)-0.14_dp*j*k)+0.3_dp*k
            end do
        end do
    end do

    ! Default-uninitialized state cannot modify a caller output.
    p = -99.0_dp
    call lombscargle_plan_into(plan,sig,p,axis=2,status=st)
    if (st /= irregular_status_invalid_input .or. any(p /= -99.0_dp)) &
        error stop 'uninitialized plan modified output'

    do mode=1,2
        call plan%initialize(t,freqs,weights=wt,floating_mean=(mode==1),status=st)
        if (st /= irregular_status_success) error stop 'plan initialization failed'
        if (.not. plan%is_initialized()) error stop 'plan not initialized'
        call lombscargle_plan_into(plan,sig,p,axis=2, &
            precenter=(mode==2),normalize='normalize',status=st)
        if (st /= irregular_status_success) error stop 'planned power failed'
        call lombscargle_into(t,sig,freqs,pref,axis=2,weights=wt, &
            floating_mean=(mode==1),precenter=(mode==2), &
            normalize='normalize',status=st)
        if (st /= irregular_status_success) error stop 'reference power failed'
        call check_real(p,pref,6.0e-13_dp)
        call lombscargle_plan_amplitude_into(plan,sig,amp,axis=-2, &
            precenter=(mode==2),status=st)
        if (st /= irregular_status_success) error stop 'planned amplitude failed'
        call lombscargle_amplitude_into(t,sig,freqs,aref,axis=2,weights=wt, &
            floating_mean=(mode==1),precenter=(mode==2),status=st)
        if (st /= irregular_status_success) error stop 'reference amplitude failed'
        if (maxval(abs(amp-aref)) > 6.0e-13_dp) error stop 'amplitude difference'
        call lombscargle_plan_into(plan,sig,p,axis=2,normalize='power',status=st)
        call lombscargle_into(t,sig,freqs,pref,axis=2,weights=wt, &
            floating_mean=(mode==1),normalize='power',status=st)
        call check_real(p,pref,3.0e-12_dp)
    end do

    ! An invalid reinitialization must preserve the previous valid plan.
    call lombscargle_plan_into(plan,sig,good,axis=2,status=st)
    call plan%initialize(t,freqs,weights=-wt,status=st)
    if (st /= irregular_status_invalid_input) error stop 'bad weights accepted'
    call lombscargle_plan_into(plan,sig,p,axis=2,status=st)
    if (st /= irregular_status_success) error stop 'invalid reinit lost previous plan'
    call check_real(p,good,0.0_dp)

    ! Rank-15, negative axes, and genuinely noncontiguous input/output arrays.
    high = 0.0_dp
    high(:,1,1,1,1,:,1,1,1,1,1,1,1,1,1) = sig(:,:,1)
    high(:,1,1,1,1,:,1,1,1,1,1,1,1,1,2) = sig(:,:,3)
    call lombscargle_plan_into(plan,high,highp,axis=6,status=st)
    if (st /= irregular_status_success) error stop 'rank-15 plan'
    call check_real(highp(:,1,1,1,1,:,1,1,1,1,1,1,1,1,1),good(:,:,1),6.0e-13_dp)
    call check_real(highp(:,1,1,1,1,:,1,1,1,1,1,1,1,1,2),good(:,:,3),6.0e-13_dp)
    long_data = -8.0_dp
    long_power = -7.0_dp
    long_data(1::2,1) = sig(1,:,1)
    long_data(1::2,2) = sig(2,:,1)
    long_data(1::2,3) = sig(1,:,2)
    long_data(1::2,4) = sig(2,:,2)
    call lombscargle_plan_into(plan,long_data(1::2,:), &
        long_power(1::2,:),axis=1,status=st)
    if (st /= irregular_status_success) error stop 'strided input and output'
    if (any(long_power(2::2,:) /= -7.0_dp)) error stop 'strided output modified gaps'
    call check_real(long_power(1::2,1),good(1,:,1),8.0e-13_dp)
    call check_real(long_power(1::2,4),good(2,:,2),8.0e-13_dp)

    ! Errors after partially processing channels never modify caller storage.
    p = 91.0_dp
    sig(2,:,3) = 0.0_dp
    call lombscargle_plan_into(plan,sig,p,axis=2,status=st)
    if (st /= irregular_status_singular) error stop 'singular channel not rejected'
    if (any(p /= 91.0_dp)) error stop 'singular channel overwrote output'
    sig(2,:,3) = nan_value
    amp = cmplx(91.0_dp,0.0_dp,dp)
    call lombscargle_plan_amplitude_into(plan,sig,amp,axis=2,status=st)
    if (st /= irregular_status_invalid_input) error stop 'nonfinite samples accepted'
    if (any(amp /= cmplx(91.0_dp,0.0_dp,dp))) error stop 'failed amplitude changed output'
    p = 91.0_dp
    call lombscargle_plan_into(plan,sig,p,axis=8,status=st)
    if (st /= irregular_status_invalid_input .or. any(p /= 91.0_dp)) &
        error stop 'invalid axis accepted'
    call lombscargle_plan_into(plan,sig,p,axis=2,normalize='bogus',status=st)
    if (st /= irregular_status_invalid_input .or. any(p /= 91.0_dp)) &
        error stop 'invalid normalization accepted'

    call plan%clear()
    if (plan%is_initialized()) error stop 'clear retained cached arrays'
    print '(a)', 'test_signal_irregular_plan: PASS'

contains
    subroutine check_real(actual, expected, tol)
        real(dp), intent(in) :: actual(..) !! measured results of a plan evaluation
        real(dp), intent(in) :: expected(..) !! independent unplanned reference results
        real(dp), intent(in) :: tol !! maximum allowed absolute difference
        real(dp), allocatable :: a(:), b(:)

        select rank (actual)
        rank(1)
            a=reshape(actual,[size(actual)])
        rank(2)
            a=reshape(actual,[size(actual)])
        rank(3)
            a=reshape(actual,[size(actual)])
        rank default
            error stop 'bad test rank'
        end select
        select rank (expected)
        rank(1)
            b=reshape(expected,[size(expected)])
        rank(2)
            b=reshape(expected,[size(expected)])
        rank(3)
            b=reshape(expected,[size(expected)])
        rank default
            error stop 'bad reference rank'
        end select
        if (size(a) /= size(b)) error stop 'test shape mismatch'
        if (maxval(abs(a-b)) > tol) then
            print '(a,es12.4)', 'max error: ', maxval(abs(a-b))
            error stop 'planned versus unplanned mismatch'
        end if
    end subroutine check_real
end program test_signal_irregular_plan
