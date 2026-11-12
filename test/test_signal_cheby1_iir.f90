! SPDX-License-Identifier: MIT
program test_signal_cheby1_iir
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, cheby1, butter_result, frequency_response, sosfreqz, freqz, &
        design_status_success, design_status_invalid_input, decimate, &
        resample_status_success, resample_status_invalid_input
    implicit none
    type(butter_result) :: filter
    real(dp), allocatable :: y(:), x(:), original(:)
    type(frequency_response) :: resp, from_sos
    complex(dp), allocatable :: z(:)
    real(dp), parameter :: tol = 4.0e-10_dp
    integer :: st, i, ord, style
    character(8), parameter :: modes(4) = ['lowpass ', 'highpass', 'bandpass', 'bandstop']

    filter = cheby1(3, 0.5_dp, 0.3_dp)
    if (filter%status /= design_status_success) error stop 'cheby1 lowpass status'
    call check(filter%b, [0.0444150219969267_dp, 0.1332450659907801_dp, &
        0.1332450659907801_dp, 0.0444150219969267_dp], tol)
    call check(filter%a, [1.0_dp,-1.3867339081961572_dp, &
        1.0542310052025194_dp,-0.3121769210309485_dp], tol)

    do style = 1,4
        do ord = 1,8
            if (style <= 2) then
                filter=cheby1(ord,0.5_dp,0.31_dp,btype=trim(modes(style)))
            else
                filter=cheby1(ord,0.5_dp,[0.18_dp,0.61_dp],btype=trim(modes(style)))
            end if
            if (filter%status /= design_status_success) error stop 'cheby1 mode and order'
            if (.not. all(ieee_is_finite(filter%b))) error stop 'nonfinite cheby1 b'
            if (.not. all(ieee_is_finite(filter%a))) error stop 'nonfinite cheby1 a'
            resp = freqz(filter%b,filter%a,w=[0.0_dp,0.08_dp,0.18_dp,0.32_dp,0.41_dp])
            from_sos = sosfreqz(filter%sos,w=[0.0_dp,0.08_dp,0.18_dp,0.32_dp,0.41_dp])
            if (maxval(abs(resp%response-from_sos%response)) > 3.0e-7_dp) &
                error stop 'Chebyshev SOS / TF mismatch'
        end do
    end do
    filter = cheby1(4, 1.0_dp, [0.2_dp,0.65_dp],btype='bandstop')
    call check(filter%b, [0.0793409967004312_dp,-0.1948620992609594_dp, &
        0.496832543289724_dp,-0.6580489907221253_dp,0.846259651924087_dp, &
        -0.6580489907221253_dp,0.496832543289724_dp, &
        -0.1948620992609593_dp,0.0793409967004312_dp],3.0e-8_dp)
    filter = cheby1(3,0.0_dp,0.3_dp)
    if (filter%status /= design_status_invalid_input) error stop 'invalid ripple'
    filter = cheby1(3,-1.0_dp,0.3_dp)
    if (filter%status /= design_status_invalid_input) error stop 'negative ripple'
    filter = cheby1(4,0.5_dp,1.0_dp)
    if (filter%status /= design_status_invalid_input) error stop 'Nyquist ripple cutoff'

    allocate(x(80))
    do i=1,size(x)
        x(i)=0.01_dp*real(i,dp)+sin(0.11_dp*real(i-1,dp))
    end do
    y=decimate(x,3,n=8,ftype='iir',status=st)
    if (st /= resample_status_success) error stop 'decimate IIR status'
    call check(y(:10),[0.0088270582823885_dp,0.3635063016442114_dp, &
        0.6793737400256057_dp,0.93244057218108_dp,1.0938345674184593_dp, &
        1.1513006201747544_dp,1.1025979889876028_dp,0.9533421476391237_dp, &
        0.7268612919634949_dp,0.4465211180851709_dp],tol)
    original=decimate(x,3,n=8,zero_phase=.false.,status=st)
    if (st /= resample_status_success) error stop 'default decimate is IIR'
    call check(original(:10), [3.4799547399084981e-7_dp,6.2454818246571804e-4_dp, &
        2.6471405873772960e-2_dp,2.0944993984184851e-1_dp, &
        5.8518746530703880e-1_dp,8.7217073735498096e-1_dp, &
        1.0420952079683687_dp,1.1557188170249828_dp, &
        1.1268629242345971_dp,1.0146721964619823_dp],tol)
    z=decimate(cmplx(x,0.15_dp*x,dp),3,n=8,ftype='iir',status=st)
    if (st /= resample_status_success) error stop 'complex IIR status'
    if (maxval(abs(z-cmplx(y,0.15_dp*y,dp)))>tol) error stop 'complex IIR values'
    y=decimate(x,3,ftype='unsupported',status=st)
    if (st /= resample_status_invalid_input .or. size(y) /= 0) error stop 'unknown filter type'
    y=decimate(x,3,ftype='iir',n=0,status=st)
    if (st /= resample_status_invalid_input .or. size(y) /= 0) error stop 'invalid order'
    print '(a)', 'test_signal_cheby1_iir passed'
contains
    subroutine check(a,b,atol)
        real(dp), intent(in) :: a(:) !! actual output
        real(dp), intent(in) :: b(:) !! reference output
        real(dp), intent(in) :: atol !! absolute error bound
        if (size(a) /= size(b)) error stop 'reference shape mismatch'
        if (size(a)>0) then
            if (maxval(abs(a-b))>atol) error stop 'reference value mismatch'
        end if
    end subroutine check
end program test_signal_cheby1_iir
