! SPDX-License-Identifier: MIT
program test_signal_tf2sos
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, tf2sos, zpk_conversion_success, &
        zpk_conversion_invalid_input, filter_sos_result, freqz, sosfreqz, &
        frequency_response
    implicit none

    type(filter_sos_result) :: ss
    type(frequency_response) :: response, expected
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: nanval
    real(dp), parameter :: freqs(5) = [0.05_dp,0.4_dp,0.8_dp,1.4_dp,2.7_dp]
    complex(dp) :: b(2), a(3), x, numer, denom, h
    integer :: j, section

    ss = tf2sos([1.0_dp,2.0_dp],[1.0_dp,2.0_dp,3.0_dp])
    call check(ss%status == zpk_conversion_success, 'digital default status')
    call check(size(ss%sos,1) == 1, 'single digital section')
    call check(maxval(abs(ss%sos(1,:) - cmplx( &
        [1.0_dp,2.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp],0.0_dp,dp))) < tol, &
        'SciPy digital coefficient reference')

    ss = tf2sos([2.0_dp],[1.0_dp,3.0_dp],analog=.true.)
    call check(ss%status == zpk_conversion_success, 'analog first order')
    call check(maxval(abs(ss%sos(1,:) - cmplx( &
        [0.0_dp,0.0_dp,2.0_dp,0.0_dp,1.0_dp,3.0_dp],0.0_dp,dp))) < tol, &
        'SciPy analog first-order reference')

    ss = tf2sos([2.0_dp],[1.0_dp],analog=.true.)
    call check(ss%status == zpk_conversion_success, 'constant analog transfer function')
    call check(maxval(abs(ss%sos(1,:) - cmplx( &
        [0.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,1.0_dp],0.0_dp,dp))) < tol, &
        'constant analog coefficients')
    ss = tf2sos([2.0_dp],[1.0_dp])
    call check(ss%status == zpk_conversion_success, 'constant digital transfer function')
    call check(maxval(abs(ss%sos(1,:) - cmplx( &
        [2.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp],0.0_dp,dp))) < tol, &
        'constant digital coefficients')

    ss = tf2sos([1.0_dp],[1.0_dp,-0.4_dp,0.12_dp])
    call check(ss%status == zpk_conversion_success, 'second-order IIR digital')
    expected = freqz([1.0_dp],[1.0_dp,-0.4_dp,0.12_dp],w=freqs)
    response = sosfreqz(ss%sos,w=freqs)
    call check(maxval(abs(response%response-expected%response)) < tol, &
        'digital frequency response')

    ss = tf2sos([1.0_dp],[1.0_dp,-0.4_dp,0.12_dp],analog=.true.)
    call check(ss%status == zpk_conversion_success,'second-order analog')
    do j = 1, size(freqs)
        x = cmplx(0.0_dp,freqs(j),dp)
        h = cmplx(1.0_dp,0.0_dp,dp)
        do section = 1, size(ss%sos,1)
            numer = (ss%sos(section,1)*x + ss%sos(section,2))*x + ss%sos(section,3)
            denom = (ss%sos(section,4)*x + ss%sos(section,5))*x + ss%sos(section,6)
            h = h * (numer/denom)
        end do
        call check(abs(h - 1.0_dp/((x-0.4_dp)*x+0.12_dp)) < tol, &
            'analog frequency response')
    end do

    ss = tf2sos([1.0_dp, 2.0_dp, 3.0_dp],[1.0_dp])
    call check(ss%status == zpk_conversion_success, 'digital improper transfer function')
    call check(maxval(abs(ss%sos(1,:) - cmplx( &
        [1.0_dp,2.0_dp,3.0_dp,1.0_dp,0.0_dp,0.0_dp],0.0_dp,dp))) < tol, &
        'digital improper polynomial')
    ! Unequal pole/zero degrees previously caused keep_odd to reject a
    ! valid transfer function if artificial origin roots could not match.
    ss = tf2sos([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[1.0_dp],pairing='keep_odd')
    call check(ss%status == zpk_conversion_success,'keep_odd improper digital')
    do j = 1, size(freqs)
        x = exp(cmplx(0.0_dp,-freqs(j),dp))
        h = cmplx(1.0_dp,0.0_dp,dp)
        do section = 1, size(ss%sos,1)
            numer = (ss%sos(section,3)*x+ss%sos(section,2))*x+ss%sos(section,1)
            denom = (ss%sos(section,6)*x+ss%sos(section,5))*x+ss%sos(section,4)
            h = h*(numer/denom)
        end do
        call check(abs(h - (1.0_dp+2.0_dp*x+3.0_dp*x*x+4.0_dp*x**3)) < tol, &
            'keep_odd improper digital response')
    end do
    ss = tf2sos([1.0_dp, 2.0_dp, 3.0_dp],[1.0_dp],analog=.true.)
    call check(ss%status == zpk_conversion_invalid_input,'analog improper rejection')
    ss = tf2sos([1.0_dp],[1.0_dp,-0.3_dp],analog=.true.,pairing='nearest')
    call check(ss%status == zpk_conversion_invalid_input,'analog nonminimal pairing')
    ss = tf2sos([1.0_dp],[1.0_dp,-0.3_dp],pairing='bogus')
    call check(ss%status == zpk_conversion_invalid_input,'invalid pairing')
    ss = tf2sos([1.0_dp],[0.0_dp,0.0_dp])
    call check(ss%status == zpk_conversion_invalid_input,'zero denominator rejection')
    ss = tf2sos([real(dp) :: ],[1.0_dp])
    call check(ss%status == zpk_conversion_invalid_input,'empty numerator rejection')
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    ss = tf2sos([nanval],[1.0_dp])
    call check(ss%status == zpk_conversion_invalid_input,'nonfinite numerator rejection')

    b = [cmplx(1.0_dp,0.6_dp,dp), cmplx(-0.35_dp,0.2_dp,dp)]
    a = [cmplx(1.0_dp,0.0_dp,dp), cmplx(-0.3_dp,0.2_dp,dp), &
         cmplx(0.12_dp,-0.04_dp,dp)]
    ss = tf2sos(b,a)
    call check(ss%status == zpk_conversion_success, 'complex polynomial support')
    do j = 1, size(freqs)
        x = exp(cmplx(0.0_dp,-freqs(j),dp))
        h = cmplx(1.0_dp,0.0_dp,dp)
        do section = 1, size(ss%sos,1)
            numer = (ss%sos(section,3)*x+ss%sos(section,2))*x+ss%sos(section,1)
            denom = (ss%sos(section,6)*x+ss%sos(section,5))*x+ss%sos(section,4)
            h = h*(numer/denom)
        end do
        call check(abs(h-(b(1)+b(2)*x)/(a(1)+a(2)*x+a(3)*x*x)) < tol, &
            'complex digital filter response')
    end do
    print '(a)', 'PASS: tf2sos digital/analog transfer-function conversion'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: label !! assertion label
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_tf2sos
