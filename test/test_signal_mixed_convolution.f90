! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Mixed real/complex 1D convolution, correlation, FFT and overlap-add tests.
program test_signal_mixed_convolution
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_convolve, only : convolve, correlate, choose_conv_method, &
        convolution_status_success, convolution_status_invalid_input
    use scifort_signal_fftconvolve, only : fftconvolve, oaconvolve
    implicit none

    integer :: i, k, state
    real(dp) :: r(500), short(200), timings(2), nanval
    complex(dp) :: c(200), long(500)
    complex(dp), allocatable :: got(:), expected(:)
    character(len=6) :: choice
    character(len=5), parameter :: modes(3) = ['full ', 'same ', 'valid']
    character(len=6), parameter :: methods(3) = ['direct', 'fft   ', 'auto  ']

    r = [(sin(0.11_dp*real(i,dp))+0.2_dp*cos(0.17_dp*real(i,dp)),i=1,500)]
    short = [(cos(0.07_dp*real(i,dp)),i=1,200)]
    c = cmplx(short, [(0.31_dp*sin(0.09_dp*real(i,dp)),i=1,200)],dp)
    long = cmplx(r, [(0.37_dp*cos(0.03_dp*real(i,dp)),i=1,500)],dp)

    do i = 1,3
        do k = 1,3
            expected = convolve(cmplx(r,0.0_dp,dp),c,trim(modes(i)),trim(methods(k)))
            got = convolve(r,c,trim(modes(i)),trim(methods(k)),status=state)
            call assert_match(got,expected,2.0e-10_dp)
            if (state /= convolution_status_success) error stop 'real-complex convolve status'
            expected = convolve(long,cmplx(short,0.0_dp,dp),trim(modes(i)),trim(methods(k)))
            got = convolve(long,short,trim(modes(i)),trim(methods(k)),status=state)
            call assert_match(got,expected,2.0e-10_dp)
            if (state /= convolution_status_success) error stop 'complex-real convolve status'
            expected = correlate(cmplx(r,0.0_dp,dp),c,trim(modes(i)),trim(methods(k)))
            got = correlate(r,c,trim(modes(i)),trim(methods(k)),status=state)
            call assert_match(got,expected,2.0e-10_dp)
            if (state /= convolution_status_success) error stop 'real-complex correlate status'
            expected = correlate(long,cmplx(short,0.0_dp,dp),trim(modes(i)),trim(methods(k)))
            got = correlate(long,short,trim(modes(i)),trim(methods(k)),status=state)
            call assert_match(got,expected,2.0e-10_dp)
            if (state /= convolution_status_success) error stop 'complex-real correlate status'
        end do
        expected = convolve(cmplx(r,0.0_dp,dp),c,trim(modes(i)),method='direct')
        got = fftconvolve(r,c,trim(modes(i)),status=state,axes=[1])
        call assert_match(got,expected,2.0e-10_dp)
        if (state /= convolution_status_success) error stop 'real-complex FFT status'
        got = oaconvolve(r,c,trim(modes(i)),status=state)
        call assert_match(got,expected,2.0e-10_dp)
        if (state /= convolution_status_success) error stop 'real-complex overlap status'
        expected = convolve(long,cmplx(short,0.0_dp,dp),trim(modes(i)),method='direct')
        got = fftconvolve(long,short,trim(modes(i)),status=state,axes=[1])
        call assert_match(got,expected,2.0e-10_dp)
        if (state /= convolution_status_success) error stop 'complex-real FFT status'
        got = oaconvolve(long,short,trim(modes(i)),status=state)
        call assert_match(got,expected,2.0e-10_dp)
        if (state /= convolution_status_success) error stop 'complex-real overlap status'
    end do

    ! The imaginary part of the second input must be conjugated in correlation.
    got = correlate([1.0_dp,2.0_dp], &
        [cmplx(0.0_dp,1.0_dp,dp),cmplx(0.0_dp,3.0_dp,dp)])
    call assert_match(got, [cmplx(0.0_dp,-3.0_dp,dp), &
        cmplx(0.0_dp,-7.0_dp,dp), cmplx(0.0_dp,-2.0_dp,dp)],0.0_dp)

    ! Noncontiguous assumed-shape inputs and both mixed directions.
    got = convolve(r(1:41:2),c(2:30:3),method='direct')
    expected = convolve(cmplx(r(1:41:2),0.0_dp,dp),c(2:30:3),method='direct')
    call assert_match(got,expected,2.0e-12_dp)
    got = fftconvolve(long(2:40:3),short(1:31:2))
    expected = convolve(long(2:40:3),cmplx(short(1:31:2),0.0_dp,dp),method='direct')
    call assert_match(got,expected,2.0e-12_dp)
    got = oaconvolve(r(2:40:2),c(1:31:2))
    expected = convolve(cmplx(r(2:40:2),0.0_dp,dp),c(1:31:2),method='direct')
    call assert_match(got,expected,2.0e-12_dp)

    choice = choose_conv_method(r,c,measure=.true.,times=timings,status=state)
    if (state /= convolution_status_success) error stop 'choose rc status'
    if (choice /= 'direct' .and. choice /= 'fft') error stop 'choose rc name'
    if (any(timings < 0.0_dp)) error stop 'choose rc timings'
    choice = choose_conv_method(long,short,measure=.true.,times=timings,status=state)
    if (state /= convolution_status_success) error stop 'choose cr status'
    if (choice /= 'direct' .and. choice /= 'fft') error stop 'choose cr name'
    if (any(timings < 0.0_dp)) error stop 'choose cr timings'

    got = convolve(r,c,mode='bad',status=state)
    call rejected(got,state)
    got = correlate(long,short,method='bogus',status=state)
    call rejected(got,state)
    got = fftconvolve(r,c,axes=[2],status=state)
    call rejected(got,state)
    got = fftconvolve(long,short,axes=[1,1],status=state)
    call rejected(got,state)
    got = fftconvolve([real(dp) ::],c,status=state)
    call rejected(got,state)
    got = oaconvolve(long,[real(dp) ::],status=state)
    call rejected(got,state)
    got = oaconvolve(r,c,mode='unsupported',status=state)
    call rejected(got,state)
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    got = fftconvolve([1.0_dp,nanval],c,status=state)
    call rejected(got,state)
    got = oaconvolve([cmplx(nanval,0.0_dp,dp)],r,status=state)
    call rejected(got,state)

    print '(a)', 'test_signal_mixed_convolution passed'
contains
    subroutine assert_match(actual,reference,tolerance)
        complex(dp), intent(in) :: actual(:) !! computed convolution output
        complex(dp), intent(in) :: reference(:) !! expected output
        real(dp), intent(in) :: tolerance !! maximum absolute error
        if (size(actual) /= size(reference)) error stop 'mixed convolution size'
        if (any(abs(actual-reference) > tolerance)) error stop 'mixed convolution values'
    end subroutine assert_match

    subroutine rejected(actual,code)
        complex(dp), intent(in) :: actual(:) !! failed operation output
        integer, intent(in) :: code !! status returned by failed operation
        if (code /= convolution_status_invalid_input) error stop 'mixed reject status'
        if (size(actual) /= 0) error stop 'mixed reject nonempty output'
    end subroutine rejected
end program test_signal_mixed_convolution
