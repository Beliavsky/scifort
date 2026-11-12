! SPDX-License-Identifier: MIT
! Deterministic tests for FFT and overlap-add convolution.
program test_signal_fftconvolve
    use scifort_signal
    implicit none
    integer :: i,j,code
    real(dp), allocatable :: a(:), b(:), got(:), expected(:)
    complex(dp), allocatable :: z(:), h(:), actual(:), reference(:)
    character(5), parameter :: mode(3) = ['full ', 'same ', 'valid']

    a = [(sin(real(i,dp)*0.17_dp)+cos(real(i,dp)*0.07_dp), i=1,501)]
    b = [(exp(-real(i,dp)*0.04_dp),i=1,19)]
    z = cmplx(a,a*0.4_dp,dp)
    h = cmplx(b,-b*0.2_dp,dp)
    do j = 1,3
        expected = convolve(a,b,trim(mode(j)))
        got = fftconvolve(a,b,trim(mode(j)),status=code)
        if (code/=convolution_status_success) error stop 'fftconvolve real status'
        call check_real(got,expected)
        got = oaconvolve(a,b,trim(mode(j)),status=code)
        if (code/=convolution_status_success) error stop 'oaconvolve real status'
        call check_real(got,expected)
        reference = convolve(z,h,trim(mode(j)))
        actual = fftconvolve(z,h,trim(mode(j)),status=code)
        if (code/=convolution_status_success) error stop 'fftconvolve complex status'
        call check_complex(actual,reference)
        actual = oaconvolve(z,h,trim(mode(j)),status=code)
        if (code/=convolution_status_success) error stop 'oaconvolve complex status'
        call check_complex(actual,reference)
    end do
    do i = 1, 12
        a = [sin(real(i,dp)),cos(real(i,dp))]
        b = [real(i,dp)*0.125_dp]
        call check_real(fftconvolve(a,b),convolve(a,b))
        call check_real(oaconvolve(b,a),convolve(b,a))
    end do
    got = fftconvolve(a,b,axes=[1],status=code)
    if (code /= convolution_status_success) error stop '1D axis selection'
    got = fftconvolve(a,b,axes=[2],status=code)
    if (code /= convolution_status_invalid_input .or. size(got) /= 0) &
        error stop 'invalid 1D axis'
    got = fftconvolve([1.0_dp,2.0_dp],[3.0_dp],mode='unsupported',status=code)
    if (code/=convolution_status_invalid_input .or. size(got)/=0) error stop 'invalid mode'
    got = fftconvolve([real(dp)::],[1.0_dp],status=code)
    if (code/=convolution_status_invalid_input .or. size(got)/=0) error stop 'empty input'
    got = oaconvolve([1.0_dp],[real(dp)::],status=code)
    if (code/=convolution_status_invalid_input .or. size(got)/=0) error stop 'empty kernel'
    print '(a)', 'test_signal_fftconvolve passed'
contains
    subroutine check_real(x,y)
        real(dp), intent(in) :: x(:) !! FFT result
        real(dp), intent(in) :: y(:) !! exact direct result
        if (size(x)/=size(y)) error stop 'real length mismatch'
        if (any(abs(x-y)>1.0e-10_dp)) error stop 'real convolution mismatch'
    end subroutine check_real
    subroutine check_complex(x,y)
        complex(dp), intent(in) :: x(:) !! FFT result
        complex(dp), intent(in) :: y(:) !! direct result
        if (size(x)/=size(y)) error stop 'complex length mismatch'
        if (any(abs(x-y)>1.0e-10_dp)) error stop 'complex convolution mismatch'
    end subroutine check_complex
end program test_signal_fftconvolve
