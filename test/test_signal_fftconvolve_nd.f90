! SPDX-License-Identifier: MIT
! Tests for multidimensional FFT convolution and axis selection.
program test_signal_fftconvolve_nd
    use scifort_signal, only : dp, fftconvolve, convolve2d, &
        convolution_status_success, convolution_status_invalid_input
    implicit none
    real(dp) :: a(4,3), b(2,2), c(3,2,2), d(2,2,2)
    real(dp), allocatable :: y(:,:), want(:,:), z(:,:,:)
    complex(dp), allocatable :: yc(:,:), mixed(:,:,:)
    integer :: i,j,k,code

    do j = 1,3
        do i = 1,4
            a(i,j) = sin(real(2*i+3*j,dp))
        end do
    end do
    b(:,1) = [1.0_dp,3.0_dp]
    b(:,2) = [-2.0_dp,0.5_dp]
    y = fftconvolve(a,b,status=code)
    if (code /= convolution_status_success) error stop 'fft rank2 status'
    want = convolve2d(a,b)
    call check2(y,want)
    y = fftconvolve(a,b,'same')
    want = convolve2d(a,b,'same')
    call check2(y,want)
    y = fftconvolve(a,b,'valid')
    want = convolve2d(a,b,'valid')
    call check2(y,want)
    yc = fftconvolve(a,cmplx(b,b,dp))
    if (maxval(abs(yc - cmplx(convolve2d(a,b),convolve2d(a,b),dp))) > 1.0e-12_dp) &
        error stop 'fft mixed rank2'

    y = fftconvolve(a,b(:,1:1),axes=[1])
    if (any(shape(y) /= [5,3])) error stop 'fft axis1 shape'
    if (maxval(abs(y)) < 1.0e-12_dp) error stop 'fft axis1 zero result'
    y = fftconvolve(a,b,axes=[1,1],status=code)
    if (code /= convolution_status_invalid_input .or. size(y) /= 0) &
        error stop 'invalid axes'
    y = fftconvolve(a,b,axes=[3],status=code)
    if (code /= convolution_status_invalid_input .or. size(y) /= 0) &
        error stop 'out-of-range axis'
    y = fftconvolve(a,b,mode='valid',axes=[1],status=code)
    if (code /= convolution_status_invalid_input .or. size(y) /= 0) &
        error stop 'invalid nonconvolution broadcast dimension'
    y = fftconvolve(a,b,axes=[integer::],status=code)
    if (code /= convolution_status_invalid_input .or. size(y) /= 0) &
        error stop 'empty axes'

    do k = 1,2
        do j = 1,2
            do i = 1,3
                c(i,j,k) = real(i+j+k,dp)
            end do
        end do
    end do
    d = 1.0_dp
    z = fftconvolve(c,d,mode='full',status=code)
    if (code /= convolution_status_success) error stop 'fft rank3 status'
    if (any(shape(z) /= [4,3,3])) error stop 'fft rank3 shape'
    if (abs(z(2,2,2) - sum(c(1:2,1:2,1:2))) > 1.0e-12_dp) &
        error stop 'fft rank3 value'
    mixed = fftconvolve(c,cmplx(d,0.3_dp*d,dp))
    if (maxval(abs(real(mixed,dp)-z)) > 1.0e-12_dp) error stop 'fft rank3 mixed'
    if (maxval(abs(aimag(mixed)-0.3_dp*z)) > 1.0e-12_dp) error stop 'fft rank3 imaginary'
    y = fftconvolve(reshape([real(dp)::],[0,3]),b,status=code)
    if (code /= convolution_status_invalid_input .or. size(y) /= 0) &
        error stop 'empty rank-two array'
    print '(a)', 'test_signal_fftconvolve_nd passed' 
contains
    subroutine check2(x,ref)
        real(dp), intent(in) :: x(:,:) !! computed convolution
        real(dp), intent(in) :: ref(:,:) !! direct reference
        if (any(shape(x) /= shape(ref))) error stop 'rank2 shape mismatch'
        if (maxval(abs(x-ref)) > 1.0e-11_dp) error stop 'rank2 mismatch'
    end subroutine check2
end program test_signal_fftconvolve_nd
