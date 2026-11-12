! SPDX-License-Identifier: MIT
! Regression coverage for cost-aware multidimensional FFT padding.
program test_signal_fft_padding_selection
    use scifort_kinds, only : dp
    use scifort_signal, only : fftconvolve, fftconvolve_nd, &
        oaconvolve_nd, convolve_nd, convolve, convolution_status_success
    implicit none
    real(dp) :: a(47,19), b(21,11), got(67,29), want(67,29)
    real(dp) :: axis_kernel(21,1), axis_got(67,19), axis_want(67,19)
    real(dp) :: a15(47,1,1,1,1,1,1,1,1,1,1,1,1,1,19)
    real(dp) :: b15(21,1,1,1,1,1,1,1,1,1,1,1,1,1,11)
    real(dp) :: y15(67,1,1,1,1,1,1,1,1,1,1,1,1,1,29)
    complex(dp) :: ca(47,19), cb(21,11), cy(67,29), cref(67,29)
    real(dp) :: long_a(901,2), short_b(43,1), oay(943,2), oaref(943,2)
    real(dp), allocatable :: x(:), h(:), y(:), ref(:)
    integer :: i, j, k, code, cases
    real(dp), parameter :: tol = 2.5e-11_dp

    cases = 0
    do j=1,19
        do i=1,47
            a(i,j) = sin(0.031_dp*real(i*j,dp)) + 0.2_dp*cos(real(i-2*j,dp))
        end do
    end do
    do j=1,11
        do i=1,21
            b(i,j) = cos(0.14_dp*real(i+j,dp))
        end do
    end do
    call convolve_nd(a,b,want,method='direct',status=code)
    if (code/=0) error stop 'direct convolution status'
    call fftconvolve_nd(a,b,got,status=code)
    if (code/=convolution_status_success) error stop 'fft convolution status'
    if (maxval(abs(got-want))>tol) error stop 'full convolution mismatch'
    cases = cases+1

    ca = cmplx(a,0.3_dp*reshape([(sin(real(k,dp)),k=1,size(a))],shape(a)),dp)
    cb = cmplx(b,-0.2_dp*b,dp)
    call convolve_nd(ca,cb,cref,method='direct',status=code)
    if (code/=0) error stop 'complex direct status'
    call fftconvolve_nd(ca,cb,cy,status=code)
    if (code/=0 .or. maxval(abs(cy-cref))>tol) error stop 'complex convolution mismatch'
    cases = cases+1

    axis_kernel(:,1) = b(:,1)
    call convolve_nd(a,axis_kernel,axis_want,axes=[1],method='direct',status=code)
    if (code/=0) error stop 'axis direct status'
    call fftconvolve_nd(a,axis_kernel,axis_got,axes=[1],status=code)
    if (code/=0 .or. maxval(abs(axis_got-axis_want))>tol) &
        error stop 'single selected axis mismatch'
    cases = cases+1

    a15 = reshape(a,shape(a15))
    b15 = reshape(b,shape(b15))
    call fftconvolve_nd(a15,b15,y15,axes=[1,15],status=code)
    if (code/=0) error stop 'rank15 fft status'
    if (maxval(abs(reshape(y15,[size(got)])-reshape(want,[size(got)])))>tol) &
        error stop 'rank15 fft results'
    cases = cases+1

    allocate(x(2000),h(81))
    do i=1,2000
        x(i) = sin(0.035_dp*real(i,dp))
    end do
    do i=1,81
        h(i) = exp(-0.019_dp*real(i,dp))
    end do
    ref = convolve(x,h,method='direct')
    do j=1,2
        long_a(:,j) = cos(real(j,dp)*[(0.014_dp*real(i,dp),i=1,901)])
    end do
    short_b(:,1) = h(1:43)
    call fftconvolve_nd(long_a,short_b,oaref,axes=[1],status=code)
    if (code/=0) error stop 'blocked reference status'
    call oaconvolve_nd(long_a,short_b,oay,axes=[1],status=code)
    if (code/=0 .or. maxval(abs(oay-oaref))>tol) error stop 'blocked OA mismatch'
    cases=cases+1

    ! High-precision 1D reference using an explicit frequency-domain path.
    y = fftconvolve(x,h,status=code)
    if (code/=0 .or. maxval(abs(y-ref))>tol) error stop 'rank1 fft mismatch'
    cases=cases+1
    print '(a,i0,a)', 'test_signal_fft_padding_selection passed: ',cases,' cases'
end program test_signal_fft_padding_selection
