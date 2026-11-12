! SPDX-License-Identifier: MIT
! Gaussian and multidimensional ndimage regression tests.
program test_ndimage_filters
    use scifort_kinds, only : dp
    use scifort_ndimage, only : gaussian_filter_nd, gaussian_filter1d_nd, &
        gaussian_laplace_nd, gaussian_gradient_magnitude_nd, &
        uniform_filter_nd, minimum_filter_nd, maximum_filter_nd, &
        sobel_nd, laplace_nd, ndimage_status_success, ndimage_status_invalid_input
    implicit none
    real(dp) :: x(5,6,3),a(5,6,3),b(5,6,3),c(5,6,3)
    real(dp) :: one(5,6,3),g(5,6,3),u(5,6,3),lo(5,6,3),hi(5,6,3)
    real(dp) :: x15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: strided(10,12,6)
    complex(dp) :: z(5,6,3),w(5,6,3)
    integer :: i,j,k,status

    do k=1,3
        do j=1,6
            do i=1,5
                x(i,j,k)=sin(0.3_dp*i)+cos(0.7_dp*j)+0.13_dp*k*i
            end do
        end do
    end do
    z=cmplx(x,0.5_dp*x+0.3_dp,dp)
    one=1.75_dp
    call gaussian_filter_nd(one,[0.7_dp,1.3_dp,0.9_dp],a,status=status)
    if (status/=ndimage_status_success) error stop 'gaussian status'
    if (maxval(abs(a-one))>5.e-14_dp) error stop 'gaussian preserves constant'
    call gaussian_filter_nd(x,0.0_dp,a,status=status)
    if (status/=ndimage_status_success .or. any(a/=x)) error stop 'zero sigma copy'
    call gaussian_filter_nd(x,[0.6_dp,1.2_dp,0.8_dp],a,status=status)
    if (status/=ndimage_status_success) error stop 'gaussian vector status'
    call gaussian_filter1d_nd(x,0.6_dp,1,b,status=status)
    if (status/=ndimage_status_success) error stop 'gaussian axis status'
    call gaussian_filter1d_nd(b,1.2_dp,2,c,status=status)
    if (status/=ndimage_status_success) error stop 'gaussian axis status 2'
    call gaussian_filter1d_nd(c,0.8_dp,3,b,status=status)
    if (maxval(abs(a-b))>1.e-14_dp) error stop 'gaussian separability'
    call gaussian_filter_nd(x,[0.6_dp,1.2_dp,0.8_dp], &
        strided(1:10:2,1:12:2,1:6:2),status=status)
    if (status/=ndimage_status_success) error stop 'strided Gaussian status'
    if (maxval(abs(a-strided(1:10:2,1:12:2,1:6:2)))>1.e-14_dp) &
        error stop 'strided Gaussian values'
    call gaussian_filter_nd(z,[0.6_dp,1.2_dp,0.8_dp],w,status=status)
    if (status/=ndimage_status_success) error stop 'complex gaussian status'
    if (maxval(abs(real(w,dp)-a))>2.e-14_dp) error stop 'complex gaussian real'
    if (maxval(abs(aimag(w)-(0.5_dp*a+0.3_dp)))>2.e-14_dp) &
        error stop 'complex gaussian imaginary'

    call gaussian_filter_nd(x,[0.8_dp,1.3_dp,0.9_dp],a, &
                            order=[1,0,2],status=status,mode='mirror')
    if (status/=ndimage_status_success) error stop 'gaussian derivatives'
    if (.not.all(abs(a)<huge(1.0_dp))) error stop 'finite derivatives'
    call gaussian_gradient_magnitude_nd(x,1.0_dp,g,status=status)
    if (status/=ndimage_status_success .or. minval(g)<0.0_dp) &
        error stop 'gradient magnitude'
    call gaussian_laplace_nd(x,[0.8_dp,1.2_dp,0.9_dp],a,status=status)
    if (status/=ndimage_status_success) error stop 'gaussian laplace status'
    call uniform_filter_nd(x,[3,4,2],u,status=status,mode='nearest')
    if (status/=ndimage_status_success) error stop 'uniform status'
    call uniform_filter_nd(z,[3,4,2],w,status=status,mode='nearest')
    if (maxval(abs(real(w,dp)-u))>1.e-14_dp) error stop 'complex uniform'
    call minimum_filter_nd(x,[3,4,2],lo,status=status)
    if (status/=ndimage_status_success) error stop 'minimum status'
    call maximum_filter_nd(x,[3,4,2],hi,status=status)
    if (status/=ndimage_status_success) error stop 'maximum status'
    if (any(lo>hi)) error stop 'minimum greater than maximum'
    call sobel_nd(one,axis=2,y=a,status=status)
    if (status/=ndimage_status_success .or. maxval(abs(a))>1.e-14_dp) &
        error stop 'sobel constant'
    call laplace_nd(one,a,status=status)
    if (status/=ndimage_status_success .or. maxval(abs(a))>1.e-14_dp) &
        error stop 'laplace constant'

    x15=2.5_dp
    call gaussian_filter_nd(x15,0.85_dp,y15,status=status)
    if (status/=ndimage_status_success) error stop 'rank-15 gaussian status'
    if (maxval(abs(x15-y15))>1.e-13_dp) error stop 'rank-15 gaussian values'
    call minimum_filter_nd(x15,3,y15,status=status)
    if (status/=ndimage_status_success .or. maxval(abs(x15-y15))>0.0_dp) &
        error stop 'rank-15 minimum'

    call gaussian_filter_nd(x,[-1.0_dp,1.0_dp,1.0_dp],a,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid sigma'
    call gaussian_filter_nd(x,[1.0_dp,1.0_dp],a,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'sigma length'
    call gaussian_filter_nd(x,1.0_dp,a,order=4,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'unsupported derivative'
    call gaussian_filter1d_nd(x,1.0_dp,4,a,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid axis'
    call gaussian_filter1d_nd(x,1.0_dp,1,a,status=status,radius=100001)
    if (status/=2) error stop 'oversized kernel guard'
    call gaussian_filter_nd(x,0.0_dp,a,order=1,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'zero sigma derivative'
    call uniform_filter_nd(x,0,a,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid uniform width'
    call minimum_filter_nd(x,[3,3],a,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid min sizes'
    call sobel_nd(x,0,a,status=status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid sobel axis'
    print '(a)', 'PASS: ndimage Gaussian, neighborhood, and gradient filters'
end program test_ndimage_filters
