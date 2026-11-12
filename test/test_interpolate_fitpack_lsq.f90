! SPDX-License-Identifier: MIT
! Banded Givens QR, FITPACK-style weighted least squares, and interpolation.
program test_interpolate_fitpack_lsq
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_interpolate, only : lsq_univariate_spline, interpolated_univariate_spline, &
        bspline_interpolator_type, bspline_array_interpolator_type, bspline_design_matrix, &
        bspline_status_singular, bspline_status_invalid_input
    implicit none
    integer, parameter :: n=83
    type(bspline_interpolator_type) :: s, i_spline, bad
    type(bspline_array_interpolator_type) :: sm
    real(dp) :: x(n), y(n), yy(n,3), w(n), t(7), query, &
        residual(n), scale, poor(3), q2(1)
    real(dp), allocatable :: design(:,:), knots(:), fit(:), mat(:,:), big_x(:), big_y(:), big_t(:)
    real(dp), allocatable :: channel(:)
    integer :: j, k, stat

    do j=1,n
        x(j)=real(j-1,dp)*0.1_dp+0.008_dp*sin(real(j,dp)*0.4_dp)
        y(j)=0.9_dp+0.3_dp*x(j)-0.05_dp*x(j)**2+0.015_dp*x(j)**3 &
            +0.04_dp*sin(1.7_dp*x(j))
        w(j)=0.7_dp+0.07_dp*real(mod(j,11),dp)
    end do
    do j=1,size(t)
        t(j)=x(1)+real(j,dp)*(x(n)-x(1))/real(size(t)+1,dp)
    end do
    yy(:,1)=y
    yy(:,2)=2.0_dp*y+0.6_dp
    yy(:,3)=3.0_dp-0.3_dp*x

    do k=1,5
        s=lsq_univariate_spline(x,y,t,w=w,k=k,status=stat)
        if (stat/=0 .or. .not. s%valid()) error stop 'LSQ scalar fit failed'
        sm=lsq_univariate_spline(x,yy,t,w=w,k=k,status=stat)
        if (stat/=0 .or. .not. sm%valid() .or. sm%channels()/=3) error stop 'LSQ matrix fit failed'
        channel=sm%evaluate(2.1_dp)
        if (abs(s%evaluate(2.1_dp)-channel(1))>1.0e-12_dp) &
            error stop 'first channel mismatch'
        do j=1,15
            query=x(1)+real(j-1,dp)*(x(n)-x(1))/14.0_dp
            channel=sm%evaluate(query)
            if (abs(channel(2)-2.0_dp*s%evaluate(query)-0.6_dp)>5.0e-12_dp) &
                error stop 'affine channel transform'
            if (abs(channel(3)-(3.0_dp-0.3_dp*query))>5.0e-12_dp) &
                error stop 'affine reproduction'
        end do
        ! Weighted residual must be orthogonal to every B-spline basis column.
        knots=s%knots()
        design=bspline_design_matrix(x,knots,k,status=stat)
        if (stat/=0) error stop 'cannot evaluate design matrix'
        fit=s%evaluate(x)
        residual=w**2*(fit-y)
        scale=max(1.0_dp,sum(abs(w*y)))
        if (maxval(abs(matmul(residual,design)))>1.0e-10_dp*scale) &
            error stop 'weighted normal equation residual'
        if (.not. ieee_is_finite(s%integrate(x(1),x(n)))) error stop 'integration'
    end do

    i_spline=interpolated_univariate_spline(x,y,k=3,status=stat)
    if (stat/=0 .or. .not. i_spline%valid()) error stop 'exact spline initialization'
    if (maxval(abs(i_spline%evaluate(x)-y))>5.0e-12_dp) &
        error stop 'exact interpolation failed'
    sm=interpolated_univariate_spline(x,yy,k=4,status=stat)
    if (stat/=0 .or. .not. sm%valid()) error stop 'exact matrix interpolation failed'
    mat=sm%evaluate(x)
    if (maxval(abs(mat-yy))>1.0e-11_dp) error stop 'matrix exact interpolation'

    bad=interpolated_univariate_spline(x,y,k=6,status=stat)
    if (stat/=bspline_status_invalid_input .or. bad%valid()) error stop 'interpolated degree rejection'
    bad=lsq_univariate_spline(x,y,t,k=0,status=stat)
    if (stat/=bspline_status_invalid_input .or. bad%valid()) error stop 'degree zero rejection'
    bad=lsq_univariate_spline(x,y,t,w=-w,status=stat)
    if (stat/=bspline_status_invalid_input .or. bad%valid()) error stop 'invalid weights'
    bad=lsq_univariate_spline(x,y,t(7:1:-1),status=stat)
    if (stat/=bspline_status_invalid_input .or. bad%valid()) error stop 'descending knots'
    bad=lsq_univariate_spline(x,y,[x(1)],status=stat)
    if (stat/=bspline_status_invalid_input .or. bad%valid()) error stop 'boundary knots'
    bad=lsq_univariate_spline(x,y,t,bbox=[x(1)+.1_dp,x(n)],status=stat)
    if (stat/=bspline_status_invalid_input .or. bad%valid()) error stop 'invalid bbox'
    poor=[x(1)+0.0001_dp,x(1)+0.0002_dp,x(1)+0.0003_dp]
    bad=lsq_univariate_spline(x,y,poor,status=stat)
    if (stat/=bspline_status_singular .or. bad%valid()) error stop 'rank-deficient knots'
    bad=lsq_univariate_spline(x,y,t,status=stat,extrapolate=.false.)
    if (stat/=0 .or. .not. ieee_is_nan(bad%evaluate(x(1)-1.0_dp))) &
        error stop 'nonextrapolating B-spline'
    q2=[0.0_dp]
    bad=lsq_univariate_spline(x,y,q2,bbox=[x(1),x(n)],k=6,status=stat)
    if (stat/=bspline_status_invalid_input) error stop 'unsupported spline degree'

    ! Large fit must not form the dense 8000-by-500 design matrix.
    allocate(big_x(8000),big_y(8000),big_t(340))
    do j=1,size(big_x)
        big_x(j)=real(j-1,dp)*0.01_dp
        big_y(j)=sin(0.1_dp*big_x(j))+0.03_dp*cos(0.7_dp*big_x(j))
    end do
    do j=1,size(big_t)
        big_t(j)=big_x(1)+real(j,dp)*(big_x(size(big_x))-big_x(1))/real(size(big_t)+1,dp)
    end do
    s=lsq_univariate_spline(big_x,big_y,big_t,status=stat)
    if (stat/=0 .or. .not. s%valid()) error stop 'large least-squares solve failed'
    if (.not. ieee_is_finite(s%evaluate(42.0_dp))) error stop 'large spline evaluation'

    print '(a)', 'FITPACK LSQ / interpolating B-spline native regression PASS'
end program test_interpolate_fitpack_lsq
