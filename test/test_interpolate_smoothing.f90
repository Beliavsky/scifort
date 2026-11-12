! SPDX-License-Identifier: MIT
! Natural cubic smoothing spline native regression, including GCV and weights.
program test_interpolate_smoothing
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_interpolate, only : make_smoothing_spline, cubic_spline, &
        piecewise_cubic_interpolator, ppoly_type
    implicit none
    integer, parameter :: n=21
    real(dp) :: x(n), y(n), y2(n,3), w(n), lam, a, b, q, &
        selected(3), approx(n), long_x(1500), long_y(1500)
    real(dp), allocatable :: co(:,:,:), vals(:,:), integ(:)
    type(ppoly_type) :: s, s2, sgcv, sbig
    type(piecewise_cubic_interpolator) :: natural
    integer :: i,j,status

    do i=1,n
        x(i)=0.04_dp*real(i-1,dp)+0.002_dp*real(i-1,dp)**1.5_dp
        y(i)=sin(3.0_dp*x(i))+0.15_dp*cos(19.0_dp*x(i))
        w(i)=0.75_dp+real(mod(i,7),dp)/9.0_dp
    end do
    s=make_smoothing_spline(x,y,lam=0.0_dp,status=status)
    natural=cubic_spline(x,y,bc_type='natural',status=j)
    if (status/=0 .or. j/=0 .or. .not. s%valid()) error stop 'zero-penalty construction'
    do i=1,101
        q=x(1)+(x(n)-x(1))*real(i-1,dp)/100.0_dp
        if (abs(s%value(q)-natural%evaluate(q))>5.0e-12_dp) error stop 'lambda=0 value'
        if (abs(s%value(q,nu=1)-natural%derivative(q))>5.0e-11_dp) error stop 'lambda=0 derivative'
    end do
    if (abs(s%integral(x(1),x(n))-natural%integrate(x(1),x(n)))>1.0e-11_dp) &
        error stop 'zero-penalty integration'

    lam=0.04_dp
    s=make_smoothing_spline(x,y,w=w,lam=lam,status=status,selected_lam=a)
    if (status/=0 .or. abs(a-lam)>1.0e-15_dp) error stop 'known lambda'
    co=s%coefficients()
    if (size(co,1)/=4 .or. size(co,2)/=n-1 .or. size(co,3)/=1) error stop 'power coefficients'
    do i=1,n
        approx(i)=s%value(x(i))
    end do
    ! The second derivative is natural at the two ends.
    if (abs(s%value(x(1),nu=2))>1.0e-9_dp .or. &
        abs(s%value(x(n),nu=2))>1.0e-9_dp) error stop 'natural boundary'
    ! The roughness penalty annihilates affine functions. Weighted residuals
    ! are therefore orthogonal to the constant and linear trends.
    if (abs(sum(w*(y-approx)))>5.0e-12_dp) error stop 'weighted constant orthogonality'
    if (abs(sum(w*x*(y-approx)))>5.0e-12_dp) error stop 'weighted linear orthogonality'
    do i=1,n-2
        if (abs(s%value(x(i+1),nu=1)-s%value(x(i+1)+1.0e-10_dp,nu=1))>2.0e-7_dp) &
            error stop 'slope continuity'
    end do

    y2(:,1)=y
    y2(:,2)=4.0_dp*y+1.2_dp
    y2(:,3)=2.0_dp-0.7_dp*x
    s2=make_smoothing_spline(x,y2,lam=lam,w=w,status=status,selected_lam=selected)
    if (status/=0 .or. s2%channels()/=3) error stop 'matrix construction'
    if (any(abs(selected-lam)>1.0e-15_dp)) error stop 'matrix penalty'
    vals=s2%evaluate(x)
    if (maxval(abs(vals(:,1)-approx))>1.0e-12_dp) error stop 'first channel'
    if (maxval(abs(vals(:,2)-4.0_dp*approx-1.2_dp))>3.0e-12_dp) error stop 'scaled channel'
    if (maxval(abs(vals(:,3)-y2(:,3)))>2.0e-12_dp) error stop 'affine reproduction'
    integ=s2%integrate(x(1),x(n))
    if (size(integ)/=3 .or. .not. all(ieee_is_finite(integ))) error stop 'matrix integrals'

    sgcv=make_smoothing_spline(x,y2,w=w,status=status,selected_lam=selected)
    if (status/=0 .or. .not. sgcv%valid()) error stop 'GCV fit'
    if (.not. all(ieee_is_finite(selected)) .or. any(selected<0.0_dp) .or. any(selected>real(n,dp))) &
        error stop 'GCV selected penalty outside [0,n]'
    vals=sgcv%evaluate(x)
    if (maxval(abs(vals(:,3)-y2(:,3)))>5.0e-12_dp) &
        error stop 'GCV affine channel'

    s=make_smoothing_spline(x(:4),y(:4),status=status)
    if (status==0 .or. s%valid()) error stop 'invalid minimum sample count'
    s=make_smoothing_spline(x,y,w=-w,status=status)
    if (status==0 .or. s%valid()) error stop 'invalid negative weights'
    s=make_smoothing_spline(x,y,lam=-1.0_dp,status=status)
    if (status==0 .or. s%valid()) error stop 'invalid negative penalty'
    s=make_smoothing_spline(x(n:1:-1),y,status=status)
    if (status==0 .or. s%valid()) error stop 'invalid knot ordering'
    s=make_smoothing_spline(x,y2(:n-1,:),status=status)
    if (status==0 .or. s%valid()) error stop 'invalid matrix shape'
    s=make_smoothing_spline(x,y2,selected_lam=selected(:2),status=status)
    if (status==0 .or. s%valid()) error stop 'invalid selected-penalty shape'

    ! O(n) banded factorization should work without allocating an n*n matrix.
    do i=1,size(long_x)
        long_x(i)=real(i-1,dp)*0.001_dp
        long_y(i)=sin(5.0_dp*long_x(i))+0.1_dp*cos(23.0_dp*long_x(i))
    end do
    sbig=make_smoothing_spline(long_x,long_y,lam=0.01_dp,status=status)
    if (status/=0 .or. .not. sbig%valid()) error stop 'large banded system'
    if (.not. ieee_is_finite(sbig%value(0.4_dp))) error stop 'large sample evaluation'
    print '(a)', 'smoothing spline native regression PASS'
end program test_interpolate_smoothing
