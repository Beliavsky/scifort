! SPDX-License-Identifier: MIT
! Native Fortran regressions for PPoly/BPoly and cubic-spline conversion.
program test_interpolate_piecewise
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_interpolate, only : ppoly_type, bpoly_type, ppoly, bpoly, &
        ppoly_from_bernstein, bpoly_from_power, bpoly_from_derivatives, &
        ppoly_from_cubic, cubic_spline, piecewise_cubic_interpolator, &
        piecewise_cubic_array_interpolator
    implicit none
    type(ppoly_type) :: p, pp, pa
    type(bpoly_type) :: b, ba
    type(piecewise_cubic_interpolator) :: cubic
    type(piecewise_cubic_array_interpolator) :: cubic_channels
    real(dp) :: x(4), y(4), y2(4,2), q, tol, derivative_samples(3,3)
    real(dp) :: query(5), yh(3,3,2)
    real(dp), allocatable :: val(:,:)
    integer :: status, i
    tol=1.0e-11_dp

    ! Global polynomial 3*x+2 on two intervals, stored in local powers.
    p=ppoly(reshape([3.0_dp,2.0_dp,3.0_dp,5.0_dp],[2,2]),[0.0_dp,1.0_dp,2.0_dp],status=status)
    if (status/=0 .or. .not. p%valid()) error stop 'power constructor'
    if (abs(p%value(1.25_dp)-5.75_dp)>tol) error stop 'power evaluation'
    if (abs(p%value(1.25_dp,nu=1)-3.0_dp)>tol) error stop 'power derivative'
    if (abs(p%integral(0.0_dp,2.0_dp)-10.0_dp)>tol) error stop 'power integral'
    if (abs(p%integral(2.0_dp,0.0_dp)+10.0_dp)>tol) error stop 'reverse integration'
    pp=p%antiderivative(2)
    pa=pp%derivative(2)
    if (abs(pa%value(1.5_dp)-p%value(1.5_dp))>tol) error stop 'power calculus'

    b=bpoly_from_power(p)
    pp=ppoly_from_bernstein(b)
    do i=1,61
        q=-1.0_dp+real(i,dp)/20.0_dp
        if (abs(pp%value(q)-p%value(q))>tol) error stop 'power/Bernstein conversion'
        if (abs(b%value(q)-p%value(q))>tol) error stop 'Bernstein basis evaluation'
        if (abs(b%value(q,nu=1)-3.0_dp)>tol) error stop 'Bernstein derivative'
    end do
    if (abs(b%integral(0.0_dp,2.0_dp)-10.0_dp)>tol) error stop 'Bernstein integration'
    ba=b%antiderivative(2)
    ba=ba%derivative(2)
    if (abs(ba%value(0.5_dp)-b%value(0.5_dp))>tol) error stop 'Bernstein calculus'

    ! Right-hand interpolation at knots, including descending coordinates.
    p=ppoly(reshape([11.0_dp,22.0_dp,33.0_dp],[1,3]),[3.0_dp,2.0_dp,1.0_dp,0.0_dp])
    if (p%value(2.0_dp)/=22.0_dp) error stop 'descending right-hand selection'
    if (p%value(1.0_dp)/=33.0_dp) error stop 'descending second knot'
    p=ppoly(reshape([11.0_dp,22.0_dp,33.0_dp],[1,3]),[3.0_dp,2.0_dp,1.0_dp,0.0_dp],periodic=.true.)
    if (p%value(3.0_dp)/=11.0_dp .or. p%value(0.0_dp)/=11.0_dp) &
        error stop 'descending periodic endpoint'
    if (abs(p%integral(0.0_dp,6.0_dp)-132.0_dp)>tol) error stop 'descending periodic integral'
    p=ppoly(reshape([1.0_dp,2.0_dp],[1,2]),[0.0_dp,1.0_dp,2.0_dp],extrapolate=.false.)
    if (.not. ieee_is_nan(p%value(-1.0_dp))) error stop 'disabled extrapolation'
    p=ppoly(reshape([1.0_dp,2.0_dp],[1,2]),[0.0_dp,2.0_dp,1.0_dp],status=status)
    if (status==0 .or. p%valid()) error stop 'invalid breakpoints accepted'

    ! Different Hermite derivative counts at each node, scalar and matrix-valued.
    derivative_samples=0.0_dp
    derivative_samples(1,1:2)=[1.0_dp,1.5_dp]
    derivative_samples(2,1)=2.0_dp
    derivative_samples(3,1:3)=[3.0_dp,-2.0_dp,0.75_dp]
    b=bpoly_from_derivatives([0.0_dp,0.8_dp,2.0_dp],derivative_samples,[2,1,3],status=status)
    if (status/=0) error stop 'Hermite construction'
    if (abs(b%value(0.0_dp)-1.0_dp)>tol) error stop 'Hermite left value'
    if (abs(b%value(0.0_dp,nu=1)-1.5_dp)>tol) error stop 'Hermite left slope'
    if (abs(b%value(2.0_dp)-3.0_dp)>tol) error stop 'Hermite right value'
    if (abs(b%value(2.0_dp,nu=1)+2.0_dp)>tol) error stop 'Hermite right slope'
    if (abs(b%value(2.0_dp,nu=2)-0.75_dp)>tol) error stop 'Hermite right curvature'
    yh(:,:,1)=derivative_samples
    yh(:,:,2)=2.0_dp*derivative_samples
    ba=bpoly_from_derivatives([0.0_dp,0.8_dp,2.0_dp],yh,[2,1,3],status=status)
    if (status/=0 .or. ba%channels()/=2) error stop 'Hermite channel construction'
    val=ba%evaluate([0.1_dp,0.5_dp,1.4_dp])
    if (any(abs(val(:,2) - 2.0_dp*b%values([0.1_dp,0.5_dp,1.4_dp]))>tol)) &
        error stop 'Hermite channels' 

    ! Connect SciFort's existing CubicSpline to the general power basis.
    x=[0.0_dp,0.75_dp,1.5_dp,3.0_dp]
    y=sin(x)
    cubic=cubic_spline(x,y,bc_type='natural',status=status)
    if (status/=0) error stop 'CubicSpline construction'
    p=ppoly_from_cubic(cubic)
    if (.not. p%valid()) error stop 'CubicSpline conversion'
    query=[-0.1_dp,0.5_dp,1.25_dp,2.2_dp,3.1_dp]
    do i=1,size(query)
        if (abs(p%value(query(i))-cubic%evaluate(query(i)))>tol) error stop 'CubicSpline values'
        if (abs(p%value(query(i),nu=2)-cubic%derivative(query(i),order=2))>tol) &
            error stop 'CubicSpline derivative'
    end do
    y2(:,1)=y
    y2(:,2)=cos(x)
    cubic_channels=cubic_spline(x,y2,bc_type='clamped',status=status)
    if (status/=0) error stop 'multichannel CubicSpline construction'
    pa=ppoly_from_cubic(cubic_channels)
    if (.not. pa%valid() .or. pa%channels()/=2) error stop 'multichannel PPoly conversion'
    val=pa%evaluate(query)
    if (any(abs(val-cubic_channels%evaluate(query))>tol)) error stop 'multichannel PPoly values'

    print '(a)', 'piecewise PPoly/BPoly native regression PASS'
end program test_interpolate_piecewise
