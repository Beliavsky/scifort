! SPDX-License-Identifier: MIT
! Native regression for roots, B-spline conversion and complex PPoly/BPoly.
program test_interpolate_piecewise_roots_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly, bpoly, ppoly_type, complex_ppoly_type, &
        complex_bpoly_type, complex_bpoly_from_power, complex_ppoly_from_bernstein, &
        bspline, bspline_interpolator_type, bspline_array_interpolator_type, ppoly_from_spline
    implicit none
    type(ppoly_type) :: p, pp
    type(complex_ppoly_type) :: cp, cpp
    type(complex_bpoly_type) :: cb, cbb
    type(bspline_interpolator_type) :: s
    type(bspline_array_interpolator_type) :: sa
    real(dp), allocatable :: roots(:), vals(:,:), knots(:)
    real(dp) :: t(11), coeff(7), q
    complex(dp) :: z(2,2), expected, integrated
    complex(dp), allocatable :: cv(:,:)
    integer :: j, status
    real(dp), parameter :: tol = 5.0e-11_dp

    ! Even-multiplicity root, with a second strictly positive polynomial.
    p = ppoly(reshape([1.0_dp,-0.5_dp,0.0625_dp, &
        1.0_dp,0.0_dp,1.0_dp],[3,2]),[0.0_dp,1.0_dp,2.0_dp],extrapolate=.false.)
    roots = p%roots()
    if (size(roots)/=1) error stop 'tangential root count'
    if (abs(roots(1)-0.25_dp)>tol) error stop 'tangential root location'

    p = ppoly(reshape([1.0_dp,-1.0_dp],[1,2]),[0.0_dp,1.0_dp,2.0_dp], &
        extrapolate=.false.)
    roots = p%roots(discontinuity=.false.)
    if (size(roots)/=0) error stop 'unwanted discontinuity root'
    roots = p%roots(discontinuity=.true.)
    if (size(roots)/=1) error stop 'discontinuity root count'
    if (roots(1)/=1.0_dp) error stop 'discontinuity root position'

    p = ppoly(reshape([0.0_dp,1.0_dp],[1,2]),[0.0_dp,1.0_dp,2.0_dp])
    roots = p%roots()
    if (size(roots)/=2) error stop 'zero section root count'
    if (roots(1)/=0.0_dp .or. .not. ieee_is_nan(roots(2))) error stop 'zero section sentinel'
    roots = p%solve(y=1.0_dp,extrapolate=.false.)
    if (size(roots)/=2) error stop 'shifted polynomial root count'
    if (abs(roots(1)-1.0_dp)>tol .or. .not. ieee_is_nan(roots(2))) &
        error stop 'shifted polynomial zero section'

    p = ppoly(reshape([1.0_dp,0.5_dp,1.0_dp,0.25_dp],[2,2]), &
        [3.0_dp,2.0_dp,1.0_dp],extrapolate=.false.)
    roots = p%roots(discontinuity=.false.)
    if (size(roots)/=2) error stop 'descending roots count'
    if (abs(roots(1)-2.5_dp)>tol .or. abs(roots(2)-1.75_dp)>tol) &
        error stop 'descending roots location'

    ! Cubic clamped B-splines, including repeated endpoints and multiple channels.
    t = [0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.5_dp,1.0_dp,1.5_dp, &
        2.0_dp,2.0_dp,2.0_dp,2.0_dp]
    coeff = [0.1_dp,0.2_dp,-0.3_dp,0.7_dp,0.1_dp,0.0_dp,0.5_dp]
    s = bspline(t,coeff,3,status=status)
    if (status/=0) error stop 'test BSpline construction'
    p = ppoly_from_spline(s,status=status)
    if (status/=0 .or. .not. p%valid()) error stop 'PPoly from scalar B-spline'
    knots = p%knots()
    if (size(knots)/=size(t) .or. any(knots/=t)) error stop 'PPoly from spline knots'
    do j=1,81
        q=-0.3_dp+real(j,dp)*0.03125_dp
        if (abs(p%value(q)-s%evaluate(q))>tol) error stop 'from spline evaluation'
        if (abs(p%value(q,nu=2)-s%derivative(q,nu=2))>tol) &
            error stop 'from spline derivative'
    end do
    sa=bspline(t,reshape([coeff,2.0_dp*coeff],[7,2]),3,status=status)
    if (status/=0) error stop 'test multichannel BSpline construction'
    pp=ppoly_from_spline(sa,status=status)
    if (status/=0 .or. pp%channels()/=2) error stop 'PPoly from multichannel BSpline'
    vals=pp%evaluate([0.2_dp,0.7_dp,1.6_dp])
    if (any(abs(vals(:,1)*2.0_dp-vals(:,2))>tol)) error stop 'from spline channel ratio'

    ! Complex piecewise polynomial and Bernstein conversion.
    z(:,1) = [cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-3.0_dp,dp)]
    z(:,2) = z(:,1)
    z(2,2) = z(2,2)+z(1,1)
    cp=ppoly(z,[0.0_dp,1.0_dp,2.0_dp],status=status)
    if (status/=0 .or. .not. cp%valid()) error stop 'complex power constructor'
    expected=z(1,1)*0.3_dp+z(2,1)
    if (abs(cp%value(0.3_dp)-expected)>tol) error stop 'complex power value'
    cpp=cp%derivative(1)
    if (abs(cpp%value(0.3_dp)-z(1,1))>tol) error stop 'complex power derivative'
    cpp=cp%antiderivative(2)
    cpp=cpp%derivative(2)
    if (abs(cpp%value(1.4_dp)-cp%value(1.4_dp))>tol) &
        error stop 'complex power repeated integration'
    cb=complex_bpoly_from_power(cp)
    cpp=complex_ppoly_from_bernstein(cb)
    cbb=bpoly(z,[0.0_dp,1.0_dp,2.0_dp],status=status)
    if (status/=0 .or. .not. cbb%valid()) error stop 'complex Bernstein constructor'
    do j=1,51
        q=-0.3_dp+real(j,dp)*0.055_dp
        if (abs(cb%value(q)-cp%value(q))>tol) error stop 'complex Bernstein conversion'
        if (abs(cpp%value(q)-cp%value(q))>tol) error stop 'complex power basis roundtrip'
    end do
    integrated=cp%integral(0.0_dp,2.0_dp)
    if (abs(integrated-cb%integral(0.0_dp,2.0_dp))>tol) &
        error stop 'complex Bernstein definite integration'
    cp=ppoly(reshape([cmplx(1.0_dp,2.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp), &
        cmplx(3.0_dp,4.0_dp,dp),cmplx(4.0_dp,5.0_dp,dp)], [1,2,2]), &
        [0.0_dp,1.0_dp,2.0_dp],status=status)
    if (status/=0 .or. cp%channels()/=2) error stop 'complex channel construction'
    cv=cp%evaluate([0.1_dp,0.5_dp,1.4_dp])
    if (size(cv,1)/=3 .or. size(cv,2)/=2) error stop 'complex channel evaluation'

    print '(a)', 'piecewise roots/complex native regression PASS'
end program test_interpolate_piecewise_roots_complex
