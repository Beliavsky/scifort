! SPDX-License-Identifier: MIT
! Independent SciPy 1.17 QR references, calculus/normal-equation identities,
! and old LSQUnivariateSpline compatibility.  No SciPy runtime needed.
program test_interpolate_make_lsq_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only: dp
    use scifort_interpolate_bspline, only: bspline_interpolator_type, &
        bspline_array_interpolator_type, make_lsq_spline, lsq_univariate_spline, &
        bspline_design_matrix, bspline_status_success, bspline_status_invalid_input, &
        bspline_status_singular
    implicit none
    integer, parameter :: n = 81, ndegree = 7
    integer, parameter :: degree(ndegree) = [0,1,2,3,5,9,15]
    real(dp), parameter :: inner(4) = [-0.6_dp,-0.2_dp,0.3_dp,0.7_dp]
    real(dp), parameter :: q(5) = [-0.97_dp,-0.31_dp,0.02_dp,0.4_dp,0.98_dp]
    real(dp) :: x(n), xx(n), y(n), yy(n,3), w(n), wpositive(n), expected(5,ndegree)
    real(dp) :: fit(n), v(3), normal_error, nanval, scale
    real(dp), allocatable :: knots(:), basis(:,:), coef(:,:), tbad(:)
    type(bspline_interpolator_type) :: spline, old, bad
    type(bspline_array_interpolator_type) :: multi, invalid_multi
    integer :: i, j, k, code

    expected(:,1) = [-0.96697129944929439_dp, -0.82280761663454594_dp, &
        0.18316647038580394_dp, 0.85737385430375435_dp, 0.80819571398100176_dp]
    expected(:,2) = [-0.81747071418445238_dp, -0.59682869821565054_dp, &
        0.13102126099007000_dp, 0.82091510276963831_dp, 0.68965228767259201_dp]
    expected(:,3) = [-0.75727485867735689_dp, -0.65150574268056138_dp, &
        0.17516050452020193_dp, 0.81936150453621803_dp, 0.68428662604850676_dp]
    expected(:,4) = [-0.76123120415951351_dp, -0.63215333710389054_dp, &
        0.16207844270109642_dp, 0.82262033611835794_dp, 0.68062611194913336_dp]
    expected(:,5) = [-0.76034065696166775_dp, -0.63656980815934605_dp, &
        0.16704051763088992_dp, 0.81924385137737721_dp, 0.68115228518271076_dp]
    expected(:,6) = [-0.76022145710772715_dp, -0.63699237702485534_dp, &
        0.16759673532940389_dp, 0.81889293934141483_dp, 0.68129464147521279_dp]
    expected(:,7) = [-0.76021970331581645_dp, -0.63699230881649715_dp, &
        0.16759835611480925_dp, 0.81889375445854551_dp, 0.68129658176505103_dp]

    do i=1,n
        x(i) = -1.0_dp + real(i-1,dp)/40.0_dp
        y(i) = sin(2.4_dp*x(i)) + 0.12_dp*cos(4.1_dp*x(i)) + 0.05_dp*x(i)**2
        w(i) = 0.7_dp + 0.03_dp*real(mod(i,13),dp)
    end do
    w(12) = 0.0_dp
    w(62) = 0.0_dp
    wpositive = w
    where (wpositive == 0.0_dp) wpositive = 1.0_dp
    yy(:,1) = y
    yy(:,2) = 2.0_dp*y + 0.3_dp
    yy(:,3) = 1.0_dp + 0.5_dp*x
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)

    do j=1,ndegree
        k = degree(j)
        knots = [spread(x(1),1,k+1),inner,spread(x(n),1,k+1)]
        spline = make_lsq_spline(x,y,knots,k=k,w=w,status=code,method='qr')
        if (code /= bspline_status_success .or. .not. spline%valid()) &
            error stop 'make_lsq_spline: scalar constructor failed'
        if (maxval(abs(spline%evaluate(q)-expected(:,j))) > 8.0e-9_dp) &
            error stop 'make_lsq_spline: SciPy reference mismatch'
        multi = make_lsq_spline(x,yy,knots,k=k,w=w,status=code)
        if (code /= bspline_status_success .or. .not. multi%valid()) &
            error stop 'make_lsq_spline: multichannel constructor failed'
        coef = multi%coefficients()
        if (any(shape(coef) /= [k+5,3])) error stop 'make_lsq_spline: coefficient shape'
        if (maxval(abs(coef(:,1)-spline%coefficients())) > 3.0e-10_dp) &
            error stop 'make_lsq_spline: shared factorization disagrees'
        do i=1,size(q)
            v = multi%evaluate(q(i))
            if (abs(v(1)-spline%evaluate(q(i))) > 1.0e-10_dp) &
                error stop 'make_lsq_spline: scalar-channel mismatch'
            if (abs(v(2)-(2.0_dp*v(1)+0.3_dp)) > 3.0e-10_dp) &
                error stop 'make_lsq_spline: affine channel mismatch'
            if (k >= 1) then
                if (abs(v(3)-(1.0_dp+0.5_dp*q(i))) > 4.0e-9_dp) &
                    error stop 'make_lsq_spline: linear reproduction'
            end if
        end do
        basis = bspline_design_matrix(x,knots,k,status=code)
        if (code /= bspline_status_success) error stop 'make_lsq_spline: design matrix'
        fit = matmul(basis,spline%coefficients())
        scale = max(1.0_dp,maxval(abs(y)))
        normal_error = maxval(abs(matmul(transpose(basis),w*w*(fit-y))))
        if (normal_error > 2.0e-9_dp*scale) &
            error stop 'make_lsq_spline: weighted least-squares stationarity'
        if (k == 3) then
            old = lsq_univariate_spline(x,y,inner,k=3,w=wpositive,status=code)
            spline = make_lsq_spline(x,y,knots,k=3,w=wpositive,status=code)
            if (code /= 0 .or. maxval(abs(old%coefficients()-spline%coefficients())) > 1.0e-10_dp) &
                error stop 'make_lsq_spline: existing LSQ API changed'
        end if
    end do

    ! Multiple identical sample locations and repeated *interior* knots are
    ! supported in the QR method, unlike legacy FITPACK strict-x fitting.
    k = 3
    knots = [spread(x(1),1,k+1),-0.3_dp,-0.3_dp,0.5_dp,spread(x(n),1,k+1)]
    xx = x
    xx(29) = xx(28)
    multi = make_lsq_spline(xx,yy,knots,k=k,w=w,status=code)
    if (code /= bspline_status_success .or. .not. multi%valid()) &
        error stop 'make_lsq_spline: repeated x and knot'
    spline = make_lsq_spline(x,y,knots,k=k,w=w,status=code,extrapolate=.false.)
    if (code /= bspline_status_success) error stop 'make_lsq_spline: no extrapolation'
    if (.not. ieee_is_nan(spline%evaluate(-1.01_dp))) &
        error stop 'make_lsq_spline: extrapolation policy'
    if (.not. ieee_is_nan(spline%evaluate(1.01_dp))) &
        error stop 'make_lsq_spline: right extrapolation policy'

    tbad = [-1.0_dp,0.101_dp,0.102_dp,1.0_dp]
    bad = make_lsq_spline(x,y,tbad,k=0,status=code)
    if (code /= bspline_status_singular .or. bad%valid()) &
        error stop 'make_lsq_spline: rank deficiency not reported'
    bad = make_lsq_spline(x,y,knots,k=3,w=-w,status=code)
    if (code /= bspline_status_invalid_input .or. bad%valid()) &
        error stop 'make_lsq_spline: negative weights accepted'
    bad = make_lsq_spline(x,y,knots,k=3,method='norm-eq',status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: unsupported solver accepted'
    bad = make_lsq_spline(x,y,knots,k=-1,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: negative degree accepted'
    bad = make_lsq_spline(x,y,knots,k=16,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: degree 16 accepted'
    bad = make_lsq_spline(x,y,knots(:size(knots)-1),k=3,status=code)
    if (code /= bspline_status_invalid_input .or. bad%valid()) &
        error stop 'make_lsq_spline: truncated knot vector accepted'
    bad = make_lsq_spline(x,y,[knots(1:5),-0.99_dp,knots(7:)],k=3,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: unsorted knots accepted'
    xx = x
    xx(30) = xx(29)-0.02_dp
    bad = make_lsq_spline(xx,y,knots,k=3,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: unsorted x accepted'
    xx = x
    xx(1) = -2.0_dp
    bad = make_lsq_spline(xx,y,knots,k=3,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: outside x accepted'
    invalid_multi = make_lsq_spline(x,yy(:n-1,:),knots,k=3,status=code)
    if (code /= bspline_status_invalid_input .or. invalid_multi%valid()) &
        error stop 'make_lsq_spline: wrong channel shape accepted'
    bad = make_lsq_spline(x,[nanval,y(2:)],knots,k=3,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: NaN values accepted'
    w(1) = nanval
    bad = make_lsq_spline(x,y,knots,k=3,w=w,status=code)
    if (code /= bspline_status_invalid_input) error stop 'make_lsq_spline: NaN weights accepted'

    print '(a)', 'PASS: make_lsq_spline explicit knots, degrees 0:15, QR, SciPy refs, validation'
end program test_interpolate_make_lsq_spline
