! SPDX-License-Identifier: MIT
program test_interpolate_periodic_bspline
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_interpolate_bspline, only : bspline_interpolator_type, &
        bspline_array_interpolator_type, make_interp_spline, &
        bspline_status_success, bspline_status_invalid_input
    implicit none
    integer, parameter :: n = 21
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: x(n), y(n,3), period, u, lhs, rhs, nanval, tolerance
    real(dp) :: v(3), v2(3)
    real(dp), allocatable :: c(:,:)
    type(bspline_interpolator_type) :: scalar, invalid, refined
    type(bspline_array_interpolator_type) :: multi, bad_multi
    integer :: i, j, k, deriv, status

    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    x(1) = -0.8_dp
    do i=2,n
        x(i) = x(i-1)+0.31_dp+0.06_dp*sin(real(i,dp))
    end do
    period = x(n)-x(1)
    do i=1,n
        u = 2.0_dp*pi*(x(i)-x(1))/period
        y(i,1) = sin(u)+0.15_dp*cos(3.0_dp*u)
        y(i,2) = cos(2.0_dp*u)+0.03_dp*sin(4.0_dp*u)
        y(i,3) = 2.0_dp
    end do
    y(n,:) = y(1,:)

    do k=2,15
        multi = make_interp_spline(x,y,k=k,bc_type='periodic',status=status)
        if (status /= bspline_status_success .or. .not. multi%valid()) &
            error stop 'periodic matrix construction failed'
        c = multi%coefficients()
        if (size(c,1) /= n+k-1 .or. size(c,2) /= 3) &
            error stop 'wrong periodic spline coefficient dimensions'
        do j=1,3
            scalar = make_interp_spline(x,y(:,j),k=k,bc_type='periodic',status=status)
            if (status /= bspline_status_success .or. .not. scalar%valid()) &
                error stop 'periodic scalar construction failed'
            if (maxval(abs(scalar%coefficients()-c(:,j))) > 5.0e-10_dp) &
                error stop 'scalar and channel coefficient mismatch'
            if (maxval(abs(scalar%evaluate(x)-y(:,j))) > 1.0e-9_dp) &
                error stop 'periodic interpolation fails at observations'
        end do
        do deriv=0,k-1
            v = multi%derivative(x(1),deriv)
            v2 = multi%derivative(x(n),deriv)
            tolerance = 2.0e-7_dp*max(1.0_dp,maxval(abs(v)),maxval(abs(v2)))
            if (maxval(abs(v-v2)) > tolerance) &
                error stop 'periodic derivative mismatch across endpoints'
        end do
        u = x(7)+0.12_dp
        do deriv=0,min(3,k)
            v = multi%derivative(u,deriv)
            v2 = multi%derivative(u+3.0_dp*period,deriv)
            tolerance = 5.0e-8_dp*max(1.0_dp,maxval(abs(v)))
            if (maxval(abs(v-v2)) > tolerance) &
                error stop 'periodic exterior evaluation'
        end do
        v = multi%integrate(x(1),x(n))
        v2 = multi%integrate(x(1)-period,x(n)+2.0_dp*period)
        if (maxval(abs(v2-4.0_dp*v)) > 3.0e-8_dp) &
            error stop 'periodic integration across cycles'
        if (k == 3) then
            scalar = make_interp_spline(x,y(:,1),k=k,bc_type='periodic')
            refined = scalar%insert_knot(x(8)+0.03_dp,m=2,status=status)
            if (status /= bspline_status_success) &
                error stop 'periodic knot insertion failed'
            do i=1,n
                if (abs(refined%evaluate(x(i))-scalar%evaluate(x(i))) > 1.0e-9_dp) &
                    error stop 'knot insertion changed periodic interior'
            end do
        end if
    end do

    scalar = make_interp_spline(x,y(:,1),bc_type='periodic',extrapolate=.false.,status=status)
    if (status /= bspline_status_success) error stop 'finite interval periodic construction'
    if (.not. ieee_is_nan(scalar%evaluate(x(1)-period))) &
        error stop 'extrapolate false must not wrap'
    if (abs(scalar%evaluate(x(2))-y(2,1)) > 2.0e-11_dp) &
        error stop 'extrapolate false changed interior'

    ! SciPy's k=0/1 construction uses ordinary pieces even with periodic BC.
    scalar = make_interp_spline(x,y(:,1),k=0,bc_type='periodic',status=status)
    if (status /= bspline_status_success) error stop 'constant periodic alias'
    if (maxval(abs(scalar%evaluate(x)-y(:,1))) > 1.0e-12_dp) &
        error stop 'degree-zero interpolation'

    ! SciPy's k=1 construction is the standard linear spline even with
    ! bc_type='periodic'; only the data endpoints are checked for equality.
    scalar = make_interp_spline(x,y(:,1),k=1,bc_type='periodic',status=status)
    if (status /= bspline_status_success) error stop 'linear periodic alias'
    lhs = scalar%evaluate(x(1)-0.1_dp)
    rhs = y(1,1)-0.1_dp*(y(2,1)-y(1,1))/(x(2)-x(1))
    if (abs(lhs-rhs)>2.0e-11_dp) error stop 'degree-one extrapolation convention'

    invalid = make_interp_spline(x,y(:,1)+[spread(0.0_dp,1,n-1),0.4_dp], &
        k=3,bc_type='periodic',status=status)
    if (status /= bspline_status_invalid_input .or. invalid%valid()) &
        error stop 'mismatched periodic endpoints accepted'
    invalid = make_interp_spline([x(1),x(1),x(3:)],y(:,1),k=3,bc_type='periodic',status=status)
    if (status /= bspline_status_invalid_input) error stop 'duplicate site accepted'
    invalid = make_interp_spline(x,y(:,1),k=3,bc_type='periodic',knots=x,status=status)
    if (status /= bspline_status_invalid_input) error stop 'explicit periodic knots accepted'
    invalid = make_interp_spline(x,y(:,1),k=-1,bc_type='periodic',status=status)
    if (status /= bspline_status_invalid_input) error stop 'negative degree accepted'
    invalid = make_interp_spline(x,y(:,1),k=16,bc_type='periodic',status=status)
    if (status /= bspline_status_invalid_input) error stop 'unsupported degree accepted'
    bad_multi = make_interp_spline(x,reshape(y(:,1),[n,1]),k=3,bc_type='periodic',status=status)
    if (status /= bspline_status_success .or. .not. bad_multi%valid()) &
        error stop 'single-channel periodic matrix rejected'
    bad_multi = make_interp_spline(x,y(:,1:2),k=3,bc_type='periodic', &
        knots=x,status=status)
    if (status /= bspline_status_invalid_input .or. bad_multi%valid()) &
        error stop 'matrix explicit knots accepted'
    bad_multi = make_interp_spline(x,reshape([nanval,y(2:n,1)],[n,1]), &
        k=3,bc_type='periodic',status=status)
    if (status /= bspline_status_invalid_input .or. bad_multi%valid()) &
        error stop 'nonfinite periodic data accepted'

    print '(a)', 'PASS: periodic B-spline construction, calculus, channels, validation'
end program test_interpolate_periodic_bspline
