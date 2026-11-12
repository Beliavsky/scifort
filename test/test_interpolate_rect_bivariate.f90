! SPDX-License-Identifier: MIT
program test_interpolate_rect_bivariate
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : rect_bivariate_spline, &
        fitpack_bivariate_tck_type, bisplev, bisplev_points, bispline_integral, &
        bispline_success, bispline_invalid_input
    implicit none
    integer :: i, j, code, kx, ky
    real(dp) :: x(9), y(11), z(9,11), xx, yy, exact, result, nanval
    real(dp) :: xp(6), yp(7), zz(6,7)
    type(fitpack_bivariate_tck_type) :: fit, bad
    real(dp), allocatable :: grid(:,:), points(:)

    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    x = [-1.3_dp,-1.0_dp,-0.7_dp,-0.1_dp,0.25_dp,0.7_dp,1.0_dp,1.65_dp,2.0_dp]
    y = [-2.0_dp,-1.4_dp,-0.9_dp,-0.4_dp,0.1_dp,0.3_dp,0.7_dp,1.1_dp,1.8_dp,2.2_dp,2.9_dp]
    do j = 1,size(y)
        do i = 1,size(x)
            z(i,j) = 1.0_dp + 2.0_dp*x(i) - 3.0_dp*y(j) + &
                x(i)*y(j)**2 + 0.4_dp*x(i)**2*y(j)
        end do
    end do
    fit = rect_bivariate_spline(x,y,z,status=code)
    if (code /= bispline_success .or. .not. fit%valid()) &
        error stop 'default degree rectangular fit'
    if (fit%kx /= 3 .or. fit%ky /= 3) error stop 'default degrees'
    do j = 1,size(y)
        do i = 1,size(x)
            if (abs(bisplev(x(i),y(j),fit)-z(i,j)) > 3.0e-11_dp) &
                error stop 'interpolates each grid observation'
        end do
    end do
    xx = 0.46_dp
    yy = -0.28_dp
    exact = 1.0_dp + 2.0_dp*xx - 3.0_dp*yy + xx*yy**2 + 0.4_dp*xx**2*yy
    if (abs(bisplev(xx,yy,fit)-exact) > 2.0e-11_dp) &
        error stop 'off-grid bicubic polynomial'
    if (abs(bisplev(xx,yy,fit,dx=1)-(2.0_dp+yy**2+0.8_dp*xx*yy)) > 2.0e-10_dp) &
        error stop 'partial derivative x'
    if (abs(bisplev(xx,yy,fit,dy=1)-(-3.0_dp+2.0_dp*xx*yy+0.4_dp*xx**2)) > 2.0e-10_dp) &
        error stop 'partial derivative y'
    if (abs(bisplev(xx,yy,fit,dx=1,dy=1)-(2.0_dp*yy+0.8_dp*xx)) > 2.0e-10_dp) &
        error stop 'mixed partial derivative'
    grid = bisplev(x(2:5),y(3:6),fit,status=code)
    if (code /= bispline_success .or. any(shape(grid) /= [4,4])) &
        error stop 'rectangular grid evaluation'
    if (maxval(abs(grid-z(2:5,3:6))) > 3.0e-11_dp) &
        error stop 'rectangular grid values'
    points = bisplev_points(x(2:5),y(3:6),fit,status=code)
    if (code /= bispline_success) error stop 'paired evaluation status'
    do i = 1,size(points)
        if (abs(points(i)-z(i+1,i+2)) > 3.0e-11_dp) &
            error stop 'paired evaluation values'
    end do

    ! Both axes integrate exactly. The rectangle limits are clipped to the
    ! spline support, and reversing either integration bound changes sign.
    exact = 0.0_dp
    ! Integral of 1+2*x-3*y+x*y**2+0.4*x**2*y over [0,1]x[-1,2].
    exact = 6.0_dp - 4.5_dp + 1.5_dp + 0.2_dp
    result = fit%integral(0.0_dp,1.0_dp,-1.0_dp,2.0_dp,status=code)
    if (code /= bispline_success .or. abs(result-exact) > 3.0e-11_dp) &
        error stop 'exact rectangular integral'
    if (abs(bispline_integral(1.0_dp,0.0_dp,-1.0_dp,2.0_dp,fit)+exact) > 3.0e-11_dp) &
        error stop 'x integral orientation'
    if (abs(fit%integral(1.0_dp,0.0_dp,2.0_dp,-1.0_dp)-exact) > 3.0e-11_dp) &
        error stop 'double orientation reversal'
    if (abs(fit%integral(-200.0_dp,0.0_dp,-100.0_dp,0.0_dp)- &
            fit%integral(x(1),0.0_dp,y(1),0.0_dp)) > 3.0e-11_dp) &
        error stop 'integral domain clipping'
    if (fit%integral(4.0_dp,5.0_dp,-1.0_dp,0.0_dp) /= 0.0_dp) &
        error stop 'no overlap integral'
    if (fit%integral(0.0_dp,0.0_dp,-1.0_dp,2.0_dp) /= 0.0_dp) &
        error stop 'zero width integral'
    result = fit%integral(0.0_dp,nanval,-1.0_dp,2.0_dp,status=code)
    if (code /= bispline_invalid_input .or. .not. ieee_is_nan(result)) &
        error stop 'invalid integral limits'
    result = bad%integral(0.0_dp,1.0_dp,-1.0_dp,2.0_dp,status=code)
    if (code /= bispline_invalid_input .or. .not. ieee_is_nan(result)) &
        error stop 'invalid TCK for integral'

    ! Exercise each degree and asymmetric degrees at minimal grid sizes.
    do kx = 1,5
        do ky = 1,5
            do i = 1,6
                xp(i) = real(i-1,dp) + 0.02_dp*real(i*i,dp)
            end do
            do j = 1,7
                yp(j) = 0.3_dp*real(j-1,dp) + 0.005_dp*real(j*j,dp)
            end do
            do j = 1,7
                do i = 1,6
                    zz(i,j) = sin(xp(i)) + cos(yp(j)) + 0.4_dp*xp(i)*yp(j)
                end do
            end do
            fit = rect_bivariate_spline(xp(:kx+1),yp(:ky+1), &
                zz(:kx+1,:ky+1),kx=kx,ky=ky,status=code)
            if (code /= bispline_success .or. .not. fit%valid()) &
                error stop 'valid degree/minimum size'
            do j = 1,ky+1
                do i = 1,kx+1
                    if (abs(bisplev(xp(i),yp(j),fit)-zz(i,j)) > 2.0e-10_dp) &
                        error stop 'minimum degree grid reconstruction'
                end do
            end do
        end do
    end do

    bad = rect_bivariate_spline(x(1:2),y(1:4),z(1:2,1:4),status=code)
    if (code /= bispline_invalid_input .or. bad%valid()) error stop 'too few x sites'
    bad = rect_bivariate_spline(x,y,z(1:8,:),status=code)
    if (code /= bispline_invalid_input .or. bad%valid()) error stop 'invalid shape'
    bad = rect_bivariate_spline(x,y,z,kx=0,status=code)
    if (code /= bispline_invalid_input .or. bad%valid()) error stop 'invalid degree'
    bad = rect_bivariate_spline([x(1),x(1),x(3:)],y,z,status=code)
    if (code /= bispline_invalid_input .or. bad%valid()) error stop 'duplicate x site'
    bad = rect_bivariate_spline(x,[y(1),nanval,y(3:)],z,status=code)
    if (code /= bispline_invalid_input .or. bad%valid()) error stop 'nonfinite y site'
    z(1,1) = nanval
    bad = rect_bivariate_spline(x,y,z,status=code)
    if (code /= bispline_invalid_input .or. bad%valid()) error stop 'nonfinite observation'

    print '(a)', 'PASS: rectangular bivariate spline interpolation and exact integration'
end program test_interpolate_rect_bivariate
