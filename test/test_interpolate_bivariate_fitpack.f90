! SPDX-License-Identifier: MIT
program test_interpolate_bivariate_fitpack
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, &
        bisplev, bisplev_points, bispline_success, bispline_invalid_input
    implicit none
    type(fitpack_bivariate_tck_type) :: surface
    integer :: ix,iy,code
    real(dp) :: xx,yy
    real(dp), allocatable :: grid(:,:), points(:)
    surface%tx = [0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp]
    surface%ty = surface%tx
    surface%kx = 2
    surface%ky = 2
    allocate(surface%c(9))
    do ix=1,3
        xx=real(ix-1,dp)/2.0_dp
        do iy=1,3
            yy=real(iy-1,dp)/2.0_dp
            surface%c((ix-1)*3+iy) = xx+2.0_dp*yy+3.0_dp*xx*yy
        end do
    end do
    if (.not.surface%valid()) error stop 'bivariate TCK validity'
    xx=0.43_dp
    yy=0.62_dp
    if (abs(bisplev(xx,yy,surface)-(xx+2.0_dp*yy+3.0_dp*xx*yy))>1.0e-12_dp) &
        error stop 'bivariate polynomial reproduction'
    if (abs(bisplev(xx,yy,surface,dx=1)-(1.0_dp+3.0_dp*yy))>1.0e-12_dp) &
        error stop 'bivariate x derivative'
    if (abs(bisplev(xx,yy,surface,dy=1)-(2.0_dp+3.0_dp*xx))>1.0e-12_dp) &
        error stop 'bivariate y derivative'
    if (abs(bisplev(xx,yy,surface,dx=1,dy=1)-3.0_dp)>1.0e-12_dp) &
        error stop 'bivariate mixed derivative'
    if (abs(bisplev(-3.0_dp,yy,surface)-(2.0_dp*yy))>1.0e-12_dp) &
        error stop 'FITPACK clamps x outside knot domain'
    grid=bisplev([0.2_dp,0.5_dp],[0.3_dp,0.6_dp,0.8_dp],surface,status=code)
    if (code/=bispline_success .or. any(shape(grid)/=[2,3])) &
        error stop 'bivariate grid interface'
    if (abs(grid(2,3)-(0.5_dp+2.0_dp*0.8_dp+3.0_dp*0.5_dp*0.8_dp))>1.0e-12_dp) &
        error stop 'bivariate grid ordering'
    points=bisplev_points([0.2_dp,0.5_dp],[0.3_dp,0.6_dp],surface,status=code)
    if (code/=bispline_success .or. maxval(abs(points-[0.98_dp,2.6_dp]))>1.0e-12_dp) &
        error stop 'paired points'
    points=bisplev_points([0.2_dp],[0.3_dp,0.5_dp],surface,status=code)
    if (code/=bispline_invalid_input .or. size(points)/=0) error stop 'paired shape rejection'
    xx=bisplev(0.2_dp,0.4_dp,surface,dx=2,status=code)
    if (code/=bispline_invalid_input) error stop 'invalid derivative rejection'
    print *, 'Bivariate FITPACK TCK native PASS'
end program test_interpolate_bivariate_fitpack
