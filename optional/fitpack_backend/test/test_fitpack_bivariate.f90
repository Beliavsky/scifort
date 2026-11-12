! SPDX-License-Identifier: MIT
program test_fitpack_bivariate
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, &
        bisplev_points, bisplev
    use scifort_optional_fitpack, only : bivariate_smoothing_spline, &
        lsq_bivariate_spline, fitpack_smooth_success, fitpack_smooth_invalid_input
    implicit none
    integer, parameter :: n=82
    real(dp) :: x(n),y(n),z(n),w(n),d,ss,xi,yi,lo(2),hi(2)
    real(dp), allocatable :: fit(:)
    integer :: i,k,status
    type(fitpack_bivariate_tck_type) :: tck
    do i=1,n
        x(i) = modulo(real(i*37,dp)*0.091_dp,2.0_dp)-1.0_dp
        y(i) = modulo(real(i*53,dp)*0.077_dp,2.0_dp)-0.5_dp
        z(i) = sin(1.7_dp*x(i))+cos(1.9_dp*y(i))+0.2_dp*x(i)*y(i) &
             + 0.02_dp*sin(7.0_dp*x(i)+3.0_dp*y(i))
        w(i) = 0.8_dp+0.004_dp*real(i,dp)
    end do
    do k=1,5
        tck=bivariate_smoothing_spline(x,y,z,w=w,k=k,s=0.2_dp, &
            status=status,residual=ss)
        if (status/=fitpack_smooth_success .or. .not.tck%valid()) &
            error stop 'adaptive bivariate smoothing failed'
        fit=bisplev_points(x,y,tck,status=status)
        if (status/=0) error stop 'bivariate fit evaluation failed'
        d=sum((w*(z-fit))**2)
        if (abs(d-ss)>1.0e-8_dp) error stop 'bivariate smoothing residual mismatch'
        if (ss>0.2003_dp) error stop 'bivariate smoothing target exceeded'
        xi=0.1_dp
        yi=0.3_dp
        if (k>1) then
            if (abs(bisplev(xi,yi,tck,dx=1,dy=1))>1.0e4_dp) &
                error stop 'invalid bivariate mixed derivative'
        end if
        lo = [minval(x),minval(y)]
        hi = [maxval(x),maxval(y)]
        tck=lsq_bivariate_spline(x,y,z, &
            [lo(1)+0.5_dp*(hi(1)-lo(1))], &
            [lo(2)+0.5_dp*(hi(2)-lo(2))], &
            w=w,kx=k,ky=k,status=status,residual=ss)
        if (status/=fitpack_smooth_success .or. .not.tck%valid()) &
            error stop 'fixed-knot bivariate fit failed'
        fit=bisplev_points(x,y,tck,status=status)
        if (status/=0) error stop 'fixed-knot bivariate evaluation failed'
        d=sum((w*(z-fit))**2)
        if (abs(d-ss)>1.0e-8_dp) error stop 'fixed-knot residual mismatch'
        if (size(tck%tx)/=2*(k+1)+1) error stop 'fixed x knot lost'
        if (size(tck%ty)/=2*(k+1)+1) error stop 'fixed y knot lost'
    end do
    tck=bivariate_smoothing_spline(x(1:n-1),y,z,status=status)
    if (status/=fitpack_smooth_invalid_input .or. tck%valid()) &
        error stop 'input size invalid handling'
    tck=lsq_bivariate_spline(x,y,z,[minval(x)],[0.5_dp],status=status)
    if (status/=fitpack_smooth_invalid_input) error stop 'bad interior x knot handling'
    print *, 'FITPACK 2.0.0 bivariate adapter PASS'
end program test_fitpack_bivariate
