! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program rect_sphere_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_sphere, only : rect_sphere_bivariate_spline, &
        sphere_bivariate_spline_type
    implicit none
    type(sphere_bivariate_spline_type) :: fit
    real(dp), allocatable :: u(:),v(:),r(:,:),qt(:),qp(:),tu(:),tv(:),c(:)
    real(dp) :: poles(2)
    integer :: mu,mv,nq,i,status,flags(6)
    read(*,*) mu,mv,nq
    allocate(u(mu),v(mv),r(mu,mv),qt(nq),qp(nq))
    read(*,*) u
    read(*,*) v
    do i = 1, mu
        read(*,*) r(i,:)
    end do
    read(*,*) flags
    read(*,*) poles
    do i = 1, nq
        read(*,*) qt(i),qp(i)
    end do
    fit = rect_sphere_bivariate_spline(u,v,r,pole_continuity=flags(1:2)/=0, &
        pole_values=poles,pole_exact=flags(3:4)/=0,pole_flat=flags(5:6)/=0,status=status)
    print '(3(i0,1x),es26.17e3)', status,fit%rank(),fit%degrees_of_freedom(),fit%residual()
    if (.not. fit%valid()) stop
    call fit%get_knots(tu,tv)
    c = fit%coefficients()
    print '(*(es26.17e3,1x))', tu
    print '(*(es26.17e3,1x))', tv
    print '(*(es26.17e3,1x))', c
    do i = 1, nq
        print '(6(es26.17e3,1x))', fit%value(qt(i),qp(i)), &
            fit%value(qt(i),qp(i),dtheta=1),fit%value(qt(i),qp(i),dphi=1), &
            fit%value(qt(i),qp(i),dtheta=2),fit%value(qt(i),qp(i),dphi=2), &
            fit%value(qt(i),qp(i),dtheta=1,dphi=1)
    end do
end program rect_sphere_driver
