! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program sphere_spline_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_sphere, only : lsq_sphere_bivariate_spline, sphere_bivariate_spline_type
    implicit none
    type(sphere_bivariate_spline_type) :: fit
    real(dp), allocatable :: theta(:),phi(:),r(:),w(:),tt(:),tp(:),qt(:),qp(:),tx(:),ty(:),c(:)
    real(dp) :: eps
    integer :: m, nt, np, nq, i, status
    read(*,*) m,nt,np,nq,eps
    allocate(theta(m),phi(m),r(m),w(m),tt(nt),tp(np),qt(nq),qp(nq))
    do i = 1, m
        read(*,*) theta(i),phi(i),r(i),w(i)
    end do
    if (nt > 0) read(*,*) tt
    read(*,*) tp
    do i = 1, nq
        read(*,*) qt(i),qp(i)
    end do
    fit = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w,eps,status)
    print '(3(i0,1x),es26.17e3)', status,fit%rank(),fit%degrees_of_freedom(),fit%residual()
    if (.not. fit%valid()) stop
    call fit%get_knots(tx,ty)
    c = fit%coefficients()
    print '(*(es26.17e3,1x))', tx
    print '(*(es26.17e3,1x))', ty
    print '(*(es26.17e3,1x))', c
    do i = 1, nq
        print '(6(es26.17e3,1x))', fit%value(qt(i),qp(i)), &
            fit%value(qt(i),qp(i),dtheta=1), fit%value(qt(i),qp(i),dphi=1), &
            fit%value(qt(i),qp(i),dtheta=2), fit%value(qt(i),qp(i),dphi=2), &
            fit%value(qt(i),qp(i),dtheta=1,dphi=1)
    end do
end program sphere_spline_driver
