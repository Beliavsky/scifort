! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program smooth_sphere_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_sphere, only : smooth_sphere_bivariate_spline, &
        sphere_bivariate_spline_type
    implicit none
    integer :: m,nq,ntest,npest,io
    do
        read(*,*,iostat=io) m,nq,ntest,npest
        if (io /= 0) exit
        block
            type(sphere_bivariate_spline_type) :: fit
            real(dp), allocatable :: t(:),p(:),r(:),w(:),qt(:),qp(:),tt(:),tp(:),c(:)
            real(dp) :: s,eps
            integer :: i,status
            allocate(t(m),p(m),r(m),w(m),qt(nq),qp(nq))
            read(*,*) s,eps
            read(*,*) t
            read(*,*) p
            read(*,*) r
            read(*,*) w
            do i = 1, nq
                read(*,*) qt(i),qp(i)
            end do
            fit = smooth_sphere_bivariate_spline(t,p,r,w=w,s=s,eps=eps, &
                ntest=ntest,npest=npest,status=status)
            print '(3(i0,1x),es26.17e3)', &
                status,fit%rank(),fit%degrees_of_freedom(),fit%residual()
            if (.not. fit%valid()) cycle
            call fit%get_knots(tt,tp)
            c = fit%coefficients()
            print '(*(es26.17e3,1x))', tt
            print '(*(es26.17e3,1x))', tp
            print '(*(es26.17e3,1x))', c
            do i = 1, nq
                print '(6(es26.17e3,1x))', fit%value(qt(i),qp(i)), &
                    fit%value(qt(i),qp(i),dtheta=1),fit%value(qt(i),qp(i),dphi=1), &
                    fit%value(qt(i),qp(i),dtheta=2),fit%value(qt(i),qp(i),dphi=2), &
                    fit%value(qt(i),qp(i),dtheta=1,dphi=1)
            end do
        end block
    end do
end program smooth_sphere_driver
