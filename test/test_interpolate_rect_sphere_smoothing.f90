! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_rect_sphere_smoothing
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_constants, only : pi => scifort_pi
    use scifort_interpolate, only : rect_sphere_bivariate_spline, &
        sphere_bivariate_spline_type, sphere_spline_numerical_failure
    implicit none
    type(sphere_bivariate_spline_type) :: fit,scaled,bad
    real(dp) :: u(12),v(16),r(12,16),targets(4),fp,previous,q,base,scale
    real(dp), allocatable :: grid(:,:),tu(:),tv(:)
    integer :: i,j,k,status,checks
    checks = 0
    do i = 1, size(u)
        u(i) = pi*real(i,dp)/13.0_dp
    end do
    do j = 1, size(v)
        v(j) = -0.7_dp+2.0_dp*pi*real(j-1,dp)/16.0_dp
        do i = 1, size(u)
            r(i,j) = cos(u(i))+0.7_dp*sin(u(i))*cos(v(j)) &
                +0.04_dp*sin(17.0_dp*u(i)+5.0_dp*v(j))
        end do
    end do
    targets = [0.01_dp,0.2_dp,2.0_dp,1000.0_dp]
    previous = -1.0_dp
    do k = 1, size(targets)
        fit = rect_sphere_bivariate_spline(u,v,r,s=targets(k), &
            pole_continuity=[.true.,.true.],pole_values=[1.0_dp,-1.0_dp], &
            pole_exact=[.true.,.true.],status=status)
        call require(status == 0 .and. fit%valid(),'adaptive fit succeeds')
        grid = fit%grid(u,v,status=status)
        call require(status == 0,'grid evaluation succeeds')
        fp = sum((grid-r)**2)
        call near(fp,fit%residual(),2e-11_dp*(1.0_dp+fp),'reported grid residual')
        call require(fp <= 1.0011_dp*targets(k),'residual target upper bound')
        call require(fp > previous,'residual grows with smoothing')
        previous = fp
        if (k < 4) call near(fp,targets(k),0.0011_dp*targets(k),'active target')
        call near(fit%value(0.0_dp,0.3_dp),1.0_dp,3e-12_dp,'exact north pole')
        call near(fit%value(pi,0.3_dp),-1.0_dp,3e-12_dp,'exact south pole')
        q = fit%value(0.8_dp,0.2_dp)
        call near(fit%value(0.8_dp,0.2_dp+4*pi),q,3e-12_dp,'periodic smoothing')
        call require(ieee_is_finite(fit%value(0.8_dp,0.2_dp,dtheta=1,dphi=1)), &
            'finite mixed derivative')
        call require(fit%rank() == -1,'rank not estimated')
        call fit%get_knots(tu,tv)
        call require(size(tu) <= size(u)+8 .and. size(tv) <= size(v)+7,'knot limits')
        if (k == 4) then
            call require(size(tu) == 8 .and. size(tv) == 8,'polynomial limit')
            call require(fit%degrees_of_freedom() == 2,'polynomial model dimension')
        end if
    end do
    fit = rect_sphere_bivariate_spline(u,v,r,s=0.2_dp,status=status)
    call require(status == 0,'default pole smoothing')
    base = fit%value(0.8_dp,0.2_dp)
    fp = fit%residual()
    do k = 1, 2
        scale = 1e150_dp
        if (k == 2) scale = 1e-150_dp
        scaled = rect_sphere_bivariate_spline(u,v,r*scale,s=0.2_dp*scale*scale,status=status)
        call require(status == 0 .and. scaled%valid(),'scaled smoothing succeeds')
        call near(scaled%value(0.8_dp,0.2_dp)/scale,base,3e-9_dp,'scaled values')
        call near((scaled%residual()/scale)/scale,fp,3e-9_dp,'scaled residual')
    end do
    fit = rect_sphere_bivariate_spline(u,v,0*r,s=0.2_dp,status=status)
    call require(status == 0 .and. fit%valid(),'zero-data smoothing')
    call near(fit%residual(),0.0_dp,0.0_dp,'zero-data residual')
    fit = rect_sphere_bivariate_spline(u,v,0*r+2.5_dp,s=huge(1.0_dp),status=status)
    call require(status == 0 .and. fit%valid(),'huge target polynomial limit')
    call near(fit%value(0.8_dp,0.2_dp),2.5_dp,3e-12_dp,'constant smoothing')
    bad = rect_sphere_bivariate_spline(u,v,r*1e200_dp,s=1e-200_dp,status=status)
    call require(status == sphere_spline_numerical_failure,'scaled target underflow guard')
    call require(.not. bad%valid(),'underflow returns invalid fit')
    bad = rect_sphere_bivariate_spline(u,v,r,s=1e-25_dp, &
        pole_continuity=[.true.,.true.],pole_values=[1.0_dp,-1.0_dp], &
        pole_exact=[.true.,.true.],status=status)
    call require(status == sphere_spline_numerical_failure,'unattainable tolerance failure')
    call require(.not. bad%valid(),'failed convergence returns invalid fit')
    call require(ieee_is_nan(bad%residual()),'invalid residual is NaN')
    print '(a,i0,a)', 'Rectangular spherical smoothing: ',checks,' checks passed.'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
    subroutine near(got,ref,tol,label)
        implicit none
        real(dp), intent(in) :: got !! Computed result.
        real(dp), intent(in) :: ref !! Reference result.
        real(dp), intent(in) :: tol !! Absolute error bound.
        character(len=*), intent(in) :: label !! Failure description.
        call require(ieee_is_finite(got),label)
        call require(abs(got-ref) <= tol,label)
    end subroutine near
end program test_interpolate_rect_sphere_smoothing
