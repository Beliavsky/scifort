! SPDX-License-Identifier: MIT
! Deterministic stdin/stdout protocol for comparison against SciPy and BSpline.
program fitpack_periodic_parametric_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_fitpack_compat, only: fitpack_tck_type, &
        fitpack_parametric_tck_type, splev
    use scifort_optional_fitpack, only: periodic_univariate_spline, &
        parametric_smoothing_spline
    implicit none
    integer :: mode,n,ndim,deg,nq,i,status
    real(dp) :: s,r
    real(dp), allocatable :: u(:),w(:),q(:),x(:,:),v(:,:),y(:)
    type(fitpack_tck_type) :: t
    type(fitpack_parametric_tck_type) :: p
    read(*,*) mode,n,ndim,deg,nq,s
    allocate(u(n),w(n),q(nq),x(n,ndim))
    read(*,*) u
    read(*,*) w
    do i=1,ndim
        read(*,*) x(:,i)
    end do
    read(*,*) q
    if (mode==0) then
        if (ndim/=1) error stop 'periodic univariate expects one channel'
        t=periodic_univariate_spline(u,x(:,1),w=w,k=deg,s=s, &
            status=status,residual=r)
        write(*,'(i0)') status
        if (status/=0) stop
        write(*,'(es27.18)') r
        write(*,'(3(i0,1x))') size(t%t),size(t%c),t%k
        write(*,'(*(es27.18,1x))') t%t
        write(*,'(*(es27.18,1x))') t%c
        write(*,'(*(es27.18,1x))') u
        y=splev(q,t)
        write(*,'(*(es27.18,1x))') y
        y=splev(q,t,der=1)
        write(*,'(*(es27.18,1x))') y
    else
        p=parametric_smoothing_spline(x,u=u,w=w,k=deg,s=s, &
            per=(mode==2),status=status,residual=r)
        write(*,'(i0)') status
        if (status/=0) stop
        write(*,'(es27.18)') r
        write(*,'(3(i0,1x))') size(p%t),size(p%c,1),p%k
        write(*,'(*(es27.18,1x))') p%t
        do i=1,ndim
            write(*,'(*(es27.18,1x))') p%c(:,i)
        end do
        write(*,'(*(es27.18,1x))') p%u
        v=splev(q,p)
        do i=1,ndim
            write(*,'(*(es27.18,1x))') v(:,i)
        end do
        v=splev(q,p,der=1)
        do i=1,ndim
            write(*,'(*(es27.18,1x))') v(:,i)
        end do
    end if
end program fitpack_periodic_parametric_driver
