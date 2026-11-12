! SPDX-License-Identifier: MIT
! Numeric regression driver for FITPACK-style parametric splines.
program fitpack_parametric_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_compat, only : fitpack_parametric_tck_type, &
        splprep, splev
    implicit none
    type(fitpack_parametric_tck_type) :: spline
    integer :: n, ndim, nq, degree, task, nt, provided_u, code, i, j, der
    real(dp), allocatable :: x(:,:), w(:), u(:), t(:), q(:), ans(:,:), point(:)
    read(*,*) n, ndim, nq, degree, task, nt, provided_u
    allocate(x(n,ndim),w(n),q(nq),u(n),t(nt))
    do i=1,n
        read(*,*) x(i,:)
    end do
    read(*,*) w
    if (provided_u /= 0) read(*,*) u
    if (nt > 0) read(*,*) t
    read(*,*) q
    if (provided_u /= 0) then
        if (task == -1) then
            spline = splprep(x,u=u,w=w,k=degree,s=0.0_dp,task=-1,t=t,status=code)
        else
            spline = splprep(x,u=u,w=w,k=degree,s=0.0_dp,status=code)
        end if
    else
        if (task == -1) then
            spline = splprep(x,w=w,k=degree,s=0.0_dp,task=-1,t=t,status=code)
        else
            spline = splprep(x,w=w,k=degree,s=0.0_dp,status=code)
        end if
    end if
    print '(i0)', code
    if (code /= 0) stop
    print '(3(i0,1x))', size(spline%t),size(spline%c,1),spline%k
    print '(*(es26.17e3,1x))',spline%u
    print '(*(es26.17e3,1x))',spline%t
    do j=1,ndim
        print '(*(es26.17e3,1x))',spline%c(:,j)
    end do
    do der=0,min(3,degree)
        ans = splev(q,spline,der=der,status=code)
        print '(i0)',code
        do i=1,nq
            print '(*(es26.17e3,1x))', ans(i,:)
        end do
    end do
    point = splev(q(1),spline,der=1,status=code)
    print '(i0)',code
    print '(*(es26.17e3,1x))',point
end program fitpack_parametric_driver
