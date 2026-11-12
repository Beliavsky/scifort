! SPDX-License-Identifier: MIT
program fitpack_tck_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_compat, only : fitpack_tck_type, splrep, splev, &
        splint, splder, splantider, spalde
    implicit none
    integer :: n,nq,k,nt,task,status,i,j,order
    real(dp) :: bbox(2),h
    real(dp), allocatable :: x(:),y(:),w(:),t(:),q(:),ds(:)
    real(dp), allocatable :: vals(:)
    type(fitpack_tck_type) :: tck,dt,at
    read(*,*) n,nq,k,task,nt
    allocate(x(n),y(n),w(n),t(nt),q(nq))
    read(*,*) bbox
    read(*,*) x
    read(*,*) y
    read(*,*) w
    if (nt > 0) read(*,*) t
    read(*,*) q
    if (task == -1) then
        tck = splrep(x,y,w=w,k=k,task=-1,t=t,bbox=bbox,status=status)
    else
        tck = splrep(x,y,w=w,k=k,s=0.0_dp,status=status)
    end if
    write(*,'(i0)') status
    if (status /= 0) stop
    write(*,'(3(i0,1x))') size(tck%t),size(tck%c),tck%k
    write(*,'(*(es27.18,1x))') tck%t
    write(*,'(*(es27.18,1x))') tck%c
    do order=0,min(k,3)
        vals=splev(q,tck,der=order,status=status)
        write(*,'(*(es27.18,1x))') vals
    end do
    do j=1,3
        vals=splev(q,tck,ext=j,status=status)
        ! ext=2 is a deliberate out-of-bounds result.
        write(*,'(i0)') status
        if (status == 0) write(*,'(*(es27.18,1x))') vals
    end do
    write(*,'(es27.18)') splint(q(1),q(nq),tck)
    h=0.5_dp*(x(2)+x(3))
    ds=spalde(h,tck,status=status)
    write(*,'(i0)') status
    if (status == 0) write(*,'(*(es27.18,1x))') ds
    dt=splder(tck,status=status)
    write(*,'(i0)') status
    if (status == 0) then
        vals=splev(q,dt)
        write(*,'(*(es27.18,1x))') vals
    end if
    at=splantider(tck,status=status)
    write(*,'(i0)') status
    if (status == 0) then
        vals=splev(q,at)
        write(*,'(*(es27.18,1x))') vals
    end if
end program fitpack_tck_driver
