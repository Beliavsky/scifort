! SPDX-License-Identifier: MIT
! Used by tools/interpolate_reference/compare_optional_fitpack_bivariate_scipy.py
program fitpack_bivariate_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, &
        bisplev_points, bisplev, bispline_success
    use scifort_optional_fitpack, only : bivariate_smoothing_spline, lsq_bivariate_spline
    implicit none
    integer :: mode, n, kx, ky, ntx, nty, i, status, nq
    real(dp) :: smooth, residual
    real(dp), allocatable :: x(:),y(:),z(:),w(:),tx(:),ty(:),qx(:),qy(:),v(:)
    type(fitpack_bivariate_tck_type) :: tck

    read(*,*) mode,n,kx,ky,smooth,ntx,nty
    allocate(x(n),y(n),z(n),w(n),tx(ntx),ty(nty))
    do i=1,n
        read(*,*) x(i),y(i),z(i),w(i)
    end do
    if (ntx > 0) read(*,*) tx
    if (nty > 0) read(*,*) ty
    if (mode==0) then
        tck = bivariate_smoothing_spline(x,y,z,w=w,k=kx,s=smooth, &
            status=status,residual=residual)
    else
        tck = lsq_bivariate_spline(x,y,z,tx,ty,w=w,kx=kx,ky=ky, &
            status=status,residual=residual)
    end if
    write(*,'(I0,1X,ES25.16)') status,residual
    if (.not.tck%valid()) stop
    write(*,'(I0,1X,I0,1X,I0,1X,I0)') size(tck%tx),size(tck%ty),tck%kx,tck%ky
    do i=1,size(tck%tx)
        write(*,'(ES25.16)') tck%tx(i)
    end do
    do i=1,size(tck%ty)
        write(*,'(ES25.16)') tck%ty(i)
    end do
    do i=1,size(tck%c)
        write(*,'(ES25.16)') tck%c(i)
    end do
    read(*,*) nq
    allocate(qx(nq),qy(nq))
    do i=1,nq
        read(*,*) qx(i),qy(i)
    end do
    v = bisplev_points(qx,qy,tck,status=status)
    if (status /= bispline_success) error stop 'bisplev_points failed'
    do i=1,nq
        write(*,'(ES25.16)') v(i)
    end do
    if (kx>1 .and. ky>1) then
        v = bisplev_points(qx,qy,tck,dx=1,dy=1,status=status)
        if (status /= bispline_success) error stop 'derivative failed'
        do i=1,nq
            write(*,'(ES25.16)') v(i)
        end do
    end if
end program fitpack_bivariate_driver
