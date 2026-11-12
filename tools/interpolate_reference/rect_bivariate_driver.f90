! SPDX-License-Identifier: MIT
! Simple stream-based reference driver (not a public API).
program rect_bivariate_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, &
        rect_bivariate_spline, bisplev, bispline_success
    implicit none
    integer :: nx, ny, kx, ky, nq, ni, status, i
    real(dp), allocatable :: x(:), y(:), z(:,:), qx(:), qy(:), bounds(:,:)
    real(dp) :: values(6)
    type(fitpack_bivariate_tck_type) :: fit

    read(*,*) nx,ny,kx,ky,nq,ni
    allocate(x(nx),y(ny),z(nx,ny),qx(nq),qy(nq),bounds(4,ni))
    read(*,*) x
    read(*,*) y
    read(*,*) z
    read(*,*) qx
    read(*,*) qy
    do i=1,ni
        read(*,*) bounds(:,i)
    end do
    fit = rect_bivariate_spline(x,y,z,kx=kx,ky=ky,status=status)
    print '(i0)',status
    if (status /= bispline_success) stop
    print '(i0,1x,i0,1x,i0)',size(fit%tx),size(fit%ty),size(fit%c)
    print '(*(es26.17,1x))',fit%tx
    print '(*(es26.17,1x))',fit%ty
    print '(*(es26.17,1x))',fit%c
    do i=1,nq
        values = 0.0_dp
        values(1) = bisplev(qx(i),qy(i),fit)
        if (kx > 1) values(2) = bisplev(qx(i),qy(i),fit,dx=1)
        if (ky > 1) values(3) = bisplev(qx(i),qy(i),fit,dy=1)
        if (kx > 2) values(4) = bisplev(qx(i),qy(i),fit,dx=2)
        if (ky > 2) values(5) = bisplev(qx(i),qy(i),fit,dy=2)
        if (kx > 1 .and. ky > 1) values(6) = bisplev(qx(i),qy(i),fit,dx=1,dy=1)
        print '(*(es26.17,1x))',values
    end do
    do i=1,ni
        print '(es26.17)',fit%integral(bounds(1,i),bounds(2,i),bounds(3,i),bounds(4,i))
    end do
end program rect_bivariate_driver
