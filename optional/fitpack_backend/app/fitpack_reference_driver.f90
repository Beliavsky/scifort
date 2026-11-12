! SPDX-License-Identifier: MIT
! Reproducible numerical reference driver for the GENUINE optional dependency.
program fitpack_reference_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_compat, only : fitpack_tck_type,splev
    use scifort_optional_fitpack, only : univariate_spline
    implicit none
    integer :: n,k,nq,code
    real(dp) :: s,residual
    real(dp), allocatable :: x(:),y(:),w(:),q(:),v(:)
    type(fitpack_tck_type) :: curve
    read(*,*) n,k,nq,s
    allocate(x(n),y(n),w(n),q(nq))
    read(*,*) x
    read(*,*) y
    read(*,*) w
    read(*,*) q
    curve=univariate_spline(x,y,w=w,k=k,s=s,status=code,residual=residual)
    write(*,'(i0)') code
    if (code /= 0) stop
    write(*,'(es27.18)') residual
    write(*,'(3(i0,1x))') size(curve%t),size(curve%c),curve%k
    write(*,'(*(es27.18,1x))') curve%t
    write(*,'(*(es27.18,1x))') curve%c
    v=splev(q,curve)
    write(*,'(*(es27.18,1x))') v
end program fitpack_reference_driver
