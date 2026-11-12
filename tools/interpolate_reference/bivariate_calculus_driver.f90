! SPDX-License-Identifier: MIT
! Reference input/output driver for bivariate TCK calculus and refinement.
program bivariate_calculus_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, &
        rect_bivariate_spline, bisplev
    implicit none
    integer :: nx,ny,kx,ky,nq,dx,dy,axis,m,i,code_fit,code_der,code_ins,code_dins
    real(dp) :: knot, bounds(4)
    real(dp), allocatable :: x(:),y(:),z(:,:),qx(:),qy(:)
    type(fitpack_bivariate_tck_type) :: fit, der, inserted, derinserted

    read(*,*) nx,ny,kx,ky,nq,dx,dy,axis,m,knot
    allocate(x(nx),y(ny),z(nx,ny),qx(nq),qy(nq))
    read(*,*) x
    read(*,*) y
    read(*,*) z
    read(*,*) qx
    read(*,*) qy
    read(*,*) bounds
    fit = rect_bivariate_spline(x,y,z,kx=kx,ky=ky,status=code_fit)
    if (code_fit /= 0) error stop 'invalid reference-grid fit'
    der = fit%partial_derivative(dx,dy,status=code_der)
    inserted = fit%insert_knot(axis,knot,m=m,status=code_ins)
    if (code_ins == 0) then
        derinserted = inserted%partial_derivative(dx,dy,status=code_dins)
    else
        code_dins = 1
    end if
    write(*,'(4(i0,1x))') code_fit,code_der,code_ins,code_dins
    if (code_der /= 0 .or. code_ins /= 0 .or. code_dins /= 0) stop
    call write_tck(fit)
    call write_tck(der)
    call write_tck(inserted)
    call write_tck(derinserted)
    do i = 1, nq
        write(*,'(6(es25.16,1x))') bisplev(qx(i),qy(i),fit), bisplev(qx(i),qy(i),der), &
            bisplev(qx(i),qy(i),inserted),bisplev(qx(i),qy(i),derinserted), &
            bisplev(qx(i),qy(i),fit,dx=dx,dy=dy), &
            bisplev(qx(i),qy(i),inserted,dx=dx,dy=dy)
    end do
    write(*,'(4(es25.16,1x))') fit%integral(bounds(1),bounds(2),bounds(3),bounds(4)), &
        der%integral(bounds(1),bounds(2),bounds(3),bounds(4)), &
        inserted%integral(bounds(1),bounds(2),bounds(3),bounds(4)), &
        derinserted%integral(bounds(1),bounds(2),bounds(3),bounds(4))
contains
    subroutine write_tck(tck)
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! tensor-product TCK
        write(*,'(5(i0,1x))') size(tck%tx),size(tck%ty),size(tck%c),tck%kx,tck%ky
        write(*,'(*(es25.16,1x))') tck%tx
        write(*,'(*(es25.16,1x))') tck%ty
        write(*,'(*(es25.16,1x))') tck%c
    end subroutine write_tck
end program bivariate_calculus_driver
