! SPDX-License-Identifier: MIT
! Line-oriented Fortran driver for independently generated SciPy root cases.
program ppoly_roots_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly_type, ppoly
    implicit none
    type(ppoly_type) :: p
    real(dp), allocatable :: x(:), c(:,:,:), roots(:)
    real(dp) :: target
    integer :: cases, i, k, nseg, nch, ext_flag, periodic_flag, disc_flag, j, code

    read(*,*) cases
    do i = 1,cases
        read(*,*) k, nseg, nch, ext_flag, periodic_flag, disc_flag, target
        allocate(x(nseg+1),c(k+1,nseg,nch))
        read(*,*) x
        read(*,*) c
        p = ppoly(c,x,extrapolate=ext_flag/=0,periodic=periodic_flag/=0,status=code)
        if (code/=0) error stop 'PPoly driver: invalid input'
        write(*,'(a,1x,i0)') 'CASE',i
        do j = 1,nch
            roots = p%solve(y=target,discontinuity=disc_flag/=0,channel=j,status=code)
            if (code/=0) error stop 'PPoly driver: solve failure'
            write(*,'(a,1x,i0,1x,i0)') 'ROOTS',j,size(roots)
            write(*,'(*(1x,es26.17e3))') roots
        end do
        deallocate(x,c)
    end do
end program ppoly_roots_driver
