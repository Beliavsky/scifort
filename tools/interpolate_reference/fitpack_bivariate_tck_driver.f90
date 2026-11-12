! SPDX-License-Identifier: MIT
! Read externally supplied FITPACK bivariate knots and coefficients, evaluate
! scattered points and mixed derivatives. Protocol used by the SciPy comparator.
program fitpack_bivariate_tck_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, bisplev
    implicit none
    type(fitpack_bivariate_tck_type) :: tck
    integer :: m, n, ny, nc, i, dx, dy, code
    real(dp) :: x, y
    read(*,*) n,ny,tck%kx,tck%ky
    allocate(tck%tx(n),tck%ty(ny))
    nc=(n-tck%kx-1)*(ny-tck%ky-1)
    allocate(tck%c(nc))
    read(*,*) tck%tx
    read(*,*) tck%ty
    read(*,*) tck%c
    read(*,*) m
    do i=1,m
        read(*,*) x,y,dx,dy
        write(*,'(ES27.18,1X,I0)') bisplev(x,y,tck,dx,dy,code),code
    end do
end program fitpack_bivariate_tck_driver
