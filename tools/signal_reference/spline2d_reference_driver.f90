! SPDX-License-Identifier: MIT
! Text protocol for independent SciPy spline2d comparison.
program spline2d_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline2d, only : cspline2d, qspline2d, spline_filter
    implicit none
    integer :: t, cases, m, n, kind_code, is_complex, i, j, status
    real(dp) :: lam, precision, xr, xi
    real(dp), allocatable :: x(:,:), y(:,:)
    complex(dp), allocatable :: xc(:,:), yc(:,:)

    read(*,*) cases
    do t = 1, cases
        read(*,*) m, n, kind_code, is_complex, lam, precision
        allocate(x(m,n), xc(m,n))
        do j = 1, n
            do i = 1, m
                read(*,*) xr, xi
                xc(i,j) = cmplx(xr, xi, dp)
                x(i,j) = xr
            end do
        end do
        if (is_complex == 1) then
            select case (kind_code)
            case (2)
                yc = qspline2d(xc, lamb=lam, precision=precision, status=status)
            case (3)
                yc = cspline2d(xc, lamb=lam, precision=precision, status=status)
            case (4)
                yc = spline_filter(xc, lmbda=lam, status=status)
            end select
        else
            select case (kind_code)
            case (2)
                y = qspline2d(x, lamb=lam, precision=precision, status=status)
            case (3)
                y = cspline2d(x, lamb=lam, precision=precision, status=status)
            case (4)
                y = spline_filter(x, lmbda=lam, status=status)
            end select
            yc = cmplx(y, 0.0_dp, dp)
        end if
        write(*,'(3(i0,1x))') status, size(yc,1), size(yc,2)
        if (status == 0) then
            do j = 1, n
                do i = 1, m
                    write(*,'(2(es25.17e3,1x))') real(yc(i,j),dp), aimag(yc(i,j))
                end do
            end do
        end if
        deallocate(x,xc)
        if (allocated(y)) deallocate(y)
        if (allocated(yc)) deallocate(yc)
    end do
end program spline2d_reference_driver
