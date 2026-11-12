! SPDX-License-Identifier: MIT
! Data-driven reference test driver; compare with scipy.interpolate in Python.
program pchip_akima_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : piecewise_cubic_interpolator, &
        pchip_interpolator, akima1d_interpolator, cubic_hermite_spline
    implicit none
    type(piecewise_cubic_interpolator) :: poly
    real(dp), allocatable :: x(:), y(:), q(:), c(:,:), dydx(:)
    real(dp) :: a, b
    integer :: count, trial, n, nq, flavor, ext, i, j, code

    read(*,*) count
    do trial = 1, count
        read(*,*) n, nq, flavor, ext
        allocate(x(n), y(n), q(nq))
        read(*,*) x
        read(*,*) y
        if (flavor == 4) then
            allocate(dydx(n))
            read(*,*) dydx
        end if
        read(*,*) q
        read(*,*) a, b
        select case (flavor)
        case (1)
            if (ext == -1) then
                poly = pchip_interpolator(x, y, status=code)
            else
                poly = pchip_interpolator(x, y, status=code, extrapolate=(ext == 1))
            end if
        case (4)
            if (ext == -1) then
                poly = cubic_hermite_spline(x, y, dydx, status=code)
            else
                poly = cubic_hermite_spline(x, y, dydx, status=code, extrapolate=(ext == 1))
            end if
        case (2,3)
            if (ext == -1) then
                if (flavor == 2) then
                    poly = akima1d_interpolator(x, y, status=code)
                else
                    poly = akima1d_interpolator(x, y, method='makima', status=code)
                end if
            else
                if (flavor == 2) then
                    poly = akima1d_interpolator(x, y, status=code, extrapolate=(ext == 1))
                else
                    poly = akima1d_interpolator(x, y, method='makima', &
                        status=code, extrapolate=(ext == 1))
                end if
            end if
        end select
        print '(a,2(1x,i0))', 'CASE', trial, code
        do i = 1, nq
            print '(5(1x,es25.16e3))', (poly%derivative(q(i),j),j=0,4)
        end do
        print '(1x,es25.16e3)', poly%integrate(a,b)
        c = poly%coefficients()
        do i = 1, n-1
            print '(4(1x,es25.16e3))', c(:,i)
        end do
        deallocate(x,y,q,c)
        if (allocated(dydx)) deallocate(dydx)
    end do
end program pchip_akima_driver
