! SPDX-License-Identifier: MIT
! Native reference driver for independent SciPy CubicSpline comparisons.
program cubic_spline_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : cubic_spline, piecewise_cubic_interpolator
    implicit none
    type(piecewise_cubic_interpolator) :: poly
    real(dp), allocatable :: x(:), y(:), q(:), c(:,:)
    real(dp) :: a, b, lv, rv
    integer :: trial, count, n, nq, flavor, ext, i, j, code

    read(*,*) count
    do trial = 1, count
        read(*,*) n, nq, flavor, ext
        allocate(x(n), y(n), q(nq))
        read(*,*) x
        read(*,*) y
        read(*,*) q
        read(*,*) a, b, lv, rv
        select case (flavor)
        case (0)
            if (ext == -1) then
                poly = cubic_spline(x,y,status=code)
            else
                poly = cubic_spline(x,y,status=code,extrapolate=(ext == 1))
            end if
        case (1)
            if (ext == -1) then
                poly = cubic_spline(x,y,'natural',status=code)
            else
                poly = cubic_spline(x,y,'natural',status=code,extrapolate=(ext == 1))
            end if
        case (2)
            if (ext == -1) then
                poly = cubic_spline(x,y,'clamped',status=code)
            else
                poly = cubic_spline(x,y,'clamped',status=code,extrapolate=(ext == 1))
            end if
        case (3)
            if (ext == -1) then
                poly = cubic_spline(x,y,'periodic',status=code)
            else
                poly = cubic_spline(x,y,'periodic',status=code,extrapolate=(ext == 1))
            end if
        case (4)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='first',bc_right='natural', &
                    left_value=lv,status=code)
            else
                poly = cubic_spline(x,y,bc_left='first',bc_right='natural', &
                    left_value=lv,status=code,extrapolate=(ext == 1))
            end if
        case (5)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_right='first',right_value=rv,status=code)
            else
                poly = cubic_spline(x,y,bc_right='first',right_value=rv,status=code,extrapolate=(ext == 1))
            end if
        case (6)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='second',bc_right='second', &
                    left_value=lv,right_value=rv,status=code)
            else
                poly = cubic_spline(x,y,bc_left='second',bc_right='second', &
                    left_value=lv,right_value=rv,status=code,extrapolate=(ext == 1))
            end if
        case (7)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='natural',bc_right='clamped',status=code)
            else
                poly = cubic_spline(x,y,bc_left='natural',bc_right='clamped', &
                    status=code,extrapolate=(ext == 1))
            end if
        case (8)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='clamped',bc_right='natural',status=code)
            else
                poly = cubic_spline(x,y,bc_left='clamped',bc_right='natural', &
                    status=code,extrapolate=(ext == 1))
            end if
        case (9)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='first',bc_right='first', &
                    left_value=lv,right_value=rv,status=code)
            else
                poly = cubic_spline(x,y,bc_left='first',bc_right='first', &
                    left_value=lv,right_value=rv,status=code,extrapolate=(ext == 1))
            end if
        case (10)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='second',left_value=lv,status=code)
            else
                poly = cubic_spline(x,y,bc_left='second',left_value=lv, &
                    status=code,extrapolate=(ext == 1))
            end if
        end select
        print '(a,2(1x,i0))', 'CASE', trial, code
        if (code /= 0) error stop 'Unexpected fitting error'
        do i = 1, nq
            print '(5(1x,es25.16e3))', (poly%derivative(q(i),j),j=0,4)
        end do
        print '(1x,es25.16e3)', poly%integrate(a,b)
        c = poly%coefficients()
        do i = 1,n-1
            print '(4(1x,es25.16e3))', c(:,i)
        end do
        deallocate(x,y,q,c)
    end do
end program cubic_spline_driver
