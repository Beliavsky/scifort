! SPDX-License-Identifier: MIT
program cubic_spline_multi_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : cubic_spline, piecewise_cubic_array_interpolator
    implicit none
    type(piecewise_cubic_array_interpolator) :: poly
    real(dp), allocatable :: x(:), y(:,:), q(:), left(:), right(:), v(:,:), c(:,:,:), intg(:)
    real(dp) :: a, b
    integer :: nt, trial, n, m, nq, flavor, ext, code, i, j, k, nu
    read(*,*) nt
    do trial = 1,nt
        read(*,*) n,m,nq,flavor,ext
        allocate(x(n),y(n,m),q(nq),left(m),right(m))
        read(*,*) x
        do k = 1,m
            read(*,*) y(:,k)
        end do
        read(*,*) q
        read(*,*) a,b
        read(*,*) left
        read(*,*) right
        select case(flavor)
        case(0)
            if (ext == -1) then
                poly = cubic_spline(x,y,status=code)
            else
                poly = cubic_spline(x,y,status=code,extrapolate=(ext==1))
            end if
        case(1)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_type='natural',status=code)
            else
                poly = cubic_spline(x,y,bc_type='natural',status=code,extrapolate=(ext==1))
            end if
        case(2)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_type='clamped',status=code)
            else
                poly = cubic_spline(x,y,bc_type='clamped',status=code,extrapolate=(ext==1))
            end if
        case(3)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_type='periodic',status=code)
            else
                poly = cubic_spline(x,y,bc_type='periodic',status=code,extrapolate=(ext==1))
            end if
        case(4)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='first',bc_right='second', &
                    left_values=left,right_values=right,status=code)
            else
                poly = cubic_spline(x,y,bc_left='first',bc_right='second', &
                    left_values=left,right_values=right,status=code,extrapolate=(ext==1))
            end if
        case(5)
            if (ext == -1) then
                poly = cubic_spline(x,y,bc_left='second',bc_right='first', &
                    left_values=left,right_values=right,status=code)
            else
                poly = cubic_spline(x,y,bc_left='second',bc_right='first', &
                    left_values=left,right_values=right,status=code,extrapolate=(ext==1))
            end if
        end select
        print '(a,2(1x,i0))', 'CASE', trial, code
        if (code /= 0) error stop 'Unexpected matrix fit error'
        do nu = 0,4
            v = poly%derivative(q,order=nu)
            do i = 1,nq
                print '(*(1x,es25.16e3))', (v(i,k),k=1,m)
            end do
        end do
        intg = poly%integrate(a,b)
        print '(*(1x,es25.16e3))', intg
        c = poly%coefficients()
        do k = 1,m
            do j = 1,n-1
                print '(4(1x,es25.16e3))', c(:,j,k)
            end do
        end do
        deallocate(x,y,q,left,right,v,c,intg)
    end do
end program cubic_spline_multi_driver
