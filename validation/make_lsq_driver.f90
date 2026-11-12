! SPDX-License-Identifier: MIT
! Harness for validation/compare_make_lsq_spline.py; reads cases on stdin.
program make_lsq_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_bspline, only: make_lsq_spline, bspline_array_interpolator_type
    implicit none
    integer :: ncases, m, deg, nch, nt, j, test, code
    real(dp), allocatable :: x(:), t(:), w(:), y(:,:), c(:,:)
    type(bspline_array_interpolator_type) :: spline

    read(*,*) ncases
    do test=1,ncases
        read(*,*) m,deg,nch,nt
        allocate(x(m),t(nt),w(m),y(m,nch))
        read(*,*) x
        read(*,*) t
        read(*,*) w
        read(*,*) y
        spline = make_lsq_spline(x,y,t,k=deg,w=w,status=code)
        write(*,'(i0)') code
        if (code == 0) then
            c = spline%coefficients()
            do j=1,nch
                write(*,'(*(es26.17e3,1x))') c(:,j)
            end do
        end if
        deallocate(x,t,w,y)
    end do
end program make_lsq_driver
