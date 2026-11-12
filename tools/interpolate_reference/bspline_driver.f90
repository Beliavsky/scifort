! SPDX-License-Identifier: MIT
! Read machine-generated interpolation cases and write results for SciPy comparison.
program bspline_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: bspline_interpolator_type, bspline_array_interpolator_type, &
        bspline, make_interp_spline, bspline_design_matrix, bspline_status_success
    implicit none
    integer :: ncases, case_id, mode, n, m, k, bc_id, extrap_id, periodic_id, nt, nq, ni
    integer :: mi, code, j, code2
    real(dp) :: insert_x
    real(dp), allocatable :: x(:), c(:,:), t(:), q(:), a(:), b(:), values(:,:), temp(:,:), dm(:,:)
    type(bspline_array_interpolator_type) :: s, si
    type(bspline_interpolator_type) :: single, deriv, anti
    character(len=16) :: bc
    read(*,*) ncases
    do case_id=1,ncases
        read(*,*) mode,n,m,k,bc_id,extrap_id,periodic_id,nt,nq,ni,mi
        allocate(x(n),c(n,m),t(nt),q(nq),a(ni),b(ni))
        read(*,*) x
        read(*,*) c
        if (nt > 0) read(*,*) t
        read(*,*) q
        read(*,*) a
        read(*,*) b
        read(*,*) insert_x
        select case (bc_id)
        case (1)
            bc = 'natural'
        case (2)
            bc = 'clamped'
        case default
            bc = 'not-a-knot'
        end select
        if (mode == 1) then
            if (nt > 0) then
                s = make_interp_spline(x,c,k=k,bc_type=trim(bc),status=code, &
                    extrapolate=extrap_id==1,knots=t)
            else
                s = make_interp_spline(x,c,k=k,bc_type=trim(bc),status=code, &
                    extrapolate=extrap_id==1)
            end if
        else
            s = bspline(t,c,k,status=code,extrapolate=extrap_id==1,periodic=periodic_id==1)
        end if
        print '(a,i0,1x,i0)', 'CASE ',case_id,code
        if (code == bspline_status_success) then
            t = s%knots()
            print '(a,i0)', 'T ', size(t)
            call print_real(t)
            c = s%coefficients()
            print '(a,i0,1x,i0)', 'C ',size(c,1),size(c,2)
            do j=1,m
                call print_real(c(:,j))
            end do
            print '(a,i0,1x,i0)', 'V ',nq,m
            values = s%evaluate(q)
            do j=1,m
                call print_real(values(:,j))
            end do
            values = s%derivative(q,1)
            do j=1,m
                call print_real(values(:,j))
            end do
            values = s%derivative(q,2)
            do j=1,m
                call print_real(values(:,j))
            end do
            print '(a,i0,1x,i0)', 'I ',ni,m
            allocate(temp(ni,m))
            do j=1,ni
                temp(j,:) = s%integrate(a(j),b(j))
            end do
            do j=1,m
                call print_real(temp(:,j))
            end do
            deallocate(temp)
            if (mi >= 0) then
                si = s%insert_knot(insert_x,m=mi,status=code2)
                print '(a,i0)', 'INS ',code2
                if (code2 == 0) then
                    t = si%knots()
                    c = si%coefficients()
                    print '(a,i0)', 'IT ',size(t)
                    call print_real(t)
                    print '(a,i0,1x,i0)', 'IC ',size(c,1),size(c,2)
                    do j=1,m
                        call print_real(c(:,j))
                    end do
                    values = si%evaluate(q)
                    do j=1,m
                        call print_real(values(:,j))
                    end do
                end if
            else
                print '(a,i0)', 'INS ', -1
            end if
            t = s%knots()
            dm = bspline_design_matrix(q,t,k,status=code2,extrapolate=.true.)
            print '(a,i0)', 'DM ',code2
            if (code2 == 0) then
                do j=1,nq
                    call print_real(dm(j,:))
                end do
            end if
            if (k >= 1) then
                c = s%coefficients()
                single = bspline(t,c(:,1),k,extrapolate=extrap_id==1,periodic=periodic_id==1)
                deriv = single%derivative_spline(status=code2)
                print '(a,i0)', 'DER ',code2
                if (code2 == 0) call print_real(deriv%evaluate(q))
            else
                print '(a,i0)', 'DER ', -1
            end if
            c = s%coefficients()
            single = bspline(t,c(:,1),k,extrapolate=extrap_id==1,periodic=periodic_id==1)
            anti = single%antiderivative(status=code2)
            print '(a,i0)', 'ANT ',code2
            if (code2 == 0) call print_real(anti%evaluate(q))
        end if
        deallocate(x,c,t,q,a,b)
    end do
contains
    subroutine print_real(arr)
        real(dp), intent(in) :: arr(:) !! values to print
        write(*,'(*(es24.16,1x))') arr
    end subroutine print_real
end program bspline_driver
