! SPDX-License-Identifier: MIT
! Test harness for SciPy comparison: periodic make_interp_spline, scalar and matrix.
program periodic_bspline_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_bspline, only : bspline_interpolator_type, &
        bspline_array_interpolator_type, make_interp_spline
    implicit none
    type(bspline_interpolator_type) :: scalar
    type(bspline_array_interpolator_type) :: spline
    real(dp), allocatable :: x(:), y(:,:), q(:), bounds(:,:), &
        knots(:), coeff(:,:), values(:,:), integrals(:,:)
    integer :: ncases, testcase, n, channels, degree, nq, ni, extrap, code, other
    integer :: i, j, order

    read(*,*) ncases
    do testcase=1,ncases
        read(*,*) n, channels, degree, nq, ni, extrap
        allocate(x(n),y(n,channels),q(nq),bounds(ni,2))
        read(*,*) x
        read(*,*) y
        read(*,*) q
        read(*,*) bounds(:,1)
        read(*,*) bounds(:,2)
        spline = make_interp_spline(x,y,k=degree,bc_type='periodic', &
            extrapolate=extrap/=0,status=code)
        write(*,'(a,1x,i0,1x,i0)') 'CASE',testcase,code
        if (code==0) then
            scalar = make_interp_spline(x,y(:,1),k=degree,bc_type='periodic', &
                extrapolate=extrap/=0,status=other)
            if (other/=0) error stop 'scalar periodic fit failed'
            knots = spline%knots()
            coeff = spline%coefficients()
            write(*,'(a,1x,i0)') 'T',size(knots)
            write(*,'(*(es26.17e3,1x))') knots
            write(*,'(a,1x,i0,1x,i0)') 'C',size(coeff,1),size(coeff,2)
            do j=1,channels
                write(*,'(*(es26.17e3,1x))') coeff(:,j)
            end do
            if (maxval(abs(scalar%coefficients()-coeff(:,1)))>1.0e-9_dp) &
                error stop 'scalar and matrix periodic coefficients differ'
            do order=0,min(3,degree)
                values = spline%derivative(q,order)
                write(*,'(a,1x,i0,1x,i0)') 'V',order,size(q)
                do j=1,channels
                    write(*,'(*(es26.17e3,1x))') values(:,j)
                end do
            end do
            allocate(integrals(ni,channels))
            do i=1,ni
                integrals(i,:) = spline%integrate(bounds(i,1),bounds(i,2))
            end do
            write(*,'(a,1x,i0)') 'I',ni
            do j=1,channels
                write(*,'(*(es26.17e3,1x))') integrals(:,j)
            end do
        end if
        deallocate(x,y,q,bounds)
        if (allocated(knots)) deallocate(knots)
        if (allocated(coeff)) deallocate(coeff)
        if (allocated(values)) deallocate(values)
        if (allocated(integrals)) deallocate(integrals)
    end do
end program periodic_bspline_driver
