! SPDX-License-Identifier: MIT
! Compare weighted B-spline least-squares / exact interpolation against SciPy.
program fitpack_lsq_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : lsq_univariate_spline, interpolated_univariate_spline, &
        bspline_array_interpolator_type
    implicit none
    type(bspline_array_interpolator_type) :: spline
    real(dp), allocatable :: x(:), y(:,:), w(:), t(:), q(:), knots(:), coeff(:,:), &
        value(:), d1(:), d2(:), d3(:), integral(:)
    real(dp) :: bbox(2)
    integer :: n, nc, nt, deg, nq, mode, status, i, j

    read(*,*) n, nc, nt, deg, nq, mode
    read(*,*) bbox
    allocate(x(n),y(n,nc),w(n),t(nt),q(nq))
    read(*,*) x
    if (nt>0) read(*,*) t
    read(*,*) w
    do i=1,n
        read(*,*) y(i,:)
    end do
    if (nq>0) read(*,*) q
    if (mode == 0) then
        spline = lsq_univariate_spline(x,y,t,w=w,k=deg,bbox=bbox,status=status)
    else
        spline = interpolated_univariate_spline(x,y,k=deg,status=status)
    end if
    print '(i0)',status
    if (status/=0) stop
    knots=spline%knots()
    coeff=spline%coefficients()
    print '(3(i0,1x))',size(knots),size(coeff,1),size(coeff,2)
    print '(*(es26.17e3,1x))',knots
    do i=1,size(coeff,1)
        print '(*(es26.17e3,1x))',coeff(i,:)
    end do
    do i=1,nq
        value=spline%evaluate(q(i))
        d1=spline%derivative(q(i),nu=1)
        d2=spline%derivative(q(i),nu=2)
        d3=spline%derivative(q(i),nu=3)
        do j=1,nc
            print '(4(es26.17e3,1x))',value(j),d1(j),d2(j),d3(j)
        end do
    end do
    integral=spline%integrate(x(1),x(n))
    print '(*(es26.17e3,1x))',integral
end program fitpack_lsq_driver
