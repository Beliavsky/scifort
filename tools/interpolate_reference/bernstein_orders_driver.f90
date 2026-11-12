! SPDX-License-Identifier: MIT
! Batched independent reference tests of BPoly.from_derivatives(orders=...).
program bernstein_orders_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : bpoly_from_derivatives, bpoly_type
    implicit none
    type(bpoly_type) :: spline
    integer :: nc, testcase, n, nderiv, nch, nq, option, order, code, i, der
    integer, allocatable :: counts(:), orders(:)
    real(dp), allocatable :: x(:), y(:,:,:), q(:), coeff(:,:,:), vals(:)

    read(*,*) nc
    do testcase=1,nc
        read(*,*) n, nderiv, nch, nq, option, order
        allocate(x(n),counts(n),orders(n-1),y(n,nderiv,nch),q(nq))
        read(*,*) x
        read(*,*) counts
        read(*,*) orders
        read(*,*) y
        read(*,*) q
        if (option==0) then
            spline=bpoly_from_derivatives(x,y,counts,orders=orders,status=code)
        else
            spline=bpoly_from_derivatives(x,y,counts,order=order,status=code)
        end if
        if (code/=0 .or. .not. spline%valid()) error stop 'Hermite orders construction failed'
        coeff=spline%coefficients()
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',testcase,size(coeff,1),size(coeff,2)
        write(*,'(*(1x,es25.17))') reshape(coeff,[size(coeff)])
        do der=0,2
            do i=1,nq
                vals=spline%evaluate(q(i),nu=der)
                write(*,'(*(1x,es25.17))') vals
            end do
        end do
        deallocate(x,counts,orders,y,q,coeff,vals)
    end do
end program bernstein_orders_driver
