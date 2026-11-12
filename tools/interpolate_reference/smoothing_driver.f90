! SPDX-License-Identifier: MIT
! Native driver for SciPy make_smoothing_spline reference comparisons.
program smoothing_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : make_smoothing_spline, ppoly_type
    implicit none
    type(ppoly_type) :: spline
    real(dp), allocatable :: x(:), y(:,:), w(:), q(:), lambdas(:), a(:,:),b(:,:),c(:,:), &
        coeff(:,:,:), val(:)
    real(dp) :: lam
    integer :: n,nc,nq,mode,j,i,status

    read (*,*) n,nc,nq,mode
    if (n<5 .or. nc<1 .or. nq<0) error stop 'invalid dimensions'
    read (*,*) lam
    allocate(x(n),y(n,nc),w(n),q(nq),lambdas(nc))
    read (*,*) x
    read (*,*) w
    do i=1,n
        read (*,*) y(i,:)
    end do
    if (nq>0) read (*,*) q
    if (mode==0) then
        spline=make_smoothing_spline(x,y,w=w,lam=lam,status=status,selected_lam=lambdas)
    else
        spline=make_smoothing_spline(x,y,w=w,status=status,selected_lam=lambdas)
    end if
    print '(i0)',status
    if (status/=0) stop
    do j=1,nc
        write (*,'(es25.16e3,1x)',advance='no') lambdas(j)
    end do
    print *
    a=spline%evaluate(q)
    b=spline%evaluate(q,nu=1)
    c=spline%evaluate(q,nu=2)
    do i=1,nq
        do j=1,nc
            print '(3(es25.16e3,1x))',a(i,j),b(i,j),c(i,j)
        end do
    end do
    coeff=spline%coefficients()
    do i=1,n-1
        do j=1,nc
            print '(4(es25.16e3,1x))',coeff(:,i,j)
        end do
    end do
    val=spline%integrate(x(1),x(n))
    do j=1,nc
        print '(es25.16e3)', val(j)
    end do
end program smoothing_driver
