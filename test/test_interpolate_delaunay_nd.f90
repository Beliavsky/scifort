! SPDX-License-Identifier: MIT
! Native, rank-generic Delaunay triangulation in dimensions four and five.
program test_interpolate_delaunay_nd
    use scifort_kinds, only : dp
    use scifort_interpolate, only : linear_nd_interpolator, linear_nd_interpolator_type, &
        complex_linear_nd_interpolator_type, delaunay_nd, &
        scattered_status_success, scattered_status_invalid_input
    implicit none
    type(linear_nd_interpolator_type) :: model
    type(complex_linear_nd_interpolator_type) :: cmodel
    real(dp), allocatable :: x(:,:),y(:,:),q(:,:),computed(:,:)
    complex(dp), allocatable :: z(:,:),cvalues(:,:)
    integer, allocatable :: simplices(:,:)
    integer :: d,n,i,j,k,code
    real(dp) :: expected(2),coeff

    do d=4,5
        n=24
        allocate(x(n,d),y(n,2),z(n,2),q(5,d))
        do j=1,d
            do i=1,n
                x(i,j)=sin(0.47_dp*i+0.73_dp*j)+cos(0.37_dp*i*j)
            end do
        end do
        x(:,1)=x(:,1)+80.0_dp
        x(:,d)=0.11_dp*x(:,d)
        do i=1,n
            y(i,1)=2.0_dp+sum([(0.13_dp*j*x(i,j),j=1,d)])
            y(i,2)=-0.8_dp+sum([(-0.21_dp*j*x(i,j),j=1,d)])
        end do
        z=cmplx(y,0.7_dp*y,dp)
        call delaunay_nd(x,simplices,code)
        if (code/=scattered_status_success .or. size(simplices,2)<n-d) &
            error stop 'native ND Delaunay construction'
        if (size(simplices,1)/=d+1 .or. any(simplices<1) .or. any(simplices>n)) &
            error stop 'native ND simplex indexing'
        model=linear_nd_interpolator(x,y,fill_value=-99.0_dp,status=code)
        if (code/=scattered_status_success .or. .not.model%valid()) &
            error stop 'native ND public real constructor'
        if (any(shape(model%simplices())/=shape(simplices))) &
            error stop 'native ND simplex access'
        q(1,:)=x(1,:)
        q(2,:)=x(n,:)
        q(3,:)=sum(x(simplices(:,1),:),dim=1)/real(d+1,dp)
        q(4,:)=sum(x,dim=1)/real(n,dp)
        q(5,:)=maxval(x,dim=1)+10.0_dp
        computed=model%evaluate(q,status=code)
        if (code/=0) error stop 'native ND real evaluation status'
        do i=1,4
            expected(1)=2.0_dp+sum([(0.13_dp*j*q(i,j),j=1,d)])
            expected(2)=-0.8_dp+sum([(-0.21_dp*j*q(i,j),j=1,d)])
            if (maxval(abs(computed(i,:)-expected))>1.0e-9_dp) &
                error stop 'native ND affine reproduction'
        end do
        if (any(computed(5,:)/=-99.0_dp)) error stop 'native ND fill value'
        cmodel=linear_nd_interpolator(x,z,rescale=.true.,fill_value=cmplx(-9.0_dp,3.0_dp,dp),status=code)
        if (code/=scattered_status_success .or. .not.cmodel%valid()) &
            error stop 'native ND complex constructor'
        cvalues=cmodel%evaluate(q,status=code)
        if (code/=0) error stop 'native ND complex evaluation status'
        if (maxval(abs(cvalues(:4,:)-cmplx(computed(:4,:),0.7_dp*computed(:4,:),dp)))>1.0e-9_dp) &
            error stop 'native ND complex affine reproduction'
        if (any(cvalues(5,:)/=cmplx(-9.0_dp,3.0_dp,dp))) error stop 'native ND complex fill'
        x(n,:)=x(n-1,:)
        model=linear_nd_interpolator(x,y,status=code)
        if (code/=scattered_status_invalid_input .or. model%valid()) &
            error stop 'native ND duplicated site validation'
        call delaunay_nd(x,simplices,code)
        if (code/=scattered_status_invalid_input) error stop 'native ND direct duplicated site rejection'
        if (allocated(simplices)) deallocate(simplices)
        deallocate(x,y,z,q)
    end do
    print '(a)', '4D/5D native Delaunay, affine real/complex, rescale and invalid inputs PASS'
end program
