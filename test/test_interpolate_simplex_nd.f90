! SPDX-License-Identifier: MIT
! Dimension-generic supplied-simplex LinearND interpolation and invalid inputs.
program test_interpolate_simplex_nd
    use scifort_kinds, only: dp
    use scifort_interpolate, only: linear_nd_interpolator, linear_nd_interpolator_type, &
        complex_linear_nd_interpolator_type, scattered_status_invalid_input, &
        scattered_status_degenerate, scattered_status_success
    implicit none
    type(linear_nd_interpolator_type) :: p
    type(complex_linear_nd_interpolator_type) :: pc
    real(dp), allocatable :: x(:,:), y(:,:), q(:,:), r(:,:), v(:,:)
    complex(dp), allocatable :: yc(:,:), vc(:,:)
    integer, allocatable :: simplex(:,:), found(:,:)
    integer :: d, k, j, code

    do d=2,10
        allocate(x(d+1,d),y(d+1,3),yc(d+1,3),q(4,d),simplex(d+1,1))
        x=0.0_dp
        do j=1,d
            x(j+1,j)=1.0_dp
        end do
        do j=1,d+1
            y(j,1)=0.4_dp+dot_product(x(j,:),[(0.25_dp*real(k,dp),k=1,d)])
            y(j,2)=1.2_dp-dot_product(x(j,:),[(0.13_dp*real(k,dp),k=1,d)])
            y(j,3)=3.0_dp+sum(x(j,:))
        end do
        yc=cmplx(y,-2.0_dp*y,dp)
        simplex(:,1)=[(j,j=1,d+1)]
        q=0.0_dp
        q(1,:)=0.0_dp
        q(2,:)=0.2_dp/real(d,dp)
        q(3,:)=0.7_dp/real(d,dp)
        q(4,:)=0.5_dp/real(d,dp)
        q(4,1)=1.5_dp   ! outside simplex
        p=linear_nd_interpolator(x,y,simplices=simplex,fill_value=-42.0_dp,status=code)
        if (code/=scattered_status_success .or. .not.p%valid()) error stop 'rank-general construct'
        if (p%dimensions()/=d .or. p%channels()/=3) error stop 'rank-general metadata'
        found=p%simplices()
        if (any(shape(found)/=shape(simplex))) error stop 'rank-general simplex shape'
        if (any(found/=simplex)) error stop 'rank-general simplex indices'
        v=p%evaluate(q,status=code)
        if (code/=0) error stop 'rank-general evaluation status'
        do j=1,3
            if (abs(v(j,1)-0.4_dp-dot_product(q(j,:),[(0.25_dp*real(k,dp),k=1,d)]))>3.0e-12_dp) &
                error stop 'rank-general affine channel 1'
            if (abs(v(j,2)-1.2_dp+dot_product(q(j,:),[(0.13_dp*real(k,dp),k=1,d)]))>3.0e-12_dp) &
                error stop 'rank-general affine channel 2'
            if (abs(v(j,3)-3.0_dp-sum(q(j,:)))>3.0e-12_dp) error stop 'rank-general affine channel 3'
        end do
        if (any(v(4,:)/=-42.0_dp)) error stop 'rank-general fill'
        pc=linear_nd_interpolator(x,yc,simplices=simplex,fill_value=cmplx(-4.0_dp,5.0_dp,dp), &
                                  rescale=.true.,status=code)
        if (code/=0 .or. .not.pc%valid()) error stop 'rank-general complex construct'
        vc=pc%evaluate(q,status=code)
        if (code/=0) error stop 'rank-general complex evaluation'
        if (maxval(abs(vc(:3,:)-cmplx(v(:3,:),-2.0_dp*v(:3,:),dp)))>3.0e-12_dp) &
            error stop 'rank-general complex channels'
        if (any(vc(4,:)/=cmplx(-4.0_dp,5.0_dp,dp))) error stop 'rank-general complex fill'
        if (d==6) then
            p=linear_nd_interpolator(x,y,status=code)
            if (code/=scattered_status_invalid_input .or. p%valid()) &
                error stop 'missing 6D tessellation should fail'
        end if
        if (d==4) then
            p=linear_nd_interpolator(x,y,status=code)
            if (code/=scattered_status_success .or. .not.p%valid()) &
                error stop 'native 4D triangulation should succeed'
            simplex(1,1)=simplex(2,1)
            p=linear_nd_interpolator(x,y,simplices=simplex,status=code)
            if (code/=scattered_status_degenerate .or. p%valid()) error stop 'duplicate simplex vertex'
            simplex(1,1)=1
            simplex(d+1,1)=d+2
            p=linear_nd_interpolator(x,y,simplices=simplex,status=code)
            if (code/=scattered_status_invalid_input .or. p%valid()) error stop 'vertex index range'
            simplex(d+1,1)=d+1
            x(d+1,:)=0.5_dp*x(2,:)+0.5_dp*x(3,:)
            p=linear_nd_interpolator(x,y,simplices=simplex,status=code)
            if (code/=scattered_status_degenerate .or. p%valid()) error stop 'singular simplex'
        end if
        deallocate(x,y,yc,q,simplex)
    end do
    print '(a)', '2D-10D external simplices, affine real/complex interpolation and errors PASS'
end program test_interpolate_simplex_nd
