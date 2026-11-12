! SPDX-License-Identifier: MIT
program test_interpolate_clough_tocher
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use scifort_kinds, only: dp
    use scifort_interpolate, only: clough_tocher_2d_interpolator, &
        clough_tocher_2d_interpolator_type, complex_clough_tocher_2d_interpolator_type, &
        ct_status_success, ct_status_invalid_input, ct_status_degenerate, griddata
    implicit none
    real(dp) :: x(9,2), v(9,2), q(5,2), expected(2), z(3,2), step(2)
    real(dp), allocatable :: out(:,:),rv(:), vv(:)
    complex(dp) :: complex_v(9,2)
    complex(dp), allocatable :: complex_out(:,:),cv(:)
    integer, allocatable :: tri(:,:)
    integer :: i,code,j,k,t1,t2,nedges
    real(dp) :: center(2),normal(2),deriv1(2),deriv2(2),d1,d2
    type(clough_tocher_2d_interpolator_type) :: p
    type(complex_clough_tocher_2d_interpolator_type) :: cp
    x(:,1)=[-1.0_dp,0.9_dp,1.1_dp,-1.1_dp,0.15_dp,-0.4_dp,0.6_dp,-0.8_dp,0.20_dp]
    x(:,2)=[-0.8_dp,-1.1_dp,1.05_dp,1.1_dp,0.23_dp,0.5_dp,-0.65_dp,-0.3_dp,-0.23_dp]
    do i=1,9
        v(i,1)=2.0_dp+3.0_dp*x(i,1)-4.0_dp*x(i,2)
        v(i,2)=1.0_dp-2.0_dp*x(i,1)+0.5_dp*x(i,2)
    end do
    complex_v=cmplx(v,2.0_dp*v,dp)
    q(1:3,:)=x(1:3,:)
    q(4,:)=[0.0_dp,0.0_dp]
    q(5,:)=[5.0_dp,5.0_dp]
    p=clough_tocher_2d_interpolator(x,v,tol=1.0e-12_dp,status=code)
    if (code/=ct_status_success .or. .not.p%valid() .or. .not.p%converged()) error stop 'CT real construction'
    if (p%channels()/=2) error stop 'CT channels'
    out=p%evaluate(q,status=code)
    if (code/=ct_status_success .or. any(shape(out)/=[5,2])) error stop 'CT batch evaluation'
    do i=1,4
        expected=[2.0_dp+3.0_dp*q(i,1)-4.0_dp*q(i,2), &
                  1.0_dp-2.0_dp*q(i,1)+0.5_dp*q(i,2)]
        if (maxval(abs(out(i,:)-expected))>1.0e-10_dp) error stop 'CT affine reproduction'
    end do
    if (.not.all(ieee_is_nan(out(5,:)))) error stop 'CT NaN outside hull'
    if (abs(p%value(q(4,:),channel=2)-expected(2))>1.0e-10_dp) error stop 'CT scalar API'
    rv=p%evaluate(q(1,:))
    if (maxval(abs(rv-v(1,:)))>1.0e-12_dp) error stop 'CT node interpolation'
    cp=clough_tocher_2d_interpolator(x,complex_v,rescale=.true., &
        fill_value=cmplx(-7.0_dp,2.0_dp,dp),tol=1.0e-12_dp,status=code)
    if (code/=ct_status_success .or. .not.cp%valid()) error stop 'CT complex fit'
    complex_out=cp%evaluate(q,status=code)
    if (code/=0) error stop 'CT complex query'
    if (maxval(abs(complex_out(:4,:)-cmplx(out(:4,:),2.0_dp*out(:4,:),dp)))>1.0e-10_dp) &
        error stop 'CT complex affine reproduction'
    if (maxval(abs(complex_out(5,:)-cmplx(-7.0_dp,2.0_dp,dp)))>1.0e-12_dp) &
        error stop 'CT complex fill'
    cv=cp%evaluate(q(1,:))
    if (maxval(abs(cv-complex_v(1,:)))>1.0e-12_dp) error stop 'CT complex scalar query'
    out=griddata(x,v,q,method='cubic',fill_value=-9.0_dp,status=code)
    if (code/=0 .or. any(abs(out(5,:)+9.0_dp)>1.0e-12_dp)) error stop 'CT griddata matrix'
    rv=griddata(x,v(:,1),q,method='cubic',status=code)
    if (code/=0 .or. maxval(abs(rv(:4)-[v(1,1),v(2,1),v(3,1),2.0_dp]))>1.0e-6_dp) &
        error stop 'CT griddata vector'
    cv=griddata(x,complex_v(:,1),q,method='cubic',status=code)
    if (code/=0 .or. maxval(abs(cv(:4)-cmplx(rv(:4),2.0_dp*rv(:4),dp)))>1.0e-8_dp) &
        error stop 'CT complex griddata vector'
    complex_out=griddata(x,complex_v,q,method='cubic',status=code)
    if (code/=0 .or. maxval(abs(complex_out(:4,:)-complex_v(1,1)))<1.0_dp) &
        error stop 'CT complex griddata matrix'
    ! Shared-interior-edge check: directional slopes from both sides agree.
    tri=p%simplices()
    nedges=0
    do i=1,size(tri,2)
        do j=1,3
            t1=tri(mod(j,3)+1,i)
            t2=tri(mod(j+1,3)+1,i)
            if (t1>t2) cycle
            do k=1,size(tri,2)
                if (k==i) cycle
                if (.not.any(tri(:,k)==t1) .or. .not.any(tri(:,k)==t2)) cycle
                center=0.5_dp*(x(t1,:)+x(t2,:))
                step=x(t2,:)-x(t1,:)
                normal=[-step(2),step(1)]
                normal=normal/sqrt(sum(normal**2))
                d1=1.0e-5_dp
                rv=p%evaluate(center)
                vv=p%evaluate(center+d1*normal)
                deriv1=(vv-rv)/d1
                vv=p%evaluate(center-d1*normal)
                deriv2=(rv-vv)/d1
                if (maxval(abs(deriv1-deriv2))>2.0e-3_dp) error stop 'CT C1 shared edge'
                nedges=nedges+1
                exit
            end do
        end do
    end do
    if (nedges<3) error stop 'not enough interior edges for C1 checks'
    p=clough_tocher_2d_interpolator(x(:2,:),v(:2,:),status=code)
    if (code/=ct_status_invalid_input .or. p%valid()) error stop 'CT two sites accepted'
    p=clough_tocher_2d_interpolator(x,v(:8,:),status=code)
    if (code/=ct_status_invalid_input .or. p%valid()) error stop 'CT values mismatch accepted'
    z(:,1)=[0.0_dp,1.0_dp,2.0_dp]
    z(:,2)=[0.0_dp,1.0_dp,2.0_dp]
    p=clough_tocher_2d_interpolator(z,[1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code/=ct_status_degenerate .or. p%valid()) error stop 'CT collinear sites accepted'
    z(:,1)=[0.0_dp,0.0_dp,2.0_dp]
    z(:,2)=[0.0_dp,0.0_dp,1.0_dp]
    p=clough_tocher_2d_interpolator(z,[1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code/=ct_status_invalid_input .or. p%valid()) error stop 'CT duplicate sites accepted'
    print '(a,i0,a)', 'Clough-Tocher real/complex/affine/C1/invalid/griddata PASS (',nedges,' shared edges)'
end program test_interpolate_clough_tocher
