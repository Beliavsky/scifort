! SPDX-License-Identifier: MIT
! Independent 3-D Delaunay/LinearND regression with affine exactness.
program test_interpolate_scattered_3d
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use scifort_kinds, only: dp
    use scifort_interpolate, only: linear_nd_interpolator, linear_nd_interpolator_type, &
        complex_linear_nd_interpolator_type, griddata, scattered_status_success, &
        scattered_status_invalid_input, scattered_status_degenerate
    implicit none
    integer, parameter :: nsite=18, nquery=9
    real(dp) :: x(nsite,3),values(nsite,2),q(nquery,3),r(2)
    real(dp), allocatable :: yy(:,:),z(:)
    complex(dp) :: cz(nsite,2)
    complex(dp), allocatable :: cc(:,:),cv(:)
    integer, allocatable :: s(:,:)
    type(linear_nd_interpolator_type) :: p
    type(complex_linear_nd_interpolator_type) :: pc
    integer :: i,code
    ! Include the exact vertices of the enclosing tetrahedron, as well as
    ! random-looking interior sites, for deterministic boundary coverage.
    x(1,:)=[0.0_dp,0.0_dp,0.0_dp]
    x(2,:)=[2.0_dp,0.0_dp,0.0_dp]
    x(3,:)=[0.0_dp,3.0_dp,0.0_dp]
    x(4,:)=[0.0_dp,0.0_dp,4.0_dp]
    do i=5,nsite
        x(i,:)=[0.05_dp+0.05_dp*real(mod(i*7,11),dp), &
                0.06_dp+0.09_dp*real(mod(i*5,13),dp), &
                0.08_dp+0.11_dp*real(mod(i*9,17),dp)]
        ! Shrink the coordinates if a point is outside the seed tetrahedron.
        if (sum(x(i,:)/[2.0_dp,3.0_dp,4.0_dp])>0.9_dp) x(i,:)=0.7_dp*x(i,:)
    end do
    do i=1,nsite
        values(i,1)=1.2_dp+2.0_dp*x(i,1)-3.0_dp*x(i,2)+0.4_dp*x(i,3)
        values(i,2)=-0.9_dp-1.0_dp*x(i,1)+0.2_dp*x(i,2)+2.1_dp*x(i,3)
    end do
    cz=cmplx(values,2.5_dp*values,dp)
    q(1:4,:)=x(1:4,:)
    q(5,:)=[0.35_dp,0.2_dp,0.4_dp]
    q(6,:)=[0.55_dp,0.5_dp,0.5_dp]
    q(7,:)=0.25_dp*sum(x(1:4,:),dim=1)
    q(8,:)=[1.0_dp,0.0_dp,0.0_dp] ! exact hull edge midpoint
    q(9,:)=[-0.3_dp,0.6_dp,1.0_dp] ! outside hull

    p=linear_nd_interpolator(x,values,fill_value=-42.0_dp,status=code)
    if (code/=scattered_status_success .or. .not.p%valid()) error stop '3D construct'
    if (p%dimensions()/=3 .or. p%channels()/=2) error stop '3D metadata'
    s=p%simplices()
    if (size(s,1)/=4 .or. size(s,2)<10) error stop '3D tetrahedra shape'
    if (any(s<1) .or. any(s>nsite)) error stop '3D tetrahedron vertex indices'
    yy=p%evaluate(q,status=code)
    if (code/=scattered_status_success) error stop '3D evaluate'
    do i=1,nquery-1
        r=[1.2_dp+2.0_dp*q(i,1)-3.0_dp*q(i,2)+0.4_dp*q(i,3), &
           -0.9_dp-q(i,1)+0.2_dp*q(i,2)+2.1_dp*q(i,3)]
        if (maxval(abs(yy(i,:)-r))>2.0e-11_dp) error stop '3D affine reproduction'
        if (abs(p%value(q(i,:),channel=2)-r(2))>2.0e-11_dp) error stop '3D scalar value'
    end do
    if (maxval(abs(yy(9,:)+42.0_dp))>1.0e-12_dp) error stop '3D out-of-hull fill'
    yy=griddata(x,values,q,method='linear',fill_value=-123.0_dp,status=code)
    if (code/=scattered_status_success) error stop '3D griddata matrix'
    if (maxval(abs(yy(9,:)+123.0_dp))>1.0e-12_dp) error stop '3D griddata fill'
    z=griddata(x,values(:,1),q,method='linear',status=code)
    if (code/=0 .or. abs(z(5)-(1.2_dp+2.0_dp*q(5,1)-3.0_dp*q(5,2)+0.4_dp*q(5,3))) &
        >2.0e-11_dp) error stop '3D griddata vector'
    if (.not.ieee_is_nan(z(9))) error stop '3D default NaN'

    pc=linear_nd_interpolator(x,cz,rescale=.true.,fill_value=cmplx(-2.0_dp,3.0_dp,dp),status=code)
    if (code/=0 .or. .not.pc%valid()) error stop '3D complex construct'
    s=pc%simplices()
    if (size(s,1)/=4) error stop '3D complex simplex shape'
    cc=pc%evaluate(q,status=code)
    if (code/=0) error stop '3D complex evaluate'
    do i=1,nquery-1
        r=[1.2_dp+2.0_dp*q(i,1)-3.0_dp*q(i,2)+0.4_dp*q(i,3), &
           -0.9_dp-q(i,1)+0.2_dp*q(i,2)+2.1_dp*q(i,3)]
        if (maxval(abs(cc(i,:)-cmplx(r,2.5_dp*r,dp)))>2.0e-11_dp) &
            error stop '3D complex affine'
    end do
    if (any(abs(cc(9,:)-cmplx(-2.0_dp,3.0_dp,dp))>1.0e-12_dp)) error stop '3D complex fill'
    cv=griddata(x,cz(:,1),q,method='linear',fill_value=cmplx(5.0_dp,4.0_dp,dp),status=code)
    if (code/=0 .or. abs(cv(9)-cmplx(5.0_dp,4.0_dp,dp))>1.0e-12_dp) &
        error stop '3D complex griddata'

    p=linear_nd_interpolator(x(:3,:),values(:3,:),status=code)
    if (code/=scattered_status_invalid_input .or. p%valid()) error stop '3D insufficient sites'
    p=linear_nd_interpolator(x,values(:nsite-1,:),status=code)
    if (code/=scattered_status_invalid_input .or. p%valid()) error stop '3D mismatched values'
    p=linear_nd_interpolator(x(1:3,1:2),values(1:3,:),status=code)
    if (code/=0) error stop '2D preserved'
    p=linear_nd_interpolator(x,values,status=code)
    if (code/=0) error stop '3D preserved'
    x(6,:)=x(5,:)
    p=linear_nd_interpolator(x,values,status=code)
    if (code/=scattered_status_invalid_input .or. p%valid()) error stop '3D duplicate sites'
    x(6,:)=[0.83_dp,0.17_dp,0.22_dp]
    x(4,1:2)=[0.21_dp,0.32_dp]
    x(:,3)=0.0_dp
    p=linear_nd_interpolator(x,values,status=code)
    if (code/=scattered_status_degenerate .or. p%valid()) error stop '3D coplanar sites'
    print '(a)', '3D Delaunay, LinearND, griddata, complex, affine and error inputs PASS'
end program test_interpolate_scattered_3d
