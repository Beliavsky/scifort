! SPDX-License-Identifier: MIT
program test_interpolate_scattered
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use scifort_kinds, only: dp
    use scifort_interpolate, only: nearest_nd_interpolator_type, complex_nearest_nd_interpolator_type, &
        linear_nd_interpolator_type, complex_linear_nd_interpolator_type, &
        nearest_nd_interpolator, linear_nd_interpolator, griddata, &
        scattered_status_success, scattered_status_invalid_input, scattered_status_degenerate
    implicit none
    real(dp) :: xy(8,2), z(8,2), queries(6,2), r, out(2)
    real(dp) :: space(5,4), vals(5,2), xq(4), fill
    complex(dp) :: cz(8,2), cv(5,2)
    real(dp), allocatable :: a(:,:),vector(:)
    complex(dp), allocatable :: ca(:,:),cvector(:)
    integer, allocatable :: tri(:,:)
    type(linear_nd_interpolator_type) :: lin
    type(complex_linear_nd_interpolator_type) :: clin
    type(nearest_nd_interpolator_type) :: nn
    type(complex_nearest_nd_interpolator_type) :: cnn
    integer :: i,code

    xy(:,1)=[-1.0_dp,0.9_dp,1.1_dp,-1.1_dp,0.15_dp,-0.4_dp,0.6_dp,-0.8_dp]
    xy(:,2)=[-0.8_dp,-1.1_dp,1.05_dp,1.1_dp,0.23_dp,0.5_dp,-0.65_dp,-0.3_dp]
    do i=1,8
        z(i,1)=2.0_dp+3.0_dp*xy(i,1)-4.0_dp*xy(i,2)
        z(i,2)=1.0_dp-2.0_dp*xy(i,1)+0.5_dp*xy(i,2)
    end do
    cz=cmplx(z,2.0_dp*z,dp)
    queries(1:3,:)=xy(1:3,:)
    queries(4,:)=[0.0_dp,0.0_dp]
    queries(5,:)=[0.2_dp,0.35_dp]
    queries(6,:)=[4.0_dp,4.0_dp]

    lin=linear_nd_interpolator(xy,z,status=code)
    if (code/=scattered_status_success .or. .not.lin%valid()) error stop 'linear constructor'
    if (lin%dimensions()/=2 .or. lin%channels()/=2) error stop 'linear metadata'
    tri=lin%simplices()
    if (size(tri,1)/=3 .or. size(tri,2)<6) error stop 'linear triangles'
    a=lin%evaluate(queries,status=code)
    if (code/=0 .or. any(shape(a)/=[6,2])) error stop 'linear batch'
    do i=1,5
        out=[2.0_dp+3.0_dp*queries(i,1)-4.0_dp*queries(i,2), &
             1.0_dp-2.0_dp*queries(i,1)+0.5_dp*queries(i,2)]
        if (maxval(abs(a(i,:)-out))>5.0e-13_dp) error stop 'linear affine reproduction'
        if (abs(lin%value(queries(i,:),channel=2)-out(2))>5.0e-13_dp) error stop 'linear value'
    end do
    if (.not.all(ieee_is_nan(a(6,:)))) error stop 'linear default NaN'
    a=griddata(xy,z,queries,method='linear',fill_value=-17.0_dp,status=code)
    if (code/=0 .or. any(abs(a(6,:)+17.0_dp)>1.0e-13_dp)) error stop 'griddata linear fill'
    if (maxval(abs(a(:5,:)-lin%evaluate(queries(:5,:))))>1.0e-13_dp) error stop 'griddata linear'
    vector=griddata(xy,z(:,1),queries,method='linear',fill_value=-9.0_dp,status=code)
    if (code/=0 .or. size(vector)/=6 .or. abs(vector(6)+9.0_dp)>1.0e-13_dp) error stop 'griddata vector'
    clin=linear_nd_interpolator(xy,cz,rescale=.true.,fill_value=cmplx(-1.0_dp,2.0_dp,dp),status=code)
    if (code/=0) error stop 'complex linear fit'
    ca=clin%evaluate(queries,status=code)
    if (code/=0 .or. maxval(abs(ca(:5,:)-cmplx(a(:5,:),2.0_dp*a(:5,:),dp)))>1.0e-12_dp) &
        error stop 'complex linear reproduction'
    if (maxval(abs(ca(6,:)-cmplx(-1.0_dp,2.0_dp,dp)))>1.0e-12_dp) error stop 'complex linear fill'
    ca=griddata(xy,cz,queries,method='linear',fill_value=cmplx(-3.0_dp,2.0_dp,dp),status=code)
    if (code/=0 .or. abs(ca(6,1)-cmplx(-3.0_dp,2.0_dp,dp))>1.0e-14_dp) error stop 'complex griddata'
    cvector=griddata(xy,cz(:,1),queries,method='linear',status=code)
    if (code/=0 .or. maxval(abs(cvector(:5)-cz(1,1)))<0.01_dp) then
        ! Ensure complex vector overload is exercised and returns varied data.
        error stop 'complex vector griddata'
    end if

    do i=1,5
        space(i,:)=[real(i,dp),10.0_dp*real(i,dp),real(i*i,dp),1.0_dp]
        vals(i,:)=[real(i,dp),3.0_dp*real(i,dp)]
    end do
    xq=space(3,:)+[0.01_dp,-0.02_dp,0.04_dp,0.0_dp]
    nn=nearest_nd_interpolator(space,vals,rescale=.true.,status=code)
    if (code/=0 .or. .not.nn%valid()) error stop 'nearest 4D fit'
    if (abs(nn%value(xq,channel=2)-9.0_dp)>1.0e-14_dp) error stop 'nearest 4D'
    a=griddata(space,vals,reshape(xq,[1,4]),method='nearest',rescale=.true.,status=code)
    if (code/=0 .or. a(1,1)/=3.0_dp) error stop 'griddata nearest 4D'
    cv=cmplx(vals,2.0_dp*vals,dp)
    cnn=nearest_nd_interpolator(space,cv,status=code)
    if (code/=0 .or. abs(cnn%value(space(5,:),channel=2)-cmplx(15.0_dp,30.0_dp,dp))>1.0e-14_dp) &
        error stop 'nearest complex 4D'
    lin=linear_nd_interpolator(xy(:3,:),z(:3,:),status=code)
    if (code/=0 .or. .not.lin%valid()) error stop 'minimal triangle'
    lin=linear_nd_interpolator(xy(:3,:),z(:2,:),status=code)
    if (code/=scattered_status_invalid_input .or. lin%valid()) error stop 'mismatched values accepted'
    lin=linear_nd_interpolator(reshape([0.0_dp,1.0_dp,2.0_dp,0.0_dp,1.0_dp,2.0_dp],[3,2]), &
        [1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code/=scattered_status_degenerate .or. lin%valid()) error stop 'collinear points accepted'
    lin=linear_nd_interpolator(reshape([0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp],[3,2]), &
        [1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code/=scattered_status_invalid_input .or. lin%valid()) error stop 'duplicate points accepted'
    a=griddata(xy,z,queries,method='cubic',fill_value=-17.0_dp,status=code)
    if (code/=0 .or. maxval(abs(a(1:5,:)-z(1,1)))<0.01_dp) error stop 'cubic griddata public dispatch'
    if (any(abs(a(6,:)+17.0_dp)>1.0e-13_dp)) error stop 'cubic griddata fill'
    print '(a)', 'NearestND, LinearND, 2D Delaunay, complex, griddata, invalid inputs PASS'
end program test_interpolate_scattered
