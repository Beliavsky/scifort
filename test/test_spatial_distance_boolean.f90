! SPDX-License-Identifier: MIT
program test_spatial_distance_boolean
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial, only : point_distance, pdist, cdist, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    character(len=14), parameter :: metrics(5) = [character(len=14) :: &
        'dice', 'rogerstanimoto', 'russellrao', 'sokalsneath', 'yule']
    real(dp), parameter :: expected(5) = [0.5_dp, 2.0_dp/3.0_dp, &
        0.75_dp, 0.8_dp, 1.0_dp]
    real(dp), parameter :: weighted_expected(5) = [2.0_dp/3.0_dp, 16.0_dp/25.0_dp, &
        15.0_dp/17.0_dp, 8.0_dp/9.0_dp, 30.0_dp/29.0_dp]
    real(dp) :: u(4), v(4), w(4), x(3,4), y(2,4), d, z(4), inf
    real(dp), allocatable :: pairwise(:), cross(:,:)
    integer :: i, code

    u = [1.0_dp,0.0_dp,1.0_dp,0.0_dp]
    v = [1.0_dp,1.0_dp,0.0_dp,0.0_dp]
    z = 0.0_dp
    w = [2.0_dp,3.0_dp,5.0_dp,7.0_dp]
    do i=1,size(metrics)
        d = point_distance(u,v,metrics(i),status=code)
        if (code/=spatial_distance_success) error stop 'Boolean scalar status'
        call assert_close('scalar Boolean', d, expected(i))
        call assert_close('weighted Boolean', &
            point_distance(u,v,metrics(i),w=w), weighted_expected(i))
        call assert_close('huge weight scale', &
            point_distance(u,v,metrics(i),w=1.0e300_dp*w), weighted_expected(i))
        call assert_close('tiny weight scale', &
            point_distance(u,v,metrics(i),w=1.0e-300_dp*w), weighted_expected(i))
        call assert_close('nonzero numeric to bool', &
            point_distance([2.0_dp,0.0_dp,-3.0_dp,0.0_dp], &
                           [9.0_dp,4.0_dp,0.0_dp,0.0_dp],metrics(i)),expected(i))
        d = point_distance(u,v,metrics(i),w=z,status=code)
        if (code/=spatial_distance_success) error stop 'zero weights status'
        if (i==5) then
            call assert_close('Yule all-zero weights',d,0.0_dp)
        else
            if (.not.ieee_is_nan(d)) error stop 'zero weights should be undefined'
        end if
        d = point_distance(u,v,metrics(i),w=[-1.0_dp,1.0_dp,1.0_dp,1.0_dp], &
            status=code)
        if (code/=spatial_distance_invalid_input .or. .not.ieee_is_nan(d)) &
            error stop 'negative Boolean weights not rejected'
        d = point_distance(u,v,metrics(i),w=[1.0_dp],status=code)
        if (code/=spatial_distance_invalid_input .or. .not.ieee_is_nan(d)) &
            error stop 'invalid Boolean weights shape not rejected'
        d = point_distance(u,v,metrics(i),p=2.0_dp,status=code)
        if (code/=spatial_distance_invalid_input .or. .not.ieee_is_nan(d)) &
            error stop 'unexpected Minkowski p accepted'
    end do

    ! SciPy's scalar sokalsneath raises on two false vectors, whereas its
    ! bulk pdist/cdist returns NaN. SciFort consistently returns NaN.
    if (.not.ieee_is_nan(point_distance(z,z,'dice'))) error stop 'empty Dice'
    if (.not.ieee_is_nan(point_distance(z,z,'sokalsneath'))) error stop 'empty Sokal-Sneath'
    call assert_close('empty Rogers-Tanimoto',point_distance(z,z,'rogerstanimoto'),0.0_dp)
    call assert_close('empty Russell-Rao',point_distance(z,z,'russellrao'),1.0_dp)
    call assert_close('empty Yule',point_distance(z,z,'yule'),0.0_dp)
    call assert_close('opposite Yule', &
        point_distance([0.0_dp,1.0_dp],[1.0_dp,0.0_dp],'yule'),2.0_dp)
    call assert_close('one-sided mismatch Yule', &
        point_distance([0.0_dp,0.0_dp],[1.0_dp,0.0_dp],'yule'),0.0_dp)
    call assert_close('zero weight ignored', &
        point_distance([1.0_dp,0.0_dp],[1.0_dp,1.0_dp],'russellrao', &
            w=[1.0_dp,0.0_dp]),0.0_dp)

    x(1,:) = u
    x(2,:) = v
    x(3,:) = z
    y(1,:) = v
    y(2,:) = z
    do i=1,size(metrics)
        pairwise=pdist(x,metrics(i),w=w,status=code)
        if (code/=spatial_distance_success .or. size(pairwise)/=3) &
            error stop 'Boolean pdist returned incorrect status or size'
        call assert_close('Boolean pdist row pair',pairwise(1),weighted_expected(i))
        cross=cdist(x,y,metrics(i),w=w,status=code)
        if (code/=spatial_distance_success .or. any(shape(cross)/=[3,2])) &
            error stop 'Boolean cdist returned incorrect status or size'
        call assert_close('Boolean cdist row pair',cross(1,1),weighted_expected(i))
        if (i==1 .or. i==4) then
            if (.not.ieee_is_nan(cross(3,2))) error stop 'undefined bulk Boolean pair'
        end if
    end do
    inf=ieee_value(0.0_dp,ieee_quiet_nan)
    d=point_distance([1.0_dp,inf],[1.0_dp,0.0_dp],'dice',status=code)
    if (code/=spatial_distance_invalid_input .or. .not.ieee_is_nan(d)) &
        error stop 'nonfinite Boolean data not rejected'
    cross=cdist(x,y,'dice',w=[1.0_dp],status=code)
    if (code/=spatial_distance_invalid_input .or. any(shape(cross)/=[0,0])) &
        error stop 'invalid cdist weight shape'
    pairwise=pdist(x,'yule',w=[0.0_dp,1.0_dp],status=code)
    if (code/=spatial_distance_invalid_input .or. size(pairwise)/=0) &
        error stop 'invalid pdist weight shape'
    print '(a)', 'test_spatial_distance_boolean PASS'

contains
    subroutine assert_close(label, actual, reference)
        character(len=*), intent(in) :: label !! identifying name of checked case
        real(dp), intent(in) :: actual !! computed spatial distance
        real(dp), intent(in) :: reference !! mathematically expected distance
        if (.not.(abs(actual-reference)<=5.0e-14_dp*max(1.0_dp,abs(reference)))) then
            print *, label, actual, reference
            error stop 'Boolean distance mismatch'
        end if
    end subroutine assert_close
end program test_spatial_distance_boolean
