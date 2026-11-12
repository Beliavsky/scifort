! SPDX-License-Identifier: MIT
program test_spatial_distance
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial, only : point_distance, cdist, pdist, squareform, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    character(len=12), parameter :: metrics(*) = [character(len=12) :: &
        'euclidean', 'sqeuclidean', 'cityblock', 'chebyshev', 'minkowski', &
        'cosine', 'correlation', 'hamming', 'jaccard', 'canberra', 'braycurtis']
    real(dp) :: x(4,4), y(3,4), u(4), v(4), w(4), got, expected
    real(dp), allocatable :: d(:,:), dv(:), dm(:,:), invalid_mat(:,:)
    integer :: code, i, j, k, m

    u = [0.0_dp,2.0_dp,0.0_dp,1.0_dp]
    v = [1.0_dp,0.0_dp,0.0_dp,3.0_dp]
    w = [1.0_dp,0.0_dp,2.0_dp,3.0_dp]
    call assert_close('Euclidean',point_distance(u,v),3.0_dp,1.0e-14_dp)
    call assert_close('squared Euclidean',point_distance(u,v,'sqeuclidean'),9.0_dp,1.0e-14_dp)
    call assert_close('Manhattan',point_distance(u,v,'cityblock'),5.0_dp,1.0e-14_dp)
    call assert_close('Chebyshev',point_distance(u,v,'chebyshev'),2.0_dp,1.0e-14_dp)
    call assert_close('Minkowski p=3',point_distance(u,v,'minkowski',p=3.0_dp), &
        17.0_dp**(1.0_dp/3.0_dp),1.0e-14_dp)
    call assert_close('weighted Euclidean',point_distance(u,v,'euclidean',w=w), &
        sqrt(13.0_dp),1.0e-14_dp)
    call assert_close('weighted Chebyshev',point_distance(u,v,'chebyshev',w=w),2.0_dp,1.0e-14_dp)
    call assert_close('Hamming',point_distance(u,v,'hamming'),0.75_dp,1.0e-14_dp)
    call assert_close('Jaccard',point_distance(u,v,'jaccard'),2.0_dp/3.0_dp,1.0e-14_dp)
    call assert_close('Canberra',point_distance(u,v,'canberra'), &
        1.0_dp+1.0_dp+0.5_dp,1.0e-14_dp)
    call assert_close('Bray-Curtis',point_distance(u,v,'braycurtis'),5.0_dp/7.0_dp,1.0e-14_dp)
    ! Avoid intermediate norm-product overflow, which can affect SciPy's
    ! implementation even though the mathematical cosine distance is finite.
    got=point_distance([1.0e200_dp,2.0e200_dp], [2.0e200_dp,-1.0e200_dp], 'cosine')
    call assert_close('extreme cosine',got,1.0_dp,1.0e-14_dp)
    got=point_distance([1.0e200_dp,1.0e200_dp], [0.0_dp,0.0_dp])
    call assert_close('extreme Euclidean',got/1.0e200_dp,sqrt(2.0_dp),1.0e-14_dp)
    got=point_distance([1.0e-200_dp,1.0e-200_dp],[0.0_dp,0.0_dp])
    call assert_close('small Euclidean',got/1.0e-200_dp,sqrt(2.0_dp),1.0e-14_dp)
    got=point_distance([1.0_dp,2.0_dp],[2.0_dp,1.0_dp],'correlation',w=[0.0_dp,1.0_dp])
    if (.not.ieee_is_nan(got)) error stop 'single support correlation'
    if (.not. ieee_is_nan(point_distance(u,u*0.0_dp,'cosine'))) error stop 'cosine zero norm'
    if (.not. ieee_is_nan(point_distance(u*0.0_dp,v*0.0_dp,'braycurtis'))) &
        error stop 'Bray-Curtis zero denominator'
    if (.not. ieee_is_nan(point_distance([1.0_dp,1.0_dp],[2.0_dp,2.0_dp],'correlation'))) &
        error stop 'correlation constant'
    call assert_close('Jaccard all-zero',point_distance(u*0.0_dp,u*0.0_dp,'jaccard'), &
        0.0_dp,1.0e-14_dp)
    call assert_close('Jaccard boolean conversion',point_distance([2.0_dp],[3.0_dp],'jaccard'), &
        0.0_dp,1.0e-14_dp)
    got=point_distance([1.0_dp],[-1.0_dp],'braycurtis')
    if (.not. got>huge(0.0_dp)) error stop 'Bray-Curtis positive infinity'

    x(1,:)=[1.0_dp,0.0_dp,2.0_dp,3.0_dp]
    x(2,:)=[0.0_dp,2.0_dp,1.0_dp,1.0_dp]
    x(3,:)=[3.0_dp,2.0_dp,1.0_dp,0.0_dp]
    x(4,:)=[2.0_dp,3.0_dp,0.0_dp,1.0_dp]
    y(1,:)=[3.0_dp,1.0_dp,3.0_dp,0.0_dp]
    y(2,:)=[1.0_dp,2.0_dp,0.0_dp,2.0_dp]
    y(3,:)=[0.0_dp,1.0_dp,3.0_dp,2.0_dp]
    do m=1,size(metrics)
        if (trim(metrics(m)) == 'minkowski') then
            d=cdist(x,y,trim(metrics(m)),p=3.4_dp,w=w,status=code)
            dv=pdist(x,trim(metrics(m)),p=3.4_dp,w=w,status=k)
        else
            d=cdist(x,y,trim(metrics(m)),w=w,status=code)
            dv=pdist(x,trim(metrics(m)),w=w,status=k)
        end if
        if (code/=spatial_distance_success .or. k/=spatial_distance_success) &
            error stop 'matrix status'
        if (any(shape(d)/=[4,3]) .or. size(dv)/=6) error stop 'distance shapes'
        do i=1,4
            do j=1,3
                if (trim(metrics(m)) == 'minkowski') then
                    expected=point_distance(x(i,:),y(j,:),trim(metrics(m)),p=3.4_dp,w=w)
                else
                    expected=point_distance(x(i,:),y(j,:),trim(metrics(m)),w=w)
                end if
                call assert_close('cdist/point agreement',d(i,j),expected,1.0e-13_dp)
            end do
        end do
        k=0
        do i=1,3
            do j=i+1,4
                k=k+1
                if (trim(metrics(m)) == 'minkowski') then
                    expected=point_distance(x(i,:),x(j,:),trim(metrics(m)),p=3.4_dp,w=w)
                else
                    expected=point_distance(x(i,:),x(j,:),trim(metrics(m)),w=w)
                end if
                call assert_close('pdist/point agreement',dv(k),expected,1.0e-13_dp)
            end do
        end do
        dm=squareform(dv,status=code)
        if (code/=0 .or. any(shape(dm)/=[4,4])) error stop 'squareform vector to matrix'
        if (any(dm(1:4:1,1:4:1)<0.0_dp)) error stop 'negative metric'
    end do

    dv=pdist(x)
    dm=squareform(dv,status=code)
    if (code/=0) error stop 'squareform distance matrix'
    block
        real(dp), allocatable :: back(:)
        back=squareform(dm,status=code)
        if (code/=0 .or. size(back)/=size(dv)) error stop 'squareform roundtrip shape'
        if (any(abs(back-dv)>1.0e-14_dp)) error stop 'squareform roundtrip value'
    end block
    invalid_mat=squareform([1.0_dp,2.0_dp],status=code)
    if (code/=spatial_distance_invalid_input .or. any(shape(invalid_mat)/=[0,0])) &
        error stop 'invalid vector length accepted'
    block
        real(dp), allocatable :: back(:)
        dm(1,1)=1.0_dp
        back=squareform(dm,status=code)
        if (code/=spatial_distance_invalid_input .or. size(back)/=0) &
            error stop 'nonzero diagonal accepted'
        back=squareform(dm,checks=.false.,status=code)
        if (code/=0 .or. size(back)/=6) error stop 'disabled squareform checks'
        dm(1,1)=0.0_dp
        dm(2,1)=dm(2,1)+1.0_dp
        back=squareform(dm,status=code)
        if (code/=spatial_distance_invalid_input) error stop 'asymmetry accepted'
    end block
    call expect_invalid('mismatched vector length',point_distance(u,[1.0_dp],status=code),code)
    got=point_distance(u,v,metric='does_not_exist',status=code)
    call expect_invalid('invalid metric',got,code)
    got=point_distance(u,v,metric='minkowski',p=0.0_dp,status=code)
    call expect_invalid('invalid exponent',got,code)
    got=point_distance(u,v,metric='euclidean',p=3.0_dp,status=code)
    call expect_invalid('inapplicable exponent',got,code)
    got=point_distance(u,v,w=[-1.0_dp,1.0_dp,1.0_dp,1.0_dp],status=code)
    call expect_invalid('negative weights',got,code)
    got=point_distance(u,v,w=[1.0_dp],status=code)
    call expect_invalid('weight length',got,code)
    got=point_distance(u,ieee_value(0.0_dp,ieee_quiet_nan)*v,status=code)
    call expect_invalid('nonfinite vector',got,code)
    d=cdist(x,y(:,1:2),status=code)
    if (code/=1 .or. any(shape(d)/=[0,0])) error stop 'mismatched columns accepted'
    dv=pdist(x(:,1:0),status=code)
    if (code/=1 .or. size(dv)/=0) error stop 'zero-dimensional points accepted'
    d=cdist(x(1:0,:),y,status=code)
    if (code/=0 .or. any(shape(d)/=[0,3])) error stop 'empty-row shape'
    dm=squareform([real(dp)::],status=code)
    if (code/=0 .or. any(shape(dm)/=[1,1]) .or. dm(1,1)/=0.0_dp) &
        error stop 'empty condensed distances'
    print '(a)', 'test_spatial_distance PASS'

contains
    subroutine assert_close(label,actual,reference,tolerance)
        character(len=*), intent(in) :: label !! failure identification
        real(dp), intent(in) :: actual !! computed distance
        real(dp), intent(in) :: reference !! expected distance
        real(dp), intent(in) :: tolerance !! absolute comparison threshold
        if (.not. abs(actual-reference)<=tolerance) then
            print *,label, actual,reference
            error stop 'metric disagreement'
        end if
    end subroutine assert_close

    subroutine expect_invalid(label,actual,status)
        character(len=*), intent(in) :: label !! failing input case
        real(dp), intent(in) :: actual !! required NaN on error
        integer, intent(in) :: status !! required invalid-input status
        if (status/=spatial_distance_invalid_input .or. .not.ieee_is_nan(actual)) then
            print *,label,actual,status
            error stop 'invalid input not rejected'
        end if
    end subroutine expect_invalid
end program test_spatial_distance
