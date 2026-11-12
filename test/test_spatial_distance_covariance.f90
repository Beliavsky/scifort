! SPDX-License-Identifier: MIT
program test_spatial_distance_covariance
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_spatial, only : point_distance, pdist, cdist, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    real(dp) :: u(2), v(2), vi(2,2), x(5,2), y(3,2), result
    real(dp), allocatable :: distances(:), matrix(:,:), invalid(:), reference(:)
    integer :: code

    u = [1.0_dp,0.0_dp]
    v = [0.0_dp,1.0_dp]
    call assert_close('seuclidean point',point_distance(u,v,'seuclidean', &
        variance=[1.0_dp,4.0_dp],status=code),sqrt(1.25_dp))
    if (code /= spatial_distance_success) error stop 'seuclidean point status'

    vi(1,:) = [2.0_dp,0.25_dp]
    vi(2,:) = [0.25_dp,1.5_dp]
    result = point_distance(u,v,'mahalanobis',inv_cov=vi,status=code)
    call assert_close('mahalanobis point',result,sqrt(3.0_dp))
    if (code /= spatial_distance_success) error stop 'mahalanobis point status'

    result = point_distance(u,v,'jensenshannon',status=code)
    call assert_close('Jensen-Shannon nat',result,sqrt(log(2.0_dp)))
    if (code /= spatial_distance_success) error stop 'JS point status'
    result = point_distance(u,v,'jensenshannon',base=2.0_dp)
    call assert_close('Jensen-Shannon base 2',result,1.0_dp)
    result = point_distance(1.0e200_dp*u,1.0e200_dp*v,'jensenshannon',base=2.0_dp)
    call assert_close('Jensen-Shannon enormous probabilities',result,1.0_dp)
    result = point_distance(1.0e-200_dp*u,1.0e-200_dp*v,'jensenshannon',base=2.0_dp)
    call assert_close('Jensen-Shannon tiny probabilities',result,1.0_dp)
    result = point_distance([1.0_dp,2.0_dp],[2.0_dp,4.0_dp],'jensenshannon')
    call assert_close('Jensen-Shannon proportional vectors',result,0.0_dp)

    x(1,:) = [0.0_dp,0.0_dp]
    x(2,:) = [3.0_dp,4.0_dp]
    x(3,:) = [1.0_dp,2.0_dp]
    x(4,:) = [7.0_dp,-1.0_dp]
    x(5,:) = [2.0_dp,8.0_dp]
    y(1,:) = [1.0_dp,3.0_dp]
    y(2,:) = [-2.0_dp,2.0_dp]
    y(3,:) = [5.0_dp,1.0_dp]

    distances = pdist(x,'seuclidean',status=code)
    if (code/=0 .or. size(distances)/=10) error stop 'auto standardized Euclidean pdist'
    ! Sample variances of x(:,1), x(:,2), respectively.
    call assert_close('auto pdist seuclidean pair',distances(1), &
        sqrt(9.0_dp/7.3_dp+16.0_dp/12.8_dp))
    matrix = cdist(x,y,'seuclidean',status=code)
    if (code/=0 .or. any(shape(matrix)/=[5,3])) error stop 'auto cdist standardized Euclidean'
    matrix = cdist(x,y,'seuclidean',variance=[2.0_dp,4.0_dp],status=code)
    if (code/=0) error stop 'specified variance cdist status'
    call assert_close('specified variance cdist pair',matrix(1,1),sqrt(0.5_dp+2.25_dp))

    distances = pdist(x,'mahalanobis',status=code)
    if (code/=0 .or. size(distances)/=10) error stop 'auto mahalanobis pdist'
    if (.not. all(distances>=0.0_dp)) error stop 'invalid auto Mahalanobis distances'
    reference = distances
    distances = pdist(reshape([x(:,1)*1.0e8_dp,x(:,2)*1.0e-8_dp],[5,2]), &
        'mahalanobis',status=code)
    if (code/=0 .or. any(abs(distances-reference)>2.0e-12_dp)) &
        error stop 'scaled Mahalanobis geometry lost or falsely singular'
    matrix = cdist(x,y,'mahalanobis',status=code)
    if (code/=0 .or. any(shape(matrix)/=[5,3])) error stop 'auto mahalanobis cdist'
    matrix = cdist(x,y,'mahalanobis',inv_cov=vi,status=code)
    if (code/=0) error stop 'specified inv_cov cdist status'
    call assert_close('specified inv_cov cdist pair',matrix(1,1), &
        point_distance(x(1,:),y(1,:),'mahalanobis',inv_cov=vi))

    result = point_distance(u,v,'seuclidean',status=code)
    if (code/=spatial_distance_invalid_input .or. .not. ieee_is_nan(result)) &
        error stop 'pointwise seuclidean requires variance'
    result = point_distance(u,v,'mahalanobis',status=code)
    if (code/=spatial_distance_invalid_input .or. .not. ieee_is_nan(result)) &
        error stop 'pointwise mahalanobis requires inverse covariance'
    result = point_distance(u,v,'seuclidean',variance=[1.0_dp],status=code)
    if (code/=spatial_distance_invalid_input) error stop 'variance shape accepted'
    result = point_distance(u,v,'seuclidean',variance=[0.0_dp,1.0_dp],status=code)
    if (code/=spatial_distance_invalid_input) error stop 'zero variance accepted'
    result = point_distance(u,v,'mahalanobis',inv_cov=reshape([1.0_dp],[1,1]),status=code)
    if (code/=spatial_distance_invalid_input) error stop 'inverse covariance shape accepted'
    result = point_distance(u,v,'jensenshannon',base=1.0_dp,status=code)
    if (code/=spatial_distance_invalid_input) error stop 'log base one accepted'
    result = point_distance([-1.0_dp,2.0_dp],v,'jensenshannon',status=code)
    if (code/=spatial_distance_invalid_input) error stop 'negative mass accepted'
    result = point_distance([0.0_dp,0.0_dp],v,'jensenshannon',status=code)
    if (code/=spatial_distance_invalid_input) error stop 'zero mass accepted'
    result = point_distance(u,v,'jensenshannon',w=[1.0_dp,1.0_dp],status=code)
    if (code/=spatial_distance_invalid_input) error stop 'JS weights accepted'
    invalid = pdist(x(1:2,:),'mahalanobis',status=code)
    if (code/=spatial_distance_invalid_input .or. size(invalid)/=0) &
        error stop 'underidentified covariance accepted'
    invalid = pdist(transpose(reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp,3.0_dp,6.0_dp], &
        [2,3])),'mahalanobis',status=code)
    if (code/=spatial_distance_invalid_input .or. size(invalid)/=0) &
        error stop 'rank-deficient covariance accepted'
    matrix = cdist(x,y,'seuclidean',variance=[1.0_dp],status=code)
    if (code/=spatial_distance_invalid_input .or. any(shape(matrix)/=[0,0])) &
        error stop 'cdist variance shape accepted'

    print '(a)', 'test_spatial_distance_covariance PASS'
contains
    subroutine assert_close(label, a, b)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        if (.not. (abs(a-b) <= 2.0e-12_dp*max(1.0_dp,abs(b)))) then
            print *, label,' got ',a,' expected ',b
            error stop 'distance mismatch'
        end if
    end subroutine assert_close
end program test_spatial_distance_covariance
