! SPDX-License-Identifier: MIT
program test_cluster_hierarchy_layout
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : linkage, inconsistent, is_isomorphic, &
        is_valid_im, is_valid_linkage, to_mlab_linkage, from_mlab_linkage, &
        dendrogram, cluster_dendrogram, leaves_list, cluster_success, &
        cluster_invalid_input
    implicit none

    real(dp) :: x(5,1), bad_im(1,4), nanval
    real(dp), allocatable :: z(:,:), mlab(:,:), again(:,:), im(:,:), chain(:,:)
    integer, allocatable :: reference(:), large1(:), large2(:)
    type(cluster_dendrogram) :: layout, sorted
    integer :: status, i, n

    if (.not. is_isomorphic([5,5,-7,100,-7],[2,2,9,8,9],status)) &
        error stop 'isomorphic relabeled assignment'
    if (status /= cluster_success) error stop 'isomorphic success status'
    if (is_isomorphic([5,5,-7,100,-7],[2,2,9,8,8])) &
        error stop 'isomorphic split and merge'
    if (is_isomorphic([1,1,2],[3,4,4])) error stop 'isomorphic mismatch'
    if (is_isomorphic([1,2],[3],status)) error stop 'isomorphic length'
    if (status /= cluster_invalid_input) error stop 'isomorphic bad size status'
    if (.not. is_isomorphic([integer::],[integer::],status)) &
        error stop 'isomorphic empty partition'
    if (status /= cluster_success) error stop 'isomorphic empty status'
    if (.not. is_isomorphic([huge(i),-huge(i),huge(i)],[-1,0,-1])) &
        error stop 'isomorphic extreme labels'
    allocate(large1(100001),large2(100001))
    large1 = [(mod(i,313),i=1,size(large1))]
    large2 = large1+100000
    if (.not. is_isomorphic(large1,large2)) error stop 'large partition isomorphism'
    large2(size(large2)) = large2(size(large2))-1
    if (is_isomorphic(large1,large2)) error stop 'large partition mismatch'

    x(:,1) = [0.0_dp,1.0_dp,2.0_dp,7.0_dp,10.0_dp]
    z = linkage(x,method='complete',status=status)
    if (status /= cluster_success .or. .not. is_valid_linkage(z)) &
        error stop 'linkage valid fixture'
    im = inconsistent(z,status=status)
    if (status /= cluster_success .or. .not. is_valid_im(im)) &
        error stop 'inconsistency validation'

    mlab = to_mlab_linkage(z,status)
    if (status /= cluster_success .or. size(mlab,1) /= 4 .or. &
        size(mlab,2) /= 3) error stop 'MATLAB conversion shape'
    if (any(mlab /= z(:,1:3))) error stop 'MATLAB conversion values'
    again = from_mlab_linkage(mlab,status)
    if (status /= cluster_success .or. any(again /= z)) &
        error stop 'MATLAB conversion roundtrip'
    mlab(1,1) = mlab(1,2)
    again = from_mlab_linkage(mlab,status)
    if (status /= cluster_invalid_input .or. size(again) /= 0) &
        error stop 'MATLAB duplicate child validation'
    mlab = to_mlab_linkage(z)
    mlab(2,1) = mlab(1,1)
    again = from_mlab_linkage(mlab,status)
    if (status /= cluster_invalid_input .or. size(again) /= 0) &
        error stop 'MATLAB repeated leaf validation'
    again = from_mlab_linkage(reshape([1.0_dp,2.0_dp,1.0_dp, &
        1.0_dp,2.0_dp,2.0_dp],[2,3]),status)
    if (status /= cluster_invalid_input .or. size(again) /= 0) &
        error stop 'MATLAB bad topology'

    layout = dendrogram(z,status)
    if (status /= cluster_success) error stop 'dendrogram default status'
    reference = leaves_list(z)
    if (any(layout%leaves /= reference)) error stop 'dendrogram leaf order'
    if (any(shape(layout%icoord) /= [4,4])) error stop 'dendrogram x shape'
    if (any(shape(layout%dcoord) /= [4,4])) error stop 'dendrogram y shape'
    if (abs(layout%icoord(4,2)-layout%icoord(4,1)) > 0.0_dp) &
        error stop 'dendrogram left U vertex'
    if (abs(layout%dcoord(4,2)-z(4,3)) > 1.0e-14_dp) &
        error stop 'dendrogram root height'
    sorted = dendrogram(z,status,count_sort='ascending')
    if (status /= cluster_success .or. size(sorted%leaves) /= 5) &
        error stop 'dendrogram count ascending'
    sorted = dendrogram(z,status,count_sort='descending')
    if (status /= cluster_success .or. size(sorted%leaves) /= 5) &
        error stop 'dendrogram count descending'
    sorted = dendrogram(z,status,distance_sort='ascending')
    if (status /= cluster_success .or. size(sorted%leaves) /= 5) &
        error stop 'dendrogram distance ascending'
    sorted = dendrogram(z,status,distance_sort='descending')
    if (status /= cluster_success .or. size(sorted%leaves) /= 5) &
        error stop 'dendrogram distance descending'
    sorted = dendrogram(z,status,count_sort='broken')
    if (status /= cluster_invalid_input .or. size(sorted%leaves) /= 0) &
        error stop 'dendrogram invalid sort'
    sorted = dendrogram(z,status,count_sort='ascending',distance_sort='ascending')
    if (status /= cluster_invalid_input .or. size(sorted%icoord) /= 0) &
        error stop 'dendrogram mixed sort'

    again = z
    again(1,4) = 3.0_dp
    if (is_valid_linkage(again)) error stop 'linkage wrong leaf count accepted'
    sorted = dendrogram(again,status)
    if (status /= cluster_invalid_input .or. size(sorted%leaves) /= 0) &
        error stop 'dendrogram rejected invalid linkage'
    mlab = to_mlab_linkage(again,status)
    if (status /= cluster_invalid_input .or. size(mlab) /= 0) &
        error stop 'MATLAB rejected invalid linkage'

    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    bad_im(1,:) = [1.0_dp,0.1_dp,1.0_dp,0.0_dp]
    if (.not. is_valid_im(bad_im)) error stop 'valid inconsistency row'
    bad_im(1,1) = -1.0_dp
    if (is_valid_im(bad_im)) error stop 'negative mean accepted'
    bad_im(1,1) = 1.0_dp
    bad_im(1,2) = -0.1_dp
    if (is_valid_im(bad_im)) error stop 'negative std accepted'
    bad_im(1,2) = 0.1_dp
    bad_im(1,3) = -1.0_dp
    if (is_valid_im(bad_im)) error stop 'negative count accepted'
    bad_im(1,3) = nanval
    if (is_valid_im(bad_im)) error stop 'nonfinite count accepted'
    if (is_valid_im(reshape([real(dp)::],[0,4]))) &
        error stop 'empty inconsistency accepted'

    ! A large chain exercises the iterative, recursion-free geometry walk.
    n = 300
    allocate(chain(n-1,4))
    chain(1,:) = [1.0_dp,2.0_dp,1.0_dp,2.0_dp]
    do i = 2,n-1
        chain(i,:) = [real(i+1,dp),real(n+i-1,dp), &
            real(i,dp),real(i+1,dp)]
    end do
    if (.not. is_valid_linkage(chain)) error stop 'deep chain invalid'
    sorted = dendrogram(chain,status)
    if (status /= cluster_success .or. size(sorted%leaves) /= n .or. &
        size(sorted%icoord,1) /= n-1) error stop 'deep chain layout'
    print *, 'cluster hierarchy validation/MATLAB/layout regression PASS'
end program test_cluster_hierarchy_layout
