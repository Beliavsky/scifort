! SPDX-License-Identifier: MIT
program test_cluster_hierarchy_flatdata
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : pdist
    use scifort_cluster_hierarchy, only : linkage, fcluster, fclusterdata, &
        fcluster_distance, fcluster_maxclust, fcluster_inconsistent, &
        fcluster_monocrit, fcluster_maxclust_monocrit, inconsistent, &
        cophenetic_correlation, cluster_success, cluster_invalid_input
    implicit none

    real(dp) :: x(4,1), four_groups(12,2), bad_matrix(2,2)
    real(dp), allocatable :: z(:,:), y(:), r(:,:), crit(:)
    integer, allocatable :: labels(:), reference(:)
    real(dp) :: corr, nanvalue
    integer :: status, k

    x(:,1) = [0.0_dp, 1.0_dp, 4.0_dp, 10.0_dp]
    z = linkage(x,status=status)
    if (status /= cluster_success) error stop 'flatdata linkage fixture'

    labels = fcluster(z,3.0_dp,criterion='distance',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fcluster distance'
    labels = fcluster(z,2.0_dp,criterion='maxclust',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fcluster maxclust'
    labels = fcluster(z,2,criterion='maxclust',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fcluster integer overload'
    labels = fcluster(z,2.9_dp,criterion='maxclust',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fcluster fractional count'
    labels = fcluster(z,0.5_dp,status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,2,3])) &
        error stop 'fcluster default inconsistent'
    labels = fclusterdata(x,0.5_dp,status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,2,3])) &
        error stop 'fclusterdata default inconsistent'

    r = inconsistent(z,depth=2,status=status)
    if (status /= cluster_success) error stop 'inconsistent fixture'
    labels = fcluster(z,0.5_dp,criterion='inconsistent',depth=2,r=r,status=status)
    if (status /= cluster_success .or. &
        any(labels /= fcluster_inconsistent(z,0.5_dp,r=r))) &
        error stop 'fcluster supplied R'
    labels = fclusterdata(x,0.5_dp,r=r,status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,2,3])) &
        error stop 'fclusterdata supplied R'

    crit = [1.0_dp,3.0_dp,6.0_dp]
    labels = fcluster(z,3.0_dp,criterion='monocrit',monocrit=crit,status=status)
    if (status /= cluster_success .or. &
        any(labels /= fcluster_monocrit(z,crit,3.0_dp))) &
        error stop 'fcluster monocrit'
    labels = fcluster(z,2.0_dp,criterion='maxclust_monocrit',monocrit=crit,status=status)
    if (status /= cluster_success .or. &
        any(labels /= fcluster_maxclust_monocrit(z,crit,2))) &
        error stop 'fcluster maxclust_monocrit'

    labels = fclusterdata(x,3.0_dp,criterion='distance',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fclusterdata distance'
    labels = fclusterdata(x,2.0_dp,criterion='maxclust',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fclusterdata maxclust'
    labels = fclusterdata(x,2,criterion='maxclust',status=status)
    if (status /= cluster_success .or. any(labels /= [1,1,1,2])) &
        error stop 'fclusterdata integer overload'

    ! The cophenetic correlation is independent of all labels and is
    ! compared to the reference SciPy result to near-double precision.
    y = pdist(x)
    corr = cophenetic_correlation(z,y,status)
    if (status /= cluster_success .or. &
        abs(corr-0.9206549208873558_dp) > 1.0e-13_dp) &
        error stop 'cophenetic correlation reference'
    corr = cophenetic_correlation(z,y*1.0e200_dp,status)
    if (status /= cluster_success .or. &
        abs(corr-0.9206549208873558_dp) > 1.0e-13_dp) &
        error stop 'cophenetic correlation scaling'

    ! Four naturally separated groups in two dimensions, matching the
    ! 12-point example in the public SciPy documentation.
    four_groups(1,:) = [0.0_dp,0.0_dp]
    four_groups(2,:) = [0.0_dp,1.0_dp]
    four_groups(3,:) = [1.0_dp,0.0_dp]
    four_groups(4,:) = [0.0_dp,4.0_dp]
    four_groups(5,:) = [0.0_dp,3.0_dp]
    four_groups(6,:) = [1.0_dp,4.0_dp]
    four_groups(7,:) = [4.0_dp,0.0_dp]
    four_groups(8,:) = [3.0_dp,0.0_dp]
    four_groups(9,:) = [4.0_dp,1.0_dp]
    four_groups(10,:) = [4.0_dp,4.0_dp]
    four_groups(11,:) = [3.0_dp,4.0_dp]
    four_groups(12,:) = [4.0_dp,3.0_dp]
    reference = [( (k-1)/3+1, k=1,12)]
    labels = fclusterdata(four_groups,1.0_dp,status=status)
    if (status /= cluster_success .or. .not. same_partition(labels,reference)) &
        error stop 'fclusterdata scipy documentation example'
    labels = fclusterdata(four_groups,1.5_dp,criterion='distance', &
        metric='cityblock',method='complete',status=status)
    z = linkage(four_groups,method='complete',metric='cityblock')
    if (status /= cluster_success .or. &
        any(labels /= fcluster_distance(z,1.5_dp))) &
        error stop 'fclusterdata cityblock complete'
    z = linkage(x)

    ! Input and option failures must report an error and no labels.
    labels = fcluster(z,2.0_dp,criterion='unrecognized',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fcluster unrecognized criterion'
    labels = fcluster(z,0.0_dp,criterion='maxclust',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fcluster zero maxclust'
    labels = fcluster(z,2.0_dp,criterion='monocrit',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fcluster missing monocrit'
    labels = fcluster(z,2.0_dp,criterion='maxclust_monocrit',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fcluster missing maxclust_monocrit'
    labels = fcluster(z,2.0_dp,criterion='monocrit', &
        monocrit=[3.0_dp,2.0_dp,1.0_dp],status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fcluster nonmonotone monocrit'
    labels = fclusterdata(x,1.0_dp,criterion='monocrit',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata unsupported monocrit without vector'
    labels = fclusterdata(x,1.0_dp,metric='not_a_metric',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata invalid metric'
    labels = fclusterdata(x,1.0_dp,method='not_a_method',status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata invalid method'
    labels = fclusterdata(x,1.0_dp,r=r(:2,:),status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata invalid R'
    labels = fclusterdata(x(:1,:),1.0_dp,status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata insufficient observations'
    labels = fclusterdata(x,1.0_dp,depth=-1,status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata negative depth'
    nanvalue = ieee_value(0.0_dp,ieee_quiet_nan)
    bad_matrix = reshape([0.0_dp,1.0_dp,2.0_dp,nanvalue],[2,2])
    labels = fclusterdata(bad_matrix,1.0_dp,status=status)
    if (status /= cluster_invalid_input .or. size(labels) /= 0) &
        error stop 'fclusterdata NaN observations'
    corr = cophenetic_correlation(z,y(:3),status)
    if (status /= cluster_invalid_input .or. .not. ieee_is_nan(corr)) &
        error stop 'cophenetic correlation wrong length'
    y(1) = nanvalue
    corr = cophenetic_correlation(z,y,status)
    if (status /= cluster_invalid_input .or. .not. ieee_is_nan(corr)) &
        error stop 'cophenetic correlation NaN input'

    y = pdist(x(:2,:))
    z = linkage(y,status=status)
    corr = cophenetic_correlation(z,y,status)
    if (status /= cluster_success .or. .not. ieee_is_nan(corr)) &
        error stop 'two leaves have undefined cophenetic correlation'
    z = linkage([0.0_dp,0.0_dp,0.0_dp],status=status)
    corr = cophenetic_correlation(z,[0.0_dp,0.0_dp,0.0_dp],status)
    if (status /= cluster_success .or. .not. ieee_is_nan(corr)) &
        error stop 'constant distances undefined correlation'
    print *, 'cluster hierarchy flatdata/cophenetic regression PASS'

contains

    pure logical function same_partition(a,b) result(ok)
        integer, intent(in) :: a(:) !! First cluster labeling.
        integer, intent(in) :: b(:) !! Second cluster labeling.
        integer :: i, j
        ok = .false.
        if (size(a) /= size(b)) return
        do i = 1,size(a)
            do j = i+1,size(a)
                if ((a(i) == a(j)) .neqv. (b(i) == b(j))) return
            end do
        end do
        ok = .true.
    end function same_partition

end program test_cluster_hierarchy_flatdata
