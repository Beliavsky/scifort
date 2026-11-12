! SPDX-License-Identifier: MIT
program test_cluster_hierarchy
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : linkage, fcluster_distance, &
        fcluster_maxclust, cophenetic_distances, leaves_list, &
        inconsistent, maxdists, maxrstat, maxinconsts, fcluster_inconsistent, &
        fcluster_monocrit, fcluster_maxclust_monocrit, cut_tree, &
        is_monotonic, num_obs_linkage, &
        cluster_success, cluster_invalid_input
    implicit none

    real(dp) :: x(4,1), z_expected(3,4), triangle(3,2)
    real(dp), allocatable :: z(:,:), y(:), r(:,:), r0(:,:), stats(:)
    integer, allocatable :: labels(:), leaves(:), cuts(:,:)
    integer :: status, j
    character(len=8) :: methods(7)

    x(:,1)=[0.0_dp,1.0_dp,4.0_dp,10.0_dp]
    z_expected(1,:)=[1.0_dp,2.0_dp,1.0_dp,2.0_dp]
    z_expected(2,:)=[3.0_dp,5.0_dp,3.0_dp,3.0_dp]
    z_expected(3,:)=[4.0_dp,6.0_dp,6.0_dp,4.0_dp]
    z=linkage(x,method='single',status=status)
    if (status/=cluster_success) error stop 'linkage observations status'
    if (any(shape(z)/=[3,4])) error stop 'linkage observations shape'
    if (any(abs(z-z_expected)>1.0e-13_dp)) error stop 'linkage single example'
    y=cophenetic_distances(z,status)
    if (status/=cluster_success) error stop 'cophenetic status'
    if (any(abs(y-[1.0_dp,3.0_dp,6.0_dp,3.0_dp,6.0_dp,6.0_dp])>1.0e-13_dp)) &
        error stop 'cophenetic values'
    labels=fcluster_distance(z,3.0_dp,status)
    if (status/=cluster_success) error stop 'fcluster distance status'
    if (any(labels/=[1,1,1,2])) error stop 'fcluster distance values'
    labels=fcluster_maxclust(z,2,status)
    if (status/=cluster_success) error stop 'fcluster maxclust status'
    if (any(labels/=[1,1,1,2])) error stop 'fcluster maxclust values'
    labels=fcluster_maxclust(z,4,status)
    if (status/=cluster_success) error stop 'fcluster maxclust singleton status'
    if (any(labels/=[1,2,3,4])) error stop 'maxclust singleton'
    leaves=leaves_list(z,status)
    if (status/=cluster_success) error stop 'leaves status'
    if (any(leaves/=[4,3,1,2])) error stop 'leaves values'

    methods=[character(len=8) :: 'single', 'complete', 'average', 'weighted', &
        'centroid', 'median', 'ward']
    do j=1,size(methods)
        z=linkage(x,method=trim(methods(j)),status=status)
        if (status/=cluster_success) error stop 'linkage methods status'
        if (any(shape(z)/=[3,4])) error stop 'linkage method shape'
        if (abs(z(3,4)-4.0_dp)>1.0e-14_dp) error stop 'wrong cluster sizes'
        y=cophenetic_distances(z,status)
        if (status/=cluster_success .or. size(y)/=6) error stop 'cophenetic methods'
        labels=fcluster_maxclust(z,1,status)
        if (status/=cluster_success) error stop 'fcluster method status'
        if (any(labels/=1)) error stop 'fcluster method not merged'
        leaves=leaves_list(z,status)
        if (status/=cluster_success .or. size(leaves)/=4) error stop 'leaves method'
    end do
    z=linkage([1.0_dp,4.0_dp,10.0_dp,3.0_dp,9.0_dp,6.0_dp], &
        method='single',status=status)
    if (status/=cluster_success) error stop 'condensed distances'
    if (any(abs(z-z_expected)>1.0e-13_dp)) error stop 'condensed matrix mismatch'


    ! Distinct observations may still produce nonmonotone centroid linkage.
    ! Cut criteria use maximum height over every descendant of a node.
    triangle(1,:)=[0.0_dp,0.0_dp]
    triangle(2,:)=[2.0_dp,0.0_dp]
    triangle(3,:)=[1.0_dp,1.8_dp]
    z=linkage(triangle,'centroid',status=status)
    if (status/=cluster_success) error stop 'centroid inversion construction'
    if (z(2,3)>=z(1,3)) error stop 'centroid inversion fixture'
    labels=fcluster_distance(z,1.9_dp,status)
    if (status/=cluster_success .or. size(labels)/=3) error stop 'inverted tree status'
    if (size(labels)/=size(pack(labels,labels>=1))) error stop 'inverted tree labels'
    if (labels(1)==labels(2) .or. labels(1)==labels(3) .or. &
        labels(2)==labels(3)) error stop 'inverted tree illegal merge'

    ! Tied distances may create fewer than the requested maximum clusters.
    z=linkage([0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], &
        method='average',status=status)
    if (status/=cluster_success) error stop 'duplicate points linkage'
    if (any(z(:,3)/=0.0_dp)) error stop 'duplicate point distances'
    labels=fcluster_maxclust(z,3,status)
    if (status/=cluster_success .or. any(labels/=1)) error stop 'tied maximum clusters'
    z=linkage(reshape([0.0_dp,1.0e200_dp,2.0e200_dp], [3,1]), &
        method='ward',status=status)
    if (status/=cluster_success) error stop 'large-scale ward'
    if (any(z(:,3)<0.0_dp)) error stop 'large-scale ward negative'

    z=linkage([1.0_dp,2.0_dp],status=status)
    if (status/=cluster_invalid_input .or. size(z)/=0) error stop 'invalid triangular input'
    z=linkage([1.0_dp,-2.0_dp,3.0_dp],status=status)
    if (status/=cluster_invalid_input .or. size(z)/=0) error stop 'negative distances'
    z=linkage([1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),3.0_dp],status=status)
    if (status/=cluster_invalid_input .or. size(z)/=0) error stop 'nonfinite distances'
    z=linkage(x,method='invalid_method',status=status)
    if (status/=cluster_invalid_input .or. size(z)/=0) error stop 'invalid method'
    z=linkage(x,metric='invalid_metric',status=status)
    if (status/=cluster_invalid_input .or. size(z)/=0) error stop 'invalid metric'

    z=linkage(x,method='complete')
    z(1,4)=3.0_dp
    labels=fcluster_distance(z,1.0_dp,status)
    if (status/=cluster_invalid_input .or. size(labels)/=0) error stop 'invalid tree counts'
    y=cophenetic_distances(z,status)
    if (status/=cluster_invalid_input .or. size(y)/=0) error stop 'invalid tree coph'
    leaves=leaves_list(z,status)
    if (status/=cluster_invalid_input .or. size(leaves)/=0) error stop 'invalid tree leaves'

    ! SciPy hierarchy diagnostics and additional flat-cut criteria.
    z = linkage(x,method='single',status=status)
    if (status/=cluster_success) error stop 'diagnostic linkage'
    if (num_obs_linkage(z,status)/=4 .or. status/=cluster_success) &
        error stop 'num_obs_linkage'
    if (.not. is_monotonic(z,status) .or. status/=cluster_success) &
        error stop 'monotonic single'
    r = inconsistent(z,status=status)
    if (status/=cluster_success .or. any(shape(r)/=[3,4])) &
        error stop 'inconsistent dimensions'
    if (maxval(abs(r(:,1)-[1.0_dp,2.0_dp,4.5_dp]))>1.0e-13_dp) &
        error stop 'inconsistent means depth 2'
    if (maxval(abs(r(:,2)-[0.0_dp,sqrt(2.0_dp),sqrt(4.5_dp)]))>1.0e-13_dp) &
        error stop 'inconsistent sample sd'
    if (any(r(:,3)/=[1.0_dp,2.0_dp,2.0_dp])) error stop 'inconsistent counts'
    if (maxval(abs(r(:,4)-[0.0_dp,sqrt(0.5_dp),sqrt(0.5_dp)]))>1.0e-13_dp) &
        error stop 'inconsistent coefficients'
    r = inconsistent(z,depth=3,status=status)
    if (status/=cluster_success) error stop 'inconsistent depth 3'
    if (abs(r(3,1)-10.0_dp/3.0_dp)>1.0e-13_dp) &
        error stop 'inconsistent depth 3 mean'
    if (abs(r(3,2)-sqrt(19.0_dp/3.0_dp))>1.0e-13_dp) &
        error stop 'inconsistent depth 3 sd'
    if (r(3,3)/=3.0_dp) error stop 'inconsistent depth 3 count'
    r0 = inconsistent(z,depth=0,status=status)
    if (status/=cluster_success .or. any(r0(:,3)/=1.0_dp)) &
        error stop 'inconsistent zero depth'
    r0 = inconsistent(z,depth=100,status=status)
    if (status/=cluster_success .or. maxval(abs(r0-r))>1.0e-13_dp) &
        error stop 'inconsistent full depth'
    stats = maxdists(z,status)
    if (status/=cluster_success .or. any(stats/=[1.0_dp,3.0_dp,6.0_dp])) &
        error stop 'maxdists'
    stats = maxrstat(z,r,1,status)
    if (status/=cluster_success .or. maxval(abs(stats-[1.0_dp,2.0_dp,10.0_dp/3.0_dp]))>1.0e-13_dp) &
        error stop 'maxRstat mean'
    stats = maxinconsts(z,r,status)
    if (status/=cluster_success .or. abs(stats(3)-r(3,4))>1.0e-13_dp) &
        error stop 'maxinconsts'
    labels = fcluster_inconsistent(z,0.8_dp,status,depth=3)
    if (status/=cluster_success .or. any(labels/=[1,1,1,2])) &
        error stop 'inconsistent flat depth 3'
    labels = fcluster_inconsistent(z,0.8_dp,status,r=r)
    if (status/=cluster_success .or. any(labels/=[1,1,1,2])) &
        error stop 'inconsistent provided R'
    labels = fcluster_inconsistent(z,0.8_dp,status)
    if (status/=cluster_success .or. any(labels/=1)) &
        error stop 'inconsistent flat default'
    stats = maxdists(z,status)
    labels = fcluster_monocrit(z,stats,3.0_dp,status)
    if (status/=cluster_success .or. any(labels/=[1,1,1,2])) &
        error stop 'monocrit flat cutoff'
    labels = fcluster_maxclust_monocrit(z,stats,2,status)
    if (status/=cluster_success .or. any(labels/=[1,1,1,2])) &
        error stop 'monocrit maxclusters two'
    labels = fcluster_maxclust_monocrit(z,stats,99,status)
    if (status/=cluster_success .or. any(labels/=[1,2,3,4])) &
        error stop 'monocrit maxclusters singleton'
    labels = fcluster_maxclust_monocrit(z,stats,1,status)
    if (status/=cluster_success .or. any(labels/=1)) &
        error stop 'monocrit maxclusters one'
    cuts = cut_tree(z,status=status)
    if (status/=cluster_success .or. any(shape(cuts)/=[4,4])) &
        error stop 'cut_tree full dimension'
    if (any(cuts(:,1)/=[1,2,3,4]) .or. any(cuts(:,2)/=[1,1,2,3]) .or. &
        any(cuts(:,3)/=[1,1,1,2]) .or. any(cuts(:,4)/=1)) &
        error stop 'cut_tree all levels'
    cuts = cut_tree(z,n_clusters=[1,4,2,3,2],status=status)
    if (status/=cluster_success .or. any(shape(cuts)/=[4,5])) &
        error stop 'cut_tree selected counts dimensions'
    if (any(cuts(:,1)/=1) .or. any(cuts(:,2)/=[1,2,3,4]) .or. &
        any(cuts(:,3)/=[1,1,1,2]) .or. any(cuts(:,4)/=[1,1,2,3]) .or. &
        any(cuts(:,5)/=[1,1,1,2])) &
        error stop 'cut_tree selected counts'
    cuts = cut_tree(z,height=[0.0_dp,1.0_dp,1.5_dp,3.0_dp,3.01_dp,7.0_dp],status=status)
    if (status/=cluster_success .or. any(shape(cuts)/=[4,6])) &
        error stop 'cut_tree heights dimensions'
    if (any(cuts(:,1)/=[1,2,3,4]) .or. any(cuts(:,2)/=[1,2,3,4]) .or. &
        any(cuts(:,3)/=[1,1,2,3]) .or. any(cuts(:,4)/=[1,1,2,3]) .or. &
        any(cuts(:,5)/=[1,1,1,2]) .or. any(cuts(:,6)/=1)) &
        error stop 'cut_tree heights strict cutoff'
    cuts = cut_tree(z,height=[7.0_dp,-1.0_dp,3.1_dp],status=status)
    if (status/=cluster_success .or. any(cuts(:,1)/=1) .or. &
        any(cuts(:,2)/=[1,2,3,4]) .or. any(cuts(:,3)/=[1,1,1,2])) &
        error stop 'cut_tree unsorted heights'
    cuts = cut_tree(z,n_clusters=[4,0],status=status)
    if (status/=cluster_invalid_input .or. size(cuts)/=0) &
        error stop 'cut_tree invalid count'
    cuts = cut_tree(z,n_clusters=[2],height=[3.0_dp],status=status)
    if (status/=cluster_invalid_input .or. size(cuts)/=0) &
        error stop 'cut_tree conflicting options'
    r0 = inconsistent(z,depth=-1,status=status)
    if (status/=cluster_invalid_input .or. size(r0)/=0) &
        error stop 'inconsistent invalid depth'
    stats = maxrstat(z,r,0,status)
    if (status/=cluster_invalid_input .or. size(stats)/=0) &
        error stop 'maxRstat invalid column'
    stats = maxrstat(z,r(:2,:),1,status)
    if (status/=cluster_invalid_input .or. size(stats)/=0) &
        error stop 'maxRstat invalid shape'
    labels = fcluster_monocrit(z,[6.0_dp,3.0_dp,1.0_dp],2.0_dp,status)
    if (status/=cluster_invalid_input .or. size(labels)/=0) &
        error stop 'monocrit invalid monotonicity'
    labels = fcluster_maxclust_monocrit(z,[1.0_dp,3.0_dp,6.0_dp],0,status)
    if (status/=cluster_invalid_input .or. size(labels)/=0) &
        error stop 'maxclust_monocrit invalid count'
    labels = fcluster_inconsistent(z,0.5_dp,status,r=r(:2,:))
    if (status/=cluster_invalid_input .or. size(labels)/=0) &
        error stop 'inconsistent flat invalid R'
    z = linkage(triangle,'centroid',status=status)
    if (is_monotonic(z,status) .or. status/=cluster_success) &
        error stop 'centroid inversion monotonicity'
    cuts = cut_tree(z,height=[1.0_dp],status=status)
    if (status/=cluster_invalid_input .or. size(cuts)/=0) &
        error stop 'cut_tree height inversion'
    cuts = cut_tree(z,n_clusters=[2],status=status)
    if (status/=cluster_success .or. any(shape(cuts)/=[3,1])) &
        error stop 'cut_tree inverted count'
    stats = maxdists(z,status)
    if (status/=cluster_success .or. stats(2)<stats(1)) &
        error stop 'maxdists nonmonotone tree'

    ! Equal merge heights must give an exact count cut but may collapse
    ! several levels at once for threshold-based fcluster.
    z = linkage([0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], &
        method='average',status=status)
    r = inconsistent(z,status=status)
    if (status/=cluster_success .or. any(r(:,2)/=0.0_dp) .or. &
        any(r(:,4)/=0.0_dp)) error stop 'inconsistent all equal'
    cuts = cut_tree(z,n_clusters=[2],status=status)
    if (status/=cluster_success .or. size(cuts,1)/=4 .or. &
        size(cuts,2)/=1 .or. size(pack(cuts(:,1),cuts(:,1)==cuts(1,1)))==4) &
        error stop 'cut_tree ties exact counts'
    cuts = cut_tree(z,height=[0.0_dp,1.0_dp],status=status)
    if (status/=cluster_success .or. any(cuts(:,1)/=[1,2,3,4]) .or. &
        any(cuts(:,2)/=1)) error stop 'cut_tree ties strict height'
    z = linkage(reshape([0.0_dp,1.0e200_dp,3.0e200_dp],[3,1]), &
        method='single',status=status)
    if (status/=cluster_success) error stop 'large diagnostic linkage'
    r = inconsistent(z,depth=3,status=status)
    if (status/=cluster_success .or. .not. all(ieee_is_finite(r))) &
        error stop 'large diagnostic stats finite'
    if (abs(r(2,4)-sqrt(0.5_dp))>1.0e-12_dp) &
        error stop 'large diagnostic coefficient'

    z = linkage([4.0_dp],status=status)
    if (status/=cluster_success) error stop 'two leaf linkage'
    r = inconsistent(z,status=status)
    if (status/=cluster_success .or. any(shape(r)/=[1,4])) &
        error stop 'two leaf diagnostics'
    if (any(r(1,:)/=[4.0_dp,0.0_dp,1.0_dp,0.0_dp])) &
        error stop 'two leaf statistics'
    cuts = cut_tree(z,status=status)
    if (status/=cluster_success .or. any(shape(cuts)/=[2,2])) &
        error stop 'two leaf cut tree shape'
    if (any(cuts(:,1)/=[1,2]) .or. any(cuts(:,2)/=1)) &
        error stop 'two leaf cut values'
    cuts = cut_tree(z,n_clusters=[integer::],status=status)
    if (status/=cluster_success .or. any(shape(cuts)/=[2,0])) &
        error stop 'empty cut selection'
    labels = fcluster_monocrit(z,[ieee_value(0.0_dp,ieee_quiet_nan)],3.0_dp,status)
    if (status/=cluster_invalid_input .or. size(labels)/=0) &
        error stop 'nonfinite monocrit'
    z(1,4) = 3.0_dp
    r = inconsistent(z,status=status)
    if (status/=cluster_invalid_input .or. size(r)/=0) &
        error stop 'invalid diagnostic linkage'
    cuts = cut_tree(z,status=status)
    if (status/=cluster_invalid_input .or. size(cuts)/=0) &
        error stop 'invalid cut tree linkage'

    print '(a)', 'cluster hierarchy regression PASS'
end program test_cluster_hierarchy
