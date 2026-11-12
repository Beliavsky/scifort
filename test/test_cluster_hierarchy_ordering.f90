! SPDX-License-Identifier: MIT
program test_cluster_hierarchy_ordering
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : linkage, leaves_list, &
        optimal_leaf_ordering, to_tree, cluster_tree_node, leaders, correspond, &
        cluster_success, cluster_invalid_input
    use scifort_spatial_distance, only : pdist
    implicit none
    real(dp) :: x(4,2), x2(3,2)
    real(dp), allocatable :: z(:,:), orig(:,:), ordered(:,:), again(:,:), &
        y(:), bad(:)
    type(cluster_tree_node), allocatable :: nodes(:)
    integer, allocatable :: leaf_order(:), pairs(:,:)
    integer :: status, j

    x(:,1) = [0.0_dp,1.0_dp,4.0_dp,10.0_dp]
    x(:,2) = [0.0_dp,0.5_dp,-1.0_dp,0.0_dp]
    y = pdist(x)
    z = linkage(y,method='single',status=status)
    if (status /= cluster_success) error stop 'ordering setup'

    ordered = optimal_leaf_ordering(z,y,status)
    if (status /= cluster_success) error stop 'ordering condensed status'
    if (any(shape(ordered) /= shape(z))) error stop 'ordering condensed shape'
    if (any(ordered(:,3:4) /= z(:,3:4))) error stop 'ordering heights/counts changed'
    do j = 1,size(z,1)
        if (any(ordered(j,1:2) /= z(j,1:2)) .and. &
            any(ordered(j,1:2) /= [z(j,2),z(j,1)])) &
            error stop 'ordering changed topology'
    end do
    leaf_order = leaves_list(ordered,status)
    if (status /= cluster_success .or. size(leaf_order) /= 4) &
        error stop 'ordered leaves list'
    if (path_cost(leaf_order,y) > path_cost(leaves_list(z),y)+1.0e-12_dp) &
        error stop 'optimal ordering increased path length'

    again = optimal_leaf_ordering(z,x,status)
    if (status /= cluster_success .or. any(again /= ordered)) &
        error stop 'ordering observations overload'
    orig = linkage(y,method='single',optimal_ordering=.true.,status=status)
    if (status /= cluster_success .or. any(orig /= ordered)) &
        error stop 'linkage condensed optimal_ordering option'
    orig = linkage(x,method='single',optimal_ordering=.true.,status=status)
    if (status /= cluster_success .or. any(orig /= ordered)) &
        error stop 'linkage observations optimal_ordering option'

    nodes = to_tree(z,status)
    if (status /= cluster_success .or. size(nodes) /= 7) &
        error stop 'to_tree shape'
    do j = 1,4
        if (nodes(j)%id /= j .or. nodes(j)%count /= 1 .or. &
            nodes(j)%left /= 0 .or. nodes(j)%right /= 0 .or. &
            nodes(j)%parent == 0) error stop 'to_tree leaf'
    end do
    if (nodes(7)%id /= 7 .or. nodes(7)%count /= 4 .or. &
        nodes(7)%parent /= 0 .or. nodes(7)%distance /= z(3,3)) &
        error stop 'to_tree root'
    do j = 5,7
        if (nodes(j)%count /= int(z(j-4,4))) error stop 'to_tree counts'
        if (nodes(j)%left /= int(z(j-4,1)) .or. &
            nodes(j)%right /= int(z(j-4,2))) error stop 'to_tree children'
        if (nodes(nodes(j)%left)%parent /= j .or. &
            nodes(nodes(j)%right)%parent /= j) error stop 'to_tree parents'
    end do

    pairs = leaders(z,[11,11,11,-7],status)
    if (status /= cluster_success .or. any(shape(pairs) /= [2,2])) &
        error stop 'leaders shape'
    if (.not. any(pairs(:,1) == 6 .and. pairs(:,2) == 11) .or. &
        .not. any(pairs(:,1) == 4 .and. pairs(:,2) == -7)) &
        error stop 'leaders subtrees'
    pairs = leaders(z,[7,7,7,7],status)
    if (status /= cluster_success .or. any(shape(pairs) /= [1,2])) &
        error stop 'leaders single cluster shape'
    if (any(pairs(1,:) /= [7,7])) error stop 'leaders root'
    pairs = leaders(z,[1,2,3,4],status)
    if (status /= cluster_success .or. any(shape(pairs) /= [4,2])) &
        error stop 'leaders singleton shape'
    pairs = leaders(z,[0,8,0,9],status)
    if (status /= cluster_invalid_input .or. size(pairs) /= 0) &
        error stop 'leaders disconnected cluster must reject'
    pairs = leaders(z,[1,2,3],status)
    if (status /= cluster_invalid_input .or. size(pairs) /= 0) &
        error stop 'leaders invalid label count'

    if (.not. correspond(z,y,status) .or. status /= cluster_success) &
        error stop 'correspond matching sizes'
    if (correspond(z,y(:4),status) .or. status /= cluster_success) &
        error stop 'correspond mismatching sizes'
    if (correspond(z,[-1.0_dp,1.0_dp],status) .or. &
        status /= cluster_invalid_input) error stop 'correspond bad distances'

    bad = [1.0_dp,-1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    ordered = optimal_leaf_ordering(z,bad,status)
    if (status /= cluster_invalid_input .or. size(ordered) /= 0) &
        error stop 'ordering negative distances'
    bad(2) = ieee_value(0.0_dp,ieee_quiet_nan)
    ordered = optimal_leaf_ordering(z,bad,status)
    if (status /= cluster_invalid_input .or. size(ordered) /= 0) &
        error stop 'ordering nonfinite distances'
    ordered = optimal_leaf_ordering(z,y(:3),status)
    if (status /= cluster_invalid_input .or. size(ordered) /= 0) &
        error stop 'ordering triangular mismatch'
    ordered = optimal_leaf_ordering(z,x(:3,:),status)
    if (status /= cluster_invalid_input .or. size(ordered) /= 0) &
        error stop 'ordering observation count mismatch'
    ordered = optimal_leaf_ordering(z,x,status,metric='bad-metric')
    if (status /= cluster_invalid_input .or. size(ordered) /= 0) &
        error stop 'ordering invalid metric'
    orig = z
    orig(1,4) = 3.0_dp
    ordered = optimal_leaf_ordering(orig,y,status)
    if (status /= cluster_invalid_input .or. size(ordered) /= 0) &
        error stop 'ordering invalid linkage'
    nodes = to_tree(orig,status)
    if (status /= cluster_invalid_input .or. size(nodes) /= 0) &
        error stop 'to_tree invalid linkage'

    z = linkage([0.0_dp,0.0_dp,0.0_dp],method='single',status=status)
    if (status /= cluster_success) error stop 'tie fixture setup'
    ordered = optimal_leaf_ordering(z,[0.0_dp,0.0_dp,0.0_dp],status)
    if (status /= cluster_success .or. any(ordered(:,3:4) /= z(:,3:4))) &
        error stop 'all-equal distances'
    z = linkage([1.0e250_dp,3.0e250_dp,2.0e250_dp],method='average')
    ordered = optimal_leaf_ordering(z,[1.0e250_dp,3.0e250_dp,2.0e250_dp],status)
    if (status /= cluster_success .or. size(ordered,1) /= 2) &
        error stop 'large finite distances'

    x2(1,:) = [0.0_dp,0.0_dp]
    x2(2,:) = [2.0_dp,0.0_dp]
    x2(3,:) = [1.0_dp,1.8_dp]
    z = linkage(x2,method='centroid',status=status)
    if (status /= cluster_success .or. z(2,3) >= z(1,3)) &
        error stop 'inverted tree setup'
    ordered = optimal_leaf_ordering(z,x2,status)
    if (status /= cluster_success .or. any(ordered(:,3:4) /= z(:,3:4))) &
        error stop 'inverted tree optimal ordering'
    print *, 'cluster hierarchy ordering/tree regression PASS'

contains
    real(dp) function path_cost(order, condensed) result(cost)
        integer, intent(in) :: order(:) !! One-based leaf sequence.
        real(dp), intent(in) :: condensed(:) !! Pairwise distances in condensed format.
        integer :: n, i, low, high, idx
        n = size(order)
        cost = 0.0_dp
        do i = 1,n-1
            low = min(order(i),order(i+1))
            high = max(order(i),order(i+1))
            idx = (low-1)*n-(low-1)*low/2 + high-low
            cost = cost + condensed(idx)
        end do
    end function path_cost
end program test_cluster_hierarchy_ordering
