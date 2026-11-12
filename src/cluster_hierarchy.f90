! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent Lance-Williams agglomerative clustering implementation.
! The linkage-table convention follows scipy.cluster.hierarchy, except that
! node identifiers are one-based to match native Fortran indexing.
module scifort_cluster_hierarchy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : pdist, spatial_distance_success
    implicit none
    private

    integer, parameter, public :: cluster_success = 0
    integer, parameter, public :: cluster_invalid_input = 1
    public :: linkage, fcluster_distance, fcluster_maxclust
    public :: cophenetic_distances, leaves_list
    public :: inconsistent, maxdists, maxrstat, maxinconsts
    public :: fcluster_inconsistent, fcluster_monocrit, fcluster_maxclust_monocrit
    public :: cut_tree, is_monotonic, num_obs_linkage
    public :: optimal_leaf_ordering, to_tree, leaders, correspond
    public :: fcluster, fclusterdata, cophenetic_correlation
    public :: is_isomorphic, is_valid_linkage, is_valid_im
    public :: to_mlab_linkage, from_mlab_linkage
    public :: dendrogram, cluster_dendrogram
    public :: cluster_tree_node
    public :: single, complete, average, weighted
    public :: centroid, median, ward

    ! Array-based analogue of SciPy's ClusterNode, with one-based node IDs.
    ! A zero parent or child means that the corresponding link is absent.
    type :: cluster_tree_node
        integer :: id = 0
        integer :: left = 0
        integer :: right = 0
        integer :: parent = 0
        integer :: count = 0
        real(dp) :: distance = 0.0_dp
    end type cluster_tree_node

    ! The numeric result of scipy.cluster.hierarchy.dendrogram(no_plot=True).
    ! Each row of icoord and dcoord is one U-shaped linkage segment in
    ! depth-first postorder; leaf indices follow SciFort's one-based scheme.
    ! No color management, figure creation, or truncation is performed.
    type :: cluster_dendrogram
        integer, allocatable :: leaves(:)
        real(dp), allocatable :: icoord(:,:)
        real(dp), allocatable :: dcoord(:,:)
    end type cluster_dendrogram

    ! Temporary dynamic-programming state for optimal leaf ordering. Only
    ! the separator choices remain allocated after a parent is processed.
    type :: ordering_state
        integer, allocatable :: leaves(:)
        real(dp), allocatable :: cost(:,:)
        integer, allocatable :: bridge_left(:,:), bridge_right(:,:)
    end type ordering_state

    interface single
        module procedure single_condensed
        module procedure single_observations
    end interface single

    interface complete
        module procedure complete_condensed
        module procedure complete_observations
    end interface complete

    interface average
        module procedure average_condensed
        module procedure average_observations
    end interface average

    interface weighted
        module procedure weighted_condensed
        module procedure weighted_observations
    end interface weighted

    interface centroid
        module procedure centroid_condensed
        module procedure centroid_observations
    end interface centroid

    interface median
        module procedure median_condensed
        module procedure median_observations
    end interface median

    interface ward
        module procedure ward_condensed
        module procedure ward_observations
    end interface ward

    interface optimal_leaf_ordering
        module procedure optimal_leaf_ordering_condensed
        module procedure optimal_leaf_ordering_observations
    end interface optimal_leaf_ordering

    interface linkage
        module procedure linkage_condensed
        module procedure linkage_observations
    end interface linkage

    interface fcluster
        module procedure fcluster_real
        module procedure fcluster_integer
    end interface fcluster

    interface fclusterdata
        module procedure fclusterdata_real
        module procedure fclusterdata_integer
    end interface fclusterdata

contains

    function single_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='single', status=status, optimal_ordering=optimal_ordering)
    end function single_condensed

    function single_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='single', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function single_observations

    function complete_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='complete', status=status, optimal_ordering=optimal_ordering)
    end function complete_condensed

    function complete_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='complete', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function complete_observations

    function average_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='average', status=status, optimal_ordering=optimal_ordering)
    end function average_condensed

    function average_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='average', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function average_observations

    function weighted_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='weighted', status=status, optimal_ordering=optimal_ordering)
    end function weighted_condensed

    function weighted_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='weighted', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function weighted_observations

    function centroid_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='centroid', status=status, optimal_ordering=optimal_ordering)
    end function centroid_condensed

    function centroid_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='centroid', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function centroid_observations

    function median_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='median', status=status, optimal_ordering=optimal_ordering)
    end function median_condensed

    function median_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='median', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function median_observations

    function ward_condensed(y, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! Condensed upper triangle pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_condensed(y, method='ward', status=status, optimal_ordering=optimal_ordering)
    end function ward_condensed

    function ward_observations(x, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Observations as (samples, features).
        character(len=*), intent(in), optional :: metric !! Distance metric, default euclidean.
        integer, intent(out), optional :: status !! Zero on success.
        logical, intent(in), optional :: optimal_ordering !! Request optimal leaf ordering.
        real(dp), allocatable :: z(:,:)

        z = linkage_observations(x, method='ward', metric=metric, &
            status=status, optimal_ordering=optimal_ordering)
    end function ward_observations

    pure logical function valid_method(method) result(ok)
        character(len=*), intent(in) :: method !! Linkage rule name.
        select case (trim(method))
        case ('single', 'complete', 'average', 'weighted', &
            'centroid', 'median', 'ward')
            ok = .true.
        case default
            ok = .false.
        end select
    end function valid_method

    pure integer function condensed_index(i, j, n) result(k)
        integer, intent(in) :: i !! Smaller one-based sample index.
        integer, intent(in) :: j !! Larger one-based sample index.
        integer, intent(in) :: n !! Total number of samples.
        k = (i-1)*n - (i-1)*i/2 + j-i
    end function condensed_index

    function linkage_observations(x, method, metric, status, optimal_ordering) result(z)
        real(dp), intent(in) :: x(:,:) !! Row-major observations (n, features).
        character(len=*), intent(in), optional :: method !! Linkage rule; default single.
        character(len=*), intent(in), optional :: metric !! Pairwise metric; default euclidean.
        integer, intent(out), optional :: status !! Zero on success, otherwise invalid input.
        logical, intent(in), optional :: optimal_ordering !! Minimize adjacent-leaf distances; default false.
        real(dp), allocatable :: z(:,:)
        real(dp), allocatable :: y(:)
        character(len=:), allocatable :: name
        integer :: code

        name = 'euclidean'
        if (present(metric)) name = metric
        if (size(x,1) < 2 .or. size(x,2) < 1) then
            allocate(z(0,4))
            if (present(status)) status = cluster_invalid_input
            return
        end if
        y = pdist(x, metric=name, status=code)
        if (code /= spatial_distance_success) then
            allocate(z(0,4))
            if (present(status)) status = cluster_invalid_input
            return
        end if
        z = linkage_condensed(y, method, status, optimal_ordering)
    end function linkage_observations

    function linkage_condensed(y, method, status, optimal_ordering) result(z)
        real(dp), intent(in) :: y(:) !! SciPy condensed pairwise distances, n*(n-1)/2.
        character(len=*), intent(in), optional :: method !! Linkage rule; default single.
        integer, intent(out), optional :: status !! Zero on success, otherwise invalid input.
        logical, intent(in), optional :: optimal_ordering !! Minimize adjacent-leaf distances; default false.
        real(dp), allocatable :: z(:,:)
        real(dp), allocatable :: ordered(:,:)
        real(dp), allocatable :: distance(:,:)
        integer, allocatable :: cardinality(:)
        logical, allocatable :: active(:)
        character(len=:), allocatable :: rule
        integer :: n, n_nodes, m, step, i, j, a, b, c, u, k, code
        real(dp) :: best, da, db, dab, updated, na, nb, nc

        if (present(status)) status = cluster_invalid_input
        allocate(z(0,4))
        m = size(y)
        if (m < 1) return
        ! Use int64 arithmetic to prevent integer overflow in triangular counts.
        n = int((1.0_dp + sqrt(1.0_dp+8.0_dp*real(m,dp)))/2.0_dp)
        if (n < 2) return
        if (int(n,int64)*int(n-1,int64)/2_int64 /= int(m,int64)) return
        if (.not. all(ieee_is_finite(y))) return
        if (any(y < 0.0_dp)) return
        rule = 'single'
        if (present(method)) rule = trim(method)
        if (.not. valid_method(rule)) return

        deallocate(z)
        allocate(z(n-1,4))
        n_nodes = 2*n-1
        allocate(distance(n_nodes,n_nodes), cardinality(n_nodes), active(n_nodes))
        distance = 0.0_dp
        cardinality = 0
        cardinality(:n) = 1
        active = .false.
        active(:n) = .true.
        k = 0
        do i=1,n-1
            do j=i+1,n
                k=k+1
                distance(i,j) = y(k)
                distance(j,i) = y(k)
            end do
        end do
        do step=1,n-1
            best = huge(1.0_dp)
            a = 0
            b = 0
            do i=1,n_nodes-1
                if (.not. active(i)) cycle
                do j=i+1,n_nodes
                    if (.not. active(j)) cycle
                    if (a == 0 .or. distance(i,j)<best) then
                        best = distance(i,j)
                        a=i
                        b=j
                    end if
                end do
            end do
            if (a == 0) then
                deallocate(z)
                allocate(z(0,4))
                return
            end if
            u = n+step
            z(step,:) = [real(a,dp), real(b,dp), best, &
                real(cardinality(a)+cardinality(b),dp)]
            dab = best
            na = real(cardinality(a),dp)
            nb = real(cardinality(b),dp)
            do c=1,u-1
                if (.not. active(c) .or. c==a .or. c==b) cycle
                da = distance(a,c)
                db = distance(b,c)
                nc = real(cardinality(c),dp)
                updated = updated_distance(rule, da, db, dab, na, nb, nc)
                if (.not. ieee_is_finite(updated) .or. updated < 0.0_dp) then
                    deallocate(z)
                    allocate(z(0,4))
                    return
                end if
                distance(u,c) = updated
                distance(c,u) = updated
            end do
            cardinality(u) = cardinality(a)+cardinality(b)
            active(a) = .false.
            active(b) = .false.
            active(u) = .true.
        end do
        if (present(optimal_ordering)) then
            if (optimal_ordering) then
                ordered = optimal_leaf_ordering_condensed(z,y,code)
                if (code /= cluster_success) then
                    deallocate(z)
                    allocate(z(0,4))
                    return
                end if
                call move_alloc(ordered,z)
            end if
        end if
        if (present(status)) status = cluster_success
    end function linkage_condensed

    pure function updated_distance(rule, da, db, dab, na, nb, nc) result(d)
        character(len=*), intent(in) :: rule !! Selected linkage rule.
        real(dp), intent(in) :: da !! Existing first-cluster distance.
        real(dp), intent(in) :: db !! Existing second-cluster distance.
        real(dp), intent(in) :: dab !! Distance between clusters being merged.
        real(dp), intent(in) :: na !! Size of first merged cluster.
        real(dp), intent(in) :: nb !! Size of second merged cluster.
        real(dp), intent(in) :: nc !! Size of remaining cluster.
        real(dp) :: d, scale, xa, xb, xab, sq, total
        d = huge(1.0_dp)
        sq = 0.0_dp
        select case (rule)
        case ('single')
            d = min(da,db)
        case ('complete')
            d = max(da,db)
        case ('average')
            d = (na/(na+nb))*da + (nb/(na+nb))*db
        case ('weighted')
            d = 0.5_dp*da + 0.5_dp*db
        case ('centroid', 'median', 'ward')
            scale = max(da,db,dab)
            if (scale == 0.0_dp) then
                d = 0.0_dp
                return
            end if
            xa = da/scale
            xb = db/scale
            xab = dab/scale
            select case (rule)
            case ('centroid')
                total = na+nb
                sq = (na/total)*xa**2 + (nb/total)*xb**2 &
                    - (na*nb/total**2)*xab**2
            case ('median')
                sq = 0.5_dp*(xa**2+xb**2)-0.25_dp*xab**2
            case ('ward')
                total = na+nb+nc
                sq = ((nc+na)/total)*xa**2 + ((nc+nb)/total)*xb**2 &
                    - (nc/total)*xab**2
            end select
            d = scale*sqrt(max(0.0_dp,sq))
        end select
    end function updated_distance

    pure logical function valid_linkage(z) result(ok)
        real(dp), intent(in) :: z(:,:) !! One-based hierarchy table, shape (n-1,4).
        integer, allocatable :: count(:), used(:)
        integer :: n, k, a, b
        ok = .false.
        if (size(z,2)/=4 .or. size(z,1)<1) return
        if (.not. all(ieee_is_finite(z))) return
        n = size(z,1)+1
        allocate(count(2*n-1), used(2*n-1))
        count=0
        used=0
        count(:n)=1
        do k=1,n-1
            if (z(k,1)<1.0_dp .or. z(k,2)<1.0_dp) return
            if (z(k,1)>real(n+k-1,dp) .or. z(k,2)>real(n+k-1,dp)) return
            if (z(k,1)/=real(int(z(k,1)),dp) .or. &
                z(k,2)/=real(int(z(k,2)),dp)) return
            a=int(z(k,1))
            b=int(z(k,2))
            if (a==b .or. z(k,3)<0.0_dp) return
            if (used(a)>0 .or. used(b)>0) return
            used(a)=used(a)+1
            used(b)=used(b)+1
            count(n+k)=count(a)+count(b)
            if (z(k,4)/=real(count(n+k),dp)) return
        end do
        if (count(2*n-1)/=n) return
        ok=.true.
    end function valid_linkage

    function cophenetic_distances(z, status) result(y)
        real(dp), intent(in) :: z(:,:) !! One-based linkage table, (n-1,4).
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: y(:)
        logical, allocatable :: members(:,:)
        integer :: n, k, a, b, i, j, low, high
        if (present(status)) status=cluster_invalid_input
        allocate(y(0))
        if (.not. valid_linkage(z)) return
        n=size(z,1)+1
        deallocate(y)
        allocate(y(n*(n-1)/2), members(n,2*n-1))
        y=0.0_dp
        members=.false.
        do i=1,n
            members(i,i)=.true.
        end do
        do k=1,n-1
            a=int(z(k,1))
            b=int(z(k,2))
            do i=1,n
                if (.not. members(i,a)) cycle
                do j=1,n
                    if (.not. members(j,b)) cycle
                    low=min(i,j)
                    high=max(i,j)
                    y(condensed_index(low,high,n))=z(k,3)
                end do
            end do
            members(:,n+k) = members(:,a) .or. members(:,b)
        end do
        if (present(status)) status=cluster_success
    end function cophenetic_distances

    function cophenetic_correlation(z, y, status) result(correlation)
        ! Pearson correlation between original condensed dissimilarities
        ! and the ultrametric/cophenetic values of a valid linkage tree.
        ! Scale each vector independently to avoid overflow in the sum of
        ! products; a zero-variance input has an undefined correlation.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix, shape (n-1,4).
        real(dp), intent(in) :: y(:) !! Original condensed pairwise distances, length n*(n-1)/2.
        integer, intent(out), optional :: status !! Zero on valid inputs, one otherwise.
        real(dp) :: correlation
        real(dp), allocatable :: c(:), yc(:), cc(:)
        real(dp) :: norm_y, norm_c, scale_y, scale_c
        integer :: n, code

        correlation = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = cluster_invalid_input
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        if (int(size(y),int64) /= int(n,int64)*int(n-1,int64)/2_int64) return
        if (.not. all(ieee_is_finite(y))) return
        if (any(y < 0.0_dp)) return
        c = cophenetic_distances(z,code)
        if (code /= cluster_success) return

        scale_y = maxval(abs(y))
        scale_c = maxval(abs(c))
        if (scale_y == 0.0_dp .or. scale_c == 0.0_dp) then
            ! Pearson r is undefined when either vector is constant.
            if (present(status)) status = cluster_success
            return
        end if
        yc = y/scale_y
        cc = c/scale_c
        yc = yc - sum(yc)/real(size(yc),dp)
        cc = cc - sum(cc)/real(size(cc),dp)
        norm_y = sqrt(sum(yc**2))
        norm_c = sqrt(sum(cc**2))
        if (norm_y > 0.0_dp .and. norm_c > 0.0_dp) then
            correlation = sum((yc/norm_y)*(cc/norm_c))
            correlation = max(-1.0_dp,min(1.0_dp,correlation))
        end if
        if (present(status)) status = cluster_success
    end function cophenetic_correlation

    function leaves_list(z, status) result(leaves)
        real(dp), intent(in) :: z(:,:) !! One-based linkage table, (n-1,4).
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: leaves(:), stack(:)
        integer :: n, top, node, pos
        if (present(status)) status=cluster_invalid_input
        allocate(leaves(0))
        if (.not. valid_linkage(z)) return
        n=size(z,1)+1
        deallocate(leaves)
        allocate(leaves(n), stack(2*n))
        top=1
        stack(top)=2*n-1
        pos=0
        do while(top>0)
            node=stack(top)
            top=top-1
            if (node<=n) then
                pos=pos+1
                leaves(pos)=node
            else
                top=top+1
                stack(top)=int(z(node-n,2))
                top=top+1
                stack(top)=int(z(node-n,1))
            end if
        end do
        if (present(status)) status=cluster_success
    end function leaves_list

    function fcluster_distance(z, threshold, status) result(labels)
        real(dp), intent(in) :: z(:,:) !! One-based linkage table, (n-1,4).
        real(dp), intent(in) :: threshold !! Maximum monotone linkage height.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: labels(:)
        integer, allocatable :: parent(:), rep(:)
        real(dp), allocatable :: height(:)
        integer :: n, i, k, a, b, root, new_label
        if (present(status)) status=cluster_invalid_input
        allocate(labels(0))
        if (.not. valid_linkage(z)) return
        if (.not. ieee_is_finite(threshold)) return
        n=size(z,1)+1
        deallocate(labels)
        allocate(labels(n), parent(n), rep(2*n-1), height(2*n-1))
        parent=[(i,i=1,n)]
        height=0.0_dp
        rep=0
        rep(:n)=[(i,i=1,n)]
        do k=1,n-1
            a=int(z(k,1))
            b=int(z(k,2))
            rep(n+k)=rep(a)
            height(n+k)=max(z(k,3),height(a),height(b))
            if (height(n+k)<=threshold) then
                call join_roots(parent,rep(a),rep(b))
            end if
        end do
        labels=0
        new_label=0
        do i=1,n
            root=find_root(parent,i)
            if (labels(root)==0) then
                new_label=new_label+1
                labels(root)=new_label
            end if
            labels(i)=labels(root)
        end do
        if (present(status)) status=cluster_success
    end function fcluster_distance

    function fcluster_maxclust(z, max_clusters, status) result(labels)
        real(dp), intent(in) :: z(:,:) !! One-based linkage table, (n-1,4).
        integer, intent(in) :: max_clusters !! Target maximum number of flat clusters >=1.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: labels(:)
        real(dp), allocatable :: heights(:)
        integer :: n, k, a, b, j, dummy
        real(dp) :: h
        if (present(status)) status=cluster_invalid_input
        allocate(labels(0))
        if (.not. valid_linkage(z)) return
        if (max_clusters<1) return
        n=size(z,1)+1
        if (max_clusters>=n) then
            deallocate(labels)
            allocate(labels(n))
            labels=[(k,k=1,n)]
            if (present(status)) status=cluster_success
            return
        end if
        allocate(heights(2*n-1))
        heights=0.0_dp
        do k=1,n-1
            a=int(z(k,1))
            b=int(z(k,2))
            heights(n+k)=max(z(k,3),heights(a),heights(b))
        end do
        ! Selection sort is adequate for the O(n^3) linkage constructor.
        do k=1,n-2
            do j=k+1,n-1
                if (heights(n+j)<heights(n+k)) then
                    h=heights(n+k)
                    heights(n+k)=heights(n+j)
                    heights(n+j)=h
                end if
            end do
        end do
        labels=fcluster_distance(z,heights(n+(n-max_clusters)),status=dummy)
        if (present(status)) status=dummy
    end function fcluster_maxclust

    function fcluster_real(z, threshold, criterion, depth, r, monocrit, status) result(labels)
        ! SciPy-style dispatcher for five flat-clustering criteria. The
        ! underlying routines remain independently callable for typed
        ! maximum-cluster counts and explicit monotone criterion vectors.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        real(dp), intent(in) :: threshold !! Distance/statistic cutoff or maximum cluster count.
        character(len=*), intent(in), optional :: criterion !! Cut rule, default 'inconsistent'.
        integer, intent(in), optional :: depth !! Inconsistency depth, default 2.
        real(dp), intent(in), optional :: r(:,:) !! Optional (n-1,4) inconsistency statistics.
        real(dp), intent(in), optional :: monocrit(:) !! Monotone node criterion for monocrit rules.
        integer, intent(out), optional :: status !! Zero on success, one on invalid input.
        integer, allocatable :: labels(:)
        character(len=:), allocatable :: rule
        integer :: code, count_clusters, n

        if (present(status)) status = cluster_invalid_input
        allocate(labels(0))
        if (.not. valid_linkage(z)) return
        if (.not. ieee_is_finite(threshold)) return
        rule = 'inconsistent'
        if (present(criterion)) rule = trim(criterion)
        select case(rule)
        case ('inconsistent')
            labels = fcluster_inconsistent(z,threshold,code,depth,r)
        case ('distance')
            labels = fcluster_distance(z,threshold,code)
        case ('maxclust','maxclust_monocrit')
            ! SciPy truncates positive floating-point counts toward zero.
            ! Reject nonpositive counts rather than silently returning n
            ! singletons; large counts are capped before conversion to int.
            if (threshold < 1.0_dp) return
            n = size(z,1)+1
            if (threshold >= real(n,dp)) then
                count_clusters = n
            else
                count_clusters = int(threshold)
            end if
            if (rule == 'maxclust') then
                labels = fcluster_maxclust(z,count_clusters,code)
            else
                if (.not. present(monocrit)) return
                labels = fcluster_maxclust_monocrit(z,monocrit,count_clusters,code)
            end if
        case ('monocrit')
            if (.not. present(monocrit)) return
            labels = fcluster_monocrit(z,monocrit,threshold,code)
        case default
            return
        end select
        if (present(status)) status = code
    end function fcluster_real

    function fclusterdata_real(x, threshold, criterion, metric, depth, method, r, status) result(labels)
        ! Observation-to-labels convenience API: pdist -> linkage -> fcluster.
        ! Like scipy.cluster.hierarchy.fclusterdata, x is a rank-two matrix.
        real(dp), intent(in) :: x(:,:) !! Observation rows, shape (n,features), n>=2.
        real(dp), intent(in) :: threshold !! Cut threshold or requested max cluster count.
        character(len=*), intent(in), optional :: criterion !! Cut rule, default inconsistent.
        character(len=*), intent(in), optional :: metric !! Pairwise metric, default euclidean.
        integer, intent(in), optional :: depth !! Inconsistency depth, default 2.
        character(len=*), intent(in), optional :: method !! Linkage rule, default single.
        real(dp), intent(in), optional :: r(:,:) !! Precomputed inconsistency matrix for this linkage.
        integer, intent(out), optional :: status !! Zero on success, one on invalid input.
        integer, allocatable :: labels(:)
        real(dp), allocatable :: z(:,:)
        integer :: code

        if (present(status)) status = cluster_invalid_input
        allocate(labels(0))
        if (size(x,1) < 2 .or. size(x,2) < 1) return
        z = linkage(x,method=method,metric=metric,status=code)
        if (code /= cluster_success) return
        labels = fcluster(z,threshold,criterion=criterion,depth=depth,r=r,status=code)
        if (present(status)) status = code
    end function fclusterdata_real

    function fcluster_integer(z, threshold, criterion, depth, r, monocrit, status) result(labels)
        ! Integer overload for natural Fortran maxclust counts and exact
        ! integer thresholds. All validation occurs in the real overload.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(in) :: threshold !! Cutoff or requested maximum number of clusters.
        character(len=*), intent(in), optional :: criterion !! Cut rule, default inconsistent.
        integer, intent(in), optional :: depth !! Inconsistency depth, default 2.
        real(dp), intent(in), optional :: r(:,:) !! Optional precomputed inconsistency matrix.
        real(dp), intent(in), optional :: monocrit(:) !! Monotone per-node criterion.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: labels(:)

        labels = fcluster_real(z,real(threshold,dp),criterion,depth,r,monocrit,status)
    end function fcluster_integer

    function fclusterdata_integer(x, threshold, criterion, metric, depth, method, r, status) result(labels)
        ! Accept ordinary integer cluster-count thresholds without requiring
        ! explicit REAL(...,dp) conversion at call sites.
        real(dp), intent(in) :: x(:,:) !! Observation matrix, shape (n,features).
        integer, intent(in) :: threshold !! Cutoff or maximum number of clusters.
        character(len=*), intent(in), optional :: criterion !! Cut rule, default inconsistent.
        character(len=*), intent(in), optional :: metric !! Pairwise metric, default euclidean.
        integer, intent(in), optional :: depth !! Inconsistency depth, default 2.
        character(len=*), intent(in), optional :: method !! Linkage rule, default single.
        real(dp), intent(in), optional :: r(:,:) !! Optional precomputed inconsistency matrix.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: labels(:)

        labels = fclusterdata_real(x,real(threshold,dp),criterion,metric,depth,method,r,status)
    end function fclusterdata_integer

    pure integer function find_root(parent, i) result(root)
        integer, intent(in) :: parent(:) !! Union-find parent links.
        integer, intent(in) :: i !! A leaf identifier.
        root=i
        do while(parent(root)/=root)
            root=parent(root)
        end do
    end function find_root

    pure subroutine join_roots(parent, a, b)
        integer, intent(inout) :: parent(:) !! Union-find parents, updated in place.
        integer, intent(in) :: a !! First leaf representative.
        integer, intent(in) :: b !! Second leaf representative.
        integer :: ra, rb
        ra=find_root(parent,a)
        rb=find_root(parent,b)
        if (ra/=rb) parent(rb)=ra
    end subroutine join_roots

    function inconsistent(z, depth, status) result(r)
        ! Inconsistency statistics for each non-singleton node.  The current
        ! merge and at most depth-1 generations of descendant merges are used.
        ! Sample standard deviations (n-1 denominator) match SciPy.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage, (n-1,4).
        integer, intent(in), optional :: depth !! Levels including current merge; default 2, >= 0.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: r(:,:)
        real(dp), allocatable :: values(:)
        integer, allocatable :: nodes(:), levels(:)
        integer :: n, d, k, top, node, level, count_links, child, i
        real(dp) :: scale, mean_value, delta, m2, sd

        if (present(status)) status = cluster_invalid_input
        allocate(r(0,4))
        if (.not. valid_linkage(z)) return
        d = 2
        if (present(depth)) d = depth
        if (d < 0) return
        d = max(1,d)  ! SciPy's depth=0 is equivalent to depth=1.
        n = size(z,1)+1
        deallocate(r)
        allocate(r(n-1,4),values(n-1),nodes(2*n-1),levels(2*n-1))
        r = 0.0_dp
        do k = 1,n-1
            top = 1
            nodes(top) = n+k
            levels(top) = 1
            count_links = 0
            do while (top > 0)
                node = nodes(top)
                level = levels(top)
                top = top-1
                if (node <= n) cycle
                count_links = count_links+1
                values(count_links) = z(node-n,3)
                if (level >= d) cycle
                do i = 1,2
                    child = int(z(node-n,i))
                    if (child <= n) cycle
                    top = top+1
                    nodes(top) = child
                    levels(top) = level+1
                end do
            end do
            scale = maxval(values(:count_links))
            mean_value = 0.0_dp
            m2 = 0.0_dp
            if (scale > 0.0_dp) then
                do i = 1,count_links
                    delta = values(i)/scale-mean_value
                    mean_value = mean_value+delta/real(i,dp)
                    m2 = m2+delta*(values(i)/scale-mean_value)
                end do
            end if
            sd = 0.0_dp
            if (count_links > 1) sd = sqrt(max(0.0_dp,m2/real(count_links-1,dp)))
            r(k,1) = scale*mean_value
            r(k,2) = scale*sd
            r(k,3) = real(count_links,dp)
            if (sd > 0.0_dp) then
                r(k,4) = (z(k,3)/scale-mean_value)/sd
            else
                r(k,4) = 0.0_dp
            end if
        end do
        if (present(status)) status = cluster_success
    end function inconsistent

    function maxdists(z, status) result(maximum)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: maximum(:)
        real(dp), allocatable :: descendant(:)
        integer :: n, k, a, b

        if (present(status)) status = cluster_invalid_input
        allocate(maximum(0))
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        allocate(descendant(2*n-1))
        descendant = 0.0_dp
        deallocate(maximum)
        allocate(maximum(n-1))
        do k = 1,n-1
            a = int(z(k,1))
            b = int(z(k,2))
            descendant(n+k) = max(z(k,3),descendant(a),descendant(b))
            maximum(k) = descendant(n+k)
        end do
        if (present(status)) status = cluster_success
    end function maxdists

    function maxrstat(z, r, column, status) result(maximum)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        real(dp), intent(in) :: r(:,:) !! Per-node statistics, shape (n-1,4).
        integer, intent(in) :: column !! One-based column of r, from 1 through 4.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: maximum(:)
        real(dp), allocatable :: descendant(:)
        integer :: n, k, a, b

        if (present(status)) status = cluster_invalid_input
        allocate(maximum(0))
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        if (any(shape(r) /= [n-1,4])) return
        if (column < 1 .or. column > 4) return
        if (.not. all(ieee_is_finite(r))) return
        allocate(descendant(2*n-1))
        descendant = -huge(1.0_dp)
        deallocate(maximum)
        allocate(maximum(n-1))
        do k = 1,n-1
            a = int(z(k,1))
            b = int(z(k,2))
            descendant(n+k) = max(r(k,column),descendant(a),descendant(b))
            maximum(k) = descendant(n+k)
        end do
        if (present(status)) status = cluster_success
    end function maxrstat

    function maxinconsts(z, r, status) result(maximum)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        real(dp), intent(in) :: r(:,:) !! Inconsistency statistics from inconsistent().
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: maximum(:)

        maximum = maxrstat(z,r,4,status)
    end function maxinconsts

    pure logical function valid_monocrit(z, criterion) result(ok)
        real(dp), intent(in) :: z(:,:) !! Already-validated linkage matrix.
        real(dp), intent(in) :: criterion(:) !! Candidate monotone node criterion.
        integer :: n, k, j, child

        n = size(z,1)+1
        ok = .false.
        if (size(criterion) /= n-1) return
        if (.not. all(ieee_is_finite(criterion))) return
        do k = 1,n-1
            do j = 1,2
                child = int(z(k,j))
                if (child <= n) cycle
                if (criterion(k) < criterion(child-n)) return
            end do
        end do
        ok = .true.
    end function valid_monocrit

    function fcluster_monocrit(z, criterion, threshold, status) result(labels)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        real(dp), intent(in) :: criterion(:) !! Monotone per-merge criterion, length n-1.
        real(dp), intent(in) :: threshold !! Cut criterion; merge if <= threshold.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: labels(:), parent(:), rep(:), groups(:)
        integer :: n, k, i, a, b, root, next_label

        if (present(status)) status = cluster_invalid_input
        allocate(labels(0))
        if (.not. valid_linkage(z)) return
        if (.not. valid_monocrit(z,criterion)) return
        if (.not. ieee_is_finite(threshold)) return
        n = size(z,1)+1
        allocate(parent(n),rep(2*n-1),groups(n))
        parent = [(i,i=1,n)]
        rep = 0
        rep(:n) = [(i,i=1,n)]
        do k = 1,n-1
            a = int(z(k,1))
            b = int(z(k,2))
            rep(n+k) = rep(a)
            if (criterion(k) <= threshold) call join_roots(parent,rep(a),rep(b))
        end do
        deallocate(labels)
        allocate(labels(n))
        groups = 0
        next_label = 0
        do i = 1,n
            root = find_root(parent,i)
            if (groups(root) == 0) then
                next_label = next_label+1
                groups(root) = next_label
            end if
            labels(i) = groups(root)
        end do
        if (present(status)) status = cluster_success
    end function fcluster_monocrit

    function fcluster_maxclust_monocrit(z, criterion, max_clusters, status) result(labels)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        real(dp), intent(in) :: criterion(:) !! Monotone per-merge criterion, length n-1.
        integer, intent(in) :: max_clusters !! Maximum number of clusters, >= 1.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: labels(:)
        real(dp), allocatable :: sorted(:)
        real(dp) :: temp
        integer :: n, i, j

        if (present(status)) status = cluster_invalid_input
        allocate(labels(0))
        if (.not. valid_linkage(z)) return
        if (.not. valid_monocrit(z,criterion)) return
        if (max_clusters < 1) return
        n = size(z,1)+1
        if (max_clusters >= n) then
            deallocate(labels)
            allocate(labels(n))
            labels = [(i,i=1,n)]
            if (present(status)) status = cluster_success
            return
        end if
        sorted = criterion
        do i = 1,n-2
            do j = i+1,n-1
                if (sorted(j) < sorted(i)) then
                    temp = sorted(i)
                    sorted(i) = sorted(j)
                    sorted(j) = temp
                end if
            end do
        end do
        labels = fcluster_monocrit(z,criterion,sorted(n-max_clusters),status)
    end function fcluster_maxclust_monocrit

    function fcluster_inconsistent(z, threshold, status, depth, r) result(labels)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        real(dp), intent(in) :: threshold !! Maximum inconsistency coefficient.
        integer, intent(out), optional :: status !! Zero on success.
        integer, intent(in), optional :: depth !! Inconsistency depth; default 2.
        real(dp), intent(in), optional :: r(:,:) !! Optional precomputed inconsistency matrix.
        integer, allocatable :: labels(:)
        real(dp), allocatable :: stats(:,:), crit(:)
        integer :: code

        if (present(status)) status = cluster_invalid_input
        allocate(labels(0))
        if (.not. valid_linkage(z)) return
        if (present(r)) then
            stats = r
        else
            stats = inconsistent(z,depth,code)
            if (code /= cluster_success) return
        end if
        crit = maxinconsts(z,stats,code)
        if (code /= cluster_success) return
        labels = fcluster_monocrit(z,crit,threshold,code)
        if (present(status)) status = code
    end function fcluster_inconsistent

    logical function is_monotonic(z, status) result(ok)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(out), optional :: status !! Zero on success.

        ok = .false.
        if (present(status)) status = cluster_invalid_input
        if (.not. valid_linkage(z)) return
        ok = all(z(2:size(z,1),3) >= z(1:size(z,1)-1,3))
        if (present(status)) status = cluster_success
    end function is_monotonic

    integer function num_obs_linkage(z, status) result(n)
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(out), optional :: status !! Zero on success.

        n = 0
        if (present(status)) status = cluster_invalid_input
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        if (present(status)) status = cluster_success
    end function num_obs_linkage

    function cut_tree(z, n_clusters, height, status) result(cuts)
        ! Full cut tree or selected exact cluster counts / height thresholds.
        ! For heights, SciPy uses a *strict* height comparison; equality does
        ! not merge that node.  One-based cluster labels are used throughout.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(in), optional :: n_clusters(:) !! Requested counts from 1 through n.
        real(dp), intent(in), optional :: height(:) !! Merge heights, only for monotone linkage.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: cuts(:,:)
        integer, allocatable :: parent(:), rep(:), groups(:)
        integer :: n, ncols, col, steps, i, k, a, b, root, next_label
        logical :: monotone

        if (present(status)) status = cluster_invalid_input
        allocate(cuts(0,0))
        if (.not. valid_linkage(z)) return
        if (present(n_clusters) .and. present(height)) return
        n = size(z,1)+1
        ncols = n
        if (present(n_clusters)) then
            ncols = size(n_clusters)
            if (any(n_clusters < 1) .or. any(n_clusters > n)) return
        end if
        if (present(height)) then
            ncols = size(height)
            if (.not. all(ieee_is_finite(height))) return
            monotone = is_monotonic(z)
            if (.not. monotone) return
        end if
        deallocate(cuts)
        allocate(cuts(n,ncols),parent(n),rep(2*n-1),groups(n))
        do col = 1,ncols
            steps = col-1
            if (present(n_clusters)) steps = n-n_clusters(col)
            if (present(height)) steps = count(z(:,3) < height(col))
            parent = [(i,i=1,n)]
            rep = 0
            rep(:n) = [(i,i=1,n)]
            do k = 1,steps
                a = int(z(k,1))
                b = int(z(k,2))
                rep(n+k) = rep(a)
                call join_roots(parent,rep(a),rep(b))
            end do
            groups = 0
            next_label = 0
            do i = 1,n
                root = find_root(parent,i)
                if (groups(root) == 0) then
                    next_label = next_label+1
                    groups(root) = next_label
                end if
                cuts(i,col) = groups(root)
            end do
        end do
        if (present(status)) status = cluster_success
    end function cut_tree


    function optimal_leaf_ordering_observations(z, x, status, metric) result(ordered)
        ! Minimize adjacent-leaf dissimilarities subject to a fixed tree.
        real(dp), intent(in) :: z(:,:) !! One-based linkage, shape (n-1,4).
        real(dp), intent(in) :: x(:,:) !! Observation rows, shape (n,features).
        integer, intent(out), optional :: status !! Zero on success.
        character(len=*), intent(in), optional :: metric !! Distance name; default euclidean.
        real(dp), allocatable :: ordered(:,:)
        real(dp), allocatable :: y(:)
        integer :: n, code
        character(len=:), allocatable :: metric_name

        if (present(status)) status = cluster_invalid_input
        allocate(ordered(0,4))
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        if (size(x,1) /= n .or. size(x,2) < 1) return
        metric_name = 'euclidean'
        if (present(metric)) metric_name = trim(metric)
        y = pdist(x,metric=metric_name,status=code)
        if (code /= spatial_distance_success) return
        deallocate(ordered)
        ordered = optimal_leaf_ordering_condensed(z,y,code)
        if (present(status)) status = code
    end function optimal_leaf_ordering_observations

    function optimal_leaf_ordering_condensed(z, y, status) result(ordered)
        ! For each subtree, compute the minimum adjacent-distance sum for
        ! each possible ordered pair of extreme leaves. Reversing either
        ! child is permitted, but no other rearrangement is performed.
        ! Separator choices use O(n^2) memory in total: each unordered
        ! leaf pair is separated by its least common ancestor exactly once.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage, shape (n-1,4).
        real(dp), intent(in) :: y(:) !! Original condensed pairwise distances.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: ordered(:,:)
        real(dp), allocatable :: dist(:,:)
        type(ordering_state), allocatable :: state(:)
        integer :: n, m, k, a, b, sa, sb, size_parent
        integer :: i, j, p, q, best_p, best_q, best_k, best_l
        integer :: left_leaf, right_leaf
        real(dp) :: scale, candidate, best, best_total

        if (present(status)) status = cluster_invalid_input
        allocate(ordered(0,4))
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        if (int(size(y),int64) /= int(n,int64)*int(n-1,int64)/2_int64) return
        if (.not. all(ieee_is_finite(y))) return
        if (any(y < 0.0_dp)) return
        scale = maxval(y)
        allocate(dist(n,n),state(2*n-1))
        dist = 0.0_dp
        m = 0
        do i = 1,n-1
            do j = i+1,n
                m = m+1
                if (scale > 0.0_dp) then
                    dist(i,j) = y(m)/scale
                    dist(j,i) = dist(i,j)
                end if
            end do
        end do
        do i = 1,n
            state(i)%leaves = [i]
            allocate(state(i)%cost(1,1))
            state(i)%cost(1,1) = 0.0_dp
        end do
        do k = 1,n-1
            a = int(z(k,1))
            b = int(z(k,2))
            sa = size(state(a)%leaves)
            sb = size(state(b)%leaves)
            size_parent = sa+sb
            state(n+k)%leaves = [state(a)%leaves,state(b)%leaves]
            allocate(state(n+k)%cost(size_parent,size_parent))
            allocate(state(n+k)%bridge_left(sa,sb),state(n+k)%bridge_right(sa,sb))
            state(n+k)%cost = huge(1.0_dp)
            do i = 1,sa
                do j = 1,sb
                    best = huge(1.0_dp)
                    best_k = 0
                    best_l = 0
                    do p = 1,sa
                        if (state(a)%cost(i,p) == huge(1.0_dp)) cycle
                        do q = 1,sb
                            if (state(b)%cost(q,j) == huge(1.0_dp)) cycle
                            candidate = state(a)%cost(i,p) + &
                                dist(state(a)%leaves(p),state(b)%leaves(q)) + &
                                state(b)%cost(q,j)
                            if (candidate < best) then
                                best = candidate
                                best_k = p
                                best_l = q
                            end if
                        end do
                    end do
                    if (best_k == 0) return
                    state(n+k)%cost(i,sa+j) = best
                    state(n+k)%cost(sa+j,i) = best
                    state(n+k)%bridge_left(i,j) = state(a)%leaves(best_k)
                    state(n+k)%bridge_right(i,j) = state(b)%leaves(best_l)
                end do
            end do
            ! Descendants need only their separator choices for traceback.
            deallocate(state(a)%cost,state(b)%cost)
        end do
        best_total = huge(1.0_dp)
        best_p = 0
        best_q = 0
        do i = 1,n
            do j = 1,n
                if (state(2*n-1)%cost(i,j) < best_total) then
                    best_total = state(2*n-1)%cost(i,j)
                    best_p = i
                    best_q = j
                end if
            end do
        end do
        if (best_p == 0) return
        left_leaf = state(2*n-1)%leaves(best_p)
        right_leaf = state(2*n-1)%leaves(best_q)
        deallocate(ordered)
        ordered = z
        call orient_ordered_tree(ordered,state,2*n-1,left_leaf,right_leaf,n)
        if (present(status)) status = cluster_success
    end function optimal_leaf_ordering_condensed

    recursive subroutine orient_ordered_tree(z, state, node, first, last, n)
        ! Trace back a minimum-cost endpoint pair, modifying only column
        ! 1 and column 2 (the two child links) of the parent linkage row.
        real(dp), intent(inout) :: z(:,:) !! Reoriented tree, updated in place.
        type(ordering_state), intent(in) :: state(:) !! Stored subtree separators.
        integer, intent(in) :: node !! Current one-based node ID.
        integer, intent(in) :: first !! Desired first leaf of this subtree.
        integer, intent(in) :: last !! Desired last leaf of this subtree.
        integer, intent(in) :: n !! Original observation count.
        integer :: a, b, sa, i, j, left_end, right_start
        real(dp) :: temp

        if (node <= n) return
        a = int(z(node-n,1))
        b = int(z(node-n,2))
        sa = size(state(a)%leaves)
        i = findloc(state(node)%leaves,first,dim=1)
        j = findloc(state(node)%leaves,last,dim=1)
        if (i <= sa) then
            left_end = state(node)%bridge_left(i,j-sa)
            right_start = state(node)%bridge_right(i,j-sa)
            call orient_ordered_tree(z,state,a,first,left_end,n)
            call orient_ordered_tree(z,state,b,right_start,last,n)
        else
            left_end = state(node)%bridge_left(j,i-sa)
            right_start = state(node)%bridge_right(j,i-sa)
            temp = z(node-n,1)
            z(node-n,1) = z(node-n,2)
            z(node-n,2) = temp
            call orient_ordered_tree(z,state,b,first,right_start,n)
            call orient_ordered_tree(z,state,a,left_end,last,n)
        end if
    end subroutine orient_ordered_tree

    function to_tree(z, status) result(nodes)
        ! Materialize an indexable tree with stable one-based node IDs.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(out), optional :: status !! Zero on success.
        type(cluster_tree_node), allocatable :: nodes(:)
        integer :: n, i, k, a, b, parent

        if (present(status)) status = cluster_invalid_input
        allocate(nodes(0))
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        deallocate(nodes)
        allocate(nodes(2*n-1))
        do i = 1,2*n-1
            nodes(i)%id = i
            nodes(i)%count = 1
        end do
        do k = 1,n-1
            a = int(z(k,1))
            b = int(z(k,2))
            parent = n+k
            nodes(parent)%left = a
            nodes(parent)%right = b
            nodes(parent)%distance = z(k,3)
            nodes(parent)%count = int(z(k,4))
            nodes(a)%parent = parent
            nodes(b)%parent = parent
        end do
        if (present(status)) status = cluster_success
    end function to_tree

    function leaders(z, labels, status) result(pairs)
        ! Return one [leader_node_id, flat_cluster_label] row per cluster.
        ! The labels must describe complete, disjoint linkage subtrees.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(in) :: labels(:) !! One integer cluster label per leaf.
        integer, intent(out), optional :: status !! Zero on success.
        integer, allocatable :: pairs(:,:)
        integer, allocatable :: same(:), parent(:), candidates(:,:)
        logical, allocatable :: uniform(:)
        integer :: n, i, k, a, b, node, nfound

        if (present(status)) status = cluster_invalid_input
        allocate(pairs(0,2))
        if (.not. valid_linkage(z)) return
        n = size(z,1)+1
        if (size(labels) /= n) return
        allocate(same(2*n-1),parent(2*n-1),uniform(2*n-1),candidates(n,2))
        same = 0
        parent = 0
        uniform = .false.
        uniform(:n) = .true.
        same(:n) = labels
        do k = 1,n-1
            node = n+k
            a = int(z(k,1))
            b = int(z(k,2))
            parent(a) = node
            parent(b) = node
            if (uniform(a) .and. uniform(b)) then
                if (same(a) == same(b)) then
                    same(node) = same(a)
                    uniform(node) = .true.
                end if
            end if
        end do
        nfound = 0
        do node = 1,2*n-1
            if (.not. uniform(node)) cycle
            if (parent(node) > 0) then
                if (uniform(parent(node))) cycle
            end if
            if (nfound > 0) then
                if (any(candidates(:nfound,2) == same(node))) return
            end if
            nfound = nfound+1
            candidates(nfound,:) = [node,same(node)]
        end do
        deallocate(pairs)
        allocate(pairs(nfound,2))
        pairs = candidates(:nfound,:)
        if (present(status)) status = cluster_success
    end function leaders

    logical function correspond(z, y, status) result(ok)
        ! Like SciPy's correspond, only check the number of observations;
        ! equality of the actual linkage distances is *not* asserted.
        real(dp), intent(in) :: z(:,:) !! One-based linkage matrix.
        real(dp), intent(in) :: y(:) !! Condensed pairwise distances.
        integer, intent(out), optional :: status !! Zero on valid inputs.
        integer :: n

        ok = .false.
        if (present(status)) status = cluster_invalid_input
        if (.not. valid_linkage(z)) return
        if (.not. all(ieee_is_finite(y))) return
        if (any(y < 0.0_dp)) return
        n = size(z,1)+1
        ok = int(size(y),int64) == int(n,int64)*int(n-1,int64)/2_int64
        if (present(status)) status = cluster_success
    end function correspond


    pure logical function is_valid_linkage(z) result(ok)
        ! Check topology, node numbers, distances, and exact leaf counts.
        real(dp), intent(in) :: z(:,:) !! One-based linkage matrix, shape (n-1,4).
        ok = valid_linkage(z)
    end function is_valid_linkage

    pure logical function is_valid_im(r) result(ok)
        ! Numeric inconsistency-matrix validator. Accepts nonnegative
        ! link-height means, standard deviations, and link counts.
        ! Unlike SciPy's permissive 1.17 validator, reject nonfinite values.
        real(dp), intent(in) :: r(:,:) !! Inconsistency statistics, shape (n-1,4).
        ok = .false.
        if (size(r,1) < 1 .or. size(r,2) /= 4) return
        if (.not. all(ieee_is_finite(r))) return
        if (any(r(:,1) < 0.0_dp)) return
        if (any(r(:,2) < 0.0_dp)) return
        if (any(r(:,3) < 0.0_dp)) return
        ok = .true.
    end function is_valid_im

    logical function is_isomorphic(labels1, labels2, status) result(equal)
        ! Equivalence of flat partitions independent of numeric label IDs.
        ! Sort indices by each label array in turn to establish a bijection.
        ! Time O(n log n), memory O(n); labels may be any Fortran integers.
        integer, intent(in) :: labels1(:) !! First flat-cluster assignment.
        integer, intent(in) :: labels2(:) !! Second flat-cluster assignment.
        integer, intent(out), optional :: status !! Zero if lengths agree.
        integer, allocatable :: order(:), work(:)
        integer :: n, i

        equal = .false.
        if (present(status)) status = cluster_invalid_input
        n = size(labels1)
        if (n /= size(labels2)) return
        if (present(status)) status = cluster_success
        allocate(order(n),work(n))
        order = [(i,i=1,n)]
        call sort_label_indices(labels1, order, work)
        do i = 2,n
            if (labels1(order(i)) == labels1(order(i-1)) .and. &
                labels2(order(i)) /= labels2(order(i-1))) return
        end do
        order = [(i,i=1,n)]
        call sort_label_indices(labels2, order, work)
        do i = 2,n
            if (labels2(order(i)) == labels2(order(i-1)) .and. &
                labels1(order(i)) /= labels1(order(i-1))) return
        end do
        equal = .true.
    end function is_isomorphic

    pure subroutine sort_label_indices(keys, order, work)
        ! Stable bottom-up mergesort; always terminates for 32-bit n.
        integer, intent(in) :: keys(:) !! Keys used for sorting.
        integer, intent(inout) :: order(:) !! Index permutation to sort.
        integer, intent(out) :: work(:) !! Workspace, same size as order.
        integer :: n, width, base, mid, last, i, j, k

        n = size(order)
        width = 1
        do while (width < n)
            base = 1
            do while (base <= n)
                mid = base + min(width,n-base+1)
                last = mid + min(width,n-mid+1) - 1
                i = base
                j = mid
                do k = base,last
                    if (i < mid) then
                        if (j > last) then
                            work(k) = order(i)
                            i = i+1
                        else if (keys(order(i)) <= keys(order(j))) then
                            work(k) = order(i)
                            i = i+1
                        else
                            work(k) = order(j)
                            j = j+1
                        end if
                    else
                        work(k) = order(j)
                        j = j+1
                    end if
                end do
                base = last+1
            end do
            order = work
            if (width > n/2) exit
            width = 2*width
        end do
    end subroutine sort_label_indices

    function to_mlab_linkage(z, status) result(mlab)
        ! The native SciFort linkage already uses MATLAB's one-based IDs;
        ! converting to MATLAB therefore just discards the leaf counts.
        real(dp), intent(in) :: z(:,:) !! Native linkage, shape (n-1,4).
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: mlab(:,:)

        if (present(status)) status = cluster_invalid_input
        allocate(mlab(0,3))
        if (.not. valid_linkage(z)) return
        mlab = z(:,1:3)
        if (present(status)) status = cluster_success
    end function to_mlab_linkage

    function from_mlab_linkage(mlab, status) result(z)
        ! Reconstruct fourth-column cluster cardinalities from MATLAB's
        ! three-column, one-based linkage representation.
        real(dp), intent(in) :: mlab(:,:) !! MATLAB linkage, shape (n-1,3).
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: z(:,:)
        integer, allocatable :: counts(:)
        integer :: n, k, a, b

        if (present(status)) status = cluster_invalid_input
        allocate(z(0,4))
        if (size(mlab,1) < 1 .or. size(mlab,2) /= 3) return
        if (.not. all(ieee_is_finite(mlab))) return
        n = size(mlab,1)+1
        allocate(counts(2*n-1))
        counts = 0
        counts(:n) = 1
        do k = 1,n-1
            if (any(mlab(k,1:2) < 1.0_dp)) return
            if (any(mlab(k,1:2) > real(n+k-1,dp))) return
            if (any(mlab(k,1:2) /= real(int(mlab(k,1:2)),dp))) return
            if (mlab(k,3) < 0.0_dp) return
            a = int(mlab(k,1))
            b = int(mlab(k,2))
            if (a == b) return
            if (counts(a) == 0 .or. counts(b) == 0) return
            counts(n+k) = counts(a)+counts(b)
            counts(a) = 0
            counts(b) = 0
        end do
        ! Rebuild node sizes using the child references, which still
        ! remain available in the original input matrix.
        counts = 0
        counts(:n) = 1
        deallocate(z)
        allocate(z(n-1,4))
        do k = 1,n-1
            a = int(mlab(k,1))
            b = int(mlab(k,2))
            counts(n+k) = counts(a)+counts(b)
            z(k,:) = [mlab(k,:),real(counts(n+k),dp)]
        end do
        if (.not. valid_linkage(z)) then
            deallocate(z)
            allocate(z(0,4))
            return
        end if
        if (present(status)) status = cluster_success
    end function from_mlab_linkage

    function dendrogram(z, status, count_sort, distance_sort) result(layout)
        ! Compute SciPy's *untruncated*, vertical dendrogram(no_plot=True)
        ! geometry without invoking a plotting package. Columns in each
        ! row describe the four vertices of one U-shaped link segment.
        ! Child sorting affects rendering but never changes the input tree.
        real(dp), intent(in) :: z(:,:) !! Valid one-based linkage matrix.
        integer, intent(out), optional :: status !! Zero on success.
        character(len=*), intent(in), optional :: count_sort !! 'none', 'ascending', or 'descending'.
        character(len=*), intent(in), optional :: distance_sort !! 'none', 'ascending', or 'descending'.
        type(cluster_dendrogram) :: layout
        integer, allocatable :: left(:), right(:), stack(:)
        real(dp), allocatable :: center(:), height(:)
        character(len=:), allocatable :: count_rule, distance_rule
        integer :: n, k, a, b, node, child, top, nleaves, nlinks

        if (present(status)) status = cluster_invalid_input
        allocate(layout%leaves(0), layout%icoord(0,4), layout%dcoord(0,4))
        if (.not. valid_linkage(z)) return
        count_rule = 'none'
        distance_rule = 'none'
        if (present(count_sort)) count_rule = trim(count_sort)
        if (present(distance_sort)) distance_rule = trim(distance_sort)
        if (.not. valid_sort_rule(count_rule) .or. &
            .not. valid_sort_rule(distance_rule)) return
        if (count_rule /= 'none' .and. distance_rule /= 'none') return

        n = size(z,1)+1
        allocate(left(2*n-1),right(2*n-1),center(2*n-1), &
            height(2*n-1),stack(4*n))
        left = 0
        right = 0
        height = 0.0_dp
        center = 0.0_dp
        do k = 1,n-1
            node = n+k
            a = int(z(k,1))
            b = int(z(k,2))
            height(node) = z(k,3)
            if (count_rule /= 'none') then
                if (count_rule == 'ascending') then
                    if (subtree_count(a,z,n) > subtree_count(b,z,n)) then
                        child = a
                        a = b
                        b = child
                    end if
                else
                    if (subtree_count(a,z,n) <= subtree_count(b,z,n)) then
                        child = a
                        a = b
                        b = child
                    end if
                end if
            else if (distance_rule /= 'none') then
                if (distance_rule == 'ascending') then
                    if (height(a) > height(b)) then
                        child = a
                        a = b
                        b = child
                    end if
                else
                    if (height(a) <= height(b)) then
                        child = a
                        a = b
                        b = child
                    end if
                end if
            end if
            left(node) = a
            right(node) = b
        end do

        deallocate(layout%leaves,layout%icoord,layout%dcoord)
        allocate(layout%leaves(n),layout%icoord(n-1,4),layout%dcoord(n-1,4))
        top = 1
        stack(top) = 2*n-1
        nleaves = 0
        nlinks = 0
        do while (top > 0)
            node = stack(top)
            top = top-1
            if (node < 0) then
                node = -node
                a = left(node)
                b = right(node)
                nlinks = nlinks+1
                layout%icoord(nlinks,:) = [center(a),center(a),center(b),center(b)]
                layout%dcoord(nlinks,:) = [height(a),height(node),height(node),height(b)]
                center(node) = 0.5_dp*(center(a)+center(b))
            else if (node <= n) then
                nleaves = nleaves+1
                layout%leaves(nleaves) = node
                center(node) = 10.0_dp*real(nleaves,dp)-5.0_dp
            else
                top = top+1
                stack(top) = -node
                top = top+1
                stack(top) = right(node)
                top = top+1
                stack(top) = left(node)
            end if
        end do
        if (present(status)) status = cluster_success
    end function dendrogram

    pure logical function valid_sort_rule(rule) result(ok)
        character(len=*), intent(in) :: rule !! Sorting convention.
        select case (rule)
        case ('none', 'ascending', 'descending')
            ok = .true.
        case default
            ok = .false.
        end select
    end function valid_sort_rule

    pure integer function subtree_count(node,z,n) result(count)
        integer, intent(in) :: node !! One-based cluster node ID.
        real(dp), intent(in) :: z(:,:) !! Valid linkage matrix.
        integer, intent(in) :: n !! Number of original observations.
        count = 1
        if (node > n) count = int(z(node-n,4))
    end function subtree_count

end module scifort_cluster_hierarchy
