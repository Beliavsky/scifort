! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Balanced k-d tree for reusable Minkowski nearest-neighbor and radius queries.
! Independent implementation.  Indices are 1-based in the original input order.
! Exact ties prefer the smaller original index, including duplicate sites.
! No BLAS, LAPACK, C bindings, or external spatial-index dependencies.
module scifort_interpolate_kdtree
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf
    use, intrinsic :: iso_fortran_env, only: int64
    use scifort_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: kdtree_success = 0
    integer, parameter, public :: kdtree_invalid_input = 1

    type, public :: kd_tree_type
        private
        real(dp), allocatable :: sites(:,:)       ! original sites (n, dimensions)
        real(dp), allocatable :: period(:)        ! box length per dimension; zero means nonperiodic
        real(dp), allocatable :: low(:,:), high(:,:) ! node bounding boxes (dimensions,n)
        integer, allocatable :: site(:), lower(:), upper(:) ! per-node site and two children
        integer, allocatable :: split_axis(:), subtree_size(:)
        integer :: root = 0
    contains
        procedure :: valid => tree_valid
        procedure :: dimensions => tree_dimensions
        procedure :: size => tree_size
        procedure :: nearest => tree_nearest
        procedure :: query => tree_query
        procedure :: query_ball_point => tree_radius
        procedure :: query_pairs => tree_query_pairs
        procedure :: query_ball_tree => tree_query_ball_tree
        procedure, private :: count_neighbors_scalar => tree_count_neighbors
        procedure, private :: count_neighbors_radii => tree_count_neighbors_radii
        generic :: count_neighbors => count_neighbors_scalar, count_neighbors_radii
        procedure, private :: weighted_count_neighbors_scalar => tree_weighted_count_neighbors
        procedure, private :: weighted_count_neighbors_radii => tree_weighted_count_neighbors_radii
        generic :: weighted_count_neighbors => weighted_count_neighbors_scalar, weighted_count_neighbors_radii
        procedure :: sparse_distance_matrix => tree_sparse_distance_matrix
    end type kd_tree_type

    public :: kd_tree

    interface kd_tree
        module procedure kd_tree_vector
        module procedure kd_tree_scalar
    end interface kd_tree

contains

    function kd_tree_vector(points,status,boxsize) result(tree)
        real(dp), intent(in) :: points(:,:) !! sites: (observations, coordinates)
        integer, intent(out), optional :: status !! 0 valid, 1 invalid dimensions or nonfinite coordinates
        real(dp), intent(in), optional :: boxsize(:) !! periodic box lengths per dimension; zero disables wrapping
        type(kd_tree_type) :: tree
        integer, allocatable :: order(:)
        integer :: n, d, i, code

        code = kdtree_invalid_input
        n = size(points,1)
        d = size(points,2)
        if (n>0 .and. d>0) then
            if (present(boxsize)) then
                if (size(boxsize)/=d) goto 900
                if (.not.all(ieee_is_finite(boxsize))) goto 900
                if (any(boxsize<0.0_dp)) goto 900
                do i=1,d
                    if (boxsize(i)==0.0_dp) cycle
                    if (any(points(:,i)<0.0_dp) .or. any(points(:,i)>=boxsize(i))) goto 900
                end do
            end if
            if (all(ieee_is_finite(points))) then
                tree%sites = points
                allocate(tree%period(d))
                tree%period=0.0_dp
                if (present(boxsize)) tree%period=boxsize
                allocate(tree%low(d,n),tree%high(d,n),tree%site(n),tree%lower(n), &
                    tree%upper(n),tree%split_axis(n),tree%subtree_size(n),order(n))
                order = [(i,i=1,n)]
                tree%lower=0
                tree%upper=0
                tree%root=build_node(tree,order,1,n,1)
                code=kdtree_success
            end if
        end if
900     if (present(status)) status=code
    end function kd_tree_vector

    function kd_tree_scalar(points,boxsize,status) result(tree)
        real(dp), intent(in) :: points(:,:) !! sites: (observations, coordinates)
        real(dp), intent(in) :: boxsize !! same periodic box length on each axis; 0 disables wrapping
        integer, intent(out), optional :: status !! 0 valid, 1 invalid dimensions or box size
        type(kd_tree_type) :: tree
        real(dp) :: periods(size(points,2))
        periods=boxsize
        tree=kd_tree_vector(points,status=status,boxsize=periods)
    end function kd_tree_scalar

    pure logical function tree_valid(self) result(ok)
        class(kd_tree_type), intent(in) :: self !! index to validate
        ok = allocated(self%sites) .and. allocated(self%site)
        if (ok) ok=self%root>0
    end function tree_valid

    pure integer function tree_dimensions(self) result(d)
        class(kd_tree_type), intent(in) :: self !! index to examine
        d=0
        if (allocated(self%sites)) d=size(self%sites,2)
    end function tree_dimensions

    pure integer function tree_size(self) result(n)
        class(kd_tree_type), intent(in) :: self !! index to examine
        n=0
        if (allocated(self%sites)) n=size(self%sites,1)
    end function tree_size

    pure logical function valid_p(p) result(ok)
        real(dp), intent(in) :: p !! Minkowski order; 1 <= p <= infinity
        ok=p>=1.0_dp .and. (ieee_is_finite(p) .or. p==ieee_value(0.0_dp,ieee_positive_inf))
    end function valid_p

    ! Evaluate a p-norm without forming unscaled powers of large or tiny
    ! coordinate gaps. Infinity gives the Chebyshev (maximum) norm.
    pure real(dp) function scaled_lp(gaps,p) result(distance)
        real(dp), intent(in) :: gaps(:) !! absolute coordinate gaps
        real(dp), intent(in) :: p !! valid Minkowski order
        real(dp) :: scale
        scale=maxval(gaps)
        distance=scale
        if (scale==0.0_dp .or. .not.ieee_is_finite(scale)) return
        if (p==ieee_value(0.0_dp,ieee_positive_inf)) return
        if (p==1.0_dp) then
            distance=scale*sum(gaps/scale)
        else if (p==2.0_dp) then
            distance=scale*sqrt(sum((gaps/scale)**2))
        else
            distance=scale*sum((gaps/scale)**p)**(1.0_dp/p)
        end if
    end function scaled_lp

    pure real(dp) function point_lp(self,a,b,p) result(distance)
        class(kd_tree_type), intent(in) :: self !! supplying tree and periodic topology
        real(dp), intent(in) :: a(:) !! first point
        real(dp), intent(in) :: b(:) !! second point
        real(dp), intent(in) :: p !! Minkowski order
        real(dp) :: gaps(size(a)),delta,wrapped
        integer :: i
        do i=1,size(a)
            if (self%period(i)>0.0_dp) then
                ! a may be an arbitrary query point, and is wrapped here.
                delta=abs(modulo(a(i),self%period(i))-modulo(b(i),self%period(i)))
                wrapped=self%period(i)-delta
                gaps(i)=min(delta,wrapped)
            else
                gaps(i)=abs(a(i)-b(i))
            end if
        end do
        distance=scaled_lp(gaps,p)
    end function point_lp

    ! Minimum/maximum possible separations from a point to an axis-aligned
    ! node box.  Periodic maxima use a conservative L/2 bound.
    pure real(dp) function point_box_lp(self,node,x,p,upper_bound) result(distance)
        class(kd_tree_type), intent(in) :: self !! tree of node boxes
        integer, intent(in) :: node !! subtree identifier
        real(dp), intent(in) :: x(:) !! query coordinates
        real(dp), intent(in) :: p !! Minkowski order
        logical, intent(in) :: upper_bound !! true for a safe maximum
        real(dp) :: gaps(size(x)),q,lo,hi,l
        integer :: i
        do i=1,size(x)
            lo=self%low(i,node)
            hi=self%high(i,node)
            l=self%period(i)
            q=x(i)
            if (l>0.0_dp) q=modulo(q,l)
            if (upper_bound) then
                gaps(i)=max(abs(lo-q),abs(hi-q))
                if (l>0.0_dp) gaps(i)=min(gaps(i),0.5_dp*l)
            else if (q<lo) then
                gaps(i)=lo-q
                if (l>0.0_dp) gaps(i)=min(gaps(i),(l-hi)+q)
            else if (q>hi) then
                gaps(i)=q-hi
                if (l>0.0_dp) gaps(i)=min(gaps(i),(l-q)+lo)
            else
                gaps(i)=0.0_dp
            end if
        end do
        distance=scaled_lp(gaps,p)
    end function point_box_lp

    pure real(dp) function box_box_lp(self,node,other,other_node,p,upper_bound) result(distance)
        class(kd_tree_type), intent(in) :: self !! first tree
        integer, intent(in) :: node !! first subtree
        class(kd_tree_type), intent(in) :: other !! second tree
        integer, intent(in) :: other_node !! second subtree
        real(dp), intent(in) :: p !! Minkowski order
        logical, intent(in) :: upper_bound !! maximum if true; minimum otherwise
        real(dp) :: gaps(size(self%period)),al,ah,bl,bh,l
        integer :: i
        do i=1,size(gaps)
            al=self%low(i,node)
            ah=self%high(i,node)
            bl=other%low(i,other_node)
            bh=other%high(i,other_node)
            l=self%period(i)
            if (upper_bound) then
                gaps(i)=max(abs(al-bh),abs(ah-bl))
                if (l>0.0_dp) gaps(i)=min(gaps(i),0.5_dp*l)
            else if (al>bh) then
                gaps(i)=al-bh
                if (l>0.0_dp) gaps(i)=min(gaps(i),(l-ah)+bl)
            else if (bl>ah) then
                gaps(i)=bl-ah
                if (l>0.0_dp) gaps(i)=min(gaps(i),(l-bh)+al)
            else
                gaps(i)=0.0_dp
            end if
        end do
        distance=scaled_lp(gaps,p)
    end function box_box_lp

    pure logical function matching_topology(self,other) result(ok)
        class(kd_tree_type), intent(in) :: self !! source tree
        class(kd_tree_type), intent(in) :: other !! destination tree
        ok=.false.
        if (.not.self%valid() .or. .not.other%valid()) return
        if (self%dimensions()/=other%dimensions()) return
        ok=all(self%period==other%period)
    end function matching_topology

    pure logical function less_site(self,i,j,axis) result(before)
        type(kd_tree_type), intent(in) :: self !! index containing coordinates
        integer, intent(in) :: i !! first input site
        integer, intent(in) :: j !! second input site
        integer, intent(in) :: axis !! active coordinate index
        before=self%sites(i,axis)<self%sites(j,axis)
        ! Exact ties are resolved by input position; the sites are finite.
        if (.not.before .and. .not.(self%sites(i,axis)>self%sites(j,axis))) before=i<j
    end function less_site

    subroutine partition_median(self,order,lo,hi,median,axis)
        type(kd_tree_type), intent(in) :: self !! immutable site coordinates
        integer, intent(inout) :: order(:) !! permuted original input indices
        integer, intent(in) :: lo !! start of active segment
        integer, intent(in) :: hi !! end of active segment
        integer, intent(in) :: median !! requested order statistic
        integer, intent(in) :: axis !! active coordinate
        integer :: left,right,pivot,i,j,t
        left=lo
        right=hi
        do while (left<right)
            pivot=order((left+right)/2)
            i=left
            j=right
            do
                do while (less_site(self,order(i),pivot,axis))
                    i=i+1
                end do
                do while (less_site(self,pivot,order(j),axis))
                    j=j-1
                end do
                if (i>j) exit
                t=order(i)
                order(i)=order(j)
                order(j)=t
                i=i+1
                j=j-1
                if (i>j) exit
            end do
            if (median<=j) then
                right=j
            else if (median>=i) then
                left=i
            else
                exit
            end if
        end do
    end subroutine partition_median

    recursive integer function build_node(self,order,lo,hi,depth) result(node)
        type(kd_tree_type), intent(inout) :: self !! index being built
        integer, intent(inout) :: order(:) !! in-place partitioned input indices
        integer, intent(in) :: lo !! first input index in partition
        integer, intent(in) :: hi !! last input index in partition
        integer, intent(in) :: depth !! one-based tree depth
        integer :: axis,child,d
        node=0
        if (lo>hi) return
        d=size(self%sites,2)
        axis=mod(depth-1,d)+1
        node=lo+(hi-lo)/2
        call partition_median(self,order,lo,hi,node,axis)
        self%site(node)=order(node)
        self%subtree_size(node)=hi-lo+1
        self%split_axis(node)=axis
        self%low(:,node)=self%sites(self%site(node),:)
        self%high(:,node)=self%low(:,node)
        if (lo<node) then
            child=build_node(self,order,lo,node-1,depth+1)
            self%lower(node)=child
            self%low(:,node)=min(self%low(:,node),self%low(:,child))
            self%high(:,node)=max(self%high(:,node),self%high(:,child))
        end if
        if (node<hi) then
            child=build_node(self,order,node+1,hi,depth+1)
            self%upper(node)=child
            self%low(:,node)=min(self%low(:,node),self%low(:,child))
            self%high(:,node)=max(self%high(:,node),self%high(:,child))
        end if
    end function build_node

    pure logical function worse(da,ia,db,ib) result(bad)
        real(dp), intent(in) :: da !! first squared distance
        integer, intent(in) :: ia !! first original site index
        real(dp), intent(in) :: db !! second squared distance
        integer, intent(in) :: ib !! second original site index
        bad=da>db .or. (.not.(da>db) .and. .not.(da<db) .and. ia>ib)
    end function worse

    subroutine heap_down(ids,dist,start,last)
        integer, intent(inout) :: ids(:) !! max-heap of original site indices
        real(dp), intent(inout) :: dist(:) !! corresponding squared distances
        integer, intent(in) :: start !! node to sift
        integer, intent(in) :: last !! last valid heap node
        integer :: parent,child,tmpid
        real(dp) :: tmpdist
        parent=start
        do while (parent<=last/2)
            child=2*parent
            if (child<last) then
                if (worse(dist(child+1),ids(child+1),dist(child),ids(child))) child=child+1
            end if
            if (.not.worse(dist(child),ids(child),dist(parent),ids(parent))) exit
            tmpid=ids(parent)
            ids(parent)=ids(child)
            ids(child)=tmpid
            tmpdist=dist(parent)
            dist(parent)=dist(child)
            dist(child)=tmpdist
            parent=child
        end do
    end subroutine heap_down

    subroutine heap_insert(ids,dist,count,k,index,d2)
        integer, intent(inout) :: ids(:) !! selected site indices
        real(dp), intent(inout) :: dist(:) !! selected squared distances
        integer, intent(inout) :: count !! current number of selected sites
        integer, intent(in) :: k !! requested neighbor count
        integer, intent(in) :: index !! candidate site index
        real(dp), intent(in) :: d2 !! candidate squared distance
        integer :: child,parent
        if (count<k) then
            count=count+1
            child=count
            do while (child>1)
                parent=child/2
                if (.not.worse(d2,index,dist(parent),ids(parent))) exit
                ids(child)=ids(parent)
                dist(child)=dist(parent)
                child=parent
            end do
            ids(child)=index
            dist(child)=d2
        else if (worse(dist(1),ids(1),d2,index)) then
            ids(1)=index
            dist(1)=d2
            call heap_down(ids,dist,1,k)
        end if
    end subroutine heap_insert

    pure real(dp) function min_distance2(self,node,x) result(distance)
        class(kd_tree_type), intent(in) :: self !! index with bounding boxes
        integer, intent(in) :: node !! index of subtree
        real(dp), intent(in) :: x(:) !! query coordinates
        distance=sum(max(max(self%low(:,node)-x,x-self%high(:,node)),0.0_dp)**2)
    end function min_distance2

    recursive subroutine search_node(self,node,x,k,ids,dist,count,visited)
        class(kd_tree_type), intent(in) :: self !! balanced index
        integer, intent(in) :: node !! subtree node identifier
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(in) :: k !! requested neighbor count
        integer, intent(inout) :: ids(:) !! candidate index max-heap
        real(dp), intent(inout) :: dist(:) !! candidate distance max-heap
        integer, intent(inout) :: count !! heap entries stored
        integer, intent(inout) :: visited !! exact site distances calculated
        integer :: a,b,t
        real(dp) :: da,db,d2
        if (node==0) return
        if (count==k) then
            if (min_distance2(self,node,x)>dist(1)) return
        end if
        d2=sum((self%sites(self%site(node),:)-x)**2)
        visited=visited+1
        call heap_insert(ids,dist,count,k,self%site(node),d2)
        a=self%lower(node)
        b=self%upper(node)
        if (a==0 .and. b==0) return
        if (a==0) then
            call search_node(self,b,x,k,ids,dist,count,visited)
        else if (b==0) then
            call search_node(self,a,x,k,ids,dist,count,visited)
        else
            da=min_distance2(self,a,x)
            db=min_distance2(self,b,x)
            if (db<da) then
                t=a
                a=b
                b=t
            end if
            call search_node(self,a,x,k,ids,dist,count,visited)
            call search_node(self,b,x,k,ids,dist,count,visited)
        end if
    end subroutine search_node

    recursive subroutine search_node_lp(self,node,x,p,k,ids,dist,count,visited)
        class(kd_tree_type), intent(in) :: self !! search tree with optional periodic topology
        integer, intent(in) :: node !! subtree to examine
        real(dp), intent(in) :: x(:) !! query coordinates
        real(dp), intent(in) :: p !! Minkowski order
        integer, intent(in) :: k !! number of nearest neighbors
        integer, intent(inout) :: ids(:) !! max heap of site identifiers
        real(dp), intent(inout) :: dist(:) !! max heap of Minkowski distances
        integer, intent(inout) :: count !! current heap size
        integer, intent(inout) :: visited !! exact point distances evaluated
        integer :: a,b,t
        real(dp) :: da,db,d
        if (node==0) return
        if (count==k) then
            if (point_box_lp(self,node,x,p,.false.)>dist(1)) return
        end if
        d=point_lp(self,x,self%sites(self%site(node),:),p)
        visited=visited+1
        call heap_insert(ids,dist,count,k,self%site(node),d)
        a=self%lower(node)
        b=self%upper(node)
        if (a==0 .and. b==0) return
        if (a==0) then
            call search_node_lp(self,b,x,p,k,ids,dist,count,visited)
        else if (b==0) then
            call search_node_lp(self,a,x,p,k,ids,dist,count,visited)
        else
            da=point_box_lp(self,a,x,p,.false.)
            db=point_box_lp(self,b,x,p,.false.)
            if (db<da) then
                t=a
                a=b
                b=t
            end if
            call search_node_lp(self,a,x,p,k,ids,dist,count,visited)
            call search_node_lp(self,b,x,p,k,ids,dist,count,visited)
        end if
    end subroutine search_node_lp

    subroutine tree_query(self,x,indices,distance2,k,status,visited,p,distance)
        class(kd_tree_type), intent(in) :: self !! index to query
        real(dp), intent(in) :: x(:) !! coordinates (dimension)
        integer, allocatable, intent(out) :: indices(:) !! original 1-based site indices in ascending distance order
        real(dp), allocatable, intent(out), optional :: distance2(:) !! squared Euclidean distances in same order
        integer, intent(in), optional :: k !! number of nearest sites (default 1)
        integer, intent(out), optional :: status !! 0 success, 1 invalid request
        integer, intent(out), optional :: visited !! number of exact site-distance computations
        real(dp), intent(in), optional :: p !! Minkowski order (default 2, may be +infinity)
        real(dp), allocatable, intent(out), optional :: distance(:) !! unsquared Minkowski distances
        real(dp), allocatable :: dist(:)
        real(dp) :: tmpdist
        real(dp) :: normp
        integer :: count,want,seen,last,tmpid,code
        logical :: generic
        code=kdtree_invalid_input
        seen=0
        allocate(indices(0))
        if (present(distance2)) allocate(distance2(0))
        if (present(distance)) allocate(distance(0))
        if (.not.self%valid()) goto 900
        if (size(x)/=self%dimensions()) goto 900
        if (.not.all(ieee_is_finite(x))) goto 900
        normp=2.0_dp
        if (present(p)) normp=p
        if (.not.valid_p(normp)) goto 900
        ! An explicit unsquared distance requests the numerically scaled
        ! path, even for a nonperiodic Euclidean tree.
        generic=any(self%period>0.0_dp) .or. present(p) .or. present(distance)
        want=1
        if (present(k)) want=k
        if (want<1 .or. want>self%size()) goto 900
        deallocate(indices)
        allocate(indices(want),dist(want))
        count=0
        if (generic) then
            call search_node_lp(self,self%root,x,normp,want,indices,dist,count,seen)
        else
            call search_node(self,self%root,x,want,indices,dist,count,seen)
        end if
        if (count/=want) error stop 'kd_tree internal error: incomplete nearest-neighbor heap'
        ! Max-heap sort yields nearest-first order; ties favor lowest input index.
        do last=want,2,-1
            tmpid=indices(1)
            indices(1)=indices(last)
            indices(last)=tmpid
            tmpdist=dist(1)
            dist(1)=dist(last)
            dist(last)=tmpdist
            call heap_down(indices,dist,1,last-1)
        end do
        if (present(distance2)) then
            if (generic) then
                distance2=dist**2
            else
                distance2=dist
            end if
        end if
        if (present(distance)) then
            if (generic) then
                distance=dist
            else
                distance=sqrt(dist)
            end if
        end if
        code=kdtree_success
900     if (present(status)) status=code
        if (present(visited)) visited=seen
    end subroutine tree_query

    integer function tree_nearest(self,x,status,visited,p) result(id)
        class(kd_tree_type), intent(in) :: self !! index to query
        real(dp), intent(in) :: x(:) !! position
        integer, intent(out), optional :: status !! 0 success, 1 invalid
        integer, intent(out), optional :: visited !! exact point evaluations
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        integer, allocatable :: ids(:)
        integer :: code,seen
        if (present(p)) then
            call self%query(x,ids,k=1,status=code,visited=seen,p=p)
        else
            call self%query(x,ids,k=1,status=code,visited=seen)
        end if
        id=0
        if (code==kdtree_success) id=ids(1)
        if (present(status)) status=code
        if (present(visited)) visited=seen
    end function tree_nearest

    subroutine tree_radius(self,x,radius,indices,status,visited,p)
        class(kd_tree_type), intent(in) :: self !! index to query
        real(dp), intent(in) :: x(:) !! position
        real(dp), intent(in) :: radius !! finite nonnegative search radius
        integer, allocatable, intent(out) :: indices(:) !! indices of sites within closed ball, ascending input order
        integer, intent(out), optional :: status !! 0 success, 1 invalid
        integer, intent(out), optional :: visited !! actual site distances evaluated
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        integer, allocatable :: buffer(:)
        integer :: count,seen,code
        real(dp) :: normp
        allocate(indices(0))
        seen=0
        code=kdtree_invalid_input
        if (.not.self%valid()) goto 900
        if (size(x)/=self%dimensions()) goto 900
        if (.not.all(ieee_is_finite(x)) .or. .not.ieee_is_finite(radius)) goto 900
        if (radius<0.0_dp) goto 900
        normp=2.0_dp
        if (present(p)) normp=p
        if (.not.valid_p(normp)) goto 900
        allocate(buffer(self%size()))
        count=0
        call radius_node(self,self%root,x,radius,normp,buffer,count,seen)
        call sort_indices(buffer,count)
        indices=buffer(:count)
        code=kdtree_success
900     if (present(status)) status=code
        if (present(visited)) visited=seen
    end subroutine tree_radius

    recursive subroutine radius_node(self,node,x,radius,p,buffer,count,visited,min_index)
        class(kd_tree_type), intent(in) :: self !! index to search
        integer, intent(in) :: node !! subtree root
        real(dp), intent(in) :: x(:) !! query point
        real(dp), intent(in) :: radius !! nonnegative inclusive Minkowski search radius
        real(dp), intent(in) :: p !! Minkowski order
        integer, intent(inout) :: buffer(:) !! matched original site indices
        integer, intent(inout) :: count !! number of matches
        integer, intent(inout) :: visited !! exact site distances evaluated
        integer, intent(in), optional :: min_index !! retain original site indices strictly greater than this
        logical :: eligible
        if (node==0) return
        if (point_box_lp(self,node,x,p,.false.)>radius) return
        eligible=.true.
        if (present(min_index)) eligible=self%site(node)>min_index
        if (eligible) then
            visited=visited+1
            if (point_lp(self,x,self%sites(self%site(node),:),p)<=radius) then
                count=count+1
                buffer(count)=self%site(node)
            end if
        end if
        call radius_node(self,self%lower(node),x,radius,p,buffer,count,visited,min_index)
        call radius_node(self,self%upper(node),x,radius,p,buffer,count,visited,min_index)
    end subroutine radius_node

    recursive subroutine count_point_node(self,node,x,radius,p,count,visited)
        class(kd_tree_type), intent(in) :: self !! tree to count points in
        integer, intent(in) :: node !! subtree root
        real(dp), intent(in) :: x(:) !! external query point
        real(dp), intent(in) :: radius !! closed-ball radius
        real(dp), intent(in) :: p !! Minkowski order
        integer(int64), intent(inout) :: count !! accumulated qualifying point pairs
        integer, intent(inout) :: visited !! actual point distance evaluations
        if (node==0) return
        if (point_box_lp(self,node,x,p,.false.)>radius) return
        if (point_box_lp(self,node,x,p,.true.)<=radius) then
            count=count+int(self%subtree_size(node),int64)
            return
        end if
        visited=visited+1
        if (point_lp(self,x,self%sites(self%site(node),:),p)<=radius) count=count+1_int64
        call count_point_node(self,self%lower(node),x,radius,p,count,visited)
        call count_point_node(self,self%upper(node),x,radius,p,count,visited)
    end subroutine count_point_node

    recursive subroutine count_tree_nodes(self,node,other,other_node,radius,p,count,visited)
        class(kd_tree_type), intent(in) :: self !! first tree to count pairs from
        integer, intent(in) :: node !! first subtree root
        class(kd_tree_type), intent(in) :: other !! second tree to count pairs from
        integer, intent(in) :: other_node !! second subtree root
        real(dp), intent(in) :: radius !! closed-ball radius
        real(dp), intent(in) :: p !! Minkowski order
        integer(int64), intent(inout) :: count !! accumulated ordered pair count
        integer, intent(inout) :: visited !! evaluated site-pair distances
        if (node==0 .or. other_node==0) return
        if (box_box_lp(self,node,other,other_node,p,.false.)>radius) return
        if (box_box_lp(self,node,other,other_node,p,.true.)<=radius) then
            count=count+int(self%subtree_size(node),int64)*int(other%subtree_size(other_node),int64)
            return
        end if
        ! Partition whichever side contains more sites; each pair belongs to
        ! exactly one singleton-vs-subtree or child-vs-subtree recursion.
        if (self%subtree_size(node)>=other%subtree_size(other_node)) then
            call count_point_node(other,other_node,self%sites(self%site(node),:),radius,p,count,visited)
            call count_tree_nodes(self,self%lower(node),other,other_node,radius,p,count,visited)
            call count_tree_nodes(self,self%upper(node),other,other_node,radius,p,count,visited)
        else
            call count_point_node(self,node,other%sites(other%site(other_node),:),radius,p,count,visited)
            call count_tree_nodes(self,node,other,other%lower(other_node),radius,p,count,visited)
            call count_tree_nodes(self,node,other,other%upper(other_node),radius,p,count,visited)
        end if
    end subroutine count_tree_nodes

    ! Dynamic results grow geometrically, avoiding O(n**2) storage for sparse queries.
    subroutine reserve_indices(buffer,needed,used)
        integer, allocatable, intent(inout) :: buffer(:) !! storage capacity to grow
        integer, intent(in) :: needed !! required number of entries
        integer, intent(in) :: used !! initialized entries to preserve
        integer, allocatable :: work(:)
        integer :: capacity
        if (needed<=size(buffer)) return
        capacity=max(64,size(buffer))
        do while (capacity<needed)
            if (capacity>huge(capacity)/2) then
                capacity=needed
                exit
            end if
            capacity=2*capacity
        end do
        allocate(work(capacity))
        if (used>0) work(:used)=buffer(:used)
        call move_alloc(work,buffer)
    end subroutine reserve_indices

    subroutine reserve_pairs(buffer,needed,used)
        integer, allocatable, intent(inout) :: buffer(:,:) !! pair rows, capacity by two
        integer, intent(in) :: needed !! required pair count
        integer, intent(in) :: used !! initialized pair rows to preserve
        integer, allocatable :: work(:,:)
        integer :: capacity
        if (needed<=size(buffer,1)) return
        capacity=max(64,size(buffer,1))
        do while (capacity<needed)
            if (capacity>huge(capacity)/2) then
                capacity=needed
                exit
            end if
            capacity=2*capacity
        end do
        allocate(work(capacity,2))
        if (used>0) work(:used,:)=buffer(:used,:)
        call move_alloc(work,buffer)
    end subroutine reserve_pairs

    subroutine tree_query_pairs(self,radius,pairs,status,visited,p)
        class(kd_tree_type), intent(in) :: self !! tree to search for unordered site pairs
        real(dp), intent(in) :: radius !! nonnegative inclusive Euclidean distance
        integer, allocatable, intent(out) :: pairs(:,:) !! (npairs,2), 1-based i<j, lexicographic order
        integer, intent(out), optional :: status !! zero on success; one on invalid request or overflow
        integer, intent(out), optional :: visited !! number of exact site-pair distances evaluated
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        integer, allocatable :: sites_in_ball(:), work(:,:)
        integer :: code,seen,i,matched,used
        real(dp) :: normp
        code=kdtree_invalid_input
        seen=0
        allocate(pairs(0,2))
        if (.not.self%valid()) goto 900
        if (.not.ieee_is_finite(radius)) goto 900
        if (radius<0.0_dp) goto 900
        normp=2.0_dp
        if (present(p)) normp=p
        if (.not.valid_p(normp)) goto 900
        allocate(sites_in_ball(self%size()),work(0,2))
        used=0
        do i=1,self%size()
            matched=0
            call radius_node(self,self%root,self%sites(i,:),radius,normp,sites_in_ball,matched,seen,min_index=i)
            if (matched==0) cycle
            if (matched>huge(used)-used) goto 900
            call sort_indices(sites_in_ball,matched)
            call reserve_pairs(work,used+matched,used)
            work(used+1:used+matched,1)=i
            work(used+1:used+matched,2)=sites_in_ball(:matched)
            used=used+matched
        end do
        pairs=work(:used,:)
        code=kdtree_success
900     if (present(status)) status=code
        if (present(visited)) visited=seen
    end subroutine tree_query_pairs

    subroutine tree_query_ball_tree(self,other,radius,row_ptr,indices,status,visited,p)
        class(kd_tree_type), intent(in) :: self !! source tree, one output row per source site
        class(kd_tree_type), intent(in) :: other !! destination tree, column site indices
        real(dp), intent(in) :: radius !! nonnegative inclusive Euclidean distance
        integer, allocatable, intent(out) :: row_ptr(:) !! one-based CSR starts; length nsource+1
        integer, allocatable, intent(out) :: indices(:) !! 1-based destination indices, sorted within rows
        integer, intent(out), optional :: status !! zero on success; one on invalid request or overflow
        integer, intent(out), optional :: visited !! actual exact site-to-site distance evaluations
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        integer, allocatable :: sites_in_ball(:), work(:)
        integer :: code,seen,i,matched,used
        real(dp) :: normp
        code=kdtree_invalid_input
        seen=0
        allocate(row_ptr(0),indices(0))
        if (.not.matching_topology(self,other)) goto 900
        if (.not.ieee_is_finite(radius)) goto 900
        if (radius<0.0_dp) goto 900
        normp=2.0_dp
        if (present(p)) normp=p
        if (.not.valid_p(normp)) goto 900
        deallocate(row_ptr)
        allocate(row_ptr(self%size()+1),sites_in_ball(other%size()),work(0))
        row_ptr(1)=1
        used=0
        do i=1,self%size()
            matched=0
            call radius_node(other,other%root,self%sites(i,:),radius,normp,sites_in_ball,matched,seen)
            if (matched>huge(used)-used-1) goto 900
            if (matched>0) then
                call sort_indices(sites_in_ball,matched)
                call reserve_indices(work,used+matched,used)
                work(used+1:used+matched)=sites_in_ball(:matched)
                used=used+matched
            end if
            row_ptr(i+1)=used+1
        end do
        indices=work(:used)
        code=kdtree_success
900     if (code/=kdtree_success) then
            if (allocated(row_ptr)) deallocate(row_ptr)
            allocate(row_ptr(0))
        end if
        if (present(status)) status=code
        if (present(visited)) visited=seen
    end subroutine tree_query_ball_tree

    integer(int64) function tree_count_neighbors(self,other,radius,status,visited,p) result(count)
        class(kd_tree_type), intent(in) :: self !! first tree (ordered pair source)
        class(kd_tree_type), intent(in) :: other !! second tree (ordered pair destination)
        real(dp), intent(in) :: radius !! nonnegative inclusive Euclidean distance
        integer, intent(out), optional :: status !! zero on success, one on invalid request
        integer, intent(out), optional :: visited !! actual evaluated site-to-site distances
        integer :: code,seen
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        real(dp) :: normp
        count=0_int64
        code=kdtree_invalid_input
        seen=0
        if (.not.matching_topology(self,other)) goto 900
        if (.not.ieee_is_finite(radius)) goto 900
        if (radius<0.0_dp) goto 900
        normp=2.0_dp
        if (present(p)) normp=p
        if (.not.valid_p(normp)) goto 900
        call count_tree_nodes(self,self%root,other,other%root,radius,normp,count,seen)
        code=kdtree_success
900     if (present(status)) status=code
        if (present(visited)) visited=seen
    end function tree_count_neighbors

    ! Vectorized counting uses a difference array across sorted radii. Each
    ! fully enclosed subtree pair updates a whole interval in O(1), while
    ! partial pairs are split only for the smaller, unresolved radii. This
    ! avoids a separate tree traversal for every histogram bin.
    pure integer function radius_lower_bound(radii,value) result(position)
        real(dp), intent(in) :: radii(:) !! sorted nondecreasing radii
        real(dp), intent(in) :: value !! distance to find; inclusive threshold
        integer :: lo,hi,middle
        lo=1
        hi=size(radii)+1
        do while (lo<hi)
            middle=lo+(hi-lo)/2
            if (radii(middle)<value) then
                lo=middle+1
            else
                hi=middle
            end if
        end do
        position=lo
    end function radius_lower_bound

    subroutine radius_sort_order(radii,order)
        real(dp), intent(in) :: radii(:) !! finite input radii, possibly unsorted
        integer, intent(out) :: order(:) !! permutation: sorted(i)=radii(order(i))
        integer, allocatable :: work(:)
        integer :: n,width,first,middle,last,i,j,k
        n=size(radii)
        order=[(i,i=1,n)]
        allocate(work(n))
        width=1
        do while (width<n)
            first=1
            do while (first<=n)
                middle=min(first+width-1,n)
                last=min(first+2*width-1,n)
                i=first
                j=middle+1
                do k=first,last
                    if (i>middle) then
                        work(k)=order(j)
                        j=j+1
                    else if (j>last) then
                        work(k)=order(i)
                        i=i+1
                    else if (radii(order(i))<=radii(order(j))) then
                        work(k)=order(i)
                        i=i+1
                    else
                        work(k)=order(j)
                        j=j+1
                    end if
                end do
                first=last+1
            end do
            order=work
            if (width>n/2) exit
            width=2*width
        end do
    end subroutine radius_sort_order

    recursive real(dp) function subtree_weight(self,node,weights,sums) result(result_sum)
        class(kd_tree_type), intent(in) :: self !! tree defining subtree topology
        integer, intent(in) :: node !! subtree root, or zero
        real(dp), intent(in) :: weights(:) !! original site weights
        real(dp), intent(inout) :: sums(:) !! output subtree weight sums
        result_sum=0.0_dp
        if (node==0) return
        result_sum=weights(self%site(node)) &
            +subtree_weight(self,self%lower(node),weights,sums) &
            +subtree_weight(self,self%upper(node),weights,sums)
        sums(node)=result_sum
    end function subtree_weight

    subroutine add_count_range(first,last,count1,count2,weight1,weight2,weighted,delta_int,delta_real)
        integer, intent(in) :: first !! first qualifying radius index
        integer, intent(in) :: last !! final radius index accounted for by this subtree
        integer, intent(in) :: count1 !! number of source points
        integer, intent(in) :: count2 !! number of destination points
        real(dp), intent(in) :: weight1 !! source subtree weight sum
        real(dp), intent(in) :: weight2 !! destination subtree weight sum
        logical, intent(in) :: weighted !! select weighted vs integer counting
        integer(int64), intent(inout) :: delta_int(:) !! integer difference array
        real(dp), intent(inout) :: delta_real(:) !! weighted difference array
        real(dp) :: product
        integer(int64) :: pair_count
        if (first>last) return
        if (weighted) then
            product=weight1*weight2
            delta_real(first)=delta_real(first)+product
            delta_real(last+1)=delta_real(last+1)-product
        else
            pair_count=int(count1,int64)*int(count2,int64)
            delta_int(first)=delta_int(first)+pair_count
            delta_int(last+1)=delta_int(last+1)-pair_count
        end if
    end subroutine add_count_range

    recursive subroutine count_radii_point_node(point_tree,point_id,other,node,radii,limit,p, &
        point_weights,other_weights,other_sums,weighted,delta_int,delta_real,visited)
        class(kd_tree_type), intent(in) :: point_tree !! tree containing the fixed point
        integer, intent(in) :: point_id !! original 1-based fixed site index
        class(kd_tree_type), intent(in) :: other !! tree being searched
        integer, intent(in) :: node !! subtree being searched
        real(dp), intent(in) :: radii(:) !! sorted radii
        integer, intent(in) :: limit !! last unresolved radius index
        real(dp), intent(in) :: p !! Minkowski norm order
        real(dp), intent(in) :: point_weights(:) !! site weights on fixed-point side
        real(dp), intent(in) :: other_weights(:) !! site weights on searched side
        real(dp), intent(in) :: other_sums(:) !! subtree weights on searched side
        logical, intent(in) :: weighted !! weighted result when true
        integer(int64), intent(inout) :: delta_int(:) !! cumulative differences
        real(dp), intent(inout) :: delta_real(:) !! weighted cumulative differences
        integer, intent(inout) :: visited !! exact distance evaluations
        real(dp) :: min_dist,max_dist,dist
        integer :: complete,remaining,k
        if (node==0 .or. limit==0) return
        min_dist=point_box_lp(other,node,point_tree%sites(point_id,:),p,.false.)
        if (min_dist>radii(limit)) return
        max_dist=point_box_lp(other,node,point_tree%sites(point_id,:),p,.true.)
        complete=radius_lower_bound(radii,max_dist)
        if (complete<=limit) then
            call add_count_range(complete,limit,1,other%subtree_size(node), &
                point_weights(point_id),other_sums(node),weighted,delta_int,delta_real)
        end if
        remaining=min(limit,complete-1)
        if (remaining==0) return
        visited=visited+1
        dist=point_lp(point_tree,point_tree%sites(point_id,:),other%sites(other%site(node),:),p)
        k=radius_lower_bound(radii,dist)
        if (k<=remaining) then
            call add_count_range(k,remaining,1,1,point_weights(point_id), &
                other_weights(other%site(node)),weighted,delta_int,delta_real)
        end if
        call count_radii_point_node(point_tree,point_id,other,other%lower(node),radii,remaining,p, &
            point_weights,other_weights,other_sums,weighted,delta_int,delta_real,visited)
        call count_radii_point_node(point_tree,point_id,other,other%upper(node),radii,remaining,p, &
            point_weights,other_weights,other_sums,weighted,delta_int,delta_real,visited)
    end subroutine count_radii_point_node

    recursive subroutine count_radii_tree_nodes(self,node,other,other_node,radii,limit,p, &
        weight_a,weight_b,sum_a,sum_b,weighted,delta_int,delta_real,visited)
        class(kd_tree_type), intent(in) :: self !! source tree
        integer, intent(in) :: node !! source subtree root
        class(kd_tree_type), intent(in) :: other !! destination tree
        integer, intent(in) :: other_node !! destination subtree root
        real(dp), intent(in) :: radii(:) !! sorted radius thresholds
        integer, intent(in) :: limit !! last unresolved radius index
        real(dp), intent(in) :: p !! Minkowski order
        real(dp), intent(in) :: weight_a(:) !! original source weights
        real(dp), intent(in) :: weight_b(:) !! original destination weights
        real(dp), intent(in) :: sum_a(:) !! source node sums
        real(dp), intent(in) :: sum_b(:) !! destination node sums
        logical, intent(in) :: weighted !! weighted result when true
        integer(int64), intent(inout) :: delta_int(:) !! integer cumulative differences
        real(dp), intent(inout) :: delta_real(:) !! weighted cumulative differences
        integer, intent(inout) :: visited !! exact point distances
        real(dp) :: min_dist,max_dist
        integer :: complete,remaining
        if (node==0 .or. other_node==0 .or. limit==0) return
        min_dist=box_box_lp(self,node,other,other_node,p,.false.)
        if (min_dist>radii(limit)) return
        max_dist=box_box_lp(self,node,other,other_node,p,.true.)
        complete=radius_lower_bound(radii,max_dist)
        if (complete<=limit) then
            call add_count_range(complete,limit,self%subtree_size(node),other%subtree_size(other_node), &
                sum_a(node),sum_b(other_node),weighted,delta_int,delta_real)
        end if
        remaining=min(limit,complete-1)
        if (remaining==0) return
        if (self%subtree_size(node)>=other%subtree_size(other_node)) then
            call count_radii_point_node(self,self%site(node),other,other_node,radii,remaining,p, &
                weight_a,weight_b,sum_b,weighted,delta_int,delta_real,visited)
            call count_radii_tree_nodes(self,self%lower(node),other,other_node,radii,remaining,p, &
                weight_a,weight_b,sum_a,sum_b,weighted,delta_int,delta_real,visited)
            call count_radii_tree_nodes(self,self%upper(node),other,other_node,radii,remaining,p, &
                weight_a,weight_b,sum_a,sum_b,weighted,delta_int,delta_real,visited)
        else
            call count_radii_point_node(other,other%site(other_node),self,node,radii,remaining,p, &
                weight_b,weight_a,sum_a,weighted,delta_int,delta_real,visited)
            call count_radii_tree_nodes(self,node,other,other%lower(other_node),radii,remaining,p, &
                weight_a,weight_b,sum_a,sum_b,weighted,delta_int,delta_real,visited)
            call count_radii_tree_nodes(self,node,other,other%upper(other_node),radii,remaining,p, &
                weight_a,weight_b,sum_a,sum_b,weighted,delta_int,delta_real,visited)
        end if
    end subroutine count_radii_tree_nodes

    ! Shared engine for weighted and unweighted vector-radius counting.
    ! No uninitialized outputs on invalid input. A zero-length radius array
    ! is valid and returns zero-length results without any traversal.
    subroutine count_radii_engine(self,other,radii,cumulative,p,weighted, &
        weights_self,weights_other,integer_result,real_result,status,visited)
        class(kd_tree_type), intent(in) :: self !! source index
        class(kd_tree_type), intent(in) :: other !! destination index
        real(dp), intent(in) :: radii(:) !! radii to evaluate
        logical, intent(in) :: cumulative !! cumulative if true; disjoint bins otherwise
        real(dp), intent(in) :: p !! Minkowski exponent
        logical, intent(in) :: weighted !! weighted vs integer count
        real(dp), intent(in), optional :: weights_self(:) !! source site weights
        real(dp), intent(in), optional :: weights_other(:) !! destination site weights
        integer(int64), allocatable, intent(out) :: integer_result(:) !! unweighted result
        real(dp), allocatable, intent(out) :: real_result(:) !! weighted result
        integer, intent(out) :: status !! zero on success
        integer, intent(out) :: visited !! number of exact distances evaluated
        integer, allocatable :: order(:)
        real(dp), allocatable :: sorted(:),wa(:),wb(:),sum_a(:),sum_b(:),delta_real(:)
        integer(int64), allocatable :: delta_int(:)
        integer :: n,i
        integer(int64) :: running_int
        real(dp) :: running_real
        status=kdtree_invalid_input
        visited=0
        allocate(integer_result(0),real_result(0))
        if (.not.matching_topology(self,other)) return
        if (.not.valid_p(p)) return
        if (.not.all(ieee_is_finite(radii))) return
        if (any(radii<0.0_dp)) return
        n=size(radii)
        if (.not.cumulative) then
            do i=2,n
                if (radii(i)<radii(i-1)) return
            end do
        end if
        if (present(weights_self)) then
            if (size(weights_self)/=self%size()) return
            if (.not.all(ieee_is_finite(weights_self))) return
        end if
        if (present(weights_other)) then
            if (size(weights_other)/=other%size()) return
            if (.not.all(ieee_is_finite(weights_other))) return
        end if
        deallocate(integer_result,real_result)
        allocate(integer_result(n),real_result(n))
        integer_result=0_int64
        real_result=0.0_dp
        status=kdtree_success
        if (n==0) return
        allocate(order(n),sorted(n),delta_int(n+1),delta_real(n+1))
        if (cumulative) then
            call radius_sort_order(radii,order)
        else
            order=[(i,i=1,n)]
        end if
        sorted=radii(order)
        allocate(wa(self%size()),wb(other%size()),sum_a(self%size()),sum_b(other%size()))
        wa=1.0_dp
        wb=1.0_dp
        if (weighted) then
            if (present(weights_self)) wa=weights_self
            if (present(weights_other)) wb=weights_other
        end if
        running_real=subtree_weight(self,self%root,wa,sum_a)
        running_real=subtree_weight(other,other%root,wb,sum_b)
        delta_int=0_int64
        delta_real=0.0_dp
        call count_radii_tree_nodes(self,self%root,other,other%root,sorted,n,p,wa,wb,sum_a,sum_b, &
            weighted,delta_int,delta_real,visited)
        if (weighted) then
            running_real=0.0_dp
            do i=1,n
                if (cumulative) then
                    running_real=running_real+delta_real(i)
                    real_result(order(i))=running_real
                else
                    real_result(i)=delta_real(i)
                end if
            end do
        else
            running_int=0_int64
            do i=1,n
                if (cumulative) then
                    running_int=running_int+delta_int(i)
                    integer_result(order(i))=running_int
                else
                    integer_result(i)=delta_int(i)
                end if
            end do
        end if
    end subroutine count_radii_engine

    function tree_count_neighbors_radii(self,other,radii,status,visited,p,cumulative) result(counts)
        class(kd_tree_type), intent(in) :: self !! source tree
        class(kd_tree_type), intent(in) :: other !! destination tree
        real(dp), intent(in) :: radii(:) !! radius thresholds, any order if cumulative
        integer, intent(out), optional :: status !! zero success; one invalid input
        integer, intent(out), optional :: visited !! actual point-distance evaluations
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        logical, intent(in), optional :: cumulative !! true cumulative, false disjoint bins
        integer(int64), allocatable :: counts(:)
        integer(int64), allocatable :: result_int(:)
        real(dp), allocatable :: result_real(:)
        real(dp) :: normp
        logical :: is_cumulative
        integer :: code,seen
        normp=2.0_dp
        if (present(p)) normp=p
        is_cumulative=.true.
        if (present(cumulative)) is_cumulative=cumulative
        call count_radii_engine(self,other,radii,is_cumulative,normp,.false., &
            integer_result=result_int,real_result=result_real,status=code,visited=seen)
        call move_alloc(result_int,counts)
        if (present(status)) status=code
        if (present(visited)) visited=seen
    end function tree_count_neighbors_radii

    real(dp) function tree_weighted_count_neighbors(self,other,radius,weights_self,weights_other, &
        status,visited,p) result(weight_sum)
        class(kd_tree_type), intent(in) :: self !! first tree
        class(kd_tree_type), intent(in) :: other !! second tree
        real(dp), intent(in) :: radius !! scalar radius threshold
        real(dp), intent(in), optional :: weights_self(:) !! source site weights
        real(dp), intent(in), optional :: weights_other(:) !! destination site weights
        integer, intent(out), optional :: status !! 0 valid, 1 invalid input
        integer, intent(out), optional :: visited !! exact point-distance evaluations
        real(dp), intent(in), optional :: p !! Minkowski order
        real(dp), allocatable :: result_real(:)
        integer(int64), allocatable :: result_int(:)
        integer :: code,seen
        real(dp) :: normp
        normp=2.0_dp
        if (present(p)) normp=p
        call count_radii_engine(self,other,[radius],.true.,normp,.true., &
            weights_self=weights_self,weights_other=weights_other, &
            integer_result=result_int,real_result=result_real,status=code,visited=seen)
        weight_sum=0.0_dp
        if (code==kdtree_success) weight_sum=result_real(1)
        if (present(status)) status=code
        if (present(visited)) visited=seen
    end function tree_weighted_count_neighbors

    function tree_weighted_count_neighbors_radii(self,other,radii,weights_self,weights_other, &
        status,visited,p,cumulative) result(weight_sums)
        class(kd_tree_type), intent(in) :: self !! source tree
        class(kd_tree_type), intent(in) :: other !! destination tree
        real(dp), intent(in) :: radii(:) !! threshold radii
        real(dp), intent(in), optional :: weights_self(:) !! source site weights
        real(dp), intent(in), optional :: weights_other(:) !! destination site weights
        integer, intent(out), optional :: status !! 0 valid, 1 invalid input
        integer, intent(out), optional :: visited !! actual point-distance evaluations
        real(dp), intent(in), optional :: p !! Minkowski order
        logical, intent(in), optional :: cumulative !! false for disjoint bins
        real(dp), allocatable :: weight_sums(:)
        real(dp), allocatable :: result_real(:)
        integer(int64), allocatable :: result_int(:)
        integer :: code,seen
        real(dp) :: normp
        logical :: is_cumulative
        normp=2.0_dp
        if (present(p)) normp=p
        is_cumulative=.true.
        if (present(cumulative)) is_cumulative=cumulative
        call count_radii_engine(self,other,radii,is_cumulative,normp,.true., &
            weights_self=weights_self,weights_other=weights_other, &
            integer_result=result_int,real_result=result_real,status=code,visited=seen)
        call move_alloc(result_real,weight_sums)
        if (present(status)) status=code
        if (present(visited)) visited=seen
    end function tree_weighted_count_neighbors_radii

    subroutine tree_sparse_distance_matrix(self,other,max_distance,row,column,distance,status,visited,p)
        class(kd_tree_type), intent(in) :: self !! first tree (matrix row sites)
        class(kd_tree_type), intent(in) :: other !! second tree (matrix column sites)
        real(dp), intent(in) :: max_distance !! inclusive Euclidean distance cutoff
        integer, allocatable, intent(out) :: row(:) !! one-based row subscripts in ascending order
        integer, allocatable, intent(out) :: column(:) !! one-based column subscripts in ascending order per row
        real(dp), allocatable, intent(out) :: distance(:) !! Euclidean distances, including explicit zeros
        integer, intent(out), optional :: status !! zero on success, one on invalid request
        integer, intent(out), optional :: visited !! exact search site-pair distances evaluated
        real(dp), intent(in), optional :: p !! Minkowski order (default 2)
        integer, allocatable :: row_ptr(:), indices(:)
        integer :: code,seen,i,k
        real(dp) :: normp
        allocate(row(0),column(0),distance(0))
        normp=2.0_dp
        if (present(p)) normp=p
        call self%query_ball_tree(other,max_distance,row_ptr,indices,status=code,visited=seen,p=normp)
        if (code==kdtree_success) then
            deallocate(row,column,distance)
            allocate(row(size(indices)),column(size(indices)),distance(size(indices)))
            column=indices
            do i=1,self%size()
                row(row_ptr(i):row_ptr(i+1)-1)=i
            end do
            do k=1,size(indices)
                distance(k)=point_lp(self,self%sites(row(k),:),other%sites(column(k),:),normp)
            end do
        end if
        if (present(status)) status=code
        if (present(visited)) visited=seen
    end subroutine tree_sparse_distance_matrix

    subroutine sort_indices(ids,n)
        integer, intent(inout) :: ids(:) !! input IDs; sort first n in ascending order
        integer, intent(in) :: n !! number of relevant IDs
        integer :: i,last,t
        ! Sorting by original index makes the radius query deterministic.
        ! Use a max heap for O(n log n), including when the ball covers all sites.
        do i=n/2,1,-1
            call sift_integer(ids,i,n)
        end do
        do last=n,2,-1
            t=ids(1)
            ids(1)=ids(last)
            ids(last)=t
            call sift_integer(ids,1,last-1)
        end do
    end subroutine sort_indices

    subroutine sift_integer(ids,first,last)
        integer, intent(inout) :: ids(:) !! integer max heap
        integer, intent(in) :: first !! node to sift
        integer, intent(in) :: last !! valid heap length
        integer :: parent,child,tmp
        parent=first
        do while (parent<=last/2)
            child=2*parent
            if (child<last) then
                if (ids(child+1)>ids(child)) child=child+1
            end if
            if (ids(child)<=ids(parent)) exit
            tmp=ids(parent)
            ids(parent)=ids(child)
            ids(child)=tmp
            parent=child
        end do
    end subroutine sift_integer

end module scifort_interpolate_kdtree
