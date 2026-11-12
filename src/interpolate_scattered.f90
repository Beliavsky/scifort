! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Scattered-data interpolation: nearest neighbor in any positive dimension and
! Delaunay/barycentric piecewise linear interpolation in native 2D-5D, or any
! geometric dimension with an externally supplied simplicial tessellation.
! Independent implementation (Bowyer-Watson); no Qhull, BLAS, or LAPACK needed.
! SciPy is used as an independent numerical test reference.
module scifort_interpolate_scattered
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: int64
    use scifort_kinds, only: dp
    use scifort_interpolate_kdtree, only: kd_tree, kd_tree_type, kdtree_success
    implicit none
    private

    integer, parameter, public :: scattered_status_success = 0
    integer, parameter, public :: scattered_status_invalid_input = 1
    integer, parameter, public :: scattered_status_degenerate = 2

    ! Balanced axis-aligned bounding-box hierarchy over 3-D tetrahedra.
    ! Leaf ranges refer to a permutation of tetrahedron indices.  This is a
    ! candidate index only: final containment uses the original extended-
    ! precision tetrahedron_barycentric calculation (no matrix-inverse cache).
    ! High-precision kind for small determinant systems (uses real128 on
    ! GNU Fortran, or real64 when the compiler does not provide it).
    integer, parameter :: geom_wanted_kind = selected_real_kind(33,4931)
    integer, parameter :: geom_kind = merge(geom_wanted_kind,dp,geom_wanted_kind>0)

    integer, parameter :: tetra_bvh_leaf_size = 8
    integer, parameter :: tetra_bvh_stack_size = 128

    type :: tetra_bvh_type
        integer :: node_count = 0
        integer, allocatable :: order(:)             ! tetra indices in leaf order
        integer, allocatable :: child_l(:), child_r(:)
        integer, allocatable :: first(:), last(:)    ! leaf ranges; zero on branches
        real(dp), allocatable :: bbox_min(:,:), bbox_max(:,:) ! (3,nodes)
        real(dp), allocatable :: tetra_min(:,:), tetra_max(:,:) ! (3,ntetra)
        real(dp) :: margin(3) = 0.0_dp              ! roundoff/barycentric boundary slack
    end type tetra_bvh_type

    ! Dimension-independent bounding-box hierarchy for supplied simplices
    ! and for natively triangulated 4D/5D clouds.
    type :: simplex_bvh_type
        integer :: node_count = 0
        integer, allocatable :: order(:), child_l(:), child_r(:), first(:), last(:)
        real(dp), allocatable :: bbox_min(:,:), bbox_max(:,:), margin(:)
    end type simplex_bvh_type

    type :: scattered_core
        type(kd_tree_type) :: spatial_index             ! reusable balanced nearest-neighbor tree
        real(dp), allocatable :: points(:,:)            ! nsite x dimension, normalized if rescale
        complex(dp), allocatable :: values(:,:)         ! nsite x channels
        real(dp), allocatable :: offset(:), scale(:)    ! query transformation
        integer, allocatable :: triangles(:,:)          ! 3 x ntriangles for 2D sites
        integer, allocatable :: tetrahedra(:,:)         ! 4 x ntetrahedra for 3D sites
        integer, allocatable :: external_simplices(:,:) ! (d+1,nsimplex), any d>=2
        real(dp), allocatable :: simplex_min(:,:), simplex_max(:,:) ! (d,nsimplex)
        type(tetra_bvh_type) :: tetra_index            ! 3D containing-tetrahedron index
        type(simplex_bvh_type) :: simplex_index         ! rank-generic containing-simplex index
        integer, allocatable :: simplex_neighbors(:,:), site_simplex(:)
        real(geom_kind), allocatable :: simplex_inv(:,:,:), simplex_base(:,:)
        logical, allocatable :: simplex_inverse_valid(:)
        logical :: linear = .false.
        complex(dp) :: fill = (0.0_dp, 0.0_dp)
    end type scattered_core

    type, public :: nearest_nd_interpolator_type
        private
        type(scattered_core) :: core
    contains
        procedure :: valid => nearest_real_valid
        procedure :: dimensions => nearest_real_dimensions
        procedure :: channels => nearest_real_channels
        procedure :: value => nearest_real_value
        procedure, private :: nearest_real_one, nearest_real_many
        generic :: evaluate => nearest_real_one, nearest_real_many
    end type nearest_nd_interpolator_type

    type, public :: complex_nearest_nd_interpolator_type
        private
        type(scattered_core) :: core
    contains
        procedure :: valid => nearest_complex_valid
        procedure :: dimensions => nearest_complex_dimensions
        procedure :: channels => nearest_complex_channels
        procedure :: value => nearest_complex_value
        procedure, private :: nearest_complex_one, nearest_complex_many
        generic :: evaluate => nearest_complex_one, nearest_complex_many
    end type complex_nearest_nd_interpolator_type

    type, public :: linear_nd_interpolator_type
        private
        type(scattered_core) :: core
    contains
        procedure :: valid => linear_real_valid
        procedure :: dimensions => linear_real_dimensions
        procedure :: channels => linear_real_channels
        procedure :: value => linear_real_value
        procedure :: simplices => linear_real_simplices
        procedure, private :: linear_real_one, linear_real_many
        generic :: evaluate => linear_real_one, linear_real_many
    end type linear_nd_interpolator_type

    type, public :: complex_linear_nd_interpolator_type
        private
        type(scattered_core) :: core
    contains
        procedure :: valid => linear_complex_valid
        procedure :: dimensions => linear_complex_dimensions
        procedure :: channels => linear_complex_channels
        procedure :: value => linear_complex_value
        procedure :: simplices => linear_complex_simplices
        procedure, private :: linear_complex_one, linear_complex_many
        generic :: evaluate => linear_complex_one, linear_complex_many
    end type complex_linear_nd_interpolator_type

    public :: nearest_nd_interpolator, linear_nd_interpolator, griddata
    public :: delaunay_2d, delaunay_3d, delaunay_nd
    interface griddata
        module procedure griddata_real_vector, griddata_real_matrix
        module procedure griddata_complex_vector, griddata_complex_matrix
    end interface
    interface nearest_nd_interpolator
        module procedure nearest_real_vector_constructor, nearest_real_matrix_constructor
        module procedure nearest_complex_vector_constructor, nearest_complex_matrix_constructor
    end interface
    interface linear_nd_interpolator
        module procedure linear_real_vector_constructor, linear_real_matrix_constructor
        module procedure linear_complex_vector_constructor, linear_complex_matrix_constructor
    end interface

contains

    pure real(dp) function nan_real() result(v)
        v = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_real

    pure logical function core_valid(p) result(ok)
        type(scattered_core), intent(in) :: p !! coordinate and ordinate data
        ok = allocated(p%points) .and. allocated(p%values) .and. allocated(p%offset) .and. allocated(p%scale)
        if (.not. ok) return
        ok = size(p%points,1)>0 .and. size(p%points,2)>0 .and. size(p%values,1)==size(p%points,1)
        if (.not. ok) return
        if (p%linear) then
            if (allocated(p%external_simplices)) then
                ok=size(p%external_simplices,1)==size(p%points,2)+1 .and. &
                   size(p%external_simplices,2)>0
                ok=ok .and. p%simplex_index%node_count>0
                return
            end if
            select case(size(p%points,2))
            case(2)
                ok = allocated(p%triangles)
            case(3)
                ok = allocated(p%tetrahedra) .and. p%tetra_index%node_count>0
            case default
                ok = .false.
            end select
        else
            ok = p%spatial_index%valid()
        end if
    end function core_valid

    pure logical function nearest_real_valid(self) result(ok)
        class(nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        ok=core_valid(self%core)
    end function nearest_real_valid

    pure logical function nearest_complex_valid(self) result(ok)
        class(complex_nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        ok=core_valid(self%core)
    end function nearest_complex_valid

    pure logical function linear_real_valid(self) result(ok)
        class(linear_nd_interpolator_type), intent(in) :: self !! interpolator
        ok=core_valid(self%core)
    end function linear_real_valid

    pure logical function linear_complex_valid(self) result(ok)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! interpolator
        ok=core_valid(self%core)
    end function linear_complex_valid

    pure integer function nearest_real_dimensions(self) result(n)
        class(nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%points)) n=size(self%core%points,2)
    end function nearest_real_dimensions

    pure integer function nearest_complex_dimensions(self) result(n)
        class(complex_nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%points)) n=size(self%core%points,2)
    end function nearest_complex_dimensions

    pure integer function linear_real_dimensions(self) result(n)
        class(linear_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%points)) n=size(self%core%points,2)
    end function linear_real_dimensions

    pure integer function linear_complex_dimensions(self) result(n)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%points)) n=size(self%core%points,2)
    end function linear_complex_dimensions

    pure integer function nearest_real_channels(self) result(n)
        class(nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%values)) n=size(self%core%values,2)
    end function nearest_real_channels

    pure integer function nearest_complex_channels(self) result(n)
        class(complex_nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%values)) n=size(self%core%values,2)
    end function nearest_complex_channels

    pure integer function linear_real_channels(self) result(n)
        class(linear_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%values)) n=size(self%core%values,2)
    end function linear_real_channels

    pure integer function linear_complex_channels(self) result(n)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%values)) n=size(self%core%values,2)
    end function linear_complex_channels

    function linear_real_simplices(self) result(v)
        class(linear_nd_interpolator_type), intent(in) :: self !! triangulation
        integer, allocatable :: v(:,:)
        allocate(v(3,0))
        if (.not. self%valid()) return
        if (allocated(self%core%external_simplices)) then
            v=self%core%external_simplices
        else if (allocated(self%core%tetrahedra)) then
            v=self%core%tetrahedra
        else
            v=self%core%triangles
        end if
    end function linear_real_simplices

    function linear_complex_simplices(self) result(v)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! triangulation
        integer, allocatable :: v(:,:)
        allocate(v(3,0))
        if (.not. self%valid()) return
        if (allocated(self%core%external_simplices)) then
            v=self%core%external_simplices
        else if (allocated(self%core%tetrahedra)) then
            v=self%core%tetrahedra
        else
            v=self%core%triangles
        end if
    end function linear_complex_simplices

    function nearest_real_vector_constructor(x,y,rescale,status) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions), at least one site
        real(dp), intent(in) :: y(:) !! value at every site
        logical, intent(in), optional :: rescale !! rescale each coordinate to unit range
        integer, intent(out), optional :: status !! status code
        type(nearest_nd_interpolator_type) :: p
        p=nearest_real_matrix_constructor(x,reshape(y,[size(y),1]),rescale,status)
    end function nearest_real_vector_constructor

    function nearest_real_matrix_constructor(x,y,rescale,status) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        real(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in), optional :: rescale !! normalize coordinate axes
        integer, intent(out), optional :: status !! status code
        type(nearest_nd_interpolator_type) :: p
        integer :: code
        call initialize(p%core,x,cmplx(y,0.0_dp,dp),.false.,rescale,code)
        if (present(status)) status=code
    end function nearest_real_matrix_constructor

    function nearest_complex_vector_constructor(x,y,rescale,status) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        complex(dp), intent(in) :: y(:) !! (sites)
        logical, intent(in), optional :: rescale !! normalize coordinate axes
        integer, intent(out), optional :: status !! status code
        type(complex_nearest_nd_interpolator_type) :: p
        p=nearest_complex_matrix_constructor(x,reshape(y,[size(y),1]),rescale,status)
    end function nearest_complex_vector_constructor

    function nearest_complex_matrix_constructor(x,y,rescale,status) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        complex(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in), optional :: rescale !! normalize coordinate axes
        integer, intent(out), optional :: status !! status code
        type(complex_nearest_nd_interpolator_type) :: p
        integer :: code
        call initialize(p%core,x,y,.false.,rescale,code)
        if (present(status)) status=code
    end function nearest_complex_matrix_constructor

    function linear_real_vector_constructor(x,y,rescale,fill_value,status,simplices) result(p)
        real(dp), intent(in) :: x(:,:) !! scattered sites (n,ndim)
        real(dp), intent(in) :: y(:) !! (sites)
        logical, intent(in), optional :: rescale !! normalize axes
        real(dp), intent(in), optional :: fill_value !! outside-convex-hull value; default NaN
        integer, intent(out), optional :: status !! status code
        integer, intent(in), optional :: simplices(:,:) !! (d+1,nsimplex) externally supplied tessellation
        type(linear_nd_interpolator_type) :: p
        p=linear_real_matrix_constructor(x,reshape(y,[size(y),1]),rescale,fill_value,status,simplices)
    end function linear_real_vector_constructor

    function linear_real_matrix_constructor(x,y,rescale,fill_value,status,simplices) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,ndim)
        real(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in), optional :: rescale !! normalize axes
        real(dp), intent(in), optional :: fill_value !! value beyond convex hull
        integer, intent(out), optional :: status !! status code
        integer, intent(in), optional :: simplices(:,:) !! (d+1,nsimplex) externally supplied tessellation
        type(linear_nd_interpolator_type) :: p
        integer :: code
        call initialize(p%core,x,cmplx(y,0.0_dp,dp),.true.,rescale,code,simplices)
        if (present(fill_value)) p%core%fill=cmplx(fill_value,0.0_dp,dp)
        if (present(status)) status=code
    end function linear_real_matrix_constructor

    function linear_complex_vector_constructor(x,y,rescale,fill_value,status,simplices) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,ndim)
        complex(dp), intent(in) :: y(:) !! (sites)
        logical, intent(in), optional :: rescale !! normalize axes
        complex(dp), intent(in), optional :: fill_value !! outside-convex-hull value
        integer, intent(out), optional :: status !! status code
        integer, intent(in), optional :: simplices(:,:) !! (d+1,nsimplex) externally supplied tessellation
        type(complex_linear_nd_interpolator_type) :: p
        p=linear_complex_matrix_constructor(x,reshape(y,[size(y),1]),rescale,fill_value,status,simplices)
    end function linear_complex_vector_constructor

    function linear_complex_matrix_constructor(x,y,rescale,fill_value,status,simplices) result(p)
        real(dp), intent(in) :: x(:,:) !! (sites,ndim)
        complex(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in), optional :: rescale !! normalize axes
        complex(dp), intent(in), optional :: fill_value !! value beyond convex hull
        integer, intent(out), optional :: status !! status code
        integer, intent(in), optional :: simplices(:,:) !! (d+1,nsimplex) externally supplied tessellation
        type(complex_linear_nd_interpolator_type) :: p
        integer :: code
        call initialize(p%core,x,y,.true.,rescale,code,simplices)
        if (present(fill_value)) p%core%fill=fill_value
        if (present(status)) status=code
    end function linear_complex_matrix_constructor

    subroutine initialize(p,x,y,linear,rescale,code,simplices)
        type(scattered_core), intent(out) :: p !! constructed interpolator core
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        complex(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in) :: linear !! request simplicial interpolation; native Delaunay in 2D-5D
        logical, intent(in), optional :: rescale !! normalize each coordinate axis
        integer, intent(out) :: code !! result status
        integer, intent(in), optional :: simplices(:,:) !! given (d+1,nsimplex), skips tessellation
        integer :: n,d,j,i
        real(dp), allocatable :: weights(:)
        logical :: nonsingular
        real(dp) :: dd
        code=scattered_status_invalid_input
        n=size(x,1)
        d=size(x,2)
        if (n<1 .or. d<1 .or. size(y,1)/=n .or. size(y,2)<1) return
        if (linear) then
            if (n<d+1) return
            if (present(simplices)) then
                if (d<2 .or. size(simplices,1)/=d+1 .or. size(simplices,2)<1) return
                if (any(simplices<1) .or. any(simplices>n)) return
            else
                if (d<2 .or. d>5) return
            end if
        end if
        if (.not.all(ieee_is_finite(x))) return
        p%linear=linear
        allocate(p%points(n,d),p%values(n,size(y,2)),p%offset(d),p%scale(d))
        p%offset=0.0_dp
        p%scale=1.0_dp
        if (present(rescale)) then
            if (rescale) then
                p%offset=minval(x,dim=1)
                p%scale=maxval(x,dim=1)-p%offset
                if (linear .and. any(p%scale==0.0_dp)) goto 900
                where(p%scale==0.0_dp) p%scale=1.0_dp
            end if
        end if
        p%points=(x-spread(p%offset,1,n))/spread(p%scale,1,n)
        p%values=y
        p%fill=cmplx(nan_real(),0.0_dp,dp)
        if (.not.linear) then
            p%spatial_index=kd_tree(p%points,status=code)
            if (code/=kdtree_success) goto 900
        end if
        if (linear) then
            ! Input must be unique; Qhull silently removes duplicates, which
            ! can obscure which ordinate should be selected at a duplicate site.
            do i=1,n
                do j=i+1,n
                    dd=sum((p%points(i,:)-p%points(j,:))**2)
                    if (dd==0.0_dp) goto 900
                end do
            end do
            if (present(simplices)) then
                p%external_simplices=simplices
                allocate(weights(d+1))
                allocate(p%simplex_min(d,size(simplices,2)),p%simplex_max(d,size(simplices,2)))
                code=scattered_status_degenerate
                do j=1,size(simplices,2)
                    ! Reject repeated IDs and affinely dependent simplices.
                    do i=1,d+1
                        if (any(simplices(i+1:,j)==simplices(i,j))) goto 900
                    end do
                    p%simplex_min(:,j)=minval(p%points(simplices(:,j),:),dim=1)
                    p%simplex_max(:,j)=maxval(p%points(simplices(:,j),:),dim=1)
                    call simplex_barycentric(p%points(simplices(:,j),:), &
                        sum(p%points(simplices(:,j),:),dim=1)/real(d+1,dp), &
                        weights,nonsingular)
                    if (.not.nonsingular) goto 900
                end do
                call build_simplex_bvh(p)
            else
                if (d==2) then
                    call delaunay_2d(p%points,p%triangles,code)
                else if (d==3) then
                    call delaunay_3d(p%points,p%tetrahedra,code)
                else
                    call delaunay_nd(p%points,p%external_simplices,code)
                    if (code==scattered_status_success) then
                        allocate(p%simplex_min(d,size(p%external_simplices,2)), &
                                 p%simplex_max(d,size(p%external_simplices,2)))
                        do j=1,size(p%external_simplices,2)
                            p%simplex_min(:,j)=minval(p%points(p%external_simplices(:,j),:),dim=1)
                            p%simplex_max(:,j)=maxval(p%points(p%external_simplices(:,j),:),dim=1)
                        end do
                        call build_simplex_bvh(p)
                    end if
                end if
                if (code/=scattered_status_success) goto 900
                if (d==3) call build_tetra_bvh(p%points,p%tetrahedra,p%tetra_index)
            end if
        end if
        code=scattered_status_success
        return
900     if (allocated(p%points)) deallocate(p%points)
        if (allocated(p%values)) deallocate(p%values)
        if (allocated(p%offset)) deallocate(p%offset)
        if (allocated(p%scale)) deallocate(p%scale)
        if (allocated(p%triangles)) deallocate(p%triangles)
        if (allocated(p%tetrahedra)) deallocate(p%tetrahedra)
        if (allocated(p%external_simplices)) deallocate(p%external_simplices)
        if (allocated(p%simplex_min)) deallocate(p%simplex_min)
        if (allocated(p%simplex_max)) deallocate(p%simplex_max)
        if (code==scattered_status_success) code=scattered_status_invalid_input
    end subroutine initialize

    ! Build a balanced rank-two bounding-volume hierarchy. The coordinates
    ! are (nsite,ndim), and simplex bounds are (ndim,nsimplex); neither rank
    ! changes with geometric dimension.
    subroutine build_simplex_bvh(p)
        type(scattered_core), intent(inout) :: p !! interpolator being constructed
        real(dp), allocatable :: centers(:,:)
        integer :: d,nt,j,root
        d=size(p%points,2)
        nt=size(p%external_simplices,2)
        if (nt<1) return
        allocate(p%simplex_index%order(nt),p%simplex_index%child_l(2*nt), &
            p%simplex_index%child_r(2*nt),p%simplex_index%first(2*nt), &
            p%simplex_index%last(2*nt),p%simplex_index%bbox_min(d,2*nt), &
            p%simplex_index%bbox_max(d,2*nt),p%simplex_index%margin(d), &
            centers(d,nt))
        p%simplex_index%order=[(j,j=1,nt)]
        p%simplex_index%child_l=0
        p%simplex_index%child_r=0
        p%simplex_index%first=0
        p%simplex_index%last=0
        p%simplex_index%node_count=0
        do j=1,nt
            centers(:,j)=sum(transpose(p%points(p%external_simplices(:,j),:))/ &
                real(d+1,dp),dim=2)
        end do
        p%simplex_index%margin=max(2.0e-13_dp, &
            64.0_dp*epsilon(1.0_dp)*maxval(abs(p%points),dim=1))
        call build_simplex_bvh_node(p%simplex_index,centers,p%simplex_min, &
                                    p%simplex_max,1,nt,root)
        call build_simplex_neighbors(p)
    end subroutine build_simplex_bvh

    recursive subroutine build_simplex_bvh_node(index,centers,boxlo,boxhi,lo,hi,node)
        type(simplex_bvh_type), intent(inout) :: index !! index to populate
        real(dp), intent(in) :: centers(:,:) !! simplex centroids, (ndim,nsimplex)
        real(dp), intent(in) :: boxlo(:,:),boxhi(:,:) !! (ndim,nsimplex)
        integer, intent(in) :: lo,hi !! range of simplex permutation
        integer, intent(out) :: node !! assigned node ID
        integer :: mid,axis,left,right,i,j
        real(dp) :: lower(size(centers,1)),upper(size(centers,1)),span(size(centers,1))
        node=index%node_count+1
        index%node_count=node
        lower=huge(0.0_dp)
        upper=-huge(0.0_dp)
        do i=lo,hi
            j=index%order(i)
            lower=min(lower,boxlo(:,j))
            upper=max(upper,boxhi(:,j))
        end do
        index%bbox_min(:,node)=lower
        index%bbox_max(:,node)=upper
        if (hi-lo+1<=8) then
            index%first(node)=lo
            index%last(node)=hi
            return
        end if
        span=maxval(centers(:,index%order(lo:hi)),dim=2)- &
             minval(centers(:,index%order(lo:hi)),dim=2)
        axis=maxloc(span,dim=1)
        mid=lo+(hi-lo)/2
        call bvh_partition_median(index%order,centers,axis,lo,hi,mid)
        call build_simplex_bvh_node(index,centers,boxlo,boxhi,lo,mid,left)
        call build_simplex_bvh_node(index,centers,boxlo,boxhi,mid+1,hi,right)
        index%child_l(node)=left
        index%child_r(node)=right
    end subroutine build_simplex_bvh_node

    ! Precompute affine barycentric transforms and neighbors across facets.
    ! In high dimensions a BVH alone visits many overlapping bounding boxes;
    ! a deterministic neighbor walk started from the nearest observed site
    ! usually locates the simplex in a small number of steps.
    subroutine build_simplex_neighbors(p)
        type(scattered_core), intent(inout) :: p !! populated simplicial interpolator
        integer :: n,d,nt,j,k,m,i,position,h,cap,other,previous
        integer, allocatable :: keys(:,:),owner(:),facet(:),count(:)
        integer :: key(size(p%points,2)),sorted(size(p%points,2)),vert
        real(geom_kind) :: a(size(p%points,2),size(p%points,2))
        logical :: valid
        integer :: code
        n=size(p%points,1)
        d=size(p%points,2)
        nt=size(p%external_simplices,2)
        allocate(p%simplex_neighbors(d+1,nt),p%site_simplex(n), &
                 p%simplex_inv(d,d,nt),p%simplex_base(d,nt),p%simplex_inverse_valid(nt))
        p%simplex_neighbors=0
        p%site_simplex=0
        do j=1,nt
            do i=1,d+1
                m=p%external_simplices(i,j)
                if (p%site_simplex(m)==0) p%site_simplex(m)=j
            end do
            m=p%external_simplices(1,j)
            p%simplex_base(:,j)=real(p%points(m,:),geom_kind)
            do k=1,d
                a(:,k)=real(p%points(p%external_simplices(k+1,j),:),geom_kind)- &
                       p%simplex_base(:,j)
            end do
            call inverse_small_matrix(a,p%simplex_inv(:,:,j),valid)
            p%simplex_inverse_valid(j)=valid
        end do
        cap=2*(d+1)*nt+11
        allocate(keys(d,cap),owner(cap),facet(cap),count(cap))
        count=0
        do j=1,nt
            do k=1,d+1
                key=pack(p%external_simplices(:,j),[(m/=k,m=1,d+1)])
                sorted=key
                do i=2,d
                    vert=sorted(i)
                    position=i-1
                    do while(position>=1)
                        if (sorted(position)<=vert) exit
                        sorted(position+1)=sorted(position)
                        position=position-1
                    end do
                    sorted(position+1)=vert
                end do
                h=0
                do i=1,d
                    h=modulo(37*h+sorted(i),cap)
                end do
                h=h+1
                do
                    if (count(h)==0) then
                        keys(:,h)=sorted
                        owner(h)=j
                        facet(h)=k
                        count(h)=1
                        exit
                    else if(all(keys(:,h)==sorted)) then
                        if (count(h)==1) then
                            other=owner(h)
                            previous=facet(h)
                            p%simplex_neighbors(k,j)=other
                            p%simplex_neighbors(previous,other)=j
                        else
                            ! Non-manifold supplied tessellation: fall back
                            ! to the guaranteed exhaustive BVH search.
                            p%simplex_neighbors(k,j)=0
                        end if
                        count(h)=count(h)+1
                        exit
                    end if
                    h=modulo(h,cap)+1
                end do
            end do
        end do
        p%spatial_index=kd_tree(p%points,status=code)
    end subroutine build_simplex_neighbors

    pure subroutine inverse_small_matrix(a,ainv,valid)
        real(geom_kind), intent(in) :: a(:,:) !! nonsingular small square matrix
        real(geom_kind), intent(out) :: ainv(:,:) !! its inverse, if valid
        logical, intent(out) :: valid !! numerical nonsingularity
        real(geom_kind) :: m(size(a,1),size(a,2)),row(size(a,1)),temp,scale
        integer :: d,i,j,pivot
        d=size(a,1)
        valid=.false.
        m=a
        ainv=0.0_geom_kind
        do i=1,d
            ainv(i,i)=1.0_geom_kind
        end do
        scale=maxval(abs(m))
        if (scale<=0.0_geom_kind) return
        do j=1,d
            pivot=j-1+maxloc(abs(m(j:d,j)),dim=1)
            if (abs(m(pivot,j))<64.0_geom_kind*real(epsilon(1.0_dp),geom_kind)*scale) return
            if (pivot/=j) then
                row=m(j,:);m(j,:)=m(pivot,:);m(pivot,:)=row
                row=ainv(j,:);ainv(j,:)=ainv(pivot,:);ainv(pivot,:)=row
            end if
            temp=m(j,j)
            m(j,:)=m(j,:)/temp
            ainv(j,:)=ainv(j,:)/temp
            do i=1,d
                if (i==j) cycle
                temp=m(i,j)
                m(i,:)=m(i,:)-temp*m(j,:)
                ainv(i,:)=ainv(i,:)-temp*ainv(j,:)
            end do
        end do
        valid=.true.
    end subroutine inverse_small_matrix

    pure subroutine cached_barycentric(p,index,query,w,ok)
        type(scattered_core), intent(in) :: p !! initialized simplex cache
        integer, intent(in) :: index !! simplex number
        real(dp), intent(in) :: query(:) !! normalized location
        real(dp), intent(out) :: w(:) !! barycentric weights
        logical, intent(out) :: ok !! numerical validity
        real(geom_kind) :: v(size(query))
        if (p%simplex_inverse_valid(index)) then
            v=matmul(p%simplex_inv(:,:,index),real(query,geom_kind)-p%simplex_base(:,index))
            w(2:)=real(v,dp)
            w(1)=real(1.0_geom_kind-sum(v),dp)
            ok=all(ieee_is_finite(w))
        else
            call simplex_barycentric(p%points(p%external_simplices(:,index),:),query,w,ok)
        end if
    end subroutine cached_barycentric

    ! Choose the lowest-numbered containing simplex, even on shared faces,
    ! preserving the deterministic behavior of the original full scan.
    subroutine find_containing_simplex(p,query,best,weights,candidates)
        type(scattered_core), intent(in) :: p !! constructed interpolator core
        real(dp), intent(in) :: query(:) !! normalized point in ndim dimensions
        integer, intent(out) :: best !! containing simplex index or zero
        real(dp), intent(out) :: weights(:) !! barycentric coordinates (ndim+1)
        integer, intent(out), optional :: candidates !! exact affine evaluations performed
        integer :: stack(128),top,node,j,k,visited,near,start,step,next,side
        real(dp) :: w(size(weights))
        logical :: nonsingular
        best=0
        weights=0.0_dp
        visited=0
        if (p%simplex_index%node_count<1) then
            if (present(candidates)) candidates=0
            return
        end if
        near=0
        if (p%spatial_index%valid()) near=p%spatial_index%nearest(query)
        start=0
        if (near>0) start=p%site_simplex(near)
        if (start>0) then
            j=start
            do step=1,min(128,size(p%external_simplices,2))
                visited=visited+1
                call cached_barycentric(p,j,query,w,nonsingular)
                if (.not.nonsingular) exit
                if (minval(w)>=-2.0e-13_dp .and. maxval(w)<=1.0_dp+2.0e-13_dp) then
                    best=j
                    weights=w
                    if (present(candidates)) candidates=visited
                    return
                end if
                side=minloc(w,dim=1)
                next=p%simplex_neighbors(side,j)
                if (next==0 .or. next==j) exit
                j=next
            end do
        end if
        ! A supplied triangulation may be disconnected or non-manifold, and
        ! near-degenerate query positions can defeat the directed walk.
        ! Fall back to the BVH and the original boundary tolerances.
        top=1
        stack(1)=1
        do while(top>0)
            node=stack(top)
            top=top-1
            if (any(query<p%simplex_index%bbox_min(:,node)-p%simplex_index%margin)) cycle
            if (any(query>p%simplex_index%bbox_max(:,node)+p%simplex_index%margin)) cycle
            if (p%simplex_index%child_l(node)>0) then
                if (top+2>size(stack)) error stop 'generic simplex BVH stack exceeded'
                stack(top+1)=p%simplex_index%child_r(node)
                stack(top+2)=p%simplex_index%child_l(node)
                top=top+2
            else
                do k=p%simplex_index%first(node),p%simplex_index%last(node)
                    j=p%simplex_index%order(k)
                    if (best>0 .and. j>=best) cycle
                    if (any(query<p%simplex_min(:,j)-p%simplex_index%margin)) cycle
                    if (any(query>p%simplex_max(:,j)+p%simplex_index%margin)) cycle
                    visited=visited+1
                    call cached_barycentric(p,j,query,w,nonsingular)
                    if (.not.nonsingular) cycle
                    if (minval(w)<-2.0e-13_dp) cycle
                    if (maxval(w)>1.0_dp+2.0e-13_dp) cycle
                    best=j
                    weights=w
                end do
            end if
        end do
        if (present(candidates)) candidates=visited
    end subroutine find_containing_simplex

    ! The 3D Delaunay implementation is an independent incremental
    ! Bowyer-Watson algorithm.  Its predicates are evaluated in extended
    ! precision where supported, but are NOT adaptive exact predicates as
    ! used by Qhull. Near-degenerate/cospherical data may be rejected.

    ! Geometry-independent Delaunay simplices in four or five dimensions.
    ! This is an incremental hull-preserving Bowyer-Watson construction.
    ! Geometric coordinates always have rank two, (sites,ndim); no
    ! dimension-specific assumed-rank SELECT RANK dispatch is needed.
    !
    ! This is not a replacement for Qhull's adaptive exact predicates.
    ! Degenerate/cospherical configurations and very large clouds require
    ! a robust external triangulation supplied with simplices=.
    subroutine delaunay_nd(x,simplices,code)
        real(dp), intent(in) :: x(:,:) !! (npoints,ndim), unique finite coordinates
        integer, allocatable, intent(out) :: simplices(:,:) !! (ndim+1,ncells), 1-based vertex indices
        integer, intent(out) :: code !! scattered_status_* code
        integer :: n,d,i,j,k,ns,capacity,nfaces,nleft,face_cap
        integer :: pick,seed_count
        integer, allocatable :: cell(:,:),boundary(:,:),ids(:),seed(:)
        logical, allocatable :: bad(:),inserted(:)
        real(dp), allocatable :: xnorm(:,:)
        real(geom_kind), allocatable :: center(:,:),radius(:)
        real(geom_kind) :: basis(size(x,2),size(x,2)),residual(size(x,2)), &
             axis(size(x,2)), norm2,best,det,scale,shift(size(x,2))
        logical :: ok
        code=scattered_status_invalid_input
        n=size(x,1)
        d=size(x,2)
        if (d==2) then
            call delaunay_2d(x,simplices,code)
            return
        else if (d==3) then
            call delaunay_3d(x,simplices,code)
            return
        end if
        if (d<4 .or. d>5 .or. n<d+1 .or. n>160) return
        if (.not.all(ieee_is_finite(x))) return
        do i=1,n
            do j=i+1,n
                if (all(x(i,:)==x(j,:))) return
            end do
        end do
        code=scattered_status_degenerate
        allocate(xnorm(n,d),seed(d+1),inserted(n))
        shift=0.5_geom_kind*(real(minval(x,dim=1),geom_kind)+ &
                               real(maxval(x,dim=1),geom_kind))
        scale=maxval(real(maxval(x,dim=1)-minval(x,dim=1),geom_kind))
        if (scale<=0.0_geom_kind) return
        xnorm=real((real(x,geom_kind)-spread(shift,1,n))/scale,dp)
        ! Greedy Gram-Schmidt chooses a well-conditioned initial simplex.
        seed=0
        seed(1)=1
        basis=0.0_geom_kind
        do seed_count=1,d
            best=0.0_geom_kind
            pick=0
            do i=1,n
                if (any(seed(:seed_count)==i)) cycle
                residual=real(xnorm(i,:)-xnorm(seed(1),:),geom_kind)
                do j=1,seed_count-1
                    residual=residual-dot_product(residual,basis(:,j))*basis(:,j)
                end do
                norm2=dot_product(residual,residual)
                if (norm2>best) then
                    best=norm2
                    pick=i
                end if
            end do
            if (pick==0 .or. best<=1.0e-24_geom_kind) return
            seed(seed_count+1)=pick
            residual=real(xnorm(pick,:)-xnorm(seed(1),:),geom_kind)
            do j=1,seed_count-1
                residual=residual-dot_product(residual,basis(:,j))*basis(:,j)
            end do
            basis(:,seed_count)=residual/sqrt(dot_product(residual,residual))
        end do
        det=nd_orientation(xnorm(seed,:))
        if (det==0.0_geom_kind) return
        if (det<0.0_geom_kind) then
            pick=seed(1);seed(1)=seed(2);seed(2)=pick
        end if
        ! Deliberately cap the initial high-dimensional implementation.
        ! A capacity overrun is returned as a status, never ERROR STOP.
        capacity=min(250000,max(2048,16*n*n))
        face_cap=(d+1)*capacity
        allocate(cell(d+1,capacity),center(d,capacity),radius(capacity), &
                 bad(capacity),boundary(d,face_cap),ids(d+1))
        cell(:,1)=seed
        call nd_circumsphere(xnorm,seed,center(:,1),radius(1),ok)
        if (.not.ok) return
        ns=1
        inserted=.false.
        inserted(seed)=.true.
        do i=1,n
            if (inserted(i)) cycle
            bad(:ns)=.false.
            nfaces=0
            do j=1,ns
                ! Exact-predicate filter and fallback both use extended
                ! precision. Cache prevents repeated matrix factorization.
                axis=real(xnorm(i,:),geom_kind)-center(:,j)
                if (dot_product(axis,axis)>=radius(j)) cycle
                bad(j)=.true.
                do k=1,d+1
                    call nd_toggle_face(boundary,nfaces,pack(cell(:,j), &
                        [(pick/=k,pick=1,d+1)]))
                end do
            end do
            call nd_visible_hull_faces(xnorm,cell,ns,xnorm(i,:),boundary,nfaces)
            ! An exactly cospherical interior point might not be in any
            ! strict cavity; split its containing cell as a fallback.
            if (.not.any(bad(:ns)) .and. nfaces==0) then
                block
                    real(dp) :: weights(d+1)
                    do j=1,ns
                        call simplex_barycentric(xnorm(cell(:,j),:),xnorm(i,:),weights,ok)
                        if (.not.ok) cycle
                        if (minval(weights)<-1.0e-13_dp) cycle
                        bad(j)=.true.
                        do k=1,d+1
                            call nd_toggle_face(boundary,nfaces,pack(cell(:,j), &
                                [(pick/=k,pick=1,d+1)]))
                        end do
                        exit
                    end do
                end block
            end if
            nleft=0
            do j=1,ns
                if (bad(j)) cycle
                nleft=nleft+1
                cell(:,nleft)=cell(:,j)
                center(:,nleft)=center(:,j)
                radius(nleft)=radius(j)
            end do
            ns=nleft
            do j=1,nfaces
                if (boundary(1,j)==0) cycle
                ids(:d)=boundary(:,j)
                ids(d+1)=i
                det=nd_orientation(xnorm(ids,:))
                if (abs(det)<=1.0e-29_geom_kind) cycle
                if (det<0.0_geom_kind) then
                    pick=ids(1);ids(1)=ids(2);ids(2)=pick
                end if
                if (ns>=capacity) return
                call nd_circumsphere(xnorm,ids,center(:,ns+1),radius(ns+1),ok)
                if (.not.ok) cycle
                ns=ns+1
                cell(:,ns)=ids
            end do
            if (ns==0) return
            inserted(i)=.true.
        end do
        allocate(simplices(d+1,ns))
        simplices=cell(:,:ns)
        code=scattered_status_success
    end subroutine delaunay_nd

    pure function nd_orientation(vertices) result(det)
        real(dp), intent(in) :: vertices(:,:) !! (ndim+1,ndim) affine simplex
        real(geom_kind) :: det
        real(geom_kind) :: a(size(vertices,2),size(vertices,2))
        real(geom_kind) :: row(size(vertices,2)),v
        integer :: d,i,j,pivot
        d=size(vertices,2)
        det=1.0_geom_kind
        do j=1,d
            a(:,j)=real(vertices(j+1,:)-vertices(1,:),geom_kind)
        end do
        do j=1,d
            pivot=j-1+maxloc(abs(a(j:d,j)),dim=1)
            if (a(pivot,j)==0.0_geom_kind) then
                det=0.0_geom_kind
                return
            end if
            if (pivot/=j) then
                row=a(j,:);a(j,:)=a(pivot,:);a(pivot,:)=row
                det=-det
            end if
            v=a(j,j)
            det=det*v
            do i=j+1,d
                a(i,j+1:d)=a(i,j+1:d)-a(i,j)/v*a(j,j+1:d)
            end do
        end do
    end function nd_orientation

    pure subroutine nd_circumsphere(points,ids,center,radius2,ok)
        real(dp), intent(in) :: points(:,:) !! (sites,ndim)
        integer, intent(in) :: ids(:) !! ndim+1 vertex indices
        real(geom_kind), intent(out) :: center(:) !! center of circumsphere
        real(geom_kind), intent(out) :: radius2 !! squared radius
        logical, intent(out) :: ok !! nonsingular cell
        real(geom_kind) :: a(size(center),size(center)),rhs(size(center)), &
             row(size(center)),origin(size(center)),delta(size(center)),v,tmp
        integer :: d,i,j,pivot
        d=size(center)
        origin=real(points(ids(1),:),geom_kind)
        do j=1,d
            delta=real(points(ids(j+1),:),geom_kind)-origin
            a(j,:)=delta
            rhs(j)=0.5_geom_kind*dot_product(delta,delta)
        end do
        ok=.false.
        do j=1,d
            pivot=j-1+maxloc(abs(a(j:d,j)),dim=1)
            if (abs(a(pivot,j))<max(1.0e-29_geom_kind, &
                    16.0_geom_kind*epsilon(1.0_geom_kind)*maxval(abs(a)))) return
            if (pivot/=j) then
                row=a(j,:);a(j,:)=a(pivot,:);a(pivot,:)=row
                tmp=rhs(j);rhs(j)=rhs(pivot);rhs(pivot)=tmp
            end if
            v=a(j,j)
            do i=j+1,d
                tmp=a(i,j)/v
                a(i,j+1:d)=a(i,j+1:d)-tmp*a(j,j+1:d)
                rhs(i)=rhs(i)-tmp*rhs(j)
            end do
        end do
        do i=d,1,-1
            rhs(i)=(rhs(i)-dot_product(a(i,i+1:d),rhs(i+1:d)))/a(i,i)
        end do
        center=origin+rhs
        radius2=dot_product(rhs,rhs)
        ok=.true.
    end subroutine nd_circumsphere

    subroutine nd_toggle_face(faces,nface,face)
        integer, intent(inout) :: faces(:,:) !! (ndim,capacity)
        integer, intent(inout) :: nface !! used entries including canceled slots
        integer, intent(in) :: face(:) !! unordered facet vertices
        integer :: j,k
        do j=1,nface
            if (faces(1,j)==0) cycle
            if (all([(any(faces(:,j)==face(k)),k=1,size(face))])) then
                faces(:,j)=0
                return
            end if
        end do
        nface=nface+1
        if (nface>size(faces,2)) error stop 'nd_toggle_face capacity exceeded'
        faces(:,nface)=face
    end subroutine nd_toggle_face

    subroutine nd_visible_hull_faces(points,cells,ns,p,faces,nfaces)
        real(dp), intent(in) :: points(:,:) !! normalized sites, (n,d)
        integer, intent(in) :: cells(:,:) !! (d+1,capacity)
        integer, intent(in) :: ns !! current cell count
        real(dp), intent(in) :: p(:) !! next site
        integer, intent(inout) :: faces(:,:) !! cavity boundary faces
        integer, intent(inout) :: nfaces !! occupied or canceled face entries
        integer, allocatable :: keys(:,:),unsorted(:,:),opposite(:),count(:)
        integer :: d,cap,j,k,h,m,idx,vert
        integer :: sorted(size(p)),face(size(p))
        real(geom_kind) :: op_side,p_side
        real(dp) :: vertices(size(p)+1,size(p))
        d=size(p)
        cap=2*(d+1)*ns+11
        allocate(keys(d,cap),unsorted(d,cap),opposite(cap),count(cap))
        count=0
        do j=1,ns
            do k=1,d+1
                face=pack(cells(:,j),[(m/=k,m=1,d+1)])
                sorted=face
                do m=2,d
                    vert=sorted(m)
                    idx=m-1
                    do while(idx>=1)
                        if (sorted(idx)<=vert) exit
                        sorted(idx+1)=sorted(idx)
                        idx=idx-1
                    end do
                    sorted(idx+1)=vert
                end do
                h=0
                do m=1,d
                    h=modulo(37*h+sorted(m),cap)
                end do
                h=h+1
                do
                    if (count(h)==0) then
                        count(h)=1
                        keys(:,h)=sorted
                        unsorted(:,h)=face
                        opposite(h)=cells(k,j)
                        exit
                    else if(all(keys(:,h)==sorted)) then
                        count(h)=count(h)+1
                        exit
                    end if
                    h=modulo(h,cap)+1
                end do
            end do
        end do
        do h=1,cap
            if (count(h)/=1) cycle
            face=unsorted(:,h)
            vertices(:d,:)=points(face,:)
            vertices(d+1,:)=points(opposite(h),:)
            op_side=nd_orientation(vertices)
            vertices(d+1,:)=p
            p_side=nd_orientation(vertices)
            if (p_side*op_side<0.0_geom_kind) call nd_toggle_face(faces,nfaces,face)
        end do
    end subroutine nd_visible_hull_faces

    pure function det3q(a,b,c) result(v)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(qkind), intent(in) :: a(3) !! first column
        real(qkind), intent(in) :: b(3) !! second column
        real(qkind), intent(in) :: c(3) !! third column
        real(qkind) :: v
        v=a(1)*(b(2)*c(3)-b(3)*c(2)) + &
          a(2)*(b(3)*c(1)-b(1)*c(3)) + &
          a(3)*(b(1)*c(2)-b(2)*c(1))
    end function det3q

    pure function orient3q(a,b,c,d) result(v)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(dp), intent(in) :: a(3) !! first vertex
        real(dp), intent(in) :: b(3) !! second vertex
        real(dp), intent(in) :: c(3) !! third vertex
        real(dp), intent(in) :: d(3) !! fourth vertex
        real(qkind) :: v
        v=det3q(real(b,qkind)-real(a,qkind), &
                real(c,qkind)-real(a,qkind), &
                real(d,qkind)-real(a,qkind))
    end function orient3q

    pure logical function insphere3d(a,b,c,d,p) result(inside)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(dp), intent(in) :: a(3) !! first vertex
        real(dp), intent(in) :: b(3) !! second vertex
        real(dp), intent(in) :: c(3) !! third vertex
        real(dp), intent(in) :: d(3) !! fourth vertex
        real(dp), intent(in) :: p(3) !! candidate site
        real(qkind) :: r(3,4),norm(4),det,orientation
        r(:,1)=real(a,qkind)-real(p,qkind)
        r(:,2)=real(b,qkind)-real(p,qkind)
        r(:,3)=real(c,qkind)-real(p,qkind)
        r(:,4)=real(d,qkind)-real(p,qkind)
        norm=sum(r*r,dim=1)
        det=-norm(1)*det3q(r(:,2),r(:,3),r(:,4)) &
            +norm(2)*det3q(r(:,1),r(:,3),r(:,4)) &
            -norm(3)*det3q(r(:,1),r(:,2),r(:,4)) &
            +norm(4)*det3q(r(:,1),r(:,2),r(:,3))
        orientation=orient3q(a,b,c,d)
        inside=det*orientation<0.0_qkind
    end function insphere3d

    ! Cache the circumsphere of a Delaunay tetrahedron.  During incremental
    ! construction this inexpensive conservative prefilter rejects most
    ! tetrahedra before the extended-precision in-sphere determinant.  The
    ! center is formed in qkind to avoid gratuitous cancellation; borderline
    ! points always go through the original determinant predicate.
    pure subroutine tetra_circumsphere(a,b,c,d,center,radius2)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(dp), intent(in) :: a(3) !! first site
        real(dp), intent(in) :: b(3) !! second site
        real(dp), intent(in) :: c(3) !! third site
        real(dp), intent(in) :: d(3) !! fourth site
        real(dp), intent(out) :: center(3) !! circumsphere center
        real(dp), intent(out) :: radius2 !! squared circumsphere radius
        real(qkind) :: u(3),v(3),w(3),origin(3),cc(3),det,bb(3)
        origin=real(a,qkind)
        u=real(b,qkind)-origin
        v=real(c,qkind)-origin
        w=real(d,qkind)-origin
        det=det3q(u,v,w)
        if (det==0.0_qkind) then
            center=0.0_dp
            radius2=huge(1.0_dp)
            return
        end if
        bb=0.5_qkind*[dot_product(u,u),dot_product(v,v),dot_product(w,w)]
        cc=origin+(bb(1)*cross3q(v,w)+bb(2)*cross3q(w,u)+ &
                   bb(3)*cross3q(u,v))/det
        center=real(cc,dp)
        radius2=real(dot_product(cc-origin,cc-origin),dp)
        if (.not.all(ieee_is_finite(center)) .or. .not.ieee_is_finite(radius2)) then
            center=0.0_dp
            radius2=huge(1.0_dp)
        end if
    end subroutine tetra_circumsphere

    pure function cross3q(a,b) result(c)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(qkind), intent(in) :: a(3) !! first vector
        real(qkind), intent(in) :: b(3) !! second vector
        real(qkind) :: c(3)
        c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    end function cross3q

    ! True only for points safely outside the approximate circumsphere.
    ! Near its boundary, never skip the exact in-sphere predicate.
    pure logical function outside_circumsphere(point,center,radius2) result(outside)
        real(dp), intent(in) :: point(3) !! candidate site
        real(dp), intent(in) :: center(3) !! cached center
        real(dp), intent(in) :: radius2 !! cached squared radius
        real(dp) :: d2,margin
        outside=.false.
        if (radius2>=huge(1.0_dp)*0.5_dp) return
        d2=sum((point-center)**2)
        if (.not.ieee_is_finite(d2)) return
        margin=1.0e-10_dp*max(1.0_dp,radius2,d2)
        outside=d2>radius2+margin
    end function outside_circumsphere

    ! Dimension-independent barycentric coordinates for a supplied simplex.
    ! Coordinates have SHAPE (ndim+1,ndim), irrespective of ndim: assumed-
    ! SHAPE rank-two arrays already support arbitrary geometric dimension.
    ! Gauss elimination is carried out in extended precision where available.
    pure subroutine simplex_barycentric(vertices,point,w,good)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(dp), intent(in) :: vertices(:,:) !! (dimension+1,dimension)
        real(dp), intent(in) :: point(:) !! (dimension)
        real(dp), intent(out) :: w(:) !! (dimension+1), barycentric weights
        logical, intent(out) :: good !! true if the affine simplex is nonsingular
        integer :: d,i,j,k,pivot
        real(qkind) :: matrix(size(point),size(point)), rhs(size(point)), &
                       row(size(point)), factor, cutoff, t, sol(size(point))
        d=size(point)
        w=0.0_dp
        good=.false.
        if (size(vertices,1)/=d+1 .or. size(vertices,2)/=d .or. size(w)/=d+1) return
        do j=1,d
            matrix(:,j)=real(vertices(j+1,:)-vertices(1,:),qkind)
        end do
        rhs=real(point-vertices(1,:),qkind)
        cutoff=64.0_qkind*real(epsilon(1.0_dp),qkind)*maxval(abs(matrix))
        do k=1,d
            pivot=k-1+maxloc(abs(matrix(k:d,k)),dim=1)
            if (abs(matrix(pivot,k))<=max(tiny(1.0_qkind),cutoff)) return
            if (pivot/=k) then
                row=matrix(k,:)
                matrix(k,:)=matrix(pivot,:)
                matrix(pivot,:)=row
                t=rhs(k)
                rhs(k)=rhs(pivot)
                rhs(pivot)=t
            end if
            do i=k+1,d
                factor=matrix(i,k)/matrix(k,k)
                matrix(i,k)=0.0_qkind
                matrix(i,k+1:d)=matrix(i,k+1:d)-factor*matrix(k,k+1:d)
                rhs(i)=rhs(i)-factor*rhs(k)
            end do
        end do
        do k=d,1,-1
            sol(k)=(rhs(k)-dot_product(matrix(k,k+1:d),sol(k+1:d)))/matrix(k,k)
        end do
        w(2:)=real(sol,dp)
        w(1)=real(1.0_qkind-sum(sol),dp)
        good=all(ieee_is_finite(w))
    end subroutine simplex_barycentric

    pure subroutine tetrahedron_barycentric(a,b,c,d,x,w)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(dp), intent(in) :: a(3) !! first vertex
        real(dp), intent(in) :: b(3) !! second vertex
        real(dp), intent(in) :: c(3) !! third vertex
        real(dp), intent(in) :: d(3) !! fourth vertex
        real(dp), intent(in) :: x(3) !! query point
        real(dp), intent(out) :: w(4) !! barycentric weights, NaN on singular tetrahedron
        real(qkind) :: ab(3),ac(3),ad(3),ax(3),det,wq(4)
        ab=real(b,qkind)-real(a,qkind)
        ac=real(c,qkind)-real(a,qkind)
        ad=real(d,qkind)-real(a,qkind)
        ax=real(x,qkind)-real(a,qkind)
        det=det3q(ab,ac,ad)
        if (det==0.0_qkind) then
            w=nan_real()
            return
        end if
        wq(2)=det3q(ax,ac,ad)/det
        wq(3)=det3q(ab,ax,ad)/det
        wq(4)=det3q(ab,ac,ax)/det
        wq(1)=1.0_qkind-sum(wq(2:4))
        w=real(wq,dp)
    end subroutine tetrahedron_barycentric

    ! The index is built once per interpolator.  Bounds are over actual
    ! tetrahedron vertices, so searching cannot miss a true interior point
    ! unless the input is already numerically ill-conditioned.  A small
    ! axis-specific margin admits the same near-boundary points as the
    ! barycentric tolerance used by the original linear scan.
    subroutine build_tetra_bvh(points,tetra,index)
        real(dp), intent(in) :: points(:,:)          !! (sites,3), normalized sites
        integer, intent(in) :: tetra(:,:)            !! (4,ntet), tetra vertex IDs
        type(tetra_bvh_type), intent(out) :: index    !! output spatial index
        real(dp), allocatable :: centers(:,:)
        real(dp) :: extent(3), magnitude(3)
        integer :: nt,j,root

        nt=size(tetra,2)
        if (nt<1) return
        allocate(index%order(nt),index%child_l(2*nt),index%child_r(2*nt), &
                 index%first(2*nt),index%last(2*nt), &
                 index%bbox_min(3,2*nt),index%bbox_max(3,2*nt), &
                 index%tetra_min(3,nt),index%tetra_max(3,nt))
        allocate(centers(3,nt))
        index%order=[(j,j=1,nt)]
        index%child_l=0
        index%child_r=0
        index%first=0
        index%last=0
        index%node_count=0
        do j=1,nt
            ! Divide before summing to avoid unnecessary overflow.
            centers(:,j)=sum(0.25_dp*transpose(points(tetra(:,j),:)),dim=2)
            index%tetra_min(:,j)=minval(points(tetra(:,j),:),dim=1)
            index%tetra_max(:,j)=maxval(points(tetra(:,j),:),dim=1)
        end do
        extent=maxval(points,dim=1)-minval(points,dim=1)
        magnitude=maxval(abs(points),dim=1)
        index%margin=max(4.0e-13_dp*max(1.0_dp,extent), &
                         64.0_dp*epsilon(1.0_dp)*magnitude)
        call build_bvh_node(index,centers,1,nt,root)
    end subroutine build_tetra_bvh

    ! Order centroid keys lexicographically by (coordinate, tetrahedron id).
    ! The tie break prevents Quickselect from stalling on coplanar centroids.
    pure logical function bvh_key_before(a,b,axis,centers) result(is_before)
        integer, intent(in) :: a !! first tetra ID
        integer, intent(in) :: b !! second tetra ID
        integer, intent(in) :: axis !! sort dimension
        real(dp), intent(in) :: centers(:,:) !! (3,ntet)
        is_before = centers(axis,a)<centers(axis,b)
        if (centers(axis,a)==centers(axis,b)) is_before=a<b
    end function bvh_key_before

    subroutine bvh_partition_median(order,centers,axis,lo,hi,median)
        integer, intent(inout) :: order(:)            !! index permutation
        real(dp), intent(in) :: centers(:,:)          !! tetra centroids
        integer, intent(in) :: axis                   !! partition dimension
        integer, intent(in) :: lo                     !! first index
        integer, intent(in) :: hi                     !! last index
        integer, intent(in) :: median                 !! desired final rank
        integer :: l,r,i,j,pivot,temp
        l=lo
        r=hi
        do while (l<r)
            i=l
            j=r
            pivot=order(l+(r-l)/2)
            do
                do while (bvh_key_before(order(i),pivot,axis,centers))
                    i=i+1
                end do
                do while (bvh_key_before(pivot,order(j),axis,centers))
                    j=j-1
                end do
                if (i>j) exit
                temp=order(i)
                order(i)=order(j)
                order(j)=temp
                i=i+1
                j=j-1
                if (i>j) exit
            end do
            if (median<=j) then
                r=j
            else if (median>=i) then
                l=i
            else
                exit
            end if
        end do
    end subroutine bvh_partition_median

    recursive subroutine build_bvh_node(index,centers,lo,hi,node)
        type(tetra_bvh_type), intent(inout) :: index !! index being populated
        real(dp), intent(in) :: centers(:,:)         !! tetra centroids (3,nt)
        integer, intent(in) :: lo                    !! initial permutation position
        integer, intent(in) :: hi                    !! final permutation position
        integer, intent(out) :: node                 !! allocated node number
        integer :: i,j,mid,axis,left,right
        real(dp) :: lower(3),upper(3),span(3)
        node=index%node_count+1
        index%node_count=node
        lower=huge(0.0_dp)
        upper=-huge(0.0_dp)
        do i=lo,hi
            j=index%order(i)
            lower=min(lower,index%tetra_min(:,j))
            upper=max(upper,index%tetra_max(:,j))
        end do
        index%bbox_min(:,node)=lower
        index%bbox_max(:,node)=upper
        if (hi-lo+1<=tetra_bvh_leaf_size) then
            index%first(node)=lo
            index%last(node)=hi
            return
        end if
        span=maxval(centers(:,index%order(lo:hi)),dim=2)- &
             minval(centers(:,index%order(lo:hi)),dim=2)
        axis=maxloc(span,dim=1)
        mid=lo+(hi-lo)/2
        call bvh_partition_median(index%order,centers,axis,lo,hi,mid)
        call build_bvh_node(index,centers,lo,mid,left)
        call build_bvh_node(index,centers,mid+1,hi,right)
        index%child_l(node)=left
        index%child_r(node)=right
    end subroutine build_bvh_node

    ! Return the lowest-numbered containing tetrahedron, matching the
    ! historical full scan on shared faces.  This preserves deterministic
    ! near-edge selection while typically testing only a handful of tets.
    subroutine find_containing_tetra(p,xyz,best,wbest,candidates)
        type(scattered_core), intent(in) :: p         !! valid 3D linear interpolator
        real(dp), intent(in) :: xyz(3)                !! normalized query coordinates
        integer, intent(out) :: best                 !! matching tetra ID; 0 outside
        real(dp), intent(out) :: wbest(4)             !! barycentric coordinates
        integer, intent(out), optional :: candidates !! exact barycentric tests
        integer :: stack(tetra_bvh_stack_size),top,node,k,j,a,b,c,d,visited
        real(dp) :: w(4)
        best=0
        visited=0
        wbest=0.0_dp
        top=1
        stack(1)=1
        do while (top>0)
            node=stack(top)
            top=top-1
            if (any(xyz<p%tetra_index%bbox_min(:,node)-p%tetra_index%margin)) cycle
            if (any(xyz>p%tetra_index%bbox_max(:,node)+p%tetra_index%margin)) cycle
            if (p%tetra_index%child_l(node)>0) then
                ! Median splits ensure depth < 128 for any practical mesh;
                ! order does not change the lowest-index match rule.
                if (top+2>tetra_bvh_stack_size) error stop 'tetra BVH stack exceeded'
                stack(top+1)=p%tetra_index%child_r(node)
                stack(top+2)=p%tetra_index%child_l(node)
                top=top+2
            else
                do k=p%tetra_index%first(node),p%tetra_index%last(node)
                    j=p%tetra_index%order(k)
                    if (best>0 .and. j>=best) cycle
                    if (any(xyz<p%tetra_index%tetra_min(:,j)-p%tetra_index%margin)) cycle
                    if (any(xyz>p%tetra_index%tetra_max(:,j)+p%tetra_index%margin)) cycle
                    a=p%tetrahedra(1,j)
                    b=p%tetrahedra(2,j)
                    c=p%tetrahedra(3,j)
                    d=p%tetrahedra(4,j)
                    visited=visited+1
                    call tetrahedron_barycentric(p%points(a,:),p%points(b,:), &
                        p%points(c,:),p%points(d,:),xyz,w)
                    if (.not.all(ieee_is_finite(w))) cycle
                    if (minval(w)<-2.0e-13_dp) cycle
                    if (maxval(w)>1.0_dp+2.0e-13_dp) cycle
                    best=j
                    wbest=w
                end do
            end if
        end do
        if (present(candidates)) candidates=visited
    end subroutine find_containing_tetra

    ! Save a face unless its opposite-oriented partner is already present.
    ! Faces are compared without regard to vertex order.
    subroutine add_cavity_face(faces,nfaces,f)
        integer, intent(inout) :: faces(:,:) !! 3 x capacity cavity faces
        integer, intent(inout) :: nfaces !! number of face entries, includes cancelled ones
        integer, intent(in) :: f(3) !! oriented boundary-face vertices
        integer :: j,k
        do j=1,nfaces
            if (faces(1,j)==0) cycle
            if (all([(any(faces(:,j)==f(k)),k=1,3)])) then
                faces(:,j)=0
                return
            end if
        end do
        nfaces=nfaces+1
        if (nfaces>size(faces,2)) error stop 'delaunay_3d cavity capacity exceeded'
        faces(:,nfaces)=f
    end subroutine add_cavity_face

    ! Store the current tetrahedral hull facets in an open-addressed hash
    ! table, then add the faces visible from the new point to the cavity.
    ! The 3D hull is never represented by artificial/super-tetra vertices:
    ! that avoids losing thin hull tetrahedra upon supertet removal.
    subroutine append_visible_hull_faces(points,tetra,nt,p,faces,nfaces)
        real(dp), intent(in) :: points(:,:) !! (sites,3), uniformly scaled
        integer, intent(in) :: tetra(:,:) !! (4,capacity), positively oriented
        integer, intent(in) :: nt !! number of current tetrahedra
        real(dp), intent(in) :: p(3) !! point about to be inserted
        integer, intent(inout) :: faces(:,:) !! cavity face list
        integer, intent(inout) :: nfaces !! occupied or canceled cavity slots
        integer, allocatable :: keys(:,:),outward(:,:),count(:)
        integer :: cap,j,k,h,u(3),f(3),v(4),tmp
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        cap=8*nt+11
        allocate(keys(3,cap),outward(3,cap),count(cap))
        keys=0
        count=0
        do j=1,nt
            v=tetra(:,j)
            do k=1,4
                select case(k)
                case(1)
                    f=[v(2),v(3),v(4)]
                case(2)
                    f=[v(1),v(4),v(3)]
                case(3)
                    f=[v(1),v(2),v(4)]
                case default
                    f=[v(1),v(3),v(2)]
                end select
                u=f
                if (u(1)>u(2)) then
                    tmp=u(1);u(1)=u(2);u(2)=tmp
                end if
                if (u(2)>u(3)) then
                    tmp=u(2);u(2)=u(3);u(3)=tmp
                end if
                if (u(1)>u(2)) then
                    tmp=u(1);u(1)=u(2);u(2)=tmp
                end if
                h=1+int(modulo(int(u(1),int64)*73856093_int64 &
                           +int(u(2),int64)*19349663_int64 &
                           +int(u(3),int64)*83492791_int64,int(cap,int64)))
                do
                    if (count(h)==0) then
                        keys(:,h)=u
                        outward(:,h)=f
                        count(h)=1
                        exit
                    else if (all(keys(:,h)==u)) then
                        count(h)=count(h)+1
                        exit
                    end if
                    h=mod(h,cap)+1
                end do
            end do
        end do
        do h=1,cap
            if (count(h)/=1) cycle
            f=outward(:,h)
            if (orient3q(points(f(1),:),points(f(2),:),points(f(3),:),p)>0.0_qkind) then
                call add_cavity_face(faces,nfaces,f)
            end if
        end do
    end subroutine append_visible_hull_faces

    subroutine delaunay_3d(x,simplices,code)
        real(dp), intent(in) :: x(:,:) !! unique finite sites, (n,3)
        integer, allocatable, intent(out) :: simplices(:,:) !! 4 x nt, 1-based site indices
        integer, intent(out) :: code !! success, invalid input or degenerate geometry
        integer :: n,i,j,nt,max_t,kept,nfaces,a,b,c,d,idx(4),seed(4),nt_final
        integer, allocatable :: tetra(:,:),faces(:,:)
        logical, allocatable :: bad(:),seen(:)
        real(dp), allocatable :: sphere_center(:,:),sphere_radius2(:)
        real(dp), allocatable :: points(:,:)
        real(dp) :: mid(3),span,measure,best,ab(3),ac(3),normal(3)
        integer, parameter :: wanted_kind=selected_real_kind(33,4931)
        integer, parameter :: qkind=merge(wanted_kind,dp,wanted_kind>0)
        real(qkind) :: volume

        code=scattered_status_invalid_input
        n=size(x,1)
        if (size(x,2)/=3 .or. n<4) return
        if (.not.all(ieee_is_finite(x))) return
        code=scattered_status_degenerate
        mid=0.5_dp*minval(x,dim=1)+0.5_dp*maxval(x,dim=1)
        span=maxval(maxval(x,dim=1)-minval(x,dim=1))
        if (span<=0.0_dp .or. .not.ieee_is_finite(span)) return
        allocate(points(n,3))
        points=(x-spread(mid,1,n))/span

        ! Choose a well-conditioned initial tetrahedron: most distant point,
        ! greatest line distance and greatest plane distance, respectively.
        seed(1)=1
        seed(2)=0
        best=0.0_dp
        do i=2,n
            measure=sum((points(i,:)-points(seed(1),:))**2)
            if (measure>best) then
                best=measure
                seed(2)=i
            end if
        end do
        if (seed(2)==0) return
        ab=points(seed(2),:)-points(seed(1),:)
        seed(3)=0
        best=0.0_dp
        do i=1,n
            if (any(seed(:2)==i)) cycle
            ac=points(i,:)-points(seed(1),:)
            normal=cross3(ab,ac)
            measure=sum(normal**2)
            if (measure>best) then
                best=measure
                seed(3)=i
            end if
        end do
        if (seed(3)==0) return
        ac=points(seed(3),:)-points(seed(1),:)
        normal=cross3(ab,ac)
        seed(4)=0
        best=0.0_dp
        do i=1,n
            if (any(seed(:3)==i)) cycle
            measure=abs(dot_product(normal,points(i,:)-points(seed(1),:)))
            if (measure>best) then
                best=measure
                seed(4)=i
            end if
        end do
        if (seed(4)==0 .or. best<=1.0e-14_dp) return
        if (orient3q(points(seed(1),:),points(seed(2),:), &
            points(seed(3),:),points(seed(4),:))<0.0_qkind) then
            seed=[seed(1),seed(3),seed(2),seed(4)]
        end if
        max_t=24*n+128
        allocate(tetra(4,max_t),bad(max_t),faces(3,4*max_t),seen(n), &
                 sphere_center(3,max_t),sphere_radius2(max_t))
        tetra(:,1)=seed
        call tetra_circumsphere(points(seed(1),:),points(seed(2),:), &
             points(seed(3),:),points(seed(4),:),sphere_center(:,1),sphere_radius2(1))
        seen=.false.
        seen(seed)=.true.
        nt=1
        do i=1,n
            if (seen(i)) cycle
            bad(:nt)=.false.
            nfaces=0
            do j=1,nt
                idx=tetra(:,j)
                if (outside_circumsphere(points(i,:),sphere_center(:,j),sphere_radius2(j))) cycle
                if (.not.insphere3d(points(idx(1),:),points(idx(2),:), &
                     points(idx(3),:),points(idx(4),:),points(i,:))) cycle
                bad(j)=.true.
                a=idx(1);b=idx(2);c=idx(3);d=idx(4)
                call add_cavity_face(faces,nfaces,[b,c,d])
                call add_cavity_face(faces,nfaces,[a,d,c])
                call add_cavity_face(faces,nfaces,[a,b,d])
                call add_cavity_face(faces,nfaces,[a,c,b])
            end do
            call append_visible_hull_faces(points,tetra,nt,points(i,:),faces,nfaces)
            if (.not.any(bad(:nt)) .and. all(faces(:,1:nfaces)==0)) then
                ! A cospherical point might lie inside a tetrahedron without
                ! strictly lying in any computed circumsphere.
                do j=1,nt
                    idx=tetra(:,j)
                    if (point_in_tetra3d(points(i,:),points(idx(1),:), &
                        points(idx(2),:),points(idx(3),:),points(idx(4),:))) exit
                end do
                if (j>nt) return
                bad(j)=.true.
                a=idx(1);b=idx(2);c=idx(3);d=idx(4)
                call add_cavity_face(faces,nfaces,[b,c,d])
                call add_cavity_face(faces,nfaces,[a,d,c])
                call add_cavity_face(faces,nfaces,[a,b,d])
                call add_cavity_face(faces,nfaces,[a,c,b])
            end if
            kept=0
            do j=1,nt
                if (bad(j)) cycle
                kept=kept+1
                tetra(:,kept)=tetra(:,j)
                sphere_center(:,kept)=sphere_center(:,j)
                sphere_radius2(kept)=sphere_radius2(j)
            end do
            nt=kept
            do j=1,nfaces
                if (faces(1,j)==0) cycle
                idx(:3)=faces(:,j)
                idx(4)=i
                volume=orient3q(points(idx(1),:),points(idx(2),:), &
                    points(idx(3),:),points(idx(4),:))
                if (volume==0.0_qkind) cycle
                if (volume<0.0_qkind) idx=[idx(1),idx(3),idx(2),idx(4)]
                nt=nt+1
                if (nt>max_t) return
                tetra(:,nt)=idx
                call tetra_circumsphere(points(idx(1),:),points(idx(2),:), &
                    points(idx(3),:),points(idx(4),:),sphere_center(:,nt),sphere_radius2(nt))
            end do
            if (nt==0) return
            seen(i)=.true.
        end do
        nt_final=nt
        allocate(simplices(4,nt_final))
        simplices=tetra(:,1:nt_final)
        code=scattered_status_success
    end subroutine delaunay_3d

    pure function cross3(a,b) result(c)
        real(dp), intent(in) :: a(3) !! vector A
        real(dp), intent(in) :: b(3) !! vector B
        real(dp) :: c(3)
        c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    end function cross3

    pure logical function point_in_tetra3d(p,a,b,c,d) result(inside)
        real(dp), intent(in) :: p(3) !! candidate point
        real(dp), intent(in) :: a(3) !! first vertex
        real(dp), intent(in) :: b(3) !! second vertex
        real(dp), intent(in) :: c(3) !! third vertex
        real(dp), intent(in) :: d(3) !! fourth vertex
        real(dp) :: w(4)
        call tetrahedron_barycentric(a,b,c,d,p,w)
        inside=all(w>=-1.0e-13_dp)
    end function point_in_tetra3d

    pure real(dp) function cross2(a,b,c) result(area)
        real(dp), intent(in) :: a(2),b(2),c(2) !! triangle vertices
        area=(b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1))
    end function cross2

    ! Use extended precision for the determinant of the in-circle test.
    ! Double precision cancellation around the enclosing supertriangle can
    ! otherwise leave thin wedges of the true convex hull untriangulated.
    pure logical function incircle(a,b,c,p) result(inside)
        real(dp), intent(in) :: a(2),b(2),c(2),p(2) !! triangle and inserted point
        integer, parameter :: wide_kind=selected_real_kind(33,4931)
        ! Fall back to double precision on compilers lacking IEEE quad.
        integer, parameter :: qkind=merge(wide_kind,dp,wide_kind>0)
        real(qkind) :: ax,ay,bx,by,cx,cy,aa,bb,cc,det,orient
        ax=real(a(1),qkind)-real(p(1),qkind)
        ay=real(a(2),qkind)-real(p(2),qkind)
        bx=real(b(1),qkind)-real(p(1),qkind)
        by=real(b(2),qkind)-real(p(2),qkind)
        cx=real(c(1),qkind)-real(p(1),qkind)
        cy=real(c(2),qkind)-real(p(2),qkind)
        aa=ax*ax+ay*ay
        bb=bx*bx+by*by
        cc=cx*cx+cy*cy
        det=aa*(bx*cy-by*cx)-bb*(ax*cy-ay*cx)+cc*(ax*by-ay*bx)
        orient=(real(b(1),qkind)-real(a(1),qkind)) * &
               (real(c(2),qkind)-real(a(2),qkind)) - &
               (real(b(2),qkind)-real(a(2),qkind)) * &
               (real(c(1),qkind)-real(a(1),qkind))
        inside=(det>0.0_qkind .and. orient>0.0_qkind) .or. &
               (det<0.0_qkind .and. orient<0.0_qkind)
    end function incircle

    subroutine delaunay_2d(x,simplices,code)
        real(dp), intent(in) :: x(:,:) !! (sites,2) or (sites,3), unique finite coordinates
        integer, allocatable, intent(out) :: simplices(:,:) !! resulting oriented triangles, 3 x nt
        integer, intent(out) :: code !! success or degenerate
        real(dp), allocatable :: q(:,:)
        integer, allocatable :: t(:,:),edges(:,:)
        logical, allocatable :: bad(:)
        integer :: n,nt,max_t,i,j,k,l,ne,kept,a,b,nt_final
        real(dp) :: cx,cy,span,area,mesh_area,hull_area
        n=size(x,1)
        code=scattered_status_degenerate
        allocate(q(n+3,2))
        q(:n,:)=x
        span=maxval(maxval(x,dim=1)-minval(x,dim=1))
        if (span<=0.0_dp .or. .not.ieee_is_finite(span)) return
        cx=0.5_dp*(minval(x(:,1))+maxval(x(:,1)))
        cy=0.5_dp*(minval(x(:,2))+maxval(x(:,2)))
        q(n+1,:)=[cx-16.0_dp*span,cy-8.0_dp*span]
        q(n+2,:)=[cx+16.0_dp*span,cy-8.0_dp*span]
        q(n+3,:)=[cx,cy+16.0_dp*span]
        max_t=8*n+32
        allocate(t(3,max_t),bad(max_t),edges(2,3*max_t))
        t(:,1)=[n+1,n+2,n+3]
        nt=1
        do i=1,n
            bad(:nt)=.false.
            ne=0
            do j=1,nt
                if (.not.incircle(q(t(1,j),:),q(t(2,j),:),q(t(3,j),:),q(i,:))) cycle
                bad(j)=.true.
                do k=1,3
                    a=t(k,j)
                    b=t(mod(k,3)+1,j)
                    ! Each internal cavity edge occurs twice with opposite
                    ! orientation. Mark paired edges and keep the boundary.
                    do l=1,ne
                        if (edges(1,l)==b .and. edges(2,l)==a) exit
                    end do
                    if (l<=ne) then
                        edges(:,l)=0
                    else
                        ne=ne+1
                        if (ne>size(edges,2)) return
                        edges(:,ne)=[a,b]
                    end if
                end do
            end do
            ! A point can lie on circumcircle boundaries of a perfectly
            ! cocircular configuration. Find the containing triangle instead.
            if (.not.any(bad(:nt))) then
                do j=1,nt
                    if (point_in_triangle(q(i,:),q(t(1,j),:),q(t(2,j),:),q(t(3,j),:))) exit
                end do
                if (j>nt) return
                bad(j)=.true.
                ne=3
                edges(:,1:3)=reshape([t(1,j),t(2,j),t(2,j),t(3,j),t(3,j),t(1,j)],[2,3])
            end if
            kept=0
            do j=1,nt
                if (bad(j)) cycle
                kept=kept+1
                t(:,kept)=t(:,j)
            end do
            nt=kept
            do k=1,ne
                if (all(edges(:,k)==0)) cycle
                a=edges(1,k)
                b=edges(2,k)
                area=cross2(q(a,:),q(b,:),q(i,:))
                if (area==0.0_dp) cycle
                nt=nt+1
                if (nt>max_t) return
                if (area>0.0_dp) then
                    t(:,nt)=[a,b,i]
                else
                    t(:,nt)=[b,a,i]
                end if
            end do
        end do
        nt_final=0
        do j=1,nt
            if (any(t(:,j)>n)) cycle
            if (cross2(q(t(1,j),:),q(t(2,j),:),q(t(3,j),:))<=0.0_dp) cycle
            nt_final=nt_final+1
            t(:,nt_final)=t(:,j)
        end do
        if (nt_final==0) return
        ! Patch boundary notches left by supertriangle circumcircle rounding.
        ! Ear clipping only fills missing regions between a true hull edge
        ! and the corresponding existing exterior boundary chain.
        call repair_hull_notches(x,t,nt_final,max_t)
        call legalize_delaunay_edges(x,t,nt_final)
        ! The supertriangle and finite-precision circumcircle predicates can
        ! occasionally leave a small gap near a very elongated convex hull.
        ! Detect it rather than silently interpolating over an incomplete mesh.
        mesh_area=0.0_dp
        do j=1,nt_final
            mesh_area=mesh_area+0.5_dp*abs(cross2(x(t(1,j),:),x(t(2,j),:),x(t(3,j),:)))
        end do
        hull_area=convex_hull_area_2d(x)
        if (hull_area<=0.0_dp) return
        if (abs(mesh_area-hull_area)>1.0e-9_dp*hull_area) return
        allocate(simplices(3,nt_final))
        simplices=t(:,1:nt_final)
        code=scattered_status_success
    end subroutine delaunay_2d

    subroutine legalize_delaunay_edges(x,tri,nt)
        real(dp), intent(in) :: x(:,:) !! finite sites (sites,2)
        integer, intent(inout) :: tri(:,:) !! (3,max_t) CCW oriented triangles
        integer, intent(in) :: nt !! active triangle count
        integer :: pass,i,j,k,common(2),nc,a,b,c,d,temp(3)
        real(dp) :: side1,side2
        logical :: flipped
        ! Lawson edge flips restore the local Delaunay condition after hull
        ! completion (and correct rare cavity updates on slender triangles).
        do pass=1,max(1,8*nt)
            flipped=.false.
            do i=1,nt-1
                do j=i+1,nt
                    nc=0
                    c=0
                    do k=1,3
                        if (any(tri(:,j)==tri(k,i))) then
                            nc=nc+1
                            if (nc<=2) common(nc)=tri(k,i)
                        else
                            c=tri(k,i)
                        end if
                    end do
                    if (nc/=2) cycle
                    d=0
                    do k=1,3
                        if (.not.any(tri(:,i)==tri(k,j))) d=tri(k,j)
                    end do
                    if (c==0 .or. d==0) cycle
                    a=common(1)
                    b=common(2)
                    side1=cross2(x(a,:),x(b,:),x(c,:))
                    side2=cross2(x(a,:),x(b,:),x(d,:))
                    if (side1*side2>=0.0_dp) cycle
                    side1=cross2(x(c,:),x(d,:),x(a,:))
                    side2=cross2(x(c,:),x(d,:),x(b,:))
                    if (side1*side2>=0.0_dp) cycle
                    if (.not.incircle(x(a,:),x(b,:),x(c,:),x(d,:))) cycle
                    temp=[c,d,a]
                    if (cross2(x(temp(1),:),x(temp(2),:),x(temp(3),:))<0.0_dp) &
                        temp=[d,c,a]
                    tri(:,i)=temp
                    temp=[d,c,b]
                    if (cross2(x(temp(1),:),x(temp(2),:),x(temp(3),:))<0.0_dp) &
                        temp=[c,d,b]
                    tri(:,j)=temp
                    flipped=.true.
                    exit
                end do
                if (flipped) exit
            end do
            if (.not.flipped) return
        end do
    end subroutine legalize_delaunay_edges

    pure integer function triangle_edge_count(tri,a,b,ntri) result(nedge)
        integer, intent(in) :: tri(:,:) !! triangle indices (3,max_t)
        integer, intent(in) :: a !! first endpoint
        integer, intent(in) :: b !! second endpoint
        integer, intent(in) :: ntri !! active triangle count
        integer :: j
        nedge=0
        do j=1,ntri
            if (any(tri(:,j)==a) .and. any(tri(:,j)==b)) nedge=nedge+1
        end do
    end function triangle_edge_count

    pure function convex_hull_indices_2d(x) result(hull)
        real(dp), intent(in) :: x(:,:) !! (sites,2) or (sites,3) planar data
        integer, allocatable :: hull(:)
        integer, allocatable :: path(:)
        integer :: n,first,p,q,j,step,nh
        real(dp) :: cross,dist1,dist2
        n=size(x,1)
        allocate(path(n+1))
        if (n<3) then
            allocate(hull(0))
            return
        end if
        first=1
        do j=2,n
            if (x(j,1)<x(first,1)) then
                first=j
            else if (x(j,1)==x(first,1) .and. x(j,2)<x(first,2)) then
                first=j
            end if
        end do
        p=first
        nh=0
        do step=1,n+1
            nh=nh+1
            path(nh)=p
            q=1
            if (q==p) q=2
            do j=1,n
                if (j==p .or. j==q) cycle
                cross=cross2(x(p,:),x(q,:),x(j,:))
                if (cross>0.0_dp) then
                    q=j
                else if (cross==0.0_dp) then
                    dist1=sum((x(j,:)-x(p,:))**2)
                    dist2=sum((x(q,:)-x(p,:))**2)
                    if (dist1>dist2) q=j
                end if
            end do
            p=q
            if (p==first) exit
        end do
        hull=path(:nh)
    end function convex_hull_indices_2d

    subroutine repair_hull_notches(x,tri,nt,max_t)
        real(dp), intent(in) :: x(:,:) !! (n,2) source points
        integer, intent(inout) :: tri(:,:) !! triangles with scratch capacity
        integer, intent(inout) :: nt !! number of completed triangles
        integer, intent(in) :: max_t !! capacity of the triangle array
        integer, allocatable :: hull(:),path(:)
        integer :: n,k,a,b,candidate,j,prev,cur,next_vertex,npath,step
        logical :: reached
        hull=convex_hull_indices_2d(x)
        n=size(x,1)
        allocate(path(n+1))
        do k=1,size(hull)
            a=hull(k)
            b=hull(mod(k,size(hull))+1)
            if (triangle_edge_count(tri,a,b,nt)>0) cycle
            do candidate=1,n
                if (candidate==a .or. candidate==b) cycle
                if (triangle_edge_count(tri,a,candidate,nt)/=1) cycle
                if (any(hull==candidate)) cycle
                path(1)=a
                path(2)=candidate
                npath=2
                prev=a
                cur=candidate
                reached=.false.
                do step=3,n+1
                    if (cur==b) then
                        reached=.true.
                        exit
                    end if
                    next_vertex=0
                    do j=1,n
                        if (j==prev .or. j==cur) cycle
                        if (triangle_edge_count(tri,cur,j,nt)==1) then
                            next_vertex=j
                            exit
                        end if
                    end do
                    if (next_vertex==0 .or. next_vertex==a) exit
                    if (next_vertex/=b .and. any(hull==next_vertex)) exit
                    npath=npath+1
                    path(npath)=next_vertex
                    prev=cur
                    cur=next_vertex
                    if (cur==b) then
                        reached=.true.
                        exit
                    end if
                end do
                if (.not.reached) cycle
                call triangulate_boundary_notch(x,path(:npath),tri,nt,max_t)
                exit
            end do
        end do
    end subroutine repair_hull_notches

    subroutine triangulate_boundary_notch(x,path,tri,nt,max_t)
        real(dp), intent(in) :: x(:,:) !! (n,2) coordinates
        integer, intent(in) :: path(:) !! polygon along concave native boundary
        integer, intent(inout) :: tri(:,:) !! CCW triangles with scratch capacity
        integer, intent(inout) :: nt !! active triangle count
        integer, intent(in) :: max_t !! capacity
        integer, allocatable :: polygon(:)
        integer :: m,j,prev,next,k,a,b,c
        real(dp) :: area,orientation
        logical :: found,contains
        polygon=path
        m=size(polygon)
        area=0.0_dp
        do j=1,m
            a=polygon(j)
            b=polygon(mod(j,m)+1)
            area=area+x(a,1)*x(b,2)-x(a,2)*x(b,1)
        end do
        if (abs(area)<=tiny(area) .or. nt+m-2>max_t) return
        orientation=sign(1.0_dp,area)
        do while(m>=3)
            found=.false.
            do j=1,m
                prev=mod(j+m-2,m)+1
                next=mod(j,m)+1
                a=polygon(prev)
                b=polygon(j)
                c=polygon(next)
                area=cross2(x(a,:),x(b,:),x(c,:))
                if (orientation*area<=tiny(area)) cycle
                contains=.false.
                do k=1,m
                    if (k==prev .or. k==j .or. k==next) cycle
                    if (point_in_triangle(x(polygon(k),:),x(a,:),x(b,:),x(c,:))) then
                        contains=.true.
                        exit
                    end if
                end do
                if (contains) cycle
                nt=nt+1
                if (area>0.0_dp) then
                    tri(:,nt)=[a,b,c]
                else
                    tri(:,nt)=[c,b,a]
                end if
                if (m==3) return
                polygon(j:m-1)=polygon(j+1:m)
                m=m-1
                found=.true.
                exit
            end do
            if (.not.found) return
        end do
    end subroutine triangulate_boundary_notch

    ! Jarvis march determines the enclosing convex hull independently from
    ! the triangulation, for detecting missing or overlapping triangles.
    pure real(dp) function convex_hull_area_2d(x) result(area)
        real(dp), intent(in) :: x(:,:) !! (sites,2) or (sites,3) planar points
        integer :: n,first,p,q,j,step
        real(dp) :: cross,dist1,dist2
        n=size(x,1)
        area=0.0_dp
        if (n<3) return
        first=1
        do j=2,n
            if (x(j,1)<x(first,1)) then
                first=j
            else if (x(j,1)==x(first,1) .and. x(j,2)<x(first,2)) then
                first=j
            end if
        end do
        p=first
        do step=1,n+1
            q=1
            if (q==p) q=2
            do j=1,n
                if (j==p .or. j==q) cycle
                cross=cross2(x(p,:),x(q,:),x(j,:))
                if (cross>0.0_dp) then
                    q=j
                else if (cross==0.0_dp) then
                    dist1=sum((x(j,:)-x(p,:))**2)
                    dist2=sum((x(q,:)-x(p,:))**2)
                    if (dist1>dist2) q=j
                end if
            end do
            area=area+x(p,1)*x(q,2)-x(q,1)*x(p,2)
            p=q
            if (p==first) exit
        end do
        area=abs(area)*0.5_dp
    end function convex_hull_area_2d

    pure logical function point_in_triangle(x,a,b,c) result(inside)
        real(dp), intent(in) :: x(2),a(2),b(2),c(2) !! query and vertices
        real(dp) :: det,wa,wb,wc
        inside=.false.
        det=cross2(a,b,c)
        if (det==0.0_dp) return
        wa=cross2(b,c,x)/det
        wb=cross2(c,a,x)/det
        wc=1.0_dp-wa-wb
        inside=min(wa,wb,wc)>=-1.0e-12_dp
    end function point_in_triangle

    function core_evaluate(p,x,code) result(v)
        type(scattered_core), intent(in) :: p !! interpolator
        real(dp), intent(in) :: x(:) !! evaluation position
        integer, intent(out) :: code !! status for malformed query
        complex(dp), allocatable :: v(:)
        integer :: j,best,a,b,c,found_tetra,k
        real(dp) :: xy(2),xyz(3),det,w1,w2,w3,w(4)
        real(dp), allocatable :: query(:), weights(:)
        logical :: nonsingular
        code=scattered_status_invalid_input
        allocate(v(0))
        if (.not.core_valid(p)) return
        deallocate(v)
        allocate(v(size(p%values,2)))
        v=cmplx(nan_real(),0.0_dp,dp)
        if (size(x)/=size(p%points,2)) return
        if (.not.all(ieee_is_finite(x))) return
        code=scattered_status_success
        if (.not.p%linear) then
            best=p%spatial_index%nearest((x-p%offset)/p%scale, status=code)
            if (code==kdtree_success) v=p%values(best,:)
            return
        end if
        v=p%fill
        if (allocated(p%external_simplices)) then
            query=(x-p%offset)/p%scale
            allocate(weights(size(query)+1))
            call find_containing_simplex(p,query,j,weights)
            if (j>0) then
                v=cmplx(0.0_dp,0.0_dp,dp)
                do k=1,size(weights)
                    v=v+weights(k)*p%values(p%external_simplices(k,j),:)
                end do
            end if
            return
        end if
        if (allocated(p%tetrahedra)) then
            xyz=(x-p%offset)/p%scale
            call find_containing_tetra(p,xyz,found_tetra,w)
            if (found_tetra>0) then
                a=p%tetrahedra(1,found_tetra)
                b=p%tetrahedra(2,found_tetra)
                c=p%tetrahedra(3,found_tetra)
                v=w(1)*p%values(a,:)+w(2)*p%values(b,:)+ &
                    w(3)*p%values(c,:)+w(4)*p%values(p%tetrahedra(4,found_tetra),:)
            end if
            return
        end if
        xy=(x-p%offset)/p%scale
        do j=1,size(p%triangles,2)
            a=p%triangles(1,j)
            b=p%triangles(2,j)
            c=p%triangles(3,j)
            det=cross2(p%points(a,:),p%points(b,:),p%points(c,:))
            if (det<=0.0_dp) cycle
            w1=cross2(p%points(b,:),p%points(c,:),xy)/det
            w2=cross2(p%points(c,:),p%points(a,:),xy)/det
            w3=1.0_dp-w1-w2
            if (min(w1,w2,w3)<-2.0e-13_dp) cycle
            v=w1*p%values(a,:)+w2*p%values(b,:)+w3*p%values(c,:)
            return
        end do
    end function core_evaluate

    function core_evaluate_many(p,x,code) result(v)
        type(scattered_core), intent(in) :: p !! interpolator
        real(dp), intent(in) :: x(:,:) !! (query points,dimensions)
        integer, intent(out) :: code !! status
        complex(dp), allocatable :: v(:,:),temp(:)
        integer :: i,stat
        code=scattered_status_invalid_input
        allocate(v(0,0))
        if (.not.core_valid(p)) return
        if (size(x,2)/=size(p%points,2)) return
        deallocate(v)
        allocate(v(size(x,1),size(p%values,2)))
        code=scattered_status_success
        do i=1,size(x,1)
            temp=core_evaluate(p,x(i,:),stat)
            v(i,:)=temp
            if (stat/=scattered_status_success) code=stat
        end do
    end function core_evaluate_many

    function nearest_real_one(self,x,status) result(v)
        class(nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! coordinates
        integer, intent(out), optional :: status !! success/invalid
        real(dp), allocatable :: v(:)
        complex(dp), allocatable :: temp(:)
        integer :: code
        temp=core_evaluate(self%core,x,code)
        v=real(temp,dp)
        if (present(status)) status=code
    end function nearest_real_one

    function nearest_real_many(self,x,status) result(v)
        class(nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:,:) !! (queries,dimensions)
        integer, intent(out), optional :: status !! success/invalid
        real(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: temp(:,:)
        integer :: code
        temp=core_evaluate_many(self%core,x,code)
        v=real(temp,dp)
        if (present(status)) status=code
    end function nearest_real_many

    function nearest_complex_one(self,x,status) result(v)
        class(complex_nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! coordinates
        integer, intent(out), optional :: status !! success/invalid
        complex(dp), allocatable :: v(:)
        integer :: code
        v=core_evaluate(self%core,x,code)
        if (present(status)) status=code
    end function nearest_complex_one

    function nearest_complex_many(self,x,status) result(v)
        class(complex_nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:,:) !! (queries,dimensions)
        integer, intent(out), optional :: status !! success/invalid
        complex(dp), allocatable :: v(:,:)
        integer :: code
        v=core_evaluate_many(self%core,x,code)
        if (present(status)) status=code
    end function nearest_complex_many

    function linear_real_one(self,x,status) result(v)
        class(linear_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! (ndim) coordinates
        integer, intent(out), optional :: status !! success/invalid
        real(dp), allocatable :: v(:)
        complex(dp), allocatable :: temp(:)
        integer :: code
        temp=core_evaluate(self%core,x,code)
        v=real(temp,dp)
        if (present(status)) status=code
    end function linear_real_one

    function linear_real_many(self,x,status) result(v)
        class(linear_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:,:) !! (queries,ndim)
        integer, intent(out), optional :: status !! success/invalid
        real(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: temp(:,:)
        integer :: code
        temp=core_evaluate_many(self%core,x,code)
        v=real(temp,dp)
        if (present(status)) status=code
    end function linear_real_many

    function linear_complex_one(self,x,status) result(v)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! (ndim) coordinates
        integer, intent(out), optional :: status !! success/invalid
        complex(dp), allocatable :: v(:)
        integer :: code
        v=core_evaluate(self%core,x,code)
        if (present(status)) status=code
    end function linear_complex_one

    function linear_complex_many(self,x,status) result(v)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:,:) !! (queries,ndim)
        integer, intent(out), optional :: status !! success/invalid
        complex(dp), allocatable :: v(:,:)
        integer :: code
        v=core_evaluate_many(self%core,x,code)
        if (present(status)) status=code
    end function linear_complex_many

    function nearest_real_value(self,x,channel,status) result(v)
        class(nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! coordinate vector
        integer, intent(in), optional :: channel !! output channel; default 1
        integer, intent(out), optional :: status !! success/invalid
        real(dp) :: v
        real(dp), allocatable :: tmp(:)
        integer :: code,c
        tmp=self%evaluate(x,status=code)
        c=1
        if (present(channel)) c=channel
        v=nan_real()
        if (c<1 .or. c>size(tmp)) code=scattered_status_invalid_input
        if (code==0) v=tmp(c)
        if (present(status)) status=code
    end function nearest_real_value

    function nearest_complex_value(self,x,channel,status) result(v)
        class(complex_nearest_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! coordinate vector
        integer, intent(in), optional :: channel !! output channel; default 1
        integer, intent(out), optional :: status !! success/invalid
        complex(dp) :: v
        complex(dp), allocatable :: tmp(:)
        integer :: code,c
        tmp=self%evaluate(x,status=code)
        c=1
        if (present(channel)) c=channel
        v=cmplx(nan_real(),0.0_dp,dp)
        if (c<1 .or. c>size(tmp)) code=scattered_status_invalid_input
        if (code==0) v=tmp(c)
        if (present(status)) status=code
    end function nearest_complex_value

    function linear_real_value(self,x,channel,status) result(v)
        class(linear_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! coordinates
        integer, intent(in), optional :: channel !! output channel; default 1
        integer, intent(out), optional :: status !! success/invalid
        real(dp) :: v
        real(dp), allocatable :: tmp(:)
        integer :: code,c
        tmp=self%evaluate(x,status=code)
        c=1
        if (present(channel)) c=channel
        v=nan_real()
        if (c<1 .or. c>size(tmp)) code=scattered_status_invalid_input
        if (code==0) v=tmp(c)
        if (present(status)) status=code
    end function linear_real_value

    function linear_complex_value(self,x,channel,status) result(v)
        class(complex_linear_nd_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! coordinates
        integer, intent(in), optional :: channel !! output channel; default 1
        integer, intent(out), optional :: status !! success/invalid
        complex(dp) :: v
        complex(dp), allocatable :: tmp(:)
        integer :: code,c
        tmp=self%evaluate(x,status=code)
        c=1
        if (present(channel)) c=channel
        v=cmplx(nan_real(),0.0_dp,dp)
        if (c<1 .or. c>size(tmp)) code=scattered_status_invalid_input
        if (code==0) v=tmp(c)
        if (present(status)) status=code
    end function linear_complex_value

    ! SciPy-style convenience interface: nearest for any dimension and native
    ! automatic Delaunay linear interpolation for 2D-5D. An explicit
    ! simplices= array enables higher-dimensional LinearND constructors.
    function griddata_real_vector(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        real(dp), intent(in) :: y(:) !! site values
        real(dp), intent(in) :: xi(:,:) !! (queries,dimensions)
        character(len=*), intent(in), optional :: method !! linear (default) or nearest
        real(dp), intent(in), optional :: fill_value !! linear outside-hull value
        logical, intent(in), optional :: rescale !! rescale to unit cube
        integer, intent(out), optional :: status !! interpolation status
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: a(:,:)
        a=griddata_real_matrix(x,reshape(y,[size(y),1]),xi,method,fill_value,rescale,status)
        if (size(a,2)==0) then
            allocate(v(0))
        else
            v=a(:,1)
        end if
    end function griddata_real_vector

    function griddata_real_matrix(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        real(dp), intent(in) :: y(:,:) !! (sites,channels)
        real(dp), intent(in) :: xi(:,:) !! (queries,dimensions)
        character(len=*), intent(in), optional :: method !! linear (default), nearest
        real(dp), intent(in), optional :: fill_value !! value beyond hull
        logical, intent(in), optional :: rescale !! unit-cube scaling
        integer, intent(out), optional :: status !! status code
        real(dp), allocatable :: v(:,:)
        type(linear_nd_interpolator_type) :: line
        type(nearest_nd_interpolator_type) :: near
        character(len=:), allocatable :: name
        integer :: code
        name='linear'
        if (present(method)) name=trim(method)
        code=scattered_status_invalid_input
        allocate(v(0,0))
        select case(name)
        case('linear')
            line=linear_nd_interpolator(x,y,rescale=rescale,fill_value=fill_value,status=code)
            if (code==0) v=line%evaluate(xi,status=code)
        case('nearest')
            near=nearest_nd_interpolator(x,y,rescale=rescale,status=code)
            if (code==0) v=near%evaluate(xi,status=code)
        end select
        if (present(status)) status=code
    end function griddata_real_matrix

    function griddata_complex_vector(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        complex(dp), intent(in) :: y(:) !! complex site values
        real(dp), intent(in) :: xi(:,:) !! query positions
        character(len=*), intent(in), optional :: method !! linear (default), nearest
        complex(dp), intent(in), optional :: fill_value !! linear outside-hull value
        logical, intent(in), optional :: rescale !! normalize coordinate axes
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: v(:)
        complex(dp), allocatable :: a(:,:)
        a=griddata_complex_matrix(x,reshape(y,[size(y),1]),xi,method,fill_value,rescale,status)
        if (size(a,2)==0) then
            allocate(v(0))
        else
            v=a(:,1)
        end if
    end function griddata_complex_vector

    function griddata_complex_matrix(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! (sites,dimensions)
        complex(dp), intent(in) :: y(:,:) !! (sites,channels)
        real(dp), intent(in) :: xi(:,:) !! query positions
        character(len=*), intent(in), optional :: method !! linear (default), nearest
        complex(dp), intent(in), optional :: fill_value !! linear outside-hull value
        logical, intent(in), optional :: rescale !! rescale to unit cube
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: v(:,:)
        type(complex_linear_nd_interpolator_type) :: line
        type(complex_nearest_nd_interpolator_type) :: near
        character(len=:), allocatable :: name
        integer :: code
        name='linear'
        if (present(method)) name=trim(method)
        code=scattered_status_invalid_input
        allocate(v(0,0))
        select case(name)
        case('linear')
            line=linear_nd_interpolator(x,y,rescale=rescale,fill_value=fill_value,status=code)
            if (code==0) v=line%evaluate(xi,status=code)
        case('nearest')
            near=nearest_nd_interpolator(x,y,rescale=rescale,status=code)
            if (code==0) v=near%evaluate(xi,status=code)
        end select
        if (present(status)) status=code
    end function griddata_complex_matrix

end module scifort_interpolate_scattered
