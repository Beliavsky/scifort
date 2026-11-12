! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! Clough-Tocher C1 cubic interpolation on a planar Delaunay triangulation.
! The 2D minimum-curvature gradient equations and affine-invariant triangular
! Bernstein-Bezier patch formulas were adapted from the published SciPy
! _interpnd.pyx implementation (Copyright (C) Pauli Virtanen, 2010; BSD-3-Clause).
! See LICENSES/scipy-interpnd-BSD-3-Clause.txt for the retained BSD notice.
! Geometric triangulation uses the independent native SciFort Bowyer-Watson
! implementation rather than SciPy's Qhull; degeneracy tolerance therefore
! differs, and cocircular datasets can select a different Delaunay diagonal.
module scifort_interpolate_clough_tocher
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only: dp
    use scifort_interpolate_scattered, only: delaunay_2d, scattered_status_success, &
        scattered_status_invalid_input, scattered_status_degenerate
    implicit none
    private

    integer, parameter, public :: ct_status_success = scattered_status_success
    integer, parameter, public :: ct_status_invalid_input = scattered_status_invalid_input
    integer, parameter, public :: ct_status_degenerate = scattered_status_degenerate
    integer, parameter, public :: ct_status_nonconverged = 3

    type :: ct_core
        real(dp), allocatable :: points(:,:)        ! (n,2), optionally normalized
        complex(dp), allocatable :: values(:,:)    ! (n,channels)
        complex(dp), allocatable :: grads(:,:,:)   ! (n,channels,2)
        integer, allocatable :: triangles(:,:)     ! (3,ntri), CCW oriented
        integer, allocatable :: neighbors(:,:)     ! (3,ntri), opposite local vertex
        complex(dp), allocatable :: controls(:,:,:) ! (19,channels,ntri), cached Bezier controls
        real(dp), allocatable :: bary_map(:,:,:)   ! (3,3,ntri), affine barycentric maps
        integer, allocatable :: bin_start(:)      ! CSR offsets of uniform spatial bins
        integer, allocatable :: bin_triangles(:)  ! CSR triangle ids, ascending by triangle
        integer :: bin_shape(2) = 0              ! (nx,ny), zero if fallback to linear search
        real(dp) :: bin_lower(2) = 0.0_dp         ! minimum coordinates
        real(dp) :: bin_step(2) = 1.0_dp          ! bin widths
        real(dp) :: bin_upper(2) = 0.0_dp         ! maximum coordinates
        real(dp) :: bin_pad(2) = 0.0_dp           ! rescaling and site-roundoff padding
        real(dp) :: query_roundoff(2) = 0.0_dp    ! normalized-coordinate roundoff bound
        real(dp) :: offset(2) = 0.0_dp
        real(dp) :: scale(2) = 1.0_dp
        complex(dp) :: fill_value = (0.0_dp,0.0_dp)
        logical :: converged = .false.
    end type ct_core

    type, public :: clough_tocher_2d_interpolator_type
        private
        type(ct_core) :: core
    contains
        procedure :: valid => real_valid
        procedure :: converged => real_converged
        procedure :: channels => real_channels
        procedure :: simplices => real_simplices
        procedure, private :: real_one, real_many
        generic :: evaluate => real_one, real_many
        procedure :: value => real_value
    end type clough_tocher_2d_interpolator_type

    type, public :: complex_clough_tocher_2d_interpolator_type
        private
        type(ct_core) :: core
    contains
        procedure :: valid => complex_valid
        procedure :: converged => complex_converged
        procedure :: channels => complex_channels
        procedure :: simplices => complex_simplices
        procedure, private :: complex_one, complex_many
        generic :: evaluate => complex_one, complex_many
        procedure :: value => complex_value
    end type complex_clough_tocher_2d_interpolator_type

    public :: clough_tocher_2d_interpolator
    interface clough_tocher_2d_interpolator
        module procedure real_vector_constructor, real_matrix_constructor
        module procedure complex_vector_constructor, complex_matrix_constructor
    end interface

contains

    pure logical function core_valid(p) result(ok)
        type(ct_core), intent(in) :: p !! triangulated interpolation data
        ok = allocated(p%points) .and. allocated(p%values) .and. allocated(p%grads) .and. &
             allocated(p%triangles) .and. allocated(p%neighbors) .and. &
             allocated(p%controls) .and. allocated(p%bary_map)
    end function core_valid

    pure logical function real_valid(self) result(ok)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        ok = core_valid(self%core)
    end function real_valid

    pure logical function complex_valid(self) result(ok)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        ok = core_valid(self%core)
    end function complex_valid

    pure logical function real_converged(self) result(ok)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        ok = self%core%converged
    end function real_converged

    pure logical function complex_converged(self) result(ok)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        ok = self%core%converged
    end function complex_converged

    pure integer function real_channels(self) result(n)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        n = 0
        if (allocated(self%core%values)) n = size(self%core%values,2)
    end function real_channels

    pure integer function complex_channels(self) result(n)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        n = 0
        if (allocated(self%core%values)) n = size(self%core%values,2)
    end function complex_channels

    function real_simplices(self) result(tri)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        integer, allocatable :: tri(:,:)
        allocate(tri(3,0))
        if (self%valid()) tri = self%core%triangles
    end function real_simplices

    function complex_simplices(self) result(tri)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        integer, allocatable :: tri(:,:)
        allocate(tri(3,0))
        if (self%valid()) tri = self%core%triangles
    end function complex_simplices

    function real_vector_constructor(x,y,rescale,fill_value,tol,maxiter,status) result(out)
        real(dp), intent(in) :: x(:,:) !! (sites,2) scattered coordinates
        real(dp), intent(in) :: y(:) !! value at each site
        logical, intent(in), optional :: rescale !! rescale axes to comparable ranges
        real(dp), intent(in), optional :: fill_value !! out-of-hull fill, default NaN
        real(dp), intent(in), optional :: tol !! iteration tolerance, default 1e-6
        integer, intent(in), optional :: maxiter !! number of gradient iterations, default 400
        integer, intent(out), optional :: status !! construction code
        type(clough_tocher_2d_interpolator_type) :: out
        out = real_matrix_constructor(x,reshape(y,[size(y),1]),rescale,fill_value,tol,maxiter,status)
    end function real_vector_constructor

    function real_matrix_constructor(x,y,rescale,fill_value,tol,maxiter,status) result(out)
        real(dp), intent(in) :: x(:,:) !! (sites,2) scattered coordinates
        real(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in), optional :: rescale !! rescale axes to comparable ranges
        real(dp), intent(in), optional :: fill_value !! out-of-hull fill, default NaN
        real(dp), intent(in), optional :: tol !! convergence tolerance, default 1e-6
        integer, intent(in), optional :: maxiter !! gradient iterations, default 400
        integer, intent(out), optional :: status !! construction code
        type(clough_tocher_2d_interpolator_type) :: out
        integer :: code
        call initialize(out%core,x,cmplx(y,0.0_dp,dp),rescale,tol,maxiter,code)
        if (present(fill_value)) out%core%fill_value=cmplx(fill_value,0.0_dp,dp)
        if (present(status)) status=code
    end function real_matrix_constructor

    function complex_vector_constructor(x,y,rescale,fill_value,tol,maxiter,status) result(out)
        real(dp), intent(in) :: x(:,:) !! (sites,2) coordinates
        complex(dp), intent(in) :: y(:) !! site ordinates
        logical, intent(in), optional :: rescale !! rescale axes
        complex(dp), intent(in), optional :: fill_value !! out-of-hull fill
        real(dp), intent(in), optional :: tol !! convergence tolerance
        integer, intent(in), optional :: maxiter !! gradient iterations
        integer, intent(out), optional :: status !! construction code
        type(complex_clough_tocher_2d_interpolator_type) :: out
        out=complex_matrix_constructor(x,reshape(y,[size(y),1]),rescale,fill_value,tol,maxiter,status)
    end function complex_vector_constructor

    function complex_matrix_constructor(x,y,rescale,fill_value,tol,maxiter,status) result(out)
        real(dp), intent(in) :: x(:,:) !! (sites,2) coordinates
        complex(dp), intent(in) :: y(:,:) !! (sites,channels)
        logical, intent(in), optional :: rescale !! rescale axes
        complex(dp), intent(in), optional :: fill_value !! out-of-hull fill
        real(dp), intent(in), optional :: tol !! convergence tolerance
        integer, intent(in), optional :: maxiter !! gradient iterations
        integer, intent(out), optional :: status !! construction code
        type(complex_clough_tocher_2d_interpolator_type) :: out
        integer :: code
        call initialize(out%core,x,y,rescale,tol,maxiter,code)
        if (present(fill_value)) out%core%fill_value=fill_value
        if (present(status)) status=code
    end function complex_matrix_constructor

    subroutine initialize(p,x,y,rescale,tol,maxiter,code)
        type(ct_core), intent(out) :: p !! initialized Clough-Tocher core
        real(dp), intent(in) :: x(:,:) !! input sites, (n,2)
        complex(dp), intent(in) :: y(:,:) !! complex site ordinates (n,channels)
        logical, intent(in), optional :: rescale !! normalize coordinates
        real(dp), intent(in), optional :: tol !! gradient relative/absolute convergence
        integer, intent(in), optional :: maxiter !! iteration limit
        integer, intent(out) :: code !! status
        real(dp) :: t, nan_value
        integer :: n, i, j, iteration_limit
        logical :: scaled
        code=ct_status_invalid_input
        n=size(x,1)
        if (size(x,2)/=2 .or. n<3 .or. size(y,1)/=n .or. size(y,2)<1) return
        if (.not.all(ieee_is_finite(x))) return
        if (.not.all(ieee_is_finite(real(y,dp))) .or. .not.all(ieee_is_finite(aimag(y)))) return
        t=1.0e-6_dp
        if (present(tol)) t=tol
        iteration_limit=400
        if (present(maxiter)) iteration_limit=maxiter
        if (.not.ieee_is_finite(t) .or. t<=0.0_dp .or. iteration_limit<1) return
        scaled=.false.
        if (present(rescale)) scaled=rescale
        p%offset=0.0_dp
        p%scale=1.0_dp
        if (scaled) then
            p%offset=sum(x,dim=1)/real(n,dp)
            p%scale=maxval(x,dim=1)-minval(x,dim=1)
            if (any(p%scale<=0.0_dp)) return
        end if
        allocate(p%points(n,2),p%values(n,size(y,2)))
        p%points=(x-spread(p%offset,1,n))/spread(p%scale,1,n)
        p%values=y
        ! Converting physical queries to normalized coordinates can lose low
        ! bits when an axis has a large origin relative to its extent.
        p%query_roundoff=16.0_dp*epsilon(1.0_dp)*(abs(p%offset)/p%scale+1.0_dp)
        nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
        p%fill_value=cmplx(nan_value,0.0_dp,dp)
        do i=1,n
            do j=i+1,n
                if (maxval(abs(p%points(i,:)-p%points(j,:)))<=tiny(1.0_dp)) then
                    call reset(p)
                    return
                end if
            end do
        end do
        call delaunay_2d(p%points,p%triangles,code)
        if (code/=ct_status_success) then
            call reset(p)
            return
        end if
        call build_topology(p,code)
        if (code/=ct_status_success) then
            call reset(p)
            return
        end if
        allocate(p%grads(n,size(y,2),2))
        call estimate_gradients(p,t,iteration_limit,p%converged,code)
        if (code/=ct_status_success .and. code/=ct_status_nonconverged) then
            call reset(p)
        else
            call prepare_patches(p)
            call build_spatial_index(p)
        end if
    end subroutine initialize

    subroutine reset(p)
        type(ct_core), intent(inout) :: p !! interpolator to invalidate on error
        if (allocated(p%points)) deallocate(p%points)
        if (allocated(p%values)) deallocate(p%values)
        if (allocated(p%grads)) deallocate(p%grads)
        if (allocated(p%triangles)) deallocate(p%triangles)
        if (allocated(p%neighbors)) deallocate(p%neighbors)
        if (allocated(p%controls)) deallocate(p%controls)
        if (allocated(p%bary_map)) deallocate(p%bary_map)
        if (allocated(p%bin_start)) deallocate(p%bin_start)
        if (allocated(p%bin_triangles)) deallocate(p%bin_triangles)
        p%bin_shape = 0
    end subroutine reset

    ! Edge sort allows mesh adjacency and per-vertex neighbors to be built in
    ! O(ntri log(ntri)) without quadratic per-site scratch arrays.
    recursive subroutine sort_edges(edges,left,right)
        integer, intent(inout) :: edges(:,:) !! 4 x nedge [endpoint1,endpoint2,triangle,opposite]
        integer, intent(in) :: left !! first sorted column
        integer, intent(in) :: right !! last sorted column
        integer :: i,j,m,tmp(4),pivot(2)
        if (left>=right) return
        i=left
        j=right
        m=(left+right)/2
        pivot=edges(1:2,m)
        do
            do while (i<=right)
                if (.not.less_edge(edges(1:2,i),pivot)) exit
                i=i+1
            end do
            do while (j>=left)
                if (.not.less_edge(pivot,edges(1:2,j))) exit
                j=j-1
            end do
            if (i>j) exit
            tmp=edges(:,i)
            edges(:,i)=edges(:,j)
            edges(:,j)=tmp
            i=i+1
            j=j-1
            if (i>j) exit
        end do
        if (left<j) call sort_edges(edges,left,j)
        if (i<right) call sort_edges(edges,i,right)
    end subroutine sort_edges

    pure logical function less_edge(a,b) result(less)
        integer, intent(in) :: a(2) !! first edge endpoints
        integer, intent(in) :: b(2) !! second edge endpoints
        less=(a(1)<b(1)) .or. (a(1)==b(1) .and. a(2)<b(2))
    end function less_edge

    subroutine build_topology(p,code)
        type(ct_core), intent(inout) :: p !! core with triangulation, augmented with neighboring triangles
        integer, intent(out) :: code !! zero if manifold triangulation and valid edge list
        integer, allocatable :: edges(:,:),degree(:)
        integer :: n, nt, ne, t,k,a,b,i,j
        n=size(p%points,1)
        nt=size(p%triangles,2)
        ne=3*nt
        allocate(edges(4,ne),p%neighbors(3,nt),degree(n))
        p%neighbors=0
        degree=0
        do t=1,nt
            do k=1,3
                a=p%triangles(mod(k,3)+1,t)
                b=p%triangles(mod(k+1,3)+1,t)
                edges(:,3*(t-1)+k)=[min(a,b),max(a,b),t,k]
            end do
        end do
        call sort_edges(edges,1,ne)
        code=ct_status_degenerate
        i=1
        do while (i<=ne)
            j=i+1
            do while (j<=ne)
                if (any(edges(1:2,j)/=edges(1:2,i))) exit
                j=j+1
            end do
            if (j-i>2) return   ! non-manifold edge
            a=edges(1,i)
            b=edges(2,i)
            degree(a)=degree(a)+1
            degree(b)=degree(b)+1
            if (j-i==2) then
                p%neighbors(edges(4,i),edges(3,i))=edges(3,i+1)
                p%neighbors(edges(4,i+1),edges(3,i+1))=edges(3,i)
            end if
            i=j
        end do
        ! A compressed neighbor array is only needed for gradient estimation.
        ! Its construction remains O(ntri log(ntri)) from the sorted edge list.
        code=ct_status_success
    end subroutine build_topology

    subroutine estimate_gradients(p,tol,maxiter,converged,code)
        type(ct_core), intent(inout) :: p !! stores estimated vertex gradients
        real(dp), intent(in) :: tol !! convergence threshold
        integer, intent(in) :: maxiter !! maximum Gauss-Seidel passes
        logical, intent(out) :: converged !! convergence indicator
        integer, intent(out) :: code !! success, degenerate, nonconvergence
        integer, allocatable :: edges(:,:),offset(:),degree(:),next(:),adj(:)
        integer :: n,nt,t,k,a,b,i,iteration,channel,j,other
        real(dp) :: q11,q12,q22,det,ex,ey,ll3,err,change
        complex(dp) :: sx,sy, proj, dx,dy, f1,f2
        n=size(p%points,1)
        nt=size(p%triangles,2)
        allocate(edges(2,3*nt),degree(n))
        degree=0
        do t=1,nt
            do k=1,3
                a=p%triangles(mod(k,3)+1,t)
                b=p%triangles(mod(k+1,3)+1,t)
                edges(:,3*(t-1)+k)=[min(a,b),max(a,b)]
            end do
        end do
        call sort_pairs(edges,1,3*nt)
        ! Build CSR adjacency with a deterministic orientation/order.
        do k=1,size(edges,2)
            if (k>1) then
                if (all(edges(:,k)==edges(:,k-1))) cycle
            end if
            degree(edges(1,k))=degree(edges(1,k))+1
            degree(edges(2,k))=degree(edges(2,k))+1
        end do
        allocate(offset(n+1),next(n))
        offset(1)=1
        do i=1,n
            offset(i+1)=offset(i)+degree(i)
        end do
        next=offset(:n)
        allocate(adj(offset(n+1)-1))
        do k=1,size(edges,2)
            if (k>1) then
                if (all(edges(:,k)==edges(:,k-1))) cycle
            end if
            a=edges(1,k)
            b=edges(2,k)
            adj(next(a))=b
            next(a)=next(a)+1
            adj(next(b))=a
            next(b)=next(b)+1
        end do
        p%grads=(0.0_dp,0.0_dp)
        code=ct_status_success
        converged=.false.
        do iteration=1,maxiter
            err=0.0_dp
            do i=1,n
                q11=0.0_dp
                q12=0.0_dp
                q22=0.0_dp
                do j=offset(i),offset(i+1)-1
                    other=adj(j)
                    ex=p%points(other,1)-p%points(i,1)
                    ey=p%points(other,2)-p%points(i,2)
                    ll3=(ex*ex+ey*ey)**1.5_dp
                    q11=q11+4.0_dp*ex*ex/ll3
                    q12=q12+4.0_dp*ex*ey/ll3
                    q22=q22+4.0_dp*ey*ey/ll3
                end do
                det=q11*q22-q12*q12
                if (.not.ieee_is_finite(det) .or. det<=epsilon(det)*max(q11*q22,1.0e-100_dp)) then
                    code=ct_status_degenerate
                    return
                end if
                do channel=1,size(p%values,2)
                    sx=(0.0_dp,0.0_dp)
                    sy=(0.0_dp,0.0_dp)
                    f1=p%values(i,channel)
                    do j=offset(i),offset(i+1)-1
                        other=adj(j)
                        ex=p%points(other,1)-p%points(i,1)
                        ey=p%points(other,2)-p%points(i,2)
                        ll3=(ex*ex+ey*ey)**1.5_dp
                        f2=p%values(other,channel)
                        proj=ex*p%grads(other,channel,1)+ey*p%grads(other,channel,2)
                        sx=sx+(6.0_dp*(f1-f2)+2.0_dp*proj)*ex/ll3
                        sy=sy+(6.0_dp*(f1-f2)+2.0_dp*proj)*ey/ll3
                    end do
                    dx=-(q22*sx-q12*sy)/det
                    dy=-(-q12*sx+q11*sy)/det
                    change=max(abs(dx-p%grads(i,channel,1)),abs(dy-p%grads(i,channel,2)))
                    err=max(err,change/max(1.0_dp,abs(dx),abs(dy)))
                    p%grads(i,channel,1)=dx
                    p%grads(i,channel,2)=dy
                end do
            end do
            if (err<tol) then
                converged=.true.
                return
            end if
        end do
        code=ct_status_nonconverged
    end subroutine estimate_gradients

    recursive subroutine sort_pairs(pairs,left,right)
        integer, intent(inout) :: pairs(:,:) !! sorted endpoint pair columns
        integer, intent(in) :: left !! lower column
        integer, intent(in) :: right !! upper column
        integer :: i,j,m,tmp(2),pivot(2)
        if (left>=right) return
        i=left
        j=right
        m=(left+right)/2
        pivot=pairs(:,m)
        do
            do while (i<=right)
                if (.not.less_edge(pairs(:,i),pivot)) exit
                i=i+1
            end do
            do while (j>=left)
                if (.not.less_edge(pivot,pairs(:,j))) exit
                j=j-1
            end do
            if (i>j) exit
            tmp=pairs(:,i)
            pairs(:,i)=pairs(:,j)
            pairs(:,j)=tmp
            i=i+1
            j=j-1
            if (i>j) exit
        end do
        if (left<j) call sort_pairs(pairs,left,j)
        if (i<right) call sort_pairs(pairs,i,right)
    end subroutine sort_pairs

    pure real(dp) function orient(a,b,c) result(o)
        real(dp), intent(in) :: a(2) !! first vertex
        real(dp), intent(in) :: b(2) !! second vertex
        real(dp), intent(in) :: c(2) !! third vertex
        o=(b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1))
    end function orient

    pure function barycentric(p,a,b,c) result(w)
        real(dp), intent(in) :: p(2) !! query
        real(dp), intent(in) :: a(2) !! vertex 1
        real(dp), intent(in) :: b(2) !! vertex 2
        real(dp), intent(in) :: c(2) !! vertex 3
        real(dp) :: w(3),den
        den=orient(a,b,c)
        w(1)=orient(b,c,p)/den
        w(2)=orient(c,a,p)/den
        w(3)=1.0_dp-w(1)-w(2)
    end function barycentric

    ! Store the cubic Bezier coefficients once for each triangle and channel.
    ! This removes all of the neighbor-centroid and vertex-gradient work from
    ! individual queries. The 19 coefficients are the nonzero basis terms on
    ! the three Clough-Tocher subtriangles (one barycentric component vanishes).
    subroutine prepare_patches(p)
        type(ct_core), intent(inout) :: p !! converged (or usable) gradient interpolation core
        real(dp) :: a(2),b(2),c(2),det
        integer :: t,v(3),j
        allocate(p%controls(19,size(p%values,2),size(p%triangles,2)))
        allocate(p%bary_map(3,3,size(p%triangles,2)))
        do t=1,size(p%triangles,2)
            v=p%triangles(:,t)
            a=p%points(v(1),:)
            b=p%points(v(2),:)
            c=p%points(v(3),:)
            det=orient(a,b,c)
            p%bary_map(1,:,t)=[(b(2)-c(2))/det,(c(1)-b(1))/det, &
                              (b(1)*c(2)-c(1)*b(2))/det]
            p%bary_map(2,:,t)=[(c(2)-a(2))/det,(a(1)-c(1))/det, &
                              (c(1)*a(2)-a(1)*c(2))/det]
            p%bary_map(3,:,t)=[-p%bary_map(1,1,t)-p%bary_map(2,1,t), &
                              -p%bary_map(1,2,t)-p%bary_map(2,2,t), &
                              1.0_dp-p%bary_map(1,3,t)-p%bary_map(2,3,t)]
            do j=1,size(p%values,2)
                call make_patch_controls(p,t,j,p%controls(:,j,t))
            end do
        end do
    end subroutine prepare_patches

    pure function patch_weights(b) result(w)
        real(dp), intent(in) :: b(3) !! triangle barycentric coordinates
        real(dp) :: w(19),u(4),m
        m=minval(b)
        u(:3)=b-m
        u(4)=3.0_dp*m
        w(1:4)=u**3
        w(5:16)=3.0_dp*[u(1)**2*u(2),u(1)**2*u(3),u(1)**2*u(4), &
                        u(2)**2*u(1),u(2)**2*u(3),u(2)**2*u(4), &
                        u(3)**2*u(1),u(3)**2*u(2),u(3)**2*u(4), &
                        u(4)**2*u(1),u(4)**2*u(2),u(4)**2*u(3)]
        w(17:19)=6.0_dp*u(4)*[u(1)*u(2),u(1)*u(3),u(2)*u(3)]
    end function patch_weights

    ! Construct a CSR spatial index over triangle bounding boxes. The number
    ! of bins scales with ntri; pathological long thin triangles may touch too
    ! many bins, in which case the exact linear search is retained instead.
    ! Triangle ids within each bin remain increasing, preserving the original
    ! deterministic choice on shared edges.
    subroutine build_spatial_index(p)
        type(ct_core), intent(inout) :: p !! core with triangulation and vertices
        integer, allocatable :: counts(:),cursor(:)
        integer :: nt,nx,ny,ncell,t,j,i,v(3),lo(2),hi(2),cell,entries
        real(dp) :: pad(2),xmin(2),xmax(2),extent(2)
        nt=size(p%triangles,2)
        nx=max(1,min(128,int(sqrt(real(nt,dp)/2.0_dp))))
        ny=nx
        ncell=nx*ny
        p%bin_lower=minval(p%points,dim=1)
        p%bin_upper=maxval(p%points,dim=1)
        extent=p%bin_upper-p%bin_lower
        if (any(extent<=0.0_dp)) return
        p%bin_step=extent/real([nx,ny],dp)
        pad=4.0e-12_dp*extent+p%query_roundoff
        p%bin_pad=pad
        allocate(counts(ncell))
        counts=0
        do t=1,nt
            v=p%triangles(:,t)
            xmin=minval(p%points(v,:),dim=1)-pad
            xmax=maxval(p%points(v,:),dim=1)+pad
            lo=max(1,min([nx,ny],floor((xmin-p%bin_lower)/p%bin_step)+1))
            hi=max(1,min([nx,ny],floor((xmax-p%bin_lower)/p%bin_step)+1))
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    counts((j-1)*nx+i)=counts((j-1)*nx+i)+1
                end do
            end do
        end do
        entries=sum(counts)
        if (entries>32*nt) return
        allocate(p%bin_start(ncell+1),p%bin_triangles(entries),cursor(ncell))
        p%bin_start(1)=1
        do i=1,ncell
            p%bin_start(i+1)=p%bin_start(i)+counts(i)
        end do
        cursor=p%bin_start(:ncell)
        do t=1,nt
            v=p%triangles(:,t)
            xmin=minval(p%points(v,:),dim=1)-pad
            xmax=maxval(p%points(v,:),dim=1)+pad
            lo=max(1,min([nx,ny],floor((xmin-p%bin_lower)/p%bin_step)+1))
            hi=max(1,min([nx,ny],floor((xmax-p%bin_lower)/p%bin_step)+1))
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    cell=(j-1)*nx+i
                    p%bin_triangles(cursor(cell))=t
                    cursor(cell)=cursor(cell)+1
                end do
            end do
        end do
        p%bin_shape=[nx,ny]
    end subroutine build_spatial_index

    pure function triangle_barycentric(p,t,q) result(b)
        type(ct_core), intent(in) :: p !! core with vertices
        integer, intent(in) :: t !! triangle index
        real(dp), intent(in) :: q(2) !! normalized query
        real(dp) :: b(3)
        integer :: v(3)
        v=p%triangles(:,t)
        b=barycentric(q,p%points(v(1),:),p%points(v(2),:),p%points(v(3),:))
    end function triangle_barycentric

    subroutine locate_triangle(p,q,t,b)
        type(ct_core), intent(in) :: p !! indexed triangulation
        real(dp), intent(in) :: q(2) !! normalized query
        integer, intent(out) :: t !! triangle id or zero outside the convex hull
        real(dp), intent(out) :: b(3) !! barycentric coordinates if inside
        integer :: cell,i,j,k,candidate
        real(dp) :: tol
        t=0
        b=0.0_dp
        if (allocated(p%bin_start)) then
            if (any(q<p%bin_lower-p%bin_pad) .or. any(q>p%bin_upper+p%bin_pad)) return
            i=max(1,min(p%bin_shape(1),floor((q(1)-p%bin_lower(1))/p%bin_step(1))+1))
            j=max(1,min(p%bin_shape(2),floor((q(2)-p%bin_lower(2))/p%bin_step(2))+1))
            cell=(j-1)*p%bin_shape(1)+i
            do k=p%bin_start(cell),p%bin_start(cell+1)-1
                candidate=p%bin_triangles(k)
                b=triangle_barycentric(p,candidate,q)
                tol=1.0e-12_dp+ &
                    p%query_roundoff(1)*maxval(abs(p%bary_map(:,1,candidate)))+ &
                    p%query_roundoff(2)*maxval(abs(p%bary_map(:,2,candidate)))
                if (minval(b)<-tol) cycle
                t=candidate
                return
            end do
        else
            do candidate=1,size(p%triangles,2)
                b=triangle_barycentric(p,candidate,q)
                tol=1.0e-12_dp+ &
                    p%query_roundoff(1)*maxval(abs(p%bary_map(:,1,candidate)))+ &
                    p%query_roundoff(2)*maxval(abs(p%bary_map(:,2,candidate)))
                if (minval(b)<-tol) cycle
                t=candidate
                return
            end do
        end if
    end subroutine locate_triangle

    function evaluate_core(p,x,status) result(out)
        type(ct_core), intent(in) :: p !! constructed interpolator
        real(dp), intent(in) :: x(:) !! query coordinates, length 2
        integer, intent(out) :: status !! success or malformed query
        complex(dp), allocatable :: out(:)
        real(dp) :: q(2),b(3),weights(19)
        integer :: t,channel
        allocate(out(0))
        status=ct_status_invalid_input
        if (.not.core_valid(p) .or. size(x)/=2) return
        if (.not.all(ieee_is_finite(x))) return
        q=(x-p%offset)/p%scale
        deallocate(out)
        allocate(out(size(p%values,2)))
        out=p%fill_value
        status=ct_status_success
        call locate_triangle(p,q,t,b)
        if (t==0) return
        weights=patch_weights(b)
        do channel=1,size(p%values,2)
            out(channel)=sum(p%controls(:,channel,t)*weights)
        end do
    end function evaluate_core

    subroutine make_patch_controls(p,t,channel,c)
        type(ct_core), intent(in) :: p !! triangulation, gradients, ordinates
        integer, intent(in) :: t !! selected triangle
        integer, intent(in) :: channel !! output column
        complex(dp), intent(out) :: c(19) !! cached tetrahedral cubic Bernstein controls
        complex(dp) :: f(3),gf(3,2),c3000,c0300,c0030,c0003,c2100,c2010,c2001, &
            c0210,c0201,c0021,c1200,c1020,c1002,c0120,c0102,c0012,c1101,c1011,c0111, &
            d12,d13,d21,d23,d31,d32
        real(dp) :: e12(2),e23(2),e31(2),g(3),centroid(2),bc(3)
        integer :: v(3),i,j,k,ntri
        v=p%triangles(:,t)
        do i=1,3
            f(i)=p%values(v(i),channel)
            gf(i,:)=p%grads(v(i),channel,:)
        end do
        e12=p%points(v(2),:)-p%points(v(1),:)
        e23=p%points(v(3),:)-p%points(v(2),:)
        e31=p%points(v(1),:)-p%points(v(3),:)
        d12=sum(gf(1,:)*e12)
        d21=-sum(gf(2,:)*e12)
        d23=sum(gf(2,:)*e23)
        d32=-sum(gf(3,:)*e23)
        d31=sum(gf(3,:)*e31)
        d13=-sum(gf(1,:)*e31)
        c3000=f(1)
        c2100=c3000+d12/3.0_dp
        c2010=c3000+d13/3.0_dp
        c0300=f(2)
        c1200=c0300+d21/3.0_dp
        c0210=c0300+d23/3.0_dp
        c0030=f(3)
        c1020=c0030+d31/3.0_dp
        c0120=c0030+d32/3.0_dp
        c2001=(c2100+c2010+c3000)/3.0_dp
        c0201=(c1200+c0300+c0210)/3.0_dp
        c0021=(c1020+c0120+c0030)/3.0_dp
        do k=1,3
            ntri=p%neighbors(k,t)
            if (ntri==0) then
                g(k)=-0.5_dp
            else
                centroid=0.0_dp
                do j=1,3
                    centroid=centroid+p%points(p%triangles(j,ntri),:)/3.0_dp
                end do
                bc=barycentric(centroid,p%points(v(1),:),p%points(v(2),:),p%points(v(3),:))
                select case(k)
                case(1)
                    g(k)=(2.0_dp*bc(3)+bc(2)-1.0_dp)/(2.0_dp-3.0_dp*bc(3)-3.0_dp*bc(2))
                case(2)
                    g(k)=(2.0_dp*bc(1)+bc(3)-1.0_dp)/(2.0_dp-3.0_dp*bc(1)-3.0_dp*bc(3))
                case(3)
                    g(k)=(2.0_dp*bc(2)+bc(1)-1.0_dp)/(2.0_dp-3.0_dp*bc(2)-3.0_dp*bc(1))
                end select
            end if
        end do
        c0111=(g(1)*(-c0300+3.0_dp*c0210-3.0_dp*c0120+c0030) &
                -c0300+2.0_dp*c0210-c0120+c0021+c0201)/2.0_dp
        c1011=(g(2)*(-c0030+3.0_dp*c1020-3.0_dp*c2010+c3000) &
                -c0030+2.0_dp*c1020-c2010+c2001+c0021)/2.0_dp
        c1101=(g(3)*(-c3000+3.0_dp*c2100-3.0_dp*c1200+c0300) &
                -c3000+2.0_dp*c2100-c1200+c2001+c0201)/2.0_dp
        c1002=(c1101+c1011+c2001)/3.0_dp
        c0102=(c1101+c0111+c0201)/3.0_dp
        c0012=(c1011+c0111+c0021)/3.0_dp
        c0003=(c1002+c0102+c0012)/3.0_dp
        c=[c3000,c0300,c0030,c0003, &
           c2100,c2010,c2001,c1200,c0210,c0201,c1020,c0120,c0021, &
           c1002,c0102,c0012,c1101,c1011,c0111]
    end subroutine make_patch_controls

    function evaluate_many(p,x,status) result(out)
        type(ct_core), intent(in) :: p !! interpolator
        real(dp), intent(in) :: x(:,:) !! (query,2)
        integer, intent(out) :: status !! status code
        complex(dp), allocatable :: out(:,:),tmp(:)
        integer :: i,code
        status=ct_status_invalid_input
        allocate(out(0,0))
        if (.not.core_valid(p) .or. size(x,2)/=2) return
        deallocate(out)
        allocate(out(size(x,1),size(p%values,2)))
        status=ct_status_success
        do i=1,size(x,1)
            tmp=evaluate_core(p,x(i,:),code)
            if (code/=ct_status_success) then
                status=code
                return
            end if
            out(i,:)=tmp
        end do
    end function evaluate_many

    function real_one(self,x,status) result(out)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        real(dp), intent(in) :: x(:) !! position, length 2
        integer, intent(out), optional :: status !! status code
        real(dp), allocatable :: out(:)
        complex(dp), allocatable :: tmp(:)
        integer :: code
        tmp=evaluate_core(self%core,x,code)
        out=real(tmp,dp)
        if (present(status)) status=code
    end function real_one

    function real_many(self,x,status) result(out)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        real(dp), intent(in) :: x(:,:) !! query array, (query,2)
        integer, intent(out), optional :: status !! status code
        real(dp), allocatable :: out(:,:)
        complex(dp), allocatable :: tmp(:,:)
        integer :: code
        tmp=evaluate_many(self%core,x,code)
        out=real(tmp,dp)
        if (present(status)) status=code
    end function real_many

    function complex_one(self,x,status) result(out)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        real(dp), intent(in) :: x(:) !! query location, length 2
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: out(:)
        integer :: code
        out=evaluate_core(self%core,x,code)
        if (present(status)) status=code
    end function complex_one

    function complex_many(self,x,status) result(out)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        real(dp), intent(in) :: x(:,:) !! (queries,2)
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: out(:,:)
        integer :: code
        out=evaluate_many(self%core,x,code)
        if (present(status)) status=code
    end function complex_many

    function real_value(self,x,channel,status) result(out)
        class(clough_tocher_2d_interpolator_type), intent(in) :: self !! real interpolator
        real(dp), intent(in) :: x(:) !! query, (2)
        integer, intent(in), optional :: channel !! output channel, default 1
        integer, intent(out), optional :: status !! status code
        real(dp) :: out
        real(dp), allocatable :: tmp(:)
        integer :: i,code
        tmp=self%evaluate(x,status=code)
        i=1
        if (present(channel)) i=channel
        out=ieee_value(0.0_dp,ieee_quiet_nan)
        if (i<1 .or. i>size(tmp)) code=ct_status_invalid_input
        if (code==ct_status_success) out=tmp(i)
        if (present(status)) status=code
    end function real_value

    function complex_value(self,x,channel,status) result(out)
        class(complex_clough_tocher_2d_interpolator_type), intent(in) :: self !! complex interpolator
        real(dp), intent(in) :: x(:) !! query, (2)
        integer, intent(in), optional :: channel !! output channel, default 1
        integer, intent(out), optional :: status !! status code
        complex(dp) :: out
        complex(dp), allocatable :: tmp(:)
        integer :: i,code
        tmp=self%evaluate(x,status=code)
        i=1
        if (present(channel)) i=channel
        out=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (i<1 .or. i>size(tmp)) code=ct_status_invalid_input
        if (code==ct_status_success) out=tmp(i)
        if (present(status)) status=code
    end function complex_value

end module scifort_interpolate_clough_tocher
