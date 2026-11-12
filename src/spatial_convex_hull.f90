! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Native 2D/3D convex hulls: planar monotone chain and validated 3D facets.
! No Qhull, BLAS/LAPACK, or other external implementation is required.
module scifort_spatial_convex_hull
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : spatial_distance_success, spatial_distance_invalid_input
    implicit none
    private

    real(dp), parameter :: orientation_tol = 128.0_dp * epsilon(1.0_dp)

    type, public :: convex_hull_type
        real(dp), allocatable :: points(:,:)    !! original sites (n,d)
        integer, allocatable :: vertices(:)     !! 1-based input indices; CCW for d=2
        integer, allocatable :: simplices(:,:) !! (d,nfacets): edges/triangles
        integer, allocatable :: neighbors(:,:) !! (d,nfacets): facet opposite local vertex
        real(dp), allocatable :: equations(:,:)!! (d+1,nfacets): outward unit normal, offset
        real(dp) :: area = 0.0_dp              !! 2D perimeter / 3D surface area
        real(dp) :: volume = 0.0_dp            !! 2D enclosed area / 3D volume
        integer :: status = spatial_distance_invalid_input
    contains
        procedure :: valid => hull_valid
    end type convex_hull_type

    public :: convex_hull
contains
    pure logical function hull_valid(self)
        class(convex_hull_type), intent(in) :: self !! hull to inspect
        hull_valid = self%status == spatial_distance_success .and. &
            allocated(self%simplices) .and. allocated(self%equations) .and. &
            allocated(self%vertices) .and. allocated(self%neighbors)
    end function hull_valid

    ! Form the complete hull of points (n,d). Only full-dimensional d=2,3
    ! geometries are supported. Invalid/degenerate inputs return status=1.
    function convex_hull(points, status) result(hull)
        real(dp), intent(in) :: points(:,:) !! input coordinates (n,d)
        integer, intent(out), optional :: status !! 0=success, 1=invalid geometry
        type(convex_hull_type) :: hull
        real(dp), allocatable :: scaled(:,:)
        real(dp) :: offset(3), span(3), scale, lo, hi
        logical :: ok
        integer :: d, n, j

        if (present(status)) status = spatial_distance_invalid_input
        n = size(points,1)
        d = size(points,2)
        if (d < 2 .or. d > 3 .or. n < d+1) return
        if (.not. all(ieee_is_finite(points))) return

        ! Centering/scaling are used ONLY internally. An isotropic scale
        ! preserves the true Euclidean geometry and facet unit normals.
        allocate(scaled(n,d))
        do j=1,d
            lo = minval(points(:,j))
            hi = maxval(points(:,j))
            offset(j) = 0.5_dp*lo + 0.5_dp*hi
            span(j) = hi-lo
        end do
        scale = maxval(span(:d))
        if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
        do j=1,d
            scaled(:,j) = (points(:,j)-offset(j))/scale
        end do
        if (.not. all(ieee_is_finite(scaled))) return

        ok = .false.
        select case(d)
        case(2)
            call build_planar(hull,scaled,points,scale,ok)
        case(3)
            call build_solid(hull,scaled,points,scale,ok)
            if (ok) call validate_solid(hull,scaled,ok)
            if (.not. ok) then
                hull = convex_hull_type()
                call retry_solid(hull,scaled,points,scale,ok)
            end if
            if (.not. ok) then
                ! A nearly coplanar horizon can invalidate an incremental
                ! triangulation. Rebuild from globally supporting planes;
                ! never return an incomplete hull as successful.
                hull = convex_hull_type()
                call build_solid_fallback(hull,scaled,points,scale,ok)
                if (ok) call validate_solid(hull,scaled,ok)
            end if
        end select
        if (.not. ok) return
        hull%points = points
        hull%status = spatial_distance_success
        if (present(status)) status = spatial_distance_success
    end function convex_hull

    pure real(dp) function cross2(a,b,c) result(z)
        real(dp), intent(in) :: a(2) !! first point
        real(dp), intent(in) :: b(2) !! second point
        real(dp), intent(in) :: c(2) !! third point
        z = (b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1))
    end function cross2

    pure function cross3(a,b) result(c)
        real(dp), intent(in) :: a(3) !! first vector
        real(dp), intent(in) :: b(3) !! second vector
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
    end function cross3

    subroutine build_planar(hull,p,original,scale,ok,tolerance)
        type(convex_hull_type), intent(inout) :: hull !! populated 2D hull
        real(dp), intent(in) :: p(:,:) !! centered/scaled points (n,2)
        real(dp), intent(in) :: original(:,:) !! original points (n,2)
        real(dp), intent(in) :: scale !! coordinate scaling factor
        logical, intent(out) :: ok !! true when full-dimensional
        real(dp), intent(in), optional :: tolerance !! minimum positive orientation
        real(dp) :: min_turn
        integer, allocatable :: order(:), sorted(:), chain(:)
        real(dp) :: direction(2), normal(2), length, twice_area
        integer :: n, i, j, k, id, t, nuniq, nh, a, b, limit

        ok = .false.
        min_turn = orientation_tol
        if (present(tolerance)) min_turn = tolerance
        n = size(p,1)
        allocate(order(n), sorted(n), chain(2*n))
        do i=1,n
            order(i) = i
        end do
        ! Stable lexicographic insertion sort for predictable hull order.
        do i=2,n
            id = order(i)
            j = i-1
            do while(j >= 1)
                t = order(j)
                if (p(t,1) < p(id,1)) exit
                if (p(t,1) == p(id,1) .and. p(t,2) <= p(id,2)) exit
                order(j+1) = order(j)
                j = j-1
            end do
            order(j+1) = id
        end do
        nuniq = 0
        do i=1,n
            id = order(i)
            if (nuniq > 0) then
                t = sorted(nuniq)
                if (all(p(t,:)==p(id,:))) cycle
            end if
            nuniq = nuniq+1
            sorted(nuniq) = id
        end do
        if (nuniq < 3) return

        k = 0
        do i=1,nuniq
            id = sorted(i)
            do while(k >= 2)
                if (cross2(p(chain(k-1),:),p(chain(k),:),p(id,:)) > min_turn) exit
                k = k-1
            end do
            k = k+1
            chain(k) = id
        end do
        limit = k+1
        do i=nuniq-1,1,-1
            id = sorted(i)
            do while(k >= limit)
                if (cross2(p(chain(k-1),:),p(chain(k),:),p(id,:)) > min_turn) exit
                k = k-1
            end do
            k = k+1
            chain(k) = id
        end do
        nh = k-1
        if (nh < 3) return
        if (cross2(p(chain(1),:),p(chain(2),:),p(chain(3),:)) <= 0.0_dp) return
        allocate(hull%vertices(nh),hull%simplices(2,nh), &
                 hull%neighbors(2,nh),hull%equations(3,nh))
        hull%vertices = chain(1:nh)
        twice_area = 0.0_dp
        do i=1,nh
            j = mod(i,nh)+1
            a = chain(i)
            b = chain(j)
            direction = p(b,:)-p(a,:)
            length = norm2(direction)
            if (length <= 0.0_dp) return
            normal = [direction(2),-direction(1)]/length
            hull%simplices(:,i) = [a,b]
            ! Neighbor opposite first endpoint touches the second endpoint.
            hull%neighbors(:,i) = [j,mod(i+nh-2,nh)+1]
            hull%equations(1:2,i) = normal
            hull%equations(3,i) = -dot_product(normal,original(a,:))
            hull%area = hull%area + length*scale
            twice_area = twice_area + cross2(p(chain(1),:),p(a,:),p(b,:))
        end do
        hull%volume = 0.5_dp*twice_area*scale**2
        ok = hull%volume > 0.0_dp .and. ieee_is_finite(hull%area) .and. &
             ieee_is_finite(hull%volume)
    end subroutine build_planar

    pure subroutine orient_triangle(tri,p,inside)
        integer, intent(inout) :: tri(3) !! facet indices, oriented outward
        real(dp), intent(in) :: p(:,:) !! scaled point coordinates
        real(dp), intent(in) :: inside(3) !! guaranteed interior point
        real(dp) :: normal(3)
        integer :: t
        normal = cross3(p(tri(2),:)-p(tri(1),:),p(tri(3),:)-p(tri(1),:))
        if (dot_product(normal,inside-p(tri(1),:)) > 0.0_dp) then
            t = tri(2)
            tri(2) = tri(3)
            tri(3) = t
        end if
    end subroutine orient_triangle

    subroutine build_solid(hull,p,original,scale,ok)
        type(convex_hull_type), intent(inout) :: hull !! populated 3D hull
        real(dp), intent(in) :: p(:,:) !! centered/scaled points (n,3)
        real(dp), intent(in) :: original(:,:) !! original coordinates (n,3)
        real(dp), intent(in) :: scale !! isotropic coordinate scale
        logical, intent(out) :: ok !! full-dimensional and watertight
        integer, allocatable :: faces(:,:), horizon(:,:), triangles(:,:)
        logical, allocatable :: active(:), visible(:), chosen(:), is_vertex(:)
        real(dp) :: inside(3), normal(3), facearea, best, magnitude, distance, volume_sum
        integer :: n, capacity, seed(4), i, f, g, nf, ne, e, t, u, a, b, c, idx, nfaces, h, v

        ok = .false.
        n = size(p,1)
        seed(1) = 1
        best = -1.0_dp
        do i=2,n
            distance = sum((p(i,:)-p(seed(1),:))**2)
            if (distance > best) then
                best = distance
                seed(2) = i
            end if
        end do
        if (best <= orientation_tol**2) return
        best = -1.0_dp
        do i=1,n
            if (i == seed(1) .or. i == seed(2)) cycle
            distance = norm2(cross3(p(seed(2),:)-p(seed(1),:),p(i,:)-p(seed(1),:)))
            if (distance > best) then
                best = distance
                seed(3) = i
            end if
        end do
        if (best <= orientation_tol) return
        normal = cross3(p(seed(2),:)-p(seed(1),:),p(seed(3),:)-p(seed(1),:))
        magnitude = norm2(normal)
        best = -1.0_dp
        do i=1,n
            if (any(seed(1:3) == i)) cycle
            distance = abs(dot_product(normal,p(i,:)-p(seed(1),:)))/magnitude
            if (distance > best) then
                best = distance
                seed(4) = i
            end if
        end do
        if (best <= orientation_tol) return
        inside = 0.25_dp*sum(p(seed,:),dim=1)

        capacity = 8*n+16
        allocate(faces(3,capacity),horizon(2,3*capacity), &
                 active(capacity),visible(capacity),chosen(n))
        active = .false.
        chosen = .false.
        chosen(seed) = .true.
        faces(:,1) = seed([1,2,3])
        faces(:,2) = seed([1,4,2])
        faces(:,3) = seed([2,4,3])
        faces(:,4) = seed([3,4,1])
        active(1:4) = .true.
        nf = 4
        do f=1,nf
            call orient_triangle(faces(:,f),p,inside)
        end do
        do i=1,n
            if (chosen(i)) cycle
            visible = .false.
            do f=1,nf
                if (.not. active(f)) cycle
                a = faces(1,f)
                normal = cross3(p(faces(2,f),:)-p(a,:),p(faces(3,f),:)-p(a,:))
                if (dot_product(normal,p(i,:)-p(a,:)) > orientation_tol*norm2(normal)) &
                    visible(f) = .true.
            end do
            if (.not. any(visible(1:nf))) cycle
            ne = 0
            do f=1,nf
                if (.not. visible(f)) cycle
                do t=1,3
                    a = faces(t,f)
                    b = faces(mod(t,3)+1,f)
                    idx = 0
                    do e=1,ne
                        if (horizon(1,e) == b .and. horizon(2,e) == a) then
                            idx = e
                            exit
                        end if
                    end do
                    if (idx > 0) then
                        horizon(:,idx) = 0
                    else
                        ne = ne+1
                        if (ne > size(horizon,2)) return
                        horizon(:,ne) = [a,b]
                    end if
                end do
            end do
            do f=1,nf
                if (visible(f)) active(f) = .false.
            end do
            h = 0
            do e=1,ne
                if (horizon(1,e) == 0) cycle
                h = h+1
                f = 0
                do g=1,nf
                    if (.not. active(g)) then
                        f = g
                        exit
                    end if
                end do
                if (f == 0) then
                    nf = nf+1
                    if (nf > capacity) return
                    f = nf
                end if
                faces(:,f) = [horizon(1,e),horizon(2,e),i]
                call orient_triangle(faces(:,f),p,inside)
                active(f) = .true.
            end do
            if (h < 3) return
        end do

        nfaces = count(active(1:nf))
        if (nfaces < 4) return
        allocate(triangles(3,nfaces))
        h = 0
        do f=1,nf
            if (.not. active(f)) cycle
            h = h + 1
            triangles(:,h) = faces(:,f)
        end do
        call assemble_solid(hull,p,original,scale,inside,triangles,ok)
    end subroutine build_solid

    ! Construct the public hull representation from an outward-oriented
    ! simplicial boundary. The same path is used for incremental and
    ! supporting-plane hulls, ensuring identical adjacency conventions.
    subroutine assemble_solid(hull,p,original,scale,inside,triangles,ok)
        type(convex_hull_type), intent(inout) :: hull !! hull to populate
        real(dp), intent(in) :: p(:,:) !! centered/scaled coordinates (n,3)
        real(dp), intent(in) :: original(:,:) !! original coordinates (n,3)
        real(dp), intent(in) :: scale !! isotropic coordinate scale
        real(dp), intent(in) :: inside(3) !! point inside the convex solid
        integer, intent(in) :: triangles(:,:) !! outward triangles (3,nfacets)
        logical, intent(out) :: ok !! true if manifold and positive volume
        logical, allocatable :: is_vertex(:)
        real(dp) :: normal(3), facearea, volume_sum
        integer :: n, nfaces, f, g, t, u, a, b, c, v, i

        ok = .false.
        n = size(p,1)
        nfaces = size(triangles,2)
        if (nfaces < 4) return
        allocate(hull%equations(4,nfaces),hull%neighbors(3,nfaces),is_vertex(n))
        is_vertex = .false.
        hull%neighbors = 0
        volume_sum = 0.0_dp
        do f=1,nfaces
            a = triangles(1,f)
            b = triangles(2,f)
            c = triangles(3,f)
            if (minval(triangles(:,f)) < 1 .or. maxval(triangles(:,f)) > n) return
            if (a==b .or. b==c .or. a==c) return
            is_vertex(triangles(:,f)) = .true.
            normal = cross3(p(b,:)-p(a,:),p(c,:)-p(a,:))
            facearea = norm2(normal)
            if (facearea <= 32.0_dp*tiny(1.0_dp)) return
            if (dot_product(normal,p(a,:)-inside) <= 0.0_dp) return
            volume_sum = volume_sum + dot_product(normal,p(a,:)-inside)
            hull%area = hull%area + 0.5_dp*facearea*scale**2
            normal = normal/facearea
            hull%equations(1:3,f) = normal
            hull%equations(4,f) = -dot_product(normal,original(a,:))
        end do
        hull%volume = volume_sum*scale**3/6.0_dp
        if (.not. ieee_is_finite(hull%area) .or. .not. ieee_is_finite(hull%volume)) return
        if (hull%volume <= 0.0_dp) return
        ! Every oriented edge must meet exactly one oppositely oriented edge.
        do f=1,nfaces
            do t=1,3
                a = triangles(mod(t,3)+1,f)
                b = triangles(mod(t+1,3)+1,f)
                do g=1,nfaces
                    if (g == f) cycle
                    do u=1,3
                        if (triangles(mod(u,3)+1,g) == b .and. &
                            triangles(mod(u+1,3)+1,g) == a) then
                            if (hull%neighbors(t,f) /= 0) return
                            hull%neighbors(t,f) = g
                        end if
                    end do
                end do
                if (hull%neighbors(t,f) == 0) return
            end do
        end do
        allocate(hull%vertices(count(is_vertex)))
        v = 0
        do i=1,n
            if (.not. is_vertex(i)) cycle
            v = v+1
            hull%vertices(v) = i
        end do
        hull%simplices = triangles
        ok = .true.
    end subroutine assemble_solid

    ! A topologically closed triangulation is not necessarily the convex hull:
    ! all input sites must satisfy every outward supporting halfspace. Test
    ! in centered/scaled coordinates to avoid subtracting large world offsets.
    subroutine validate_solid(hull,p,ok)
        type(convex_hull_type), intent(in) :: hull !! candidate hull
        real(dp), intent(in) :: p(:,:) !! centered/scaled coordinates (n,3)
        logical, intent(out) :: ok !! true if complete convex polytope
        real(dp) :: normal(3), unit(3), origin(3), tol
        integer :: nfaces, f, a, b, c

        ok = .false.
        if (.not. allocated(hull%simplices)) return
        if (.not. allocated(hull%vertices)) return
        if (.not. allocated(hull%neighbors)) return
        nfaces = size(hull%simplices,2)
        if (nfaces /= 2*size(hull%vertices)-4) return
        tol = 512.0_dp*epsilon(1.0_dp)
        do f=1,nfaces
            a = hull%simplices(1,f)
            b = hull%simplices(2,f)
            c = hull%simplices(3,f)
            origin = p(a,:)
            normal = cross3(p(b,:)-origin,p(c,:)-origin)
            if (norm2(normal) <= tol*tol) return
            unit = normal/norm2(normal)
            if (maxval(matmul(p-spread(origin,1,size(p,1)),unit)) > tol) return
        end do
        ok = .true.
    end subroutine validate_solid


    ! A convex hull is independent of insertion order, but floating-point
    ! near-coplanarity is not. Repeat the fast construction with deterministic
    ! shuffled orders before invoking the slower supporting-plane fallback.
    ! All results are remapped to the caller's original point indices and
    ! validated against *every* original point, never against jittered data.
    subroutine retry_solid(hull,p,original,scale,ok)
        use, intrinsic :: iso_fortran_env, only : int64
        type(convex_hull_type), intent(inout) :: hull !! output hull
        real(dp), intent(in) :: p(:,:) !! centered/scaled coordinates (n,3)
        real(dp), intent(in) :: original(:,:) !! original coordinates (n,3)
        real(dp), intent(in) :: scale !! isotropic scale
        logical, intent(out) :: ok !! true if a complete hull is recovered
        type(convex_hull_type) :: candidate
        integer, allocatable :: order(:), triangles(:,:)
        real(dp) :: inside(3)
        integer(int64) :: state
        integer :: n, i, j, attempt, temp, f
        logical :: found

        ok = .false.
        n = size(p,1)
        allocate(order(n))
        inside = sum(p,dim=1)/real(n,dp)
        do attempt=1,24
            order = [(i,i=1,n)]
            state = 1_int64 + 7919_int64*int(attempt,int64) + int(n,int64)
            do i=n,2,-1
                state = mod(48271_int64*state,2147483647_int64)
                j = 1+int(mod(state,int(i,int64)))
                temp = order(i)
                order(i) = order(j)
                order(j) = temp
            end do
            candidate = convex_hull_type()
            call build_solid(candidate,p(order,:),original(order,:),scale,found)
            if (.not. found) cycle
            call validate_solid(candidate,p(order,:),found)
            if (.not. found) cycle
            allocate(triangles(3,size(candidate%simplices,2)))
            do f=1,size(candidate%simplices,2)
                triangles(:,f) = order(candidate%simplices(:,f))
            end do
            hull = convex_hull_type()
            call assemble_solid(hull,p,original,scale,inside,triangles,found)
            deallocate(triangles)
            if (.not. found) cycle
            call validate_solid(hull,p,found)
            if (.not. found) cycle
            ok = .true.
            return
        end do
    end subroutine retry_solid

    ! Robust but more expensive fallback. Every candidate triangle determines
    ! a plane; only globally supporting planes are retained. All coplanar
    ! sites are projected into that plane and the planar hull is triangulated.
    ! This handles multiple coplanar facets and difficult insertion horizons.
    ! O(n^4) worst case; only called if the fast incremental hull fails its
    ! manifold or geometric containment checks.
    subroutine build_solid_fallback(hull,p,original,scale,ok)
        type(convex_hull_type), intent(inout) :: hull !! populated fallback hull
        real(dp), intent(in) :: p(:,:) !! centered/scaled coordinates (n,3)
        real(dp), intent(in) :: original(:,:) !! original coordinates (n,3)
        real(dp), intent(in) :: scale !! isotropic coordinate scale
        logical, intent(out) :: ok !! true for complete fallback hull
        integer, allocatable :: triangles(:,:), ids(:)
        real(dp), allocatable :: plane_normals(:,:), plane_offsets(:), projected(:,:)
        type(convex_hull_type) :: face_hull
        real(dp) :: normal(3), unit(3), basis_u(3), basis_v(3), origin(3), inside(3)
        real(dp) :: signed, vmin, vmax, mag, offset, plane_tol, match_tol
        logical :: is_new, polygon_ok
        integer :: n, capacity, np, nt, nc, i, j, k, l, f, t, m, a, b, c

        ok = .false.
        n = size(p,1)
        capacity = 2*n + 8
        plane_tol = 1.0e-11_dp
        match_tol = 1.0e-10_dp
        allocate(triangles(3,capacity),plane_normals(3,capacity), &
                 plane_offsets(capacity),ids(n),projected(n,2))
        np = 0
        nt = 0
        inside = sum(p,dim=1)/real(n,dp)
        do i=1,n-2
            do j=i+1,n-1
                do k=j+1,n
                    normal = cross3(p(j,:)-p(i,:),p(k,:)-p(i,:))
                    mag = norm2(normal)
                    if (mag <= plane_tol) cycle
                    unit = normal/mag
                    vmin = huge(1.0_dp)
                    vmax = -huge(1.0_dp)
                    do l=1,n
                        signed = dot_product(unit,p(l,:)-p(i,:))
                        vmin = min(vmin,signed)
                        vmax = max(vmax,signed)
                        if (vmin < -plane_tol .and. vmax > plane_tol) exit
                    end do
                    if (vmin < -plane_tol .and. vmax > plane_tol) cycle
                    if (vmax > plane_tol) unit = -unit
                    offset = -dot_product(unit,p(i,:))
                    is_new = .true.
                    do f=1,np
                        if (dot_product(unit,plane_normals(:,f)) >= 1.0_dp-match_tol .and. &
                            abs(offset-plane_offsets(f)) <= match_tol) then
                            is_new = .false.
                            exit
                        end if
                    end do
                    if (.not. is_new) cycle
                    if (np >= capacity) return
                    ! Project coplanar sites onto an oriented local 2D basis.
                    nc = 0
                    origin = p(i,:)
                    basis_u = (p(j,:)-origin)/norm2(p(j,:)-origin)
                    basis_v = cross3(unit,basis_u)
                    do l=1,n
                        signed = dot_product(unit,p(l,:)-origin)
                        if (abs(signed) > plane_tol) cycle
                        nc = nc + 1
                        ids(nc) = l
                        projected(nc,1) = dot_product(p(l,:)-origin,basis_u)
                        projected(nc,2) = dot_product(p(l,:)-origin,basis_v)
                    end do
                    if (nc < 3) cycle
                    face_hull = convex_hull_type()
                    call build_planar(face_hull,projected(:nc,:),projected(:nc,:),1.0_dp, &
                                      polygon_ok,tolerance=1.0e-12_dp)
                    if (.not. polygon_ok) cycle
                    m = size(face_hull%vertices)
                    if (nt + m-2 > capacity) return
                    a = ids(face_hull%vertices(1))
                    do t=2,m-1
                        b = ids(face_hull%vertices(t))
                        c = ids(face_hull%vertices(t+1))
                        nt = nt + 1
                        triangles(:,nt) = [a,b,c]
                    end do
                    np = np+1
                    plane_normals(:,np) = unit
                    plane_offsets(np) = offset
                end do
            end do
        end do
        if (nt < 4) return
        call assemble_solid(hull,p,original,scale,inside,triangles(:,:nt),ok)
    end subroutine build_solid_fallback
end module scifort_spatial_convex_hull
