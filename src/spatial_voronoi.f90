! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Spherical Voronoi diagrams on circles (2D) and spheres (3D).
! Original Fortran convex-hull construction; no external geometry library.
module scifort_spatial_voronoi
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : spatial_distance_success, spatial_distance_invalid_input
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: twopi = 2.0_dp*pi

    type, public :: spherical_voronoi_type
        real(dp), allocatable :: points(:,:)       !! supplied generator coordinates (n,d)
        real(dp), allocatable :: vertices(:,:)     !! circumcenters on sphere (m,d)
        real(dp), allocatable :: center(:)         !! sphere center, size d
        real(dp) :: radius = 1.0_dp                !! sphere radius (>0)
        integer, allocatable :: simplices(:,:)    !! outward convex hull facets (d,m), 1-based
        integer, allocatable :: region_offsets(:) !! CSR region starts; n+1 entries, 1-based
        integer, allocatable :: region_vertices(:)!! indices into vertices, 1-based
        integer :: status = spatial_distance_invalid_input
    contains
        procedure :: valid => voronoi_valid
        procedure :: region => voronoi_region
        procedure :: calculate_areas => voronoi_areas
        procedure :: region_measures => voronoi_region_measures
    end type spherical_voronoi_type

    public :: spherical_voronoi
contains
    pure function cross3(a,b) result(c)
        real(dp), intent(in) :: a(3) !! first vector
        real(dp), intent(in) :: b(3) !! second vector
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
    end function cross3

    pure logical function voronoi_valid(self)
        class(spherical_voronoi_type), intent(in) :: self !! diagram
        voronoi_valid = self%status == spatial_distance_success .and. &
            allocated(self%vertices) .and. allocated(self%region_offsets) .and. &
            allocated(self%region_vertices)
    end function voronoi_valid

    function voronoi_region(self, generator) result(indices)
        class(spherical_voronoi_type), intent(in) :: self !! diagram
        integer, intent(in) :: generator !! 1-based generator number
        integer, allocatable :: indices(:)
        if (.not. self%valid()) then
            allocate(indices(0))
        else if (generator < 1 .or. generator > size(self%points,1)) then
            allocate(indices(0))
        else
            indices = self%region_vertices(self%region_offsets(generator):self%region_offsets(generator+1)-1)
        end if
    end function voronoi_region

    ! Returns an object with status=1 and no vertices if validation or hull construction fails.
    function spherical_voronoi(points, radius, center, threshold, status) result(diagram)
        real(dp), intent(in) :: points(:,:) !! generators (n,d), d=2 or 3
        real(dp), intent(in), optional :: radius !! sphere radius, default 1
        real(dp), intent(in), optional :: center(:) !! center (d), default zeros
        real(dp), intent(in), optional :: threshold !! relative duplicate/on-sphere tolerance, default 1e-6
        integer, intent(out), optional :: status !! zero success; one invalid geometry
        type(spherical_voronoi_type) :: diagram
        real(dp), allocatable :: p(:,:)
        real(dp) :: tolerance, r
        integer :: n, d, i, j
        logical :: good

        n = size(points,1)
        d = size(points,2)
        if (present(status)) status = spatial_distance_invalid_input
        if (d < 2 .or. d > 3 .or. n < d+1) return
        if (.not. all(ieee_is_finite(points))) return
        r = 1.0_dp
        if (present(radius)) r = radius
        if (.not. ieee_is_finite(r) .or. r <= 0.0_dp) return
        tolerance = 1.0e-6_dp
        if (present(threshold)) tolerance = threshold
        if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp) return
        allocate(diagram%center(d))
        diagram%center = 0.0_dp
        if (present(center)) then
            if (size(center) /= d) return
            if (.not. all(ieee_is_finite(center))) return
            diagram%center = center
        end if
        allocate(p(n,d))
        do j=1,d
            p(:,j) = (points(:,j)-diagram%center(j))/r
        end do
        if (.not. all(ieee_is_finite(p))) return
        do i=1,n
            if (abs(norm2(p(i,:))-1.0_dp) >= tolerance) return
            do j=1,i-1
                if (norm2(p(i,:)-p(j,:)) < tolerance) return
            end do
        end do
        diagram%radius = r
        good = .false.
        select case(d)
        case(2)
            call build_circle(diagram,p,good)
        case(3)
            call build_sphere(diagram,p,good)
        end select
        if (.not. good) then
            if (allocated(diagram%vertices)) deallocate(diagram%vertices)
            if (allocated(diagram%region_offsets)) deallocate(diagram%region_offsets)
            if (allocated(diagram%region_vertices)) deallocate(diagram%region_vertices)
            if (allocated(diagram%simplices)) deallocate(diagram%simplices)
            return
        end if
        diagram%points = points
        diagram%status = spatial_distance_success
        if (present(status)) status = spatial_distance_success
    end function spherical_voronoi

    subroutine build_circle(diagram, p, ok)
        type(spherical_voronoi_type), intent(inout) :: diagram !! circle diagram to populate
        real(dp), intent(in) :: p(:,:) !! normalized generators (n,2)
        logical, intent(out) :: ok !! valid diagram produced
        real(dp), allocatable :: angle(:)
        integer, allocatable :: order(:), inverse(:)
        real(dp) :: gap, midpoint
        integer :: n, i, j, tmp, current, following, previous

        n = size(p,1)
        ok = .false.
        allocate(order(n),inverse(n),angle(n))
        do i=1,n
            order(i) = i
            angle(i) = modulo(atan2(p(i,2),p(i,1)),twopi)
        end do
        ! Stable insertion sort, sufficient for moderately sized circular diagrams.
        do i=2,n
            tmp = order(i)
            j = i-1
            do while(j >= 1)
                if (angle(order(j)) <= angle(tmp)) exit
                order(j+1) = order(j)
                j = j-1
            end do
            order(j+1) = tmp
        end do
        allocate(diagram%vertices(n,2),diagram%simplices(2,n))
        do i=1,n
            current = order(i)
            following = order(mod(i,n)+1)
            gap = modulo(angle(following)-angle(current),twopi)
            if (gap < 1.0e-14_dp) return
            midpoint = angle(current)+0.5_dp*gap
            diagram%vertices(i,:) = diagram%center + &
                diagram%radius*[cos(midpoint),sin(midpoint)]
            diagram%simplices(:,i) = [current,following]
            inverse(current) = i
        end do
        allocate(diagram%region_offsets(n+1),diagram%region_vertices(2*n))
        do i=1,n
            diagram%region_offsets(i) = 2*(i-1)+1
            current = inverse(i)
            previous = mod(current+n-2,n)+1
            diagram%region_vertices(2*i-1:2*i) = [previous,current]
        end do
        diagram%region_offsets(n+1) = 2*n+1
        ok = .true.
    end subroutine build_circle

    ! Incremental triangular convex hull. All input sites lie on the same
    ! sphere; facets are stored in consistent outward orientation. The seed
    ! tetrahedron centroid is strictly inside every subsequent hull.
    subroutine build_sphere(diagram, p, ok)
        type(spherical_voronoi_type), intent(inout) :: diagram !! diagram being constructed
        real(dp), intent(in) :: p(:,:) !! normalized generators (n,3)
        logical, intent(out) :: ok !! true if closed triangular hull was constructed
        integer, allocatable :: faces(:,:), hull(:,:), edges(:,:), neighbor(:,:), offsets(:), region(:)
        logical, allocatable :: alive(:), seed(:), visible(:)
        real(dp) :: centroid(3), normal(3), distance, best
        integer :: n, capacity, ids(4), i, k, f, g, nf, nh, e, ne, a, b, t, m, pos
        integer :: first, nextface, site, nincident, index, h

        ok = .false.
        n = size(p,1)
        ids(1) = 1
        ! A long edge, a wide triangle, and a high-volume fourth vertex
        ! reduce sensitivity to the input order and near-degenerate seeds.
        best = -1.0_dp
        do i=2,n
            distance = sum((p(i,:)-p(ids(1),:))**2)
            if (distance > best) then
                best = distance
                ids(2) = i
            end if
        end do
        best = -1.0_dp
        do i=1,n
            if (i == ids(1) .or. i == ids(2)) cycle
            distance = norm2(cross3(p(ids(2),:)-p(ids(1),:),p(i,:)-p(ids(1),:)))
            if (distance > best) then
                best = distance
                ids(3) = i
            end if
        end do
        if (best < 1.0e-12_dp) return
        normal = cross3(p(ids(2),:)-p(ids(1),:),p(ids(3),:)-p(ids(1),:))
        best = -1.0_dp
        do i=1,n
            if (any(ids(1:3)==i)) cycle
            distance = abs(dot_product(normal,p(i,:)-p(ids(1),:)))/norm2(normal)
            if (distance > best) then
                best = distance
                ids(4) = i
            end if
        end do
        if (best < 1.0e-10_dp) return
        centroid = 0.25_dp*sum(p(ids,:),dim=1)
        capacity = 8*n+16
        allocate(faces(3,capacity),alive(capacity),seed(n),visible(capacity),edges(2,3*capacity))
        alive = .false.
        seed = .false.
        seed(ids) = .true.
        faces(:,1) = ids([1,2,3])
        faces(:,2) = ids([1,4,2])
        faces(:,3) = ids([2,4,3])
        faces(:,4) = ids([3,4,1])
        alive(1:4) = .true.
        nf = 4
        do f=1,4
            call orient_face(faces(:,f),p,centroid)
        end do

        do i=1,n
            if (seed(i)) cycle
            visible = .false.
            do f=1,nf
                if (.not. alive(f)) cycle
                a = faces(1,f)
                normal = cross3(p(faces(2,f),:)-p(a,:),p(faces(3,f),:)-p(a,:))
                distance = dot_product(normal,p(i,:)-p(a,:))
                if (distance > 128.0_dp*epsilon(1.0_dp)*norm2(normal)) visible(f) = .true.
            end do
            if (.not. any(visible)) cycle
            ! Each visible face contributes oriented edges; internal edges
            ! cancel in opposite pairs, leaving the horizon boundary.
            ne = 0
            do f=1,nf
                if (.not. visible(f)) cycle
                do t=1,3
                    a = faces(t,f)
                    b = faces(mod(t,3)+1,f)
                    index = 0
                    do e=1,ne
                        if (edges(1,e) == b .and. edges(2,e) == a) then
                            index = e
                            exit
                        end if
                    end do
                    if (index > 0) then
                        edges(:,index) = 0
                    else
                        ne = ne+1
                        edges(:,ne) = [a,b]
                    end if
                end do
            end do
            do f=1,nf
                if (visible(f)) alive(f) = .false.
            end do
            nh = 0
            do e=1,ne
                if (edges(1,e) == 0) cycle
                nh = nh+1
                ! Reuse dead facets to avoid excessive temporary growth.
                f = 0
                do g=1,nf
                    if (.not. alive(g)) then
                        f = g
                        exit
                    end if
                end do
                if (f == 0) then
                    nf = nf+1
                    if (nf > capacity) return
                    f = nf
                end if
                faces(:,f) = [edges(1,e),edges(2,e),i]
                call orient_face(faces(:,f),p,centroid)
                alive(f) = .true.
            end do
            if (nh < 3) return
        end do

        m = count(alive(1:nf))
        if (m < 4) return
        allocate(hull(3,m),diagram%vertices(m,3))
        k = 0
        do f=1,nf
            if (.not. alive(f)) cycle
            k = k+1
            hull(:,k) = faces(:,f)
            a = faces(1,f)
            normal = cross3(p(faces(2,f),:)-p(a,:),p(faces(3,f),:)-p(a,:))
            distance = norm2(normal)
            if (distance < 1.0e-14_dp) return
            diagram%vertices(k,:) = diagram%center + diagram%radius*(normal/distance)
        end do
        diagram%simplices = hull
        ! A triangle edge is shared by exactly two hull facets, in opposite
        ! orientations. These links give a numerically stable cyclic ordering
        ! of the dual vertices around each generator.
        allocate(neighbor(3,m))
        neighbor = 0
        do f=1,m
            do t=1,3
                a = hull(t,f)
                b = hull(mod(t,3)+1,f)
                do g=1,m
                    if (g == f) cycle
                    do e=1,3
                        if (hull(e,g) == b .and. hull(mod(e,3)+1,g) == a) then
                            neighbor(t,f) = g
                            exit
                        end if
                    end do
                    if (neighbor(t,f) > 0) exit
                end do
                if (neighbor(t,f) == 0) return
            end do
        end do
        allocate(offsets(n+1),region(3*m))
        index = 1
        do site=1,n
            offsets(site) = index
            first = 0
            nincident = 0
            do f=1,m
                if (any(hull(:,f) == site)) then
                    nincident = nincident+1
                    if (first == 0) first = f
                end if
            end do
            if (nincident < 3) return
            f = first
            do h=1,nincident
                if (index > size(region)) return
                region(index) = f
                index = index+1
                pos = 0
                do t=1,3
                    if (hull(t,f) == site) pos = t
                end do
                if (pos == 0) return
                nextface = neighbor(pos,f)
                if (h < nincident .and. nextface == first) return
                f = nextface
            end do
            if (f /= first) return
        end do
        offsets(n+1) = index
        if (index /= size(region)+1) return
        diagram%region_offsets = offsets
        diagram%region_vertices = region
        ok = .true.
    end subroutine build_sphere

    pure subroutine orient_face(face, p, interior)
        integer, intent(inout) :: face(3) !! triangle indices, returned outward
        real(dp), intent(in) :: p(:,:) !! normalized sites
        real(dp), intent(in) :: interior(3) !! interior point
        real(dp) :: n(3)
        integer :: tmp
        n = cross3(p(face(2),:)-p(face(1),:),p(face(3),:)-p(face(1),:))
        if (dot_product(n,interior-p(face(1),:)) > 0.0_dp) then
            tmp = face(2)
            face(2) = face(3)
            face(3) = tmp
        end if
    end subroutine orient_face

    ! In 3D, spherical polygon area is the sum of oriented-ring triangle
    ! solid angles, each computed using atan2 for stable small angles.
    ! In 2D, match SciPy: return the circular arc associated with each hull
    ! edge (Voronoi vertex), not the individual generator-region arc.
    function voronoi_areas(self) result(areas)
        class(spherical_voronoi_type), intent(in) :: self !! finished diagram
        real(dp), allocatable :: areas(:)
        real(dp) :: a(3), b(3), c(3), den, det, alpha, beta
        integer :: n, d, i, j, k, first, last, v1, v2
        if (.not. self%valid()) then
            allocate(areas(0))
            return
        end if
        n = size(self%points,1)
        d = size(self%points,2)
        allocate(areas(n))
        areas = 0.0_dp
        if (d == 2) then
            do i=1,n
                v1 = self%simplices(1,i)
                v2 = self%simplices(2,i)
                alpha = atan2(self%points(v1,2)-self%center(2), &
                              self%points(v1,1)-self%center(1))
                beta = atan2(self%points(v2,2)-self%center(2), &
                             self%points(v2,1)-self%center(1))
                areas(i) = self%radius*modulo(beta-alpha,twopi)
            end do
        else
            do i=1,n
                a = (self%points(i,:)-self%center)/self%radius
                first = self%region_offsets(i)
                last = self%region_offsets(i+1)-1
                do j=first,last
                    k = j+1
                    if (j == last) k = first
                    v1 = self%region_vertices(j)
                    v2 = self%region_vertices(k)
                    b = (self%vertices(v1,:)-self%center)/self%radius
                    c = (self%vertices(v2,:)-self%center)/self%radius
                    det = abs(dot_product(a,cross3(b,c)))
                    den = 1.0_dp + dot_product(a,b) + dot_product(b,c) + dot_product(c,a)
                    areas(i) = areas(i) + 2.0_dp*atan2(det,den)*self%radius**2
                end do
            end do
        end if
    end function voronoi_areas
    ! Per-generator Voronoi region measures (2D arc length or 3D spherical area).
    ! In 2D, SciPy's calculate_areas instead returns the hull-edge arc lengths.
    function voronoi_region_measures(self) result(areas)
        class(spherical_voronoi_type), intent(in) :: self !! finished diagram
        real(dp), allocatable :: areas(:)
        real(dp) :: alpha, beta
        integer :: n, i, v1, v2
        if (.not. self%valid()) then
            allocate(areas(0))
            return
        end if
        if (size(self%points,2) == 3) then
            areas = self%calculate_areas()
            return
        end if
        n = size(self%points,1)
        allocate(areas(n))
        do i=1,n
            v1 = self%region_vertices(self%region_offsets(i))
            v2 = self%region_vertices(self%region_offsets(i)+1)
            alpha = atan2(self%vertices(v1,2)-self%center(2), &
                          self%vertices(v1,1)-self%center(1))
            beta = atan2(self%vertices(v2,2)-self%center(2), &
                         self%vertices(v2,1)-self%center(1))
            areas(i) = self%radius*modulo(beta-alpha,twopi)
        end do
    end function voronoi_region_measures
end module scifort_spatial_voronoi
