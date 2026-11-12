! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bounded 2D/3D intersections of closed affine halfspaces. Candidate vertices
! are computed directly from plane intersections and tested against every
! inequality. Polar faces are triangulated via the native planar hull solver.
! No Qhull, BLAS/LAPACK, or external geometry library is required.
module scifort_spatial_halfspace
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : spatial_distance_success, spatial_distance_invalid_input
    use scifort_spatial_convex_hull, only : convex_hull, convex_hull_type
    implicit none
    private

    type, public :: halfspace_intersection_type
        real(dp), allocatable :: halfspaces(:,:)     !! input inequalities (n,d+1), A*x+b<=0
        real(dp), allocatable :: interior_point(:)   !! strictly feasible point (d)
        real(dp), allocatable :: intersections(:,:)  !! unique primal vertices (nv,d)
        real(dp), allocatable :: dual_points(:,:)    !! polar points -A/(A*x0+b) (n,d)
        integer, allocatable :: dual_facets(:,:)     !! simplicial facets (d,nf), 1-based
        integer, allocatable :: dual_vertices(:)    !! dual-hull vertex indices, 1-based
        real(dp), allocatable :: dual_equations(:,:) !! outward dual planes (d+1,nf)
        real(dp) :: dual_area = 0.0_dp               !! 2D perimeter or 3D surface area
        real(dp) :: dual_volume = 0.0_dp             !! 2D enclosed area or 3D volume
        integer :: status = spatial_distance_invalid_input
    contains
        procedure :: valid => halfspace_valid
    end type halfspace_intersection_type

    public :: halfspace_intersection
contains
    pure logical function halfspace_valid(self) result(ok)
        class(halfspace_intersection_type), intent(in) :: self !! polytope to inspect
        ok = self%status == spatial_distance_success .and. &
             allocated(self%intersections) .and. allocated(self%dual_points) .and. &
             allocated(self%dual_facets) .and. allocated(self%dual_vertices)
    end function halfspace_valid

    pure function cross3(a,b) result(c)
        real(dp), intent(in) :: a(3) !! first 3D vector
        real(dp), intent(in) :: b(3) !! second 3D vector
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
    end function cross3

    pure logical function recession_ray(q,direction) result(recedes)
        real(dp), intent(in) :: q(:,:) !! scaled polar constraints (n,d)
        real(dp), intent(in) :: direction(:) !! trial ray (d)
        real(dp) :: mag
        recedes = .false.
        mag = norm2(direction)
        if (mag < 128.0_dp*epsilon(1.0_dp)) return
        recedes = maxval(matmul(q,direction/mag)) <= 128.0_dp*epsilon(1.0_dp)
    end function recession_ray

    ! A bounded nonempty 2D/3D polytope has no recession direction.
    ! The extreme rays of a polyhedral cone in 3D occur on the crossing
    ! of two boundary planes (or within a nontrivial common nullspace).
    pure logical function is_unbounded(q) result(unbounded)
        real(dp), intent(in) :: q(:,:) !! scaled polar points (n,d)
        real(dp) :: ray(3)
        integer :: i,j,k,d,n
        d = size(q,2)
        n = size(q,1)
        unbounded = .false.
        do k=1,d
            ray = 0.0_dp
            ray(k) = 1.0_dp
            if (recession_ray(q,ray(:d)) .or. recession_ray(q,-ray(:d))) then
                unbounded = .true.
                return
            end if
        end do
        if (d==2) then
            do i=1,n
                ray(1:2) = [-q(i,2),q(i,1)]
                if (recession_ray(q,ray(:2)) .or. recession_ray(q,-ray(:2))) then
                    unbounded = .true.
                    return
                end if
            end do
        else
            do i=1,n-1
                do j=i+1,n
                    ray = cross3(q(i,:),q(j,:))
                    if (recession_ray(q,ray) .or. recession_ray(q,-ray)) then
                        unbounded = .true.
                        return
                    end if
                end do
            end do
        end if
    end function is_unbounded

    subroutine add_vertex(q,z,vertices,nvertices)
        real(dp), intent(in) :: q(:,:) !! polar points (n,d), q*x<=1
        real(dp), intent(in) :: z(:) !! candidate normalized primal vertex (d)
        real(dp), intent(inout) :: vertices(:,:) !! accepted vertices (capacity,d)
        integer, intent(inout) :: nvertices !! number of accepted vertices
        real(dp) :: tol
        integer :: j
        if (.not. all(ieee_is_finite(z))) return
        tol = 1.0e-8_dp * max(1.0_dp,maxval(abs(z)))
        if (maxval(matmul(q,z)) > 1.0_dp+tol) return
        do j=1,nvertices
            if (maxval(abs(vertices(j,:)-z)) <= tol) return
        end do
        if (nvertices >= size(vertices,1)) return
        nvertices = nvertices+1
        vertices(nvertices,:) = z
    end subroutine add_vertex

    ! Compute vertices of a bounded, full-dimensional polytope. Halfspaces
    ! follow SciPy's row convention sum(a*x)+b<=0. The supplied point must
    ! satisfy every inequality strictly. Status=1 denotes invalid, unbounded,
    ! or numerically degenerate/unsupported geometry.
    function halfspace_intersection(halfspaces, interior_point, status) result(poly)
        real(dp), intent(in) :: halfspaces(:,:) !! inequalities (n,d+1), A*x+b<=0
        real(dp), intent(in) :: interior_point(:) !! strictly feasible point (d)
        integer, intent(out), optional :: status !! 0=success, 1=invalid geometry
        type(halfspace_intersection_type) :: poly
        type(convex_hull_type) :: face_hull, dual_hull
        real(dp), allocatable :: q(:,:), original_dual(:,:), unit(:,:), &
                                  slacks(:), vertices(:,:), projected(:,:), &
                                  equations(:,:)
        integer, allocatable :: faces(:,:), active(:), order(:)
        logical, allocatable :: used(:)
        real(dp) :: scale, slack, length, det, det_tol, c(3), z(3), &
                    qnormal(3), direction(3), u(3), v(3), t, &
                    min_t, max_t, tol, area, volume, edge(3), normal_length
        integer :: n, d, i, j, k, l, nf, nv, na, ia, ib, nf_old, nface, &
                   idx_min, idx_max, code

        if (present(status)) status = spatial_distance_invalid_input
        d = size(interior_point)
        n = size(halfspaces,1)
        if (d < 2 .or. d > 3 .or. n < d+1) return
        if (size(halfspaces,2) /= d+1) return
        if (.not. all(ieee_is_finite(halfspaces))) return
        if (.not. all(ieee_is_finite(interior_point))) return

        allocate(unit(n,d), slacks(n), q(n,d), original_dual(n,d))
        unit = 0.0_dp
        slacks = huge(1.0_dp)
        original_dual = 0.0_dp
        scale = huge(1.0_dp)
        do i=1,n
            length = norm2(halfspaces(i,1:d))
            if (.not. ieee_is_finite(length)) return
            if (length == 0.0_dp) then
                if (halfspaces(i,d+1) >= 0.0_dp) return
                cycle
            end if
            unit(i,:) = halfspaces(i,1:d)/length
            slack = -(dot_product(unit(i,:),interior_point) + halfspaces(i,d+1)/length)
            if (.not. ieee_is_finite(slack) .or. slack <= 0.0_dp) return
            slacks(i) = slack
            scale = min(scale,slack)
            original_dual(i,:) = unit(i,:)/slack
        end do
        if (scale <= 0.0_dp .or. scale == huge(1.0_dp)) return
        if (.not. all(ieee_is_finite(original_dual))) return
        do i=1,n
            if (slacks(i)==huge(1.0_dp)) then
                q(i,:) = 0.0_dp
            else
                q(i,:) = unit(i,:)*(scale/slacks(i))
            end if
        end do
        if (is_unbounded(q)) return

        ! In 3D a bounded polyhedron of n halfspaces has at most 2*n-4
        ! vertices. Extra capacity tolerates intermediate roundoff splits.
        allocate(vertices(4*n,d))
        nv = 0
        if (d==2) then
            do i=1,n-1
                do j=i+1,n
                    det = q(i,1)*q(j,2)-q(i,2)*q(j,1)
                    det_tol = 64.0_dp*epsilon(1.0_dp)*norm2(q(i,:))*norm2(q(j,:))
                    if (abs(det) <= det_tol) cycle
                    z(1:2) = [q(j,2)-q(i,2),q(i,1)-q(j,1)]/det
                    call add_vertex(q,z(:2),vertices,nv)
                end do
            end do
        else
            do i=1,n-2
                do j=i+1,n-1
                    do k=j+1,n
                        c = cross3(q(j,:),q(k,:))
                        det = dot_product(q(i,:),c)
                        det_tol = 64.0_dp*epsilon(1.0_dp)* &
                                  norm2(q(i,:))*norm2(q(j,:))*norm2(q(k,:))
                        if (abs(det) <= det_tol) cycle
                        z = (c+cross3(q(k,:),q(i,:))+cross3(q(i,:),q(j,:)))/det
                        call add_vertex(q,z,vertices,nv)
                    end do
                end do
            end do
        end if
        if (nv < d+1) return

        ! Each primal vertex determines one face of the polar hull:
        ! q_i*z = 1 for all active (supporting) constraints. Triangulate
        ! coplanar dual polygons in their tangent plane for 3D.
        allocate(faces(d,4*n), equations(d+1,4*n), used(n),active(n))
        faces = 0
        used = .false.
        nf = 0
        area = 0.0_dp
        volume = 0.0_dp
        do l=1,nv
            z(:d) = vertices(l,:)
            normal_length = norm2(z(:d))
            if (normal_length <= 0.0_dp) return
            qnormal(:d) = z(:d)/normal_length
            tol = 1.0e-8_dp*max(1.0_dp,normal_length)
            na = 0
            do i=1,n
                if (abs(dot_product(q(i,:),z(:d))-1.0_dp) > tol) cycle
                na = na+1
                active(na) = i
            end do
            if (na < d) return
            nf_old = nf
            if (d==2) then
                direction(:2) = [-qnormal(2),qnormal(1)]
                min_t = huge(1.0_dp)
                max_t = -huge(1.0_dp)
                idx_min = 0
                idx_max = 0
                do i=1,na
                    t = dot_product(q(active(i),:),direction(:2))
                    if (t<min_t) then
                        min_t = t
                        idx_min = active(i)
                    end if
                    if (t>max_t) then
                        max_t = t
                        idx_max = active(i)
                    end if
                end do
                if (idx_min==idx_max) return
                nf = nf+1
                faces(:,nf) = [idx_min,idx_max]
                area = area+norm2(q(idx_max,:)-q(idx_min,:))
                volume = volume+0.5_dp*(q(idx_min,1)*q(idx_max,2)- &
                                        q(idx_min,2)*q(idx_max,1))
            else
                ! Form right-handed tangent basis (u,v,qnormal).
                direction = [1.0_dp,0.0_dp,0.0_dp]
                if (abs(qnormal(1))>0.8_dp) direction = [0.0_dp,1.0_dp,0.0_dp]
                u = cross3(direction,qnormal)
                u = u/norm2(u)
                v = cross3(qnormal,u)
                allocate(projected(na,2))
                do i=1,na
                    projected(i,1) = dot_product(q(active(i),:),u)
                    projected(i,2) = dot_product(q(active(i),:),v)
                end do
                face_hull = convex_hull(projected,status=code)
                deallocate(projected)
                if (code/=spatial_distance_success .or. .not. face_hull%valid()) return
                order = active(face_hull%vertices)
                nface = size(order)
                if (nface<3 .or. nf+nface-2>size(faces,2)) return
                do i=2,nface-1
                    nf = nf+1
                    faces(:,nf) = [order(1),order(i),order(i+1)]
                    ia = faces(1,nf)
                    ib = faces(2,nf)
                    k = faces(3,nf)
                    edge = cross3(q(ib,:)-q(ia,:),q(k,:)-q(ia,:))
                    area = area+0.5_dp*norm2(edge)
                    volume = volume+dot_product(q(ia,:),cross3(q(ib,:),q(k,:)))/6.0_dp
                end do
                deallocate(order)
            end if
            do i=nf_old+1,nf
                used(faces(:,i)) = .true.
            end do
            equations(1:d,nf_old+1:nf) = spread(qnormal(:d),dim=2,ncopies=nf-nf_old)
            equations(d+1,nf_old+1:nf) = -1.0_dp/normal_length/scale
        end do
        if (nf < d+1 .or. volume<=0.0_dp) return

        poly%halfspaces = halfspaces
        poly%interior_point = interior_point
        poly%intersections = spread(interior_point,dim=1,ncopies=nv) + &
                             scale*vertices(:nv,:)
        if (.not. all(ieee_is_finite(poly%intersections))) return
        poly%dual_points = original_dual
        poly%dual_facets = faces(:,1:nf)
        poly%dual_equations = equations(:,1:nf)
        if (d==2) then
            dual_hull = convex_hull(q,status=code)
            if (code/=spatial_distance_success) return
            poly%dual_vertices = dual_hull%vertices
        else
            poly%dual_vertices = pack([(i,i=1,n)],used)
        end if
        poly%dual_area = area/scale**(d-1)
        poly%dual_volume = volume/scale**d
        if (.not. ieee_is_finite(poly%dual_area) .or. &
            .not. ieee_is_finite(poly%dual_volume)) return
        poly%status = spatial_distance_success
        if (present(status)) status = spatial_distance_success
    end function halfspace_intersection
end module scifort_spatial_halfspace
