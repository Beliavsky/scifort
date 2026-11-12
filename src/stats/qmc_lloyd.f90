! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Lloyd post-processing follows the relaxation and stopping convention of
! SciPy 1.17.0 scipy/stats/_qmc.py.  The finite Voronoi vertices are computed
! by an original SciFort N-dimensional Bowyer-Watson Delaunay implementation
! rather than Qhull.  See CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.

module scifort_qmc_lloyd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter :: status_success = 0
    integer, parameter :: status_invalid_shape = 1
    integer, parameter :: status_invalid_parameter = 2
    integer, parameter :: status_out_of_bounds = 3

    public :: qmc_lloyd

contains

    subroutine qmc_lloyd(sample, status, tol, maxiter, n_iterations)
        real(dp), intent(inout) :: sample(:, :) !! unit-cube sample optimized in place, d >= 2
        integer, intent(out) :: status !! QMC-compatible status code
        real(dp), intent(in), optional :: tol !! L1-spacing change tolerance, default 1e-5
        integer, intent(in), optional :: maxiter !! maximum Lloyd iterations, default 10
        integer, intent(out), optional :: n_iterations !! iterations actually performed

        integer :: i
        integer :: iterations
        integer :: maxiter_local
        real(dp) :: decay
        real(dp) :: l1_new
        real(dp) :: l1_old
        real(dp) :: root
        real(dp) :: tol_local

        status = status_success
        iterations = 0
        if (present(n_iterations)) n_iterations = 0
        tol_local = 1.0e-5_dp
        if (present(tol)) tol_local = tol
        maxiter_local = 10
        if (present(maxiter)) maxiter_local = maxiter

        if (size(sample,1) < 1) then
            status = status_invalid_shape
            return
        end if
        if (size(sample,2) < 2) then
            status = status_invalid_parameter
            return
        end if
        if (.not. all(ieee_is_finite(sample))) then
            status = status_invalid_parameter
            return
        end if
        if (any(sample < 0.0_dp) .or. any(sample > 1.0_dp)) then
            status = status_out_of_bounds
            return
        end if
        if (size(sample,1) < size(sample,2) + 1) then
            status = status_invalid_parameter
            return
        end if
        if (.not. ieee_is_finite(tol_local) .or. tol_local < 0.0_dp .or. &
                maxiter_local < 0) then
            status = status_invalid_parameter
            return
        end if
        if (maxiter_local == 0) return

        root = -real(maxiter_local,dp) / log(0.1_dp)
        l1_old = minimum_cityblock_distance(sample)
        do i = 0, maxiter_local - 1
            decay = exp(-real(i,dp) / root) + 0.9_dp
            call lloyd_iteration(sample, decay, status)
            if (status /= status_success) return
            iterations = iterations + 1
            l1_new = minimum_cityblock_distance(sample)
            if (abs(l1_new - l1_old) < tol_local) exit
            l1_old = l1_new
        end do
        if (present(n_iterations)) n_iterations = iterations
    end subroutine qmc_lloyd

    subroutine lloyd_iteration(sample, decay, status)
        real(dp), intent(inout) :: sample(:, :) !! current unit-cube sites
        real(dp), intent(in) :: decay !! positive relaxation multiplier
        integer, intent(out) :: status !! QMC-compatible status code

        integer :: d
        integer :: i
        integer :: j
        integer :: n
        integer :: n_vertices
        integer, allocatable :: counts(:)
        integer, allocatable :: simplex(:, :)
        logical, allocatable :: seen(:, :)
        real(dp), allocatable :: candidate(:)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: sums(:, :)
        real(dp), allocatable :: unique_vertices(:, :)

        status = status_success
        n = size(sample,1)
        d = size(sample,2)
        call delaunay_simplices(sample, simplex, centers, status)
        if (status /= status_success) return
        if (size(simplex,1) == 0) then
            status = status_invalid_parameter
            return
        end if

        allocate(sums(n,d), counts(n))
        allocate(unique_vertices(size(simplex,1),d))
        allocate(seen(n,size(simplex,1)))
        allocate(candidate(d))
        sums = 0.0_dp
        counts = 0
        seen = .false.
        n_vertices = 0

        do i = 1, size(simplex,1)
            j = find_matching_vertex(unique_vertices, n_vertices, centers(i,:))
            if (j == 0) then
                n_vertices = n_vertices + 1
                unique_vertices(n_vertices,:) = centers(i,:)
                j = n_vertices
            end if
            call accumulate_simplex_vertex(simplex(i,:), centers(i,:), j, seen, sums, counts)
        end do

        do i = 1, n
            if (counts(i) == 0) cycle
            candidate = sample(i,:) + &
                (sums(i,:) / real(counts(i),dp) - sample(i,:)) * decay
            if (all(ieee_is_finite(candidate)) .and. all(candidate >= 0.0_dp) .and. &
                    all(candidate <= 1.0_dp)) then
                sample(i,:) = candidate
            end if
        end do
    end subroutine lloyd_iteration

    subroutine delaunay_simplices(sample, simplex_out, centers_out, status)
        real(dp), intent(in) :: sample(:, :) !! sites whose finite Delaunay simplices are required
        integer, allocatable, intent(out) :: simplex_out(:, :) !! site indices of final simplices
        real(dp), allocatable, intent(out) :: centers_out(:, :) !! circumcenter of each simplex
        integer, intent(out) :: status !! QMC-compatible status code

        integer :: boundary_count
        integer :: capacity
        integer :: d
        integer :: facet_count
        integer :: i
        integer :: j
        integer :: k
        integer :: keep
        integer :: n
        integer :: n_bad
        integer :: nsimp
        integer :: p
        integer, allocatable :: facet(:)
        integer, allocatable :: facets(:, :)
        integer, allocatable :: facet_multiplicity(:)
        integer, allocatable :: simplex(:, :)
        integer, allocatable :: new_simplex(:)
        logical, allocatable :: bad(:)
        real(dp) :: lscale
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: points(:, :)
        real(dp), allocatable :: radius_squared(:)

        status = status_success
        n = size(sample,1)
        d = size(sample,2)
        allocate(points(n+d+1,d))
        points(1:n,:) = sample
        lscale = 4.0_dp * real(d + 1,dp)
        points(n+1,:) = 0.5_dp - lscale
        do i = 1, d
            points(n+1+i,:) = 0.5_dp
            points(n+1+i,i) = 0.5_dp + lscale
        end do

        capacity = max(16,4*(d+1))
        allocate(simplex(capacity,d+1), centers(capacity,d), radius_squared(capacity))
        allocate(new_simplex(d+1), facet(d))
        nsimp = 0
        do i = 1, d + 1
            new_simplex(i) = n + i
        end do
        call append_simplex(points, new_simplex, simplex, centers, radius_squared, &
            nsimp, capacity, status)
        if (status /= status_success) return

        do p = 1, n
            allocate(bad(nsimp))
            bad = .false.
            n_bad = 0
            do i = 1, nsimp
                if (point_in_circumsphere(points(p,:), centers(i,:), radius_squared(i))) then
                    bad(i) = .true.
                    n_bad = n_bad + 1
                end if
            end do
            if (n_bad == 0) then
                deallocate(bad)
                status = status_invalid_parameter
                return
            end if

            allocate(facets(n_bad*(d+1),d))
            allocate(facet_multiplicity(n_bad*(d+1)))
            facet_count = 0
            facet_multiplicity = 0
            do i = 1, nsimp
                if (.not. bad(i)) cycle
                do j = 1, d + 1
                    call simplex_facet(simplex(i,:), j, facet)
                    call sort_integers(facet)
                    k = locate_facet(facets, facet_count, facet)
                    if (k == 0) then
                        facet_count = facet_count + 1
                        facets(facet_count,:) = facet
                        facet_multiplicity(facet_count) = 1
                    else
                        facet_multiplicity(k) = facet_multiplicity(k) + 1
                    end if
                end do
            end do

            keep = 0
            do i = 1, nsimp
                if (bad(i)) cycle
                keep = keep + 1
                if (keep /= i) then
                    simplex(keep,:) = simplex(i,:)
                    centers(keep,:) = centers(i,:)
                    radius_squared(keep) = radius_squared(i)
                end if
            end do
            nsimp = keep

            boundary_count = 0
            do i = 1, facet_count
                if (facet_multiplicity(i) /= 1) cycle
                boundary_count = boundary_count + 1
                new_simplex(1:d) = facets(i,:)
                new_simplex(d+1) = p
                call append_simplex(points, new_simplex, simplex, centers, radius_squared, &
                    nsimp, capacity, status)
                if (status /= status_success) then
                    deallocate(bad, facets, facet_multiplicity)
                    return
                end if
            end do
            deallocate(bad, facets, facet_multiplicity)
            if (boundary_count == 0) then
                status = status_invalid_parameter
                return
            end if
        end do

        keep = 0
        do i = 1, nsimp
            if (any(simplex(i,:) > n)) cycle
            keep = keep + 1
            if (keep /= i) then
                simplex(keep,:) = simplex(i,:)
                centers(keep,:) = centers(i,:)
            end if
        end do
        allocate(simplex_out(keep,d+1), centers_out(keep,d))
        if (keep > 0) then
            simplex_out = simplex(1:keep,:)
            centers_out = centers(1:keep,:)
        end if
    end subroutine delaunay_simplices

    subroutine append_simplex(points, indices, simplex, centers, radius_squared, &
            nsimp, capacity, status)
        real(dp), intent(in) :: points(:, :) !! all ordinary and super-simplex vertices
        integer, intent(in) :: indices(:) !! d+1 vertex indices for the new simplex
        integer, allocatable, intent(inout) :: simplex(:, :) !! simplex storage
        real(dp), allocatable, intent(inout) :: centers(:, :) !! circumcenter storage
        real(dp), allocatable, intent(inout) :: radius_squared(:) !! squared circumradii
        integer, intent(inout) :: nsimp !! number of stored simplices
        integer, intent(inout) :: capacity !! allocated simplex capacity
        integer, intent(out) :: status !! QMC-compatible status code

        integer :: info
        real(dp) :: center(size(points,2))
        real(dp) :: radius2

        call simplex_circumsphere(points(indices,:), center, radius2, info)
        if (info /= 0) then
            status = status_success
            return
        end if
        if (nsimp >= capacity) then
            call grow_simplex_storage(simplex, centers, radius_squared, capacity)
        end if
        nsimp = nsimp + 1
        simplex(nsimp,:) = indices
        centers(nsimp,:) = center
        radius_squared(nsimp) = radius2
        status = status_success
    end subroutine append_simplex

    subroutine grow_simplex_storage(simplex, centers, radius_squared, capacity)
        integer, allocatable, intent(inout) :: simplex(:, :) !! simplex storage to grow
        real(dp), allocatable, intent(inout) :: centers(:, :) !! circumcenter storage to grow
        real(dp), allocatable, intent(inout) :: radius_squared(:) !! radius storage to grow
        integer, intent(inout) :: capacity !! current capacity, doubled on return

        integer, allocatable :: simplex_new(:, :)
        real(dp), allocatable :: centers_new(:, :)
        real(dp), allocatable :: radius_new(:)
        integer :: old_capacity

        old_capacity = capacity
        capacity = max(2*capacity,1)
        allocate(simplex_new(capacity,size(simplex,2)))
        allocate(centers_new(capacity,size(centers,2)))
        allocate(radius_new(capacity))
        simplex_new = 0
        centers_new = 0.0_dp
        radius_new = 0.0_dp
        simplex_new(1:old_capacity,:) = simplex
        centers_new(1:old_capacity,:) = centers
        radius_new(1:old_capacity) = radius_squared
        call move_alloc(simplex_new,simplex)
        call move_alloc(centers_new,centers)
        call move_alloc(radius_new,radius_squared)
    end subroutine grow_simplex_storage

    subroutine simplex_circumsphere(vertices, center, radius_squared, info)
        real(dp), intent(in) :: vertices(:, :) !! d+1 simplex vertices in d dimensions
        real(dp), intent(out) :: center(:) !! circumcenter
        real(dp), intent(out) :: radius_squared !! squared circumradius
        integer, intent(out) :: info !! zero on success, one for singular/nonfinite geometry

        integer :: d
        integer :: i
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)

        d = size(vertices,2)
        info = 0
        center = 0.0_dp
        radius_squared = 0.0_dp
        if (size(vertices,1) /= d + 1 .or. size(center) /= d) then
            info = 1
            return
        end if
        allocate(a(d,d), b(d))
        do i = 1, d
            a(i,:) = 2.0_dp * (vertices(i+1,:) - vertices(1,:))
            b(i) = sum(vertices(i+1,:)**2) - sum(vertices(1,:)**2)
        end do
        call solve_linear_system(a, b, center, info)
        if (info /= 0 .or. .not. all(ieee_is_finite(center))) then
            info = 1
            return
        end if
        radius_squared = sum((center - vertices(1,:))**2)
        if (.not. ieee_is_finite(radius_squared) .or. radius_squared < 0.0_dp) info = 1
    end subroutine simplex_circumsphere

    subroutine solve_linear_system(a, b, x, info)
        real(dp), intent(in) :: a(:, :) !! square coefficient matrix
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! solution
        integer, intent(out) :: info !! zero on success, one when singular

        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot
        real(dp) :: factor
        real(dp) :: pivot_tol
        real(dp) :: tmp
        real(dp), allocatable :: aa(:, :)
        real(dp), allocatable :: bb(:)
        real(dp), allocatable :: row_tmp(:)

        n = size(b)
        info = 0
        x = 0.0_dp
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
            info = 1
            return
        end if
        allocate(aa(n,n), bb(n), row_tmp(n))
        aa = a
        bb = b
        pivot_tol = 256.0_dp * epsilon(1.0_dp) * max(1.0_dp,maxval(abs(aa)))
        do k = 1, n - 1
            pivot = k
            do i = k + 1, n
                if (abs(aa(i,k)) > abs(aa(pivot,k))) pivot = i
            end do
            if (abs(aa(pivot,k)) <= pivot_tol) then
                info = 1
                return
            end if
            if (pivot /= k) then
                row_tmp = aa(k,:)
                aa(k,:) = aa(pivot,:)
                aa(pivot,:) = row_tmp
                tmp = bb(k)
                bb(k) = bb(pivot)
                bb(pivot) = tmp
            end if
            do i = k + 1, n
                factor = aa(i,k) / aa(k,k)
                aa(i,k) = 0.0_dp
                do j = k + 1, n
                    aa(i,j) = aa(i,j) - factor * aa(k,j)
                end do
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do
        if (abs(aa(n,n)) <= pivot_tol) then
            info = 1
            return
        end if
        do i = n, 1, -1
            tmp = bb(i)
            do j = i + 1, n
                tmp = tmp - aa(i,j) * x(j)
            end do
            if (abs(aa(i,i)) <= pivot_tol) then
                info = 1
                return
            end if
            x(i) = tmp / aa(i,i)
        end do
    end subroutine solve_linear_system

    pure function point_in_circumsphere(point, center, radius_squared) result(inside)
        real(dp), intent(in) :: point(:) !! point tested against a simplex circumsphere
        real(dp), intent(in) :: center(:) !! circumsphere center
        real(dp), intent(in) :: radius_squared !! squared circumsphere radius
        logical :: inside

        real(dp) :: distance_squared
        real(dp) :: tolerance

        distance_squared = sum((point - center)**2)
        tolerance = 512.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp,abs(radius_squared),abs(distance_squared))
        inside = distance_squared <= radius_squared + tolerance
    end function point_in_circumsphere

    subroutine simplex_facet(simplex, omitted, facet)
        integer, intent(in) :: simplex(:) !! d+1 simplex vertex indices
        integer, intent(in) :: omitted !! one-based vertex omitted from the facet
        integer, intent(out) :: facet(:) !! remaining d vertex indices

        integer :: i
        integer :: k

        k = 0
        do i = 1, size(simplex)
            if (i == omitted) cycle
            k = k + 1
            facet(k) = simplex(i)
        end do
    end subroutine simplex_facet

    subroutine sort_integers(values)
        integer, intent(inout) :: values(:) !! integer keys sorted in ascending order

        integer :: i
        integer :: j
        integer :: key

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = key
        end do
    end subroutine sort_integers

    function locate_facet(facets, nfacets, target) result(index)
        integer, intent(in) :: facets(:, :) !! existing sorted facet keys
        integer, intent(in) :: nfacets !! number of populated facet rows
        integer, intent(in) :: target(:) !! sorted facet sought
        integer :: index

        integer :: i

        index = 0
        do i = 1, nfacets
            if (all(facets(i,:) == target)) then
                index = i
                return
            end if
        end do
    end function locate_facet

    function find_matching_vertex(vertices, n_vertices, target) result(index)
        real(dp), intent(in) :: vertices(:, :) !! existing unique Voronoi vertices
        integer, intent(in) :: n_vertices !! number of populated vertex rows
        real(dp), intent(in) :: target(:) !! candidate vertex
        integer :: index

        integer :: i
        real(dp) :: scale_value
        real(dp) :: tolerance

        index = 0
        scale_value = max(1.0_dp,maxval(abs(target)))
        tolerance = 1024.0_dp * epsilon(1.0_dp) * scale_value
        do i = 1, n_vertices
            if (maxval(abs(vertices(i,:) - target)) <= tolerance) then
                index = i
                return
            end if
        end do
    end function find_matching_vertex

    subroutine accumulate_simplex_vertex(simplex, center, vertex_id, seen, sums, counts)
        integer, intent(in) :: simplex(:) !! site indices touching this Voronoi vertex
        real(dp), intent(in) :: center(:) !! Voronoi vertex coordinates
        integer, intent(in) :: vertex_id !! deduplicated one-based Voronoi vertex id
        logical, intent(inout) :: seen(:, :) !! site-by-vertex incidence already accumulated
        real(dp), intent(inout) :: sums(:, :) !! accumulated finite vertices for each site
        integer, intent(inout) :: counts(:) !! unique finite-vertex count for each site

        integer :: i
        integer :: site

        do i = 1, size(simplex)
            site = simplex(i)
            if (site < 1 .or. site > size(sums,1)) cycle
            if (seen(site,vertex_id)) cycle
            seen(site,vertex_id) = .true.
            sums(site,:) = sums(site,:) + center
            counts(site) = counts(site) + 1
        end do
    end subroutine accumulate_simplex_vertex

    function minimum_cityblock_distance(sample) result(value)
        real(dp), intent(in) :: sample(:, :) !! sample whose minimum pairwise L1 spacing is returned
        real(dp) :: value

        integer :: i
        integer :: j

        value = huge(1.0_dp)
        do i = 1, size(sample,1) - 1
            do j = i + 1, size(sample,1)
                value = min(value,sum(abs(sample(i,:) - sample(j,:))))
            end do
        end do
    end function minimum_cityblock_distance

end module scifort_qmc_lloyd
