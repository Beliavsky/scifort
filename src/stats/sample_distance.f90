! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_sample_distance
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: wasserstein_distance
    public :: wasserstein_distance_nd
    public :: energy_distance

    interface wasserstein_distance_nd
        module procedure wasserstein_distance_nd_vector
        module procedure wasserstein_distance_nd_matrix
    end interface wasserstein_distance_nd

contains

    pure function wasserstein_distance(u_values, v_values, u_weights, v_weights) result(distance)
        real(dp), intent(in) :: u_values(:) !! first one-dimensional support/sample
        real(dp), intent(in) :: v_values(:) !! second one-dimensional support/sample
        real(dp), intent(in), optional :: u_weights(:) !! nonnegative first-sample weights
        real(dp), intent(in), optional :: v_weights(:) !! nonnegative second-sample weights
        real(dp) :: distance

        distance = cdf_distance(1, u_values, v_values, u_weights, v_weights)
    end function wasserstein_distance

    pure function wasserstein_distance_nd_vector(u_values, v_values, u_weights, v_weights) result(distance)
        real(dp), intent(in) :: u_values(:) !! first one-dimensional support/sample
        real(dp), intent(in) :: v_values(:) !! second one-dimensional support/sample
        real(dp), intent(in), optional :: u_weights(:) !! nonnegative first-sample weights
        real(dp), intent(in), optional :: v_weights(:) !! nonnegative second-sample weights
        real(dp) :: distance

        distance = wasserstein_distance(u_values, v_values, u_weights, v_weights)
    end function wasserstein_distance_nd_vector

    pure function wasserstein_distance_nd_matrix(u_values, v_values, u_weights, v_weights) result(distance)
        real(dp), intent(in) :: u_values(:, :) !! first multidimensional support/sample, observations by rows
        real(dp), intent(in) :: v_values(:, :) !! second multidimensional support/sample, observations by rows
        real(dp), intent(in), optional :: u_weights(:) !! nonnegative first-sample weights
        real(dp), intent(in), optional :: v_weights(:) !! nonnegative second-sample weights
        real(dp) :: distance

        real(dp), allocatable :: cost(:, :)
        real(dp), allocatable :: wu(:)
        real(dp), allocatable :: wv(:)
        real(dp) :: positive_inf
        logical :: u_has_inf
        logical :: v_has_inf
        integer :: i
        integer :: j

        if (size(u_values, 1) == 0 .or. size(v_values, 1) == 0) then
            error stop 'scifort: sample distance distributions must be nonempty', quiet=.true.
        end if
        if (size(u_values, 2) /= size(v_values, 2)) then
            error stop 'scifort: multidimensional samples must have equal dimension', quiet=.true.
        end if
        if (any(ieee_is_nan(u_values)) .or. any(ieee_is_nan(v_values))) then
            error stop 'scifort: multidimensional sample values must not contain NaN', quiet=.true.
        end if

        positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
        u_has_inf = any(.not. ieee_is_finite(u_values))
        v_has_inf = any(.not. ieee_is_finite(v_values))
        if (u_has_inf .neqv. v_has_inf) then
            distance = positive_inf
            return
        else if (u_has_inf .and. v_has_inf) then
            distance = quiet_nan(0.0_dp)
            return
        end if

        allocate(wu(size(u_values, 1)), wv(size(v_values, 1)))
        call prepare_weights(u_weights, size(wu), wu)
        call prepare_weights(v_weights, size(wv), wv)
        call force_unit_sum(wu)
        call force_unit_sum(wv)

        allocate(cost(size(wu), size(wv)))
        do j = 1, size(wv)
            do i = 1, size(wu)
                cost(i, j) = safe_euclidean_distance(u_values(i, :), v_values(j, :))
            end do
        end do
        distance = minimum_transport_cost(cost, wu, wv)
    end function wasserstein_distance_nd_matrix

    pure function energy_distance(u_values, v_values, u_weights, v_weights) result(distance)
        real(dp), intent(in) :: u_values(:) !! first one-dimensional support/sample
        real(dp), intent(in) :: v_values(:) !! second one-dimensional support/sample
        real(dp), intent(in), optional :: u_weights(:) !! nonnegative first-sample weights
        real(dp), intent(in), optional :: v_weights(:) !! nonnegative second-sample weights
        real(dp) :: distance
        real(dp) :: l2

        l2 = cdf_distance(2, u_values, v_values, u_weights, v_weights)
        distance = sqrt(2.0_dp) * l2
    end function energy_distance

    pure function cdf_distance(power, u_values, v_values, u_weights, v_weights) result(distance)
        integer, intent(in) :: power
        real(dp), intent(in) :: u_values(:)
        real(dp), intent(in) :: v_values(:)
        real(dp), intent(in), optional :: u_weights(:)
        real(dp), intent(in), optional :: v_weights(:)
        real(dp) :: distance

        real(dp), allocatable :: u(:)
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: wu(:)
        real(dp), allocatable :: wv(:)
        real(dp) :: cdf_u
        real(dp) :: cdf_v
        real(dp) :: current
        real(dp) :: delta
        real(dp) :: diff
        real(dp) :: next_value
        real(dp) :: total
        integer :: i
        integer :: j
        integer :: nneg_inf
        integer :: npos_inf

        if (size(u_values) == 0 .or. size(v_values) == 0) then
            error stop 'scifort: sample distance distributions must be nonempty', quiet=.true.
        end if

        if (any(ieee_is_nan(u_values)) .or. any(ieee_is_nan(v_values))) then
            distance = quiet_nan(0.0_dp)
            return
        end if

        nneg_inf = count_negative_infinity(u_values) + count_negative_infinity(v_values)
        npos_inf = count_positive_infinity(u_values) + count_positive_infinity(v_values)
        if (nneg_inf > 1 .or. npos_inf > 1) then
            distance = quiet_nan(0.0_dp)
            return
        end if

        allocate(u(size(u_values)), v(size(v_values)))
        allocate(wu(size(u_values)), wv(size(v_values)))
        u = u_values
        v = v_values
        call prepare_weights(u_weights, size(u), wu)
        call prepare_weights(v_weights, size(v), wv)
        call sort_pairs(u, wu)
        call sort_pairs(v, wv)

        i = 1
        j = 1
        cdf_u = 0.0_dp
        cdf_v = 0.0_dp
        total = 0.0_dp
        current = min(u(1), v(1))

        do
            do while (i <= size(u))
                if (.not. same_value(u(i), current)) exit
                cdf_u = cdf_u + wu(i)
                i = i + 1
            end do
            do while (j <= size(v))
                if (.not. same_value(v(j), current)) exit
                cdf_v = cdf_v + wv(j)
                j = j + 1
            end do

            if (i > size(u) .and. j > size(v)) exit
            if (i > size(u)) then
                next_value = v(j)
            else if (j > size(v)) then
                next_value = u(i)
            else
                next_value = min(u(i), v(j))
            end if

            diff = abs(cdf_u - cdf_v)
            if (.not. ieee_is_finite(current) .or. .not. ieee_is_finite(next_value)) then
                if (diff <= 0.0_dp) then
                    distance = quiet_nan(0.0_dp)
                    return
                end if
                distance = ieee_value(0.0_dp, ieee_positive_inf)
                return
            end if
            delta = next_value - current
            if (power == 1) then
                total = total + diff * delta
            else
                total = total + diff * diff * delta
            end if
            current = next_value
        end do

        if (power == 1) then
            distance = total
        else
            distance = sqrt(total)
        end if
    end function cdf_distance

    pure function minimum_transport_cost(cost, source_weights, target_weights) result(distance)
        real(dp), intent(in) :: cost(:, :)
        real(dp), intent(in) :: source_weights(:)
        real(dp), intent(in) :: target_weights(:)
        real(dp) :: distance

        real(dp), allocatable :: capacity(:, :)
        real(dp), allocatable :: edge_cost(:, :)
        real(dp), allocatable :: shortest(:)
        real(dp) :: augment
        real(dp) :: candidate
        real(dp) :: flow_tol
        real(dp) :: positive_inf
        real(dp) :: relax_tol
        real(dp) :: sent
        integer, allocatable :: parent(:)
        integer :: from_node
        integer :: i
        integer :: iteration
        integer :: j
        integer :: node
        integer :: nodes
        integer :: pass
        integer :: sink
        integer :: source
        integer :: to_node
        logical :: changed

        source = 1
        sink = size(source_weights) + size(target_weights) + 2
        nodes = sink
        positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
        allocate(capacity(nodes, nodes), edge_cost(nodes, nodes))
        allocate(shortest(nodes), parent(nodes))
        capacity = 0.0_dp
        edge_cost = 0.0_dp

        do i = 1, size(source_weights)
            capacity(source, 1 + i) = source_weights(i)
        end do
        do j = 1, size(target_weights)
            capacity(1 + size(source_weights) + j, sink) = target_weights(j)
        end do
        do j = 1, size(target_weights)
            to_node = 1 + size(source_weights) + j
            do i = 1, size(source_weights)
                from_node = 1 + i
                if (.not. ieee_is_finite(cost(i, j))) cycle
                capacity(from_node, to_node) = 1.0_dp
                edge_cost(from_node, to_node) = cost(i, j)
                edge_cost(to_node, from_node) = -cost(i, j)
            end do
        end do

        sent = 0.0_dp
        distance = 0.0_dp
        iteration = 0
        flow_tol = 64.0_dp * epsilon(1.0_dp)
        do while (sent < 1.0_dp - flow_tol)
            iteration = iteration + 1
            if (iteration > 1000000) then
                error stop 'scifort: multidimensional transport solver did not converge', quiet=.true.
            end if

            shortest = positive_inf
            shortest(source) = 0.0_dp
            parent = 0
            do pass = 1, nodes - 1
                changed = .false.
                do from_node = 1, nodes
                    if (.not. ieee_is_finite(shortest(from_node))) cycle
                    do to_node = 1, nodes
                        if (capacity(from_node, to_node) <= 0.0_dp) cycle
                        candidate = shortest(from_node) + edge_cost(from_node, to_node)
                        if (.not. ieee_is_finite(shortest(to_node))) then
                            shortest(to_node) = candidate
                            parent(to_node) = from_node
                            changed = .true.
                        else
                            relax_tol = 64.0_dp * epsilon(1.0_dp) * &
                                max(1.0_dp, abs(shortest(to_node)), abs(candidate))
                            if (candidate < shortest(to_node) - relax_tol) then
                                shortest(to_node) = candidate
                                parent(to_node) = from_node
                                changed = .true.
                            end if
                        end if
                    end do
                end do
                if (.not. changed) exit
            end do

            if (parent(sink) == 0) then
                if (1.0_dp - sent <= flow_tol) exit
                distance = positive_inf
                return
            end if

            augment = 1.0_dp - sent
            node = sink
            do while (node /= source)
                from_node = parent(node)
                if (from_node == 0) then
                    error stop 'scifort: invalid multidimensional transport path', quiet=.true.
                end if
                augment = min(augment, capacity(from_node, node))
                node = from_node
            end do
            if (augment <= 0.0_dp) then
                error stop 'scifort: nonpositive multidimensional transport augmentation', quiet=.true.
            end if

            node = sink
            do while (node /= source)
                from_node = parent(node)
                capacity(from_node, node) = capacity(from_node, node) - augment
                if (capacity(from_node, node) < 0.0_dp) capacity(from_node, node) = 0.0_dp
                capacity(node, from_node) = capacity(node, from_node) + augment
                node = from_node
            end do
            sent = sent + augment
            if (1.0_dp - sent <= flow_tol) sent = 1.0_dp
            distance = distance + augment * shortest(sink)
        end do

        if (distance < 0.0_dp) then
            relax_tol = 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(distance))
            if (distance >= -relax_tol) distance = 0.0_dp
        end if
    end function minimum_transport_cost

    pure function safe_euclidean_distance(a, b) result(distance)
        real(dp), intent(in) :: a(:)
        real(dp), intent(in) :: b(:)
        real(dp) :: distance

        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: diff
        real(dp) :: factor
        real(dp) :: positive_inf
        real(dp) :: ratio
        real(dp) :: scale
        real(dp) :: sumsq
        integer :: i

        positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
        scale = 0.0_dp
        sumsq = 1.0_dp
        do i = 1, size(a)
            aa = abs(a(i))
            bb = abs(b(i))
            if ((a(i) > 0.0_dp .and. b(i) < 0.0_dp) .or. &
                (a(i) < 0.0_dp .and. b(i) > 0.0_dp)) then
                if (aa > huge(0.0_dp) - bb) then
                    distance = positive_inf
                    return
                end if
                diff = aa + bb
            else
                diff = abs(a(i) - b(i))
            end if

            if (diff > 0.0_dp) then
                if (scale < diff) then
                    if (scale <= 0.0_dp) then
                        sumsq = 1.0_dp
                    else
                        ratio = scale / diff
                        sumsq = 1.0_dp + sumsq * ratio * ratio
                    end if
                    scale = diff
                else
                    ratio = diff / scale
                    sumsq = sumsq + ratio * ratio
                end if
            end if
        end do

        if (scale <= 0.0_dp) then
            distance = 0.0_dp
            return
        end if
        factor = sqrt(sumsq)
        if (scale > huge(0.0_dp) / factor) then
            distance = positive_inf
        else
            distance = scale * factor
        end if
    end function safe_euclidean_distance

    pure subroutine force_unit_sum(weights)
        real(dp), intent(inout) :: weights(:)

        real(dp) :: residual
        integer :: i
        integer :: imax

        if (size(weights) == 1) then
            weights(1) = 1.0_dp
            return
        end if
        imax = 1
        do i = 2, size(weights)
            if (weights(i) > weights(imax)) imax = i
        end do
        residual = 1.0_dp - sum(weights)
        weights(imax) = weights(imax) + residual
    end subroutine force_unit_sum

    pure subroutine prepare_weights(input_weights, n, weights)
        real(dp), intent(in), optional :: input_weights(:)
        integer, intent(in) :: n
        real(dp), intent(out) :: weights(n)

        real(dp) :: total
        integer :: i

        if (.not. present(input_weights)) then
            weights = 1.0_dp / real(n, dp)
            return
        end if
        if (size(input_weights) /= n) then
            error stop 'scifort: value and weight arrays must have the same size', quiet=.true.
        end if

        total = 0.0_dp
        do i = 1, n
            if (ieee_is_nan(input_weights(i)) .or. .not. ieee_is_finite(input_weights(i))) then
                error stop 'scifort: weight sum must be positive and finite', quiet=.true.
            end if
            if (input_weights(i) < 0.0_dp) then
                error stop 'scifort: sample-distance weights must be nonnegative', quiet=.true.
            end if
            if (input_weights(i) > huge(total) - total) then
                error stop 'scifort: weight sum must be positive and finite', quiet=.true.
            end if
            total = total + input_weights(i)
        end do
        if (total <= 0.0_dp) then
            error stop 'scifort: weight sum must be positive and finite', quiet=.true.
        end if
        weights = input_weights / total
    end subroutine prepare_weights

    pure subroutine sort_pairs(values, weights)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout) :: weights(:)

        real(dp), allocatable :: values_work(:)
        real(dp), allocatable :: weights_work(:)

        if (size(values) <= 1) return
        allocate(values_work(size(values)), weights_work(size(weights)))
        call merge_sort_pairs(values, weights, values_work, weights_work, 1, size(values))
    end subroutine sort_pairs

    pure recursive subroutine merge_sort_pairs(values, weights, values_work, weights_work, left, right)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(inout) :: weights(:)
        real(dp), intent(inout) :: values_work(:)
        real(dp), intent(inout) :: weights_work(:)
        integer, intent(in) :: left
        integer, intent(in) :: right

        integer :: i
        integer :: j
        integer :: k
        integer :: mid

        if (left >= right) return
        mid = (left + right) / 2
        call merge_sort_pairs(values, weights, values_work, weights_work, left, mid)
        call merge_sort_pairs(values, weights, values_work, weights_work, mid + 1, right)

        i = left
        j = mid + 1
        do k = left, right
            if (i > mid) then
                values_work(k) = values(j)
                weights_work(k) = weights(j)
                j = j + 1
            else if (j > right) then
                values_work(k) = values(i)
                weights_work(k) = weights(i)
                i = i + 1
            else if (values(i) <= values(j)) then
                values_work(k) = values(i)
                weights_work(k) = weights(i)
                i = i + 1
            else
                values_work(k) = values(j)
                weights_work(k) = weights(j)
                j = j + 1
            end if
        end do
        values(left:right) = values_work(left:right)
        weights(left:right) = weights_work(left:right)
    end subroutine merge_sort_pairs

    pure logical function same_value(a, b) result(equal)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b

        equal = (a <= b .and. b <= a)
    end function same_value

    pure integer function count_positive_infinity(values) result(n)
        real(dp), intent(in) :: values(:)
        integer :: i

        n = 0
        do i = 1, size(values)
            if (.not. ieee_is_finite(values(i)) .and. .not. ieee_is_nan(values(i))) then
                if (values(i) > 0.0_dp) n = n + 1
            end if
        end do
    end function count_positive_infinity

    pure integer function count_negative_infinity(values) result(n)
        real(dp), intent(in) :: values(:)
        integer :: i

        n = 0
        do i = 1, size(values)
            if (.not. ieee_is_finite(values(i)) .and. .not. ieee_is_nan(values(i))) then
                if (values(i) < 0.0_dp) n = n + 1
            end if
        end do
    end function count_negative_infinity

end module scifort_sample_distance
