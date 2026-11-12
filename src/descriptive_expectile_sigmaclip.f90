! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_expectile_sigmaclip
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: expectile
    public :: sigmaclip

contains

    pure function expectile(a, alpha, weights) result(value)
        real(dp), intent(in) :: a(..) !! observations; all elements are used
        real(dp), intent(in), optional :: alpha !! expectile level in [0,1], default 0.5
        real(dp), intent(in), optional :: weights(..) !! NumPy-style broadcastable case weights
        real(dp) :: value

        real(dp) :: alpha_value
        integer, allocatable :: dims_a(:)
        real(dp), allocatable :: flat_a(:)
        real(dp), allocatable :: flat_w(:)
        logical :: ok

        alpha_value = 0.5_dp
        if (present(alpha)) alpha_value = alpha
        if (ieee_is_nan(alpha_value) .or. alpha_value < 0.0_dp .or. alpha_value > 1.0_dp) then
            value = quiet_nan(0.0_dp)
            return
        end if

        call flatten_assumed_real(a, flat_a, dims_a)
        if (size(dims_a) > 15 .or. size(flat_a) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (any(ieee_is_nan(flat_a))) then
            value = quiet_nan(0.0_dp)
            return
        end if

        if (present(weights)) then
            call broadcast_weights(weights, dims_a, flat_w, ok)
            if (.not. ok) then
                value = quiet_nan(0.0_dp)
                return
            end if
            if (any(ieee_is_nan(flat_w))) then
                value = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(flat_w(size(flat_a)))
            flat_w = 1.0_dp
        end if

        value = expectile_vector(flat_a, flat_w, alpha_value)
    end function expectile

    pure subroutine sigmaclip(a, clipped, lower, upper, low, high)
        real(dp), intent(in) :: a(..) !! observations; input is flattened before clipping
        real(dp), allocatable, intent(out) :: clipped(:) !! retained observations in input order
        real(dp), intent(out) :: lower !! final lower clipping threshold
        real(dp), intent(out) :: upper !! final upper clipping threshold
        real(dp), intent(in), optional :: low !! lower standard-deviation factor, default 4
        real(dp), intent(in), optional :: high !! upper standard-deviation factor, default 4

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next(:)
        integer, allocatable :: dims(:)
        real(dp) :: high_value
        real(dp) :: low_value
        real(dp) :: mean_value
        real(dp) :: std_value

        low_value = 4.0_dp
        high_value = 4.0_dp
        if (present(low)) low_value = low
        if (present(high)) high_value = high

        call flatten_assumed_real(a, current, dims)
        if (size(dims) > 15) then
            allocate(clipped(0))
            lower = quiet_nan(0.0_dp)
            upper = quiet_nan(0.0_dp)
            return
        end if

        do
            if (size(current) == 0) then
                lower = quiet_nan(0.0_dp)
                upper = quiet_nan(0.0_dp)
                clipped = current
                return
            end if
            if (any(.not. ieee_is_finite(current))) then
                allocate(clipped(0))
                lower = quiet_nan(0.0_dp)
                upper = quiet_nan(0.0_dp)
                return
            end if

            call mean_std_population(current, mean_value, std_value)
            if (.not. ieee_is_finite(mean_value) .or. .not. ieee_is_finite(std_value) .or. &
                    .not. ieee_is_finite(low_value) .or. .not. ieee_is_finite(high_value)) then
                allocate(clipped(0))
                lower = quiet_nan(0.0_dp)
                upper = quiet_nan(0.0_dp)
                return
            end if
            lower = mean_value - std_value * low_value
            upper = mean_value + std_value * high_value
            if (ieee_is_nan(lower) .or. ieee_is_nan(upper)) then
                allocate(clipped(0))
                lower = quiet_nan(0.0_dp)
                upper = quiet_nan(0.0_dp)
                return
            end if

            next = pack(current, current >= lower .and. current <= upper)
            if (size(next) == size(current)) then
                clipped = next
                return
            end if
            current = next
        end do
    end subroutine sigmaclip

    pure function expectile_vector(x, w, alpha) result(value)
        real(dp), intent(in) :: x(:) !! flattened observations
        real(dp), intent(in) :: w(:) !! flattened broadcast weights
        real(dp), intent(in) :: alpha !! validated expectile level
        real(dp) :: value

        real(dp), allocatable :: sx(:)
        real(dp), allocatable :: sw(:)
        real(dp) :: denominator
        real(dp) :: lower
        real(dp) :: numerator
        real(dp) :: sum_low
        real(dp) :: sum_total
        real(dp) :: upper
        real(dp) :: weight_low
        real(dp) :: weight_total
        integer :: k

        if (size(x) /= size(w) .or. size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        weight_total = sum(w)
        if (same_real(weight_total, 0.0_dp) .or. ieee_is_nan(weight_total)) then
            value = quiet_nan(0.0_dp)
            return
        end if

        if (same_real(alpha, 0.0_dp)) then
            value = minval(x)
            return
        else if (same_real(alpha, 1.0_dp)) then
            value = maxval(x)
            return
        end if

        if (any(.not. ieee_is_finite(x))) then
            value = weighted_mean_nonfinite(x, w)
            return
        end if

        if (same_real(alpha, 0.5_dp)) then
            value = sum(w * x) / weight_total
            return
        end if

        if (any(w < 0.0_dp)) then
            value = expectile_secant(x, w, alpha)
            return
        end if

        sx = x
        sw = w
        call sort_pairs(sx, sw)
        weight_total = sum(sw)
        sum_total = sum(sw * sx)
        weight_low = 0.0_dp
        sum_low = 0.0_dp

        do k = 0, size(sx)
            denominator = (1.0_dp - alpha) * weight_low + alpha * (weight_total - weight_low)
            if (.not. same_real(denominator, 0.0_dp)) then
                numerator = (1.0_dp - alpha) * sum_low + alpha * (sum_total - sum_low)
                value = numerator / denominator
                if (k == 0) then
                    lower = -huge(0.0_dp)
                else
                    lower = sx(k)
                end if
                if (k == size(sx)) then
                    upper = huge(0.0_dp)
                else
                    upper = sx(k + 1)
                end if
                if (value >= lower .and. value <= upper) return
            end if
            if (k < size(sx)) then
                weight_low = weight_low + sw(k + 1)
                sum_low = sum_low + sw(k + 1) * sx(k + 1)
            end if
        end do

        value = expectile_secant(x, w, alpha)
    end function expectile_vector

    pure function expectile_secant(x, w, alpha) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: w(:) !! weights with nonzero total
        real(dp), intent(in) :: alpha !! expectile level in (0,1)
        real(dp) :: value

        real(dp) :: p
        real(dp) :: p0
        real(dp) :: p1
        real(dp) :: q0
        real(dp) :: q1
        real(dp) :: ratio
        real(dp) :: swap_value
        integer :: iteration

        if (alpha >= 0.5_dp) then
            p0 = sum(w * x) / sum(w)
            p1 = maxval(x)
        else
            p0 = minval(x)
            p1 = sum(w * x) / sum(w)
        end if
        if (same_real(p0, p1)) then
            value = p0
            return
        end if

        q0 = expectile_equation(p0, x, w, alpha)
        q1 = expectile_equation(p1, x, w, alpha)
        if (abs(q1) < abs(q0)) then
            swap_value = p0
            p0 = p1
            p1 = swap_value
            swap_value = q0
            q0 = q1
            q1 = swap_value
        end if

        do iteration = 1, 50
            if (same_real(q1, q0)) then
                value = 0.5_dp * (p1 + p0)
                return
            end if
            if (abs(q1) > abs(q0)) then
                ratio = q0 / q1
                p = (-ratio * p1 + p0) / (1.0_dp - ratio)
            else
                ratio = q1 / q0
                p = (-ratio * p0 + p1) / (1.0_dp - ratio)
            end if
            if (.not. ieee_is_finite(p)) then
                value = quiet_nan(0.0_dp)
                return
            end if
            if (abs(p - p1) <= 1.48e-8_dp) then
                value = p
                return
            end if
            p0 = p1
            q0 = q1
            p1 = p
            q1 = expectile_equation(p1, x, w, alpha)
        end do
        value = p1
    end function expectile_secant

    pure function expectile_equation(t, x, w, alpha) result(value)
        real(dp), intent(in) :: t !! trial expectile
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: w(:) !! case weights
        real(dp), intent(in) :: alpha !! expectile level
        real(dp) :: value

        real(dp) :: coefficient
        integer :: i

        value = 0.0_dp
        do i = 1, size(x)
            if (x(i) <= t) then
                coefficient = 1.0_dp - alpha
            else
                coefficient = alpha
            end if
            value = value + w(i) * coefficient * (t - x(i))
        end do
    end function expectile_equation

    pure function weighted_mean_nonfinite(x, w) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: w(:) !! weights
        real(dp) :: value

        integer :: i
        logical :: negative_infinity_present
        logical :: positive_infinity_present
        real(dp) :: weighted_sum
        real(dp) :: weight_total

        negative_infinity_present = .false.
        positive_infinity_present = .false.
        weighted_sum = 0.0_dp
        weight_total = sum(w)
        do i = 1, size(x)
            if (same_real(w(i), 0.0_dp) .and. .not. ieee_is_finite(x(i))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            if (.not. ieee_is_finite(x(i))) then
                if (x(i) > 0.0_dp .and. .not. same_real(w(i), 0.0_dp)) positive_infinity_present = .true.
                if (x(i) < 0.0_dp .and. .not. same_real(w(i), 0.0_dp)) negative_infinity_present = .true.
            else
                weighted_sum = weighted_sum + w(i) * x(i)
            end if
        end do
        if (positive_infinity_present .and. negative_infinity_present) then
            value = quiet_nan(0.0_dp)
        else if (positive_infinity_present) then
            value = positive_infinity(0.0_dp)
        else if (negative_infinity_present) then
            value = negative_infinity(0.0_dp)
        else
            value = weighted_sum / weight_total
        end if
    end function weighted_mean_nonfinite

    pure subroutine broadcast_weights(weights, target_dims, expanded, ok)
        real(dp), intent(in) :: weights(..) !! scalar or broadcastable weight array
        integer, intent(in) :: target_dims(:) !! target extents
        real(dp), allocatable, intent(out) :: expanded(:) !! weights broadcast to target shape
        logical, intent(out) :: ok !! true when broadcasting succeeds

        integer :: coordinate
        integer :: d
        integer, allocatable :: dims_w(:)
        real(dp), allocatable :: flat_w(:)
        integer :: linear
        integer :: offset
        integer :: rank_w
        integer :: source_coordinate
        integer :: source_index
        integer, allocatable :: strides_w(:)
        integer :: target_index
        integer :: target_size

        call flatten_assumed_real(weights, flat_w, dims_w)
        rank_w = size(dims_w)
        if (rank_w > size(target_dims)) then
            allocate(expanded(0))
            ok = .false.
            return
        end if
        offset = size(target_dims) - rank_w
        do d = 1, rank_w
            if (dims_w(d) /= 1 .and. dims_w(d) /= target_dims(offset + d)) then
                allocate(expanded(0))
                ok = .false.
                return
            end if
        end do

        target_size = product(target_dims)
        if (size(target_dims) == 0) target_size = 1
        allocate(expanded(target_size))
        allocate(strides_w(rank_w))
        if (rank_w > 0) then
            strides_w(1) = 1
            do d = 2, rank_w
                strides_w(d) = strides_w(d - 1) * dims_w(d - 1)
            end do
        end if

        do target_index = 0, target_size - 1
            linear = target_index
            source_index = 1
            do d = 1, size(target_dims)
                coordinate = mod(linear, target_dims(d))
                linear = linear / target_dims(d)
                if (d > offset) then
                    source_coordinate = coordinate
                    if (dims_w(d - offset) == 1) source_coordinate = 0
                    source_index = source_index + source_coordinate * strides_w(d - offset)
                end if
            end do
            expanded(target_index + 1) = flat_w(source_index)
        end do
        ok = .true.
    end subroutine broadcast_weights

    pure subroutine flatten_assumed_real(a, flat, dims)
        real(dp), intent(in) :: a(..) !! scalar or array to flatten in Fortran element order
        real(dp), allocatable, intent(out) :: flat(:) !! flattened values
        integer, allocatable, intent(out) :: dims(:) !! source extents

        select rank (source => a)
        rank (0)
            allocate(dims(0), flat(1))
            flat(1) = source
        rank (1)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (2)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (3)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (4)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (5)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (6)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (7)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (8)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (9)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (10)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (11)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (12)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (13)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (14)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        rank (15)
            dims = shape(source)
            flat = reshape(source, [size(source)])
        end select
    end subroutine flatten_assumed_real

    pure subroutine mean_std_population(x, mean_value, std_value)
        real(dp), intent(in) :: x(:) !! finite observations
        real(dp), intent(out) :: mean_value !! arithmetic mean
        real(dp), intent(out) :: std_value !! population standard deviation

        real(dp) :: delta
        real(dp) :: mean_running
        real(dp) :: m2
        integer :: i

        mean_running = 0.0_dp
        m2 = 0.0_dp
        do i = 1, size(x)
            delta = x(i) - mean_running
            mean_running = mean_running + delta / real(i, dp)
            m2 = m2 + delta * (x(i) - mean_running)
        end do
        mean_value = mean_running
        if (m2 < 0.0_dp .and. abs(m2) <= 16.0_dp * epsilon(1.0_dp) * &
                max(1.0_dp, abs(mean_value))) m2 = 0.0_dp
        std_value = sqrt(max(0.0_dp, m2 / real(size(x), dp)))
    end subroutine mean_std_population

    pure subroutine sort_pairs(values, weights)
        real(dp), intent(inout) :: values(:) !! values to sort ascending
        real(dp), intent(inout) :: weights(:) !! weights permuted with values

        integer :: i
        integer :: j
        real(dp) :: key_value
        real(dp) :: key_weight

        do i = 2, size(values)
            key_value = values(i)
            key_weight = weights(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key_value) exit
                values(j + 1) = values(j)
                weights(j + 1) = weights(j)
                j = j - 1
            end do
            values(j + 1) = key_value
            weights(j + 1) = key_weight
        end do
    end subroutine sort_pairs

    pure elemental function same_real(a, b) result(equal)
        real(dp), intent(in) :: a !! first value
        real(dp), intent(in) :: b !! second value
        logical :: equal

        equal = a <= b .and. b <= a
    end function same_real

end module scifort_descriptive_expectile_sigmaclip
