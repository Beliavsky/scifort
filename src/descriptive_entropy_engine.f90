! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_entropy_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_special_elementary, only : entr, rel_entr
    implicit none
    private

    public :: differential_entropy_all_assumed_rank
    public :: differential_entropy_reduce_assumed_rank
    public :: differential_entropy_vector
    public :: entropy_all_assumed_rank
    public :: entropy_reduce_assumed_rank
    public :: entropy_vector

contains

    pure function entropy_vector(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:) !! probability weights in one reduction slice
        real(dp), intent(in), optional :: qk(:) !! comparison probability weights
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: p(:)
        real(dp), allocatable :: q(:)
        logical, allocatable :: keep(:)
        real(dp) :: p_sum
        character(len=16) :: policy
        real(dp) :: q_sum

        if (present(base)) then
            if (.not. ieee_is_nan(base)) then
                if (base <= 0.0_dp) then
                    error stop 'scifort: entropy base must be positive', quiet=.true.
                end if
            end if
        end if
        if (present(qk)) then
            if (size(qk) /= size(pk)) then
                error stop 'scifort: entropy paired slices must have equal size', quiet=.true.
            end if
        end if

        policy = 'propagate'
        if (present(nan_policy)) policy = lowercase(trim(adjustl(nan_policy)))
        if (.not. valid_nan_policy(policy)) then
            error stop 'scifort: invalid nan_policy', quiet=.true.
        end if

        if (present(qk)) then
            if (policy == 'raise') then
                if (any(ieee_is_nan(pk)) .or. any(ieee_is_nan(qk))) then
                    error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
                end if
                p = pk
                q = qk
            else if (policy == 'propagate') then
                if (any(ieee_is_nan(pk)) .or. any(ieee_is_nan(qk))) then
                    value = quiet_nan(0.0_dp)
                    return
                end if
                p = pk
                q = qk
            else
                allocate(keep(size(pk)))
                keep = .not. ieee_is_nan(pk) .and. .not. ieee_is_nan(qk)
                p = pack(pk, keep)
                q = pack(qk, keep)
            end if
        else
            if (policy == 'raise') then
                if (any(ieee_is_nan(pk))) then
                    error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
                end if
                p = pk
            else if (policy == 'propagate') then
                if (any(ieee_is_nan(pk))) then
                    value = quiet_nan(0.0_dp)
                    return
                end if
                p = pk
            else
                p = pack(pk, .not. ieee_is_nan(pk))
            end if
        end if

        if (size(p) == 0) then
            value = 0.0_dp
            if (present(base)) value = divide_by_log_base(value, base)
            return
        end if

        p_sum = sum(p)
        call normalize_probability_vector(p, p_sum)
        if (present(qk)) then
            q_sum = sum(q)
            call normalize_probability_vector(q, q_sum)
            value = sum(rel_entr(p, q))
        else
            value = sum(entr(p))
        end if
        if (present(base)) value = divide_by_log_base(value, base)
    end function entropy_vector

    pure function entropy_all_assumed_rank(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(..) !! probability array, rank 1 through 15
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer, allocatable :: joint_dims(:)
        integer, allocatable :: p_dims(:)
        real(dp), allocatable :: p_flat(:)
        real(dp), allocatable :: p_joint(:)
        integer, allocatable :: q_dims(:)
        real(dp), allocatable :: q_flat(:)
        real(dp), allocatable :: q_joint(:)
        logical :: ok

        call flatten_assumed_real(pk, p_flat, p_dims)
        if (size(p_dims) < 1 .or. size(p_dims) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if

        if (present(qk)) then
            call flatten_assumed_real(qk, q_flat, q_dims)
            if (size(q_dims) > 15) then
                value = quiet_nan(0.0_dp)
                return
            end if
            call broadcast_pair(p_flat, p_dims, q_flat, q_dims, p_joint, q_joint, &
                joint_dims, ok)
            if (.not. ok) then
                error stop 'scifort: entropy pk and qk are not broadcastable', quiet=.true.
            end if
            value = entropy_vector(p_joint, q_joint, base, nan_policy)
        else
            value = entropy_vector(p_flat, base=base, nan_policy=nan_policy)
        end if
    end function entropy_all_assumed_rank

    pure subroutine entropy_reduce_assumed_rank(pk, qk, axis, base, nan_policy, &
            reduced, out_shape, ok)
        real(dp), intent(in) :: pk(..) !! probability array, rank 1 through 15
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction output
        integer, intent(out) :: out_shape(:) !! joint broadcast shape without reduced axis
        logical, intent(out) :: ok !! true when rank/axis/broadcasting are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer :: i
        integer, allocatable :: joint_dims(:)
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        integer, allocatable :: p_dims(:)
        real(dp), allocatable :: p_flat(:)
        real(dp), allocatable :: p_joint(:)
        real(dp), allocatable :: p_slice(:)
        integer, allocatable :: q_dims(:)
        real(dp), allocatable :: q_flat(:)
        real(dp), allocatable :: q_joint(:)
        real(dp), allocatable :: q_slice(:)
        integer, allocatable :: strides(:)

        ok = .false.
        call flatten_assumed_real(pk, p_flat, p_dims)
        if (size(p_dims) < 1 .or. size(p_dims) > 15) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (size(out_shape) /= size(p_dims) - 1) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > size(p_dims)) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if

        if (present(qk)) then
            call flatten_assumed_real(qk, q_flat, q_dims)
            if (size(q_dims) > size(p_dims)) then
                error stop 'scifort: entropy qk rank may not exceed pk rank for axis reduction', &
                    quiet=.true.
            end if
            call broadcast_pair_to_rank(p_flat, p_dims, q_flat, q_dims, p_joint, q_joint, &
                joint_dims, ok)
            if (.not. ok) then
                error stop 'scifort: entropy pk and qk are not broadcastable', quiet=.true.
            end if
        else
            joint_dims = p_dims
            p_joint = p_flat
        end if

        allocate(strides(size(joint_dims)))
        strides(1) = 1
        do d = 2, size(joint_dims)
            strides(d) = strides(d - 1) * joint_dims(d - 1)
        end do

        out_dim = 0
        do d = 1, size(joint_dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = joint_dims(d)
        end do
        n_reduce = joint_dims(axis)
        n_out = product(out_shape)
        allocate(reduced(n_out))
        allocate(p_slice(n_reduce))
        if (present(qk)) allocate(q_slice(n_reduce))

        do i = 0, n_out - 1
            linear = i
            base_index = 1
            out_dim = 0
            do d = 1, size(joint_dims)
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, out_shape(out_dim))
                linear = linear / out_shape(out_dim)
                base_index = base_index + coordinate * strides(d)
            end do
            do d = 1, n_reduce
                p_slice(d) = p_joint(base_index + (d - 1) * strides(axis))
                if (present(qk)) then
                    q_slice(d) = q_joint(base_index + (d - 1) * strides(axis))
                end if
            end do
            if (present(qk)) then
                reduced(i + 1) = entropy_vector(p_slice, q_slice, base, nan_policy)
            else
                reduced(i + 1) = entropy_vector(p_slice, base=base, nan_policy=nan_policy)
            end if
        end do
        ok = .true.
    end subroutine entropy_reduce_assumed_rank

    pure function differential_entropy_vector(values, window_length, base, method, &
            nan_policy) result(value)
        real(dp), intent(in) :: values(:) !! sample values in one reduction slice
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: clean(:)
        integer :: m
        character(len=16) :: method_value
        integer :: n
        character(len=16) :: policy

        if (present(base)) then
            if (.not. ieee_is_nan(base)) then
                if (base <= 0.0_dp) then
                    error stop 'scifort: differential_entropy base must be positive', quiet=.true.
                end if
            end if
        end if

        policy = 'propagate'
        if (present(nan_policy)) policy = lowercase(trim(adjustl(nan_policy)))
        if (.not. valid_nan_policy(policy)) then
            error stop 'scifort: invalid nan_policy', quiet=.true.
        end if

        if (policy == 'raise') then
            if (any(ieee_is_nan(values))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            clean = values
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(values))) then
                value = quiet_nan(0.0_dp)
                return
            end if
            clean = values
        else
            clean = pack(values, .not. ieee_is_nan(values))
        end if

        n = size(clean)
        if (present(window_length)) then
            m = window_length
        else
            m = int(floor(sqrt(real(n, dp)) + 0.5_dp))
        end if
        if (m < 1 .or. 2 * m >= n) then
            value = quiet_nan(0.0_dp)
            return
        end if

        method_value = 'auto'
        if (present(method)) method_value = lowercase(trim(adjustl(method)))
        if (.not. valid_entropy_method(method_value)) then
            error stop 'scifort: invalid differential_entropy method', quiet=.true.
        end if
        if (method_value == 'auto') then
            if (n <= 10) then
                method_value = 'van es'
            else if (n <= 1000) then
                method_value = 'ebrahimi'
            else
                method_value = 'vasicek'
            end if
        end if

        call sort_values(clean)
        select case (method_value)
        case ('vasicek')
            value = vasicek_entropy(clean, m)
        case ('van es')
            value = van_es_entropy(clean, m)
        case ('ebrahimi')
            value = ebrahimi_entropy(clean, m)
        case ('correa')
            value = correa_entropy(clean, m)
        end select
        if (present(base)) value = divide_by_log_base(value, base)
    end function differential_entropy_vector

    pure function differential_entropy_all_assumed_rank(values, window_length, base, &
            method, nan_policy) result(value)
        real(dp), intent(in) :: values(..) !! sample array, rank 1 through 15
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! estimator selector
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        call flatten_assumed_real(values, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = differential_entropy_vector(flat, window_length, base, method, nan_policy)
    end function differential_entropy_all_assumed_rank

    pure subroutine differential_entropy_reduce_assumed_rank(values, axis, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        real(dp), intent(in) :: values(..) !! sample array, rank 1 through 15
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! estimator selector
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction output
        integer, intent(out) :: out_shape(:) !! input shape without reduced axis
        logical, intent(out) :: ok !! true when rank and axis are valid

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)

        ok = .false.
        call flatten_assumed_real(values, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (size(out_shape) /= size(dims) - 1) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > size(dims)) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if

        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do
        out_dim = 0
        do d = 1, size(dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = dims(d)
        end do
        n_reduce = dims(axis)
        n_out = product(out_shape)
        allocate(reduced(n_out))
        allocate(slice(n_reduce))

        do i = 0, n_out - 1
            linear = i
            base_index = 1
            out_dim = 0
            do d = 1, size(dims)
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, out_shape(out_dim))
                linear = linear / out_shape(out_dim)
                base_index = base_index + coordinate * strides(d)
            end do
            do d = 1, n_reduce
                slice(d) = flat(base_index + (d - 1) * strides(axis))
            end do
            reduced(i + 1) = differential_entropy_vector(slice, window_length, base, &
                method, nan_policy)
        end do
        ok = .true.
    end subroutine differential_entropy_reduce_assumed_rank

    pure subroutine normalize_probability_vector(values, total)
        real(dp), intent(inout) :: values(:) !! probability weights to normalize
        real(dp), intent(in) :: total !! sum of probability weights

        integer :: i

        if (ieee_is_nan(total)) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (abs(total) > 0.0_dp) then
            values = values / total
            return
        end if
        do i = 1, size(values)
            if (values(i) > 0.0_dp) then
                values(i) = positive_infinity(0.0_dp)
            else if (values(i) < 0.0_dp) then
                values(i) = negative_infinity(0.0_dp)
            else
                values(i) = quiet_nan(0.0_dp)
            end if
        end do
    end subroutine normalize_probability_vector

    pure function divide_by_log_base(value, base) result(scaled)
        real(dp), intent(in) :: value !! entropy in natural-log units
        real(dp), intent(in) :: base !! requested logarithm base
        real(dp) :: scaled

        real(dp) :: log_base

        if (ieee_is_nan(value) .or. ieee_is_nan(base)) then
            scaled = quiet_nan(0.0_dp)
            return
        end if
        log_base = log(base)
        if (abs(log_base) <= 0.0_dp) then
            if (abs(value) <= 0.0_dp) then
                scaled = quiet_nan(0.0_dp)
            else if (value > 0.0_dp) then
                scaled = positive_infinity(0.0_dp)
            else
                scaled = negative_infinity(0.0_dp)
            end if
            return
        end if
        if (.not. ieee_is_finite(log_base)) then
            if (.not. ieee_is_finite(value)) then
                scaled = quiet_nan(0.0_dp)
            else
                scaled = 0.0_dp
            end if
            return
        end if
        scaled = value / log_base
    end function divide_by_log_base

    pure function vasicek_entropy(x, m) result(value)
        real(dp), intent(in) :: x(:) !! sorted sample
        integer, intent(in) :: m !! spacing window length
        real(dp) :: value

        real(dp) :: difference
        integer :: i
        integer :: left
        integer :: n
        integer :: right
        real(dp) :: total

        n = size(x)
        total = 0.0_dp
        do i = 1, n
            left = max(1, i - m)
            right = min(n, i + m)
            difference = x(right) - x(left)
            total = total + log_nonnegative(real(n, dp) * difference / real(2 * m, dp))
        end do
        value = total / real(n, dp)
    end function vasicek_entropy

    pure function van_es_entropy(x, m) result(value)
        real(dp), intent(in) :: x(:) !! sorted sample
        integer, intent(in) :: m !! spacing window length
        real(dp) :: value

        real(dp) :: difference
        integer :: i
        integer :: k
        integer :: n
        real(dp) :: term1
        real(dp) :: harmonic

        n = size(x)
        term1 = 0.0_dp
        do i = 1, n - m
            difference = x(i + m) - x(i)
            term1 = term1 + log_nonnegative(real(n + 1, dp) * difference / real(m, dp))
        end do
        term1 = term1 / real(n - m, dp)
        harmonic = 0.0_dp
        do k = m, n
            harmonic = harmonic + 1.0_dp / real(k, dp)
        end do
        value = term1 + harmonic + log(real(m, dp)) - log(real(n + 1, dp))
    end function van_es_entropy

    pure function ebrahimi_entropy(x, m) result(value)
        real(dp), intent(in) :: x(:) !! sorted sample
        integer, intent(in) :: m !! spacing window length
        real(dp) :: value

        real(dp) :: ci
        real(dp) :: difference
        integer :: i
        integer :: left
        integer :: n
        integer :: right
        real(dp) :: total

        n = size(x)
        total = 0.0_dp
        do i = 1, n
            left = max(1, i - m)
            right = min(n, i + m)
            difference = x(right) - x(left)
            if (i <= m) then
                ci = 1.0_dp + real(i - 1, dp) / real(m, dp)
            else if (i >= n - m + 1) then
                ci = 1.0_dp + real(n - i, dp) / real(m, dp)
            else
                ci = 2.0_dp
            end if
            total = total + log_nonnegative(real(n, dp) * difference / (ci * real(m, dp)))
        end do
        value = total / real(n, dp)
    end function ebrahimi_entropy

    pure function correa_entropy(x, m) result(value)
        real(dp), intent(in) :: x(:) !! sorted sample
        integer, intent(in) :: m !! spacing window length
        real(dp) :: value

        real(dp) :: denominator
        real(dp) :: difference
        integer :: i
        integer :: j
        integer :: jj
        integer :: n
        real(dp) :: numerator
        real(dp) :: ratio
        real(dp) :: total
        real(dp) :: xbar

        n = size(x)
        total = 0.0_dp
        do i = 1, n
            xbar = 0.0_dp
            do j = -m, m
                jj = min(n, max(1, i + j))
                xbar = xbar + x(jj)
            end do
            xbar = xbar / real(2 * m + 1, dp)

            numerator = 0.0_dp
            denominator = 0.0_dp
            do j = -m, m
                jj = min(n, max(1, i + j))
                difference = x(jj) - xbar
                numerator = numerator + difference * real(j, dp)
                denominator = denominator + difference * difference
            end do
            denominator = real(n, dp) * denominator
            ratio = safe_ratio(numerator, denominator)
            total = total + log_nonnegative(ratio)
        end do
        value = -total / real(n, dp)
    end function correa_entropy

    pure function safe_ratio(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! ratio numerator
        real(dp), intent(in) :: denominator !! ratio denominator
        real(dp) :: value

        if (ieee_is_nan(numerator) .or. ieee_is_nan(denominator)) then
            value = quiet_nan(0.0_dp)
        else if (abs(denominator) > 0.0_dp) then
            value = numerator / denominator
        else if (abs(numerator) <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else if ((numerator > 0.0_dp .and. denominator >= 0.0_dp) .or. &
                (numerator < 0.0_dp .and. denominator < 0.0_dp)) then
            value = positive_infinity(0.0_dp)
        else
            value = negative_infinity(0.0_dp)
        end if
    end function safe_ratio

    pure function log_nonnegative(x) result(value)
        real(dp), intent(in) :: x !! value whose logarithm is requested
        real(dp) :: value

        if (ieee_is_nan(x)) then
            value = quiet_nan(0.0_dp)
        else if (x > 0.0_dp) then
            value = log(x)
        else if (x < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = negative_infinity(0.0_dp)
        end if
    end function log_nonnegative

    pure logical function valid_nan_policy(policy) result(valid)
        character(len=*), intent(in) :: policy !! requested NaN handling policy

        select case (trim(adjustl(policy)))
        case ('propagate', 'omit', 'raise')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_nan_policy

    pure logical function valid_entropy_method(method) result(valid)
        character(len=*), intent(in) :: method !! requested differential-entropy estimator

        select case (trim(adjustl(method)))
        case ('vasicek', 'van es', 'ebrahimi', 'correa', 'auto')
            valid = .true.
        case default
            valid = .false.
        end select
    end function valid_entropy_method

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text !! input text
        character(len=len(text)) :: lower

        integer :: i
        integer :: code

        lower = text
        do i = 1, len(text)
            code = iachar(lower(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

    pure subroutine broadcast_pair(p_flat, p_dims, q_flat, q_dims, p_joint, q_joint, &
            joint_dims, ok)
        real(dp), intent(in) :: p_flat(:) !! flattened first operand
        integer, intent(in) :: p_dims(:) !! first operand shape
        real(dp), intent(in) :: q_flat(:) !! flattened second operand
        integer, intent(in) :: q_dims(:) !! second operand shape
        real(dp), allocatable, intent(out) :: p_joint(:) !! broadcast first operand
        real(dp), allocatable, intent(out) :: q_joint(:) !! broadcast second operand
        integer, allocatable, intent(out) :: joint_dims(:) !! joint broadcast shape
        logical, intent(out) :: ok !! true when broadcasting succeeds

        integer :: joint_rank

        joint_rank = max(size(p_dims), size(q_dims))
        allocate(joint_dims(joint_rank))
        call joint_broadcast_shape(p_dims, q_dims, joint_dims, ok)
        if (.not. ok) then
            allocate(p_joint(0), q_joint(0))
            return
        end if
        call broadcast_flat(p_flat, p_dims, joint_dims, p_joint)
        call broadcast_flat(q_flat, q_dims, joint_dims, q_joint)
    end subroutine broadcast_pair

    pure subroutine broadcast_pair_to_rank(p_flat, p_dims, q_flat, q_dims, p_joint, &
            q_joint, joint_dims, ok)
        real(dp), intent(in) :: p_flat(:) !! flattened first operand
        integer, intent(in) :: p_dims(:) !! first operand shape
        real(dp), intent(in) :: q_flat(:) !! flattened second operand
        integer, intent(in) :: q_dims(:) !! second operand shape
        real(dp), allocatable, intent(out) :: p_joint(:) !! broadcast first operand
        real(dp), allocatable, intent(out) :: q_joint(:) !! broadcast second operand
        integer, allocatable, intent(out) :: joint_dims(:) !! joint shape at rank(p)
        logical, intent(out) :: ok !! true when broadcasting succeeds

        allocate(joint_dims(size(p_dims)))
        call joint_broadcast_shape(p_dims, q_dims, joint_dims, ok)
        if (.not. ok) then
            allocate(p_joint(0), q_joint(0))
            return
        end if
        call broadcast_flat(p_flat, p_dims, joint_dims, p_joint)
        call broadcast_flat(q_flat, q_dims, joint_dims, q_joint)
    end subroutine broadcast_pair_to_rank

    pure subroutine joint_broadcast_shape(a_dims, b_dims, joint_dims, ok)
        integer, intent(in) :: a_dims(:) !! first shape
        integer, intent(in) :: b_dims(:) !! second shape
        integer, intent(out) :: joint_dims(:) !! trailing-aligned broadcast shape
        logical, intent(out) :: ok !! true when shapes are broadcast-compatible

        integer :: a_extent
        integer :: a_offset
        integer :: b_extent
        integer :: b_offset
        integer :: d

        ok = .false.
        if (size(joint_dims) /= max(size(a_dims), size(b_dims))) return
        a_offset = size(joint_dims) - size(a_dims)
        b_offset = size(joint_dims) - size(b_dims)
        do d = 1, size(joint_dims)
            a_extent = 1
            if (d > a_offset) a_extent = a_dims(d - a_offset)
            b_extent = 1
            if (d > b_offset) b_extent = b_dims(d - b_offset)
            if (a_extent /= b_extent .and. a_extent /= 1 .and. b_extent /= 1) return
            joint_dims(d) = max(a_extent, b_extent)
        end do
        ok = .true.
    end subroutine joint_broadcast_shape

    pure subroutine broadcast_flat(source, source_dims, joint_dims, target)
        real(dp), intent(in) :: source(:) !! source values in Fortran order
        integer, intent(in) :: source_dims(:) !! source shape
        integer, intent(in) :: joint_dims(:) !! target broadcast shape
        real(dp), allocatable, intent(out) :: target(:) !! broadcast values

        integer :: coordinate
        integer :: d
        integer :: joint_linear
        integer :: linear
        integer :: offset
        integer :: source_d
        integer :: source_index
        integer, allocatable :: source_strides(:)

        allocate(target(product(joint_dims)))
        if (size(source_dims) == 0) then
            if (size(source) == 1) target = source(1)
            return
        end if
        allocate(source_strides(size(source_dims)))
        source_strides(1) = 1
        do d = 2, size(source_dims)
            source_strides(d) = source_strides(d - 1) * source_dims(d - 1)
        end do
        offset = size(joint_dims) - size(source_dims)
        do joint_linear = 0, size(target) - 1
            linear = joint_linear
            source_index = 1
            do d = 1, size(joint_dims)
                coordinate = mod(linear, joint_dims(d))
                linear = linear / joint_dims(d)
                if (d <= offset) cycle
                source_d = d - offset
                if (source_dims(source_d) /= 1) then
                    source_index = source_index + coordinate * source_strides(source_d)
                end if
            end do
            target(joint_linear + 1) = source(source_index)
        end do
    end subroutine broadcast_flat

    pure subroutine sort_values(values)
        real(dp), intent(inout) :: values(:) !! values sorted in ascending order

        real(dp), allocatable :: temporary(:)

        if (size(values) <= 1) return
        allocate(temporary(size(values)))
        call merge_sort(values, temporary, 1, size(values))
    end subroutine sort_values

    recursive pure subroutine merge_sort(values, temporary, left, right)
        real(dp), intent(inout) :: values(:) !! values being sorted
        real(dp), intent(inout) :: temporary(:) !! merge workspace
        integer, intent(in) :: left !! first inclusive index
        integer, intent(in) :: right !! last inclusive index

        integer :: i
        integer :: j
        integer :: k
        integer :: middle

        if (left >= right) return
        middle = left + (right - left) / 2
        call merge_sort(values, temporary, left, middle)
        call merge_sort(values, temporary, middle + 1, right)
        i = left
        j = middle + 1
        do k = left, right
            if (i > middle) then
                temporary(k) = values(j)
                j = j + 1
            else if (j > right) then
                temporary(k) = values(i)
                i = i + 1
            else if (values(i) <= values(j)) then
                temporary(k) = values(i)
                i = i + 1
            else
                temporary(k) = values(j)
                j = j + 1
            end if
        end do
        values(left:right) = temporary(left:right)
    end subroutine merge_sort

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..) !! real scalar/array of rank 0 through 15
        real(dp), allocatable, intent(out) :: flat(:) !! Fortran-order flattened values
        integer, allocatable, intent(out) :: dims(:) !! input shape

        select rank (x)
        rank (0)
            allocate(dims(0), flat(1))
            flat(1) = x
        rank (1)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (2)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (3)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (4)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (5)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (6)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (7)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (8)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (9)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (10)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (11)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (12)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (13)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (14)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank (15)
            dims = shape(x)
            flat = reshape(x, [size(x)])
        rank default
            allocate(dims(0), flat(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_descriptive_entropy_engine
