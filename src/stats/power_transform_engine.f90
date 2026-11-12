! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_power_transform_engine
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan
    use scifort_normal, only : normal_ppf
    implicit none
    private

    public :: boxcox
    public :: boxcox_fit
    public :: boxcox_llf_all_assumed_rank
    public :: boxcox_llf_reduce_assumed_rank
    public :: boxcox_normmax
    public :: boxcox_normmax_all
    public :: boxcox_llf_vector
    public :: yeojohnson
    public :: yeojohnson_fit
    public :: yeojohnson_llf_all_assumed_rank
    public :: yeojohnson_llf_reduce_assumed_rank
    public :: yeojohnson_normmax
    public :: yeojohnson_llf_vector

contains

    pure elemental function boxcox(x, lmbda) result(y)
        real(dp), intent(in) :: x !! input value
        real(dp), intent(in) :: lmbda !! Box-Cox parameter
        real(dp) :: y

        if (ieee_is_nan(x) .or. ieee_is_nan(lmbda)) then
            y = quiet_nan(0.0_dp)
        else if (x < 0.0_dp) then
            y = quiet_nan(0.0_dp)
        else if (x <= 0.0_dp .and. x >= 0.0_dp) then
            if (lmbda > 0.0_dp) then
                y = -1.0_dp / lmbda
            else
                y = negative_infinity(0.0_dp)
            end if
        else if (abs(lmbda) <= epsilon(1.0_dp)) then
            y = log(x)
        else
            y = expm1_safe(lmbda * log(x)) / lmbda
        end if
    end function boxcox

    pure elemental function yeojohnson(x, lmbda) result(y)
        real(dp), intent(in) :: x !! input value
        real(dp), intent(in) :: lmbda !! Yeo-Johnson parameter
        real(dp) :: y

        real(dp) :: exponent
        real(dp) :: log_term
        real(dp), parameter :: eps_lambda = epsilon(1.0_dp)
        real(dp), parameter :: log_huge = log(huge(1.0_dp))

        if (ieee_is_nan(x) .or. ieee_is_nan(lmbda)) then
            y = quiet_nan(0.0_dp)
            return
        end if

        if (x >= 0.0_dp) then
            log_term = log1p_safe(x)
            if (abs(lmbda) <= eps_lambda) then
                y = log_term
            else
                exponent = lmbda * log_term
                if (exponent > log_huge) then
                    y = positive_infinity(0.0_dp)
                else
                    y = expm1_safe(exponent) / lmbda
                end if
            end if
        else
            log_term = log1p_safe(-x)
            if (abs(lmbda - 2.0_dp) <= eps_lambda) then
                y = -log_term
            else
                exponent = (2.0_dp - lmbda) * log_term
                if (exponent > log_huge) then
                    y = negative_infinity(0.0_dp)
                else
                    y = -expm1_safe(exponent) / (2.0_dp - lmbda)
                end if
            end if
        end if
    end function yeojohnson

    pure function boxcox_llf_vector(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda !! Box-Cox parameter
        real(dp), intent(in) :: data(:) !! one reduction slice
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp) :: logsum
        real(dp) :: var_y
        integer :: i
        integer :: n
        character(len=16) :: policy

        call prepare_sample(data, nan_policy, x, policy, value)
        if (ieee_is_nan(value)) return
        n = size(x)
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (any(x < 0.0_dp)) then
            value = quiet_nan(0.0_dp)
            return
        end if
        if (any((x <= 0.0_dp) .and. (x >= 0.0_dp))) then
            if (lmbda > 0.0_dp .and. lmbda < 1.0_dp) then
                value = positive_infinity(0.0_dp)
            else if (lmbda > 1.0_dp) then
                y = boxcox(x, lmbda)
                var_y = population_variance(y)
                if (var_y <= 0.0_dp .and. var_y >= 0.0_dp) then
                    value = quiet_nan(0.0_dp)
                else
                    value = negative_infinity(0.0_dp)
                end if
            else
                value = quiet_nan(0.0_dp)
            end if
            return
        end if

        allocate(y(n))
        logsum = 0.0_dp
        do i = 1, n
            logsum = logsum + log(x(i))
            y(i) = boxcox(x(i), lmbda)
        end do
        if (any(.not. ieee_is_finite(y))) then
            value = quiet_nan(0.0_dp)
            return
        end if
        var_y = population_variance(y)
        if (var_y <= 0.0_dp .and. var_y >= 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else if (.not. ieee_is_finite(var_y) .or. var_y < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = (lmbda - 1.0_dp) * logsum - 0.5_dp * real(n, dp) * log(var_y)
        end if
    end function boxcox_llf_vector

    pure function yeojohnson_llf_vector(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda !! Yeo-Johnson parameter
        real(dp), intent(in) :: data(:) !! one reduction slice
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp) :: jacobian_sum
        real(dp) :: var_y
        integer :: i
        integer :: n
        character(len=16) :: policy

        call prepare_sample(data, nan_policy, x, policy, value)
        if (ieee_is_nan(value)) return
        n = size(x)
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        allocate(y(n))
        jacobian_sum = 0.0_dp
        do i = 1, n
            y(i) = yeojohnson(x(i), lmbda)
            if (x(i) > 0.0_dp) then
                jacobian_sum = jacobian_sum + log1p_safe(x(i))
            else if (x(i) < 0.0_dp) then
                jacobian_sum = jacobian_sum - log1p_safe(-x(i))
            end if
        end do
        if (any(.not. ieee_is_finite(y))) then
            value = negative_infinity(0.0_dp)
            return
        end if
        var_y = population_variance(y)
        if (var_y <= 0.0_dp .and. var_y >= 0.0_dp) then
            value = positive_infinity(0.0_dp)
        else if (.not. ieee_is_finite(var_y) .or. var_y < 0.0_dp) then
            value = negative_infinity(0.0_dp)
        else
            value = -0.5_dp * real(n, dp) * log(var_y) + (lmbda - 1.0_dp) * jacobian_sum
        end if
    end function yeojohnson_llf_vector

    subroutine boxcox_fit(x, y, lmbda)
        real(dp), intent(in) :: x(:) !! positive finite one-dimensional sample
        real(dp), allocatable, intent(out) :: y(:) !! transformed sample
        real(dp), intent(out) :: lmbda !! maximum-likelihood transform parameter

        if (size(x) == 0) then
            allocate(y(0))
            lmbda = quiet_nan(0.0_dp)
            return
        end if
        if (any(.not. ieee_is_finite(x)) .or. any(x <= 0.0_dp)) then
            error stop 'scifort: boxcox_fit data must be positive and finite', quiet=.true.
        end if
        if (all((x <= x(1)) .and. (x >= x(1)))) then
            error stop 'scifort: boxcox_fit data must not be constant', quiet=.true.
        end if
        lmbda = boxcox_normmax(x, method='mle')
        y = boxcox(x, lmbda)
    end subroutine boxcox_fit

    subroutine yeojohnson_fit(x, y, lmbda)
        real(dp), intent(in) :: x(:) !! finite one-dimensional sample
        real(dp), allocatable, intent(out) :: y(:) !! transformed sample
        real(dp), intent(out) :: lmbda !! maximum-likelihood transform parameter

        if (size(x) == 0) then
            allocate(y(0))
            lmbda = quiet_nan(0.0_dp)
            return
        end if
        if (any(.not. ieee_is_finite(x))) then
            error stop 'scifort: yeojohnson_fit data must be finite', quiet=.true.
        end if
        lmbda = yeojohnson_normmax(x)
        y = yeojohnson(x, lmbda)
    end subroutine yeojohnson_fit

    function boxcox_normmax(x, method, brack) result(lmbda)
        real(dp), intent(in) :: x(:) !! positive finite sample
        character(len=*), intent(in), optional :: method !! pearsonr or mle
        real(dp), intent(in), optional :: brack(2) !! initial search interval
        real(dp) :: lmbda

        character(len=16) :: method_value
        real(dp) :: interval(2)

        if (size(x) == 0 .or. any(.not. ieee_is_finite(x)) .or. any(x <= 0.0_dp)) then
            error stop 'scifort: boxcox_normmax data must be positive and finite', quiet=.true.
        end if
        method_value = 'pearsonr'
        if (present(method)) method_value = lowercase(trim(adjustl(method)))
        if (method_value /= 'pearsonr' .and. method_value /= 'mle') then
            error stop "scifort: boxcox_normmax method must be 'pearsonr' or 'mle'", quiet=.true.
        end if
        interval = [-2.0_dp, 2.0_dp]
        if (present(brack)) interval = brack
        if (interval(1) >= interval(2)) then
            error stop 'scifort: boxcox_normmax bracket must be increasing', quiet=.true.
        end if
        lmbda = optimize_boxcox(x, method_value, interval)
    end function boxcox_normmax

    function boxcox_normmax_all(x, brack) result(values)
        real(dp), intent(in) :: x(:) !! positive finite sample
        real(dp), intent(in), optional :: brack(2) !! initial search interval
        real(dp) :: values(2)

        values(1) = boxcox_normmax(x, method='pearsonr', brack=brack)
        values(2) = boxcox_normmax(x, method='mle', brack=brack)
    end function boxcox_normmax_all

    function yeojohnson_normmax(x, brack) result(lmbda)
        real(dp), intent(in) :: x(:) !! finite sample
        real(dp), intent(in), optional :: brack(2) !! initial search interval
        real(dp) :: lmbda

        real(dp) :: interval(2)

        if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) then
            error stop 'scifort: yeojohnson_normmax data must be finite and nonempty', quiet=.true.
        end if
        if (all(x <= 0.0_dp) .and. all(x >= 0.0_dp)) then
            lmbda = 1.0_dp
            return
        end if
        interval = [-2.0_dp, 2.0_dp]
        if (present(brack)) interval = brack
        if (interval(1) >= interval(2)) then
            error stop 'scifort: yeojohnson_normmax bracket must be increasing', quiet=.true.
        end if
        lmbda = optimize_yeojohnson(x, interval)
    end function yeojohnson_normmax

    pure function boxcox_llf_all_assumed_rank(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(..)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        call flatten_assumed_real(data, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = boxcox_llf_vector(lmbda, flat, nan_policy)
    end function boxcox_llf_all_assumed_rank

    pure function yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(..)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        call flatten_assumed_real(data, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = yeojohnson_llf_vector(lmbda, flat, nan_policy)
    end function yeojohnson_llf_all_assumed_rank

    pure subroutine boxcox_llf_reduce_assumed_rank(lmbda, data, axis, nan_policy, reduced, out_shape, ok)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(..)
        integer, intent(in) :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable, intent(out) :: reduced(:)
        integer, intent(out) :: out_shape(:)
        logical, intent(out) :: ok

        call reduce_llf_assumed_rank(lmbda, data, axis, nan_policy, reduced, out_shape, ok, .false.)
    end subroutine boxcox_llf_reduce_assumed_rank

    pure subroutine yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis, nan_policy, reduced, out_shape, ok)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(..)
        integer, intent(in) :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable, intent(out) :: reduced(:)
        integer, intent(out) :: out_shape(:)
        logical, intent(out) :: ok

        call reduce_llf_assumed_rank(lmbda, data, axis, nan_policy, reduced, out_shape, ok, .true.)
    end subroutine yeojohnson_llf_reduce_assumed_rank

    pure subroutine reduce_llf_assumed_rank(lmbda, data, axis, nan_policy, reduced, out_shape, ok, is_yj)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(..)
        integer, intent(in) :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable, intent(out) :: reduced(:)
        integer, intent(out) :: out_shape(:)
        logical, intent(out) :: ok
        logical, intent(in) :: is_yj

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
        call flatten_assumed_real(data, flat, dims)
        if (size(dims) < 1 .or. size(dims) > 15) then
            allocate(reduced(0))
            if (size(out_shape) > 0) out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > size(dims) .or. size(out_shape) /= size(dims) - 1) then
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
        if (size(out_shape) == 0) n_out = 1
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
            if (is_yj) then
                reduced(i + 1) = yeojohnson_llf_vector(lmbda, slice, nan_policy)
            else
                reduced(i + 1) = boxcox_llf_vector(lmbda, slice, nan_policy)
            end if
        end do
        ok = .true.
    end subroutine reduce_llf_assumed_rank

    pure subroutine prepare_sample(data, nan_policy, x, policy, sentinel)
        real(dp), intent(in) :: data(:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable, intent(out) :: x(:)
        character(len=16), intent(out) :: policy
        real(dp), intent(out) :: sentinel

        policy = 'propagate'
        if (present(nan_policy)) policy = lowercase(trim(adjustl(nan_policy)))
        if (.not. valid_nan_policy(policy)) then
            error stop 'scifort: invalid nan_policy', quiet=.true.
        end if
        sentinel = 0.0_dp
        if (policy == 'raise') then
            if (any(ieee_is_nan(data))) then
                error stop "scifort: nan_policy='raise' found NaN", quiet=.true.
            end if
            x = data
        else if (policy == 'propagate') then
            if (any(ieee_is_nan(data))) then
                sentinel = quiet_nan(0.0_dp)
                allocate(x(0))
                return
            end if
            x = data
        else
            x = pack(data, .not. ieee_is_nan(data))
        end if
    end subroutine prepare_sample

    pure function population_variance(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        real(dp) :: delta
        real(dp) :: mean_value
        real(dp) :: m2
        integer :: i

        if (size(x) == 0) then
            value = quiet_nan(0.0_dp)
            return
        end if
        mean_value = 0.0_dp
        m2 = 0.0_dp
        do i = 1, size(x)
            delta = x(i) - mean_value
            mean_value = mean_value + delta / real(i, dp)
            m2 = m2 + delta * (x(i) - mean_value)
        end do
        value = m2 / real(size(x), dp)
    end function population_variance

    function optimize_boxcox(x, method, initial) result(best_x)
        real(dp), intent(in) :: x(:)
        character(len=*), intent(in) :: method
        real(dp), intent(in) :: initial(2)
        real(dp) :: best_x

        best_x = grid_golden_search(x, initial, method, .false.)
    end function optimize_boxcox

    function optimize_yeojohnson(x, initial) result(best_x)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: initial(2)
        real(dp) :: best_x

        best_x = grid_golden_search(x, initial, 'mle', .true.)
    end function optimize_yeojohnson

    function grid_golden_search(x, initial, method, is_yj) result(best_x)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: initial(2)
        character(len=*), intent(in) :: method
        logical, intent(in) :: is_yj
        real(dp) :: best_x

        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fc
        real(dp) :: fd
        real(dp) :: left
        real(dp) :: right
        real(dp) :: step
        real(dp) :: width
        integer :: expand
        integer :: i
        integer :: ibest
        integer, parameter :: ngrid = 81
        real(dp), parameter :: phi = 0.6180339887498948482_dp
        real(dp), parameter :: tol = 1.0e-11_dp

        left = initial(1)
        right = initial(2)
        do expand = 0, 8
            step = (right - left) / real(ngrid - 1, dp)
            ibest = 1
            fb = objective_value(left, x, method, is_yj)
            do i = 2, ngrid
                a = left + real(i - 1, dp) * step
                fa = objective_value(a, x, method, is_yj)
                if (fa > fb) then
                    fb = fa
                    ibest = i
                end if
            end do
            if (ibest > 1 .and. ibest < ngrid) exit
            width = right - left
            if (ibest == 1) then
                right = left + step
                left = left - 2.0_dp * width
            else
                left = right - step
                right = right + 2.0_dp * width
            end if
        end do
        step = (right - left) / real(ngrid - 1, dp)
        a = max(left, left + real(ibest - 2, dp) * step)
        b = min(right, left + real(ibest, dp) * step)
        if (ibest == 1) a = left
        if (ibest == ngrid) b = right

        c = b - phi * (b - a)
        d = a + phi * (b - a)
        fc = objective_value(c, x, method, is_yj)
        fd = objective_value(d, x, method, is_yj)
        do i = 1, 200
            if (abs(b - a) <= tol * (1.0_dp + abs(c) + abs(d))) exit
            if (fc > fd) then
                b = d
                d = c
                fd = fc
                c = b - phi * (b - a)
                fc = objective_value(c, x, method, is_yj)
            else
                a = c
                c = d
                fc = fd
                d = a + phi * (b - a)
                fd = objective_value(d, x, method, is_yj)
            end if
        end do
        if (fc >= fd) then
            best_x = c
        else
            best_x = d
        end if
    end function grid_golden_search

    function objective_value(lmbda, x, method, is_yj) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: x(:)
        character(len=*), intent(in) :: method
        logical, intent(in) :: is_yj
        real(dp) :: value

        if (is_yj) then
            value = yeojohnson_llf_vector(lmbda, x)
            if (.not. ieee_is_finite(value)) value = negative_infinity(0.0_dp)
        else if (method == 'mle') then
            value = boxcox_llf_vector(lmbda, x)
            if (.not. ieee_is_finite(value)) value = negative_infinity(0.0_dp)
        else
            value = boxcox_pearson_objective(lmbda, x)
        end if
    end function objective_value

    function boxcox_pearson_objective(lmbda, x) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: x(:)
        real(dp) :: value

        integer :: i
        integer :: n
        real(dp), allocatable :: expected(:)
        real(dp), allocatable :: transformed(:)
        real(dp) :: p

        n = size(x)
        if (n < 2) then
            value = negative_infinity(0.0_dp)
            return
        end if
        allocate(expected(n), transformed(n))
        transformed = boxcox(x, lmbda)
        if (any(.not. ieee_is_finite(transformed))) then
            value = negative_infinity(0.0_dp)
            return
        end if
        call sort_real(transformed)
        do i = 1, n
            if (i == 1) then
                p = 1.0_dp - 0.5_dp ** (1.0_dp / real(n, dp))
            else if (i == n) then
                p = 0.5_dp ** (1.0_dp / real(n, dp))
            else
                p = (real(i, dp) - 0.3175_dp) / (real(n, dp) + 0.365_dp)
            end if
            expected(i) = normal_ppf(p)
        end do
        value = correlation(expected, transformed)
        if (.not. ieee_is_finite(value)) value = negative_infinity(0.0_dp)
    end function boxcox_pearson_objective

    pure function correlation(x, y) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        real(dp) :: value
        real(dp) :: mx
        real(dp) :: my
        real(dp) :: sxx
        real(dp) :: sxy
        real(dp) :: syy

        mx = sum(x) / real(size(x), dp)
        my = sum(y) / real(size(y), dp)
        sxx = sum((x - mx) * (x - mx))
        syy = sum((y - my) * (y - my))
        sxy = sum((x - mx) * (y - my))
        if (sxx <= 0.0_dp .or. syy <= 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = sxy / sqrt(sxx * syy)
        end if
    end function correlation

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:)
        if (size(values) > 1) call quicksort_real(values, 1, size(values))
    end subroutine sort_real

    recursive subroutine quicksort_real(values, left, right)
        real(dp), intent(inout) :: values(:)
        integer, intent(in) :: left
        integer, intent(in) :: right
        integer :: i
        integer :: j
        real(dp) :: pivot
        real(dp) :: temp

        i = left
        j = right
        pivot = values((left + right) / 2)
        do
            do while (values(i) < pivot)
                i = i + 1
            end do
            do while (values(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                temp = values(i)
                values(i) = values(j)
                values(j) = temp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (left < j) call quicksort_real(values, left, j)
        if (i < right) call quicksort_real(values, i, right)
    end subroutine quicksort_real

    pure function valid_nan_policy(policy) result(ok)
        character(len=*), intent(in) :: policy
        logical :: ok
        ok = policy == 'propagate' .or. policy == 'omit' .or. policy == 'raise'
    end function valid_nan_policy

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: lower
        integer :: c
        integer :: i

        lower = text
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c + 32)
        end do
    end function lowercase

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..)
        real(dp), allocatable, intent(out) :: flat(:)
        integer, allocatable, intent(out) :: dims(:)

        select rank (x)
        rank (1); dims = shape(x); flat = reshape(x, [size(x)])
        rank (2); dims = shape(x); flat = reshape(x, [size(x)])
        rank (3); dims = shape(x); flat = reshape(x, [size(x)])
        rank (4); dims = shape(x); flat = reshape(x, [size(x)])
        rank (5); dims = shape(x); flat = reshape(x, [size(x)])
        rank (6); dims = shape(x); flat = reshape(x, [size(x)])
        rank (7); dims = shape(x); flat = reshape(x, [size(x)])
        rank (8); dims = shape(x); flat = reshape(x, [size(x)])
        rank (9); dims = shape(x); flat = reshape(x, [size(x)])
        rank (10); dims = shape(x); flat = reshape(x, [size(x)])
        rank (11); dims = shape(x); flat = reshape(x, [size(x)])
        rank (12); dims = shape(x); flat = reshape(x, [size(x)])
        rank (13); dims = shape(x); flat = reshape(x, [size(x)])
        rank (14); dims = shape(x); flat = reshape(x, [size(x)])
        rank (15); dims = shape(x); flat = reshape(x, [size(x)])
        rank default
            allocate(dims(0), flat(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_power_transform_engine
