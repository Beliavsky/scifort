! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_highrank_stats_engine
    use scifort_descriptive_lowrank, only : sem_low => sem, variation_low => variation, &
        iqr_low => iqr, mad_low => median_abs_deviation, moment_low => moment, &
        skew_low => skew, kurtosis_low => kurtosis
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: stat_kind_sem = 1
    integer, parameter, public :: stat_kind_variation = 2
    integer, parameter, public :: stat_kind_iqr = 3
    integer, parameter, public :: stat_kind_mad = 4
    integer, parameter, public :: stat_kind_moment = 5
    integer, parameter, public :: stat_kind_skew = 6
    integer, parameter, public :: stat_kind_kurtosis = 7

    public :: reduce_assumed_rank_stat
    public :: iqr_assumed_rank_all

contains

    pure subroutine reduce_assumed_rank_stat(x, axis, stat_kind, reduced, out_shape, ok, &
            ddof, rng, scale, order, center, bias, fisher, nan_policy)
        real(dp), intent(in) :: x(..) !! observations with rank 5 to 15
        integer, intent(in) :: axis !! Fortran axis to reduce
        integer, intent(in) :: stat_kind !! internal descriptive-statistic selector
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! result shape with reduced axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid
        integer, intent(in), optional :: ddof !! divisor adjustment for sem or variation
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints for IQR
        real(dp), intent(in), optional :: scale !! divisor for IQR or MAD
        integer, intent(in), optional :: order !! nonnegative moment order
        real(dp), intent(in), optional :: center !! fixed center for moment
        logical, intent(in), optional :: bias !! biased skew/kurtosis estimator
        logical, intent(in), optional :: fisher !! Fisher kurtosis when true
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise policy

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
        if (rank(x) < 5 .or. rank(x) > 15) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (size(out_shape) /= rank(x) - 1) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > rank(x)) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if

        call flatten_assumed_real(x, flat, dims)
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
            reduced(i + 1) = evaluate_stat(slice, stat_kind, ddof, rng, scale, order, &
                center, bias, fisher, nan_policy)
        end do
        ok = .true.
    end subroutine reduce_assumed_rank_stat

    pure function iqr_assumed_rank_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(..) !! observations with rank 5 to 15
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to percentile range
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise policy
        real(dp) :: value

        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)

        if (rank(x) < 5 .or. rank(x) > 15) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call flatten_assumed_real(x, flat, dims)
        if (present(nan_policy)) then
            if (present(rng)) then
                if (present(scale)) then
                    value = iqr_low(flat, rng=rng, scale=scale, nan_policy=nan_policy)
                else
                    value = iqr_low(flat, rng=rng, nan_policy=nan_policy)
                end if
            else if (present(scale)) then
                value = iqr_low(flat, scale=scale, nan_policy=nan_policy)
            else
                value = iqr_low(flat, nan_policy=nan_policy)
            end if
        else
            if (present(rng)) then
                if (present(scale)) then
                    value = iqr_low(flat, rng=rng, scale=scale)
                else
                    value = iqr_low(flat, rng=rng)
                end if
            else if (present(scale)) then
                value = iqr_low(flat, scale=scale)
            else
                value = iqr_low(flat)
            end if
        end if
    end function iqr_assumed_rank_all

    pure function evaluate_stat(x, stat_kind, ddof, rng, scale, order, center, bias, &
            fisher, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! one-dimensional reduction slice
        integer, intent(in) :: stat_kind !! internal descriptive-statistic selector
        integer, intent(in), optional :: ddof !! divisor adjustment
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints
        real(dp), intent(in), optional :: scale !! scale divisor
        integer, intent(in), optional :: order !! moment order
        real(dp), intent(in), optional :: center !! moment center
        logical, intent(in), optional :: bias !! biased estimator selector
        logical, intent(in), optional :: fisher !! Fisher kurtosis selector
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise policy
        real(dp) :: value

        logical :: bias_value
        integer :: ddof_value
        logical :: fisher_value
        integer :: order_value
        real(dp) :: rng_value(2)
        real(dp) :: scale_value

        select case (stat_kind)
        case (stat_kind_sem)
            ddof_value = 1
            if (present(ddof)) ddof_value = ddof
            if (present(nan_policy)) then
                value = sem_low(x, ddof=ddof_value, nan_policy=nan_policy)
            else
                value = sem_low(x, ddof=ddof_value)
            end if
        case (stat_kind_variation)
            ddof_value = 0
            if (present(ddof)) ddof_value = ddof
            if (present(nan_policy)) then
                value = variation_low(x, ddof=ddof_value, nan_policy=nan_policy)
            else
                value = variation_low(x, ddof=ddof_value)
            end if
        case (stat_kind_iqr)
            rng_value = [25.0_dp, 75.0_dp]
            if (present(rng)) rng_value = rng
            scale_value = 1.0_dp
            if (present(scale)) scale_value = scale
            if (present(nan_policy)) then
                value = iqr_low(x, rng=rng_value, scale=scale_value, &
                    nan_policy=nan_policy)
            else
                value = iqr_low(x, rng=rng_value, scale=scale_value)
            end if
        case (stat_kind_mad)
            scale_value = 1.0_dp
            if (present(scale)) scale_value = scale
            if (present(nan_policy)) then
                value = mad_low(x, scale=scale_value, nan_policy=nan_policy)
            else
                value = mad_low(x, scale=scale_value)
            end if
        case (stat_kind_moment)
            order_value = 1
            if (present(order)) order_value = order
            if (present(center)) then
                if (present(nan_policy)) then
                    value = moment_low(x, order=order_value, center=center, &
                        nan_policy=nan_policy)
                else
                    value = moment_low(x, order=order_value, center=center)
                end if
            else
                if (present(nan_policy)) then
                    value = moment_low(x, order=order_value, nan_policy=nan_policy)
                else
                    value = moment_low(x, order=order_value)
                end if
            end if
        case (stat_kind_skew)
            bias_value = .true.
            if (present(bias)) bias_value = bias
            if (present(nan_policy)) then
                value = skew_low(x, bias=bias_value, nan_policy=nan_policy)
            else
                value = skew_low(x, bias=bias_value)
            end if
        case (stat_kind_kurtosis)
            fisher_value = .true.
            if (present(fisher)) fisher_value = fisher
            bias_value = .true.
            if (present(bias)) bias_value = bias
            if (present(nan_policy)) then
                value = kurtosis_low(x, fisher=fisher_value, bias=bias_value, &
                    nan_policy=nan_policy)
            else
                value = kurtosis_low(x, fisher=fisher_value, bias=bias_value)
            end if
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function evaluate_stat

    pure subroutine flatten_assumed_real(x, flat, dims)
        real(dp), intent(in) :: x(..) !! array to flatten in Fortran element order
        real(dp), allocatable, intent(out) :: flat(:) !! flattened observations
        integer, allocatable, intent(out) :: dims(:) !! source extents

        select rank (a => x)
        rank (5)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (6)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (7)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (8)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (9)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (10)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (11)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (12)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (13)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (14)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank (15)
            dims = shape(a)
            flat = reshape(a, [size(a)])
        rank default
            allocate(dims(0))
            allocate(flat(0))
        end select
    end subroutine flatten_assumed_real

end module scifort_descriptive_highrank_stats_engine
