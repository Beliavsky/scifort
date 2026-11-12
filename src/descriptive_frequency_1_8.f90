! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_frequency_1_8
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_descriptive_frequency_engine, only : percentileofscore_scalar, &
        percentileofscore_vector, scoreatpercentile_all_assumed_rank, &
        scoreatpercentile_reduce_assumed_rank, scoreatpercentile_vector, &
        scoreatpercentile_vector_many
    implicit none
    private

    interface percentileofscore
        module procedure percentileofscore_score_scalar
        module procedure percentileofscore_score_rank1
        module procedure percentileofscore_score_rank2
        module procedure percentileofscore_score_rank3
        module procedure percentileofscore_score_rank4
        module procedure percentileofscore_score_rank5
        module procedure percentileofscore_score_rank6
        module procedure percentileofscore_score_rank7
        module procedure percentileofscore_score_rank8
    end interface percentileofscore

    interface scoreatpercentile
        module procedure scoreatpercentile_all_scalar_rank1
        module procedure scoreatpercentile_axis_scalar_rank1
        module procedure scoreatpercentile_all_vector_rank1
        module procedure scoreatpercentile_axis_vector_rank1
        module procedure scoreatpercentile_all_scalar_rank2
        module procedure scoreatpercentile_axis_scalar_rank2
        module procedure scoreatpercentile_all_vector_rank2
        module procedure scoreatpercentile_axis_vector_rank2
        module procedure scoreatpercentile_all_scalar_rank3
        module procedure scoreatpercentile_axis_scalar_rank3
        module procedure scoreatpercentile_all_vector_rank3
        module procedure scoreatpercentile_axis_vector_rank3
        module procedure scoreatpercentile_all_scalar_rank4
        module procedure scoreatpercentile_axis_scalar_rank4
        module procedure scoreatpercentile_all_vector_rank4
        module procedure scoreatpercentile_axis_vector_rank4
        module procedure scoreatpercentile_all_scalar_rank5
        module procedure scoreatpercentile_axis_scalar_rank5
        module procedure scoreatpercentile_all_vector_rank5
        module procedure scoreatpercentile_axis_vector_rank5
        module procedure scoreatpercentile_all_scalar_rank6
        module procedure scoreatpercentile_axis_scalar_rank6
        module procedure scoreatpercentile_all_vector_rank6
        module procedure scoreatpercentile_axis_vector_rank6
        module procedure scoreatpercentile_all_scalar_rank7
        module procedure scoreatpercentile_axis_scalar_rank7
        module procedure scoreatpercentile_all_vector_rank7
        module procedure scoreatpercentile_axis_vector_rank7
        module procedure scoreatpercentile_all_scalar_rank8
        module procedure scoreatpercentile_axis_scalar_rank8
        module procedure scoreatpercentile_all_vector_rank8
        module procedure scoreatpercentile_axis_vector_rank8
    end interface scoreatpercentile

    interface scoreatpercentile_all
        module procedure scoreatpercentile_all_scalar_rank1
        module procedure scoreatpercentile_all_vector_rank1
        module procedure scoreatpercentile_all_scalar_rank2
        module procedure scoreatpercentile_all_vector_rank2
        module procedure scoreatpercentile_all_scalar_rank3
        module procedure scoreatpercentile_all_vector_rank3
        module procedure scoreatpercentile_all_scalar_rank4
        module procedure scoreatpercentile_all_vector_rank4
        module procedure scoreatpercentile_all_scalar_rank5
        module procedure scoreatpercentile_all_vector_rank5
        module procedure scoreatpercentile_all_scalar_rank6
        module procedure scoreatpercentile_all_vector_rank6
        module procedure scoreatpercentile_all_scalar_rank7
        module procedure scoreatpercentile_all_vector_rank7
        module procedure scoreatpercentile_all_scalar_rank8
        module procedure scoreatpercentile_all_vector_rank8
    end interface scoreatpercentile_all

    public :: percentileofscore
    public :: scoreatpercentile
    public :: scoreatpercentile_all

contains

    pure function percentileofscore_score_scalar(a, score, kind, nan_policy) result(value)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score !! score whose percentile rank is requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = percentileofscore_scalar(a, score, kind, nan_policy)
    end function percentileofscore_score_scalar

    pure function percentileofscore_score_rank1(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank1

    pure function scoreatpercentile_all_scalar_rank1(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank1

    pure function scoreatpercentile_all_vector_rank1(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank1

    pure function scoreatpercentile_axis_scalar_rank1(x, per, axis, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        if (axis /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = scoreatpercentile_vector(x, per, limit, interpolation_method)
    end function scoreatpercentile_axis_scalar_rank1

    pure function scoreatpercentile_axis_vector_rank1(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        if (axis /= 1) then
            allocate(values(0))
            return
        end if
        values = scoreatpercentile_vector_many(x, per, limit, interpolation_method)
    end function scoreatpercentile_axis_vector_rank1

    pure function percentileofscore_score_rank2(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank2

    pure function scoreatpercentile_all_scalar_rank2(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank2

    pure function scoreatpercentile_all_vector_rank2(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank2

    pure function scoreatpercentile_axis_scalar_rank2(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        integer :: out_shape(1)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank2

    pure function scoreatpercentile_axis_vector_rank2(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :)

        integer :: out_shape(1)
        integer :: result_shape(2)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank2

    pure function percentileofscore_score_rank3(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank3

    pure function scoreatpercentile_all_scalar_rank3(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank3

    pure function scoreatpercentile_all_vector_rank3(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank3

    pure function scoreatpercentile_axis_scalar_rank3(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :)

        integer :: out_shape(2)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank3

    pure function scoreatpercentile_axis_vector_rank3(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :)

        integer :: out_shape(2)
        integer :: result_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank3

    pure function percentileofscore_score_rank4(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank4

    pure function scoreatpercentile_all_scalar_rank4(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank4

    pure function scoreatpercentile_all_vector_rank4(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank4

    pure function scoreatpercentile_axis_scalar_rank4(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :)

        integer :: out_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank4

    pure function scoreatpercentile_axis_vector_rank4(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :)

        integer :: out_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank4

    pure function percentileofscore_score_rank5(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank5

    pure function scoreatpercentile_all_scalar_rank5(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank5

    pure function scoreatpercentile_all_vector_rank5(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank5

    pure function scoreatpercentile_axis_scalar_rank5(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :)

        integer :: out_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank5

    pure function scoreatpercentile_axis_vector_rank5(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: out_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank5

    pure function percentileofscore_score_rank6(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank6

    pure function scoreatpercentile_all_scalar_rank6(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank6

    pure function scoreatpercentile_all_vector_rank6(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank6

    pure function scoreatpercentile_axis_scalar_rank6(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: out_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank6

    pure function scoreatpercentile_axis_vector_rank6(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: out_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank6

    pure function percentileofscore_score_rank7(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank7

    pure function scoreatpercentile_all_scalar_rank7(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank7

    pure function scoreatpercentile_all_vector_rank7(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank7

    pure function scoreatpercentile_axis_scalar_rank7(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: out_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank7

    pure function scoreatpercentile_axis_vector_rank7(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: out_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank7

    pure function percentileofscore_score_rank8(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank8

    pure function scoreatpercentile_all_scalar_rank8(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank8

    pure function scoreatpercentile_all_vector_rank8(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank8

    pure function scoreatpercentile_axis_scalar_rank8(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: out_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank8

    pure function scoreatpercentile_axis_vector_rank8(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: out_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank8

end module scifort_descriptive_frequency_1_8
