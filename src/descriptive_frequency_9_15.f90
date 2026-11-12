! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_frequency_9_15
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_descriptive_frequency_engine, only : percentileofscore_scalar, &
        percentileofscore_vector, scoreatpercentile_all_assumed_rank, &
        scoreatpercentile_reduce_assumed_rank, scoreatpercentile_vector, &
        scoreatpercentile_vector_many
    implicit none
    private

    interface percentileofscore
        module procedure percentileofscore_score_rank9
        module procedure percentileofscore_score_rank10
        module procedure percentileofscore_score_rank11
        module procedure percentileofscore_score_rank12
        module procedure percentileofscore_score_rank13
        module procedure percentileofscore_score_rank14
        module procedure percentileofscore_score_rank15
    end interface percentileofscore

    interface scoreatpercentile
        module procedure scoreatpercentile_all_scalar_rank9
        module procedure scoreatpercentile_axis_scalar_rank9
        module procedure scoreatpercentile_all_vector_rank9
        module procedure scoreatpercentile_axis_vector_rank9
        module procedure scoreatpercentile_all_scalar_rank10
        module procedure scoreatpercentile_axis_scalar_rank10
        module procedure scoreatpercentile_all_vector_rank10
        module procedure scoreatpercentile_axis_vector_rank10
        module procedure scoreatpercentile_all_scalar_rank11
        module procedure scoreatpercentile_axis_scalar_rank11
        module procedure scoreatpercentile_all_vector_rank11
        module procedure scoreatpercentile_axis_vector_rank11
        module procedure scoreatpercentile_all_scalar_rank12
        module procedure scoreatpercentile_axis_scalar_rank12
        module procedure scoreatpercentile_all_vector_rank12
        module procedure scoreatpercentile_axis_vector_rank12
        module procedure scoreatpercentile_all_scalar_rank13
        module procedure scoreatpercentile_axis_scalar_rank13
        module procedure scoreatpercentile_all_vector_rank13
        module procedure scoreatpercentile_axis_vector_rank13
        module procedure scoreatpercentile_all_scalar_rank14
        module procedure scoreatpercentile_axis_scalar_rank14
        module procedure scoreatpercentile_all_vector_rank14
        module procedure scoreatpercentile_axis_vector_rank14
        module procedure scoreatpercentile_all_scalar_rank15
        module procedure scoreatpercentile_axis_scalar_rank15
        module procedure scoreatpercentile_all_vector_rank15
        module procedure scoreatpercentile_axis_vector_rank15
    end interface scoreatpercentile

    interface scoreatpercentile_all
        module procedure scoreatpercentile_all_scalar_rank9
        module procedure scoreatpercentile_all_vector_rank9
        module procedure scoreatpercentile_all_scalar_rank10
        module procedure scoreatpercentile_all_vector_rank10
        module procedure scoreatpercentile_all_scalar_rank11
        module procedure scoreatpercentile_all_vector_rank11
        module procedure scoreatpercentile_all_scalar_rank12
        module procedure scoreatpercentile_all_vector_rank12
        module procedure scoreatpercentile_all_scalar_rank13
        module procedure scoreatpercentile_all_vector_rank13
        module procedure scoreatpercentile_all_scalar_rank14
        module procedure scoreatpercentile_all_vector_rank14
        module procedure scoreatpercentile_all_scalar_rank15
        module procedure scoreatpercentile_all_vector_rank15
    end interface scoreatpercentile_all

    public :: percentileofscore
    public :: scoreatpercentile
    public :: scoreatpercentile_all

contains

    pure function percentileofscore_score_rank9(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank9

    pure function scoreatpercentile_all_scalar_rank9(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank9

    pure function scoreatpercentile_all_vector_rank9(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank9

    pure function scoreatpercentile_axis_scalar_rank9(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: out_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank9

    pure function scoreatpercentile_axis_vector_rank9(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: out_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank9

    pure function percentileofscore_score_rank10(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank10

    pure function scoreatpercentile_all_scalar_rank10(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank10

    pure function scoreatpercentile_all_vector_rank10(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank10

    pure function scoreatpercentile_axis_scalar_rank10(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: out_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank10

    pure function scoreatpercentile_axis_vector_rank10(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank10

    pure function percentileofscore_score_rank11(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank11

    pure function scoreatpercentile_all_scalar_rank11(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank11

    pure function scoreatpercentile_all_vector_rank11(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank11

    pure function scoreatpercentile_axis_scalar_rank11(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank11

    pure function scoreatpercentile_axis_vector_rank11(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank11

    pure function percentileofscore_score_rank12(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank12

    pure function scoreatpercentile_all_scalar_rank12(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank12

    pure function scoreatpercentile_all_vector_rank12(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank12

    pure function scoreatpercentile_axis_scalar_rank12(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank12

    pure function scoreatpercentile_axis_vector_rank12(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank12

    pure function percentileofscore_score_rank13(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank13

    pure function scoreatpercentile_all_scalar_rank13(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank13

    pure function scoreatpercentile_all_vector_rank13(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank13

    pure function scoreatpercentile_axis_scalar_rank13(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank13

    pure function scoreatpercentile_axis_vector_rank13(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank13

    pure function percentileofscore_score_rank14(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank14

    pure function scoreatpercentile_all_scalar_rank14(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank14

    pure function scoreatpercentile_all_vector_rank14(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank14

    pure function scoreatpercentile_axis_scalar_rank14(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank14

    pure function scoreatpercentile_axis_vector_rank14(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank14

    pure function percentileofscore_score_rank15(a, score, kind, nan_policy) result(values)
        real(dp), intent(in) :: a(:) !! one-dimensional reference sample
        real(dp), intent(in) :: score(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! scores whose percentile ranks are requested
        character(len=*), intent(in), optional :: kind !! rank/weak/strict/mean
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        real(dp), allocatable :: flat_score(:)
        real(dp), allocatable :: flat_values(:)

        flat_score = reshape(score, [size(score)])
        flat_values = percentileofscore_vector(a, flat_score, kind, nan_policy)
        values = reshape(flat_values, shape(score))
    end function percentileofscore_score_rank15

    pure function scoreatpercentile_all_scalar_rank15(x, per, limit, interpolation_method) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp) :: value

        real(dp), allocatable :: temporary(:)

        temporary = scoreatpercentile_all_assumed_rank(x, [per], limit, interpolation_method)
        value = temporary(1)
    end function scoreatpercentile_all_scalar_rank15

    pure function scoreatpercentile_all_vector_rank15(x, per, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:)

        values = scoreatpercentile_all_assumed_rank(x, per, limit, interpolation_method)
    end function scoreatpercentile_all_vector_rank15

    pure function scoreatpercentile_axis_scalar_rank15(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per !! percentile in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, [per], reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) then
            out_shape = 0
        end if
        values = reshape(reduced, out_shape)
    end function scoreatpercentile_axis_scalar_rank15

    pure function scoreatpercentile_axis_vector_rank15(x, per, axis, limit, interpolation_method) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: per(:) !! percentiles in [0, 100]
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: limit(2) !! closed lower/upper filtering interval
        character(len=*), intent(in), optional :: interpolation_method !! fraction/lower/higher
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call scoreatpercentile_reduce_assumed_rank(x, axis, per, reduced, out_shape, ok, &
            limit, interpolation_method)
        if (.not. ok) out_shape = 0
        result_shape(1) = size(per)
        result_shape(2:) = out_shape
        values = reshape(reduced, result_shape)
    end function scoreatpercentile_axis_vector_rank15

end module scifort_descriptive_frequency_9_15
