! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions_high_rank_weighted_11_15
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_special_reductions_high_rank_weighted_support, only : &
        complex_log_from_phase, complex_nan, drop_dimension_shape, &
        flatten_complex_weights, flatten_real_weights, &
        weighted_complex_prepared, weighted_real_prepared
    implicit none
    private

    integer, parameter :: axes_11(11) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    integer, parameter :: axes_12(12) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    integer, parameter :: axes_13(13) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    integer, parameter :: axes_14(14) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    integer, parameter :: axes_15(15) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

    public :: logsumexp
    public :: logsumexp_keepdims
    public :: logsumexp_signed
    public :: logsumexp_signed_keepdims

    interface logsumexp
        module procedure logsumexp_11d_all_weighted
        module procedure logsumexp_11d_dim_weighted
        module procedure logsumexp_real_11d_all_complex_b
        module procedure logsumexp_real_11d_dim_complex_b
        module procedure logsumexp_complex_11d_all_weighted
        module procedure logsumexp_complex_11d_dim_weighted
        module procedure logsumexp_complex_11d_all_real_b
        module procedure logsumexp_complex_11d_dim_real_b
        module procedure logsumexp_12d_all_weighted
        module procedure logsumexp_12d_dim_weighted
        module procedure logsumexp_real_12d_all_complex_b
        module procedure logsumexp_real_12d_dim_complex_b
        module procedure logsumexp_complex_12d_all_weighted
        module procedure logsumexp_complex_12d_dim_weighted
        module procedure logsumexp_complex_12d_all_real_b
        module procedure logsumexp_complex_12d_dim_real_b
        module procedure logsumexp_13d_all_weighted
        module procedure logsumexp_13d_dim_weighted
        module procedure logsumexp_real_13d_all_complex_b
        module procedure logsumexp_real_13d_dim_complex_b
        module procedure logsumexp_complex_13d_all_weighted
        module procedure logsumexp_complex_13d_dim_weighted
        module procedure logsumexp_complex_13d_all_real_b
        module procedure logsumexp_complex_13d_dim_real_b
        module procedure logsumexp_14d_all_weighted
        module procedure logsumexp_14d_dim_weighted
        module procedure logsumexp_real_14d_all_complex_b
        module procedure logsumexp_real_14d_dim_complex_b
        module procedure logsumexp_complex_14d_all_weighted
        module procedure logsumexp_complex_14d_dim_weighted
        module procedure logsumexp_complex_14d_all_real_b
        module procedure logsumexp_complex_14d_dim_real_b
        module procedure logsumexp_15d_all_weighted
        module procedure logsumexp_15d_dim_weighted
        module procedure logsumexp_real_15d_all_complex_b
        module procedure logsumexp_real_15d_dim_complex_b
        module procedure logsumexp_complex_15d_all_weighted
        module procedure logsumexp_complex_15d_dim_weighted
        module procedure logsumexp_complex_15d_all_real_b
        module procedure logsumexp_complex_15d_dim_real_b
    end interface logsumexp

    interface logsumexp_signed
        module procedure logsumexp_signed_11d_all_weighted
        module procedure logsumexp_signed_11d_dim_weighted
        module procedure logsumexp_signed_real_11d_all_complex_b
        module procedure logsumexp_signed_real_11d_dim_complex_b
        module procedure logsumexp_signed_complex_11d_all_weighted
        module procedure logsumexp_signed_complex_11d_dim_weighted
        module procedure logsumexp_signed_complex_11d_all_real_b
        module procedure logsumexp_signed_complex_11d_dim_real_b
        module procedure logsumexp_signed_12d_all_weighted
        module procedure logsumexp_signed_12d_dim_weighted
        module procedure logsumexp_signed_real_12d_all_complex_b
        module procedure logsumexp_signed_real_12d_dim_complex_b
        module procedure logsumexp_signed_complex_12d_all_weighted
        module procedure logsumexp_signed_complex_12d_dim_weighted
        module procedure logsumexp_signed_complex_12d_all_real_b
        module procedure logsumexp_signed_complex_12d_dim_real_b
        module procedure logsumexp_signed_13d_all_weighted
        module procedure logsumexp_signed_13d_dim_weighted
        module procedure logsumexp_signed_real_13d_all_complex_b
        module procedure logsumexp_signed_real_13d_dim_complex_b
        module procedure logsumexp_signed_complex_13d_all_weighted
        module procedure logsumexp_signed_complex_13d_dim_weighted
        module procedure logsumexp_signed_complex_13d_all_real_b
        module procedure logsumexp_signed_complex_13d_dim_real_b
        module procedure logsumexp_signed_14d_all_weighted
        module procedure logsumexp_signed_14d_dim_weighted
        module procedure logsumexp_signed_real_14d_all_complex_b
        module procedure logsumexp_signed_real_14d_dim_complex_b
        module procedure logsumexp_signed_complex_14d_all_weighted
        module procedure logsumexp_signed_complex_14d_dim_weighted
        module procedure logsumexp_signed_complex_14d_all_real_b
        module procedure logsumexp_signed_complex_14d_dim_real_b
        module procedure logsumexp_signed_15d_all_weighted
        module procedure logsumexp_signed_15d_dim_weighted
        module procedure logsumexp_signed_real_15d_all_complex_b
        module procedure logsumexp_signed_real_15d_dim_complex_b
        module procedure logsumexp_signed_complex_15d_all_weighted
        module procedure logsumexp_signed_complex_15d_dim_weighted
        module procedure logsumexp_signed_complex_15d_all_real_b
        module procedure logsumexp_signed_complex_15d_dim_real_b
    end interface logsumexp_signed

    interface logsumexp_keepdims
        module procedure logsumexp_keepdims_11d_weighted
        module procedure logsumexp_keepdims_real_11d_complex_b
        module procedure logsumexp_keepdims_complex_11d_weighted
        module procedure logsumexp_keepdims_complex_11d_real_b
        module procedure logsumexp_keepdims_12d_weighted
        module procedure logsumexp_keepdims_real_12d_complex_b
        module procedure logsumexp_keepdims_complex_12d_weighted
        module procedure logsumexp_keepdims_complex_12d_real_b
        module procedure logsumexp_keepdims_13d_weighted
        module procedure logsumexp_keepdims_real_13d_complex_b
        module procedure logsumexp_keepdims_complex_13d_weighted
        module procedure logsumexp_keepdims_complex_13d_real_b
        module procedure logsumexp_keepdims_14d_weighted
        module procedure logsumexp_keepdims_real_14d_complex_b
        module procedure logsumexp_keepdims_complex_14d_weighted
        module procedure logsumexp_keepdims_complex_14d_real_b
        module procedure logsumexp_keepdims_15d_weighted
        module procedure logsumexp_keepdims_real_15d_complex_b
        module procedure logsumexp_keepdims_complex_15d_weighted
        module procedure logsumexp_keepdims_complex_15d_real_b
    end interface logsumexp_keepdims

    interface logsumexp_signed_keepdims
        module procedure logsumexp_signed_keepdims_11d_weighted
        module procedure logsumexp_signed_keepdims_real_11d_complex_b
        module procedure logsumexp_signed_keepdims_complex_11d_weighted
        module procedure logsumexp_signed_keepdims_complex_11d_real_b
        module procedure logsumexp_signed_keepdims_12d_weighted
        module procedure logsumexp_signed_keepdims_real_12d_complex_b
        module procedure logsumexp_signed_keepdims_complex_12d_weighted
        module procedure logsumexp_signed_keepdims_complex_12d_real_b
        module procedure logsumexp_signed_keepdims_13d_weighted
        module procedure logsumexp_signed_keepdims_real_13d_complex_b
        module procedure logsumexp_signed_keepdims_complex_13d_weighted
        module procedure logsumexp_signed_keepdims_complex_13d_real_b
        module procedure logsumexp_signed_keepdims_14d_weighted
        module procedure logsumexp_signed_keepdims_real_14d_complex_b
        module procedure logsumexp_signed_keepdims_complex_14d_weighted
        module procedure logsumexp_signed_keepdims_complex_14d_real_b
        module procedure logsumexp_signed_keepdims_15d_weighted
        module procedure logsumexp_signed_keepdims_real_15d_complex_b
        module procedure logsumexp_signed_keepdims_complex_15d_weighted
        module procedure logsumexp_signed_keepdims_complex_15d_real_b
    end interface logsumexp_signed_keepdims

contains

    ! Weighted rank-11 wrappers.
    pure subroutine logsumexp_signed_keepdims_11d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(11)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_11d_weighted

    pure function logsumexp_keepdims_11d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_11d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_11d_weighted

    pure subroutine logsumexp_signed_11d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_11d_weighted(x, axes_11, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_11d_all_weighted

    pure function logsumexp_11d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_11d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_11d_all_weighted

    pure subroutine logsumexp_signed_11d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        call logsumexp_signed_keepdims_11d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_11d_dim_weighted

    pure function logsumexp_11d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_11d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_11d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_11d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(11)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_11d_complex_b

    pure function logsumexp_keepdims_real_11d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_real_11d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_11d_complex_b

    pure subroutine logsumexp_signed_real_11d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_11d_complex_b(x, axes_11, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_11d_all_complex_b

    pure function logsumexp_real_11d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_11d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_11d_all_complex_b

    pure subroutine logsumexp_signed_real_11d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_11d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_11d_dim_complex_b

    pure function logsumexp_real_11d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_real_11d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_11d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_11d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(11)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_11d_weighted

    pure function logsumexp_keepdims_complex_11d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_11d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_11d_weighted

    pure subroutine logsumexp_signed_complex_11d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_11d_weighted(x, axes_11, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_11d_all_weighted

    pure function logsumexp_complex_11d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_11d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_11d_all_weighted

    pure subroutine logsumexp_signed_complex_11d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_11d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_11d_dim_weighted

    pure function logsumexp_complex_11d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_complex_11d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_11d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_11d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(11)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_11d_real_b

    pure function logsumexp_keepdims_complex_11d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_11d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_11d_real_b

    pure subroutine logsumexp_signed_complex_11d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_11d_real_b(x, axes_11, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_11d_all_real_b

    pure function logsumexp_complex_11d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_11d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_11d_all_real_b

    pure subroutine logsumexp_signed_complex_11d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :)
        integer :: outshape(10)

        if (dim < 1 .or. dim > 11) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_11d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_11d_dim_real_b

    pure function logsumexp_complex_11d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 11
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_complex_11d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_11d_dim_real_b

    ! Weighted rank-12 wrappers.
    pure subroutine logsumexp_signed_keepdims_12d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(12)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_12d_weighted

    pure function logsumexp_keepdims_12d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_12d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_12d_weighted

    pure subroutine logsumexp_signed_12d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_12d_weighted(x, axes_12, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_12d_all_weighted

    pure function logsumexp_12d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_12d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_12d_all_weighted

    pure subroutine logsumexp_signed_12d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_12d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_12d_dim_weighted

    pure function logsumexp_12d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_12d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_12d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_12d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(12)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_12d_complex_b

    pure function logsumexp_keepdims_real_12d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_12d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_12d_complex_b

    pure subroutine logsumexp_signed_real_12d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_12d_complex_b(x, axes_12, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_12d_all_complex_b

    pure function logsumexp_real_12d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_12d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_12d_all_complex_b

    pure subroutine logsumexp_signed_real_12d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_12d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_12d_dim_complex_b

    pure function logsumexp_real_12d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_real_12d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_12d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_12d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(12)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_12d_weighted

    pure function logsumexp_keepdims_complex_12d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_12d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_12d_weighted

    pure subroutine logsumexp_signed_complex_12d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_12d_weighted(x, axes_12, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_12d_all_weighted

    pure function logsumexp_complex_12d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_12d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_12d_all_weighted

    pure subroutine logsumexp_signed_complex_12d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_12d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_12d_dim_weighted

    pure function logsumexp_complex_12d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_complex_12d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_12d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_12d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(12)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_12d_real_b

    pure function logsumexp_keepdims_complex_12d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_12d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_12d_real_b

    pure subroutine logsumexp_signed_complex_12d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_12d_real_b(x, axes_12, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_12d_all_real_b

    pure function logsumexp_complex_12d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_12d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_12d_all_real_b

    pure subroutine logsumexp_signed_complex_12d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :)
        integer :: outshape(11)

        if (dim < 1 .or. dim > 12) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_12d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_12d_dim_real_b

    pure function logsumexp_complex_12d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 12
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_complex_12d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_12d_dim_real_b

    ! Weighted rank-13 wrappers.
    pure subroutine logsumexp_signed_keepdims_13d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(13)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_13d_weighted

    pure function logsumexp_keepdims_13d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_13d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_13d_weighted

    pure subroutine logsumexp_signed_13d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_13d_weighted(x, axes_13, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_13d_all_weighted

    pure function logsumexp_13d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_13d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_13d_all_weighted

    pure subroutine logsumexp_signed_13d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_13d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_13d_dim_weighted

    pure function logsumexp_13d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_13d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_13d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_13d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(13)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_13d_complex_b

    pure function logsumexp_keepdims_real_13d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_13d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_13d_complex_b

    pure subroutine logsumexp_signed_real_13d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_13d_complex_b(x, axes_13, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_13d_all_complex_b

    pure function logsumexp_real_13d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_13d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_13d_all_complex_b

    pure subroutine logsumexp_signed_real_13d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        call logsumexp_signed_keepdims_real_13d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_13d_dim_complex_b

    pure function logsumexp_real_13d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_real_13d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_13d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_13d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(13)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_13d_weighted

    pure function logsumexp_keepdims_complex_13d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_13d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_13d_weighted

    pure subroutine logsumexp_signed_complex_13d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_13d_weighted(x, axes_13, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_13d_all_weighted

    pure function logsumexp_complex_13d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_13d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_13d_all_weighted

    pure subroutine logsumexp_signed_complex_13d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_13d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_13d_dim_weighted

    pure function logsumexp_complex_13d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_complex_13d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_13d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_13d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(13)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_13d_real_b

    pure function logsumexp_keepdims_complex_13d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_13d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_13d_real_b

    pure subroutine logsumexp_signed_complex_13d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_13d_real_b(x, axes_13, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_13d_all_real_b

    pure function logsumexp_complex_13d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_13d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_13d_all_real_b

    pure subroutine logsumexp_signed_complex_13d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        integer :: outshape(12)

        if (dim < 1 .or. dim > 13) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_13d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_13d_dim_real_b

    pure function logsumexp_complex_13d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 13
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_complex_13d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_13d_dim_real_b

    ! Weighted rank-14 wrappers.
    pure subroutine logsumexp_signed_keepdims_14d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(14)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_14d_weighted

    pure function logsumexp_keepdims_14d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_14d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_14d_weighted

    pure subroutine logsumexp_signed_14d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_14d_weighted(x, axes_14, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_14d_all_weighted

    pure function logsumexp_14d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_14d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_14d_all_weighted

    pure subroutine logsumexp_signed_14d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_14d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_14d_dim_weighted

    pure function logsumexp_14d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_14d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_14d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_14d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(14)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_14d_complex_b

    pure function logsumexp_keepdims_real_14d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_14d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_14d_complex_b

    pure subroutine logsumexp_signed_real_14d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_14d_complex_b(x, axes_14, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_14d_all_complex_b

    pure function logsumexp_real_14d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_14d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_14d_all_complex_b

    pure subroutine logsumexp_signed_real_14d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_14d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_14d_dim_complex_b

    pure function logsumexp_real_14d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_real_14d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_14d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_14d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(14)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_14d_weighted

    pure function logsumexp_keepdims_complex_14d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_14d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_14d_weighted

    pure subroutine logsumexp_signed_complex_14d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_14d_weighted(x, axes_14, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_14d_all_weighted

    pure function logsumexp_complex_14d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_14d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_14d_all_weighted

    pure subroutine logsumexp_signed_complex_14d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_14d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_14d_dim_weighted

    pure function logsumexp_complex_14d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_complex_14d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_14d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_14d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(14)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_14d_real_b

    pure function logsumexp_keepdims_complex_14d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_14d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_14d_real_b

    pure subroutine logsumexp_signed_complex_14d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_14d_real_b(x, axes_14, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_14d_all_real_b

    pure function logsumexp_complex_14d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_14d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_14d_all_real_b

    pure subroutine logsumexp_signed_complex_14d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(13)

        if (dim < 1 .or. dim > 14) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_14d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_14d_dim_real_b

    pure function logsumexp_complex_14d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 14
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :)

        call logsumexp_signed_complex_14d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_14d_dim_real_b

    ! Weighted rank-15 wrappers.
    pure subroutine logsumexp_signed_keepdims_15d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(15)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_15d_weighted

    pure function logsumexp_keepdims_15d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_15d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14), &
            & size(logabs, 15)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_15d_weighted

    pure subroutine logsumexp_signed_15d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_15d_weighted(x, axes_15, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_15d_all_weighted

    pure function logsumexp_15d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_15d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_15d_all_weighted

    pure subroutine logsumexp_signed_15d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_15d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_15d_dim_weighted

    pure function logsumexp_15d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_15d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_15d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_15d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(15)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_15d_complex_b

    pure function logsumexp_keepdims_real_15d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_15d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14), &
            & size(logabs, 15)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_15d_complex_b

    pure subroutine logsumexp_signed_real_15d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_15d_complex_b(x, axes_15, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_15d_all_complex_b

    pure function logsumexp_real_15d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_15d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_15d_all_complex_b

    pure subroutine logsumexp_signed_real_15d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_15d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_15d_dim_complex_b

    pure function logsumexp_real_15d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_real_15d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_15d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_15d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(15)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_15d_weighted

    pure function logsumexp_keepdims_complex_15d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_15d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14), &
            & size(logabs, 15)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_15d_weighted

    pure subroutine logsumexp_signed_complex_15d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_15d_weighted(x, axes_15, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_15d_all_weighted

    pure function logsumexp_complex_15d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_15d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_15d_all_weighted

    pure subroutine logsumexp_signed_complex_15d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_15d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_15d_dim_weighted

    pure function logsumexp_complex_15d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_complex_15d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_15d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_15d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(15)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14), &
            & outshape(15)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_15d_real_b

    pure function logsumexp_keepdims_complex_15d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_15d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14), &
            & size(logabs, 15)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_15d_real_b

    pure subroutine logsumexp_signed_complex_15d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_15d_real_b(x, axes_15, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_15d_all_real_b

    pure function logsumexp_complex_15d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_15d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_15d_all_real_b

    pure subroutine logsumexp_signed_complex_15d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :, :, :, :, :, :)
        integer :: outshape(14)

        if (dim < 1 .or. dim > 15) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                & 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_15d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10), outshape(11), &
            & outshape(12), outshape(13), outshape(14)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_15d_dim_real_b

    pure function logsumexp_complex_15d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 15
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :, :, :, :, :)

        call logsumexp_signed_complex_15d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10), &
            & size(logabs, 11), size(logabs, 12), &
            & size(logabs, 13), size(logabs, 14)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_15d_dim_real_b

end module scifort_special_reductions_high_rank_weighted_11_15
