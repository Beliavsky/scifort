! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions_high_rank_weighted_6_10
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_special_reductions_high_rank_weighted_support, only : &
        complex_log_from_phase, complex_nan, drop_dimension_shape, &
        flatten_complex_weights, flatten_real_weights, &
        weighted_complex_prepared, weighted_real_prepared
    implicit none
    private

    integer, parameter :: axes_6(6) = [1, 2, 3, 4, 5, 6]
    integer, parameter :: axes_7(7) = [1, 2, 3, 4, 5, 6, 7]
    integer, parameter :: axes_8(8) = [1, 2, 3, 4, 5, 6, 7, 8]
    integer, parameter :: axes_9(9) = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    integer, parameter :: axes_10(10) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    public :: logsumexp
    public :: logsumexp_keepdims
    public :: logsumexp_signed
    public :: logsumexp_signed_keepdims

    interface logsumexp
        module procedure logsumexp_6d_all_weighted
        module procedure logsumexp_6d_dim_weighted
        module procedure logsumexp_real_6d_all_complex_b
        module procedure logsumexp_real_6d_dim_complex_b
        module procedure logsumexp_complex_6d_all_weighted
        module procedure logsumexp_complex_6d_dim_weighted
        module procedure logsumexp_complex_6d_all_real_b
        module procedure logsumexp_complex_6d_dim_real_b
        module procedure logsumexp_7d_all_weighted
        module procedure logsumexp_7d_dim_weighted
        module procedure logsumexp_real_7d_all_complex_b
        module procedure logsumexp_real_7d_dim_complex_b
        module procedure logsumexp_complex_7d_all_weighted
        module procedure logsumexp_complex_7d_dim_weighted
        module procedure logsumexp_complex_7d_all_real_b
        module procedure logsumexp_complex_7d_dim_real_b
        module procedure logsumexp_8d_all_weighted
        module procedure logsumexp_8d_dim_weighted
        module procedure logsumexp_real_8d_all_complex_b
        module procedure logsumexp_real_8d_dim_complex_b
        module procedure logsumexp_complex_8d_all_weighted
        module procedure logsumexp_complex_8d_dim_weighted
        module procedure logsumexp_complex_8d_all_real_b
        module procedure logsumexp_complex_8d_dim_real_b
        module procedure logsumexp_9d_all_weighted
        module procedure logsumexp_9d_dim_weighted
        module procedure logsumexp_real_9d_all_complex_b
        module procedure logsumexp_real_9d_dim_complex_b
        module procedure logsumexp_complex_9d_all_weighted
        module procedure logsumexp_complex_9d_dim_weighted
        module procedure logsumexp_complex_9d_all_real_b
        module procedure logsumexp_complex_9d_dim_real_b
        module procedure logsumexp_10d_all_weighted
        module procedure logsumexp_10d_dim_weighted
        module procedure logsumexp_real_10d_all_complex_b
        module procedure logsumexp_real_10d_dim_complex_b
        module procedure logsumexp_complex_10d_all_weighted
        module procedure logsumexp_complex_10d_dim_weighted
        module procedure logsumexp_complex_10d_all_real_b
        module procedure logsumexp_complex_10d_dim_real_b
    end interface logsumexp

    interface logsumexp_signed
        module procedure logsumexp_signed_6d_all_weighted
        module procedure logsumexp_signed_6d_dim_weighted
        module procedure logsumexp_signed_real_6d_all_complex_b
        module procedure logsumexp_signed_real_6d_dim_complex_b
        module procedure logsumexp_signed_complex_6d_all_weighted
        module procedure logsumexp_signed_complex_6d_dim_weighted
        module procedure logsumexp_signed_complex_6d_all_real_b
        module procedure logsumexp_signed_complex_6d_dim_real_b
        module procedure logsumexp_signed_7d_all_weighted
        module procedure logsumexp_signed_7d_dim_weighted
        module procedure logsumexp_signed_real_7d_all_complex_b
        module procedure logsumexp_signed_real_7d_dim_complex_b
        module procedure logsumexp_signed_complex_7d_all_weighted
        module procedure logsumexp_signed_complex_7d_dim_weighted
        module procedure logsumexp_signed_complex_7d_all_real_b
        module procedure logsumexp_signed_complex_7d_dim_real_b
        module procedure logsumexp_signed_8d_all_weighted
        module procedure logsumexp_signed_8d_dim_weighted
        module procedure logsumexp_signed_real_8d_all_complex_b
        module procedure logsumexp_signed_real_8d_dim_complex_b
        module procedure logsumexp_signed_complex_8d_all_weighted
        module procedure logsumexp_signed_complex_8d_dim_weighted
        module procedure logsumexp_signed_complex_8d_all_real_b
        module procedure logsumexp_signed_complex_8d_dim_real_b
        module procedure logsumexp_signed_9d_all_weighted
        module procedure logsumexp_signed_9d_dim_weighted
        module procedure logsumexp_signed_real_9d_all_complex_b
        module procedure logsumexp_signed_real_9d_dim_complex_b
        module procedure logsumexp_signed_complex_9d_all_weighted
        module procedure logsumexp_signed_complex_9d_dim_weighted
        module procedure logsumexp_signed_complex_9d_all_real_b
        module procedure logsumexp_signed_complex_9d_dim_real_b
        module procedure logsumexp_signed_10d_all_weighted
        module procedure logsumexp_signed_10d_dim_weighted
        module procedure logsumexp_signed_real_10d_all_complex_b
        module procedure logsumexp_signed_real_10d_dim_complex_b
        module procedure logsumexp_signed_complex_10d_all_weighted
        module procedure logsumexp_signed_complex_10d_dim_weighted
        module procedure logsumexp_signed_complex_10d_all_real_b
        module procedure logsumexp_signed_complex_10d_dim_real_b
    end interface logsumexp_signed

    interface logsumexp_keepdims
        module procedure logsumexp_keepdims_6d_weighted
        module procedure logsumexp_keepdims_real_6d_complex_b
        module procedure logsumexp_keepdims_complex_6d_weighted
        module procedure logsumexp_keepdims_complex_6d_real_b
        module procedure logsumexp_keepdims_7d_weighted
        module procedure logsumexp_keepdims_real_7d_complex_b
        module procedure logsumexp_keepdims_complex_7d_weighted
        module procedure logsumexp_keepdims_complex_7d_real_b
        module procedure logsumexp_keepdims_8d_weighted
        module procedure logsumexp_keepdims_real_8d_complex_b
        module procedure logsumexp_keepdims_complex_8d_weighted
        module procedure logsumexp_keepdims_complex_8d_real_b
        module procedure logsumexp_keepdims_9d_weighted
        module procedure logsumexp_keepdims_real_9d_complex_b
        module procedure logsumexp_keepdims_complex_9d_weighted
        module procedure logsumexp_keepdims_complex_9d_real_b
        module procedure logsumexp_keepdims_10d_weighted
        module procedure logsumexp_keepdims_real_10d_complex_b
        module procedure logsumexp_keepdims_complex_10d_weighted
        module procedure logsumexp_keepdims_complex_10d_real_b
    end interface logsumexp_keepdims

    interface logsumexp_signed_keepdims
        module procedure logsumexp_signed_keepdims_6d_weighted
        module procedure logsumexp_signed_keepdims_real_6d_complex_b
        module procedure logsumexp_signed_keepdims_complex_6d_weighted
        module procedure logsumexp_signed_keepdims_complex_6d_real_b
        module procedure logsumexp_signed_keepdims_7d_weighted
        module procedure logsumexp_signed_keepdims_real_7d_complex_b
        module procedure logsumexp_signed_keepdims_complex_7d_weighted
        module procedure logsumexp_signed_keepdims_complex_7d_real_b
        module procedure logsumexp_signed_keepdims_8d_weighted
        module procedure logsumexp_signed_keepdims_real_8d_complex_b
        module procedure logsumexp_signed_keepdims_complex_8d_weighted
        module procedure logsumexp_signed_keepdims_complex_8d_real_b
        module procedure logsumexp_signed_keepdims_9d_weighted
        module procedure logsumexp_signed_keepdims_real_9d_complex_b
        module procedure logsumexp_signed_keepdims_complex_9d_weighted
        module procedure logsumexp_signed_keepdims_complex_9d_real_b
        module procedure logsumexp_signed_keepdims_10d_weighted
        module procedure logsumexp_signed_keepdims_real_10d_complex_b
        module procedure logsumexp_signed_keepdims_complex_10d_weighted
        module procedure logsumexp_signed_keepdims_complex_10d_real_b
    end interface logsumexp_signed_keepdims

contains

    ! Weighted rank-6 wrappers.
    pure subroutine logsumexp_signed_keepdims_6d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(6)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_6d_weighted

    pure function logsumexp_keepdims_6d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :)

        call logsumexp_signed_keepdims_6d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_6d_weighted

    pure subroutine logsumexp_signed_6d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_6d_weighted(x, axes_6, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_6d_all_weighted

    pure function logsumexp_6d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_6d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_6d_all_weighted

    pure subroutine logsumexp_signed_6d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :)
        integer :: outshape(5)

        if (dim < 1 .or. dim > 6) then
            allocate(logabs(0, 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_6d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_6d_dim_weighted

    pure function logsumexp_6d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :)

        call logsumexp_signed_6d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_6d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_6d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(6)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_6d_complex_b

    pure function logsumexp_keepdims_real_6d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_6d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_6d_complex_b

    pure subroutine logsumexp_signed_real_6d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_real_6d_complex_b(x, axes_6, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_6d_all_complex_b

    pure function logsumexp_real_6d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_6d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_6d_all_complex_b

    pure subroutine logsumexp_signed_real_6d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :)
        integer :: outshape(5)

        if (dim < 1 .or. dim > 6) then
            allocate(logabs(0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_6d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_6d_dim_complex_b

    pure function logsumexp_real_6d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :)

        call logsumexp_signed_real_6d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_6d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_6d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(6)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_6d_weighted

    pure function logsumexp_keepdims_complex_6d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_6d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_6d_weighted

    pure subroutine logsumexp_signed_complex_6d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_complex_6d_weighted(x, axes_6, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_6d_all_weighted

    pure function logsumexp_complex_6d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_6d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_6d_all_weighted

    pure subroutine logsumexp_signed_complex_6d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :)
        integer :: outshape(5)

        if (dim < 1 .or. dim > 6) then
            allocate(logabs(0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_6d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_6d_dim_weighted

    pure function logsumexp_complex_6d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :)

        call logsumexp_signed_complex_6d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_6d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_6d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(6)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_6d_real_b

    pure function logsumexp_keepdims_complex_6d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_6d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_6d_real_b

    pure subroutine logsumexp_signed_complex_6d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_complex_6d_real_b(x, axes_6, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_6d_all_real_b

    pure function logsumexp_complex_6d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_6d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_6d_all_real_b

    pure subroutine logsumexp_signed_complex_6d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :)
        integer :: outshape(5)

        if (dim < 1 .or. dim > 6) then
            allocate(logabs(0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_6d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_6d_dim_real_b

    pure function logsumexp_complex_6d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :)

        call logsumexp_signed_complex_6d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_6d_dim_real_b

    ! Weighted rank-7 wrappers.
    pure subroutine logsumexp_signed_keepdims_7d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(7)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_7d_weighted

    pure function logsumexp_keepdims_7d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_7d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_7d_weighted

    pure subroutine logsumexp_signed_7d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_7d_weighted(x, axes_7, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_7d_all_weighted

    pure function logsumexp_7d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_7d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_7d_all_weighted

    pure subroutine logsumexp_signed_7d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
            allocate(logabs(0, 0, 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_7d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_7d_dim_weighted

    pure function logsumexp_7d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :)

        call logsumexp_signed_7d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_7d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_7d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(7)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_7d_complex_b

    pure function logsumexp_keepdims_real_7d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :)

        call logsumexp_signed_keepdims_real_7d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_7d_complex_b

    pure subroutine logsumexp_signed_real_7d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_real_7d_complex_b(x, axes_7, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_7d_all_complex_b

    pure function logsumexp_real_7d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_7d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_7d_all_complex_b

    pure subroutine logsumexp_signed_real_7d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
            allocate(logabs(0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_7d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_7d_dim_complex_b

    pure function logsumexp_real_7d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        call logsumexp_signed_real_7d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_7d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_7d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(7)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_7d_weighted

    pure function logsumexp_keepdims_complex_7d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :)

        call logsumexp_signed_keepdims_complex_7d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_7d_weighted

    pure subroutine logsumexp_signed_complex_7d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_complex_7d_weighted(x, axes_7, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_7d_all_weighted

    pure function logsumexp_complex_7d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_7d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_7d_all_weighted

    pure subroutine logsumexp_signed_complex_7d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
            allocate(logabs(0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_7d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_7d_dim_weighted

    pure function logsumexp_complex_7d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        call logsumexp_signed_complex_7d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_7d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_7d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(7)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_7d_real_b

    pure function logsumexp_keepdims_complex_7d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :)

        call logsumexp_signed_keepdims_complex_7d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_7d_real_b

    pure subroutine logsumexp_signed_complex_7d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_complex_7d_real_b(x, axes_7, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_7d_all_real_b

    pure function logsumexp_complex_7d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_7d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_7d_all_real_b

    pure subroutine logsumexp_signed_complex_7d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :)
        integer :: outshape(6)

        if (dim < 1 .or. dim > 7) then
            allocate(logabs(0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_7d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_7d_dim_real_b

    pure function logsumexp_complex_7d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 7
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :)

        call logsumexp_signed_complex_7d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_7d_dim_real_b

    ! Weighted rank-8 wrappers.
    pure subroutine logsumexp_signed_keepdims_8d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(8)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_8d_weighted

    pure function logsumexp_keepdims_8d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_8d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_8d_weighted

    pure subroutine logsumexp_signed_8d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_8d_weighted(x, axes_8, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_8d_all_weighted

    pure function logsumexp_8d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_8d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_8d_all_weighted

    pure subroutine logsumexp_signed_8d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_8d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_8d_dim_weighted

    pure function logsumexp_8d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :)

        call logsumexp_signed_8d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_8d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_8d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(8)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_8d_complex_b

    pure function logsumexp_keepdims_real_8d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_real_8d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_8d_complex_b

    pure subroutine logsumexp_signed_real_8d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_real_8d_complex_b(x, axes_8, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_8d_all_complex_b

    pure function logsumexp_real_8d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_8d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_8d_all_complex_b

    pure subroutine logsumexp_signed_real_8d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_8d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_8d_dim_complex_b

    pure function logsumexp_real_8d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :)

        call logsumexp_signed_real_8d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_8d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_8d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(8)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_8d_weighted

    pure function logsumexp_keepdims_complex_8d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_complex_8d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_8d_weighted

    pure subroutine logsumexp_signed_complex_8d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_complex_8d_weighted(x, axes_8, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_8d_all_weighted

    pure function logsumexp_complex_8d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_8d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_8d_all_weighted

    pure subroutine logsumexp_signed_complex_8d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_8d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_8d_dim_weighted

    pure function logsumexp_complex_8d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :)

        call logsumexp_signed_complex_8d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_8d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_8d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(8)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_8d_real_b

    pure function logsumexp_keepdims_complex_8d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :)

        call logsumexp_signed_keepdims_complex_8d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_8d_real_b

    pure subroutine logsumexp_signed_complex_8d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_complex_8d_real_b(x, axes_8, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_8d_all_real_b

    pure function logsumexp_complex_8d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_8d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_8d_all_real_b

    pure subroutine logsumexp_signed_complex_8d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :)
        integer :: outshape(7)

        if (dim < 1 .or. dim > 8) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_8d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_8d_dim_real_b

    pure function logsumexp_complex_8d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 8
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :)

        call logsumexp_signed_complex_8d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_8d_dim_real_b

    ! Weighted rank-9 wrappers.
    pure subroutine logsumexp_signed_keepdims_9d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(9)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_9d_weighted

    pure function logsumexp_keepdims_9d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_9d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_9d_weighted

    pure subroutine logsumexp_signed_9d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_9d_weighted(x, axes_9, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_9d_all_weighted

    pure function logsumexp_9d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_9d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_9d_all_weighted

    pure subroutine logsumexp_signed_9d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_9d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_9d_dim_weighted

    pure function logsumexp_9d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_9d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_9d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_9d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(9)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_9d_complex_b

    pure function logsumexp_keepdims_real_9d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_real_9d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_9d_complex_b

    pure subroutine logsumexp_signed_real_9d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_real_9d_complex_b(x, axes_9, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_9d_all_complex_b

    pure function logsumexp_real_9d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_9d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_9d_all_complex_b

    pure subroutine logsumexp_signed_real_9d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_9d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_9d_dim_complex_b

    pure function logsumexp_real_9d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :)

        call logsumexp_signed_real_9d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_9d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_9d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(9)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_9d_weighted

    pure function logsumexp_keepdims_complex_9d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_complex_9d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_9d_weighted

    pure subroutine logsumexp_signed_complex_9d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_9d_weighted(x, axes_9, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_9d_all_weighted

    pure function logsumexp_complex_9d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_9d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_9d_all_weighted

    pure subroutine logsumexp_signed_complex_9d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_9d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_9d_dim_weighted

    pure function logsumexp_complex_9d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :)

        call logsumexp_signed_complex_9d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_9d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_9d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(9)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_9d_real_b

    pure function logsumexp_keepdims_complex_9d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_keepdims_complex_9d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_9d_real_b

    pure subroutine logsumexp_signed_complex_9d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_9d_real_b(x, axes_9, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_9d_all_real_b

    pure function logsumexp_complex_9d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_9d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_9d_all_real_b

    pure subroutine logsumexp_signed_complex_9d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :)
        integer :: outshape(8)

        if (dim < 1 .or. dim > 9) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_9d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_9d_dim_real_b

    pure function logsumexp_complex_9d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 9
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :)

        call logsumexp_signed_complex_9d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_9d_dim_real_b

    ! Weighted rank-10 wrappers.
    pure subroutine logsumexp_signed_keepdims_10d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :, :) !! signs of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        real(dp), allocatable :: flat_sign(:)
        logical :: valid
        integer :: outshape(10)

        call flatten_real_weights(b, bflat, bshape)
        call weighted_real_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_sign, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        sign_value = reshape(flat_sign, shape(sign_value))
    end subroutine logsumexp_signed_keepdims_10d_weighted

    pure function logsumexp_keepdims_10d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :, :)

        call logsumexp_signed_keepdims_10d_weighted(x, dims, logabs, sign_value, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_10d_weighted

    pure subroutine logsumexp_signed_10d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_10d_weighted(x, axes_10, kept_logabs, kept_sign_value, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        logabs = sum(kept_logabs)
        sign_value = sum(kept_sign_value)
    end subroutine logsumexp_signed_10d_all_weighted

    pure function logsumexp_10d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_10d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_10d_all_weighted

    pure subroutine logsumexp_signed_10d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, &
            & :, :, :, :, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        real(dp), allocatable :: kept_sign_value(:, :, :, :, &
            & :, :, :, :, :, :)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(sign_value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_10d_weighted(x, [dim], kept_logabs, kept_sign_value, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(sign_value(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign_value, shape(sign_value))
    end subroutine logsumexp_signed_10d_dim_weighted

    pure function logsumexp_10d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        real(dp), allocatable :: values(:, :, :, :, :, :, :, &
            & :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        real(dp), allocatable :: sign_value(:, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_10d_dim_weighted(x, dim, logabs, sign_value, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        if (size(values) == 0) return
        values = logabs
        where (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_10d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_real_10d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_x(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(10)

        call flatten_complex_weights(b, bflat, bshape)
        allocate(complex_x(size(x)))
        if (size(x) > 0) complex_x = cmplx(reshape(x, [size(x)]), 0.0_dp, kind=dp)
        call weighted_complex_prepared(complex_x, shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_real_10d_complex_b

    pure function logsumexp_keepdims_real_10d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_real_10d_complex_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_real_10d_complex_b

    pure subroutine logsumexp_signed_real_10d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_real_10d_complex_b(x, axes_10, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_real_10d_all_complex_b

    pure function logsumexp_real_10d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_real_10d_all_complex_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_10d_all_complex_b

    pure subroutine logsumexp_signed_real_10d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_real_10d_complex_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_real_10d_dim_complex_b

    pure function logsumexp_real_10d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, &
            & :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_real_10d_dim_complex_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_real_10d_dim_complex_b

    pure subroutine logsumexp_signed_keepdims_complex_10d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        complex(dp), allocatable :: bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(10)

        call flatten_complex_weights(b, bflat, bshape)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, bflat, bshape, &
            flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_10d_weighted

    pure function logsumexp_keepdims_complex_10d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_complex_10d_weighted(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_10d_weighted

    pure subroutine logsumexp_signed_complex_10d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_10d_weighted(x, axes_10, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_10d_all_weighted

    pure function logsumexp_complex_10d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_10d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_10d_all_weighted

    pure subroutine logsumexp_signed_complex_10d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_10d_weighted(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_10d_dim_weighted

    pure function logsumexp_complex_10d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        complex(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_complex_10d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_10d_dim_weighted

    pure subroutine logsumexp_signed_keepdims_complex_10d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :, :) !! log-magnitudes of weighted exponential sums
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :, :) !! unit phases of weighted exponential sums
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        integer, allocatable :: bshape(:)
        real(dp), allocatable :: bflat(:)
        complex(dp), allocatable :: complex_bflat(:)
        real(dp), allocatable :: flat_logabs(:)
        complex(dp), allocatable :: flat_phase(:)
        logical :: valid
        integer :: outshape(10)

        call flatten_real_weights(b, bflat, bshape)
        allocate(complex_bflat(size(bflat)))
        if (size(bflat) > 0) complex_bflat = cmplx(bflat, 0.0_dp, kind=dp)
        call weighted_complex_prepared(reshape(x, [size(x)]), shape(x), dims, &
            complex_bflat, bshape, flat_logabs, flat_phase, outshape, valid)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9), outshape(10)))
        if (.not. valid .or. size(logabs) == 0) return
        logabs = reshape(flat_logabs, shape(logabs))
        phase = reshape(flat_phase, shape(phase))
    end subroutine logsumexp_signed_keepdims_complex_10d_real_b

    pure function logsumexp_keepdims_complex_10d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :, :) !! weighted log-sum-exp with reduced dimensions retained

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :, :)

        call logsumexp_signed_keepdims_complex_10d_real_b(x, dims, logabs, phase, b)
        if (size(logabs) == 0) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9), size(logabs, 10)))
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_keepdims_complex_10d_real_b

    pure subroutine logsumexp_signed_complex_10d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_10d_real_b(x, axes_10, kept_logabs, kept_phase, b)
        if (size(kept_logabs) == 0) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        logabs = sum(kept_logabs)
        phase = sum(kept_phase)
    end subroutine logsumexp_signed_complex_10d_all_real_b

    pure function logsumexp_complex_10d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_10d_all_real_b(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_10d_all_real_b

    pure subroutine logsumexp_signed_complex_10d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), allocatable, intent(out) :: logabs(:, :, &
            & :, :, :, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, &
            & :, :, :, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, &
            & :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, &
            & :, :, :, :, :, :)
        integer :: outshape(9)

        if (dim < 1 .or. dim > 10) then
            allocate(logabs(0, 0, 0, 0, 0, 0, 0, 0, 0))
            allocate(phase(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_10d_real_b(x, [dim], kept_logabs, kept_phase, b)
        call drop_dimension_shape(shape(kept_logabs), dim, outshape)
        allocate(logabs(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        allocate(phase(outshape(1), outshape(2), &
            & outshape(3), outshape(4), outshape(5), &
            & outshape(6), outshape(7), outshape(8), &
            & outshape(9)))
        if (size(logabs) == 0) return
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_10d_dim_real_b

    pure function logsumexp_complex_10d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :, :, :, &
            & :, :) !! logarithmic terms to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 10
        real(dp), intent(in) :: b(..) !! scalar or trailing weights broadcast jointly with x
        complex(dp), allocatable :: values(:, :, :, :, :, :, &
            & :, :, :) !! weighted log-sum-exp after dropping the reduced dimension

        real(dp), allocatable :: logabs(:, :, :, :, :, :, :, &
            & :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :, :, &
            & :, :, :)

        call logsumexp_signed_complex_10d_dim_real_b(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), &
            & size(logabs, 3), size(logabs, 4), &
            & size(logabs, 5), size(logabs, 6), &
            & size(logabs, 7), size(logabs, 8), &
            & size(logabs, 9)))
        if (size(values) == 0) return
        values = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_10d_dim_real_b

end module scifort_special_reductions_high_rank_weighted_6_10
