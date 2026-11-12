! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions_low
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: log_softmax
    public :: logsumexp
    public :: logsumexp_signed
    public :: logsumexp_keepdims
    public :: logsumexp_signed_keepdims
    public :: softmax

    interface logsumexp
        module procedure logsumexp_1d
        module procedure logsumexp_2d_all
        module procedure logsumexp_2d_dim
        module procedure logsumexp_3d_all
        module procedure logsumexp_3d_dim
        module procedure logsumexp_4d_all
        module procedure logsumexp_4d_dim
        module procedure logsumexp_5d_all
        module procedure logsumexp_5d_dim
        module procedure logsumexp_6d_all
        module procedure logsumexp_6d_dim
        module procedure logsumexp_5d_all_weighted
        module procedure logsumexp_5d_dim_weighted
        module procedure logsumexp_real_5d_all_complex_b
        module procedure logsumexp_real_5d_dim_complex_b
        module procedure logsumexp_4d_all_weighted
        module procedure logsumexp_4d_dim_weighted
        module procedure logsumexp_real_4d_all_complex_b
        module procedure logsumexp_real_4d_dim_complex_b
        module procedure logsumexp_real_1d_complex_b
        module procedure logsumexp_real_2d_all_complex_b
        module procedure logsumexp_real_2d_dim_complex_b
        module procedure logsumexp_real_3d_all_complex_b
        module procedure logsumexp_real_3d_dim_complex_b
        module procedure logsumexp_complex_1d
        module procedure logsumexp_complex_1d_real_b
        module procedure logsumexp_complex_2d_all
        module procedure logsumexp_complex_2d_all_real_b
        module procedure logsumexp_complex_2d_dim
        module procedure logsumexp_complex_2d_dim_real_b
        module procedure logsumexp_complex_3d_all
        module procedure logsumexp_complex_3d_all_real_b
        module procedure logsumexp_complex_3d_dim
        module procedure logsumexp_complex_3d_dim_real_b
        module procedure logsumexp_complex_4d_all
        module procedure logsumexp_complex_4d_dim
        module procedure logsumexp_complex_5d_all
        module procedure logsumexp_complex_5d_dim
        module procedure logsumexp_complex_6d_all
        module procedure logsumexp_complex_6d_dim
        module procedure logsumexp_complex_5d_all_weighted
        module procedure logsumexp_complex_5d_dim_weighted
        module procedure logsumexp_complex_5d_all_real_b
        module procedure logsumexp_complex_5d_dim_real_b
        module procedure logsumexp_complex_4d_all_weighted
        module procedure logsumexp_complex_4d_dim_weighted
        module procedure logsumexp_complex_4d_all_real_b
        module procedure logsumexp_complex_4d_dim_real_b
    end interface logsumexp

    interface logsumexp_signed
        module procedure logsumexp_signed_1d
        module procedure logsumexp_signed_2d_all
        module procedure logsumexp_signed_2d_dim
        module procedure logsumexp_signed_3d_all
        module procedure logsumexp_signed_3d_dim
        module procedure logsumexp_signed_4d_all
        module procedure logsumexp_signed_4d_dim
        module procedure logsumexp_signed_5d_all
        module procedure logsumexp_signed_5d_dim
        module procedure logsumexp_signed_6d_all
        module procedure logsumexp_signed_6d_dim
        module procedure logsumexp_signed_5d_all_weighted
        module procedure logsumexp_signed_5d_dim_weighted
        module procedure logsumexp_signed_real_5d_all_complex_b
        module procedure logsumexp_signed_real_5d_dim_complex_b
        module procedure logsumexp_signed_4d_all_weighted
        module procedure logsumexp_signed_4d_dim_weighted
        module procedure logsumexp_signed_real_4d_all_complex_b
        module procedure logsumexp_signed_real_4d_dim_complex_b
        module procedure logsumexp_signed_real_1d_complex_b
        module procedure logsumexp_signed_real_2d_all_complex_b
        module procedure logsumexp_signed_real_2d_dim_complex_b
        module procedure logsumexp_signed_real_3d_all_complex_b
        module procedure logsumexp_signed_real_3d_dim_complex_b
        module procedure logsumexp_signed_complex_1d
        module procedure logsumexp_signed_complex_1d_real_b
        module procedure logsumexp_signed_complex_2d_all
        module procedure logsumexp_signed_complex_2d_all_real_b
        module procedure logsumexp_signed_complex_2d_dim
        module procedure logsumexp_signed_complex_2d_dim_real_b
        module procedure logsumexp_signed_complex_3d_all
        module procedure logsumexp_signed_complex_3d_all_real_b
        module procedure logsumexp_signed_complex_3d_dim
        module procedure logsumexp_signed_complex_3d_dim_real_b
        module procedure logsumexp_signed_complex_4d_all
        module procedure logsumexp_signed_complex_4d_dim
        module procedure logsumexp_signed_complex_5d_all
        module procedure logsumexp_signed_complex_5d_dim
        module procedure logsumexp_signed_complex_6d_all
        module procedure logsumexp_signed_complex_6d_dim
        module procedure logsumexp_signed_complex_5d_all_weighted
        module procedure logsumexp_signed_complex_5d_dim_weighted
        module procedure logsumexp_signed_complex_5d_all_real_b
        module procedure logsumexp_signed_complex_5d_dim_real_b
        module procedure logsumexp_signed_complex_4d_all_weighted
        module procedure logsumexp_signed_complex_4d_dim_weighted
        module procedure logsumexp_signed_complex_4d_all_real_b
        module procedure logsumexp_signed_complex_4d_dim_real_b
    end interface logsumexp_signed

    interface logsumexp_keepdims
        module procedure logsumexp_keepdims_2d
        module procedure logsumexp_keepdims_3d
        module procedure logsumexp_keepdims_4d
        module procedure logsumexp_keepdims_5d
        module procedure logsumexp_keepdims_6d
        module procedure logsumexp_keepdims_5d_weighted
        module procedure logsumexp_keepdims_real_5d_complex_b
        module procedure logsumexp_keepdims_4d_weighted
        module procedure logsumexp_keepdims_real_4d_complex_b
        module procedure logsumexp_keepdims_real_2d_complex_b
        module procedure logsumexp_keepdims_real_3d_complex_b
        module procedure logsumexp_keepdims_complex_2d
        module procedure logsumexp_keepdims_complex_2d_real_b
        module procedure logsumexp_keepdims_complex_3d
        module procedure logsumexp_keepdims_complex_3d_real_b
        module procedure logsumexp_keepdims_complex_4d
        module procedure logsumexp_keepdims_complex_5d
        module procedure logsumexp_keepdims_complex_6d
        module procedure logsumexp_keepdims_complex_5d_weighted
        module procedure logsumexp_keepdims_complex_5d_real_b
        module procedure logsumexp_keepdims_complex_4d_weighted
        module procedure logsumexp_keepdims_complex_4d_real_b
    end interface logsumexp_keepdims

    interface logsumexp_signed_keepdims
        module procedure logsumexp_signed_keepdims_2d
        module procedure logsumexp_signed_keepdims_3d
        module procedure logsumexp_signed_keepdims_4d
        module procedure logsumexp_signed_keepdims_5d
        module procedure logsumexp_signed_keepdims_6d
        module procedure logsumexp_signed_keepdims_5d_weighted
        module procedure logsumexp_signed_keepdims_real_5d_complex_b
        module procedure logsumexp_signed_keepdims_4d_weighted
        module procedure logsumexp_signed_keepdims_real_4d_complex_b
        module procedure logsumexp_signed_keepdims_real_2d_complex_b
        module procedure logsumexp_signed_keepdims_real_3d_complex_b
        module procedure logsumexp_signed_keepdims_complex_2d
        module procedure logsumexp_signed_keepdims_complex_2d_real_b
        module procedure logsumexp_signed_keepdims_complex_3d
        module procedure logsumexp_signed_keepdims_complex_3d_real_b
        module procedure logsumexp_signed_keepdims_complex_4d
        module procedure logsumexp_signed_keepdims_complex_5d
        module procedure logsumexp_signed_keepdims_complex_6d
        module procedure logsumexp_signed_keepdims_complex_5d_weighted
        module procedure logsumexp_signed_keepdims_complex_5d_real_b
        module procedure logsumexp_signed_keepdims_complex_4d_weighted
        module procedure logsumexp_signed_keepdims_complex_4d_real_b
    end interface logsumexp_signed_keepdims

    interface softmax
        module procedure softmax_1d
        module procedure softmax_2d
        module procedure softmax_2d_dims
        module procedure softmax_3d
        module procedure softmax_3d_dims
        module procedure softmax_4d
        module procedure softmax_4d_dims
        module procedure softmax_5d
        module procedure softmax_5d_dims
        module procedure softmax_6d
        module procedure softmax_6d_dims
        module procedure softmax_complex_1d
        module procedure softmax_complex_2d
        module procedure softmax_complex_2d_dims
        module procedure softmax_complex_3d
        module procedure softmax_complex_3d_dims
        module procedure softmax_complex_4d
        module procedure softmax_complex_4d_dims
        module procedure softmax_complex_5d
        module procedure softmax_complex_5d_dims
        module procedure softmax_complex_6d
        module procedure softmax_complex_6d_dims
    end interface softmax

    interface log_softmax
        module procedure log_softmax_1d
        module procedure log_softmax_2d
        module procedure log_softmax_2d_dims
        module procedure log_softmax_3d
        module procedure log_softmax_3d_dims
        module procedure log_softmax_4d
        module procedure log_softmax_4d_dims
        module procedure log_softmax_5d
        module procedure log_softmax_5d_dims
        module procedure log_softmax_6d
        module procedure log_softmax_6d_dims
        module procedure log_softmax_complex_1d
        module procedure log_softmax_complex_2d
        module procedure log_softmax_complex_2d_dims
        module procedure log_softmax_complex_3d
        module procedure log_softmax_complex_3d_dims
        module procedure log_softmax_complex_4d
        module procedure log_softmax_complex_4d_dims
        module procedure log_softmax_complex_5d
        module procedure log_softmax_complex_5d_dims
        module procedure log_softmax_complex_6d
        module procedure log_softmax_complex_6d_dims
    end interface log_softmax

contains

    pure function logsumexp_1d(x, b) result(value)
        real(dp), intent(in) :: x(:) !! one-dimensional values to reduce
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)
        real(dp) :: value

        real(dp), allocatable :: terms(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: logabs
        real(dp) :: sign_value
        logical :: ok

        if (present(b)) then
            call broadcast_real_terms_1d(x, b, terms, weights, ok)
            if (.not. ok) then
                value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_core(terms, logabs, sign_value, weights)
        else
            call logsumexp_core(x, logabs, sign_value)
        end if

        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = logabs
        end if
    end function logsumexp_1d

    pure function logsumexp_2d_all(x, b) result(value)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values reduced over both dimensions
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)
        real(dp) :: value

        real(dp), allocatable :: terms(:, :)
        real(dp), allocatable :: weights(:, :)
        real(dp) :: logabs
        real(dp) :: sign_value
        logical :: ok

        if (present(b)) then
            call broadcast_real_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
        end if

        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = logabs
        end if
    end function logsumexp_2d_all

    pure function logsumexp_2d_dim(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)
        real(dp), allocatable :: values(:)

        integer :: i
        logical :: ok
        real(dp), allocatable :: terms(:, :)
        real(dp), allocatable :: weights(:, :)
        real(dp) :: logabs
        real(dp) :: sign_value

        if (dim /= 1 .and. dim /= 2) then
            allocate(values(0))
            return
        end if

        if (present(b)) then
            call broadcast_real_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                if (dim == 1) then
                    allocate(values(size(x, 2)))
                else
                    allocate(values(size(x, 1)))
                end if
                values = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2)))
            terms = x
        end if

        if (dim == 1) then
            allocate(values(size(terms, 2)))
            do i = 1, size(terms, 2)
                if (present(b)) then
                    call logsumexp_core(terms(:, i), logabs, sign_value, weights(:, i))
                else
                    call logsumexp_core(terms(:, i), logabs, sign_value)
                end if
                if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                    values(i) = quiet_nan(0.0_dp)
                else
                    values(i) = logabs
                end if
            end do
        else
            allocate(values(size(terms, 1)))
            do i = 1, size(terms, 1)
                if (present(b)) then
                    call logsumexp_core(terms(i, :), logabs, sign_value, weights(i, :))
                else
                    call logsumexp_core(terms(i, :), logabs, sign_value)
                end if
                if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                    values(i) = quiet_nan(0.0_dp)
                else
                    values(i) = logabs
                end if
            end do
        end if
    end function logsumexp_2d_dim

    pure subroutine logsumexp_signed_1d(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:) !! one-dimensional values to reduce
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of the weighted exponential sum: -1, 0, or +1
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)

        logical :: ok
        real(dp), allocatable :: terms(:)
        real(dp), allocatable :: weights(:)

        if (present(b)) then
            call broadcast_real_terms_1d(x, b, terms, weights, ok)
            if (.not. ok) then
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_core(terms, logabs, sign_value, weights)
        else
            call logsumexp_core(x, logabs, sign_value)
        end if
    end subroutine logsumexp_signed_1d

    pure subroutine logsumexp_signed_2d_all(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values reduced over both dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of the weighted exponential sum: -1, 0, or +1
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)

        logical :: ok
        real(dp), allocatable :: terms(:, :)
        real(dp), allocatable :: weights(:, :)

        if (present(b)) then
            call broadcast_real_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
        end if
    end subroutine logsumexp_signed_2d_all

    pure subroutine logsumexp_signed_2d_dim(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        real(dp), allocatable, intent(out) :: logabs(:) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:) !! signs along the reduced dimension
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)

        integer :: i
        logical :: ok
        real(dp), allocatable :: terms(:, :)
        real(dp), allocatable :: weights(:, :)

        if (dim /= 1 .and. dim /= 2) then
            allocate(logabs(0), sign_value(0))
            return
        end if

        if (present(b)) then
            call broadcast_real_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                if (dim == 1) then
                    allocate(logabs(size(x, 2)), sign_value(size(x, 2)))
                else
                    allocate(logabs(size(x, 1)), sign_value(size(x, 1)))
                end if
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2)))
            terms = x
        end if

        if (dim == 1) then
            allocate(logabs(size(terms, 2)), sign_value(size(terms, 2)))
            do i = 1, size(terms, 2)
                if (present(b)) then
                    call logsumexp_core(terms(:, i), logabs(i), sign_value(i), weights(:, i))
                else
                    call logsumexp_core(terms(:, i), logabs(i), sign_value(i))
                end if
            end do
        else
            allocate(logabs(size(terms, 1)), sign_value(size(terms, 1)))
            do i = 1, size(terms, 1)
                if (present(b)) then
                    call logsumexp_core(terms(i, :), logabs(i), sign_value(i), weights(i, :))
                else
                    call logsumexp_core(terms(i, :), logabs(i), sign_value(i))
                end if
            end do
        end if
    end subroutine logsumexp_signed_2d_dim

    pure function logsumexp_keepdims_2d(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)
        real(dp), allocatable :: values(:, :)

        integer :: i
        integer :: j
        integer :: mode
        logical :: ok
        real(dp), allocatable :: terms(:, :)
        real(dp), allocatable :: weights(:, :)
        real(dp) :: logabs
        real(dp) :: sign_value

        mode = reduction_dims_mode(dims)
        if (mode < 0) then
            allocate(values(0, 0))
            return
        end if

        if (present(b)) then
            call broadcast_real_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_real(values, size(x, 1), size(x, 2), mode)
                values = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2)))
            terms = x
        end if

        select case (mode)
        case (0)
            allocate(values(size(terms, 1), size(terms, 2)))
            if (.not. present(b)) then
                values = terms
            else
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_core([terms(i, j)], logabs, sign_value, [weights(i, j)])
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, j) = quiet_nan(0.0_dp)
                        else
                            values(i, j) = logabs
                        end if
                    end do
                end do
            end if
        case (1)
            allocate(values(1, size(terms, 2)))
            if (present(b)) then
                values(1, :) = logsumexp_2d_dim(terms, 1, weights)
            else
                values(1, :) = logsumexp_2d_dim(terms, 1)
            end if
        case (2)
            allocate(values(size(terms, 1), 1))
            if (present(b)) then
                values(:, 1) = logsumexp_2d_dim(terms, 2, weights)
            else
                values(:, 1) = logsumexp_2d_dim(terms, 2)
            end if
        case (3)
            allocate(values(1, 1))
            if (present(b)) then
                values(1, 1) = logsumexp_2d_all(terms, weights)
            else
                values(1, 1) = logsumexp_2d_all(x)
            end if
        end select
    end function logsumexp_keepdims_2d

    pure subroutine logsumexp_signed_keepdims_2d(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :) !! log-absolute sums with reduced dimensions retained
        real(dp), allocatable, intent(out) :: sign_value(:, :) !! signs with reduced dimensions retained
        real(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable weights multiplying exp(x)

        integer :: i
        integer :: j
        integer :: mode
        logical :: ok
        real(dp), allocatable :: reduced_logabs(:)
        real(dp), allocatable :: reduced_sign(:)
        real(dp), allocatable :: terms(:, :)
        real(dp), allocatable :: weights(:, :)

        mode = reduction_dims_mode(dims)
        if (mode < 0) then
            allocate(logabs(0, 0), sign_value(0, 0))
            return
        end if

        if (present(b)) then
            call broadcast_real_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_real(logabs, size(x, 1), size(x, 2), mode)
                call allocate_keepdims_real(sign_value, size(x, 1), size(x, 2), mode)
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2)))
            terms = x
        end if

        select case (mode)
        case (0)
            allocate(logabs(size(terms, 1), size(terms, 2)))
            allocate(sign_value(size(terms, 1), size(terms, 2)))
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_core([terms(i, j)], logabs(i, j), sign_value(i, j), &
                            [weights(i, j)])
                    else
                        call logsumexp_core([terms(i, j)], logabs(i, j), sign_value(i, j))
                    end if
                end do
            end do
        case (1)
            allocate(logabs(1, size(terms, 2)), sign_value(1, size(terms, 2)))
            if (present(b)) then
                call logsumexp_signed_2d_dim(terms, 1, reduced_logabs, reduced_sign, weights)
            else
                call logsumexp_signed_2d_dim(terms, 1, reduced_logabs, reduced_sign)
            end if
            logabs(1, :) = reduced_logabs
            sign_value(1, :) = reduced_sign
        case (2)
            allocate(logabs(size(terms, 1), 1), sign_value(size(terms, 1), 1))
            if (present(b)) then
                call logsumexp_signed_2d_dim(terms, 2, reduced_logabs, reduced_sign, weights)
            else
                call logsumexp_signed_2d_dim(terms, 2, reduced_logabs, reduced_sign)
            end if
            logabs(:, 1) = reduced_logabs
            sign_value(:, 1) = reduced_sign
        case (3)
            allocate(logabs(1, 1), sign_value(1, 1))
            if (present(b)) then
                call logsumexp_signed_2d_all(terms, logabs(1, 1), sign_value(1, 1), weights)
            else
                call logsumexp_signed_2d_all(terms, logabs(1, 1), sign_value(1, 1))
            end if
        end select
    end subroutine logsumexp_signed_keepdims_2d

    pure subroutine broadcast_real_terms_1d(x, b, terms, weights, ok)
        real(dp), intent(in) :: x(:) !! one-dimensional logarithmic terms to broadcast
        real(dp), intent(in) :: b(..) !! scalar or rank-1 weights broadcast jointly with x
        real(dp), allocatable, intent(out) :: terms(:) !! input terms broadcast to the joint shape
        real(dp), allocatable, intent(out) :: weights(:) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: i
        integer :: ib
        integer :: ix
        integer :: n

        ok = .true.
        select rank (b)
        rank (0)
            n = size(x)
            allocate(terms(n), weights(n))
            terms = x
            weights = b
        rank (1)
            if (size(x) /= size(b) .and. size(x) /= 1 .and. size(b) /= 1) then
                ok = .false.
                allocate(terms(0), weights(0))
                return
            end if
            n = size(x)
            if (n == 1) n = size(b)
            allocate(terms(n), weights(n))
            do i = 1, n
                ix = merge(1, i, size(x) == 1)
                ib = merge(1, i, size(b) == 1)
                terms(i) = x(ix)
                weights(i) = b(ib)
            end do
        rank default
            ok = .false.
            allocate(terms(0), weights(0))
        end select
    end subroutine broadcast_real_terms_1d

    pure subroutine broadcast_real_terms_2d(x, b, terms, weights, ok)
        real(dp), intent(in) :: x(:, :) !! two-dimensional logarithmic terms to broadcast
        real(dp), intent(in) :: b(..) !! scalar, trailing-vector, or rank-2 weights
        real(dp), allocatable, intent(out) :: terms(:, :) !! input terms broadcast to the joint shape
        real(dp), allocatable, intent(out) :: weights(:, :) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: n1
        integer :: n2

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2)), weights(size(x, 1), size(x, 2)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = size(b)
        rank (2)
            b1 = size(b, 1)
            b2 = size(b, 2)
        rank default
            ok = .false.
            allocate(terms(0, 0), weights(0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1)) then
            ok = .false.
            allocate(terms(0, 0), weights(0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        allocate(terms(n1, n2), weights(n1, n2))
        select rank (b)
        rank (1)
            do j = 1, n2
                jx = merge(1, j, size(x, 2) == 1)
                jb = merge(1, j, size(b) == 1)
                do i = 1, n1
                    ix = merge(1, i, size(x, 1) == 1)
                    terms(i, j) = x(ix, jx)
                    weights(i, j) = b(jb)
                end do
            end do
        rank (2)
            do j = 1, n2
                jx = merge(1, j, size(x, 2) == 1)
                jb = merge(1, j, size(b, 2) == 1)
                do i = 1, n1
                    ix = merge(1, i, size(x, 1) == 1)
                    ib = merge(1, i, size(b, 1) == 1)
                    terms(i, j) = x(ix, jx)
                    weights(i, j) = b(ib, jb)
                end do
            end do
        end select
    end subroutine broadcast_real_terms_2d

    pure function reduction_dims_mode(dims) result(mode)
        integer, intent(in) :: dims(:) !! Fortran dimensions to reduce
        integer :: mode

        if (size(dims) == 0) then
            mode = 0
        else if (size(dims) == 1) then
            if (dims(1) == 1) then
                mode = 1
            else if (dims(1) == 2) then
                mode = 2
            else
                mode = -1
            end if
        else if (size(dims) == 2) then
            if ((dims(1) == 1 .and. dims(2) == 2) .or. &
                    (dims(1) == 2 .and. dims(2) == 1)) then
                mode = 3
            else
                mode = -1
            end if
        else
            mode = -1
        end if
    end function reduction_dims_mode

    pure subroutine allocate_keepdims_real(values, n1, n2, mode)
        real(dp), allocatable, intent(out) :: values(:, :) !! allocated keepdims result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: mode !! reduction mode from reduction_dims_mode

        if (mode < 0) then
            allocate(values(0, 0))
        else if (mode == 0) then
            allocate(values(n1, n2))
        else if (mode == 1) then
            allocate(values(1, n2))
        else if (mode == 2) then
            allocate(values(n1, 1))
        else
            allocate(values(1, 1))
        end if
    end subroutine allocate_keepdims_real

    pure subroutine logsumexp_core(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:) !! one-dimensional logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of the weighted exponential sum
        real(dp), intent(in), optional :: b(:) !! optional weights multiplying exp(x)

        integer :: i
        logical :: active
        logical :: has_negative_infinity_term
        logical :: has_positive_infinity_term
        real(dp) :: negative_compensation
        real(dp) :: negative_sum
        real(dp) :: positive_compensation
        real(dp) :: positive_sum
        real(dp) :: maximum_term
        real(dp) :: scaled_term
        real(dp) :: term_log
        real(dp) :: weight

        if (size(x) == 0) then
            logabs = negative_infinity(0.0_dp)
            sign_value = 0.0_dp
            return
        end if

        maximum_term = negative_infinity(0.0_dp)
        do i = 1, size(x)
            if (present(b)) then
                weight = b(i)
                active = weight /= 0.0_dp
            else
                weight = 1.0_dp
                active = .true.
            end if
            if (.not. active) cycle
            if (ieee_is_nan(x(i)) .or. ieee_is_nan(weight)) then
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
            if (present(b)) then
                term_log = x(i) + log(abs(weight))
            else
                term_log = x(i)
            end if
            if (ieee_is_nan(term_log)) then
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
            maximum_term = max(maximum_term, term_log)
        end do

        if (maximum_term == negative_infinity(0.0_dp)) then
            logabs = maximum_term
            sign_value = 0.0_dp
            return
        end if

        if (.not. ieee_is_finite(maximum_term)) then
            has_positive_infinity_term = .false.
            has_negative_infinity_term = .false.
            do i = 1, size(x)
                if (present(b)) then
                    weight = b(i)
                    active = weight /= 0.0_dp
                else
                    weight = 1.0_dp
                    active = .true.
                end if
                if (.not. active) cycle
                if (present(b)) then
                    term_log = x(i) + log(abs(weight))
                else
                    term_log = x(i)
                end if
                if (term_log == positive_infinity(0.0_dp)) then
                    if (weight > 0.0_dp) then
                        has_positive_infinity_term = .true.
                    else
                        has_negative_infinity_term = .true.
                    end if
                end if
            end do
            if (has_positive_infinity_term .and. has_negative_infinity_term) then
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
            else if (has_positive_infinity_term) then
                logabs = positive_infinity(0.0_dp)
                sign_value = 1.0_dp
            else if (has_negative_infinity_term) then
                logabs = positive_infinity(0.0_dp)
                sign_value = -1.0_dp
            else
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
            end if
            return
        end if

        positive_sum = 0.0_dp
        negative_sum = 0.0_dp
        positive_compensation = 0.0_dp
        negative_compensation = 0.0_dp
        do i = 1, size(x)
            if (present(b)) then
                weight = b(i)
                active = weight /= 0.0_dp
            else
                weight = 1.0_dp
                active = .true.
            end if
            if (.not. active) cycle
            if (present(b)) then
                term_log = x(i) + log(abs(weight))
            else
                term_log = x(i)
            end if
            scaled_term = exp(term_log - maximum_term)
            if (weight > 0.0_dp) then
                call compensated_add(positive_sum, positive_compensation, scaled_term)
            else
                call compensated_add(negative_sum, negative_compensation, scaled_term)
            end if
        end do

        if (positive_sum > negative_sum) then
            sign_value = 1.0_dp
            logabs = maximum_term + log(positive_sum) + log1p_safe(-negative_sum / positive_sum)
        else if (negative_sum > positive_sum) then
            sign_value = -1.0_dp
            logabs = maximum_term + log(negative_sum) + log1p_safe(-positive_sum / negative_sum)
        else
            sign_value = 0.0_dp
            logabs = negative_infinity(0.0_dp)
        end if
    end subroutine logsumexp_core


    pure subroutine compensated_add(total, compensation, value)
        real(dp), intent(inout) :: total !! running compensated sum
        real(dp), intent(inout) :: compensation !! correction term carried between additions
        real(dp), intent(in) :: value !! nonnegative value added to the running sum

        real(dp) :: adjusted
        real(dp) :: updated

        adjusted = value - compensation
        updated = total + adjusted
        compensation = (updated - total) - adjusted
        total = updated
    end subroutine compensated_add

    pure function softmax_1d(x) result(values)
        real(dp), intent(in) :: x(:) !! one-dimensional values to normalize
        real(dp), allocatable :: values(:)

        real(dp) :: maximum_value
        real(dp) :: total

        allocate(values(size(x)))
        if (size(x) == 0) return
        if (any(ieee_is_nan(x))) then
            values = quiet_nan(0.0_dp)
            return
        end if

        maximum_value = maxval(x)
        if (.not. ieee_is_finite(maximum_value)) then
            values = quiet_nan(0.0_dp)
            return
        end if

        values = exp(x - maximum_value)
        total = sum(values)
        values = values / total
    end function softmax_1d

    pure function softmax_2d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        real(dp), allocatable :: values(:, :)

        integer :: i
        real(dp) :: maximum_value
        real(dp) :: total

        allocate(values(size(x, 1), size(x, 2)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            if (any(ieee_is_nan(x))) then
                values = quiet_nan(0.0_dp)
                return
            end if
            maximum_value = maxval(x)
            if (.not. ieee_is_finite(maximum_value)) then
                values = quiet_nan(0.0_dp)
                return
            end if
            values = exp(x - maximum_value)
            total = sum(values)
            values = values / total
        else if (dim == 1) then
            do i = 1, size(x, 2)
                values(:, i) = softmax_1d(x(:, i))
            end do
        else if (dim == 2) then
            do i = 1, size(x, 1)
                values(i, :) = softmax_1d(x(i, :))
            end do
        else
            values = quiet_nan(0.0_dp)
        end if
    end function softmax_2d

    pure function softmax_2d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp), allocatable :: values(:, :)

        integer :: mode

        mode = reduction_dims_mode(dims)
        if (mode == 0) then
            allocate(values(size(x, 1), size(x, 2)))
            where (ieee_is_finite(x))
                values = 1.0_dp
            elsewhere
                values = quiet_nan(0.0_dp)
            end where
        else if (mode == 1) then
            values = softmax_2d(x, 1)
        else if (mode == 2) then
            values = softmax_2d(x, 2)
        else if (mode == 3) then
            values = softmax_2d(x)
        else
            allocate(values(size(x, 1), size(x, 2)))
            values = quiet_nan(0.0_dp)
        end if
    end function softmax_2d_dims

    pure function log_softmax_1d(x) result(values)
        real(dp), intent(in) :: x(:) !! one-dimensional values to log-normalize
        real(dp), allocatable :: values(:)

        integer :: i
        real(dp) :: maximum_value
        real(dp) :: normalizer

        allocate(values(size(x)))
        if (size(x) == 0) return
        if (any(ieee_is_nan(x))) then
            values = quiet_nan(0.0_dp)
            return
        end if

        maximum_value = maxval(x)
        if (.not. ieee_is_finite(maximum_value)) then
            if (maximum_value > 0.0_dp) then
                do i = 1, size(x)
                    if (x(i) > 0.0_dp .and. .not. ieee_is_finite(x(i))) then
                        values(i) = quiet_nan(0.0_dp)
                    else
                        values(i) = negative_infinity(0.0_dp)
                    end if
                end do
            else
                values = quiet_nan(0.0_dp)
            end if
            return
        end if

        values = x - maximum_value
        normalizer = log(sum(exp(values)))
        values = values - normalizer
    end function log_softmax_1d

    pure function log_softmax_2d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        real(dp), allocatable :: values(:, :)

        integer :: i
        integer :: j
        real(dp) :: maximum_value
        real(dp) :: normalizer

        allocate(values(size(x, 1), size(x, 2)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            if (any(ieee_is_nan(x))) then
                values = quiet_nan(0.0_dp)
                return
            end if
            maximum_value = maxval(x)
            if (.not. ieee_is_finite(maximum_value)) then
                if (maximum_value > 0.0_dp) then
                    do j = 1, size(x, 2)
                        do i = 1, size(x, 1)
                            if (x(i, j) > 0.0_dp .and. .not. ieee_is_finite(x(i, j))) then
                                values(i, j) = quiet_nan(0.0_dp)
                            else
                                values(i, j) = negative_infinity(0.0_dp)
                            end if
                        end do
                    end do
                else
                    values = quiet_nan(0.0_dp)
                end if
                return
            end if
            values = x - maximum_value
            normalizer = log(sum(exp(values)))
            values = values - normalizer
        else if (dim == 1) then
            do i = 1, size(x, 2)
                values(:, i) = log_softmax_1d(x(:, i))
            end do
        else if (dim == 2) then
            do i = 1, size(x, 1)
                values(i, :) = log_softmax_1d(x(i, :))
            end do
        else
            values = quiet_nan(0.0_dp)
        end if
    end function log_softmax_2d


    pure function log_softmax_2d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp), allocatable :: values(:, :)

        integer :: mode

        mode = reduction_dims_mode(dims)
        if (mode == 0) then
            allocate(values(size(x, 1), size(x, 2)))
            where (ieee_is_finite(x))
                values = 0.0_dp
            elsewhere
                values = quiet_nan(0.0_dp)
            end where
        else if (mode == 1) then
            values = log_softmax_2d(x, 1)
        else if (mode == 2) then
            values = log_softmax_2d(x, 2)
        else if (mode == 3) then
            values = log_softmax_2d(x)
        else
            allocate(values(size(x, 1), size(x, 2)))
            values = quiet_nan(0.0_dp)
        end if
    end function log_softmax_2d_dims


    pure function logsumexp_real_1d_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:) !! one-dimensional real values to reduce
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        value = logsumexp_complex_1d(cmplx(x, 0.0_dp, kind=dp), b)
    end function logsumexp_real_1d_complex_b

    pure function logsumexp_real_2d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :) !! two-dimensional real values reduced over both dimensions
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        value = logsumexp_complex_2d_all(cmplx(x, 0.0_dp, kind=dp), b)
    end function logsumexp_real_2d_all_complex_b

    pure function logsumexp_real_2d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional real values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:)

        values = logsumexp_complex_2d_dim(cmplx(x, 0.0_dp, kind=dp), dim, b)
    end function logsumexp_real_2d_dim_complex_b

    pure subroutine logsumexp_signed_real_1d_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:) !! one-dimensional real values to reduce
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_1d(cmplx(x, 0.0_dp, kind=dp), logabs, phase, b)
    end subroutine logsumexp_signed_real_1d_complex_b

    pure subroutine logsumexp_signed_real_2d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :) !! two-dimensional real values reduced over both dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_2d_all(cmplx(x, 0.0_dp, kind=dp), logabs, phase, b)
    end subroutine logsumexp_signed_real_2d_all_complex_b

    pure subroutine logsumexp_signed_real_2d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :) !! two-dimensional real values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        real(dp), allocatable, intent(out) :: logabs(:) !! log-magnitudes along the reduced dimension
        complex(dp), allocatable, intent(out) :: phase(:) !! unit phases along the reduced dimension
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_2d_dim(cmplx(x, 0.0_dp, kind=dp), dim, logabs, phase, b)
    end subroutine logsumexp_signed_real_2d_dim_complex_b

    pure function logsumexp_keepdims_real_2d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :) !! two-dimensional real values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :)

        values = logsumexp_keepdims_complex_2d(cmplx(x, 0.0_dp, kind=dp), dims, b)
    end function logsumexp_keepdims_real_2d_complex_b

    pure subroutine logsumexp_signed_keepdims_real_2d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :) !! two-dimensional real values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :) !! log-magnitudes with reduced dimensions retained
        complex(dp), allocatable, intent(out) :: phase(:, :) !! unit phases with reduced dimensions retained
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_keepdims_complex_2d(cmplx(x, 0.0_dp, kind=dp), dims, logabs, phase, b)
    end subroutine logsumexp_signed_keepdims_real_2d_complex_b

    pure function logsumexp_complex_1d(x, b) result(value)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex values to reduce
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights
        complex(dp) :: value

        complex(dp), allocatable :: terms(:)
        complex(dp), allocatable :: weights(:)
        complex(dp) :: phase
        logical :: ok
        real(dp) :: logabs

        if (present(b)) then
            call broadcast_complex_terms_1d(x, b, terms, weights, ok)
            if (.not. ok) then
                value = complex_nan()
                return
            end if
            call logsumexp_complex_core(terms, logabs, phase, weights)
        else
            call logsumexp_complex_core(x, logabs, phase)
        end if
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_1d

    pure function logsumexp_complex_1d_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex values to reduce
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: complex_b(:)
        complex(dp) :: value

        select rank (b)
        rank (0)
            value = logsumexp_complex_1d(x, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b(size(b)))
            complex_b = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_1d(x, complex_b)
        rank default
            value = complex_nan()
        end select
    end function logsumexp_complex_1d_real_b

    pure function logsumexp_complex_2d_all(x, b) result(value)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values reduced over both dimensions
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights
        complex(dp) :: value

        complex(dp), allocatable :: terms(:, :)
        complex(dp), allocatable :: weights(:, :)
        complex(dp) :: phase
        logical :: ok
        real(dp) :: logabs

        if (present(b)) then
            call broadcast_complex_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                value = complex_nan()
                return
            end if
            call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        end if
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_2d_all

    pure function logsumexp_complex_2d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values reduced over both dimensions
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp) :: value

        select rank (b)
        rank (0)
            value = logsumexp_complex_2d_all(x, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_2d_all(x, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_2d_all(x, complex_b2)
        rank default
            value = complex_nan()
        end select
    end function logsumexp_complex_2d_all_real_b

    pure function logsumexp_complex_2d_dim(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights
        complex(dp), allocatable :: values(:)

        complex(dp), allocatable :: terms(:, :)
        complex(dp), allocatable :: weights(:, :)
        complex(dp) :: phase
        integer :: i
        logical :: ok
        real(dp) :: logabs

        if (dim /= 1 .and. dim /= 2) then
            allocate(values(0))
            return
        end if

        if (present(b)) then
            call broadcast_complex_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                if (dim == 1) then
                    allocate(values(size(x, 2)))
                else
                    allocate(values(size(x, 1)))
                end if
                values = complex_nan()
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2)))
            terms = x
        end if

        if (dim == 1) then
            allocate(values(size(terms, 2)))
            do i = 1, size(terms, 2)
                if (present(b)) then
                    call logsumexp_complex_core(terms(:, i), logabs, phase, weights(:, i))
                else
                    call logsumexp_complex_core(terms(:, i), logabs, phase)
                end if
                values(i) = complex_log_from_phase(logabs, phase)
            end do
        else
            allocate(values(size(terms, 1)))
            do i = 1, size(terms, 1)
                if (present(b)) then
                    call logsumexp_complex_core(terms(i, :), logabs, phase, weights(i, :))
                else
                    call logsumexp_complex_core(terms(i, :), logabs, phase)
                end if
                values(i) = complex_log_from_phase(logabs, phase)
            end do
        end if
    end function logsumexp_complex_2d_dim

    pure function logsumexp_complex_2d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: values(:)

        if (dim /= 1 .and. dim /= 2) then
            allocate(values(0))
            return
        end if
        select rank (b)
        rank (0)
            values = logsumexp_complex_2d_dim(x, dim, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_2d_dim(x, dim, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_2d_dim(x, dim, complex_b2)
        rank default
            if (dim == 1) then
                allocate(values(size(x, 2)))
            else
                allocate(values(size(x, 1)))
            end if
            values = complex_nan()
        end select
    end function logsumexp_complex_2d_dim_real_b

    pure subroutine logsumexp_signed_complex_1d(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex values to reduce
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:)
        complex(dp), allocatable :: weights(:)
        logical :: ok

        if (present(b)) then
            call broadcast_complex_terms_1d(x, b, terms, weights, ok)
            if (.not. ok) then
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
            call logsumexp_complex_core(terms, logabs, phase, weights)
        else
            call logsumexp_complex_core(x, logabs, phase)
        end if
    end subroutine logsumexp_signed_complex_1d

    pure subroutine logsumexp_signed_complex_1d_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex values to reduce
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_1d(x, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_1d(x, logabs, phase, complex_b1)
        rank default
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_1d_real_b

    pure subroutine logsumexp_signed_complex_2d_all(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values reduced over both dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:, :)
        complex(dp), allocatable :: weights(:, :)
        logical :: ok

        if (present(b)) then
            call broadcast_complex_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
            call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        end if
    end subroutine logsumexp_signed_complex_2d_all

    pure subroutine logsumexp_signed_complex_2d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values reduced over both dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_2d_all(x, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_2d_all(x, logabs, phase, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_2d_all(x, logabs, phase, complex_b2)
        rank default
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_2d_all_real_b

    pure subroutine logsumexp_signed_complex_2d_dim(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        real(dp), allocatable, intent(out) :: logabs(:) !! log-magnitudes along the reduced dimension
        complex(dp), allocatable, intent(out) :: phase(:) !! unit phases along the reduced dimension
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:, :)
        complex(dp), allocatable :: weights(:, :)
        integer :: i
        logical :: ok

        if (dim /= 1 .and. dim /= 2) then
            allocate(logabs(0), phase(0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                if (dim == 1) then
                    allocate(logabs(size(x, 2)), phase(size(x, 2)))
                else
                    allocate(logabs(size(x, 1)), phase(size(x, 1)))
                end if
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2)))
            terms = x
        end if

        if (dim == 1) then
            allocate(logabs(size(terms, 2)), phase(size(terms, 2)))
            do i = 1, size(terms, 2)
                if (present(b)) then
                    call logsumexp_complex_core(terms(:, i), logabs(i), phase(i), weights(:, i))
                else
                    call logsumexp_complex_core(terms(:, i), logabs(i), phase(i))
                end if
            end do
        else
            allocate(logabs(size(terms, 1)), phase(size(terms, 1)))
            do i = 1, size(terms, 1)
                if (present(b)) then
                    call logsumexp_complex_core(terms(i, :), logabs(i), phase(i), weights(i, :))
                else
                    call logsumexp_complex_core(terms(i, :), logabs(i), phase(i))
                end if
            end do
        end if
    end subroutine logsumexp_signed_complex_2d_dim

    pure subroutine logsumexp_signed_complex_2d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, either 1 or 2
        real(dp), allocatable, intent(out) :: logabs(:) !! log-magnitudes along the reduced dimension
        complex(dp), allocatable, intent(out) :: phase(:) !! unit phases along the reduced dimension
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)

        if (dim /= 1 .and. dim /= 2) then
            allocate(logabs(0), phase(0))
            return
        end if
        select rank (b)
        rank (0)
            call logsumexp_signed_complex_2d_dim(x, dim, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_2d_dim(x, dim, logabs, phase, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_2d_dim(x, dim, logabs, phase, complex_b2)
        rank default
            if (dim == 1) then
                allocate(logabs(size(x, 2)), phase(size(x, 2)))
            else
                allocate(logabs(size(x, 1)), phase(size(x, 1)))
            end if
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_2d_dim_real_b

    pure function logsumexp_keepdims_complex_2d(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights
        complex(dp), allocatable :: values(:, :)

        complex(dp), allocatable :: terms(:, :)
        complex(dp), allocatable :: weights(:, :)
        logical :: ok

        if (reduction_dims_mode(dims) < 0) then
            allocate(values(0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_complex(values, size(x, 1), size(x, 2), &
                    reduction_dims_mode(dims))
                if (size(values) > 0) values = complex_nan()
                return
            end if
            values = logsumexp_keepdims_complex_2d_full(terms, dims, weights)
        else
            values = logsumexp_keepdims_complex_2d_full(x, dims)
        end if
    end function logsumexp_keepdims_complex_2d

    pure function logsumexp_keepdims_complex_2d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: values(:, :)

        select rank (b)
        rank (0)
            values = logsumexp_keepdims_complex_2d(x, dims, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_2d(x, dims, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_2d(x, dims, complex_b2)
        rank default
            call allocate_keepdims_complex(values, size(x, 1), size(x, 2), &
                reduction_dims_mode(dims))
            if (size(values) > 0) values = complex_nan()
        end select
    end function logsumexp_keepdims_complex_2d_real_b

    pure function logsumexp_keepdims_complex_2d_full(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in), optional :: b(:, :) !! optional fully broadcast complex weights
        complex(dp), allocatable :: values(:, :)

        complex(dp) :: phase
        integer :: i
        integer :: j
        integer :: mode
        real(dp) :: logabs

        mode = reduction_dims_mode(dims)
        if (mode < 0) then
            allocate(values(0, 0))
            return
        end if

        select case (mode)
        case (0)
            allocate(values(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    if (present(b)) then
                        call logsumexp_complex_core([x(i, j)], logabs, phase, [b(i, j)])
                    else
                        call logsumexp_complex_core([x(i, j)], logabs, phase)
                    end if
                    values(i, j) = complex_log_from_phase(logabs, phase)
                end do
            end do
        case (1)
            allocate(values(1, size(x, 2)))
            if (present(b)) then
                values(1, :) = logsumexp_complex_2d_dim(x, 1, b)
            else
                values(1, :) = logsumexp_complex_2d_dim(x, 1)
            end if
        case (2)
            allocate(values(size(x, 1), 1))
            if (present(b)) then
                values(:, 1) = logsumexp_complex_2d_dim(x, 2, b)
            else
                values(:, 1) = logsumexp_complex_2d_dim(x, 2)
            end if
        case (3)
            allocate(values(1, 1))
            if (present(b)) then
                values(1, 1) = logsumexp_complex_2d_all(x, b)
            else
                values(1, 1) = logsumexp_complex_2d_all(x)
            end if
        end select
    end function logsumexp_keepdims_complex_2d_full

    pure subroutine logsumexp_signed_keepdims_complex_2d(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :) !! log-magnitudes with reduced dimensions retained
        complex(dp), allocatable, intent(out) :: phase(:, :) !! unit phases with reduced dimensions retained
        complex(dp), intent(in), optional :: b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:, :)
        complex(dp), allocatable :: weights(:, :)
        logical :: ok

        if (reduction_dims_mode(dims) < 0) then
            allocate(logabs(0, 0), phase(0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_2d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_real(logabs, size(x, 1), size(x, 2), reduction_dims_mode(dims))
                call allocate_keepdims_complex(phase, size(x, 1), size(x, 2), reduction_dims_mode(dims))
                if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
                if (size(phase) > 0) phase = complex_nan()
                return
            end if
            call logsumexp_signed_keepdims_complex_2d_full(terms, dims, logabs, phase, weights)
        else
            call logsumexp_signed_keepdims_complex_2d_full(x, dims, logabs, phase)
        end if
    end subroutine logsumexp_signed_keepdims_complex_2d

    pure subroutine logsumexp_signed_keepdims_complex_2d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :) !! log-magnitudes with reduced dimensions retained
        complex(dp), allocatable, intent(out) :: phase(:, :) !! unit phases with reduced dimensions retained
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_keepdims_complex_2d(x, dims, logabs, phase, &
                cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_2d(x, dims, logabs, phase, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_2d(x, dims, logabs, phase, complex_b2)
        rank default
            call allocate_keepdims_real(logabs, size(x, 1), size(x, 2), reduction_dims_mode(dims))
            call allocate_keepdims_complex(phase, size(x, 1), size(x, 2), reduction_dims_mode(dims))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
        end select
    end subroutine logsumexp_signed_keepdims_complex_2d_real_b

    pure subroutine logsumexp_signed_keepdims_complex_2d_full(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :) !! log-magnitudes with reduced dimensions retained
        complex(dp), allocatable, intent(out) :: phase(:, :) !! unit phases with reduced dimensions retained
        complex(dp), intent(in), optional :: b(:, :) !! optional fully broadcast complex weights

        complex(dp), allocatable :: reduced_phase(:)
        integer :: i
        integer :: j
        integer :: mode
        real(dp), allocatable :: reduced_logabs(:)

        mode = reduction_dims_mode(dims)
        if (mode < 0) then
            allocate(logabs(0, 0), phase(0, 0))
            return
        end if

        select case (mode)
        case (0)
            allocate(logabs(size(x, 1), size(x, 2)))
            allocate(phase(size(x, 1), size(x, 2)))
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    if (present(b)) then
                        call logsumexp_complex_core([x(i, j)], logabs(i, j), phase(i, j), [b(i, j)])
                    else
                        call logsumexp_complex_core([x(i, j)], logabs(i, j), phase(i, j))
                    end if
                end do
            end do
        case (1)
            allocate(logabs(1, size(x, 2)), phase(1, size(x, 2)))
            if (present(b)) then
                call logsumexp_signed_complex_2d_dim(x, 1, reduced_logabs, reduced_phase, b)
            else
                call logsumexp_signed_complex_2d_dim(x, 1, reduced_logabs, reduced_phase)
            end if
            logabs(1, :) = reduced_logabs
            phase(1, :) = reduced_phase
        case (2)
            allocate(logabs(size(x, 1), 1), phase(size(x, 1), 1))
            if (present(b)) then
                call logsumexp_signed_complex_2d_dim(x, 2, reduced_logabs, reduced_phase, b)
            else
                call logsumexp_signed_complex_2d_dim(x, 2, reduced_logabs, reduced_phase)
            end if
            logabs(:, 1) = reduced_logabs
            phase(:, 1) = reduced_phase
        case (3)
            allocate(logabs(1, 1), phase(1, 1))
            if (present(b)) then
                call logsumexp_signed_complex_2d_all(x, logabs(1, 1), phase(1, 1), b)
            else
                call logsumexp_signed_complex_2d_all(x, logabs(1, 1), phase(1, 1))
            end if
        end select
    end subroutine logsumexp_signed_keepdims_complex_2d_full

    pure subroutine broadcast_complex_terms_1d(x, b, terms, weights, ok)
        complex(dp), intent(in) :: x(:) !! one-dimensional logarithmic terms to broadcast
        complex(dp), intent(in) :: b(..) !! scalar or rank-1 weights broadcast jointly with x
        complex(dp), allocatable, intent(out) :: terms(:) !! input terms broadcast to the joint shape
        complex(dp), allocatable, intent(out) :: weights(:) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: i
        integer :: ib
        integer :: ix
        integer :: n

        ok = .true.
        select rank (b)
        rank (0)
            n = size(x)
            allocate(terms(n), weights(n))
            terms = x
            weights = b
        rank (1)
            if (size(x) /= size(b) .and. size(x) /= 1 .and. size(b) /= 1) then
                ok = .false.
                allocate(terms(0), weights(0))
                return
            end if
            n = size(x)
            if (n == 1) n = size(b)
            allocate(terms(n), weights(n))
            do i = 1, n
                ix = merge(1, i, size(x) == 1)
                ib = merge(1, i, size(b) == 1)
                terms(i) = x(ix)
                weights(i) = b(ib)
            end do
        rank default
            ok = .false.
            allocate(terms(0), weights(0))
        end select
    end subroutine broadcast_complex_terms_1d

    pure subroutine broadcast_complex_terms_2d(x, b, terms, weights, ok)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional logarithmic terms to broadcast
        complex(dp), intent(in) :: b(..) !! scalar, trailing-vector, or rank-2 weights
        complex(dp), allocatable, intent(out) :: terms(:, :) !! input terms broadcast to the joint shape
        complex(dp), allocatable, intent(out) :: weights(:, :) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: n1
        integer :: n2

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2)), weights(size(x, 1), size(x, 2)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = size(b)
        rank (2)
            b1 = size(b, 1)
            b2 = size(b, 2)
        rank default
            ok = .false.
            allocate(terms(0, 0), weights(0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1)) then
            ok = .false.
            allocate(terms(0, 0), weights(0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        allocate(terms(n1, n2), weights(n1, n2))
        select rank (b)
        rank (1)
            do j = 1, n2
                jx = merge(1, j, size(x, 2) == 1)
                jb = merge(1, j, size(b) == 1)
                do i = 1, n1
                    ix = merge(1, i, size(x, 1) == 1)
                    terms(i, j) = x(ix, jx)
                    weights(i, j) = b(jb)
                end do
            end do
        rank (2)
            do j = 1, n2
                jx = merge(1, j, size(x, 2) == 1)
                jb = merge(1, j, size(b, 2) == 1)
                do i = 1, n1
                    ix = merge(1, i, size(x, 1) == 1)
                    ib = merge(1, i, size(b, 1) == 1)
                    terms(i, j) = x(ix, jx)
                    weights(i, j) = b(ib, jb)
                end do
            end do
        end select
    end subroutine broadcast_complex_terms_2d

    pure subroutine allocate_keepdims_complex(values, n1, n2, mode)
        complex(dp), allocatable, intent(out) :: values(:, :) !! allocated keepdims result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: mode !! reduction mode from reduction_dims_mode

        if (mode < 0) then
            allocate(values(0, 0))
        else if (mode == 0) then
            allocate(values(n1, n2))
        else if (mode == 1) then
            allocate(values(1, n2))
        else if (mode == 2) then
            allocate(values(n1, 1))
        else
            allocate(values(1, 1))
        end if
    end subroutine allocate_keepdims_complex

    pure subroutine logsumexp_complex_core(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in), optional :: b(:) !! optional complex weights multiplying exp(x)

        integer :: i
        logical :: active
        real(dp) :: imag_compensation
        real(dp) :: imag_sum
        real(dp) :: maximum_term
        real(dp) :: magnitude
        real(dp) :: real_compensation
        real(dp) :: real_sum
        real(dp) :: scaled_magnitude
        real(dp) :: term_log
        complex(dp) :: term_phase
        complex(dp) :: weight

        if (size(x) == 0) then
            logabs = negative_infinity(0.0_dp)
            phase = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        maximum_term = negative_infinity(0.0_dp)
        do i = 1, size(x)
            if (present(b)) then
                weight = b(i)
                if (complex_has_nan(weight)) then
                    logabs = quiet_nan(0.0_dp)
                    phase = complex_nan()
                    return
                end if
                active = abs(weight) /= 0.0_dp
            else
                weight = cmplx(1.0_dp, 0.0_dp, kind=dp)
                active = .true.
            end if
            if (.not. active) cycle
            if (complex_has_nan(x(i)) .or. .not. ieee_is_finite(aimag(x(i)))) then
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
            if (present(b)) then
                magnitude = abs(weight)
                if (.not. ieee_is_finite(magnitude)) then
                    logabs = quiet_nan(0.0_dp)
                    phase = complex_nan()
                    return
                end if
                term_log = real(x(i), dp) + log(magnitude)
            else
                term_log = real(x(i), dp)
            end if
            if (ieee_is_nan(term_log)) then
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
            maximum_term = max(maximum_term, term_log)
        end do

        if (maximum_term == negative_infinity(0.0_dp)) then
            logabs = maximum_term
            phase = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if
        if (.not. ieee_is_finite(maximum_term)) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if

        real_sum = 0.0_dp
        imag_sum = 0.0_dp
        real_compensation = 0.0_dp
        imag_compensation = 0.0_dp
        do i = 1, size(x)
            if (present(b)) then
                weight = b(i)
                if (abs(weight) == 0.0_dp) cycle
                magnitude = abs(weight)
                term_log = real(x(i), dp) + log(magnitude)
                term_phase = weight / magnitude
            else
                term_log = real(x(i), dp)
                term_phase = cmplx(1.0_dp, 0.0_dp, kind=dp)
            end if
            term_phase = term_phase * exp(cmplx(0.0_dp, aimag(x(i)), kind=dp))
            scaled_magnitude = exp(term_log - maximum_term)
            call compensated_add(real_sum, real_compensation, &
                scaled_magnitude * real(term_phase, dp))
            call compensated_add(imag_sum, imag_compensation, &
                scaled_magnitude * aimag(term_phase))
        end do

        magnitude = hypot(real_sum, imag_sum)
        if (magnitude == 0.0_dp) then
            logabs = negative_infinity(0.0_dp)
            phase = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else if (.not. ieee_is_finite(magnitude)) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        else
            logabs = maximum_term + log(magnitude)
            phase = cmplx(real_sum / magnitude, imag_sum / magnitude, kind=dp)
        end if
    end subroutine logsumexp_complex_core

    pure function complex_log_from_phase(logabs, phase) result(value)
        real(dp), intent(in) :: logabs !! logarithm of a complex result magnitude
        complex(dp), intent(in) :: phase !! unit phase, or zero for an exact zero result
        complex(dp) :: value

        if (ieee_is_nan(logabs) .or. complex_has_nan(phase)) then
            value = complex_nan()
        else if (abs(phase) == 0.0_dp) then
            value = cmplx(logabs, 0.0_dp, kind=dp)
        else
            value = cmplx(logabs, atan2(aimag(phase), real(phase, dp)), kind=dp)
        end if
    end function complex_log_from_phase

    pure function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

    pure function complex_has_nan(value) result(has_nan)
        complex(dp), intent(in) :: value !! complex value inspected componentwise
        logical :: has_nan

        has_nan = ieee_is_nan(real(value, dp)) .or. ieee_is_nan(aimag(value))
    end function complex_has_nan

    pure function softmax_complex_1d(x) result(values)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex values to normalize
        complex(dp), allocatable :: values(:)

        complex(dp) :: normalizer

        allocate(values(size(x)))
        if (size(x) == 0) return
        normalizer = logsumexp_complex_1d(x)
        if (complex_has_nan(normalizer)) then
            values = complex_nan()
            return
        end if
        values = exp(x - normalizer)
    end function softmax_complex_1d

    pure function softmax_complex_2d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        complex(dp), allocatable :: values(:, :)

        integer :: i
        complex(dp) :: normalizer

        allocate(values(size(x, 1), size(x, 2)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            normalizer = logsumexp_complex_2d_all(x)
            if (complex_has_nan(normalizer)) then
                values = complex_nan()
            else
                values = exp(x - normalizer)
            end if
        else if (dim == 1) then
            do i = 1, size(x, 2)
                values(:, i) = softmax_complex_1d(x(:, i))
            end do
        else if (dim == 2) then
            do i = 1, size(x, 1)
                values(i, :) = softmax_complex_1d(x(i, :))
            end do
        else
            values = complex_nan()
        end if
    end function softmax_complex_2d

    pure function softmax_complex_2d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp), allocatable :: values(:, :)

        integer :: mode

        mode = reduction_dims_mode(dims)
        if (mode == 0) then
            allocate(values(size(x, 1), size(x, 2)))
            where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                values = cmplx(1.0_dp, 0.0_dp, kind=dp)
            elsewhere
                values = complex_nan()
            end where
        else if (mode == 1) then
            values = softmax_complex_2d(x, 1)
        else if (mode == 2) then
            values = softmax_complex_2d(x, 2)
        else if (mode == 3) then
            values = softmax_complex_2d(x)
        else
            allocate(values(size(x, 1), size(x, 2)))
            values = complex_nan()
        end if
    end function softmax_complex_2d_dims

    ! Complex max-shift branch semantics follow SciPy 1.17.0 _logsumexp.py.
    ! SciPy's BSD-3-Clause notice is reproduced in THIRD_PARTY_LICENSES.md.
    pure function log_softmax_complex_1d(x) result(values)
        complex(dp), intent(in) :: x(:) !! one-dimensional complex values to log-normalize
        complex(dp), allocatable :: values(:)

        integer :: i
        complex(dp) :: maximum_value
        complex(dp) :: normalizer
        complex(dp), allocatable :: shifted(:)

        allocate(values(size(x)))
        if (size(x) == 0) return

        maximum_value = x(1)
        do i = 2, size(x)
            if (complex_has_nan(x(i))) then
                maximum_value = complex_nan()
                exit
            end if
            if (real(x(i), dp) > real(maximum_value, dp)) then
                maximum_value = x(i)
            else if (real(x(i), dp) < real(maximum_value, dp)) then
                cycle
            else if (aimag(x(i)) > aimag(maximum_value)) then
                maximum_value = x(i)
            end if
        end do
        if (.not. ieee_is_finite(real(maximum_value, dp)) .or. &
                .not. ieee_is_finite(aimag(maximum_value))) then
            maximum_value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end if

        allocate(shifted(size(x)))
        shifted = x - maximum_value
        normalizer = logsumexp_complex_1d(shifted)
        if (complex_has_nan(normalizer)) then
            values = complex_nan()
            return
        end if
        values = shifted - normalizer
    end function log_softmax_complex_1d

    pure function log_softmax_complex_2d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to log-normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        complex(dp), allocatable :: values(:, :)

        integer :: i
        complex(dp) :: normalizer

        allocate(values(size(x, 1), size(x, 2)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            normalizer = logsumexp_complex_2d_all(x)
            if (complex_has_nan(normalizer)) then
                values = complex_nan()
            else
                values = x - normalizer
            end if
        else if (dim == 1) then
            do i = 1, size(x, 2)
                values(:, i) = log_softmax_complex_1d(x(:, i))
            end do
        else if (dim == 2) then
            do i = 1, size(x, 1)
                values(i, :) = log_softmax_complex_1d(x(i, :))
            end do
        else
            values = complex_nan()
        end if
    end function log_softmax_complex_2d


    pure function log_softmax_complex_2d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :) !! two-dimensional complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp), allocatable :: values(:, :)

        integer :: mode

        mode = reduction_dims_mode(dims)
        if (mode == 0) then
            allocate(values(size(x, 1), size(x, 2)))
            where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                values = cmplx(0.0_dp, 0.0_dp, kind=dp)
            elsewhere
                values = complex_nan()
            end where
        else if (mode == 1) then
            values = log_softmax_complex_2d(x, 1)
        else if (mode == 2) then
            values = log_softmax_complex_2d(x, 2)
        else if (mode == 3) then
            values = log_softmax_complex_2d(x)
        else
            allocate(values(size(x, 1), size(x, 2)))
            values = complex_nan()
        end if
    end function log_softmax_complex_2d_dims


    pure function logsumexp_3d_all(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values reduced over all dimensions
        real(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable exponential weights
        real(dp) :: value

        real(dp), allocatable :: terms(:, :, :)
        real(dp), allocatable :: weights(:, :, :)
        real(dp) :: logabs
        real(dp) :: sign_value
        logical :: ok

        if (present(b)) then
            call broadcast_real_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
        end if

        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = logabs
        end if
    end function logsumexp_3d_all

    pure function logsumexp_3d_dim(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        real(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable exponential weights
        real(dp), allocatable :: values(:, :)

        integer :: i
        integer :: j
        integer :: k
        logical :: ok
        real(dp), allocatable :: terms(:, :, :)
        real(dp), allocatable :: weights(:, :, :)
        real(dp) :: logabs
        real(dp) :: sign_value

        if (dim < 1 .or. dim > 3) then
            allocate(values(0, 0))
            return
        end if

        if (present(b)) then
            call broadcast_real_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                select case (dim)
                case (1)
                    allocate(values(size(x, 2), size(x, 3)))
                case (2)
                    allocate(values(size(x, 1), size(x, 3)))
                case (3)
                    allocate(values(size(x, 1), size(x, 2)))
                end select
                values = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2), size(x, 3)))
            terms = x
        end if

        select case (dim)
        case (1)
            allocate(values(size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    if (present(b)) then
                        call logsumexp_core(terms(:, j, k), logabs, sign_value, weights(:, j, k))
                    else
                        call logsumexp_core(terms(:, j, k), logabs, sign_value)
                    end if
                    if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                        values(j, k) = quiet_nan(0.0_dp)
                    else
                        values(j, k) = logabs
                    end if
                end do
            end do
        case (2)
            allocate(values(size(terms, 1), size(terms, 3)))
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_core(terms(i, :, k), logabs, sign_value, weights(i, :, k))
                    else
                        call logsumexp_core(terms(i, :, k), logabs, sign_value)
                    end if
                    if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                        values(i, k) = quiet_nan(0.0_dp)
                    else
                        values(i, k) = logabs
                    end if
                end do
            end do
        case (3)
            allocate(values(size(terms, 1), size(terms, 2)))
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_core(terms(i, j, :), logabs, sign_value, weights(i, j, :))
                    else
                        call logsumexp_core(terms(i, j, :), logabs, sign_value)
                    end if
                    if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                        values(i, j) = quiet_nan(0.0_dp)
                    else
                        values(i, j) = logabs
                    end if
                end do
            end do
        end select
    end function logsumexp_3d_dim

    pure subroutine logsumexp_signed_3d_all(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of the weighted exponential sum: -1, 0, or +1
        real(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable exponential weights

        logical :: ok
        real(dp), allocatable :: terms(:, :, :)
        real(dp), allocatable :: weights(:, :, :)

        if (present(b)) then
            call broadcast_real_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
        end if
    end subroutine logsumexp_signed_3d_all

    pure subroutine logsumexp_signed_3d_dim(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        real(dp), allocatable, intent(out) :: logabs(:, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :) !! signs along the reduced dimension
        real(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable exponential weights

        integer :: i
        integer :: j
        integer :: k
        logical :: ok
        real(dp), allocatable :: terms(:, :, :)
        real(dp), allocatable :: weights(:, :, :)

        if (dim < 1 .or. dim > 3) then
            allocate(logabs(0, 0), sign_value(0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_real_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                select case (dim)
                case (1)
                    allocate(logabs(size(x, 2), size(x, 3)), sign_value(size(x, 2), size(x, 3)))
                case (2)
                    allocate(logabs(size(x, 1), size(x, 3)), sign_value(size(x, 1), size(x, 3)))
                case (3)
                    allocate(logabs(size(x, 1), size(x, 2)), sign_value(size(x, 1), size(x, 2)))
                end select
                logabs = quiet_nan(0.0_dp)
                sign_value = quiet_nan(0.0_dp)
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2), size(x, 3)))
            terms = x
        end if

        select case (dim)
        case (1)
            allocate(logabs(size(terms, 2), size(terms, 3)), sign_value(size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    if (present(b)) then
                        call logsumexp_core(terms(:, j, k), logabs(j, k), sign_value(j, k), &
                            weights(:, j, k))
                    else
                        call logsumexp_core(terms(:, j, k), logabs(j, k), sign_value(j, k))
                    end if
                end do
            end do
        case (2)
            allocate(logabs(size(terms, 1), size(terms, 3)), sign_value(size(terms, 1), size(terms, 3)))
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_core(terms(i, :, k), logabs(i, k), sign_value(i, k), &
                            weights(i, :, k))
                    else
                        call logsumexp_core(terms(i, :, k), logabs(i, k), sign_value(i, k))
                    end if
                end do
            end do
        case (3)
            allocate(logabs(size(terms, 1), size(terms, 2)), sign_value(size(terms, 1), size(terms, 2)))
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_core(terms(i, j, :), logabs(i, j), sign_value(i, j), &
                            weights(i, j, :))
                    else
                        call logsumexp_core(terms(i, j, :), logabs(i, j), sign_value(i, j))
                    end if
                end do
            end do
        end select
    end subroutine logsumexp_signed_3d_dim

    pure function logsumexp_keepdims_3d(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable exponential weights
        real(dp), allocatable :: values(:, :, :)

        logical :: ok
        real(dp), allocatable :: terms(:, :, :)
        real(dp), allocatable :: weights(:, :, :)

        if (reduction_dims_mode_3d(dims) < 0) then
            allocate(values(0, 0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_real_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_real_3d(values, size(x, 1), size(x, 2), size(x, 3), &
                    reduction_dims_mode_3d(dims))
                if (size(values) > 0) values = quiet_nan(0.0_dp)
                return
            end if
            values = logsumexp_keepdims_3d_full(terms, dims, weights)
        else
            values = logsumexp_keepdims_3d_full(x, dims)
        end if
    end function logsumexp_keepdims_3d

    pure function logsumexp_keepdims_3d_full(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in), optional :: b(:, :, :) !! optional fully broadcast exponential weights
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: mode
        real(dp) :: logabs
        real(dp) :: sign_value

        mode = reduction_dims_mode_3d(dims)
        if (mode < 0) then
            allocate(values(0, 0, 0))
            return
        end if
        call allocate_keepdims_real_3d(values, size(x, 1), size(x, 2), size(x, 3), mode)

        select case (mode)
        case (0)
            if (.not. present(b)) then
                values = x
            else
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        do i = 1, size(x, 1)
                            call logsumexp_core([x(i, j, k)], logabs, sign_value, [b(i, j, k)])
                            if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                                values(i, j, k) = quiet_nan(0.0_dp)
                            else
                                values(i, j, k) = logabs
                            end if
                        end do
                    end do
                end do
            end if
        case (1, 2, 4)
            if (mode == 1) then
                if (present(b)) then
                    values(1, :, :) = logsumexp_3d_dim(x, 1, b)
                else
                    values(1, :, :) = logsumexp_3d_dim(x, 1)
                end if
            else if (mode == 2) then
                if (present(b)) then
                    values(:, 1, :) = logsumexp_3d_dim(x, 2, b)
                else
                    values(:, 1, :) = logsumexp_3d_dim(x, 2)
                end if
            else
                if (present(b)) then
                    values(:, :, 1) = logsumexp_3d_dim(x, 3, b)
                else
                    values(:, :, 1) = logsumexp_3d_dim(x, 3)
                end if
            end if
        case (3)
            do k = 1, size(x, 3)
                if (present(b)) then
                    values(1, 1, k) = logsumexp_2d_all(x(:, :, k), b(:, :, k))
                else
                    values(1, 1, k) = logsumexp_2d_all(x(:, :, k))
                end if
            end do
        case (5)
            do j = 1, size(x, 2)
                if (present(b)) then
                    values(1, j, 1) = logsumexp_2d_all(x(:, j, :), b(:, j, :))
                else
                    values(1, j, 1) = logsumexp_2d_all(x(:, j, :))
                end if
            end do
        case (6)
            do i = 1, size(x, 1)
                if (present(b)) then
                    values(i, 1, 1) = logsumexp_2d_all(x(i, :, :), b(i, :, :))
                else
                    values(i, 1, 1) = logsumexp_2d_all(x(i, :, :))
                end if
            end do
        case (7)
            if (present(b)) then
                values(1, 1, 1) = logsumexp_3d_all(x, b)
            else
                values(1, 1, 1) = logsumexp_3d_all(x)
            end if
        end select
    end function logsumexp_keepdims_3d_full

    pure subroutine logsumexp_signed_keepdims_3d(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: &
            logabs(:, :, :) !! log-absolute sums with reduced dimensions retained
        real(dp), allocatable, intent(out) :: sign_value(:, :, :) !! signs with kept dimensions
        real(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable exponential weights

        logical :: ok
        real(dp), allocatable :: terms(:, :, :)
        real(dp), allocatable :: weights(:, :, :)

        if (reduction_dims_mode_3d(dims) < 0) then
            allocate(logabs(0, 0, 0), sign_value(0, 0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_real_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_real_3d(logabs, size(x, 1), size(x, 2), size(x, 3), &
                    reduction_dims_mode_3d(dims))
                call allocate_keepdims_real_3d(sign_value, size(x, 1), size(x, 2), size(x, 3), &
                    reduction_dims_mode_3d(dims))
                if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
                if (size(sign_value) > 0) sign_value = quiet_nan(0.0_dp)
                return
            end if
            call logsumexp_signed_keepdims_3d_full(terms, dims, logabs, sign_value, weights)
        else
            call logsumexp_signed_keepdims_3d_full(x, dims, logabs, sign_value)
        end if
    end subroutine logsumexp_signed_keepdims_3d

    pure subroutine logsumexp_signed_keepdims_3d_full(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: &
            logabs(:, :, :) !! log-absolute sums with reduced dimensions retained
        real(dp), allocatable, intent(out) :: sign_value(:, :, :) !! signs with kept dimensions
        real(dp), intent(in), optional :: b(:, :, :) !! optional fully broadcast exponential weights

        integer :: i
        integer :: j
        integer :: k
        integer :: mode
        real(dp), allocatable :: reduced_logabs(:, :)
        real(dp), allocatable :: reduced_sign(:, :)

        mode = reduction_dims_mode_3d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0), sign_value(0, 0, 0))
            return
        end if
        call allocate_keepdims_real_3d(logabs, size(x, 1), size(x, 2), size(x, 3), mode)
        call allocate_keepdims_real_3d(sign_value, size(x, 1), size(x, 2), size(x, 3), mode)

        select case (mode)
        case (0)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        if (present(b)) then
                            call logsumexp_core([x(i, j, k)], logabs(i, j, k), &
                                sign_value(i, j, k), [b(i, j, k)])
                        else
                            call logsumexp_core([x(i, j, k)], logabs(i, j, k), &
                                sign_value(i, j, k))
                        end if
                    end do
                end do
            end do
        case (1, 2, 4)
            if (mode == 1) then
                if (present(b)) then
                    call logsumexp_signed_3d_dim(x, 1, reduced_logabs, reduced_sign, b)
                else
                    call logsumexp_signed_3d_dim(x, 1, reduced_logabs, reduced_sign)
                end if
                logabs(1, :, :) = reduced_logabs
                sign_value(1, :, :) = reduced_sign
            else if (mode == 2) then
                if (present(b)) then
                    call logsumexp_signed_3d_dim(x, 2, reduced_logabs, reduced_sign, b)
                else
                    call logsumexp_signed_3d_dim(x, 2, reduced_logabs, reduced_sign)
                end if
                logabs(:, 1, :) = reduced_logabs
                sign_value(:, 1, :) = reduced_sign
            else
                if (present(b)) then
                    call logsumexp_signed_3d_dim(x, 3, reduced_logabs, reduced_sign, b)
                else
                    call logsumexp_signed_3d_dim(x, 3, reduced_logabs, reduced_sign)
                end if
                logabs(:, :, 1) = reduced_logabs
                sign_value(:, :, 1) = reduced_sign
            end if
        case (3)
            do k = 1, size(x, 3)
                if (present(b)) then
                    call logsumexp_signed_2d_all(x(:, :, k), logabs(1, 1, k), &
                        sign_value(1, 1, k), b(:, :, k))
                else
                    call logsumexp_signed_2d_all(x(:, :, k), logabs(1, 1, k), &
                        sign_value(1, 1, k))
                end if
            end do
        case (5)
            do j = 1, size(x, 2)
                if (present(b)) then
                    call logsumexp_signed_2d_all(x(:, j, :), logabs(1, j, 1), &
                        sign_value(1, j, 1), b(:, j, :))
                else
                    call logsumexp_signed_2d_all(x(:, j, :), logabs(1, j, 1), &
                        sign_value(1, j, 1))
                end if
            end do
        case (6)
            do i = 1, size(x, 1)
                if (present(b)) then
                    call logsumexp_signed_2d_all(x(i, :, :), logabs(i, 1, 1), &
                        sign_value(i, 1, 1), b(i, :, :))
                else
                    call logsumexp_signed_2d_all(x(i, :, :), logabs(i, 1, 1), &
                        sign_value(i, 1, 1))
                end if
            end do
        case (7)
            if (present(b)) then
                call logsumexp_signed_3d_all(x, logabs(1, 1, 1), sign_value(1, 1, 1), b)
            else
                call logsumexp_signed_3d_all(x, logabs(1, 1, 1), sign_value(1, 1, 1))
            end if
        end select
    end subroutine logsumexp_signed_keepdims_3d_full

    pure subroutine broadcast_real_terms_3d(x, b, terms, weights, ok)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional logarithmic terms to broadcast
        real(dp), intent(in) :: b(..) !! scalar or trailing-array weights broadcast jointly with x
        real(dp), allocatable, intent(out) :: terms(:, :, :) !! input terms broadcast to the joint shape
        real(dp), allocatable, intent(out) :: weights(:, :, :) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: b3
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: k
        integer :: kb
        integer :: kx
        integer :: n1
        integer :: n2
        integer :: n3

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2), size(x, 3)), &
                weights(size(x, 1), size(x, 2), size(x, 3)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = 1
            b3 = size(b)
        rank (2)
            b1 = 1
            b2 = size(b, 1)
            b3 = size(b, 2)
        rank (3)
            b1 = size(b, 1)
            b2 = size(b, 2)
            b3 = size(b, 3)
        rank default
            ok = .false.
            allocate(terms(0, 0, 0), weights(0, 0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1) .or. &
                (size(x, 3) /= b3 .and. size(x, 3) /= 1 .and. b3 /= 1)) then
            ok = .false.
            allocate(terms(0, 0, 0), weights(0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        if (n3 == 1) n3 = b3
        allocate(terms(n1, n2, n3), weights(n1, n2, n3))
        select rank (b)
        rank (1)
            do k = 1, n3
                kx = merge(1, k, size(x, 3) == 1)
                kb = merge(1, k, size(b) == 1)
                do j = 1, n2
                    jx = merge(1, j, size(x, 2) == 1)
                    do i = 1, n1
                        ix = merge(1, i, size(x, 1) == 1)
                        terms(i, j, k) = x(ix, jx, kx)
                        weights(i, j, k) = b(kb)
                    end do
                end do
            end do
        rank (2)
            do k = 1, n3
                kx = merge(1, k, size(x, 3) == 1)
                kb = merge(1, k, size(b, 2) == 1)
                do j = 1, n2
                    jx = merge(1, j, size(x, 2) == 1)
                    jb = merge(1, j, size(b, 1) == 1)
                    do i = 1, n1
                        ix = merge(1, i, size(x, 1) == 1)
                        terms(i, j, k) = x(ix, jx, kx)
                        weights(i, j, k) = b(jb, kb)
                    end do
                end do
            end do
        rank (3)
            do k = 1, n3
                kx = merge(1, k, size(x, 3) == 1)
                kb = merge(1, k, size(b, 3) == 1)
                do j = 1, n2
                    jx = merge(1, j, size(x, 2) == 1)
                    jb = merge(1, j, size(b, 2) == 1)
                    do i = 1, n1
                        ix = merge(1, i, size(x, 1) == 1)
                        ib = merge(1, i, size(b, 1) == 1)
                        terms(i, j, k) = x(ix, jx, kx)
                        weights(i, j, k) = b(ib, jb, kb)
                    end do
                end do
            end do
        end select
    end subroutine broadcast_real_terms_3d

    pure function reduction_dims_mode_3d(dims) result(mode)
        integer, intent(in) :: dims(:) !! Fortran dimensions to reduce
        integer :: mode

        integer :: i
        logical :: seen(3)

        mode = 0
        seen = .false.
        if (size(dims) > 3) then
            mode = -1
            return
        end if
        do i = 1, size(dims)
            if (dims(i) < 1 .or. dims(i) > 3) then
                mode = -1
                return
            end if
            if (seen(dims(i))) then
                mode = -1
                return
            end if
            seen(dims(i)) = .true.
            mode = mode + 2 ** (dims(i) - 1)
        end do
    end function reduction_dims_mode_3d

    pure subroutine allocate_keepdims_real_3d(values, n1, n2, n3, mode)
        real(dp), allocatable, intent(out) :: values(:, :, :) !! allocated rank-3 keepdims result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: n3 !! unreduced third-dimension extent
        integer, intent(in) :: mode !! reduction bit mask from reduction_dims_mode_3d

        integer :: m1
        integer :: m2
        integer :: m3

        if (mode < 0) then
            allocate(values(0, 0, 0))
            return
        end if
        m1 = merge(1, n1, iand(mode, 1) /= 0)
        m2 = merge(1, n2, iand(mode, 2) /= 0)
        m3 = merge(1, n3, iand(mode, 4) /= 0)
        allocate(values(m1, m2, m3))
    end subroutine allocate_keepdims_real_3d

    pure function softmax_3d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(softmax_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k) = softmax_1d(x(:, j, k))
                end do
            end do
        else if (dim == 2) then
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k) = softmax_1d(x(i, :, k))
                end do
            end do
        else if (dim == 3) then
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :) = softmax_1d(x(i, j, :))
                end do
            end do
        else
            values = quiet_nan(0.0_dp)
        end if
    end function softmax_3d

    pure function softmax_3d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        mode = reduction_dims_mode_3d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(x))
                values = 1.0_dp
            elsewhere
                values = quiet_nan(0.0_dp)
            end where
        case (1, 2, 4)
            if (mode == 1) values = softmax_3d(x, 1)
            if (mode == 2) values = softmax_3d(x, 2)
            if (mode == 4) values = softmax_3d(x, 3)
        case (3)
            do k = 1, size(x, 3)
                values(:, :, k) = softmax_2d(x(:, :, k))
            end do
        case (5)
            do j = 1, size(x, 2)
                values(:, j, :) = softmax_2d(x(:, j, :))
            end do
        case (6)
            do i = 1, size(x, 1)
                values(i, :, :) = softmax_2d(x(i, :, :))
            end do
        case (7)
            values = softmax_3d(x)
        case default
            values = quiet_nan(0.0_dp)
        end select
    end function softmax_3d_dims

    pure function log_softmax_3d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(log_softmax_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k) = log_softmax_1d(x(:, j, k))
                end do
            end do
        else if (dim == 2) then
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k) = log_softmax_1d(x(i, :, k))
                end do
            end do
        else if (dim == 3) then
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :) = log_softmax_1d(x(i, j, :))
                end do
            end do
        else
            values = quiet_nan(0.0_dp)
        end if
    end function log_softmax_3d

    pure function log_softmax_3d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        mode = reduction_dims_mode_3d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(x))
                values = 0.0_dp
            elsewhere
                values = quiet_nan(0.0_dp)
            end where
        case (1, 2, 4)
            if (mode == 1) values = log_softmax_3d(x, 1)
            if (mode == 2) values = log_softmax_3d(x, 2)
            if (mode == 4) values = log_softmax_3d(x, 3)
        case (3)
            do k = 1, size(x, 3)
                values(:, :, k) = log_softmax_2d(x(:, :, k))
            end do
        case (5)
            do j = 1, size(x, 2)
                values(:, j, :) = log_softmax_2d(x(:, j, :))
            end do
        case (6)
            do i = 1, size(x, 1)
                values(i, :, :) = log_softmax_2d(x(i, :, :))
            end do
        case (7)
            values = log_softmax_3d(x)
        case default
            values = quiet_nan(0.0_dp)
        end select
    end function log_softmax_3d_dims


    pure function logsumexp_real_3d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :) !! rank-3 real values reduced over all dimensions
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        value = logsumexp_complex_3d_all(cmplx(x, 0.0_dp, kind=dp), b)
    end function logsumexp_real_3d_all_complex_b

    pure function logsumexp_real_3d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional real values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :)

        values = logsumexp_complex_3d_dim(cmplx(x, 0.0_dp, kind=dp), dim, b)
    end function logsumexp_real_3d_dim_complex_b

    pure subroutine logsumexp_signed_real_3d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :) !! rank-3 real values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_3d_all(cmplx(x, 0.0_dp, kind=dp), logabs, phase, b)
    end subroutine logsumexp_signed_real_3d_all_complex_b

    pure subroutine logsumexp_signed_real_3d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional real values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        real(dp), allocatable, intent(out) :: logabs(:, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_3d_dim(cmplx(x, 0.0_dp, kind=dp), dim, logabs, phase, b)
    end subroutine logsumexp_signed_real_3d_dim_complex_b

    pure function logsumexp_keepdims_real_3d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional real values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :)

        values = logsumexp_keepdims_complex_3d(cmplx(x, 0.0_dp, kind=dp), dims, b)
    end function logsumexp_keepdims_real_3d_complex_b

    pure subroutine logsumexp_signed_keepdims_real_3d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :) !! three-dimensional real values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! log-magnitudes with kept dimensions
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! phases with kept dimensions
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_keepdims_complex_3d(cmplx(x, 0.0_dp, kind=dp), dims, &
            logabs, phase, b)
    end subroutine logsumexp_signed_keepdims_real_3d_complex_b

    pure function logsumexp_complex_3d_all(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :) !! complex values reduced over all dimensions
        complex(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable complex weights
        complex(dp) :: value

        complex(dp), allocatable :: terms(:, :, :)
        complex(dp), allocatable :: weights(:, :, :)
        complex(dp) :: phase
        logical :: ok
        real(dp) :: logabs

        if (present(b)) then
            call broadcast_complex_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                value = complex_nan()
                return
            end if
            call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        end if
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_3d_all

    pure function logsumexp_complex_3d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :) !! complex values reduced over all dimensions
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp) :: value

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: complex_b3(:, :, :)

        select rank (b)
        rank (0)
            value = logsumexp_complex_3d_all(x, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_3d_all(x, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_3d_all(x, complex_b2)
        rank (3)
            allocate(complex_b3(size(b, 1), size(b, 2), size(b, 3)))
            complex_b3 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_3d_all(x, complex_b3)
        rank default
            value = complex_nan()
        end select
    end function logsumexp_complex_3d_all_real_b

    pure function logsumexp_complex_3d_dim(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        complex(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable complex weights
        complex(dp), allocatable :: values(:, :)

        complex(dp), allocatable :: terms(:, :, :)
        complex(dp), allocatable :: weights(:, :, :)
        complex(dp) :: phase
        integer :: i
        integer :: j
        integer :: k
        logical :: ok
        real(dp) :: logabs

        if (dim < 1 .or. dim > 3) then
            allocate(values(0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                select case (dim)
                case (1)
                    allocate(values(size(x, 2), size(x, 3)))
                case (2)
                    allocate(values(size(x, 1), size(x, 3)))
                case (3)
                    allocate(values(size(x, 1), size(x, 2)))
                end select
                values = complex_nan()
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2), size(x, 3)))
            terms = x
        end if

        select case (dim)
        case (1)
            allocate(values(size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    if (present(b)) then
                        call logsumexp_complex_core(terms(:, j, k), logabs, phase, weights(:, j, k))
                    else
                        call logsumexp_complex_core(terms(:, j, k), logabs, phase)
                    end if
                    values(j, k) = complex_log_from_phase(logabs, phase)
                end do
            end do
        case (2)
            allocate(values(size(terms, 1), size(terms, 3)))
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_complex_core(terms(i, :, k), logabs, phase, weights(i, :, k))
                    else
                        call logsumexp_complex_core(terms(i, :, k), logabs, phase)
                    end if
                    values(i, k) = complex_log_from_phase(logabs, phase)
                end do
            end do
        case (3)
            allocate(values(size(terms, 1), size(terms, 2)))
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_complex_core(terms(i, j, :), logabs, phase, weights(i, j, :))
                    else
                        call logsumexp_complex_core(terms(i, j, :), logabs, phase)
                    end if
                    values(i, j) = complex_log_from_phase(logabs, phase)
                end do
            end do
        end select
    end function logsumexp_complex_3d_dim

    pure function logsumexp_complex_3d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: values(:, :)

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: complex_b3(:, :, :)

        if (dim < 1 .or. dim > 3) then
            allocate(values(0, 0))
            return
        end if
        select rank (b)
        rank (0)
            values = logsumexp_complex_3d_dim(x, dim, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_3d_dim(x, dim, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_3d_dim(x, dim, complex_b2)
        rank (3)
            allocate(complex_b3(size(b, 1), size(b, 2), size(b, 3)))
            complex_b3 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_3d_dim(x, dim, complex_b3)
        rank default
            select case (dim)
            case (1)
                allocate(values(size(x, 2), size(x, 3)))
            case (2)
                allocate(values(size(x, 1), size(x, 3)))
            case (3)
                allocate(values(size(x, 1), size(x, 2)))
            end select
            values = complex_nan()
        end select
    end function logsumexp_complex_3d_dim_real_b

    pure subroutine logsumexp_signed_complex_3d_all(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:, :, :)
        complex(dp), allocatable :: weights(:, :, :)
        logical :: ok

        if (present(b)) then
            call broadcast_complex_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
            call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
                reshape(weights, [size(weights)]))
        else
            call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        end if
    end subroutine logsumexp_signed_complex_3d_all

    pure subroutine logsumexp_signed_complex_3d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! complex values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: complex_b3(:, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_3d_all(x, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_3d_all(x, logabs, phase, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_3d_all(x, logabs, phase, complex_b2)
        rank (3)
            allocate(complex_b3(size(b, 1), size(b, 2), size(b, 3)))
            complex_b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_3d_all(x, logabs, phase, complex_b3)
        rank default
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_3d_all_real_b

    pure subroutine logsumexp_signed_complex_3d_dim(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        real(dp), allocatable, intent(out) :: logabs(:, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :) !! reduced unit phases
        complex(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:, :, :)
        complex(dp), allocatable :: weights(:, :, :)
        integer :: i
        integer :: j
        integer :: k
        logical :: ok

        if (dim < 1 .or. dim > 3) then
            allocate(logabs(0, 0), phase(0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                select case (dim)
                case (1)
                    allocate(logabs(size(x, 2), size(x, 3)), phase(size(x, 2), size(x, 3)))
                case (2)
                    allocate(logabs(size(x, 1), size(x, 3)), phase(size(x, 1), size(x, 3)))
                case (3)
                    allocate(logabs(size(x, 1), size(x, 2)), phase(size(x, 1), size(x, 2)))
                end select
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
                return
            end if
        else
            allocate(terms(size(x, 1), size(x, 2), size(x, 3)))
            terms = x
        end if

        select case (dim)
        case (1)
            allocate(logabs(size(terms, 2), size(terms, 3)), phase(size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    if (present(b)) then
                        call logsumexp_complex_core(terms(:, j, k), logabs(j, k), phase(j, k), &
                            weights(:, j, k))
                    else
                        call logsumexp_complex_core(terms(:, j, k), logabs(j, k), phase(j, k))
                    end if
                end do
            end do
        case (2)
            allocate(logabs(size(terms, 1), size(terms, 3)), phase(size(terms, 1), size(terms, 3)))
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_complex_core(terms(i, :, k), logabs(i, k), phase(i, k), &
                            weights(i, :, k))
                    else
                        call logsumexp_complex_core(terms(i, :, k), logabs(i, k), phase(i, k))
                    end if
                end do
            end do
        case (3)
            allocate(logabs(size(terms, 1), size(terms, 2)), phase(size(terms, 1), size(terms, 2)))
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    if (present(b)) then
                        call logsumexp_complex_core(terms(i, j, :), logabs(i, j), phase(i, j), &
                            weights(i, j, :))
                    else
                        call logsumexp_complex_core(terms(i, j, :), logabs(i, j), phase(i, j))
                    end if
                end do
            end do
        end select
    end subroutine logsumexp_signed_complex_3d_dim

    pure subroutine logsumexp_signed_complex_3d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 3
        real(dp), allocatable, intent(out) :: logabs(:, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: complex_b3(:, :, :)

        if (dim < 1 .or. dim > 3) then
            allocate(logabs(0, 0), phase(0, 0))
            return
        end if
        select rank (b)
        rank (0)
            call logsumexp_signed_complex_3d_dim(x, dim, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_3d_dim(x, dim, logabs, phase, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_3d_dim(x, dim, logabs, phase, complex_b2)
        rank (3)
            allocate(complex_b3(size(b, 1), size(b, 2), size(b, 3)))
            complex_b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_3d_dim(x, dim, logabs, phase, complex_b3)
        rank default
            select case (dim)
            case (1)
                allocate(logabs(size(x, 2), size(x, 3)), phase(size(x, 2), size(x, 3)))
            case (2)
                allocate(logabs(size(x, 1), size(x, 3)), phase(size(x, 1), size(x, 3)))
            case (3)
                allocate(logabs(size(x, 1), size(x, 2)), phase(size(x, 1), size(x, 2)))
            end select
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_3d_dim_real_b

    pure function logsumexp_keepdims_complex_3d(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable complex weights
        complex(dp), allocatable :: values(:, :, :)

        complex(dp), allocatable :: terms(:, :, :)
        complex(dp), allocatable :: weights(:, :, :)
        logical :: ok

        if (reduction_dims_mode_3d(dims) < 0) then
            allocate(values(0, 0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_complex_3d(values, size(x, 1), size(x, 2), size(x, 3), &
                    reduction_dims_mode_3d(dims))
                if (size(values) > 0) values = complex_nan()
                return
            end if
            values = logsumexp_keepdims_complex_3d_full(terms, dims, weights)
        else
            values = logsumexp_keepdims_complex_3d_full(x, dims)
        end if
    end function logsumexp_keepdims_complex_3d

    pure function logsumexp_keepdims_complex_3d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: values(:, :, :)

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: complex_b3(:, :, :)

        select rank (b)
        rank (0)
            values = logsumexp_keepdims_complex_3d(x, dims, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_3d(x, dims, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_3d(x, dims, complex_b2)
        rank (3)
            allocate(complex_b3(size(b, 1), size(b, 2), size(b, 3)))
            complex_b3 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_3d(x, dims, complex_b3)
        rank default
            call allocate_keepdims_complex_3d(values, size(x, 1), size(x, 2), size(x, 3), &
                reduction_dims_mode_3d(dims))
            if (size(values) > 0) values = complex_nan()
        end select
    end function logsumexp_keepdims_complex_3d_real_b

    pure function logsumexp_keepdims_complex_3d_full(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in), optional :: b(:, :, :) !! optional fully broadcast complex weights
        complex(dp), allocatable :: values(:, :, :)

        complex(dp) :: phase
        integer :: i
        integer :: j
        integer :: k
        integer :: mode
        real(dp) :: logabs

        mode = reduction_dims_mode_3d(dims)
        if (mode < 0) then
            allocate(values(0, 0, 0))
            return
        end if
        call allocate_keepdims_complex_3d(values, size(x, 1), size(x, 2), size(x, 3), mode)

        select case (mode)
        case (0)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        if (present(b)) then
                            call logsumexp_complex_core([x(i, j, k)], logabs, phase, [b(i, j, k)])
                        else
                            call logsumexp_complex_core([x(i, j, k)], logabs, phase)
                        end if
                        values(i, j, k) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (1, 2, 4)
            if (mode == 1) then
                if (present(b)) then
                    values(1, :, :) = logsumexp_complex_3d_dim(x, 1, b)
                else
                    values(1, :, :) = logsumexp_complex_3d_dim(x, 1)
                end if
            else if (mode == 2) then
                if (present(b)) then
                    values(:, 1, :) = logsumexp_complex_3d_dim(x, 2, b)
                else
                    values(:, 1, :) = logsumexp_complex_3d_dim(x, 2)
                end if
            else
                if (present(b)) then
                    values(:, :, 1) = logsumexp_complex_3d_dim(x, 3, b)
                else
                    values(:, :, 1) = logsumexp_complex_3d_dim(x, 3)
                end if
            end if
        case (3)
            do k = 1, size(x, 3)
                if (present(b)) then
                    values(1, 1, k) = logsumexp_complex_2d_all(x(:, :, k), b(:, :, k))
                else
                    values(1, 1, k) = logsumexp_complex_2d_all(x(:, :, k))
                end if
            end do
        case (5)
            do j = 1, size(x, 2)
                if (present(b)) then
                    values(1, j, 1) = logsumexp_complex_2d_all(x(:, j, :), b(:, j, :))
                else
                    values(1, j, 1) = logsumexp_complex_2d_all(x(:, j, :))
                end if
            end do
        case (6)
            do i = 1, size(x, 1)
                if (present(b)) then
                    values(i, 1, 1) = logsumexp_complex_2d_all(x(i, :, :), b(i, :, :))
                else
                    values(i, 1, 1) = logsumexp_complex_2d_all(x(i, :, :))
                end if
            end do
        case (7)
            if (present(b)) then
                values(1, 1, 1) = logsumexp_complex_3d_all(x, b)
            else
                values(1, 1, 1) = logsumexp_complex_3d_all(x)
            end if
        end select
    end function logsumexp_keepdims_complex_3d_full

    pure subroutine logsumexp_signed_keepdims_complex_3d(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! log-magnitudes with kept dimensions
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! unit phases with kept dimensions
        complex(dp), intent(in), optional :: &
            b(..) !! optional scalar or broadcastable complex weights

        complex(dp), allocatable :: terms(:, :, :)
        complex(dp), allocatable :: weights(:, :, :)
        logical :: ok

        if (reduction_dims_mode_3d(dims) < 0) then
            allocate(logabs(0, 0, 0), phase(0, 0, 0))
            return
        end if
        if (present(b)) then
            call broadcast_complex_terms_3d(x, b, terms, weights, ok)
            if (.not. ok) then
                call allocate_keepdims_real_3d(logabs, size(x, 1), size(x, 2), size(x, 3), &
                    reduction_dims_mode_3d(dims))
                call allocate_keepdims_complex_3d(phase, size(x, 1), size(x, 2), size(x, 3), &
                    reduction_dims_mode_3d(dims))
                if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
                if (size(phase) > 0) phase = complex_nan()
                return
            end if
            call logsumexp_signed_keepdims_complex_3d_full(terms, dims, logabs, phase, weights)
        else
            call logsumexp_signed_keepdims_complex_3d_full(x, dims, logabs, phase)
        end if
    end subroutine logsumexp_signed_keepdims_complex_3d

    pure subroutine logsumexp_signed_keepdims_complex_3d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! log-magnitudes with kept dimensions
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! unit phases with kept dimensions
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: complex_b1(:)
        complex(dp), allocatable :: complex_b2(:, :)
        complex(dp), allocatable :: complex_b3(:, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_keepdims_complex_3d(x, dims, logabs, phase, &
                cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(complex_b1(size(b)))
            complex_b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_3d(x, dims, logabs, phase, complex_b1)
        rank (2)
            allocate(complex_b2(size(b, 1), size(b, 2)))
            complex_b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_3d(x, dims, logabs, phase, complex_b2)
        rank (3)
            allocate(complex_b3(size(b, 1), size(b, 2), size(b, 3)))
            complex_b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_3d(x, dims, logabs, phase, complex_b3)
        rank default
            call allocate_keepdims_real_3d(logabs, size(x, 1), size(x, 2), size(x, 3), &
                reduction_dims_mode_3d(dims))
            call allocate_keepdims_complex_3d(phase, size(x, 1), size(x, 2), size(x, 3), &
                reduction_dims_mode_3d(dims))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
        end select
    end subroutine logsumexp_signed_keepdims_complex_3d_real_b

    pure subroutine logsumexp_signed_keepdims_complex_3d_full(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! log-magnitudes with kept dimensions
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! unit phases with kept dimensions
        complex(dp), intent(in), optional :: b(:, :, :) !! optional fully broadcast complex weights

        complex(dp), allocatable :: reduced_phase(:, :)
        integer :: i
        integer :: j
        integer :: k
        integer :: mode
        real(dp), allocatable :: reduced_logabs(:, :)

        mode = reduction_dims_mode_3d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0), phase(0, 0, 0))
            return
        end if
        call allocate_keepdims_real_3d(logabs, size(x, 1), size(x, 2), size(x, 3), mode)
        call allocate_keepdims_complex_3d(phase, size(x, 1), size(x, 2), size(x, 3), mode)

        select case (mode)
        case (0)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        if (present(b)) then
                            call logsumexp_complex_core([x(i, j, k)], logabs(i, j, k), &
                                phase(i, j, k), [b(i, j, k)])
                        else
                            call logsumexp_complex_core([x(i, j, k)], logabs(i, j, k), &
                                phase(i, j, k))
                        end if
                    end do
                end do
            end do
        case (1, 2, 4)
            if (mode == 1) then
                if (present(b)) then
                    call logsumexp_signed_complex_3d_dim(x, 1, reduced_logabs, reduced_phase, b)
                else
                    call logsumexp_signed_complex_3d_dim(x, 1, reduced_logabs, reduced_phase)
                end if
                logabs(1, :, :) = reduced_logabs
                phase(1, :, :) = reduced_phase
            else if (mode == 2) then
                if (present(b)) then
                    call logsumexp_signed_complex_3d_dim(x, 2, reduced_logabs, reduced_phase, b)
                else
                    call logsumexp_signed_complex_3d_dim(x, 2, reduced_logabs, reduced_phase)
                end if
                logabs(:, 1, :) = reduced_logabs
                phase(:, 1, :) = reduced_phase
            else
                if (present(b)) then
                    call logsumexp_signed_complex_3d_dim(x, 3, reduced_logabs, reduced_phase, b)
                else
                    call logsumexp_signed_complex_3d_dim(x, 3, reduced_logabs, reduced_phase)
                end if
                logabs(:, :, 1) = reduced_logabs
                phase(:, :, 1) = reduced_phase
            end if
        case (3)
            do k = 1, size(x, 3)
                if (present(b)) then
                    call logsumexp_signed_complex_2d_all(x(:, :, k), logabs(1, 1, k), &
                        phase(1, 1, k), b(:, :, k))
                else
                    call logsumexp_signed_complex_2d_all(x(:, :, k), logabs(1, 1, k), &
                        phase(1, 1, k))
                end if
            end do
        case (5)
            do j = 1, size(x, 2)
                if (present(b)) then
                    call logsumexp_signed_complex_2d_all(x(:, j, :), logabs(1, j, 1), &
                        phase(1, j, 1), b(:, j, :))
                else
                    call logsumexp_signed_complex_2d_all(x(:, j, :), logabs(1, j, 1), &
                        phase(1, j, 1))
                end if
            end do
        case (6)
            do i = 1, size(x, 1)
                if (present(b)) then
                    call logsumexp_signed_complex_2d_all(x(i, :, :), logabs(i, 1, 1), &
                        phase(i, 1, 1), b(i, :, :))
                else
                    call logsumexp_signed_complex_2d_all(x(i, :, :), logabs(i, 1, 1), &
                        phase(i, 1, 1))
                end if
            end do
        case (7)
            if (present(b)) then
                call logsumexp_signed_complex_3d_all(x, logabs(1, 1, 1), phase(1, 1, 1), b)
            else
                call logsumexp_signed_complex_3d_all(x, logabs(1, 1, 1), phase(1, 1, 1))
            end if
        end select
    end subroutine logsumexp_signed_keepdims_complex_3d_full

    pure subroutine broadcast_complex_terms_3d(x, b, terms, weights, ok)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional logarithmic terms to broadcast
        complex(dp), intent(in) :: b(..) !! scalar or trailing-array weights broadcast jointly with x
        complex(dp), allocatable, intent(out) :: terms(:, :, :) !! input terms broadcast to the joint shape
        complex(dp), allocatable, intent(out) :: weights(:, :, :) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: b3
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: k
        integer :: kb
        integer :: kx
        integer :: n1
        integer :: n2
        integer :: n3

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2), size(x, 3)), &
                weights(size(x, 1), size(x, 2), size(x, 3)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = 1
            b3 = size(b)
        rank (2)
            b1 = 1
            b2 = size(b, 1)
            b3 = size(b, 2)
        rank (3)
            b1 = size(b, 1)
            b2 = size(b, 2)
            b3 = size(b, 3)
        rank default
            ok = .false.
            allocate(terms(0, 0, 0), weights(0, 0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1) .or. &
                (size(x, 3) /= b3 .and. size(x, 3) /= 1 .and. b3 /= 1)) then
            ok = .false.
            allocate(terms(0, 0, 0), weights(0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        if (n3 == 1) n3 = b3
        allocate(terms(n1, n2, n3), weights(n1, n2, n3))
        select rank (b)
        rank (1)
            do k = 1, n3
                kx = merge(1, k, size(x, 3) == 1)
                kb = merge(1, k, size(b) == 1)
                do j = 1, n2
                    jx = merge(1, j, size(x, 2) == 1)
                    do i = 1, n1
                        ix = merge(1, i, size(x, 1) == 1)
                        terms(i, j, k) = x(ix, jx, kx)
                        weights(i, j, k) = b(kb)
                    end do
                end do
            end do
        rank (2)
            do k = 1, n3
                kx = merge(1, k, size(x, 3) == 1)
                kb = merge(1, k, size(b, 2) == 1)
                do j = 1, n2
                    jx = merge(1, j, size(x, 2) == 1)
                    jb = merge(1, j, size(b, 1) == 1)
                    do i = 1, n1
                        ix = merge(1, i, size(x, 1) == 1)
                        terms(i, j, k) = x(ix, jx, kx)
                        weights(i, j, k) = b(jb, kb)
                    end do
                end do
            end do
        rank (3)
            do k = 1, n3
                kx = merge(1, k, size(x, 3) == 1)
                kb = merge(1, k, size(b, 3) == 1)
                do j = 1, n2
                    jx = merge(1, j, size(x, 2) == 1)
                    jb = merge(1, j, size(b, 2) == 1)
                    do i = 1, n1
                        ix = merge(1, i, size(x, 1) == 1)
                        ib = merge(1, i, size(b, 1) == 1)
                        terms(i, j, k) = x(ix, jx, kx)
                        weights(i, j, k) = b(ib, jb, kb)
                    end do
                end do
            end do
        end select
    end subroutine broadcast_complex_terms_3d

    pure subroutine allocate_keepdims_complex_3d(values, n1, n2, n3, mode)
        complex(dp), allocatable, intent(out) :: values(:, :, :) !! allocated rank-3 keepdims result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: n3 !! unreduced third-dimension extent
        integer, intent(in) :: mode !! reduction bit mask from reduction_dims_mode_3d

        integer :: m1
        integer :: m2
        integer :: m3

        if (mode < 0) then
            allocate(values(0, 0, 0))
            return
        end if
        m1 = merge(1, n1, iand(mode, 1) /= 0)
        m2 = merge(1, n2, iand(mode, 2) /= 0)
        m3 = merge(1, n3, iand(mode, 4) /= 0)
        allocate(values(m1, m2, m3))
    end subroutine allocate_keepdims_complex_3d

    pure function softmax_complex_3d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        complex(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(softmax_complex_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k) = softmax_complex_1d(x(:, j, k))
                end do
            end do
        else if (dim == 2) then
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k) = softmax_complex_1d(x(i, :, k))
                end do
            end do
        else if (dim == 3) then
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :) = softmax_complex_1d(x(i, j, :))
                end do
            end do
        else
            values = complex_nan()
        end if
    end function softmax_complex_3d

    pure function softmax_complex_3d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        mode = reduction_dims_mode_3d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                values = cmplx(1.0_dp, 0.0_dp, kind=dp)
            elsewhere
                values = complex_nan()
            end where
        case (1, 2, 4)
            if (mode == 1) values = softmax_complex_3d(x, 1)
            if (mode == 2) values = softmax_complex_3d(x, 2)
            if (mode == 4) values = softmax_complex_3d(x, 3)
        case (3)
            do k = 1, size(x, 3)
                values(:, :, k) = softmax_complex_2d(x(:, :, k))
            end do
        case (5)
            do j = 1, size(x, 2)
                values(:, j, :) = softmax_complex_2d(x(:, j, :))
            end do
        case (6)
            do i = 1, size(x, 1)
                values(i, :, :) = softmax_complex_2d(x(i, :, :))
            end do
        case (7)
            values = softmax_complex_3d(x)
        case default
            values = complex_nan()
        end select
    end function softmax_complex_3d_dims

    pure function log_softmax_complex_3d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to log-normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        complex(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(log_softmax_complex_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k) = log_softmax_complex_1d(x(:, j, k))
                end do
            end do
        else if (dim == 2) then
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k) = log_softmax_complex_1d(x(i, :, k))
                end do
            end do
        else if (dim == 3) then
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :) = log_softmax_complex_1d(x(i, j, :))
                end do
            end do
        else
            values = complex_nan()
        end if
    end function log_softmax_complex_3d

    pure function log_softmax_complex_3d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :) !! three-dimensional complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3)))
        mode = reduction_dims_mode_3d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                values = cmplx(0.0_dp, 0.0_dp, kind=dp)
            elsewhere
                values = complex_nan()
            end where
        case (1, 2, 4)
            if (mode == 1) values = log_softmax_complex_3d(x, 1)
            if (mode == 2) values = log_softmax_complex_3d(x, 2)
            if (mode == 4) values = log_softmax_complex_3d(x, 3)
        case (3)
            do k = 1, size(x, 3)
                values(:, :, k) = log_softmax_complex_2d(x(:, :, k))
            end do
        case (5)
            do j = 1, size(x, 2)
                values(:, j, :) = log_softmax_complex_2d(x(:, j, :))
            end do
        case (6)
            do i = 1, size(x, 1)
                values(i, :, :) = log_softmax_complex_2d(x(i, :, :))
            end do
        case (7)
            values = log_softmax_complex_3d(x)
        case default
            values = complex_nan()
        end select
    end function log_softmax_complex_3d_dims

    pure function logsumexp_4d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_core(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_4d_all

    pure function logsumexp_4d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        real(dp) :: logabs
        real(dp) :: sign_value

        if (dim < 1 .or. dim > 4) then
            allocate(values(0, 0, 0))
            return
        end if

        select case (dim)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        call logsumexp_core(x(:, j, k, l), logabs, sign_value)
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(j, k, l) = quiet_nan(0.0_dp)
                        else
                            values(j, k, l) = logabs
                        end if
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        call logsumexp_core(x(i, :, k, l), logabs, sign_value)
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, k, l) = quiet_nan(0.0_dp)
                        else
                            values(i, k, l) = logabs
                        end if
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_core(x(i, j, :, l), logabs, sign_value)
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, j, l) = quiet_nan(0.0_dp)
                        else
                            values(i, j, l) = logabs
                        end if
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_core(x(i, j, k, :), logabs, sign_value)
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, j, k) = quiet_nan(0.0_dp)
                        else
                            values(i, j, k) = logabs
                        end if
                    end do
                end do
            end do
        end select
    end function logsumexp_4d_dim

    pure subroutine logsumexp_signed_4d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum: zero or +1

        call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_4d_all

    pure subroutine logsumexp_signed_4d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :) !! signs along reduced dimension

        integer :: i
        integer :: j
        integer :: k
        integer :: l

        if (dim < 1 .or. dim > 4) then
            allocate(logabs(0, 0, 0), sign_value(0, 0, 0))
            return
        end if

        select case (dim)
        case (1)
            allocate(logabs(size(x, 2), size(x, 3), size(x, 4)))
            allocate(sign_value(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        call logsumexp_core(x(:, j, k, l), logabs(j, k, l), &
                            sign_value(j, k, l))
                    end do
                end do
            end do
        case (2)
            allocate(logabs(size(x, 1), size(x, 3), size(x, 4)))
            allocate(sign_value(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        call logsumexp_core(x(i, :, k, l), logabs(i, k, l), &
                            sign_value(i, k, l))
                    end do
                end do
            end do
        case (3)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 4)))
            allocate(sign_value(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_core(x(i, j, :, l), logabs(i, j, l), &
                            sign_value(i, j, l))
                    end do
                end do
            end do
        case (4)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3)))
            allocate(sign_value(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_core(x(i, j, k, :), logabs(i, j, k), &
                            sign_value(i, j, k))
                    end do
                end do
            end do
        end select
    end subroutine logsumexp_signed_4d_dim

    pure function reduction_dims_mode_4d(dims) result(mode)
        integer, intent(in) :: dims(:) !! Fortran dimensions to reduce
        integer :: mode

        integer :: i
        logical :: seen(4)

        mode = 0
        seen = .false.
        if (size(dims) > 4) then
            mode = -1
            return
        end if
        do i = 1, size(dims)
            if (dims(i) < 1 .or. dims(i) > 4) then
                mode = -1
                return
            end if
            if (seen(dims(i))) then
                mode = -1
                return
            end if
            seen(dims(i)) = .true.
            mode = mode + 2 ** (dims(i) - 1)
        end do
    end function reduction_dims_mode_4d

    pure subroutine allocate_keepdims_real_4d(values, n1, n2, n3, n4, mode)
        real(dp), allocatable, intent(out) :: values(:, :, :, :) !! allocated rank-4 result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: n3 !! unreduced third-dimension extent
        integer, intent(in) :: n4 !! unreduced fourth-dimension extent
        integer, intent(in) :: mode !! reduction bit mask from reduction_dims_mode_4d

        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4

        if (mode < 0) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        m1 = merge(1, n1, iand(mode, 1) /= 0)
        m2 = merge(1, n2, iand(mode, 2) /= 0)
        m3 = merge(1, n3, iand(mode, 4) /= 0)
        m4 = merge(1, n4, iand(mode, 8) /= 0)
        allocate(values(m1, m2, m3, m4))
    end subroutine allocate_keepdims_real_4d

    pure function logsumexp_keepdims_4d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode

        mode = reduction_dims_mode_4d(dims)
        call allocate_keepdims_real_4d(values, size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), mode)
        if (mode < 0) return

        select case (mode)
        case (0)
            values = x
        case (1, 2, 4, 8)
            if (mode == 1) values = reshape(logsumexp_4d_dim(x, 1), shape(values))
            if (mode == 2) values = reshape(logsumexp_4d_dim(x, 2), shape(values))
            if (mode == 4) values = reshape(logsumexp_4d_dim(x, 3), shape(values))
            if (mode == 8) values = reshape(logsumexp_4d_dim(x, 4), shape(values))
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    values(1, 1, k, l) = logsumexp_2d_all(x(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    values(1, j, 1, l) = logsumexp_2d_all(x(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    values(i, 1, 1, l) = logsumexp_2d_all(x(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(1, j, k, 1) = logsumexp_2d_all(x(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, 1, k, 1) = logsumexp_2d_all(x(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, 1, 1) = logsumexp_2d_all(x(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                values(1, 1, 1, l) = logsumexp_3d_all(x(:, :, :, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                values(1, 1, k, 1) = logsumexp_3d_all(x(:, :, k, :))
            end do
        case (13)
            do j = 1, size(x, 2)
                values(1, j, 1, 1) = logsumexp_3d_all(x(:, j, :, :))
            end do
        case (14)
            do i = 1, size(x, 1)
                values(i, 1, 1, 1) = logsumexp_3d_all(x(i, :, :, :))
            end do
        case (15)
            values(1, 1, 1, 1) = logsumexp_4d_all(x)
        end select
    end function logsumexp_keepdims_4d

    pure subroutine logsumexp_signed_keepdims_4d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! log-absolute sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :) !! signs of reduced sums

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        real(dp), allocatable :: reduced_logabs(:, :, :)
        real(dp), allocatable :: reduced_sign(:, :, :)

        mode = reduction_dims_mode_4d(dims)
        call allocate_keepdims_real_4d(logabs, size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), mode)
        call allocate_keepdims_real_4d(sign_value, size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), mode)
        if (mode < 0) return

        select case (mode)
        case (0)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        do i = 1, size(x, 1)
                            call logsumexp_core([x(i, j, k, l)], logabs(i, j, k, l), &
                                sign_value(i, j, k, l))
                        end do
                    end do
                end do
            end do
        case (1, 2, 4, 8)
            if (mode == 1) call logsumexp_signed_4d_dim(x, 1, reduced_logabs, reduced_sign)
            if (mode == 2) call logsumexp_signed_4d_dim(x, 2, reduced_logabs, reduced_sign)
            if (mode == 4) call logsumexp_signed_4d_dim(x, 3, reduced_logabs, reduced_sign)
            if (mode == 8) call logsumexp_signed_4d_dim(x, 4, reduced_logabs, reduced_sign)
            logabs = reshape(reduced_logabs, shape(logabs))
            sign_value = reshape(reduced_sign, shape(sign_value))
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    call logsumexp_signed_2d_all(x(:, :, k, l), logabs(1, 1, k, l), &
                        sign_value(1, 1, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    call logsumexp_signed_2d_all(x(:, j, :, l), logabs(1, j, 1, l), &
                        sign_value(1, j, 1, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    call logsumexp_signed_2d_all(x(i, :, :, l), logabs(i, 1, 1, l), &
                        sign_value(i, 1, 1, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    call logsumexp_signed_2d_all(x(:, j, k, :), logabs(1, j, k, 1), &
                        sign_value(1, j, k, 1))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    call logsumexp_signed_2d_all(x(i, :, k, :), logabs(i, 1, k, 1), &
                        sign_value(i, 1, k, 1))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    call logsumexp_signed_2d_all(x(i, j, :, :), logabs(i, j, 1, 1), &
                        sign_value(i, j, 1, 1))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                call logsumexp_signed_3d_all(x(:, :, :, l), logabs(1, 1, 1, l), &
                    sign_value(1, 1, 1, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                call logsumexp_signed_3d_all(x(:, :, k, :), logabs(1, 1, k, 1), &
                    sign_value(1, 1, k, 1))
            end do
        case (13)
            do j = 1, size(x, 2)
                call logsumexp_signed_3d_all(x(:, j, :, :), logabs(1, j, 1, 1), &
                    sign_value(1, j, 1, 1))
            end do
        case (14)
            do i = 1, size(x, 1)
                call logsumexp_signed_3d_all(x(i, :, :, :), logabs(i, 1, 1, 1), &
                    sign_value(i, 1, 1, 1))
            end do
        case (15)
            call logsumexp_signed_4d_all(x, logabs(1, 1, 1, 1), &
                sign_value(1, 1, 1, 1))
        end select
    end subroutine logsumexp_signed_keepdims_4d

    pure function softmax_4d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        real(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(softmax_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(:, j, k, l) = softmax_1d(x(:, j, k, l))
                    end do
                end do
            end do
        else if (dim == 2) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, :, k, l) = softmax_1d(x(i, :, k, l))
                    end do
                end do
            end do
        else if (dim == 3) then
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, :, l) = softmax_1d(x(i, j, :, l))
                    end do
                end do
            end do
        else if (dim == 4) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k, :) = softmax_1d(x(i, j, k, :))
                    end do
                end do
            end do
        else
            values = quiet_nan(0.0_dp)
        end if
    end function softmax_4d

    pure function softmax_4d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        mode = reduction_dims_mode_4d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(x))
                values = 1.0_dp
            elsewhere
                values = quiet_nan(0.0_dp)
            end where
        case (1, 2, 4, 8)
            if (mode == 1) values = softmax_4d(x, 1)
            if (mode == 2) values = softmax_4d(x, 2)
            if (mode == 4) values = softmax_4d(x, 3)
            if (mode == 8) values = softmax_4d(x, 4)
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    values(:, :, k, l) = softmax_2d(x(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    values(:, j, :, l) = softmax_2d(x(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    values(i, :, :, l) = softmax_2d(x(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k, :) = softmax_2d(x(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k, :) = softmax_2d(x(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :, :) = softmax_2d(x(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                values(:, :, :, l) = softmax_3d(x(:, :, :, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                values(:, :, k, :) = softmax_3d(x(:, :, k, :))
            end do
        case (13)
            do j = 1, size(x, 2)
                values(:, j, :, :) = softmax_3d(x(:, j, :, :))
            end do
        case (14)
            do i = 1, size(x, 1)
                values(i, :, :, :) = softmax_3d(x(i, :, :, :))
            end do
        case (15)
            values = softmax_4d(x)
        case default
            values = quiet_nan(0.0_dp)
        end select
    end function softmax_4d_dims

    pure function log_softmax_4d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        real(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(log_softmax_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(:, j, k, l) = log_softmax_1d(x(:, j, k, l))
                    end do
                end do
            end do
        else if (dim == 2) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, :, k, l) = log_softmax_1d(x(i, :, k, l))
                    end do
                end do
            end do
        else if (dim == 3) then
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, :, l) = log_softmax_1d(x(i, j, :, l))
                    end do
                end do
            end do
        else if (dim == 4) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k, :) = log_softmax_1d(x(i, j, k, :))
                    end do
                end do
            end do
        else
            values = quiet_nan(0.0_dp)
        end if
    end function log_softmax_4d

    pure function log_softmax_4d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! four-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        mode = reduction_dims_mode_4d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(x))
                values = 0.0_dp
            elsewhere
                values = quiet_nan(0.0_dp)
            end where
        case (1, 2, 4, 8)
            if (mode == 1) values = log_softmax_4d(x, 1)
            if (mode == 2) values = log_softmax_4d(x, 2)
            if (mode == 4) values = log_softmax_4d(x, 3)
            if (mode == 8) values = log_softmax_4d(x, 4)
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    values(:, :, k, l) = log_softmax_2d(x(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    values(:, j, :, l) = log_softmax_2d(x(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    values(i, :, :, l) = log_softmax_2d(x(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k, :) = log_softmax_2d(x(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k, :) = log_softmax_2d(x(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :, :) = log_softmax_2d(x(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                values(:, :, :, l) = log_softmax_3d(x(:, :, :, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                values(:, :, k, :) = log_softmax_3d(x(:, :, k, :))
            end do
        case (13)
            do j = 1, size(x, 2)
                values(:, j, :, :) = log_softmax_3d(x(:, j, :, :))
            end do
        case (14)
            do i = 1, size(x, 1)
                values(i, :, :, :) = log_softmax_3d(x(i, :, :, :))
            end do
        case (15)
            values = log_softmax_4d(x)
        case default
            values = quiet_nan(0.0_dp)
        end select
    end function log_softmax_4d_dims

    pure function logsumexp_complex_4d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 values reduced over all dimensions
        complex(dp) :: value

        complex(dp) :: phase
        real(dp) :: logabs

        call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_4d_all

    pure function logsumexp_complex_4d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        complex(dp), allocatable :: values(:, :, :)

        complex(dp) :: phase
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        real(dp) :: logabs

        if (dim < 1 .or. dim > 4) then
            allocate(values(0, 0, 0))
            return
        end if

        select case (dim)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        call logsumexp_complex_core(x(:, j, k, l), logabs, phase)
                        values(j, k, l) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        call logsumexp_complex_core(x(i, :, k, l), logabs, phase)
                        values(i, k, l) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_complex_core(x(i, j, :, l), logabs, phase)
                        values(i, j, l) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_complex_core(x(i, j, k, :), logabs, phase)
                        values(i, j, k) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        end select
    end function logsumexp_complex_4d_dim

    pure subroutine logsumexp_signed_complex_4d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the exponential sum

        call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_4d_all

    pure subroutine logsumexp_signed_complex_4d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! reduced unit phases

        integer :: i
        integer :: j
        integer :: k
        integer :: l

        if (dim < 1 .or. dim > 4) then
            allocate(logabs(0, 0, 0), phase(0, 0, 0))
            return
        end if

        select case (dim)
        case (1)
            allocate(logabs(size(x, 2), size(x, 3), size(x, 4)))
            allocate(phase(size(x, 2), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        call logsumexp_complex_core(x(:, j, k, l), logabs(j, k, l), &
                            phase(j, k, l))
                    end do
                end do
            end do
        case (2)
            allocate(logabs(size(x, 1), size(x, 3), size(x, 4)))
            allocate(phase(size(x, 1), size(x, 3), size(x, 4)))
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        call logsumexp_complex_core(x(i, :, k, l), logabs(i, k, l), &
                            phase(i, k, l))
                    end do
                end do
            end do
        case (3)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 4)))
            allocate(phase(size(x, 1), size(x, 2), size(x, 4)))
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_complex_core(x(i, j, :, l), logabs(i, j, l), &
                            phase(i, j, l))
                    end do
                end do
            end do
        case (4)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3)))
            allocate(phase(size(x, 1), size(x, 2), size(x, 3)))
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        call logsumexp_complex_core(x(i, j, k, :), logabs(i, j, k), &
                            phase(i, j, k))
                    end do
                end do
            end do
        end select
    end subroutine logsumexp_signed_complex_4d_dim

    pure subroutine allocate_keepdims_complex_4d(values, n1, n2, n3, n4, mode)
        complex(dp), allocatable, intent(out) :: values(:, :, :, :) !! allocated rank-4 result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: n3 !! unreduced third-dimension extent
        integer, intent(in) :: n4 !! unreduced fourth-dimension extent
        integer, intent(in) :: mode !! reduction bit mask from reduction_dims_mode_4d

        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4

        if (mode < 0) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        m1 = merge(1, n1, iand(mode, 1) /= 0)
        m2 = merge(1, n2, iand(mode, 2) /= 0)
        m3 = merge(1, n3, iand(mode, 4) /= 0)
        m4 = merge(1, n4, iand(mode, 8) /= 0)
        allocate(values(m1, m2, m3, m4))
    end subroutine allocate_keepdims_complex_4d

    pure function logsumexp_keepdims_complex_4d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :)

        complex(dp) :: phase
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        real(dp) :: logabs

        mode = reduction_dims_mode_4d(dims)
        call allocate_keepdims_complex_4d(values, size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), mode)
        if (mode < 0) return

        select case (mode)
        case (0)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        do i = 1, size(x, 1)
                            call logsumexp_complex_core([x(i, j, k, l)], logabs, phase)
                            values(i, j, k, l) = complex_log_from_phase(logabs, phase)
                        end do
                    end do
                end do
            end do
        case (1, 2, 4, 8)
            if (mode == 1) values = reshape(logsumexp_complex_4d_dim(x, 1), shape(values))
            if (mode == 2) values = reshape(logsumexp_complex_4d_dim(x, 2), shape(values))
            if (mode == 4) values = reshape(logsumexp_complex_4d_dim(x, 3), shape(values))
            if (mode == 8) values = reshape(logsumexp_complex_4d_dim(x, 4), shape(values))
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    values(1, 1, k, l) = logsumexp_complex_2d_all(x(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    values(1, j, 1, l) = logsumexp_complex_2d_all(x(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    values(i, 1, 1, l) = logsumexp_complex_2d_all(x(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(1, j, k, 1) = logsumexp_complex_2d_all(x(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, 1, k, 1) = logsumexp_complex_2d_all(x(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, 1, 1) = logsumexp_complex_2d_all(x(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                values(1, 1, 1, l) = logsumexp_complex_3d_all(x(:, :, :, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                values(1, 1, k, 1) = logsumexp_complex_3d_all(x(:, :, k, :))
            end do
        case (13)
            do j = 1, size(x, 2)
                values(1, j, 1, 1) = logsumexp_complex_3d_all(x(:, j, :, :))
            end do
        case (14)
            do i = 1, size(x, 1)
                values(i, 1, 1, 1) = logsumexp_complex_3d_all(x(i, :, :, :))
            end do
        case (15)
            values(1, 1, 1, 1) = logsumexp_complex_4d_all(x)
        end select
    end function logsumexp_keepdims_complex_4d

    pure subroutine logsumexp_signed_keepdims_complex_4d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced unit phases

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        real(dp), allocatable :: reduced_logabs(:, :, :)
        complex(dp), allocatable :: reduced_phase(:, :, :)

        mode = reduction_dims_mode_4d(dims)
        call allocate_keepdims_real_4d(logabs, size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), mode)
        call allocate_keepdims_complex_4d(phase, size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), mode)
        if (mode < 0) return

        select case (mode)
        case (0)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        do i = 1, size(x, 1)
                            call logsumexp_complex_core([x(i, j, k, l)], &
                                logabs(i, j, k, l), phase(i, j, k, l))
                        end do
                    end do
                end do
            end do
        case (1, 2, 4, 8)
            if (mode == 1) call logsumexp_signed_complex_4d_dim(x, 1, reduced_logabs, reduced_phase)
            if (mode == 2) call logsumexp_signed_complex_4d_dim(x, 2, reduced_logabs, reduced_phase)
            if (mode == 4) call logsumexp_signed_complex_4d_dim(x, 3, reduced_logabs, reduced_phase)
            if (mode == 8) call logsumexp_signed_complex_4d_dim(x, 4, reduced_logabs, reduced_phase)
            logabs = reshape(reduced_logabs, shape(logabs))
            phase = reshape(reduced_phase, shape(phase))
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    call logsumexp_signed_complex_2d_all(x(:, :, k, l), logabs(1, 1, k, l), &
                        phase(1, 1, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    call logsumexp_signed_complex_2d_all(x(:, j, :, l), logabs(1, j, 1, l), &
                        phase(1, j, 1, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    call logsumexp_signed_complex_2d_all(x(i, :, :, l), logabs(i, 1, 1, l), &
                        phase(i, 1, 1, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    call logsumexp_signed_complex_2d_all(x(:, j, k, :), logabs(1, j, k, 1), &
                        phase(1, j, k, 1))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    call logsumexp_signed_complex_2d_all(x(i, :, k, :), logabs(i, 1, k, 1), &
                        phase(i, 1, k, 1))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    call logsumexp_signed_complex_2d_all(x(i, j, :, :), logabs(i, j, 1, 1), &
                        phase(i, j, 1, 1))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                call logsumexp_signed_complex_3d_all(x(:, :, :, l), logabs(1, 1, 1, l), &
                    phase(1, 1, 1, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                call logsumexp_signed_complex_3d_all(x(:, :, k, :), logabs(1, 1, k, 1), &
                    phase(1, 1, k, 1))
            end do
        case (13)
            do j = 1, size(x, 2)
                call logsumexp_signed_complex_3d_all(x(:, j, :, :), logabs(1, j, 1, 1), &
                    phase(1, j, 1, 1))
            end do
        case (14)
            do i = 1, size(x, 1)
                call logsumexp_signed_complex_3d_all(x(i, :, :, :), logabs(i, 1, 1, 1), &
                    phase(i, 1, 1, 1))
            end do
        case (15)
            call logsumexp_signed_complex_4d_all(x, logabs(1, 1, 1, 1), phase(1, 1, 1, 1))
        end select
    end subroutine logsumexp_signed_keepdims_complex_4d

    pure function softmax_complex_4d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        complex(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(softmax_complex_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(:, j, k, l) = softmax_complex_1d(x(:, j, k, l))
                    end do
                end do
            end do
        else if (dim == 2) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, :, k, l) = softmax_complex_1d(x(i, :, k, l))
                    end do
                end do
            end do
        else if (dim == 3) then
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, :, l) = softmax_complex_1d(x(i, j, :, l))
                    end do
                end do
            end do
        else if (dim == 4) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k, :) = softmax_complex_1d(x(i, j, k, :))
                    end do
                end do
            end do
        else
            values = complex_nan()
        end if
    end function softmax_complex_4d

    pure function softmax_complex_4d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        mode = reduction_dims_mode_4d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                values = cmplx(1.0_dp, 0.0_dp, kind=dp)
            elsewhere
                values = complex_nan()
            end where
        case (1, 2, 4, 8)
            if (mode == 1) values = softmax_complex_4d(x, 1)
            if (mode == 2) values = softmax_complex_4d(x, 2)
            if (mode == 4) values = softmax_complex_4d(x, 3)
            if (mode == 8) values = softmax_complex_4d(x, 4)
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    values(:, :, k, l) = softmax_complex_2d(x(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    values(:, j, :, l) = softmax_complex_2d(x(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    values(i, :, :, l) = softmax_complex_2d(x(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k, :) = softmax_complex_2d(x(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k, :) = softmax_complex_2d(x(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :, :) = softmax_complex_2d(x(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                values(:, :, :, l) = softmax_complex_3d(x(:, :, :, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                values(:, :, k, :) = softmax_complex_3d(x(:, :, k, :))
            end do
        case (13)
            do j = 1, size(x, 2)
                values(:, j, :, :) = softmax_complex_3d(x(:, j, :, :))
            end do
        case (14)
            do i = 1, size(x, 1)
                values(i, :, :, :) = softmax_complex_3d(x(i, :, :, :))
            end do
        case (15)
            values = softmax_complex_4d(x)
        case default
            values = complex_nan()
        end select
    end function softmax_complex_4d_dims

    pure function log_softmax_complex_4d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to log-normalize
        integer, intent(in), optional :: dim !! reduced dimension; absent normalizes all values
        complex(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        if (size(x) == 0) return

        if (.not. present(dim)) then
            values = reshape(log_softmax_complex_1d(reshape(x, [size(x)])), shape(x))
        else if (dim == 1) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do j = 1, size(x, 2)
                        values(:, j, k, l) = log_softmax_complex_1d(x(:, j, k, l))
                    end do
                end do
            end do
        else if (dim == 2) then
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    do i = 1, size(x, 1)
                        values(i, :, k, l) = log_softmax_complex_1d(x(i, :, k, l))
                    end do
                end do
            end do
        else if (dim == 3) then
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, :, l) = log_softmax_complex_1d(x(i, j, :, l))
                    end do
                end do
            end do
        else if (dim == 4) then
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    do i = 1, size(x, 1)
                        values(i, j, k, :) = log_softmax_complex_1d(x(i, j, k, :))
                    end do
                end do
            end do
        else
            values = complex_nan()
        end if
    end function log_softmax_complex_4d

    pure function log_softmax_complex_4d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! four-dimensional complex values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode

        allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        mode = reduction_dims_mode_4d(dims)
        select case (mode)
        case (0)
            where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                values = cmplx(0.0_dp, 0.0_dp, kind=dp)
            elsewhere
                values = complex_nan()
            end where
        case (1, 2, 4, 8)
            if (mode == 1) values = log_softmax_complex_4d(x, 1)
            if (mode == 2) values = log_softmax_complex_4d(x, 2)
            if (mode == 4) values = log_softmax_complex_4d(x, 3)
            if (mode == 8) values = log_softmax_complex_4d(x, 4)
        case (3)
            do l = 1, size(x, 4)
                do k = 1, size(x, 3)
                    values(:, :, k, l) = log_softmax_complex_2d(x(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(x, 4)
                do j = 1, size(x, 2)
                    values(:, j, :, l) = log_softmax_complex_2d(x(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(x, 4)
                do i = 1, size(x, 1)
                    values(i, :, :, l) = log_softmax_complex_2d(x(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(x, 3)
                do j = 1, size(x, 2)
                    values(:, j, k, :) = log_softmax_complex_2d(x(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(x, 3)
                do i = 1, size(x, 1)
                    values(i, :, k, :) = log_softmax_complex_2d(x(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(x, 2)
                do i = 1, size(x, 1)
                    values(i, j, :, :) = log_softmax_complex_2d(x(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(x, 4)
                values(:, :, :, l) = log_softmax_complex_3d(x(:, :, :, l))
            end do
        case (11)
            do k = 1, size(x, 3)
                values(:, :, k, :) = log_softmax_complex_3d(x(:, :, k, :))
            end do
        case (13)
            do j = 1, size(x, 2)
                values(:, j, :, :) = log_softmax_complex_3d(x(:, j, :, :))
            end do
        case (14)
            do i = 1, size(x, 1)
                values(i, :, :, :) = log_softmax_complex_3d(x(i, :, :, :))
            end do
        case (15)
            values = log_softmax_complex_4d(x)
        case default
            values = complex_nan()
        end select
    end function log_softmax_complex_4d_dims

    pure function logsumexp_4d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        real(dp) :: value

        logical :: ok
        real(dp) :: logabs
        real(dp) :: sign_value
        real(dp), allocatable :: terms(:, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :)

        call broadcast_real_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            value = quiet_nan(0.0_dp)
            return
        end if
        call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
            reshape(weights, [size(weights)]))
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
            value = quiet_nan(0.0_dp)
        else
            value = logabs
        end if
    end function logsumexp_4d_all_weighted

    pure function logsumexp_4d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        real(dp), allocatable :: values(:, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        logical :: ok
        real(dp) :: logabs
        real(dp) :: sign_value
        real(dp), allocatable :: terms(:, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :)

        if (dim < 1 .or. dim > 4) then
            allocate(values(0, 0, 0))
            return
        end if
        call broadcast_real_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            select case (dim)
            case (1)
                allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            case (2)
                allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            case (3)
                allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            case (4)
                allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            end select
            values = quiet_nan(0.0_dp)
            return
        end if

        select case (dim)
        case (1)
            allocate(values(size(terms, 2), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        call logsumexp_core(terms(:, j, k, l), logabs, sign_value, &
                            weights(:, j, k, l))
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(j, k, l) = quiet_nan(0.0_dp)
                        else
                            values(j, k, l) = logabs
                        end if
                    end do
                end do
            end do
        case (2)
            allocate(values(size(terms, 1), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do i = 1, size(terms, 1)
                        call logsumexp_core(terms(i, :, k, l), logabs, sign_value, &
                            weights(i, :, k, l))
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, k, l) = quiet_nan(0.0_dp)
                        else
                            values(i, k, l) = logabs
                        end if
                    end do
                end do
            end do
        case (3)
            allocate(values(size(terms, 1), size(terms, 2), size(terms, 4)))
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_core(terms(i, j, :, l), logabs, sign_value, &
                            weights(i, j, :, l))
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, j, l) = quiet_nan(0.0_dp)
                        else
                            values(i, j, l) = logabs
                        end if
                    end do
                end do
            end do
        case (4)
            allocate(values(size(terms, 1), size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_core(terms(i, j, k, :), logabs, sign_value, &
                            weights(i, j, k, :))
                        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                            values(i, j, k) = quiet_nan(0.0_dp)
                        else
                            values(i, j, k) = logabs
                        end if
                    end do
                end do
            end do
        end select
    end function logsumexp_4d_dim_weighted

    pure subroutine logsumexp_signed_4d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        logical :: ok
        real(dp), allocatable :: terms(:, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :)

        call broadcast_real_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
            reshape(weights, [size(weights)]))
    end subroutine logsumexp_signed_4d_all_weighted

    pure subroutine logsumexp_signed_4d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        logical :: ok
        real(dp), allocatable :: terms(:, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :)

        if (dim < 1 .or. dim > 4) then
            allocate(logabs(0, 0, 0), sign_value(0, 0, 0))
            return
        end if
        call broadcast_real_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            select case (dim)
            case (1)
                allocate(logabs(size(x, 2), size(x, 3), size(x, 4)), &
                    sign_value(size(x, 2), size(x, 3), size(x, 4)))
            case (2)
                allocate(logabs(size(x, 1), size(x, 3), size(x, 4)), &
                    sign_value(size(x, 1), size(x, 3), size(x, 4)))
            case (3)
                allocate(logabs(size(x, 1), size(x, 2), size(x, 4)), &
                    sign_value(size(x, 1), size(x, 2), size(x, 4)))
            case (4)
                allocate(logabs(size(x, 1), size(x, 2), size(x, 3)), &
                    sign_value(size(x, 1), size(x, 2), size(x, 3)))
            end select
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if

        select case (dim)
        case (1)
            allocate(logabs(size(terms, 2), size(terms, 3), size(terms, 4)), &
                sign_value(size(terms, 2), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        call logsumexp_core(terms(:, j, k, l), logabs(j, k, l), &
                            sign_value(j, k, l), weights(:, j, k, l))
                    end do
                end do
            end do
        case (2)
            allocate(logabs(size(terms, 1), size(terms, 3), size(terms, 4)), &
                sign_value(size(terms, 1), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do i = 1, size(terms, 1)
                        call logsumexp_core(terms(i, :, k, l), logabs(i, k, l), &
                            sign_value(i, k, l), weights(i, :, k, l))
                    end do
                end do
            end do
        case (3)
            allocate(logabs(size(terms, 1), size(terms, 2), size(terms, 4)), &
                sign_value(size(terms, 1), size(terms, 2), size(terms, 4)))
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_core(terms(i, j, :, l), logabs(i, j, l), &
                            sign_value(i, j, l), weights(i, j, :, l))
                    end do
                end do
            end do
        case (4)
            allocate(logabs(size(terms, 1), size(terms, 2), size(terms, 3)), &
                sign_value(size(terms, 1), size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_core(terms(i, j, k, :), logabs(i, j, k), &
                            sign_value(i, j, k), weights(i, j, k, :))
                    end do
                end do
            end do
        end select
    end subroutine logsumexp_signed_4d_dim_weighted

    pure function logsumexp_keepdims_4d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        logical :: ok
        real(dp) :: logabs
        real(dp) :: sign_value
        real(dp), allocatable :: terms(:, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :)

        mode = reduction_dims_mode_4d(dims)
        if (mode < 0) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        call broadcast_real_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            call allocate_keepdims_real_4d(values, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), mode)
            if (size(values) > 0) values = quiet_nan(0.0_dp)
            return
        end if
        call allocate_keepdims_real_4d(values, size(terms, 1), size(terms, 2), &
            size(terms, 3), size(terms, 4), mode)

        select case (mode)
        case (0)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        do i = 1, size(terms, 1)
                            call logsumexp_core([terms(i, j, k, l)], logabs, sign_value, &
                                [weights(i, j, k, l)])
                            if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) then
                                values(i, j, k, l) = quiet_nan(0.0_dp)
                            else
                                values(i, j, k, l) = logabs
                            end if
                        end do
                    end do
                end do
            end do
        case (1)
            values(1, :, :, :) = logsumexp_4d_dim_weighted(terms, 1, weights)
        case (2)
            values(:, 1, :, :) = logsumexp_4d_dim_weighted(terms, 2, weights)
        case (4)
            values(:, :, 1, :) = logsumexp_4d_dim_weighted(terms, 3, weights)
        case (8)
            values(:, :, :, 1) = logsumexp_4d_dim_weighted(terms, 4, weights)
        case (3)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    values(1, 1, k, l) = logsumexp_2d_all(terms(:, :, k, l), weights(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    values(1, j, 1, l) = logsumexp_2d_all(terms(:, j, :, l), weights(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(terms, 4)
                do i = 1, size(terms, 1)
                    values(i, 1, 1, l) = logsumexp_2d_all(terms(i, :, :, l), weights(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    values(1, j, k, 1) = logsumexp_2d_all(terms(:, j, k, :), weights(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    values(i, 1, k, 1) = logsumexp_2d_all(terms(i, :, k, :), weights(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    values(i, j, 1, 1) = logsumexp_2d_all(terms(i, j, :, :), weights(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(terms, 4)
                values(1, 1, 1, l) = logsumexp_3d_all(terms(:, :, :, l), weights(:, :, :, l))
            end do
        case (11)
            do k = 1, size(terms, 3)
                values(1, 1, k, 1) = logsumexp_3d_all(terms(:, :, k, :), weights(:, :, k, :))
            end do
        case (13)
            do j = 1, size(terms, 2)
                values(1, j, 1, 1) = logsumexp_3d_all(terms(:, j, :, :), weights(:, j, :, :))
            end do
        case (14)
            do i = 1, size(terms, 1)
                values(i, 1, 1, 1) = logsumexp_3d_all(terms(i, :, :, :), weights(i, :, :, :))
            end do
        case (15)
            values(1, 1, 1, 1) = logsumexp_4d_all_weighted(terms, weights)
        end select
    end function logsumexp_keepdims_4d_weighted

    pure subroutine logsumexp_signed_keepdims_4d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        logical :: ok
        real(dp), allocatable :: reduced_logabs(:, :, :)
        real(dp), allocatable :: reduced_sign(:, :, :)
        real(dp), allocatable :: terms(:, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :)

        mode = reduction_dims_mode_4d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0), sign_value(0, 0, 0, 0))
            return
        end if
        call broadcast_real_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            call allocate_keepdims_real_4d(logabs, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), mode)
            call allocate_keepdims_real_4d(sign_value, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), mode)
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(sign_value) > 0) sign_value = quiet_nan(0.0_dp)
            return
        end if
        call allocate_keepdims_real_4d(logabs, size(terms, 1), size(terms, 2), &
            size(terms, 3), size(terms, 4), mode)
        call allocate_keepdims_real_4d(sign_value, size(terms, 1), size(terms, 2), &
            size(terms, 3), size(terms, 4), mode)

        select case (mode)
        case (0)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        do i = 1, size(terms, 1)
                            call logsumexp_core([terms(i, j, k, l)], logabs(i, j, k, l), &
                                sign_value(i, j, k, l), [weights(i, j, k, l)])
                        end do
                    end do
                end do
            end do
        case (1, 2, 4, 8)
            if (mode == 1) call logsumexp_signed_4d_dim_weighted(terms, 1, reduced_logabs, reduced_sign, weights)
            if (mode == 2) call logsumexp_signed_4d_dim_weighted(terms, 2, reduced_logabs, reduced_sign, weights)
            if (mode == 4) call logsumexp_signed_4d_dim_weighted(terms, 3, reduced_logabs, reduced_sign, weights)
            if (mode == 8) call logsumexp_signed_4d_dim_weighted(terms, 4, reduced_logabs, reduced_sign, weights)
            logabs = reshape(reduced_logabs, shape(logabs))
            sign_value = reshape(reduced_sign, shape(sign_value))
        case (3)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    call logsumexp_signed_2d_all(terms(:, :, k, l), logabs(1, 1, k, l), &
                        sign_value(1, 1, k, l), weights(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    call logsumexp_signed_2d_all(terms(:, j, :, l), logabs(1, j, 1, l), &
                        sign_value(1, j, 1, l), weights(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(terms, 4)
                do i = 1, size(terms, 1)
                    call logsumexp_signed_2d_all(terms(i, :, :, l), logabs(i, 1, 1, l), &
                        sign_value(i, 1, 1, l), weights(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    call logsumexp_signed_2d_all(terms(:, j, k, :), logabs(1, j, k, 1), &
                        sign_value(1, j, k, 1), weights(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    call logsumexp_signed_2d_all(terms(i, :, k, :), logabs(i, 1, k, 1), &
                        sign_value(i, 1, k, 1), weights(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    call logsumexp_signed_2d_all(terms(i, j, :, :), logabs(i, j, 1, 1), &
                        sign_value(i, j, 1, 1), weights(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(terms, 4)
                call logsumexp_signed_3d_all(terms(:, :, :, l), logabs(1, 1, 1, l), &
                    sign_value(1, 1, 1, l), weights(:, :, :, l))
            end do
        case (11)
            do k = 1, size(terms, 3)
                call logsumexp_signed_3d_all(terms(:, :, k, :), logabs(1, 1, k, 1), &
                    sign_value(1, 1, k, 1), weights(:, :, k, :))
            end do
        case (13)
            do j = 1, size(terms, 2)
                call logsumexp_signed_3d_all(terms(:, j, :, :), logabs(1, j, 1, 1), &
                    sign_value(1, j, 1, 1), weights(:, j, :, :))
            end do
        case (14)
            do i = 1, size(terms, 1)
                call logsumexp_signed_3d_all(terms(i, :, :, :), logabs(i, 1, 1, 1), &
                    sign_value(i, 1, 1, 1), weights(i, :, :, :))
            end do
        case (15)
            call logsumexp_signed_4d_all_weighted(terms, logabs(1, 1, 1, 1), &
                sign_value(1, 1, 1, 1), weights)
        end select
    end subroutine logsumexp_signed_keepdims_4d_weighted

    pure function logsumexp_real_4d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        value = logsumexp_complex_4d_all_weighted(cmplx(x, 0.0_dp, kind=dp), b)
    end function logsumexp_real_4d_all_complex_b

    pure function logsumexp_real_4d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :)

        values = logsumexp_complex_4d_dim_weighted(cmplx(x, 0.0_dp, kind=dp), dim, b)
    end function logsumexp_real_4d_dim_complex_b

    pure subroutine logsumexp_signed_real_4d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! complex unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_4d_all_weighted(cmplx(x, 0.0_dp, kind=dp), &
            logabs, phase, b)
    end subroutine logsumexp_signed_real_4d_all_complex_b

    pure subroutine logsumexp_signed_real_4d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! reduced phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_4d_dim_weighted(cmplx(x, 0.0_dp, kind=dp), dim, &
            logabs, phase, b)
    end subroutine logsumexp_signed_real_4d_dim_complex_b

    pure function logsumexp_keepdims_real_4d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :, :)

        values = logsumexp_keepdims_complex_4d_weighted(cmplx(x, 0.0_dp, kind=dp), dims, b)
    end function logsumexp_keepdims_real_4d_complex_b

    pure subroutine logsumexp_signed_keepdims_real_4d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_keepdims_complex_4d_weighted(cmplx(x, 0.0_dp, kind=dp), &
            dims, logabs, phase, b)
    end subroutine logsumexp_signed_keepdims_real_4d_complex_b

    pure function logsumexp_complex_4d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        complex(dp) :: phase
        complex(dp), allocatable :: terms(:, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :)
        logical :: ok
        real(dp) :: logabs

        call broadcast_complex_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            value = complex_nan()
            return
        end if
        call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
            reshape(weights, [size(weights)]))
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_4d_all_weighted

    pure function logsumexp_complex_4d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :)

        complex(dp) :: phase
        complex(dp), allocatable :: terms(:, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :)
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        logical :: ok
        real(dp) :: logabs

        if (dim < 1 .or. dim > 4) then
            allocate(values(0, 0, 0))
            return
        end if
        call broadcast_complex_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            select case (dim)
            case (1)
                allocate(values(size(x, 2), size(x, 3), size(x, 4)))
            case (2)
                allocate(values(size(x, 1), size(x, 3), size(x, 4)))
            case (3)
                allocate(values(size(x, 1), size(x, 2), size(x, 4)))
            case (4)
                allocate(values(size(x, 1), size(x, 2), size(x, 3)))
            end select
            values = complex_nan()
            return
        end if

        select case (dim)
        case (1)
            allocate(values(size(terms, 2), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        call logsumexp_complex_core(terms(:, j, k, l), logabs, phase, &
                            weights(:, j, k, l))
                        values(j, k, l) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (2)
            allocate(values(size(terms, 1), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do i = 1, size(terms, 1)
                        call logsumexp_complex_core(terms(i, :, k, l), logabs, phase, &
                            weights(i, :, k, l))
                        values(i, k, l) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (3)
            allocate(values(size(terms, 1), size(terms, 2), size(terms, 4)))
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_complex_core(terms(i, j, :, l), logabs, phase, &
                            weights(i, j, :, l))
                        values(i, j, l) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        case (4)
            allocate(values(size(terms, 1), size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_complex_core(terms(i, j, k, :), logabs, phase, &
                            weights(i, j, k, :))
                        values(i, j, k) = complex_log_from_phase(logabs, phase)
                    end do
                end do
            end do
        end select
    end function logsumexp_complex_4d_dim_weighted

    pure subroutine logsumexp_signed_complex_4d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        complex(dp), allocatable :: terms(:, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :)
        logical :: ok

        call broadcast_complex_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
            reshape(weights, [size(weights)]))
    end subroutine logsumexp_signed_complex_4d_all_weighted

    pure subroutine logsumexp_signed_complex_4d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        complex(dp), allocatable :: terms(:, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :)
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        logical :: ok

        if (dim < 1 .or. dim > 4) then
            allocate(logabs(0, 0, 0), phase(0, 0, 0))
            return
        end if
        call broadcast_complex_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            select case (dim)
            case (1)
                allocate(logabs(size(x, 2), size(x, 3), size(x, 4)), &
                    phase(size(x, 2), size(x, 3), size(x, 4)))
            case (2)
                allocate(logabs(size(x, 1), size(x, 3), size(x, 4)), &
                    phase(size(x, 1), size(x, 3), size(x, 4)))
            case (3)
                allocate(logabs(size(x, 1), size(x, 2), size(x, 4)), &
                    phase(size(x, 1), size(x, 2), size(x, 4)))
            case (4)
                allocate(logabs(size(x, 1), size(x, 2), size(x, 3)), &
                    phase(size(x, 1), size(x, 2), size(x, 3)))
            end select
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if

        select case (dim)
        case (1)
            allocate(logabs(size(terms, 2), size(terms, 3), size(terms, 4)), &
                phase(size(terms, 2), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        call logsumexp_complex_core(terms(:, j, k, l), logabs(j, k, l), &
                            phase(j, k, l), weights(:, j, k, l))
                    end do
                end do
            end do
        case (2)
            allocate(logabs(size(terms, 1), size(terms, 3), size(terms, 4)), &
                phase(size(terms, 1), size(terms, 3), size(terms, 4)))
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do i = 1, size(terms, 1)
                        call logsumexp_complex_core(terms(i, :, k, l), logabs(i, k, l), &
                            phase(i, k, l), weights(i, :, k, l))
                    end do
                end do
            end do
        case (3)
            allocate(logabs(size(terms, 1), size(terms, 2), size(terms, 4)), &
                phase(size(terms, 1), size(terms, 2), size(terms, 4)))
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_complex_core(terms(i, j, :, l), logabs(i, j, l), &
                            phase(i, j, l), weights(i, j, :, l))
                    end do
                end do
            end do
        case (4)
            allocate(logabs(size(terms, 1), size(terms, 2), size(terms, 3)), &
                phase(size(terms, 1), size(terms, 2), size(terms, 3)))
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    do i = 1, size(terms, 1)
                        call logsumexp_complex_core(terms(i, j, k, :), logabs(i, j, k), &
                            phase(i, j, k), weights(i, j, k, :))
                    end do
                end do
            end do
        end select
    end subroutine logsumexp_signed_complex_4d_dim_weighted

    pure function logsumexp_keepdims_complex_4d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :, :)

        complex(dp) :: phase_value
        complex(dp), allocatable :: terms(:, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :)
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        logical :: ok
        real(dp) :: logabs_value

        mode = reduction_dims_mode_4d(dims)
        if (mode < 0) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        call broadcast_complex_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            call allocate_keepdims_complex_4d(values, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), mode)
            if (size(values) > 0) values = complex_nan()
            return
        end if
        call allocate_keepdims_complex_4d(values, size(terms, 1), size(terms, 2), &
            size(terms, 3), size(terms, 4), mode)

        select case (mode)
        case (0)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        do i = 1, size(terms, 1)
                            call logsumexp_complex_core([terms(i, j, k, l)], logabs_value, &
                                phase_value, [weights(i, j, k, l)])
                            values(i, j, k, l) = complex_log_from_phase(logabs_value, phase_value)
                        end do
                    end do
                end do
            end do
        case (1)
            values(1, :, :, :) = logsumexp_complex_4d_dim_weighted(terms, 1, weights)
        case (2)
            values(:, 1, :, :) = logsumexp_complex_4d_dim_weighted(terms, 2, weights)
        case (4)
            values(:, :, 1, :) = logsumexp_complex_4d_dim_weighted(terms, 3, weights)
        case (8)
            values(:, :, :, 1) = logsumexp_complex_4d_dim_weighted(terms, 4, weights)
        case (3)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    values(1, 1, k, l) = logsumexp_complex_2d_all(terms(:, :, k, l), weights(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    values(1, j, 1, l) = logsumexp_complex_2d_all(terms(:, j, :, l), weights(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(terms, 4)
                do i = 1, size(terms, 1)
                    values(i, 1, 1, l) = logsumexp_complex_2d_all(terms(i, :, :, l), weights(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    values(1, j, k, 1) = logsumexp_complex_2d_all(terms(:, j, k, :), weights(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    values(i, 1, k, 1) = logsumexp_complex_2d_all(terms(i, :, k, :), weights(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    values(i, j, 1, 1) = logsumexp_complex_2d_all(terms(i, j, :, :), weights(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(terms, 4)
                values(1, 1, 1, l) = logsumexp_complex_3d_all(terms(:, :, :, l), weights(:, :, :, l))
            end do
        case (11)
            do k = 1, size(terms, 3)
                values(1, 1, k, 1) = logsumexp_complex_3d_all(terms(:, :, k, :), weights(:, :, k, :))
            end do
        case (13)
            do j = 1, size(terms, 2)
                values(1, j, 1, 1) = logsumexp_complex_3d_all(terms(:, j, :, :), weights(:, j, :, :))
            end do
        case (14)
            do i = 1, size(terms, 1)
                values(i, 1, 1, 1) = logsumexp_complex_3d_all(terms(i, :, :, :), weights(i, :, :, :))
            end do
        case (15)
            values(1, 1, 1, 1) = logsumexp_complex_4d_all_weighted(terms, weights)
        end select
    end function logsumexp_keepdims_complex_4d_weighted

    pure subroutine logsumexp_signed_keepdims_complex_4d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        complex(dp), allocatable :: reduced_phase(:, :, :)
        complex(dp), allocatable :: terms(:, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :)
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: mode
        logical :: ok
        real(dp), allocatable :: reduced_logabs(:, :, :)

        mode = reduction_dims_mode_4d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0), phase(0, 0, 0, 0))
            return
        end if
        call broadcast_complex_terms_4d(x, b, terms, weights, ok)
        if (.not. ok) then
            call allocate_keepdims_real_4d(logabs, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), mode)
            call allocate_keepdims_complex_4d(phase, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), mode)
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
            return
        end if
        call allocate_keepdims_real_4d(logabs, size(terms, 1), size(terms, 2), &
            size(terms, 3), size(terms, 4), mode)
        call allocate_keepdims_complex_4d(phase, size(terms, 1), size(terms, 2), &
            size(terms, 3), size(terms, 4), mode)

        select case (mode)
        case (0)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    do j = 1, size(terms, 2)
                        do i = 1, size(terms, 1)
                            call logsumexp_complex_core([terms(i, j, k, l)], logabs(i, j, k, l), &
                                phase(i, j, k, l), [weights(i, j, k, l)])
                        end do
                    end do
                end do
            end do
        case (1, 2, 4, 8)
            if (mode == 1) call logsumexp_signed_complex_4d_dim_weighted(terms, 1, reduced_logabs, reduced_phase, weights)
            if (mode == 2) call logsumexp_signed_complex_4d_dim_weighted(terms, 2, reduced_logabs, reduced_phase, weights)
            if (mode == 4) call logsumexp_signed_complex_4d_dim_weighted(terms, 3, reduced_logabs, reduced_phase, weights)
            if (mode == 8) call logsumexp_signed_complex_4d_dim_weighted(terms, 4, reduced_logabs, reduced_phase, weights)
            logabs = reshape(reduced_logabs, shape(logabs))
            phase = reshape(reduced_phase, shape(phase))
        case (3)
            do l = 1, size(terms, 4)
                do k = 1, size(terms, 3)
                    call logsumexp_signed_complex_2d_all(terms(:, :, k, l), logabs(1, 1, k, l), &
                        phase(1, 1, k, l), weights(:, :, k, l))
                end do
            end do
        case (5)
            do l = 1, size(terms, 4)
                do j = 1, size(terms, 2)
                    call logsumexp_signed_complex_2d_all(terms(:, j, :, l), logabs(1, j, 1, l), &
                        phase(1, j, 1, l), weights(:, j, :, l))
                end do
            end do
        case (6)
            do l = 1, size(terms, 4)
                do i = 1, size(terms, 1)
                    call logsumexp_signed_complex_2d_all(terms(i, :, :, l), logabs(i, 1, 1, l), &
                        phase(i, 1, 1, l), weights(i, :, :, l))
                end do
            end do
        case (9)
            do k = 1, size(terms, 3)
                do j = 1, size(terms, 2)
                    call logsumexp_signed_complex_2d_all(terms(:, j, k, :), logabs(1, j, k, 1), &
                        phase(1, j, k, 1), weights(:, j, k, :))
                end do
            end do
        case (10)
            do k = 1, size(terms, 3)
                do i = 1, size(terms, 1)
                    call logsumexp_signed_complex_2d_all(terms(i, :, k, :), logabs(i, 1, k, 1), &
                        phase(i, 1, k, 1), weights(i, :, k, :))
                end do
            end do
        case (12)
            do j = 1, size(terms, 2)
                do i = 1, size(terms, 1)
                    call logsumexp_signed_complex_2d_all(terms(i, j, :, :), logabs(i, j, 1, 1), &
                        phase(i, j, 1, 1), weights(i, j, :, :))
                end do
            end do
        case (7)
            do l = 1, size(terms, 4)
                call logsumexp_signed_complex_3d_all(terms(:, :, :, l), logabs(1, 1, 1, l), &
                    phase(1, 1, 1, l), weights(:, :, :, l))
            end do
        case (11)
            do k = 1, size(terms, 3)
                call logsumexp_signed_complex_3d_all(terms(:, :, k, :), logabs(1, 1, k, 1), &
                    phase(1, 1, k, 1), weights(:, :, k, :))
            end do
        case (13)
            do j = 1, size(terms, 2)
                call logsumexp_signed_complex_3d_all(terms(:, j, :, :), logabs(1, j, 1, 1), &
                    phase(1, j, 1, 1), weights(:, j, :, :))
            end do
        case (14)
            do i = 1, size(terms, 1)
                call logsumexp_signed_complex_3d_all(terms(i, :, :, :), logabs(i, 1, 1, 1), &
                    phase(i, 1, 1, 1), weights(i, :, :, :))
            end do
        case (15)
            call logsumexp_signed_complex_4d_all_weighted(terms, logabs(1, 1, 1, 1), &
                phase(1, 1, 1, 1), weights)
        end select
    end subroutine logsumexp_signed_keepdims_complex_4d_weighted

    pure function logsumexp_complex_4d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp) :: value

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)

        select rank (b)
        rank (0)
            value = logsumexp_complex_4d_all_weighted(x, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_4d_all_weighted(x, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_4d_all_weighted(x, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_4d_all_weighted(x, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_4d_all_weighted(x, b4)
        rank default
            value = complex_nan()
        end select
    end function logsumexp_complex_4d_all_real_b

    pure function logsumexp_complex_4d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: values(:, :, :)

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)

        select rank (b)
        rank (0)
            values = logsumexp_complex_4d_dim_weighted(x, dim, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_4d_dim_weighted(x, dim, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_4d_dim_weighted(x, dim, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_4d_dim_weighted(x, dim, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_4d_dim_weighted(x, dim, b4)
        rank default
            if (dim < 1 .or. dim > 4) then
                allocate(values(0, 0, 0))
            else
                select case (dim)
                case (1)
                    allocate(values(size(x, 2), size(x, 3), size(x, 4)))
                case (2)
                    allocate(values(size(x, 1), size(x, 3), size(x, 4)))
                case (3)
                    allocate(values(size(x, 1), size(x, 2), size(x, 4)))
                case (4)
                    allocate(values(size(x, 1), size(x, 2), size(x, 3)))
                end select
                values = complex_nan()
            end if
        end select
    end function logsumexp_complex_4d_dim_real_b

    pure subroutine logsumexp_signed_complex_4d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        complex(dp), intent(out) :: phase !! unit phase of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_4d_all_weighted(x, logabs, phase, &
                cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_all_weighted(x, logabs, phase, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_all_weighted(x, logabs, phase, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_all_weighted(x, logabs, phase, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_all_weighted(x, logabs, phase, b4)
        rank default
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_4d_all_real_b

    pure subroutine logsumexp_signed_complex_4d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 4
        real(dp), allocatable, intent(out) :: logabs(:, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_4d_dim_weighted(x, dim, logabs, phase, &
                cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_dim_weighted(x, dim, logabs, phase, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_dim_weighted(x, dim, logabs, phase, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_dim_weighted(x, dim, logabs, phase, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_4d_dim_weighted(x, dim, logabs, phase, b4)
        rank default
            if (dim < 1 .or. dim > 4) then
                allocate(logabs(0, 0, 0), phase(0, 0, 0))
            else
                select case (dim)
                case (1)
                    allocate(logabs(size(x, 2), size(x, 3), size(x, 4)), &
                        phase(size(x, 2), size(x, 3), size(x, 4)))
                case (2)
                    allocate(logabs(size(x, 1), size(x, 3), size(x, 4)), &
                        phase(size(x, 1), size(x, 3), size(x, 4)))
                case (3)
                    allocate(logabs(size(x, 1), size(x, 2), size(x, 4)), &
                        phase(size(x, 1), size(x, 2), size(x, 4)))
                case (4)
                    allocate(logabs(size(x, 1), size(x, 2), size(x, 3)), &
                        phase(size(x, 1), size(x, 2), size(x, 3)))
                end select
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
            end if
        end select
    end subroutine logsumexp_signed_complex_4d_dim_real_b

    pure function logsumexp_keepdims_complex_4d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: values(:, :, :, :)

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)

        select rank (b)
        rank (0)
            values = logsumexp_keepdims_complex_4d_weighted(x, dims, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_4d_weighted(x, dims, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_4d_weighted(x, dims, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_4d_weighted(x, dims, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_4d_weighted(x, dims, b4)
        rank default
            call allocate_keepdims_complex_4d(values, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), reduction_dims_mode_4d(dims))
            if (size(values) > 0) values = complex_nan()
        end select
    end function logsumexp_keepdims_complex_4d_real_b

    pure subroutine logsumexp_signed_keepdims_complex_4d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_keepdims_complex_4d_weighted(x, dims, logabs, phase, &
                cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_4d_weighted(x, dims, logabs, phase, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_4d_weighted(x, dims, logabs, phase, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_4d_weighted(x, dims, logabs, phase, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_4d_weighted(x, dims, logabs, phase, b4)
        rank default
            call allocate_keepdims_real_4d(logabs, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), reduction_dims_mode_4d(dims))
            call allocate_keepdims_complex_4d(phase, size(x, 1), size(x, 2), size(x, 3), &
                size(x, 4), reduction_dims_mode_4d(dims))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
        end select
    end subroutine logsumexp_signed_keepdims_complex_4d_real_b

    pure subroutine broadcast_real_terms_4d(x, b, terms, weights, ok)
        real(dp), intent(in) :: x(:, :, :, :) !! rank-4 real logarithmic terms to broadcast
        real(dp), intent(in) :: b(..) !! scalar or trailing-array weights broadcast jointly with x
        real(dp), allocatable, intent(out) :: terms(:, :, :, :) !! terms broadcast to the joint shape
        real(dp), allocatable, intent(out) :: weights(:, :, :, :) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: b3
        integer :: b4
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: k
        integer :: kb
        integer :: kx
        integer :: l
        integer :: lb
        integer :: lx
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2), size(x, 3), size(x, 4)), &
                weights(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = 1
            b3 = 1
            b4 = size(b)
        rank (2)
            b1 = 1
            b2 = 1
            b3 = size(b, 1)
            b4 = size(b, 2)
        rank (3)
            b1 = 1
            b2 = size(b, 1)
            b3 = size(b, 2)
            b4 = size(b, 3)
        rank (4)
            b1 = size(b, 1)
            b2 = size(b, 2)
            b3 = size(b, 3)
            b4 = size(b, 4)
        rank default
            ok = .false.
            allocate(terms(0, 0, 0, 0), weights(0, 0, 0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1) .or. &
                (size(x, 3) /= b3 .and. size(x, 3) /= 1 .and. b3 /= 1) .or. &
                (size(x, 4) /= b4 .and. size(x, 4) /= 1 .and. b4 /= 1)) then
            ok = .false.
            allocate(terms(0, 0, 0, 0), weights(0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        if (n3 == 1) n3 = b3
        if (n4 == 1) n4 = b4
        allocate(terms(n1, n2, n3, n4), weights(n1, n2, n3, n4))

        select rank (b)
        rank (1)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(lb)
                        end do
                    end do
                end do
            end do
        rank (2)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b, 2) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    kb = merge(1, k, size(b, 1) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(kb, lb)
                        end do
                    end do
                end do
            end do
        rank (3)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b, 3) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    kb = merge(1, k, size(b, 2) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        jb = merge(1, j, size(b, 1) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(jb, kb, lb)
                        end do
                    end do
                end do
            end do
        rank (4)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b, 4) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    kb = merge(1, k, size(b, 3) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        jb = merge(1, j, size(b, 2) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            ib = merge(1, i, size(b, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(ib, jb, kb, lb)
                        end do
                    end do
                end do
            end do
        end select
    end subroutine broadcast_real_terms_4d

    pure subroutine broadcast_complex_terms_4d(x, b, terms, weights, ok)
        complex(dp), intent(in) :: x(:, :, :, :) !! rank-4 complex logarithmic terms to broadcast
        complex(dp), intent(in) :: b(..) !! scalar or trailing-array weights broadcast jointly with x
        complex(dp), allocatable, intent(out) :: terms(:, :, :, :) !! terms broadcast to the joint shape
        complex(dp), allocatable, intent(out) :: weights(:, :, :, :) !! weights broadcast to the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: b3
        integer :: b4
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: k
        integer :: kb
        integer :: kx
        integer :: l
        integer :: lb
        integer :: lx
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2), size(x, 3), size(x, 4)), &
                weights(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = 1
            b3 = 1
            b4 = size(b)
        rank (2)
            b1 = 1
            b2 = 1
            b3 = size(b, 1)
            b4 = size(b, 2)
        rank (3)
            b1 = 1
            b2 = size(b, 1)
            b3 = size(b, 2)
            b4 = size(b, 3)
        rank (4)
            b1 = size(b, 1)
            b2 = size(b, 2)
            b3 = size(b, 3)
            b4 = size(b, 4)
        rank default
            ok = .false.
            allocate(terms(0, 0, 0, 0), weights(0, 0, 0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1) .or. &
                (size(x, 3) /= b3 .and. size(x, 3) /= 1 .and. b3 /= 1) .or. &
                (size(x, 4) /= b4 .and. size(x, 4) /= 1 .and. b4 /= 1)) then
            ok = .false.
            allocate(terms(0, 0, 0, 0), weights(0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        if (n3 == 1) n3 = b3
        if (n4 == 1) n4 = b4
        allocate(terms(n1, n2, n3, n4), weights(n1, n2, n3, n4))

        select rank (b)
        rank (1)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(lb)
                        end do
                    end do
                end do
            end do
        rank (2)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b, 2) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    kb = merge(1, k, size(b, 1) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(kb, lb)
                        end do
                    end do
                end do
            end do
        rank (3)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b, 3) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    kb = merge(1, k, size(b, 2) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        jb = merge(1, j, size(b, 1) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(jb, kb, lb)
                        end do
                    end do
                end do
            end do
        rank (4)
            do l = 1, n4
                lx = merge(1, l, size(x, 4) == 1)
                lb = merge(1, l, size(b, 4) == 1)
                do k = 1, n3
                    kx = merge(1, k, size(x, 3) == 1)
                    kb = merge(1, k, size(b, 3) == 1)
                    do j = 1, n2
                        jx = merge(1, j, size(x, 2) == 1)
                        jb = merge(1, j, size(b, 2) == 1)
                        do i = 1, n1
                            ix = merge(1, i, size(x, 1) == 1)
                            ib = merge(1, i, size(b, 1) == 1)
                            terms(i, j, k, l) = x(ix, jx, kx, lx)
                            weights(i, j, k, l) = b(ib, jb, kb, lb)
                        end do
                    end do
                end do
            end do
        end select
    end subroutine broadcast_complex_terms_4d


    pure subroutine broadcast_real_terms_5d(x, b, terms, weights, ok)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 logarithmic terms to broadcast
        real(dp), intent(in) :: b(..) !! scalar or trailing-array weights broadcast jointly with x
        real(dp), allocatable, intent(out) :: terms(:, :, :, :, :) !! terms at the joint shape
        real(dp), allocatable, intent(out) :: weights(:, :, :, :, :) !! weights at the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: b3
        integer :: b4
        integer :: b5
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: k
        integer :: kb
        integer :: kx
        integer :: l
        integer :: lb
        integer :: lx
        integer :: m
        integer :: mb
        integer :: mx
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)), &
                weights(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = 1
            b3 = 1
            b4 = 1
            b5 = size(b, 1)
        rank (2)
            b1 = 1
            b2 = 1
            b3 = 1
            b4 = size(b, 1)
            b5 = size(b, 2)
        rank (3)
            b1 = 1
            b2 = 1
            b3 = size(b, 1)
            b4 = size(b, 2)
            b5 = size(b, 3)
        rank (4)
            b1 = 1
            b2 = size(b, 1)
            b3 = size(b, 2)
            b4 = size(b, 3)
            b5 = size(b, 4)
        rank (5)
            b1 = size(b, 1)
            b2 = size(b, 2)
            b3 = size(b, 3)
            b4 = size(b, 4)
            b5 = size(b, 5)
        rank default
            ok = .false.
            allocate(terms(0, 0, 0, 0, 0), weights(0, 0, 0, 0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1) .or. &
                (size(x, 3) /= b3 .and. size(x, 3) /= 1 .and. b3 /= 1) .or. &
                (size(x, 4) /= b4 .and. size(x, 4) /= 1 .and. b4 /= 1) .or. &
                (size(x, 5) /= b5 .and. size(x, 5) /= 1 .and. b5 /= 1)) then
            ok = .false.
            allocate(terms(0, 0, 0, 0, 0), weights(0, 0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        if (n3 == 1) n3 = b3
        if (n4 == 1) n4 = b4
        if (n5 == 1) n5 = b5
        allocate(terms(n1, n2, n3, n4, n5), weights(n1, n2, n3, n4, n5))

        select rank (b)
        rank (1)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 1) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (2)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 2) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 1) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (3)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 3) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 2) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        kb = merge(1, k, size(b, 1) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(kb, lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (4)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 4) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 3) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        kb = merge(1, k, size(b, 2) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            jb = merge(1, j, size(b, 1) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(jb, kb, lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (5)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 5) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 4) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        kb = merge(1, k, size(b, 3) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            jb = merge(1, j, size(b, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                ib = merge(1, i, size(b, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(ib, jb, kb, lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        end select
    end subroutine broadcast_real_terms_5d

    pure subroutine broadcast_complex_terms_5d(x, b, terms, weights, ok)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 logarithmic terms to broadcast
        complex(dp), intent(in) :: b(..) !! scalar or trailing-array weights broadcast jointly with x
        complex(dp), allocatable, intent(out) :: terms(:, :, :, :, :) !! terms at the joint shape
        complex(dp), allocatable, intent(out) :: weights(:, :, :, :, :) !! weights at the joint shape
        logical, intent(out) :: ok !! true when x and b have a NumPy-compatible joint shape

        integer :: b1
        integer :: b2
        integer :: b3
        integer :: b4
        integer :: b5
        integer :: i
        integer :: ib
        integer :: ix
        integer :: j
        integer :: jb
        integer :: jx
        integer :: k
        integer :: kb
        integer :: kx
        integer :: l
        integer :: lb
        integer :: lx
        integer :: m
        integer :: mb
        integer :: mx
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5

        ok = .true.
        select rank (b)
        rank (0)
            allocate(terms(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)), &
                weights(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
            terms = x
            weights = b
            return
        rank (1)
            b1 = 1
            b2 = 1
            b3 = 1
            b4 = 1
            b5 = size(b, 1)
        rank (2)
            b1 = 1
            b2 = 1
            b3 = 1
            b4 = size(b, 1)
            b5 = size(b, 2)
        rank (3)
            b1 = 1
            b2 = 1
            b3 = size(b, 1)
            b4 = size(b, 2)
            b5 = size(b, 3)
        rank (4)
            b1 = 1
            b2 = size(b, 1)
            b3 = size(b, 2)
            b4 = size(b, 3)
            b5 = size(b, 4)
        rank (5)
            b1 = size(b, 1)
            b2 = size(b, 2)
            b3 = size(b, 3)
            b4 = size(b, 4)
            b5 = size(b, 5)
        rank default
            ok = .false.
            allocate(terms(0, 0, 0, 0, 0), weights(0, 0, 0, 0, 0))
            return
        end select

        if ((size(x, 1) /= b1 .and. size(x, 1) /= 1 .and. b1 /= 1) .or. &
                (size(x, 2) /= b2 .and. size(x, 2) /= 1 .and. b2 /= 1) .or. &
                (size(x, 3) /= b3 .and. size(x, 3) /= 1 .and. b3 /= 1) .or. &
                (size(x, 4) /= b4 .and. size(x, 4) /= 1 .and. b4 /= 1) .or. &
                (size(x, 5) /= b5 .and. size(x, 5) /= 1 .and. b5 /= 1)) then
            ok = .false.
            allocate(terms(0, 0, 0, 0, 0), weights(0, 0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        if (n1 == 1) n1 = b1
        if (n2 == 1) n2 = b2
        if (n3 == 1) n3 = b3
        if (n4 == 1) n4 = b4
        if (n5 == 1) n5 = b5
        allocate(terms(n1, n2, n3, n4, n5), weights(n1, n2, n3, n4, n5))

        select rank (b)
        rank (1)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 1) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (2)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 2) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 1) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (3)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 3) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 2) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        kb = merge(1, k, size(b, 1) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(kb, lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (4)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 4) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 3) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        kb = merge(1, k, size(b, 2) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            jb = merge(1, j, size(b, 1) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(jb, kb, lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        rank (5)
            do m = 1, n5
                mx = merge(1, m, size(x, 5) == 1)
                mb = merge(1, m, size(b, 5) == 1)
                do l = 1, n4
                    lx = merge(1, l, size(x, 4) == 1)
                    lb = merge(1, l, size(b, 4) == 1)
                    do k = 1, n3
                        kx = merge(1, k, size(x, 3) == 1)
                        kb = merge(1, k, size(b, 3) == 1)
                        do j = 1, n2
                            jx = merge(1, j, size(x, 2) == 1)
                            jb = merge(1, j, size(b, 2) == 1)
                            do i = 1, n1
                                ix = merge(1, i, size(x, 1) == 1)
                                ib = merge(1, i, size(b, 1) == 1)
                                terms(i, j, k, l, m) = x(ix, jx, kx, lx, mx)
                                weights(i, j, k, l, m) = b(ib, jb, kb, lb, mb)
                            end do
                        end do
                    end do
                end do
            end do
        end select
    end subroutine broadcast_complex_terms_5d


    pure function reduction_dims_mode_5d(dims) result(mode)
        integer, intent(in) :: dims(:) !! Fortran dimensions to reduce
        integer :: mode

        integer :: i
        logical :: seen(5)

        mode = 0
        seen = .false.
        if (size(dims) > 5) then
            mode = -1
            return
        end if
        do i = 1, size(dims)
            if (dims(i) < 1 .or. dims(i) > 5) then
                mode = -1
                return
            end if
            if (seen(dims(i))) then
                mode = -1
                return
            end if
            seen(dims(i)) = .true.
            mode = mode + 2 ** (dims(i) - 1)
        end do
    end function reduction_dims_mode_5d

    pure subroutine allocate_keepdims_real_5d(values, n1, n2, n3, n4, n5, mode)
        real(dp), allocatable, intent(out) :: values(:, :, :, :, :) !! allocated keepdims result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: n3 !! unreduced third-dimension extent
        integer, intent(in) :: n4 !! unreduced fourth-dimension extent
        integer, intent(in) :: n5 !! unreduced fifth-dimension extent
        integer, intent(in) :: mode !! rank-5 reduction bit mask, or negative when invalid

        if (mode < 0) then
            allocate(values(0, 0, 0, 0, 0))
        else
            allocate(values(merge(1, n1, btest(mode, 0)), &
                merge(1, n2, btest(mode, 1)), merge(1, n3, btest(mode, 2)), &
                merge(1, n4, btest(mode, 3)), merge(1, n5, btest(mode, 4))))
        end if
    end subroutine allocate_keepdims_real_5d

    pure subroutine allocate_keepdims_complex_5d(values, n1, n2, n3, n4, n5, mode)
        complex(dp), allocatable, intent(out) :: values(:, :, :, :, :) !! allocated keepdims result
        integer, intent(in) :: n1 !! unreduced first-dimension extent
        integer, intent(in) :: n2 !! unreduced second-dimension extent
        integer, intent(in) :: n3 !! unreduced third-dimension extent
        integer, intent(in) :: n4 !! unreduced fourth-dimension extent
        integer, intent(in) :: n5 !! unreduced fifth-dimension extent
        integer, intent(in) :: mode !! rank-5 reduction bit mask, or negative when invalid

        if (mode < 0) then
            allocate(values(0, 0, 0, 0, 0))
        else
            allocate(values(merge(1, n1, btest(mode, 0)), &
                merge(1, n2, btest(mode, 1)), merge(1, n3, btest(mode, 2)), &
                merge(1, n4, btest(mode, 3)), merge(1, n5, btest(mode, 4))))
        end if
    end subroutine allocate_keepdims_complex_5d

    pure subroutine logsumexp_signed_keepdims_5d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! log-absolute sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :) !! signs of reduced sums

        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: p
        integer :: reduction_size
        real(dp), allocatable :: buffer(:)

        mode = reduction_dims_mode_5d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0, 0), sign_value(0, 0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        allocate(logabs(m1, m2, m3, m4, m5), sign_value(m1, m2, m3, m4, m5))

        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        allocate(buffer(max(1, reduction_size)))

        do o5 = 1, m5
            lo5 = merge(1, o5, btest(mode, 4))
            hi5 = merge(n5, o5, btest(mode, 4))
            do o4 = 1, m4
                lo4 = merge(1, o4, btest(mode, 3))
                hi4 = merge(n4, o4, btest(mode, 3))
                do o3 = 1, m3
                    lo3 = merge(1, o3, btest(mode, 2))
                    hi3 = merge(n3, o3, btest(mode, 2))
                    do o2 = 1, m2
                        lo2 = merge(1, o2, btest(mode, 1))
                        hi2 = merge(n2, o2, btest(mode, 1))
                        do o1 = 1, m1
                            lo1 = merge(1, o1, btest(mode, 0))
                            hi1 = merge(n1, o1, btest(mode, 0))
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                buffer(p) = x(i1, i2, i3, i4, i5)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            call logsumexp_core(buffer(1:p), logabs(o1, o2, o3, o4, o5), &
                                sign_value(o1, o2, o3, o4, o5))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine logsumexp_signed_keepdims_5d

    pure function logsumexp_keepdims_5d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :)

        call logsumexp_signed_keepdims_5d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_5d

    pure function logsumexp_5d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_core(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_5d_all

    pure function logsumexp_5d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable :: values(:, :, :, :)

        real(dp), allocatable :: kept(:, :, :, :, :)

        if (dim < 1 .or. dim > 5) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
        case (5)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        end select
        kept = logsumexp_keepdims_5d(x, [dim])
        values = reshape(kept, shape(values))
    end function logsumexp_5d_dim

    pure subroutine logsumexp_signed_5d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum: zero or +1

        call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_5d_all

    pure subroutine logsumexp_signed_5d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: kept_logabs(:, :, :, :, :)
        real(dp), allocatable :: kept_sign(:, :, :, :, :)

        if (dim < 1 .or. dim > 5) then
            allocate(logabs(0, 0, 0, 0), sign_value(0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(logabs(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
            allocate(sign_value(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        case (2)
            allocate(logabs(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
            allocate(sign_value(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
        case (3)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
            allocate(sign_value(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
        case (4)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
            allocate(sign_value(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
        case (5)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
            allocate(sign_value(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        end select
        call logsumexp_signed_keepdims_5d(x, [dim], kept_logabs, kept_sign)
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign, shape(sign_value))
    end subroutine logsumexp_signed_5d_dim

    pure subroutine logsumexp_signed_keepdims_5d_weighted(x, dims, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        real(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: p
        integer :: reduction_size
        logical :: ok
        real(dp), allocatable :: terms(:, :, :, :, :)
        real(dp), allocatable :: weight_buffer(:)
        real(dp), allocatable :: weights(:, :, :, :, :)

        mode = reduction_dims_mode_5d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0, 0), sign_value(0, 0, 0, 0, 0))
            return
        end if

        call broadcast_real_terms_5d(x, b, terms, weights, ok)
        if (.not. ok) then
            n1 = size(x, 1)
            n2 = size(x, 2)
            n3 = size(x, 3)
            n4 = size(x, 4)
            n5 = size(x, 5)
            m1 = merge(1, n1, btest(mode, 0))
            m2 = merge(1, n2, btest(mode, 1))
            m3 = merge(1, n3, btest(mode, 2))
            m4 = merge(1, n4, btest(mode, 3))
            m5 = merge(1, n5, btest(mode, 4))
            allocate(logabs(m1, m2, m3, m4, m5), sign_value(m1, m2, m3, m4, m5))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(sign_value) > 0) sign_value = quiet_nan(0.0_dp)
            return
        end if

        n1 = size(terms, 1)
        n2 = size(terms, 2)
        n3 = size(terms, 3)
        n4 = size(terms, 4)
        n5 = size(terms, 5)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        allocate(logabs(m1, m2, m3, m4, m5), sign_value(m1, m2, m3, m4, m5))

        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        allocate(buffer(max(1, reduction_size)), weight_buffer(max(1, reduction_size)))

        do o5 = 1, m5
            lo5 = merge(1, o5, btest(mode, 4))
            hi5 = merge(n5, o5, btest(mode, 4))
            do o4 = 1, m4
                lo4 = merge(1, o4, btest(mode, 3))
                hi4 = merge(n4, o4, btest(mode, 3))
                do o3 = 1, m3
                    lo3 = merge(1, o3, btest(mode, 2))
                    hi3 = merge(n3, o3, btest(mode, 2))
                    do o2 = 1, m2
                        lo2 = merge(1, o2, btest(mode, 1))
                        hi2 = merge(n2, o2, btest(mode, 1))
                        do o1 = 1, m1
                            lo1 = merge(1, o1, btest(mode, 0))
                            hi1 = merge(n1, o1, btest(mode, 0))
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                buffer(p) = terms(i1, i2, i3, i4, i5)
                                                weight_buffer(p) = weights(i1, i2, i3, i4, i5)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            call logsumexp_core(buffer(1:p), logabs(o1, o2, o3, o4, o5), &
                                sign_value(o1, o2, o3, o4, o5), weight_buffer(1:p))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine logsumexp_signed_keepdims_5d_weighted

    pure function logsumexp_keepdims_5d_weighted(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        real(dp), allocatable :: values(:, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :)

        call logsumexp_signed_keepdims_5d_weighted(x, dims, values, signs, b)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_5d_weighted

    pure subroutine logsumexp_signed_5d_all_weighted(x, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the absolute weighted exponential sum
        real(dp), intent(out) :: sign_value !! sign of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        logical :: ok
        real(dp), allocatable :: terms(:, :, :, :, :)
        real(dp), allocatable :: weights(:, :, :, :, :)

        call broadcast_real_terms_5d(x, b, terms, weights, ok)
        if (.not. ok) then
            logabs = quiet_nan(0.0_dp)
            sign_value = quiet_nan(0.0_dp)
            return
        end if
        call logsumexp_core(reshape(terms, [size(terms)]), logabs, sign_value, &
            reshape(weights, [size(weights)]))
    end subroutine logsumexp_signed_5d_all_weighted

    pure function logsumexp_5d_all_weighted(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_signed_5d_all_weighted(x, value, sign_value, b)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_5d_all_weighted

    pure subroutine logsumexp_signed_5d_dim_weighted(x, dim, logabs, sign_value, b)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :) !! reduced signs
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        real(dp), allocatable :: kept_logabs(:, :, :, :, :)
        real(dp), allocatable :: kept_sign(:, :, :, :, :)

        if (dim < 1 .or. dim > 5) then
            allocate(logabs(0, 0, 0, 0), sign_value(0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_5d_weighted(x, [dim], kept_logabs, kept_sign, b)
        select case (dim)
        case (1)
            allocate(logabs(size(kept_logabs, 2), size(kept_logabs, 3), &
                size(kept_logabs, 4), size(kept_logabs, 5)))
            allocate(sign_value(size(kept_sign, 2), size(kept_sign, 3), &
                size(kept_sign, 4), size(kept_sign, 5)))
        case (2)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 3), &
                size(kept_logabs, 4), size(kept_logabs, 5)))
            allocate(sign_value(size(kept_sign, 1), size(kept_sign, 3), &
                size(kept_sign, 4), size(kept_sign, 5)))
        case (3)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 2), &
                size(kept_logabs, 4), size(kept_logabs, 5)))
            allocate(sign_value(size(kept_sign, 1), size(kept_sign, 2), &
                size(kept_sign, 4), size(kept_sign, 5)))
        case (4)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 2), &
                size(kept_logabs, 3), size(kept_logabs, 5)))
            allocate(sign_value(size(kept_sign, 1), size(kept_sign, 2), &
                size(kept_sign, 3), size(kept_sign, 5)))
        case (5)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 2), &
                size(kept_logabs, 3), size(kept_logabs, 4)))
            allocate(sign_value(size(kept_sign, 1), size(kept_sign, 2), &
                size(kept_sign, 3), size(kept_sign, 4)))
        end select
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign, shape(sign_value))
    end subroutine logsumexp_signed_5d_dim_weighted

    pure function logsumexp_5d_dim_weighted(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        real(dp), allocatable :: values(:, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :)

        call logsumexp_signed_5d_dim_weighted(x, dim, values, signs, b)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_5d_dim_weighted

    pure function softmax_5d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        if (present(dim)) then
            values = softmax_5d_dims(x, [dim])
        else
            values = softmax_5d_dims(x, [1, 2, 3, 4, 5])
        end if
    end function softmax_5d

    pure function softmax_5d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        call normalize_real_5d_dims(x, dims, .false., values)
    end function softmax_5d_dims

    pure function log_softmax_5d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        if (present(dim)) then
            values = log_softmax_5d_dims(x, [dim])
        else
            values = log_softmax_5d_dims(x, [1, 2, 3, 4, 5])
        end if
    end function log_softmax_5d

    pure function log_softmax_5d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        call normalize_real_5d_dims(x, dims, .true., values)
    end function log_softmax_5d_dims

    pure subroutine logsumexp_signed_keepdims_complex_5d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! log-magnitudes of sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :) !! unit phases of sums

        complex(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: p
        integer :: reduction_size

        mode = reduction_dims_mode_5d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0, 0), phase(0, 0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        allocate(logabs(m1, m2, m3, m4, m5), phase(m1, m2, m3, m4, m5))

        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        allocate(buffer(max(1, reduction_size)))

        do o5 = 1, m5
            lo5 = merge(1, o5, btest(mode, 4))
            hi5 = merge(n5, o5, btest(mode, 4))
            do o4 = 1, m4
                lo4 = merge(1, o4, btest(mode, 3))
                hi4 = merge(n4, o4, btest(mode, 3))
                do o3 = 1, m3
                    lo3 = merge(1, o3, btest(mode, 2))
                    hi3 = merge(n3, o3, btest(mode, 2))
                    do o2 = 1, m2
                        lo2 = merge(1, o2, btest(mode, 1))
                        hi2 = merge(n2, o2, btest(mode, 1))
                        do o1 = 1, m1
                            lo1 = merge(1, o1, btest(mode, 0))
                            hi1 = merge(n1, o1, btest(mode, 0))
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                buffer(p) = x(i1, i2, i3, i4, i5)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            call logsumexp_complex_core(buffer(1:p), &
                                logabs(o1, o2, o3, o4, o5), phase(o1, o2, o3, o4, o5))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine logsumexp_signed_keepdims_complex_5d

    pure function logsumexp_keepdims_complex_5d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :)

        complex(dp), allocatable :: phase(:, :, :, :, :)
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        real(dp), allocatable :: logabs(:, :, :, :, :)

        call logsumexp_signed_keepdims_complex_5d(x, dims, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), &
            size(logabs, 4), size(logabs, 5)))
        do i5 = 1, size(values, 5)
            do i4 = 1, size(values, 4)
                do i3 = 1, size(values, 3)
                    do i2 = 1, size(values, 2)
                        do i1 = 1, size(values, 1)
                            values(i1, i2, i3, i4, i5) = complex_log_from_phase( &
                                logabs(i1, i2, i3, i4, i5), phase(i1, i2, i3, i4, i5))
                        end do
                    end do
                end do
            end do
        end do
    end function logsumexp_keepdims_complex_5d

    pure function logsumexp_complex_5d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 values reduced over all dimensions
        complex(dp) :: value

        complex(dp) :: phase
        real(dp) :: logabs

        call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_5d_all

    pure function logsumexp_complex_5d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        complex(dp), allocatable :: values(:, :, :, :)

        complex(dp), allocatable :: kept(:, :, :, :, :)

        if (dim < 1 .or. dim > 5) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
        case (5)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        end select
        kept = logsumexp_keepdims_complex_5d(x, [dim])
        values = reshape(kept, shape(values))
    end function logsumexp_complex_5d_dim

    pure subroutine logsumexp_signed_complex_5d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the exponential-sum magnitude
        complex(dp), intent(out) :: phase !! unit complex phase of the exponential sum

        call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_5d_all

    pure subroutine logsumexp_signed_complex_5d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! log-magnitudes of reduced sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! unit phases of reduced sums

        real(dp), allocatable :: kept_logabs(:, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, :)

        if (dim < 1 .or. dim > 5) then
            allocate(logabs(0, 0, 0, 0), phase(0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(logabs(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
            allocate(phase(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        case (2)
            allocate(logabs(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
            allocate(phase(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
        case (3)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
            allocate(phase(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
        case (4)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
            allocate(phase(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
        case (5)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
            allocate(phase(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
        end select
        call logsumexp_signed_keepdims_complex_5d(x, [dim], kept_logabs, kept_phase)
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_5d_dim

    pure subroutine logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        complex(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: p
        integer :: reduction_size
        logical :: ok
        complex(dp), allocatable :: terms(:, :, :, :, :)
        complex(dp), allocatable :: weight_buffer(:)
        complex(dp), allocatable :: weights(:, :, :, :, :)

        mode = reduction_dims_mode_5d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0, 0), phase(0, 0, 0, 0, 0))
            return
        end if

        call broadcast_complex_terms_5d(x, b, terms, weights, ok)
        if (.not. ok) then
            n1 = size(x, 1)
            n2 = size(x, 2)
            n3 = size(x, 3)
            n4 = size(x, 4)
            n5 = size(x, 5)
            m1 = merge(1, n1, btest(mode, 0))
            m2 = merge(1, n2, btest(mode, 1))
            m3 = merge(1, n3, btest(mode, 2))
            m4 = merge(1, n4, btest(mode, 3))
            m5 = merge(1, n5, btest(mode, 4))
            allocate(logabs(m1, m2, m3, m4, m5), phase(m1, m2, m3, m4, m5))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
            return
        end if

        n1 = size(terms, 1)
        n2 = size(terms, 2)
        n3 = size(terms, 3)
        n4 = size(terms, 4)
        n5 = size(terms, 5)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        allocate(logabs(m1, m2, m3, m4, m5), phase(m1, m2, m3, m4, m5))

        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        allocate(buffer(max(1, reduction_size)), weight_buffer(max(1, reduction_size)))

        do o5 = 1, m5
            lo5 = merge(1, o5, btest(mode, 4))
            hi5 = merge(n5, o5, btest(mode, 4))
            do o4 = 1, m4
                lo4 = merge(1, o4, btest(mode, 3))
                hi4 = merge(n4, o4, btest(mode, 3))
                do o3 = 1, m3
                    lo3 = merge(1, o3, btest(mode, 2))
                    hi3 = merge(n3, o3, btest(mode, 2))
                    do o2 = 1, m2
                        lo2 = merge(1, o2, btest(mode, 1))
                        hi2 = merge(n2, o2, btest(mode, 1))
                        do o1 = 1, m1
                            lo1 = merge(1, o1, btest(mode, 0))
                            hi1 = merge(n1, o1, btest(mode, 0))
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                buffer(p) = terms(i1, i2, i3, i4, i5)
                                                weight_buffer(p) = weights(i1, i2, i3, i4, i5)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            call logsumexp_complex_core(buffer(1:p), &
                                logabs(o1, o2, o3, o4, o5), phase(o1, o2, o3, o4, o5), &
                                weight_buffer(1:p))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine logsumexp_signed_keepdims_complex_5d_weighted

    pure function logsumexp_keepdims_complex_5d_weighted(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :, :)
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5

        call logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), &
            size(logabs, 4), size(logabs, 5)))
        do i5 = 1, size(values, 5)
            do i4 = 1, size(values, 4)
                do i3 = 1, size(values, 3)
                    do i2 = 1, size(values, 2)
                        do i1 = 1, size(values, 1)
                            values(i1, i2, i3, i4, i5) = complex_log_from_phase( &
                                logabs(i1, i2, i3, i4, i5), phase(i1, i2, i3, i4, i5))
                        end do
                    end do
                end do
            end do
        end do
    end function logsumexp_keepdims_complex_5d_weighted

    pure subroutine logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the weighted sum magnitude
        complex(dp), intent(out) :: phase !! unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        logical :: ok
        complex(dp), allocatable :: terms(:, :, :, :, :)
        complex(dp), allocatable :: weights(:, :, :, :, :)

        call broadcast_complex_terms_5d(x, b, terms, weights, ok)
        if (.not. ok) then
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
            return
        end if
        call logsumexp_complex_core(reshape(terms, [size(terms)]), logabs, phase, &
            reshape(weights, [size(weights)]))
    end subroutine logsumexp_signed_complex_5d_all_weighted

    pure function logsumexp_complex_5d_all_weighted(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        real(dp) :: logabs
        complex(dp) :: phase

        call logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_5d_all_weighted

    pure subroutine logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        real(dp), allocatable :: kept_logabs(:, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, :)

        if (dim < 1 .or. dim > 5) then
            allocate(logabs(0, 0, 0, 0), phase(0, 0, 0, 0))
            return
        end if
        call logsumexp_signed_keepdims_complex_5d_weighted( &
            x, [dim], kept_logabs, kept_phase, b)
        select case (dim)
        case (1)
            allocate(logabs(size(kept_logabs, 2), size(kept_logabs, 3), &
                size(kept_logabs, 4), size(kept_logabs, 5)))
            allocate(phase(size(kept_phase, 2), size(kept_phase, 3), &
                size(kept_phase, 4), size(kept_phase, 5)))
        case (2)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 3), &
                size(kept_logabs, 4), size(kept_logabs, 5)))
            allocate(phase(size(kept_phase, 1), size(kept_phase, 3), &
                size(kept_phase, 4), size(kept_phase, 5)))
        case (3)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 2), &
                size(kept_logabs, 4), size(kept_logabs, 5)))
            allocate(phase(size(kept_phase, 1), size(kept_phase, 2), &
                size(kept_phase, 4), size(kept_phase, 5)))
        case (4)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 2), &
                size(kept_logabs, 3), size(kept_logabs, 5)))
            allocate(phase(size(kept_phase, 1), size(kept_phase, 2), &
                size(kept_phase, 3), size(kept_phase, 5)))
        case (5)
            allocate(logabs(size(kept_logabs, 1), size(kept_logabs, 2), &
                size(kept_logabs, 3), size(kept_logabs, 4)))
            allocate(phase(size(kept_phase, 1), size(kept_phase, 2), &
                size(kept_phase, 3), size(kept_phase, 4)))
        end select
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_5d_dim_weighted

    pure function logsumexp_complex_5d_dim_weighted(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :, :)

        real(dp), allocatable :: logabs(:, :, :, :)
        complex(dp), allocatable :: phase(:, :, :, :)
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4

        call logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), size(logabs, 4)))
        do i4 = 1, size(values, 4)
            do i3 = 1, size(values, 3)
                do i2 = 1, size(values, 2)
                    do i1 = 1, size(values, 1)
                        values(i1, i2, i3, i4) = complex_log_from_phase( &
                            logabs(i1, i2, i3, i4), phase(i1, i2, i3, i4))
                    end do
                end do
            end do
        end do
    end function logsumexp_complex_5d_dim_weighted

    pure function logsumexp_real_5d_all_complex_b(x, b) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp) :: value

        value = logsumexp_complex_5d_all_weighted(cmplx(x, 0.0_dp, kind=dp), b)
    end function logsumexp_real_5d_all_complex_b

    pure function logsumexp_real_5d_dim_complex_b(x, dim, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :, :)

        values = logsumexp_complex_5d_dim_weighted(cmplx(x, 0.0_dp, kind=dp), dim, b)
    end function logsumexp_real_5d_dim_complex_b

    pure subroutine logsumexp_signed_real_5d_all_complex_b(x, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the weighted sum magnitude
        complex(dp), intent(out) :: phase !! unit phase of the weighted exponential sum
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_5d_all_weighted( &
            cmplx(x, 0.0_dp, kind=dp), logabs, phase, b)
    end subroutine logsumexp_signed_real_5d_all_complex_b

    pure subroutine logsumexp_signed_real_5d_dim_complex_b(x, dim, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_complex_5d_dim_weighted( &
            cmplx(x, 0.0_dp, kind=dp), dim, logabs, phase, b)
    end subroutine logsumexp_signed_real_5d_dim_complex_b

    pure function logsumexp_keepdims_real_5d_complex_b(x, dims, b) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights
        complex(dp), allocatable :: values(:, :, :, :, :)

        values = logsumexp_keepdims_complex_5d_weighted(cmplx(x, 0.0_dp, kind=dp), dims, b)
    end function logsumexp_keepdims_real_5d_complex_b

    pure subroutine logsumexp_signed_keepdims_real_5d_complex_b(x, dims, logabs, phase, b)
        real(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 real logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :) !! reduced unit phases
        complex(dp), intent(in) :: b(..) !! scalar or broadcastable complex exponential weights

        call logsumexp_signed_keepdims_complex_5d_weighted( &
            cmplx(x, 0.0_dp, kind=dp), dims, logabs, phase, b)
    end subroutine logsumexp_signed_keepdims_real_5d_complex_b

    pure function logsumexp_complex_5d_all_real_b(x, b) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp) :: value

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)
        complex(dp), allocatable :: b5(:, :, :, :, :)

        select rank (b)
        rank (0)
            value = logsumexp_complex_5d_all_weighted(x, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_5d_all_weighted(x, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_5d_all_weighted(x, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_5d_all_weighted(x, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_5d_all_weighted(x, b4)
        rank (5)
            allocate(b5(size(b, 1), size(b, 2), size(b, 3), size(b, 4), size(b, 5)))
            b5 = cmplx(b, 0.0_dp, kind=dp)
            value = logsumexp_complex_5d_all_weighted(x, b5)
        rank default
            value = complex_nan()
        end select
    end function logsumexp_complex_5d_all_real_b

    pure function logsumexp_complex_5d_dim_real_b(x, dim, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: values(:, :, :, :)

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)
        complex(dp), allocatable :: b5(:, :, :, :, :)

        select rank (b)
        rank (0)
            values = logsumexp_complex_5d_dim_weighted(x, dim, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_5d_dim_weighted(x, dim, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_5d_dim_weighted(x, dim, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_5d_dim_weighted(x, dim, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_5d_dim_weighted(x, dim, b4)
        rank (5)
            allocate(b5(size(b, 1), size(b, 2), size(b, 3), size(b, 4), size(b, 5)))
            b5 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_complex_5d_dim_weighted(x, dim, b5)
        rank default
            if (dim < 1 .or. dim > 5) then
                allocate(values(0, 0, 0, 0))
            else
                select case (dim)
                case (1)
                    allocate(values(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
                case (2)
                    allocate(values(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
                case (3)
                    allocate(values(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
                case (4)
                    allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
                case (5)
                    allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
                end select
                values = complex_nan()
            end if
        end select
    end function logsumexp_complex_5d_dim_real_b

    pure subroutine logsumexp_signed_complex_5d_all_real_b(x, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        real(dp), intent(out) :: logabs !! logarithm of the weighted sum magnitude
        complex(dp), intent(out) :: phase !! unit phase of the weighted exponential sum
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)
        complex(dp), allocatable :: b5(:, :, :, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_5d_all_weighted( &
                x, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b4)
        rank (5)
            allocate(b5(size(b, 1), size(b, 2), size(b, 3), size(b, 4), size(b, 5)))
            b5 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_all_weighted(x, logabs, phase, b5)
        rank default
            logabs = quiet_nan(0.0_dp)
            phase = complex_nan()
        end select
    end subroutine logsumexp_signed_complex_5d_all_real_b

    pure subroutine logsumexp_signed_complex_5d_dim_real_b(x, dim, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 5
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)
        complex(dp), allocatable :: b5(:, :, :, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_complex_5d_dim_weighted( &
                x, dim, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b4)
        rank (5)
            allocate(b5(size(b, 1), size(b, 2), size(b, 3), size(b, 4), size(b, 5)))
            b5 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_complex_5d_dim_weighted(x, dim, logabs, phase, b5)
        rank default
            if (dim < 1 .or. dim > 5) then
                allocate(logabs(0, 0, 0, 0), phase(0, 0, 0, 0))
            else
                select case (dim)
                case (1)
                    allocate(logabs(size(x, 2), size(x, 3), size(x, 4), size(x, 5)), &
                        phase(size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
                case (2)
                    allocate(logabs(size(x, 1), size(x, 3), size(x, 4), size(x, 5)), &
                        phase(size(x, 1), size(x, 3), size(x, 4), size(x, 5)))
                case (3)
                    allocate(logabs(size(x, 1), size(x, 2), size(x, 4), size(x, 5)), &
                        phase(size(x, 1), size(x, 2), size(x, 4), size(x, 5)))
                case (4)
                    allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 5)), &
                        phase(size(x, 1), size(x, 2), size(x, 3), size(x, 5)))
                case (5)
                    allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4)), &
                        phase(size(x, 1), size(x, 2), size(x, 3), size(x, 4)))
                end select
                logabs = quiet_nan(0.0_dp)
                phase = complex_nan()
            end if
        end select
    end subroutine logsumexp_signed_complex_5d_dim_real_b

    pure function logsumexp_keepdims_complex_5d_real_b(x, dims, b) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights
        complex(dp), allocatable :: values(:, :, :, :, :)

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)
        complex(dp), allocatable :: b5(:, :, :, :, :)

        select rank (b)
        rank (0)
            values = logsumexp_keepdims_complex_5d_weighted( &
                x, dims, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_5d_weighted(x, dims, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_5d_weighted(x, dims, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_5d_weighted(x, dims, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_5d_weighted(x, dims, b4)
        rank (5)
            allocate(b5(size(b, 1), size(b, 2), size(b, 3), size(b, 4), size(b, 5)))
            b5 = cmplx(b, 0.0_dp, kind=dp)
            values = logsumexp_keepdims_complex_5d_weighted(x, dims, b5)
        rank default
            call allocate_keepdims_complex_5d(values, size(x, 1), size(x, 2), &
                size(x, 3), size(x, 4), size(x, 5), reduction_dims_mode_5d(dims))
            if (size(values) > 0) values = complex_nan()
        end select
    end function logsumexp_keepdims_complex_5d_real_b

    pure subroutine logsumexp_signed_keepdims_complex_5d_real_b(x, dims, logabs, phase, b)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! rank-5 complex logarithmic terms
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :) !! reduced unit phases
        real(dp), intent(in) :: b(..) !! scalar or broadcastable real exponential weights

        complex(dp), allocatable :: b1(:)
        complex(dp), allocatable :: b2(:, :)
        complex(dp), allocatable :: b3(:, :, :)
        complex(dp), allocatable :: b4(:, :, :, :)
        complex(dp), allocatable :: b5(:, :, :, :, :)

        select rank (b)
        rank (0)
            call logsumexp_signed_keepdims_complex_5d_weighted( &
                x, dims, logabs, phase, cmplx(b, 0.0_dp, kind=dp))
        rank (1)
            allocate(b1(size(b)))
            b1 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b1)
        rank (2)
            allocate(b2(size(b, 1), size(b, 2)))
            b2 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b2)
        rank (3)
            allocate(b3(size(b, 1), size(b, 2), size(b, 3)))
            b3 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b3)
        rank (4)
            allocate(b4(size(b, 1), size(b, 2), size(b, 3), size(b, 4)))
            b4 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b4)
        rank (5)
            allocate(b5(size(b, 1), size(b, 2), size(b, 3), size(b, 4), size(b, 5)))
            b5 = cmplx(b, 0.0_dp, kind=dp)
            call logsumexp_signed_keepdims_complex_5d_weighted(x, dims, logabs, phase, b5)
        rank default
            call allocate_keepdims_real_5d(logabs, size(x, 1), size(x, 2), &
                size(x, 3), size(x, 4), size(x, 5), reduction_dims_mode_5d(dims))
            call allocate_keepdims_complex_5d(phase, size(x, 1), size(x, 2), &
                size(x, 3), size(x, 4), size(x, 5), reduction_dims_mode_5d(dims))
            if (size(logabs) > 0) logabs = quiet_nan(0.0_dp)
            if (size(phase) > 0) phase = complex_nan()
        end select
    end subroutine logsumexp_signed_keepdims_complex_5d_real_b

    pure function softmax_complex_5d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        if (present(dim)) then
            values = softmax_complex_5d_dims(x, [dim])
        else
            values = softmax_complex_5d_dims(x, [1, 2, 3, 4, 5])
        end if
    end function softmax_complex_5d

    pure function softmax_complex_5d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        call normalize_complex_5d_dims(x, dims, .false., values)
    end function softmax_complex_5d_dims

    pure function log_softmax_complex_5d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        if (present(dim)) then
            values = log_softmax_complex_5d_dims(x, [dim])
        else
            values = log_softmax_complex_5d_dims(x, [1, 2, 3, 4, 5])
        end if
    end function log_softmax_complex_5d

    pure function log_softmax_complex_5d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5))

        call normalize_complex_5d_dims(x, dims, .true., values)
    end function log_softmax_complex_5d_dims


    pure subroutine normalize_real_5d_dims(x, dims, logarithmic, values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        logical, intent(in) :: logarithmic !! return log-softmax when true, softmax otherwise
        real(dp), intent(out) :: values(:, :, :, :, :) !! normalized values with the input shape

        real(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: p
        integer :: reduction_size
        real(dp), allocatable :: reduced(:)

        mode = reduction_dims_mode_5d(dims)
        if (mode < 0) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (mode == 0) then
            if (logarithmic) then
                where (ieee_is_finite(x))
                    values = 0.0_dp
                elsewhere
                    values = quiet_nan(0.0_dp)
                end where
            else
                where (ieee_is_finite(x))
                    values = 1.0_dp
                elsewhere
                    values = quiet_nan(0.0_dp)
                end where
            end if
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        allocate(buffer(max(1, reduction_size)))

        do o5 = 1, m5
            lo5 = merge(1, o5, btest(mode, 4))
            hi5 = merge(n5, o5, btest(mode, 4))
            do o4 = 1, m4
                lo4 = merge(1, o4, btest(mode, 3))
                hi4 = merge(n4, o4, btest(mode, 3))
                do o3 = 1, m3
                    lo3 = merge(1, o3, btest(mode, 2))
                    hi3 = merge(n3, o3, btest(mode, 2))
                    do o2 = 1, m2
                        lo2 = merge(1, o2, btest(mode, 1))
                        hi2 = merge(n2, o2, btest(mode, 1))
                        do o1 = 1, m1
                            lo1 = merge(1, o1, btest(mode, 0))
                            hi1 = merge(n1, o1, btest(mode, 0))
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                buffer(p) = x(i1, i2, i3, i4, i5)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            if (logarithmic) then
                                reduced = log_softmax_1d(buffer(1:p))
                            else
                                reduced = softmax_1d(buffer(1:p))
                            end if
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                values(i1, i2, i3, i4, i5) = reduced(p)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine normalize_real_5d_dims

    pure subroutine normalize_complex_5d_dims(x, dims, logarithmic, values)
        complex(dp), intent(in) :: x(:, :, :, :, :) !! five-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        logical, intent(in) :: logarithmic !! return log-softmax when true, softmax otherwise
        complex(dp), intent(out) :: values(:, :, :, :, :) !! normalized values with input shape

        complex(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: p
        integer :: reduction_size
        complex(dp), allocatable :: reduced(:)

        mode = reduction_dims_mode_5d(dims)
        if (mode < 0) then
            values = complex_nan()
            return
        end if
        if (mode == 0) then
            if (logarithmic) then
                where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                    values = cmplx(0.0_dp, 0.0_dp, kind=dp)
                elsewhere
                    values = complex_nan()
                end where
            else
                where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                    values = cmplx(1.0_dp, 0.0_dp, kind=dp)
                elsewhere
                    values = complex_nan()
                end where
            end if
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        allocate(buffer(max(1, reduction_size)))

        do o5 = 1, m5
            lo5 = merge(1, o5, btest(mode, 4))
            hi5 = merge(n5, o5, btest(mode, 4))
            do o4 = 1, m4
                lo4 = merge(1, o4, btest(mode, 3))
                hi4 = merge(n4, o4, btest(mode, 3))
                do o3 = 1, m3
                    lo3 = merge(1, o3, btest(mode, 2))
                    hi3 = merge(n3, o3, btest(mode, 2))
                    do o2 = 1, m2
                        lo2 = merge(1, o2, btest(mode, 1))
                        hi2 = merge(n2, o2, btest(mode, 1))
                        do o1 = 1, m1
                            lo1 = merge(1, o1, btest(mode, 0))
                            hi1 = merge(n1, o1, btest(mode, 0))
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                buffer(p) = x(i1, i2, i3, i4, i5)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            if (logarithmic) then
                                reduced = log_softmax_complex_1d(buffer(1:p))
                            else
                                reduced = softmax_complex_1d(buffer(1:p))
                            end if
                            p = 0
                            do i5 = lo5, hi5
                                do i4 = lo4, hi4
                                    do i3 = lo3, hi3
                                        do i2 = lo2, hi2
                                            do i1 = lo1, hi1
                                                p = p + 1
                                                values(i1, i2, i3, i4, i5) = reduced(p)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine normalize_complex_5d_dims

    pure function reduction_dims_mode_6d(dims) result(mode)
        integer, intent(in) :: dims(:) !! Fortran dimensions to reduce
        integer :: mode

        integer :: i
        logical :: seen(6)

        mode = 0
        seen = .false.
        if (size(dims) > 6) then
            mode = -1
            return
        end if
        do i = 1, size(dims)
            if (dims(i) < 1 .or. dims(i) > 6) then
                mode = -1
                return
            end if
            if (seen(dims(i))) then
                mode = -1
                return
            end if
            seen(dims(i)) = .true.
            mode = mode + 2 ** (dims(i) - 1)
        end do
    end function reduction_dims_mode_6d

    pure subroutine logsumexp_signed_keepdims_6d(x, dims, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :) !! log-absolute sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :, :) !! signs of reduced sums

        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: hi6
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: i6
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: lo6
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: n6
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: o6
        integer :: p
        integer :: reduction_size
        real(dp), allocatable :: buffer(:)

        mode = reduction_dims_mode_6d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0, 0, 0), sign_value(0, 0, 0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        n6 = size(x, 6)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        m6 = merge(1, n6, btest(mode, 5))
        allocate(logabs(m1, m2, m3, m4, m5, m6), sign_value(m1, m2, m3, m4, m5, m6))

        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        if (btest(mode, 5)) reduction_size = reduction_size * n6
        allocate(buffer(max(1, reduction_size)))

        do o6 = 1, m6
            lo6 = merge(1, o6, btest(mode, 5))
            hi6 = merge(n6, o6, btest(mode, 5))
            do o5 = 1, m5
                lo5 = merge(1, o5, btest(mode, 4))
                hi5 = merge(n5, o5, btest(mode, 4))
                do o4 = 1, m4
                    lo4 = merge(1, o4, btest(mode, 3))
                    hi4 = merge(n4, o4, btest(mode, 3))
                    do o3 = 1, m3
                        lo3 = merge(1, o3, btest(mode, 2))
                        hi3 = merge(n3, o3, btest(mode, 2))
                        do o2 = 1, m2
                            lo2 = merge(1, o2, btest(mode, 1))
                            hi2 = merge(n2, o2, btest(mode, 1))
                            do o1 = 1, m1
                                lo1 = merge(1, o1, btest(mode, 0))
                                hi1 = merge(n1, o1, btest(mode, 0))
                                p = 0
                                do i6 = lo6, hi6
                                    do i5 = lo5, hi5
                                        do i4 = lo4, hi4
                                            do i3 = lo3, hi3
                                                do i2 = lo2, hi2
                                                    do i1 = lo1, hi1
                                                        p = p + 1
                                                        buffer(p) = x(i1, i2, i3, i4, i5, i6)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                                call logsumexp_core(buffer(1:p), logabs(o1, o2, o3, o4, o5, o6), &
                                    sign_value(o1, o2, o3, o4, o5, o6))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine logsumexp_signed_keepdims_6d

    pure function logsumexp_keepdims_6d(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable :: values(:, :, :, :, :, :)

        real(dp), allocatable :: signs(:, :, :, :, :, :)

        call logsumexp_signed_keepdims_6d(x, dims, values, signs)
        if (size(values) == 0) return
        where (ieee_is_nan(signs) .or. signs < 0.0_dp)
            values = quiet_nan(0.0_dp)
        end where
    end function logsumexp_keepdims_6d

    pure function logsumexp_6d_all(x) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! rank-6 values reduced over all dimensions
        real(dp) :: value

        real(dp) :: sign_value

        call logsumexp_core(reshape(x, [size(x)]), value, sign_value)
        if (ieee_is_nan(sign_value) .or. sign_value < 0.0_dp) value = quiet_nan(0.0_dp)
    end function logsumexp_6d_all

    pure function logsumexp_6d_dim(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable :: values(:, :, :, :, :)

        real(dp), allocatable :: kept(:, :, :, :, :, :)

        if (dim < 1 .or. dim > 6) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4), size(x, 5), size(x, 6)))
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 5), size(x, 6)))
        case (5)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 6)))
        case (6)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        end select
        kept = logsumexp_keepdims_6d(x, [dim])
        values = reshape(kept, shape(values))
    end function logsumexp_6d_dim

    pure subroutine logsumexp_signed_6d_all(x, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! rank-6 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the absolute exponential sum
        real(dp), intent(out) :: sign_value !! sign of the exponential sum: zero or +1

        call logsumexp_core(reshape(x, [size(x)]), logabs, sign_value)
    end subroutine logsumexp_signed_6d_all

    pure subroutine logsumexp_signed_6d_dim(x, dim, logabs, sign_value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! log-absolute reduced sums
        real(dp), allocatable, intent(out) :: sign_value(:, :, :, :, :) !! signs of reduced sums

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, :)
        real(dp), allocatable :: kept_sign(:, :, :, :, :, :)

        if (dim < 1 .or. dim > 6) then
            allocate(logabs(0, 0, 0, 0, 0), sign_value(0, 0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(logabs(size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6)), &
                sign_value(size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (2)
            allocate(logabs(size(x, 1), size(x, 3), size(x, 4), size(x, 5), size(x, 6)), &
                sign_value(size(x, 1), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (3)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 4), size(x, 5), size(x, 6)), &
                sign_value(size(x, 1), size(x, 2), size(x, 4), size(x, 5), size(x, 6)))
        case (4)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 5), size(x, 6)), &
                sign_value(size(x, 1), size(x, 2), size(x, 3), size(x, 5), size(x, 6)))
        case (5)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 6)), &
                sign_value(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 6)))
        case (6)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)), &
                sign_value(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        end select
        call logsumexp_signed_keepdims_6d(x, [dim], kept_logabs, kept_sign)
        logabs = reshape(kept_logabs, shape(logabs))
        sign_value = reshape(kept_sign, shape(sign_value))
    end subroutine logsumexp_signed_6d_dim

    pure function softmax_6d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6))

        if (present(dim)) then
            values = softmax_6d_dims(x, [dim])
        else
            values = softmax_6d_dims(x, [1, 2, 3, 4, 5, 6])
        end if
    end function softmax_6d

    pure function softmax_6d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6))

        call normalize_real_6d_dims(x, dims, .false., values)
    end function softmax_6d_dims

    pure function log_softmax_6d(x, dim) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6))

        if (present(dim)) then
            values = log_softmax_6d_dims(x, [dim])
        else
            values = log_softmax_6d_dims(x, [1, 2, 3, 4, 5, 6])
        end if
    end function log_softmax_6d

    pure function log_softmax_6d_dims(x, dims) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        real(dp) :: values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6))

        call normalize_real_6d_dims(x, dims, .true., values)
    end function log_softmax_6d_dims

    pure subroutine logsumexp_signed_keepdims_complex_6d(x, dims, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :, :) !! log-magnitudes of sums
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :, :) !! unit phases of sums

        complex(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: hi6
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: i6
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: lo6
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: n6
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: o6
        integer :: p
        integer :: reduction_size

        mode = reduction_dims_mode_6d(dims)
        if (mode < 0) then
            allocate(logabs(0, 0, 0, 0, 0, 0), phase(0, 0, 0, 0, 0, 0))
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        n6 = size(x, 6)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        m6 = merge(1, n6, btest(mode, 5))
        allocate(logabs(m1, m2, m3, m4, m5, m6), phase(m1, m2, m3, m4, m5, m6))

        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        if (btest(mode, 5)) reduction_size = reduction_size * n6
        allocate(buffer(max(1, reduction_size)))

        do o6 = 1, m6
            lo6 = merge(1, o6, btest(mode, 5))
            hi6 = merge(n6, o6, btest(mode, 5))
            do o5 = 1, m5
                lo5 = merge(1, o5, btest(mode, 4))
                hi5 = merge(n5, o5, btest(mode, 4))
                do o4 = 1, m4
                    lo4 = merge(1, o4, btest(mode, 3))
                    hi4 = merge(n4, o4, btest(mode, 3))
                    do o3 = 1, m3
                        lo3 = merge(1, o3, btest(mode, 2))
                        hi3 = merge(n3, o3, btest(mode, 2))
                        do o2 = 1, m2
                            lo2 = merge(1, o2, btest(mode, 1))
                            hi2 = merge(n2, o2, btest(mode, 1))
                            do o1 = 1, m1
                                lo1 = merge(1, o1, btest(mode, 0))
                                hi1 = merge(n1, o1, btest(mode, 0))
                                p = 0
                                do i6 = lo6, hi6
                                    do i5 = lo5, hi5
                                        do i4 = lo4, hi4
                                            do i3 = lo3, hi3
                                                do i2 = lo2, hi2
                                                    do i1 = lo1, hi1
                                                        p = p + 1
                                                        buffer(p) = x(i1, i2, i3, i4, i5, i6)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                                call logsumexp_complex_core(buffer(1:p), &
                                    logabs(o1, o2, o3, o4, o5, o6), &
                                    phase(o1, o2, o3, o4, o5, o6))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine logsumexp_signed_keepdims_complex_6d

    pure function logsumexp_keepdims_complex_6d(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dims(:) !! reduced dimensions retained with extent one
        complex(dp), allocatable :: values(:, :, :, :, :, :)

        complex(dp), allocatable :: phase(:, :, :, :, :, :)
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: i6
        real(dp), allocatable :: logabs(:, :, :, :, :, :)

        call logsumexp_signed_keepdims_complex_6d(x, dims, logabs, phase)
        allocate(values(size(logabs, 1), size(logabs, 2), size(logabs, 3), &
            size(logabs, 4), size(logabs, 5), size(logabs, 6)))
        do i6 = 1, size(values, 6)
            do i5 = 1, size(values, 5)
                do i4 = 1, size(values, 4)
                    do i3 = 1, size(values, 3)
                        do i2 = 1, size(values, 2)
                            do i1 = 1, size(values, 1)
                                values(i1, i2, i3, i4, i5, i6) = complex_log_from_phase( &
                                    logabs(i1, i2, i3, i4, i5, i6), phase(i1, i2, i3, i4, i5, i6))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end function logsumexp_keepdims_complex_6d

    pure function logsumexp_complex_6d_all(x) result(value)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! rank-6 values reduced over all dimensions
        complex(dp) :: value

        complex(dp) :: phase
        real(dp) :: logabs

        call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
        value = complex_log_from_phase(logabs, phase)
    end function logsumexp_complex_6d_all

    pure function logsumexp_complex_6d_dim(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        complex(dp), allocatable :: values(:, :, :, :, :)

        complex(dp), allocatable :: kept(:, :, :, :, :, :)

        if (dim < 1 .or. dim > 6) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(values(size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (2)
            allocate(values(size(x, 1), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (3)
            allocate(values(size(x, 1), size(x, 2), size(x, 4), size(x, 5), size(x, 6)))
        case (4)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 5), size(x, 6)))
        case (5)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 6)))
        case (6)
            allocate(values(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        end select
        kept = logsumexp_keepdims_complex_6d(x, [dim])
        values = reshape(kept, shape(values))
    end function logsumexp_complex_6d_dim

    pure subroutine logsumexp_signed_complex_6d_all(x, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! rank-6 values reduced over all dimensions
        real(dp), intent(out) :: logabs !! logarithm of the exponential-sum magnitude
        complex(dp), intent(out) :: phase !! unit complex phase of the exponential sum

        call logsumexp_complex_core(reshape(x, [size(x)]), logabs, phase)
    end subroutine logsumexp_signed_complex_6d_all

    pure subroutine logsumexp_signed_complex_6d_dim(x, dim, logabs, phase)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to reduce
        integer, intent(in) :: dim !! Fortran dimension to reduce, from 1 through 6
        real(dp), allocatable, intent(out) :: logabs(:, :, :, :, :) !! reduced log-magnitudes
        complex(dp), allocatable, intent(out) :: phase(:, :, :, :, :) !! reduced unit phases

        real(dp), allocatable :: kept_logabs(:, :, :, :, :, :)
        complex(dp), allocatable :: kept_phase(:, :, :, :, :, :)

        if (dim < 1 .or. dim > 6) then
            allocate(logabs(0, 0, 0, 0, 0), phase(0, 0, 0, 0, 0))
            return
        end if
        select case (dim)
        case (1)
            allocate(logabs(size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6)), &
                phase(size(x, 2), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (2)
            allocate(logabs(size(x, 1), size(x, 3), size(x, 4), size(x, 5), size(x, 6)), &
                phase(size(x, 1), size(x, 3), size(x, 4), size(x, 5), size(x, 6)))
        case (3)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 4), size(x, 5), size(x, 6)), &
                phase(size(x, 1), size(x, 2), size(x, 4), size(x, 5), size(x, 6)))
        case (4)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 5), size(x, 6)), &
                phase(size(x, 1), size(x, 2), size(x, 3), size(x, 5), size(x, 6)))
        case (5)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 6)), &
                phase(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 6)))
        case (6)
            allocate(logabs(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)), &
                phase(size(x, 1), size(x, 2), size(x, 3), size(x, 4), size(x, 5)))
        end select
        call logsumexp_signed_keepdims_complex_6d(x, [dim], kept_logabs, kept_phase)
        logabs = reshape(kept_logabs, shape(logabs))
        phase = reshape(kept_phase, shape(phase))
    end subroutine logsumexp_signed_complex_6d_dim

    pure function softmax_complex_6d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), size(x, 5), size(x, 6))

        if (present(dim)) then
            values = softmax_complex_6d_dims(x, [dim])
        else
            values = softmax_complex_6d_dims(x, [1, 2, 3, 4, 5, 6])
        end if
    end function softmax_complex_6d

    pure function softmax_complex_6d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), size(x, 5), size(x, 6))

        call normalize_complex_6d_dims(x, dims, .false., values)
    end function softmax_complex_6d_dims

    pure function log_softmax_complex_6d(x, dim) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to log-normalize
        integer, intent(in), optional :: dim !! optional Fortran dimension to normalize
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), size(x, 5), size(x, 6))

        if (present(dim)) then
            values = log_softmax_complex_6d_dims(x, [dim])
        else
            values = log_softmax_complex_6d_dims(x, [1, 2, 3, 4, 5, 6])
        end if
    end function log_softmax_complex_6d

    pure function log_softmax_complex_6d_dims(x, dims) result(values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to log-normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        complex(dp) :: values(size(x, 1), size(x, 2), size(x, 3), &
            size(x, 4), size(x, 5), size(x, 6))

        call normalize_complex_6d_dims(x, dims, .true., values)
    end function log_softmax_complex_6d_dims

    pure subroutine normalize_real_6d_dims(x, dims, logarithmic, values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        logical, intent(in) :: logarithmic !! return log-softmax when true, softmax otherwise
        real(dp), intent(out) :: values(:, :, :, :, :, :) !! normalized values with input shape

        real(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: hi6
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: i6
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: lo6
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: n6
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: o6
        integer :: p
        integer :: reduction_size
        real(dp), allocatable :: reduced(:)

        mode = reduction_dims_mode_6d(dims)
        if (mode < 0) then
            values = quiet_nan(0.0_dp)
            return
        end if
        if (mode == 0) then
            if (logarithmic) then
                where (ieee_is_finite(x))
                    values = 0.0_dp
                elsewhere
                    values = quiet_nan(0.0_dp)
                end where
            else
                where (ieee_is_finite(x))
                    values = 1.0_dp
                elsewhere
                    values = quiet_nan(0.0_dp)
                end where
            end if
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        n6 = size(x, 6)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        m6 = merge(1, n6, btest(mode, 5))
        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        if (btest(mode, 5)) reduction_size = reduction_size * n6
        allocate(buffer(max(1, reduction_size)))

        do o6 = 1, m6
            lo6 = merge(1, o6, btest(mode, 5))
            hi6 = merge(n6, o6, btest(mode, 5))
            do o5 = 1, m5
                lo5 = merge(1, o5, btest(mode, 4))
                hi5 = merge(n5, o5, btest(mode, 4))
                do o4 = 1, m4
                    lo4 = merge(1, o4, btest(mode, 3))
                    hi4 = merge(n4, o4, btest(mode, 3))
                    do o3 = 1, m3
                        lo3 = merge(1, o3, btest(mode, 2))
                        hi3 = merge(n3, o3, btest(mode, 2))
                        do o2 = 1, m2
                            lo2 = merge(1, o2, btest(mode, 1))
                            hi2 = merge(n2, o2, btest(mode, 1))
                            do o1 = 1, m1
                                lo1 = merge(1, o1, btest(mode, 0))
                                hi1 = merge(n1, o1, btest(mode, 0))
                                p = 0
                                do i6 = lo6, hi6
                                    do i5 = lo5, hi5
                                        do i4 = lo4, hi4
                                            do i3 = lo3, hi3
                                                do i2 = lo2, hi2
                                                    do i1 = lo1, hi1
                                                        p = p + 1
                                                        buffer(p) = x(i1, i2, i3, i4, i5, i6)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                                if (logarithmic) then
                                    reduced = log_softmax_1d(buffer(1:p))
                                else
                                    reduced = softmax_1d(buffer(1:p))
                                end if
                                p = 0
                                do i6 = lo6, hi6
                                    do i5 = lo5, hi5
                                        do i4 = lo4, hi4
                                            do i3 = lo3, hi3
                                                do i2 = lo2, hi2
                                                    do i1 = lo1, hi1
                                                        p = p + 1
                                                        values(i1, i2, i3, i4, i5, i6) = reduced(p)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine normalize_real_6d_dims

    pure subroutine normalize_complex_6d_dims(x, dims, logarithmic, values)
        complex(dp), intent(in) :: x(:, :, :, :, :, :) !! six-dimensional values to normalize
        integer, intent(in) :: dims(:) !! Fortran dimensions normalized jointly
        logical, intent(in) :: logarithmic !! return log-softmax when true, softmax otherwise
        complex(dp), intent(out) :: values(:, :, :, :, :, :) !! normalized values with input shape

        complex(dp), allocatable :: buffer(:)
        integer :: hi1
        integer :: hi2
        integer :: hi3
        integer :: hi4
        integer :: hi5
        integer :: hi6
        integer :: i1
        integer :: i2
        integer :: i3
        integer :: i4
        integer :: i5
        integer :: i6
        integer :: lo1
        integer :: lo2
        integer :: lo3
        integer :: lo4
        integer :: lo5
        integer :: lo6
        integer :: m1
        integer :: m2
        integer :: m3
        integer :: m4
        integer :: m5
        integer :: m6
        integer :: mode
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: n4
        integer :: n5
        integer :: n6
        integer :: o1
        integer :: o2
        integer :: o3
        integer :: o4
        integer :: o5
        integer :: o6
        integer :: p
        integer :: reduction_size
        complex(dp), allocatable :: reduced(:)

        mode = reduction_dims_mode_6d(dims)
        if (mode < 0) then
            values = complex_nan()
            return
        end if
        if (mode == 0) then
            if (logarithmic) then
                where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                    values = cmplx(0.0_dp, 0.0_dp, kind=dp)
                elsewhere
                    values = complex_nan()
                end where
            else
                where (ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x)))
                    values = cmplx(1.0_dp, 0.0_dp, kind=dp)
                elsewhere
                    values = complex_nan()
                end where
            end if
            return
        end if

        n1 = size(x, 1)
        n2 = size(x, 2)
        n3 = size(x, 3)
        n4 = size(x, 4)
        n5 = size(x, 5)
        n6 = size(x, 6)
        m1 = merge(1, n1, btest(mode, 0))
        m2 = merge(1, n2, btest(mode, 1))
        m3 = merge(1, n3, btest(mode, 2))
        m4 = merge(1, n4, btest(mode, 3))
        m5 = merge(1, n5, btest(mode, 4))
        m6 = merge(1, n6, btest(mode, 5))
        reduction_size = 1
        if (btest(mode, 0)) reduction_size = reduction_size * n1
        if (btest(mode, 1)) reduction_size = reduction_size * n2
        if (btest(mode, 2)) reduction_size = reduction_size * n3
        if (btest(mode, 3)) reduction_size = reduction_size * n4
        if (btest(mode, 4)) reduction_size = reduction_size * n5
        if (btest(mode, 5)) reduction_size = reduction_size * n6
        allocate(buffer(max(1, reduction_size)))

        do o6 = 1, m6
            lo6 = merge(1, o6, btest(mode, 5))
            hi6 = merge(n6, o6, btest(mode, 5))
            do o5 = 1, m5
                lo5 = merge(1, o5, btest(mode, 4))
                hi5 = merge(n5, o5, btest(mode, 4))
                do o4 = 1, m4
                    lo4 = merge(1, o4, btest(mode, 3))
                    hi4 = merge(n4, o4, btest(mode, 3))
                    do o3 = 1, m3
                        lo3 = merge(1, o3, btest(mode, 2))
                        hi3 = merge(n3, o3, btest(mode, 2))
                        do o2 = 1, m2
                            lo2 = merge(1, o2, btest(mode, 1))
                            hi2 = merge(n2, o2, btest(mode, 1))
                            do o1 = 1, m1
                                lo1 = merge(1, o1, btest(mode, 0))
                                hi1 = merge(n1, o1, btest(mode, 0))
                                p = 0
                                do i6 = lo6, hi6
                                    do i5 = lo5, hi5
                                        do i4 = lo4, hi4
                                            do i3 = lo3, hi3
                                                do i2 = lo2, hi2
                                                    do i1 = lo1, hi1
                                                        p = p + 1
                                                        buffer(p) = x(i1, i2, i3, i4, i5, i6)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                                if (logarithmic) then
                                    reduced = log_softmax_complex_1d(buffer(1:p))
                                else
                                    reduced = softmax_complex_1d(buffer(1:p))
                                end if
                                p = 0
                                do i6 = lo6, hi6
                                    do i5 = lo5, hi5
                                        do i4 = lo4, hi4
                                            do i3 = lo3, hi3
                                                do i2 = lo2, hi2
                                                    do i1 = lo1, hi1
                                                        p = p + 1
                                                        values(i1, i2, i3, i4, i5, i6) = reduced(p)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine normalize_complex_6d_dims

end module scifort_special_reductions_low
