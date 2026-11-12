! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Clarkson-Woodruff / CountSketch transform (scipy.linalg.clarkson_woodruff_transform).
! One random +/-1 entry per input row, with a randomly selected sketch row.
! The same immutable plan can be reused for several matrices and batches.
! Original implementation; no SciPy source has been copied.
module scifort_linalg_clarkson_woodruff
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    implicit none
    private

    ! Public components allow callers to use externally generated random hashes,
    ! including hashes produced with NumPy for exact cross-language checks.
    type, public :: countsketch_plan
        integer :: sketch_rows = 0 !! number of output rows
        integer, allocatable :: bucket(:) !! one-based output row for each input row
        integer, allocatable :: sign(:) !! +1 or -1 for each input row
    end type countsketch_plan

    public :: make_countsketch_plan, clarkson_woodruff_transform, clarkson_woodruff_csr

    interface clarkson_woodruff_transform
        module procedure sketch_real
        module procedure sketch_complex
        module procedure sketch_real_batch
        module procedure sketch_complex_batch
    end interface clarkson_woodruff_transform

    interface clarkson_woodruff_csr
        module procedure sketch_csr_real
        module procedure sketch_csr_complex
    end interface clarkson_woodruff_csr

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64

contains

    subroutine make_countsketch_plan(n_rows, sketch_rows, plan, status, seed)
        integer, intent(in) :: n_rows !! number of input rows, >= 0
        integer, intent(in) :: sketch_rows !! number of sketch rows, >= 1 unless n_rows is zero
        type(countsketch_plan), intent(out) :: plan !! generated hash bucket/sign plan
        integer, intent(out) :: status !! linalg_status_* result
        integer(int64), intent(in), optional :: seed !! local deterministic Park-Miller RNG seed
        integer(int64) :: rng_state
        real(dp) :: draw
        integer :: i

        status = linalg_status_invalid_shape
        if (n_rows < 0 .or. sketch_rows < 0) return
        if (n_rows > 0 .and. sketch_rows == 0) return
        plan%sketch_rows = sketch_rows
        allocate(plan%bucket(n_rows), plan%sign(n_rows))
        if (present(seed)) then
            ! Integer products remain below 2^63 even for maximal default-int
            ! sketch sizes. An explicit seed never touches Fortran's global RNG.
            rng_state = modulo(seed, rng_modulus - 1_int64) + 1_int64
            do i = 1, n_rows
                rng_state = modulo(rng_state * rng_multiplier, rng_modulus)
                plan%bucket(i) = 1 + int(((rng_state - 1_int64) * &
                    int(sketch_rows, int64)) / (rng_modulus - 1_int64))
                rng_state = modulo(rng_state * rng_multiplier, rng_modulus)
                if (rng_state <= (rng_modulus - 1_int64) / 2_int64) then
                    plan%sign(i) = 1
                else
                    plan%sign(i) = -1
                end if
            end do
        else
            do i = 1, n_rows
                call random_number(draw)
                plan%bucket(i) = 1 + min(sketch_rows - 1, int(draw * real(sketch_rows, dp)))
                call random_number(draw)
                plan%sign(i) = 1
                if (draw >= 0.5_dp) plan%sign(i) = -1
            end do
        end if
        status = linalg_status_success
    end subroutine make_countsketch_plan

    pure logical function plan_is_valid(plan, n_rows, sketch_rows) result(ok)
        type(countsketch_plan), intent(in) :: plan !! proposed sketch plan
        integer, intent(in) :: n_rows !! number of source rows
        integer, intent(in) :: sketch_rows !! number of result rows

        ok = .false.
        if (.not. allocated(plan%bucket) .or. .not. allocated(plan%sign)) return
        if (sketch_rows < 0 .or. n_rows < 0) return
        if (n_rows > 0 .and. sketch_rows < 1) return
        if (plan%sketch_rows /= sketch_rows) return
        if (size(plan%bucket) /= n_rows .or. size(plan%sign) /= n_rows) return
        if (any(plan%bucket < 1) .or. any(plan%bucket > sketch_rows)) return
        if (any(plan%sign /= 1 .and. plan%sign /= -1)) return
        ok = .true.
    end function plan_is_valid

    subroutine resolve_plan(n_rows, sketch_rows, supplied, seed, work, status)
        integer, intent(in) :: n_rows !! number of input rows
        integer, intent(in) :: sketch_rows !! requested sketch rows
        type(countsketch_plan), intent(in), optional :: supplied !! previously constructed plan
        integer(int64), intent(in), optional :: seed !! seed if supplied is omitted
        type(countsketch_plan), intent(out) :: work !! validated or generated plan
        integer, intent(out) :: status !! linalg_status_* result

        status = linalg_status_invalid_shape
        if (present(supplied)) then
            if (present(seed)) return
            if (.not. plan_is_valid(supplied, n_rows, sketch_rows)) return
            work = supplied
            status = linalg_status_success
        else
            call make_countsketch_plan(n_rows, sketch_rows, work, status, seed)
        end if
    end subroutine resolve_plan

    subroutine sketch_real(a, sketch_rows, result, status, plan, seed)
        real(dp), intent(in) :: a(:, :) !! n-by-d dense input matrix
        integer, intent(in) :: sketch_rows !! requested sketch dimension, >= 1 for n > 0
        real(dp), intent(out) :: result(:, :) !! caller-allocated sketch_rows-by-d result
        integer, intent(out) :: status !! linalg_status_* result
        type(countsketch_plan), intent(in), optional :: plan !! reusable hash/sign plan
        integer(int64), intent(in), optional :: seed !! local reproducibility seed if plan omitted
        type(countsketch_plan) :: work
        integer :: i, n, d

        result = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a, 1)
        d = size(a, 2)
        if (size(result, 1) /= sketch_rows .or. size(result, 2) /= d) return
        call resolve_plan(n, sketch_rows, plan, seed, work, status)
        if (status /= linalg_status_success) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        do i = 1, n
            result(work%bucket(i), :) = result(work%bucket(i), :) + real(work%sign(i), dp) * a(i, :)
        end do
        if (.not. all(ieee_is_finite(result))) then
            result = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine sketch_real

    subroutine sketch_complex(a, sketch_rows, result, status, plan, seed)
        complex(dp), intent(in) :: a(:, :) !! n-by-d complex input matrix
        integer, intent(in) :: sketch_rows !! requested sketch dimension, >= 1 for n > 0
        complex(dp), intent(out) :: result(:, :) !! caller-allocated sketch_rows-by-d result
        integer, intent(out) :: status !! linalg_status_* result
        type(countsketch_plan), intent(in), optional :: plan !! reusable hash/sign plan
        integer(int64), intent(in), optional :: seed !! local reproducibility seed if plan omitted
        type(countsketch_plan) :: work
        integer :: i, n, d

        result = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 1)
        d = size(a, 2)
        if (size(result, 1) /= sketch_rows .or. size(result, 2) /= d) return
        call resolve_plan(n, sketch_rows, plan, seed, work, status)
        if (status /= linalg_status_success) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(a,dp))) .or. .not. all(ieee_is_finite(aimag(a)))) return
        do i = 1, n
            result(work%bucket(i), :) = result(work%bucket(i), :) + real(work%sign(i), dp) * a(i, :)
        end do
        if (.not. all(ieee_is_finite(real(result,dp))) .or. &
            .not. all(ieee_is_finite(aimag(result)))) then
            result = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        status = linalg_status_success
    end subroutine sketch_complex

    subroutine sketch_real_batch(a, sketch_rows, result, status, plan, seed)
        real(dp), intent(in) :: a(:, :, :) !! batch-by-n-by-d input matrices
        integer, intent(in) :: sketch_rows !! requested sketch dimension
        real(dp), intent(out) :: result(:, :, :) !! batch-by-sketch_rows-by-d output
        integer, intent(out) :: status !! linalg_status_* result
        type(countsketch_plan), intent(in), optional :: plan !! shared sketch plan across batches
        integer(int64), intent(in), optional :: seed !! local reproducibility seed if plan omitted
        type(countsketch_plan) :: work
        integer :: i, b, n

        result = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a, 2)
        if (size(result, 1) /= size(a, 1) .or. size(result, 2) /= sketch_rows .or. &
            size(result, 3) /= size(a, 3)) return
        call resolve_plan(n, sketch_rows, plan, seed, work, status)
        if (status /= linalg_status_success) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        do b = 1, size(a, 1)
            do i = 1, n
                result(b, work%bucket(i), :) = result(b, work%bucket(i), :) + &
                    real(work%sign(i), dp) * a(b, i, :)
            end do
        end do
        if (.not. all(ieee_is_finite(result))) then
            result = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine sketch_real_batch

    subroutine sketch_complex_batch(a, sketch_rows, result, status, plan, seed)
        complex(dp), intent(in) :: a(:, :, :) !! batch-by-n-by-d complex matrices
        integer, intent(in) :: sketch_rows !! requested sketch dimension
        complex(dp), intent(out) :: result(:, :, :) !! batch-by-sketch_rows-by-d output
        integer, intent(out) :: status !! linalg_status_* result
        type(countsketch_plan), intent(in), optional :: plan !! shared sketch plan across batches
        integer(int64), intent(in), optional :: seed !! local reproducibility seed if plan omitted
        type(countsketch_plan) :: work
        integer :: i, b, n

        result = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 2)
        if (size(result, 1) /= size(a, 1) .or. size(result, 2) /= sketch_rows .or. &
            size(result, 3) /= size(a, 3)) return
        call resolve_plan(n, sketch_rows, plan, seed, work, status)
        if (status /= linalg_status_success) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(a,dp))) .or. .not. all(ieee_is_finite(aimag(a)))) return
        do b = 1, size(a, 1)
            do i = 1, n
                result(b, work%bucket(i), :) = result(b, work%bucket(i), :) + &
                    real(work%sign(i), dp) * a(b, i, :)
            end do
        end do
        if (.not. all(ieee_is_finite(real(result,dp))) .or. &
            .not. all(ieee_is_finite(aimag(result)))) then
            result = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        status = linalg_status_success
    end subroutine sketch_complex_batch

    pure logical function valid_csr(row_ptr, col_ind, nnz, ncols) result(ok)
        integer, intent(in) :: row_ptr(:) !! one-based CSR offsets, length nrows+1
        integer, intent(in) :: col_ind(:) !! one-based CSR column indices
        integer, intent(in) :: nnz !! number of stored entries
        integer, intent(in) :: ncols !! number of input columns

        ok = .false.
        if (ncols < 0 .or. size(row_ptr) < 1) return
        if (size(col_ind) /= nnz) return
        if (row_ptr(1) /= 1 .or. row_ptr(size(row_ptr)) /= nnz + 1) return
        if (any(row_ptr < 1) .or. any(row_ptr > nnz + 1)) return
        if (any(row_ptr(2:) < row_ptr(:size(row_ptr)-1))) return
        if (any(col_ind < 1) .or. any(col_ind > ncols)) return
        ok = .true.
    end function valid_csr

    subroutine sketch_csr_real(row_ptr, col_ind, values, ncols, sketch_rows, result, status, plan, seed)
        integer, intent(in) :: row_ptr(:) !! one-based CSR row offsets, length nrows+1
        integer, intent(in) :: col_ind(:) !! one-based columns, length nnz
        real(dp), intent(in) :: values(:) !! finite values, length nnz
        integer, intent(in) :: ncols !! number of columns of sparse input
        integer, intent(in) :: sketch_rows !! number of sketch rows
        real(dp), intent(out) :: result(:, :) !! sketch_rows-by-ncols result
        integer, intent(out) :: status !! linalg_status_* result
        type(countsketch_plan), intent(in), optional :: plan !! input-row mapping
        integer(int64), intent(in), optional :: seed !! local RNG seed if plan omitted
        type(countsketch_plan) :: work
        integer :: i, j, n

        result = 0.0_dp
        status = linalg_status_invalid_shape
        if (.not. valid_csr(row_ptr,col_ind,size(values),ncols)) return
        if (size(result, 1) /= sketch_rows .or. size(result, 2) /= ncols) return
        n = size(row_ptr) - 1
        call resolve_plan(n, sketch_rows, plan, seed, work, status)
        if (status /= linalg_status_success) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(values))) return
        do i = 1, n
            do j = row_ptr(i), row_ptr(i+1)-1
                result(work%bucket(i), col_ind(j)) = result(work%bucket(i), col_ind(j)) + &
                    real(work%sign(i),dp) * values(j)
            end do
        end do
        if (.not. all(ieee_is_finite(result))) then
            result = 0.0_dp
            return
        end if
        status = linalg_status_success
    end subroutine sketch_csr_real

    subroutine sketch_csr_complex(row_ptr, col_ind, values, ncols, sketch_rows, result, status, plan, seed)
        integer, intent(in) :: row_ptr(:) !! one-based CSR row offsets, length nrows+1
        integer, intent(in) :: col_ind(:) !! one-based columns, length nnz
        complex(dp), intent(in) :: values(:) !! finite complex entries, length nnz
        integer, intent(in) :: ncols !! number of columns of sparse input
        integer, intent(in) :: sketch_rows !! number of sketch rows
        complex(dp), intent(out) :: result(:, :) !! sketch_rows-by-ncols result
        integer, intent(out) :: status !! linalg_status_* result
        type(countsketch_plan), intent(in), optional :: plan !! input-row mapping
        integer(int64), intent(in), optional :: seed !! local RNG seed if plan omitted
        type(countsketch_plan) :: work
        integer :: i, j, n

        result = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (.not. valid_csr(row_ptr,col_ind,size(values),ncols)) return
        if (size(result, 1) /= sketch_rows .or. size(result, 2) /= ncols) return
        n = size(row_ptr) - 1
        call resolve_plan(n, sketch_rows, plan, seed, work, status)
        if (status /= linalg_status_success) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(values,dp))) .or. &
            .not. all(ieee_is_finite(aimag(values)))) return
        do i = 1, n
            do j = row_ptr(i), row_ptr(i+1)-1
                result(work%bucket(i), col_ind(j)) = result(work%bucket(i), col_ind(j)) + &
                    real(work%sign(i),dp) * values(j)
            end do
        end do
        if (.not. all(ieee_is_finite(real(result,dp))) .or. &
            .not. all(ieee_is_finite(aimag(result)))) then
            result = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        status = linalg_status_success
    end subroutine sketch_csr_complex

end module scifort_linalg_clarkson_woodruff
