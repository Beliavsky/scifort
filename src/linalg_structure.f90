! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Scalar and single-batch structural matrix inspection.
! Independent implementation of the defining diagonal and symmetry tests.
module scifort_linalg_structure
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private
    public :: bandwidth, issymmetric, ishermitian

    interface bandwidth
        module procedure bandwidth_i32_2
        module procedure bandwidth_i32_3
        module procedure bandwidth_i64_2
        module procedure bandwidth_i64_3
        module procedure bandwidth_r_2
        module procedure bandwidth_r_3
        module procedure bandwidth_c_2
        module procedure bandwidth_c_3
    end interface bandwidth

    interface issymmetric
        module procedure issymmetric_i32_2
        module procedure issymmetric_i32_3
        module procedure issymmetric_i64_2
        module procedure issymmetric_i64_3
        module procedure issymmetric_r_2
        module procedure issymmetric_r_3
        module procedure issymmetric_c_2
        module procedure issymmetric_c_3
    end interface issymmetric

    interface ishermitian
        module procedure ishermitian_i32_2
        module procedure ishermitian_i32_3
        module procedure ishermitian_i64_2
        module procedure ishermitian_i64_3
        module procedure ishermitian_r_2
        module procedure ishermitian_r_3
        module procedure ishermitian_c_2
        module procedure ishermitian_c_3
    end interface ishermitian

contains

    pure logical function valid_tolerances(atol, rtol) result(ok)
        real(dp), intent(in), optional :: atol !! absolute comparison tolerance, finite >=0
        real(dp), intent(in), optional :: rtol !! relative comparison tolerance, finite >=0

        ok = .true.
        if (present(atol)) then
            if (.not. ieee_is_finite(atol)) ok = .false.
            if (atol < 0.0_dp) ok = .false.
        end if
        if (present(rtol)) then
            if (.not. ieee_is_finite(rtol)) ok = .false.
            if (rtol < 0.0_dp) ok = .false.
        end if
    end function valid_tolerances

    pure logical function close_i32(a, b) result(ok)
        integer(int32), intent(in) :: a !! first entry
        integer(int32), intent(in) :: b !! reference entry

        ! SciPy 1.17 uses exact comparison for integer arrays, irrespective
        ! of atol/rtol. Preserve all int32 bits without promotion to float.
        ok = a == b
    end function close_i32

    pure logical function close_i64(a, b) result(ok)
        integer(int64), intent(in) :: a !! first entry
        integer(int64), intent(in) :: b !! reference entry

        ! No floating-point conversion: int64 comparisons are always exact.
        ok = a == b
    end function close_i64

    pure logical function close_r(a, b, at, rt, approximate) result(ok)
        real(dp), intent(in) :: a !! first entry
        real(dp), intent(in) :: b !! reference entry
        real(dp), intent(in) :: at !! nonnegative absolute tolerance
        real(dp), intent(in) :: rt !! nonnegative relative tolerance
        logical, intent(in) :: approximate !! whether nonzero tolerances were given
        real(dp) :: scale_value

        ok = .false.
        if (ieee_is_nan(a) .or. ieee_is_nan(b)) return
        if (a == b) then
            ok = .true.
            return
        end if
        if (.not. approximate) return
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) return
        ! Scaling prevents overflow in a-b when the input magnitudes are huge.
        scale_value = max(abs(a), abs(b), tiny(a))
        if (at >= 2.0_dp*min(scale_value, huge(at)/2.0_dp)) then
            ok = .true.
        else
            ok = abs(a/scale_value-b/scale_value) <= &
                at/scale_value + rt*(abs(b)/scale_value)
        end if
    end function close_r

    pure logical function close_c(a, b, at, rt, approximate) result(ok)
        complex(dp), intent(in) :: a !! first entry
        complex(dp), intent(in) :: b !! reference entry
        real(dp), intent(in) :: at !! nonnegative absolute tolerance
        real(dp), intent(in) :: rt !! nonnegative relative tolerance
        logical, intent(in) :: approximate !! whether nonzero tolerances were given
        real(dp) :: scale_value

        ok = .false.
        if (ieee_is_nan(real(a, dp)) .or. ieee_is_nan(aimag(a))) return
        if (ieee_is_nan(real(b, dp)) .or. ieee_is_nan(aimag(b))) return
        if (a == b) then
            ok = .true.
            return
        end if
        if (.not. approximate) return
        if (.not. ieee_is_finite(real(a, dp)) .or. .not. ieee_is_finite(aimag(a))) return
        if (.not. ieee_is_finite(real(b, dp)) .or. .not. ieee_is_finite(aimag(b))) return
        scale_value = max(abs(real(a, dp)), abs(aimag(a)), &
            abs(real(b, dp)), abs(aimag(b)), tiny(at))
        if (at >= 3.0_dp*min(scale_value, huge(at)/3.0_dp)) then
            ok = .true.
        else
            ok = abs(a/scale_value-b/scale_value) <= &
                at/scale_value + rt*abs(b/scale_value)
        end if
    end function close_c

    pure logical function inspect_i32(a) result(ok)
        integer(int32), intent(in) :: a(:, :) !! square input matrix
        integer :: i, j, n

        ok = .false.
        n = size(a, 1)
        ! Integer matrices have real diagonals; they need no diagonal scan.
        do j = 1, n
            do i = j + 1, n
                if (.not. close_i32(a(i, j), a(j, i))) return
            end do
        end do
        ok = .true.
    end function inspect_i32

    pure logical function inspect_i64(a) result(ok)
        integer(int64), intent(in) :: a(:, :) !! square input matrix
        integer :: i, j, n

        ok = .false.
        n = size(a, 1)
        ! Integer matrices have real diagonals; they need no diagonal scan.
        do j = 1, n
            do i = j + 1, n
                if (.not. close_i64(a(i, j), a(j, i))) return
            end do
        end do
        ok = .true.
    end function inspect_i64

    pure logical function inspect_r(a, at, rt, approximate) result(ok)
        real(dp), intent(in) :: a(:, :) !! square input matrix
        real(dp), intent(in) :: at !! nonnegative absolute tolerance
        real(dp), intent(in) :: rt !! nonnegative relative tolerance
        logical, intent(in) :: approximate !! true if any tolerance is positive
        integer :: i, j, n

        ok = .false.
        n = size(a, 1)
        ! Exact comparison ignores real diagonal entries (including NaN).
        ! Approximate comparison checks diagonal NaNs, as numpy.allclose does.
        if (approximate) then
            do i = 1, n
                if (.not. close_r(a(i, i), a(i, i), at, rt, approximate)) return
            end do
        end if
        do j = 1, n
            do i = j + 1, n
                if (.not. close_r(a(i, j), a(j, i), at, rt, approximate)) return
                if (approximate) then
                    if (.not. close_r(a(j, i), a(i, j), at, rt, approximate)) return
                end if
            end do
        end do
        ok = .true.
    end function inspect_r

    pure logical function inspect_c(a, hermitian, at, rt, approximate) result(ok)
        complex(dp), intent(in) :: a(:, :) !! square input matrix
        logical, intent(in) :: hermitian !! conjugate transpose instead of transpose
        real(dp), intent(in) :: at !! nonnegative absolute tolerance
        real(dp), intent(in) :: rt !! nonnegative relative tolerance
        logical, intent(in) :: approximate !! true if any tolerance is positive
        integer :: i, j, n

        ok = .false.
        n = size(a, 1)
        do i = 1, n
            if (hermitian) then
                if (.not. close_c(a(i, i), conjg(a(i, i)), at, rt, approximate)) return
            else if (approximate) then
                if (.not. close_c(a(i, i), a(i, i), at, rt, approximate)) return
            end if
        end do
        do j = 1, n
            do i = j + 1, n
                if (hermitian) then
                    if (.not. close_c(a(i, j), conjg(a(j, i)), at, rt, approximate)) return
                    if (approximate) then
                        if (.not. close_c(a(j, i), conjg(a(i, j)), at, rt, approximate)) return
                    end if
                else
                    if (.not. close_c(a(i, j), a(j, i), at, rt, approximate)) return
                    if (approximate) then
                        if (.not. close_c(a(j, i), a(i, j), at, rt, approximate)) return
                    end if
                end if
            end do
        end do
        ok = .true.
    end function inspect_c

    pure subroutine bandwidth_i32_2(a, lower, upper, status)
        integer(int32), intent(in) :: a(:, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower !! lower bandwidth(s), >=0
        integer, intent(out) :: upper !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        lower = 0
        upper = 0
        status = linalg_status_success
        ! SciPy returns (0,0) for any matrix with a singleton dimension.
        if (min(size(a, 1), size(a, 2)) <= 1) return
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (a(i, j) /= 0_int32) then
                    ! SciPy 1.17 excludes below-diagonal entries from its
                    ! last column, and above-diagonal entries from its last row.
                    ! The difference matters for rectangular matrices.
                    if (i > j .and. j == size(a, 2)) cycle
                    if (j > i .and. i == size(a, 1)) cycle
                    lower = max(lower, i-j)
                    upper = max(upper, j-i)
                end if
            end do
        end do
    end subroutine bandwidth_i32_2

    pure subroutine bandwidth_i32_3(a, lower, upper, status)
        integer(int32), intent(in) :: a(:, :, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower(:) !! lower bandwidth(s), >=0
        integer, intent(out) :: upper(:) !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: batch, b, status_one

        lower = 0
        upper = 0
        status = linalg_status_invalid_shape
        batch = size(a, 1)
        if (batch == 0) return
        if (size(lower) /= batch .or. size(upper) /= batch) return
        do b = 1, batch
            call bandwidth_i32_2(a(b, :, :), lower(b), upper(b), status_one)
        end do
        status = linalg_status_success
    end subroutine bandwidth_i32_3

    pure subroutine bandwidth_i64_2(a, lower, upper, status)
        integer(int64), intent(in) :: a(:, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower !! lower bandwidth(s), >=0
        integer, intent(out) :: upper !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        lower = 0
        upper = 0
        status = linalg_status_success
        ! SciPy returns (0,0) for any matrix with a singleton dimension.
        if (min(size(a, 1), size(a, 2)) <= 1) return
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (a(i, j) /= 0_int64) then
                    ! SciPy 1.17 excludes below-diagonal entries from its
                    ! last column, and above-diagonal entries from its last row.
                    ! The difference matters for rectangular matrices.
                    if (i > j .and. j == size(a, 2)) cycle
                    if (j > i .and. i == size(a, 1)) cycle
                    lower = max(lower, i-j)
                    upper = max(upper, j-i)
                end if
            end do
        end do
    end subroutine bandwidth_i64_2

    pure subroutine bandwidth_i64_3(a, lower, upper, status)
        integer(int64), intent(in) :: a(:, :, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower(:) !! lower bandwidth(s), >=0
        integer, intent(out) :: upper(:) !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: batch, b, status_one

        lower = 0
        upper = 0
        status = linalg_status_invalid_shape
        batch = size(a, 1)
        if (batch == 0) return
        if (size(lower) /= batch .or. size(upper) /= batch) return
        do b = 1, batch
            call bandwidth_i64_2(a(b, :, :), lower(b), upper(b), status_one)
        end do
        status = linalg_status_success
    end subroutine bandwidth_i64_3

    pure subroutine bandwidth_r_2(a, lower, upper, status)
        real(dp), intent(in) :: a(:, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower !! lower bandwidth(s), >=0
        integer, intent(out) :: upper !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        lower = 0
        upper = 0
        status = linalg_status_success
        ! SciPy returns (0,0) for any matrix with a singleton dimension.
        if (min(size(a, 1), size(a, 2)) <= 1) return
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (a(i, j) /= 0.0_dp) then
                    ! SciPy 1.17 excludes below-diagonal entries from its
                    ! last column, and above-diagonal entries from its last row.
                    ! The difference matters for rectangular matrices.
                    if (i > j .and. j == size(a, 2)) cycle
                    if (j > i .and. i == size(a, 1)) cycle
                    lower = max(lower, i-j)
                    upper = max(upper, j-i)
                end if
            end do
        end do
    end subroutine bandwidth_r_2

    pure subroutine bandwidth_r_3(a, lower, upper, status)
        real(dp), intent(in) :: a(:, :, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower(:) !! lower bandwidth(s), >=0
        integer, intent(out) :: upper(:) !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: batch, b, status_one

        lower = 0
        upper = 0
        status = linalg_status_invalid_shape
        batch = size(a, 1)
        if (batch == 0) return
        if (size(lower) /= batch .or. size(upper) /= batch) return
        do b = 1, batch
            call bandwidth_r_2(a(b, :, :), lower(b), upper(b), status_one)
        end do
        status = linalg_status_success
    end subroutine bandwidth_r_3

    pure subroutine bandwidth_c_2(a, lower, upper, status)
        complex(dp), intent(in) :: a(:, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower !! lower bandwidth(s), >=0
        integer, intent(out) :: upper !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: i, j

        lower = 0
        upper = 0
        status = linalg_status_success
        ! SciPy returns (0,0) for any matrix with a singleton dimension.
        if (min(size(a, 1), size(a, 2)) <= 1) return
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (a(i, j) /= cmplx(0.0_dp, 0.0_dp, dp)) then
                    ! SciPy 1.17 excludes below-diagonal entries from its
                    ! last column, and above-diagonal entries from its last row.
                    ! The difference matters for rectangular matrices.
                    if (i > j .and. j == size(a, 2)) cycle
                    if (j > i .and. i == size(a, 1)) cycle
                    lower = max(lower, i-j)
                    upper = max(upper, j-i)
                end if
            end do
        end do
    end subroutine bandwidth_c_2

    pure subroutine bandwidth_c_3(a, lower, upper, status)
        complex(dp), intent(in) :: a(:, :, :) !! matrix or one-batch matrix array
        integer, intent(out) :: lower(:) !! lower bandwidth(s), >=0
        integer, intent(out) :: upper(:) !! upper bandwidth(s), >=0
        integer, intent(out) :: status !! success or invalid_shape
        integer :: batch, b, status_one

        lower = 0
        upper = 0
        status = linalg_status_invalid_shape
        batch = size(a, 1)
        if (batch == 0) return
        if (size(lower) /= batch .or. size(upper) /= batch) return
        do b = 1, batch
            call bandwidth_c_2(a(b, :, :), lower(b), upper(b), status_one)
        end do
        status = linalg_status_success
    end subroutine bandwidth_c_3

    pure subroutine issymmetric_i32_2(a, result_value, status, atol, rtol)
        integer(int32), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        result_value = inspect_i32(a)
        status = linalg_status_success
    end subroutine issymmetric_i32_2

    pure subroutine issymmetric_i32_3(a, result_value, status, atol, rtol)
        integer(int32), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        do b = 1, batch
            result_value(b) = inspect_i32(a(b, :, :))
        end do
        status = linalg_status_success
    end subroutine issymmetric_i32_3

    pure subroutine issymmetric_i64_2(a, result_value, status, atol, rtol)
        integer(int64), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        result_value = inspect_i64(a)
        status = linalg_status_success
    end subroutine issymmetric_i64_2

    pure subroutine issymmetric_i64_3(a, result_value, status, atol, rtol)
        integer(int64), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        do b = 1, batch
            result_value(b) = inspect_i64(a(b, :, :))
        end do
        status = linalg_status_success
    end subroutine issymmetric_i64_3

    pure subroutine issymmetric_r_2(a, result_value, status, atol, rtol)
        real(dp), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        result_value = inspect_r(a, at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        status = linalg_status_success
    end subroutine issymmetric_r_2

    pure subroutine issymmetric_r_3(a, result_value, status, atol, rtol)
        real(dp), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        do b = 1, batch
            result_value(b) = inspect_r(a(b, :, :), at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        end do
        status = linalg_status_success
    end subroutine issymmetric_r_3

    pure subroutine issymmetric_c_2(a, result_value, status, atol, rtol)
        complex(dp), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        result_value = inspect_c(a, .false., at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        status = linalg_status_success
    end subroutine issymmetric_c_2

    pure subroutine issymmetric_c_3(a, result_value, status, atol, rtol)
        complex(dp), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        do b = 1, batch
            result_value(b) = inspect_c(a(b, :, :), .false., at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        end do
        status = linalg_status_success
    end subroutine issymmetric_c_3

    pure subroutine ishermitian_i32_2(a, result_value, status, atol, rtol)
        integer(int32), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        result_value = inspect_i32(a)
        status = linalg_status_success
    end subroutine ishermitian_i32_2

    pure subroutine ishermitian_i32_3(a, result_value, status, atol, rtol)
        integer(int32), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        do b = 1, batch
            result_value(b) = inspect_i32(a(b, :, :))
        end do
        status = linalg_status_success
    end subroutine ishermitian_i32_3

    pure subroutine ishermitian_i64_2(a, result_value, status, atol, rtol)
        integer(int64), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        result_value = inspect_i64(a)
        status = linalg_status_success
    end subroutine ishermitian_i64_2

    pure subroutine ishermitian_i64_3(a, result_value, status, atol, rtol)
        integer(int64), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        do b = 1, batch
            result_value(b) = inspect_i64(a(b, :, :))
        end do
        status = linalg_status_success
    end subroutine ishermitian_i64_3

    pure subroutine ishermitian_r_2(a, result_value, status, atol, rtol)
        real(dp), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        result_value = inspect_r(a, at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        status = linalg_status_success
    end subroutine ishermitian_r_2

    pure subroutine ishermitian_r_3(a, result_value, status, atol, rtol)
        real(dp), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        do b = 1, batch
            result_value(b) = inspect_r(a(b, :, :), at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        end do
        status = linalg_status_success
    end subroutine ishermitian_r_3

    pure subroutine ishermitian_c_2(a, result_value, status, atol, rtol)
        complex(dp), intent(in) :: a(:, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        if (size(a, 1) /= size(a, 2)) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        result_value = inspect_c(a, .true., at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        status = linalg_status_success
    end subroutine ishermitian_c_2

    pure subroutine ishermitian_c_3(a, result_value, status, atol, rtol)
        complex(dp), intent(in) :: a(:, :, :) !! square matrix or batch of square matrices
        logical, intent(out) :: result_value(:) !! symmetry/Hermitian result(s)
        integer, intent(out) :: status !! success or invalid_shape
        real(dp), intent(in), optional :: atol !! finite absolute tolerance, >=0
        real(dp), intent(in), optional :: rtol !! finite relative tolerance, >=0
        real(dp) :: at, rt
        integer :: b, batch

        result_value = .false.
        status = linalg_status_invalid_shape
        if (.not. valid_tolerances(atol, rtol)) return
        batch = size(a, 1)
        if (batch == 0) return
        if (size(a, 2) /= size(a, 3)) return
        if (size(result_value) /= batch) return
        at = 0.0_dp
        rt = 0.0_dp
        if (present(atol)) at = atol
        if (present(rtol)) rt = rtol
        do b = 1, batch
            result_value(b) = inspect_c(a(b, :, :), .true., at, rt, at > 0.0_dp .or. rt > 0.0_dp)
        end do
        status = linalg_status_success
    end subroutine ishermitian_c_3

end module scifort_linalg_structure
