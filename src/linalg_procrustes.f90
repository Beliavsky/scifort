! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Orthogonal/unitary Procrustes alignment via the SVD of A^H B.
! Independent implementation of the standard Schoenemann (1966) formula.
module scifort_linalg_procrustes
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_svd, only : svd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: orthogonal_procrustes

    interface orthogonal_procrustes
        module procedure procrustes_real
        module procedure procrustes_complex
        module procedure procrustes_real_complex
        module procedure procrustes_complex_real
    end interface orthogonal_procrustes

contains

    pure logical function finite_complex(a) result(valid)
        complex(dp), intent(in) :: a(:, :) !! matrix whose components are checked
        valid = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    pure subroutine scaled_product3(x, y, z, result_value, valid)
        real(dp), intent(in) :: x !! nonnegative finite multiplier
        real(dp), intent(in) :: y !! nonnegative finite multiplier
        real(dp), intent(in) :: z !! nonnegative finite multiplier
        real(dp), intent(out) :: result_value !! nonnegative product, if representable
        logical, intent(out) :: valid !! true when the product is representable
        real(dp) :: mantissa
        integer :: power

        result_value = 0.0_dp
        valid = .false.
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y) .or. &
            .not. ieee_is_finite(z)) return
        if (min(x, y, z) < 0.0_dp) return
        if (x == 0.0_dp .or. y == 0.0_dp .or. z == 0.0_dp) then
            valid = .true.
            return
        end if
        ! Combine normalized binary fractions before changing the exponent.
        ! This avoids premature over/underflow for e.g. 1e200 * 1e-200.
        mantissa = fraction(x)*fraction(y)*fraction(z)
        power = exponent(x)+exponent(y)+exponent(z)+exponent(mantissa)
        mantissa = fraction(mantissa)
        if (power > maxexponent(result_value)) return
        if (power == maxexponent(result_value)) then
            if (mantissa > fraction(huge(result_value))) return
        end if
        if (power < minexponent(result_value)-digits(result_value)) then
            valid = .true.
            return
        end if
        result_value = scale(mantissa, power)
        valid = ieee_is_finite(result_value)
    end subroutine scaled_product3

    subroutine procrustes_real(a, b, r, scale_factor, status, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! source matrix (m,n), finite or empty
        real(dp), intent(in) :: b(:, :) !! target matrix, exactly shape (m,n)
        real(dp), intent(out) :: r(:, :) !! optimal orthogonal factor, shape (n,n)
        real(dp), intent(out) :: scale_factor !! nuclear norm of transpose(A)*B
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! positive SVD Jacobi sweep limit; default 120
        real(dp), allocatable :: cross(:, :), u(:, :), vh(:, :), s(:)
        real(dp) :: amax, bmax
        logical :: valid
        integer :: m, n, i

        r = 0.0_dp
        scale_factor = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a, 1)
        n = size(a, 2)
        if (size(b, 1) /= m .or. size(b, 2) /= n) return
        if (size(r, 1) /= n .or. size(r, 2) /= n) return
        if (present(max_sweeps)) then
            if (max_sweeps < 1) return
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        if (m == 0 .or. n == 0) then
            do i = 1, n
                r(i, i) = 1.0_dp
            end do
            status = linalg_status_success
            return
        end if
        amax = maxval(abs(a))
        bmax = maxval(abs(b))
        if (amax == 0.0_dp .or. bmax == 0.0_dp) then
            do i = 1, n
                r(i, i) = 1.0_dp
            end do
            status = linalg_status_success
            return
        end if
        allocate(cross(n, n), u(n, n), vh(n, n), s(n))
        cross = matmul(transpose(a/amax), b/bmax)
        call svd(cross, u, s, vh, status, max_sweeps=max_sweeps)
        if (status /= linalg_status_success) return
        call scaled_product3(sum(s), amax, bmax, scale_factor, valid)
        if (.not. valid) then
            status = linalg_status_no_convergence
            return
        end if
        r = matmul(u, vh)
    end subroutine procrustes_real

    subroutine procrustes_complex(a, b, r, scale_factor, status, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! source matrix (m,n), finite or empty
        complex(dp), intent(in) :: b(:, :) !! target matrix, exactly shape (m,n)
        complex(dp), intent(out) :: r(:, :) !! optimal unitary factor, shape (n,n)
        real(dp), intent(out) :: scale_factor !! nuclear norm of conjugate-transpose(A)*B
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! positive SVD Jacobi sweep limit; default 120
        complex(dp), allocatable :: cross(:, :), u(:, :), vh(:, :)
        real(dp), allocatable :: s(:)
        real(dp) :: amax, bmax
        logical :: valid
        integer :: m, n, i

        r = cmplx(0.0_dp, 0.0_dp, dp)
        scale_factor = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(a, 1)
        n = size(a, 2)
        if (size(b, 1) /= m .or. size(b, 2) /= n) return
        if (size(r, 1) /= n .or. size(r, 2) /= n) return
        if (present(max_sweeps)) then
            if (max_sweeps < 1) return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        if (m == 0 .or. n == 0) then
            do i = 1, n
                r(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            status = linalg_status_success
            return
        end if
        ! Maximize over the components to avoid an overflowing complex abs().
        amax = max(maxval(abs(real(a, dp))), maxval(abs(aimag(a))))
        bmax = max(maxval(abs(real(b, dp))), maxval(abs(aimag(b))))
        if (amax == 0.0_dp .or. bmax == 0.0_dp) then
            do i = 1, n
                r(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            status = linalg_status_success
            return
        end if
        allocate(cross(n, n), u(n, n), vh(n, n), s(n))
        cross = matmul(transpose(conjg(a/amax)), b/bmax)
        call svd(cross, u, s, vh, status, max_sweeps=max_sweeps)
        if (status /= linalg_status_success) return
        call scaled_product3(sum(s), amax, bmax, scale_factor, valid)
        if (.not. valid) then
            status = linalg_status_no_convergence
            return
        end if
        r = matmul(u, vh)
    end subroutine procrustes_complex

    subroutine procrustes_real_complex(a, b, r, scale_factor, status, max_sweeps)
        real(dp), intent(in) :: a(:, :) !! real source matrix (m,n)
        complex(dp), intent(in) :: b(:, :) !! complex target matrix (m,n)
        complex(dp), intent(out) :: r(:, :) !! unitary matrix (n,n)
        real(dp), intent(out) :: scale_factor !! nuclear norm of transpose(A)*B
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! positive SVD Jacobi sweep limit

        call procrustes_complex(cmplx(a, 0.0_dp, dp), b, r, scale_factor, status, max_sweeps)
    end subroutine procrustes_real_complex

    subroutine procrustes_complex_real(a, b, r, scale_factor, status, max_sweeps)
        complex(dp), intent(in) :: a(:, :) !! complex source matrix (m,n)
        real(dp), intent(in) :: b(:, :) !! real target matrix (m,n)
        complex(dp), intent(out) :: r(:, :) !! unitary matrix (n,n)
        real(dp), intent(out) :: scale_factor !! nuclear norm of conjugate-transpose(A)*B
        integer, intent(out) :: status !! linalg_status_* outcome
        integer, intent(in), optional :: max_sweeps !! positive SVD Jacobi sweep limit

        call procrustes_complex(a, cmplx(b, 0.0_dp, dp), r, scale_factor, status, max_sweeps)
    end subroutine procrustes_complex_real

end module scifort_linalg_procrustes
