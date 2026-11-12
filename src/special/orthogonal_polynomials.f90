! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Classical orthogonal-polynomial evaluators.
!
! The nonnegative-integer-degree recurrences follow NIST DLMF 18.9.
! Shifted families are formed by the standard affine transformations in
! DLMF Chapter 18. This module is an independent implementation.

module scifort_orthogonal_polynomials
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: assoc_laguerre
    public :: eval_chebyc
    public :: eval_chebys
    public :: eval_chebyt
    public :: eval_chebyu
    public :: eval_genlaguerre
    public :: eval_gegenbauer
    public :: eval_hermite
    public :: eval_hermitenorm
    public :: eval_jacobi
    public :: eval_laguerre
    public :: eval_legendre
    public :: eval_sh_chebyt
    public :: eval_sh_chebyu
    public :: eval_sh_jacobi
    public :: eval_sh_legendre

    interface assoc_laguerre
        module procedure assoc_laguerre_real
        module procedure assoc_laguerre_complex
    end interface assoc_laguerre

    interface eval_chebyc
        module procedure eval_chebyc_real
        module procedure eval_chebyc_complex
    end interface eval_chebyc

    interface eval_chebys
        module procedure eval_chebys_real
        module procedure eval_chebys_complex
    end interface eval_chebys

    interface eval_chebyt
        module procedure eval_chebyt_real
        module procedure eval_chebyt_complex
    end interface eval_chebyt

    interface eval_chebyu
        module procedure eval_chebyu_real
        module procedure eval_chebyu_complex
    end interface eval_chebyu

    interface eval_genlaguerre
        module procedure eval_genlaguerre_real
        module procedure eval_genlaguerre_complex
    end interface eval_genlaguerre

    interface eval_gegenbauer
        module procedure eval_gegenbauer_real
        module procedure eval_gegenbauer_complex
    end interface eval_gegenbauer

    interface eval_jacobi
        module procedure eval_jacobi_real
        module procedure eval_jacobi_complex
    end interface eval_jacobi

    interface eval_laguerre
        module procedure eval_laguerre_real
        module procedure eval_laguerre_complex
    end interface eval_laguerre

    interface eval_legendre
        module procedure eval_legendre_real
        module procedure eval_legendre_complex
    end interface eval_legendre

    interface eval_sh_chebyt
        module procedure eval_sh_chebyt_real
        module procedure eval_sh_chebyt_complex
    end interface eval_sh_chebyt

    interface eval_sh_chebyu
        module procedure eval_sh_chebyu_real
        module procedure eval_sh_chebyu_complex
    end interface eval_sh_chebyu

    interface eval_sh_jacobi
        module procedure eval_sh_jacobi_real
        module procedure eval_sh_jacobi_complex
    end interface eval_sh_jacobi

    interface eval_sh_legendre
        module procedure eval_sh_legendre_real
        module procedure eval_sh_legendre_complex
    end interface eval_sh_legendre

contains

    pure elemental function assoc_laguerre_real(x, n, k) result(value)
        real(dp), intent(in) :: x !! finite real evaluation argument
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in), optional :: k !! generalized-Laguerre order, > -1; default 0
        real(dp) :: value

        real(dp) :: alpha

        alpha = 0.0_dp
        if (present(k)) alpha = k
        value = eval_genlaguerre_real(n, alpha, x)
    end function assoc_laguerre_real

    pure elemental function assoc_laguerre_complex(x, n, k) result(value)
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in), optional :: k !! generalized-Laguerre order, > -1; default 0
        complex(dp) :: value

        real(dp) :: alpha

        alpha = 0.0_dp
        if (present(k)) alpha = k
        value = eval_genlaguerre_complex(n, alpha, x)
    end function assoc_laguerre_complex

    pure elemental function eval_legendre_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_legendre_complex(n, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_legendre_real

    pure elemental function eval_legendre_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        complex(dp) :: current
        complex(dp) :: previous
        complex(dp) :: next_value
        integer :: k

        if (.not. valid_degree_argument(n, x)) then
            value = complex_nan()
            return
        end if
        if (n == 0) then
            value = cmplx(1.0_dp, 0.0_dp, dp)
            return
        end if

        previous = cmplx(1.0_dp, 0.0_dp, dp)
        current = x
        if (n == 1) then
            value = current
            return
        end if

        do k = 1, n - 1
            next_value = ((2.0_dp * real(k, dp) + 1.0_dp) * x * current - &
                real(k, dp) * previous) / real(k + 1, dp)
            previous = current
            current = next_value
        end do
        value = current
    end function eval_legendre_complex

    pure elemental function eval_jacobi_real(n, alpha, beta, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: alpha !! Jacobi alpha parameter, > -1
        real(dp), intent(in) :: beta !! Jacobi beta parameter, > -1
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_jacobi_complex(n, alpha, beta, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_jacobi_real

    pure elemental function eval_jacobi_complex(n, alpha, beta, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: alpha !! Jacobi alpha parameter, > -1
        real(dp), intent(in) :: beta !! Jacobi beta parameter, > -1
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        real(dp) :: a
        real(dp) :: ab
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        complex(dp) :: current
        complex(dp) :: previous
        complex(dp) :: next_value
        integer :: k

        if (.not. valid_degree_argument(n, x) .or. .not. ieee_is_finite(alpha) .or. &
                .not. ieee_is_finite(beta) .or. alpha <= -1.0_dp .or. beta <= -1.0_dp) then
            value = complex_nan()
            return
        end if
        if (n == 0) then
            value = cmplx(1.0_dp, 0.0_dp, dp)
            return
        end if

        previous = cmplx(1.0_dp, 0.0_dp, dp)
        current = 0.5_dp * ((alpha - beta) + (alpha + beta + 2.0_dp) * x)
        if (n == 1) then
            value = current
            return
        end if

        ab = alpha + beta
        do k = 1, n - 1
            a = 2.0_dp * real(k + 1, dp) * (real(k + 1, dp) + ab) * &
                (2.0_dp * real(k, dp) + ab)
            b = (2.0_dp * real(k, dp) + ab + 1.0_dp) * &
                (2.0_dp * real(k, dp) + ab) * (2.0_dp * real(k, dp) + ab + 2.0_dp)
            c = (2.0_dp * real(k, dp) + ab + 1.0_dp) * (alpha * alpha - beta * beta)
            d = 2.0_dp * (real(k, dp) + alpha) * (real(k, dp) + beta) * &
                (2.0_dp * real(k, dp) + ab + 2.0_dp)
            next_value = ((b * x + c) * current - d * previous) / a
            previous = current
            current = next_value
        end do
        value = current
    end function eval_jacobi_complex

    pure elemental function eval_gegenbauer_real(n, alpha, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: alpha !! Gegenbauer parameter, > -1/2
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_gegenbauer_complex(n, alpha, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_gegenbauer_real

    pure elemental function eval_gegenbauer_complex(n, alpha, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: alpha !! Gegenbauer parameter, > -1/2
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        complex(dp) :: current
        complex(dp) :: previous
        complex(dp) :: next_value
        integer :: k

        if (.not. valid_degree_argument(n, x) .or. .not. ieee_is_finite(alpha) .or. &
                alpha <= -0.5_dp) then
            value = complex_nan()
            return
        end if
        if (n == 0) then
            value = cmplx(1.0_dp, 0.0_dp, dp)
            return
        end if

        previous = cmplx(1.0_dp, 0.0_dp, dp)
        current = 2.0_dp * alpha * x
        if (n == 1) then
            value = current
            return
        end if

        do k = 1, n - 1
            next_value = (2.0_dp * (real(k, dp) + alpha) * x * current - &
                (real(k, dp) + 2.0_dp * alpha - 1.0_dp) * previous) / real(k + 1, dp)
            previous = current
            current = next_value
        end do
        value = current
    end function eval_gegenbauer_complex

    pure elemental function eval_laguerre_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_genlaguerre_complex(n, 0.0_dp, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_laguerre_real

    pure elemental function eval_laguerre_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        value = eval_genlaguerre_complex(n, 0.0_dp, x)
    end function eval_laguerre_complex

    pure elemental function eval_genlaguerre_real(n, alpha, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: alpha !! generalized-Laguerre parameter, > -1
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_genlaguerre_complex(n, alpha, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_genlaguerre_real

    pure elemental function eval_genlaguerre_complex(n, alpha, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: alpha !! generalized-Laguerre parameter, > -1
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        complex(dp) :: current
        complex(dp) :: previous
        complex(dp) :: next_value
        integer :: k

        if (.not. valid_degree_argument(n, x) .or. .not. ieee_is_finite(alpha) .or. &
                alpha <= -1.0_dp) then
            value = complex_nan()
            return
        end if
        if (n == 0) then
            value = cmplx(1.0_dp, 0.0_dp, dp)
            return
        end if

        previous = cmplx(1.0_dp, 0.0_dp, dp)
        current = alpha + 1.0_dp - x
        if (n == 1) then
            value = current
            return
        end if

        do k = 1, n - 1
            next_value = ((2.0_dp * real(k, dp) + alpha + 1.0_dp - x) * current - &
                (real(k, dp) + alpha) * previous) / real(k + 1, dp)
            previous = current
            current = next_value
        end do
        value = current
    end function eval_genlaguerre_complex

    pure elemental function eval_hermite(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        real(dp) :: current
        real(dp) :: previous
        real(dp) :: next_value
        integer :: k

        if (n < 0 .or. ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if
        if (n == 0) then
            value = 1.0_dp
            return
        end if

        previous = 1.0_dp
        current = 2.0_dp * x
        if (n == 1) then
            value = current
            return
        end if
        do k = 1, n - 1
            next_value = 2.0_dp * x * current - 2.0_dp * real(k, dp) * previous
            previous = current
            current = next_value
        end do
        value = current
    end function eval_hermite

    pure elemental function eval_hermitenorm(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        real(dp) :: current
        real(dp) :: previous
        real(dp) :: next_value
        integer :: k

        if (n < 0 .or. ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            return
        end if
        if (n == 0) then
            value = 1.0_dp
            return
        end if

        previous = 1.0_dp
        current = x
        if (n == 1) then
            value = current
            return
        end if
        do k = 1, n - 1
            next_value = x * current - real(k, dp) * previous
            previous = current
            current = next_value
        end do
        value = current
    end function eval_hermitenorm

    pure elemental function eval_chebyt_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_chebyt_complex(n, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_chebyt_real

    pure elemental function eval_chebyt_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        value = chebyshev_recurrence(n, x, .true.)
    end function eval_chebyt_complex

    pure elemental function eval_chebyu_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = real(eval_chebyu_complex(n, cmplx(x, 0.0_dp, dp)), dp)
    end function eval_chebyu_real

    pure elemental function eval_chebyu_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        value = chebyshev_recurrence(n, x, .false.)
    end function eval_chebyu_complex

    pure elemental function eval_chebyc_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = 2.0_dp * eval_chebyt_real(n, 0.5_dp * x)
    end function eval_chebyc_real

    pure elemental function eval_chebyc_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        value = 2.0_dp * eval_chebyt_complex(n, 0.5_dp * x)
    end function eval_chebyc_complex

    pure elemental function eval_chebys_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp) :: value

        value = eval_chebyu_real(n, 0.5_dp * x)
    end function eval_chebys_real

    pure elemental function eval_chebys_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        complex(dp) :: value

        value = eval_chebyu_complex(n, 0.5_dp * x)
    end function eval_chebys_complex

    pure elemental function eval_sh_legendre_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite shifted argument
        real(dp) :: value

        value = eval_legendre_real(n, 2.0_dp * x - 1.0_dp)
    end function eval_sh_legendre_real

    pure elemental function eval_sh_legendre_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite shifted complex argument
        complex(dp) :: value

        value = eval_legendre_complex(n, 2.0_dp * x - 1.0_dp)
    end function eval_sh_legendre_complex

    pure elemental function eval_sh_chebyt_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite shifted argument
        real(dp) :: value

        value = eval_chebyt_real(n, 2.0_dp * x - 1.0_dp)
    end function eval_sh_chebyt_real

    pure elemental function eval_sh_chebyt_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite shifted complex argument
        complex(dp) :: value

        value = eval_chebyt_complex(n, 2.0_dp * x - 1.0_dp)
    end function eval_sh_chebyt_complex

    pure elemental function eval_sh_chebyu_real(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: x !! finite shifted argument
        real(dp) :: value

        value = eval_chebyu_real(n, 2.0_dp * x - 1.0_dp)
    end function eval_sh_chebyu_real

    pure elemental function eval_sh_chebyu_complex(n, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite shifted complex argument
        complex(dp) :: value

        value = eval_chebyu_complex(n, 2.0_dp * x - 1.0_dp)
    end function eval_sh_chebyu_complex

    pure elemental function eval_sh_jacobi_real(n, p, q, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: p !! shifted-Jacobi parameter satisfying p - q > -1
        real(dp), intent(in) :: q !! shifted-Jacobi parameter, > 0
        real(dp), intent(in) :: x !! finite shifted argument
        real(dp) :: value

        real(dp) :: scale

        if (.not. valid_shifted_jacobi_parameters(p, q)) then
            value = quiet_nan(x)
            return
        end if
        scale = jacobi_shift_scale(n, p)
        if (ieee_is_nan(scale)) then
            value = quiet_nan(x)
        else
            value = eval_jacobi_real(n, p - q, q - 1.0_dp, 2.0_dp * x - 1.0_dp) / scale
        end if
    end function eval_sh_jacobi_real

    pure elemental function eval_sh_jacobi_complex(n, p, q, x) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: p !! shifted-Jacobi parameter satisfying p - q > -1
        real(dp), intent(in) :: q !! shifted-Jacobi parameter, > 0
        complex(dp), intent(in) :: x !! finite shifted complex argument
        complex(dp) :: value

        real(dp) :: scale

        if (.not. valid_shifted_jacobi_parameters(p, q)) then
            value = complex_nan()
            return
        end if
        scale = jacobi_shift_scale(n, p)
        if (ieee_is_nan(scale)) then
            value = complex_nan()
        else
            value = eval_jacobi_complex(n, p - q, q - 1.0_dp, 2.0_dp * x - 1.0_dp) / scale
        end if
    end function eval_sh_jacobi_complex

    pure elemental function chebyshev_recurrence(n, x, first_kind) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: x !! finite complex evaluation argument
        logical, intent(in) :: first_kind !! true for T_n, false for U_n
        complex(dp) :: value

        complex(dp) :: current
        complex(dp) :: previous
        complex(dp) :: next_value
        integer :: k

        if (.not. valid_degree_argument(n, x)) then
            value = complex_nan()
            return
        end if
        if (n == 0) then
            value = cmplx(1.0_dp, 0.0_dp, dp)
            return
        end if

        previous = cmplx(1.0_dp, 0.0_dp, dp)
        if (first_kind) then
            current = x
        else
            current = 2.0_dp * x
        end if
        if (n == 1) then
            value = current
            return
        end if

        do k = 1, n - 1
            next_value = 2.0_dp * x * current - previous
            previous = current
            current = next_value
        end do
        value = current
    end function chebyshev_recurrence

    pure elemental function jacobi_shift_scale(n, p) result(scale)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: p !! finite shifted-Jacobi p parameter
        real(dp) :: scale

        integer :: k

        if (n < 0 .or. .not. ieee_is_finite(p)) then
            scale = quiet_nan(p)
            return
        end if
        scale = 1.0_dp
        do k = 1, n
            scale = scale * (real(n + k, dp) + p - 1.0_dp) / real(k, dp)
        end do
    end function jacobi_shift_scale

    pure elemental logical function valid_degree_argument(n, x) result(valid)
        integer, intent(in) :: n !! candidate nonnegative degree
        complex(dp), intent(in) :: x !! candidate finite complex argument

        valid = n >= 0 .and. ieee_is_finite(real(x, dp)) .and. ieee_is_finite(aimag(x))
    end function valid_degree_argument

    pure elemental logical function valid_shifted_jacobi_parameters(p, q) result(valid)
        real(dp), intent(in) :: p !! shifted-Jacobi p parameter
        real(dp), intent(in) :: q !! shifted-Jacobi q parameter

        valid = ieee_is_finite(p) .and. ieee_is_finite(q) .and. q > 0.0_dp .and. p - q > -1.0_dp
    end function valid_shifted_jacobi_parameters

    pure function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), dp)
    end function complex_nan

end module scifort_orthogonal_polynomials
