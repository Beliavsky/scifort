! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Coefficient constructors for classical orthogonal polynomials.
!
! The recurrence relations follow NIST DLMF Chapter 18. Public coefficient
! arrays use descending powers, matching scipy.special orthopoly1d.coeffs.
! This module is an independent SciFort implementation.

module scifort_orthogonal_constructors
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    public :: chebyc
    public :: chebys
    public :: chebyt
    public :: chebyu
    public :: gegenbauer
    public :: genlaguerre
    public :: hermite
    public :: hermitenorm
    public :: jacobi
    public :: laguerre
    public :: legendre
    public :: sh_chebyt
    public :: sh_chebyu
    public :: sh_jacobi
    public :: sh_legendre

contains

    pure function legendre(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        integer :: k
        logical :: make_monic

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        make_monic = optional_monic(monic)
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        current = [0.0_dp, 1.0_dp]
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
            return
        end if

        do k = 1, n - 1
            next_value = recurrence_step(previous, current, 2.0_dp * real(k, dp) + 1.0_dp, &
                0.0_dp, real(k, dp), real(k + 1, dp))
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function legendre

    pure function jacobi(n, alpha, beta, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        real(dp), intent(in) :: alpha !! Jacobi alpha parameter, > -1
        real(dp), intent(in) :: beta !! Jacobi beta parameter, > -1
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        real(dp) :: a
        real(dp) :: ab
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        integer :: k
        logical :: make_monic

        if (n < 0 .or. .not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(beta) .or. &
                alpha <= -1.0_dp .or. beta <= -1.0_dp) then
            allocate(coefficients(0))
            return
        end if
        make_monic = optional_monic(monic)
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        current = [0.5_dp * (alpha - beta), 0.5_dp * (alpha + beta + 2.0_dp)]
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
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
            next_value = recurrence_step(previous, current, b, c, d, a)
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function jacobi

    pure function gegenbauer(n, alpha, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        real(dp), intent(in) :: alpha !! Gegenbauer parameter, > -1/2
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        integer :: k
        logical :: make_monic

        if (n < 0 .or. .not. ieee_is_finite(alpha) .or. alpha <= -0.5_dp) then
            allocate(coefficients(0))
            return
        end if
        make_monic = optional_monic(monic)
        if (make_monic .and. n > 0) then
            coefficients = monic_gegenbauer_coefficients(n, alpha)
            return
        end if

        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        current = [0.0_dp, 2.0_dp * alpha]
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
            return
        end if

        do k = 1, n - 1
            next_value = recurrence_step(previous, current, &
                2.0_dp * (real(k, dp) + alpha), 0.0_dp, &
                real(k, dp) + 2.0_dp * alpha - 1.0_dp, real(k + 1, dp))
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function gegenbauer

    pure function laguerre(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = genlaguerre(n, 0.0_dp, optional_monic(monic))
    end function laguerre

    pure function genlaguerre(n, alpha, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        real(dp), intent(in) :: alpha !! generalized-Laguerre parameter, > -1
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        integer :: k
        logical :: make_monic

        if (n < 0 .or. .not. ieee_is_finite(alpha) .or. alpha <= -1.0_dp) then
            allocate(coefficients(0))
            return
        end if
        make_monic = optional_monic(monic)
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        current = [alpha + 1.0_dp, -1.0_dp]
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
            return
        end if

        do k = 1, n - 1
            next_value = recurrence_step(previous, current, -1.0_dp, &
                2.0_dp * real(k, dp) + alpha + 1.0_dp, real(k, dp) + alpha, &
                real(k + 1, dp))
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function genlaguerre

    pure function hermite(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = hermite_family(n, .false., optional_monic(monic))
    end function hermite

    pure function hermitenorm(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = hermite_family(n, .true., optional_monic(monic))
    end function hermitenorm

    pure function chebyt(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = chebyshev_family(n, .true., optional_monic(monic))
    end function chebyt

    pure function chebyu(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = chebyshev_family(n, .false., optional_monic(monic))
    end function chebyu

    pure function chebyc(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: base(:)

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        base = chebyt(n)
        coefficients = affine_transform(base, 0.5_dp, 0.0_dp)
        coefficients = 2.0_dp * coefficients
        if (optional_monic(monic)) call normalize_monic(coefficients)
    end function chebyc

    pure function chebys(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: base(:)

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        base = chebyu(n)
        coefficients = affine_transform(base, 0.5_dp, 0.0_dp)
        if (optional_monic(monic)) call normalize_monic(coefficients)
    end function chebys

    pure function sh_legendre(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        integer :: k

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, optional_monic(monic))
            return
        end if
        current = [-1.0_dp, 2.0_dp]
        if (n == 1) then
            coefficients = finalize_coefficients(current, optional_monic(monic))
            return
        end if
        do k = 1, n - 1
            next_value = recurrence_step(previous, current, &
                2.0_dp * (2.0_dp * real(k, dp) + 1.0_dp), &
                -(2.0_dp * real(k, dp) + 1.0_dp), real(k, dp), real(k + 1, dp))
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, optional_monic(monic))
    end function sh_legendre

    pure function sh_jacobi(n, p, q, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        real(dp), intent(in) :: p !! shifted-Jacobi parameter satisfying p - q > -1
        real(dp), intent(in) :: q !! shifted-Jacobi parameter, > 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        real(dp) :: a
        real(dp) :: ab
        real(dp) :: alpha
        real(dp) :: b
        real(dp) :: beta
        real(dp) :: c
        real(dp) :: d
        real(dp) :: scale
        integer :: k

        if (n < 0 .or. .not. ieee_is_finite(p) .or. .not. ieee_is_finite(q) .or. &
                q <= 0.0_dp .or. p - q <= -1.0_dp) then
            allocate(coefficients(0))
            return
        end if
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, optional_monic(monic))
            return
        end if
        current = [-q, p + 1.0_dp]
        if (n == 1) then
            scale = p + 1.0_dp
            coefficients = finalize_coefficients(current / scale, optional_monic(monic))
            return
        end if

        alpha = p - q
        beta = q - 1.0_dp
        ab = alpha + beta
        do k = 1, n - 1
            a = 2.0_dp * real(k + 1, dp) * (real(k + 1, dp) + ab) * &
                (2.0_dp * real(k, dp) + ab)
            b = (2.0_dp * real(k, dp) + ab + 1.0_dp) * &
                (2.0_dp * real(k, dp) + ab) * (2.0_dp * real(k, dp) + ab + 2.0_dp)
            c = (2.0_dp * real(k, dp) + ab + 1.0_dp) * (alpha * alpha - beta * beta)
            d = 2.0_dp * (real(k, dp) + alpha) * (real(k, dp) + beta) * &
                (2.0_dp * real(k, dp) + ab + 2.0_dp)
            next_value = recurrence_step(previous, current, 2.0_dp * b, c - b, d, a)
            previous = current
            current = next_value
        end do
        scale = shifted_jacobi_scale(n, p)
        coefficients = finalize_coefficients(current / scale, optional_monic(monic))
    end function sh_jacobi

    pure function sh_chebyt(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = shifted_chebyshev_family(n, .true., optional_monic(monic))
    end function sh_chebyt

    pure function sh_chebyu(n, monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in), optional :: monic !! normalize leading coefficient to one; default false
        real(dp), allocatable :: coefficients(:)

        coefficients = shifted_chebyshev_family(n, .false., optional_monic(monic))
    end function sh_chebyu

    pure function hermite_family(n, probabilists, make_monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in) :: probabilists !! true for probabilists' Hermite normalization
        logical, intent(in) :: make_monic !! normalize leading coefficient to one
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        real(dp) :: xscale
        integer :: k

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        if (probabilists) then
            xscale = 1.0_dp
        else
            xscale = 2.0_dp
        end if
        current = [0.0_dp, xscale]
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
            return
        end if

        do k = 1, n - 1
            if (probabilists) then
                next_value = recurrence_step(previous, current, 1.0_dp, 0.0_dp, &
                    real(k, dp), 1.0_dp)
            else
                next_value = recurrence_step(previous, current, 2.0_dp, 0.0_dp, &
                    2.0_dp * real(k, dp), 1.0_dp)
            end if
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function hermite_family

    pure function chebyshev_family(n, first_kind, make_monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in) :: first_kind !! true for T_n, false for U_n
        logical, intent(in) :: make_monic !! normalize leading coefficient to one
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        integer :: k

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        if (first_kind) then
            current = [0.0_dp, 1.0_dp]
        else
            current = [0.0_dp, 2.0_dp]
        end if
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
            return
        end if

        do k = 1, n - 1
            next_value = recurrence_step(previous, current, 2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp)
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function chebyshev_family

    pure function monic_gegenbauer_coefficients(n, alpha) result(coefficients)
        integer, intent(in) :: n !! positive polynomial degree
        real(dp), intent(in) :: alpha !! Gegenbauer parameter, > -1/2
        real(dp), allocatable :: coefficients(:)

        real(dp) :: denominator
        real(dp) :: factor
        integer :: j
        integer :: k
        integer :: power

        allocate(coefficients(n + 1), source=0.0_dp)
        coefficients(1) = 1.0_dp
        factor = 1.0_dp
        do k = 1, n / 2
            factor = factor * real((n - 2 * k + 2) * (n - 2 * k + 1), dp) / &
                (4.0_dp * real(k, dp))
            denominator = 1.0_dp
            do j = n - k, n - 1
                denominator = denominator * (alpha + real(j, dp))
            end do
            power = n - 2 * k
            coefficients(n - power + 1) = (-1.0_dp) ** k * factor / denominator
        end do
    end function monic_gegenbauer_coefficients

    pure function shifted_chebyshev_family(n, first_kind, make_monic) result(coefficients)
        integer, intent(in) :: n !! polynomial degree, >= 0
        logical, intent(in) :: first_kind !! true for shifted T_n, false for shifted U_n
        logical, intent(in) :: make_monic !! normalize leading coefficient to one
        real(dp), allocatable :: coefficients(:)

        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: previous(:)
        integer :: k

        if (n < 0) then
            allocate(coefficients(0))
            return
        end if
        previous = [1.0_dp]
        if (n == 0) then
            coefficients = finalize_coefficients(previous, make_monic)
            return
        end if
        if (first_kind) then
            current = [-1.0_dp, 2.0_dp]
        else
            current = [-2.0_dp, 4.0_dp]
        end if
        if (n == 1) then
            coefficients = finalize_coefficients(current, make_monic)
            return
        end if
        do k = 1, n - 1
            next_value = recurrence_step(previous, current, 4.0_dp, -2.0_dp, 1.0_dp, 1.0_dp)
            previous = current
            current = next_value
        end do
        coefficients = finalize_coefficients(current, make_monic)
    end function shifted_chebyshev_family

    pure function recurrence_step(previous, current, xscale, constant, prevscale, denominator) &
            result(next_value)
        real(dp), intent(in) :: previous(:) !! ascending coefficients of degree k-1 polynomial
        real(dp), intent(in) :: current(:) !! ascending coefficients of degree k polynomial
        real(dp), intent(in) :: xscale !! multiplier of x times current polynomial
        real(dp), intent(in) :: constant !! constant multiplier of current polynomial
        real(dp), intent(in) :: prevscale !! subtracted multiplier of previous polynomial
        real(dp), intent(in) :: denominator !! nonzero recurrence denominator
        real(dp), allocatable :: next_value(:)

        allocate(next_value(size(current) + 1), source=0.0_dp)
        next_value(1:size(current)) = next_value(1:size(current)) + constant * current
        next_value(2:size(current) + 1) = next_value(2:size(current) + 1) + xscale * current
        next_value(1:size(previous)) = next_value(1:size(previous)) - prevscale * previous
        next_value = next_value / denominator
    end function recurrence_step

    pure function finalize_coefficients(ascending, make_monic) result(coefficients)
        real(dp), intent(in) :: ascending(:) !! coefficients ordered from constant to leading power
        logical, intent(in) :: make_monic !! normalize leading coefficient to one
        real(dp), allocatable :: coefficients(:)

        integer :: j
        integer :: ncoeff

        ncoeff = size(ascending)
        allocate(coefficients(ncoeff))
        do j = 1, ncoeff
            coefficients(j) = ascending(ncoeff - j + 1)
        end do
        if (make_monic) call normalize_monic(coefficients)
    end function finalize_coefficients

    pure function affine_transform(coefficients, slope, intercept) result(transformed)
        real(dp), intent(in) :: coefficients(:) !! descending coefficients of p(x)
        real(dp), intent(in) :: slope !! affine argument slope in p(slope*x+intercept)
        real(dp), intent(in) :: intercept !! affine argument intercept in p(slope*x+intercept)
        real(dp), allocatable :: transformed(:)

        real(dp), allocatable :: next_value(:)
        real(dp), allocatable :: work(:)
        integer :: j
        integer :: m

        if (size(coefficients) == 0) then
            allocate(transformed(0))
            return
        end if
        work = [coefficients(1)]
        do j = 2, size(coefficients)
            m = size(work)
            allocate(next_value(m + 1), source=0.0_dp)
            next_value(1:m) = next_value(1:m) + slope * work
            next_value(2:m + 1) = next_value(2:m + 1) + intercept * work
            next_value(m + 1) = next_value(m + 1) + coefficients(j)
            call move_alloc(next_value, work)
        end do
        transformed = work
    end function affine_transform

    pure subroutine normalize_monic(coefficients)
        real(dp), intent(inout) :: coefficients(:) !! descending coefficients to normalize in place

        real(dp) :: leading

        if (size(coefficients) == 0) return
        leading = coefficients(1)
        if (.not. ieee_is_finite(leading) .or. abs(leading) <= 0.0_dp) then
            coefficients = 0.0_dp
            return
        end if
        coefficients = coefficients / leading
    end subroutine normalize_monic

    pure elemental logical function optional_monic(monic) result(make_monic)
        logical, intent(in), optional :: monic !! optional monic flag; absent means false

        make_monic = .false.
        if (present(monic)) make_monic = monic
    end function optional_monic

    pure function shifted_jacobi_scale(n, p) result(scale)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: p !! finite shifted-Jacobi p parameter
        real(dp) :: scale

        integer :: k

        scale = 1.0_dp
        do k = 1, n
            scale = scale * (real(n + k, dp) + p - 1.0_dp) / real(k, dp)
        end do
    end function shifted_jacobi_scale

end module scifort_orthogonal_constructors
