! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real factorial and finite combinatorial helpers matching SciPy's floating
! special-function interfaces on the real axis.

module scifort_factorials
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_gamma_combinatorics, only : binom, poch
    use scifort_gamma_family, only : gamma
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    public :: comb
    public :: factorial
    public :: factorial2
    public :: factorialk
    public :: perm

    interface factorial
        module procedure factorial_integer
        module procedure factorial_real
    end interface factorial

    interface comb
        module procedure comb_integer
        module procedure comb_real
    end interface comb

    interface perm
        module procedure perm_integer
        module procedure perm_real
    end interface perm

contains

    pure elemental function factorial_real(n) result(y)
        real(dp), intent(in) :: n !! real factorial argument; negative values use zero extension
        real(dp) :: y

        integer :: factor_index
        integer :: integer_part
        real(dp) :: factor
        real(dp) :: fractional_part

        if (ieee_is_nan(n)) then
            y = n
            return
        end if
        if (n < 0.0_dp) then
            y = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(n)) then
            y = positive_infinity(n)
            return
        end if
        if (n >= 171.0_dp) then
            y = positive_infinity(n)
            return
        end if

        ! Form Gamma(n + 1) by recurrence from an argument in [1, 2).
        ! Near the binary64 overflow boundary, exp(log(Gamma(n + 1))) can
        ! magnify a one-ulp log-gamma error into a relative error of order
        ! 1e-13.  The short positive product avoids that amplification and
        ! is also independent of the platform libm's large-argument log-gamma.
        integer_part = int(floor(n))
        fractional_part = n - real(integer_part, dp)
        y = gamma(1.0_dp + fractional_part)
        do factor_index = 1, integer_part
            factor = fractional_part + real(factor_index, dp)
            if (y > huge(1.0_dp) / factor) then
                y = positive_infinity(n)
                return
            end if
            y = y * factor
        end do
    end function factorial_real

    pure elemental function factorial_integer(n) result(y)
        integer, intent(in) :: n !! integer factorial argument; negative values use zero extension
        real(dp) :: y

        y = factorialk(n, 1)
    end function factorial_integer

    pure elemental function factorial2(n) result(y)
        integer, intent(in) :: n !! integer double-factorial argument; negative values return zero
        real(dp) :: y

        y = factorialk(n, 2)
    end function factorial2

    pure elemental function factorialk(n, k) result(y)
        integer, intent(in) :: n !! integer multifactorial argument; negative values return zero
        integer, intent(in) :: k !! positive multifactorial step/order
        real(dp) :: y

        integer :: factor
        real(dp) :: factor_real

        if (k <= 0) then
            y = quiet_nan(real(n, dp))
            return
        end if
        if (n < 0) then
            y = 0.0_dp
            return
        end if
        if (n == 0) then
            y = 1.0_dp
            return
        end if

        y = 1.0_dp
        factor = n
        do while (factor > 0)
            factor_real = real(factor, dp)
            if (y > huge(1.0_dp) / factor_real) then
                y = positive_infinity(factor_real)
                return
            end if
            y = y * factor_real
            factor = factor - k
        end do
    end function factorialk

    pure elemental function comb_real(n, k, repetition) result(y)
        real(dp), intent(in) :: n !! number of available objects in the floating SciPy convention
        real(dp), intent(in) :: k !! number selected in the floating SciPy convention
        logical, intent(in), optional :: repetition !! allow combinations with repetition; default false
        real(dp) :: y

        logical :: use_repetition
        real(dp) :: upper

        use_repetition = .false.
        if (present(repetition)) use_repetition = repetition

        if (use_repetition) then
            if (ieee_is_nan(n)) then
                y = quiet_nan(n)
                return
            end if
            if (ieee_is_nan(k)) then
                y = quiet_nan(k)
                return
            end if
            if (k == 0.0_dp .and. n >= 0.0_dp) then
                y = 1.0_dp
                return
            end if
            upper = add_minus_one(n, k)
            y = binom(upper, k)
            return
        end if

        if (ieee_is_nan(n) .or. ieee_is_nan(k)) then
            y = 0.0_dp
            return
        end if
        if (k > n) then
            y = 0.0_dp
            return
        end if
        if (n < 0.0_dp) then
            y = 0.0_dp
            return
        end if
        if (k < 0.0_dp) then
            y = 0.0_dp
            return
        end if

        y = binom(n, k)
    end function comb_real

    pure elemental function comb_integer(n, k, repetition) result(y)
        integer, intent(in) :: n !! integer number of available objects
        integer, intent(in) :: k !! integer number selected
        logical, intent(in), optional :: repetition !! allow combinations with repetition; default false
        real(dp) :: y

        if (present(repetition)) then
            y = comb_real(real(n, dp), real(k, dp), repetition)
        else
            y = comb_real(real(n, dp), real(k, dp))
        end if
    end function comb_integer

    pure elemental function perm_real(n, k) result(y)
        real(dp), intent(in) :: n !! number of available objects in the floating SciPy convention
        real(dp), intent(in) :: k !! number selected in the floating SciPy convention
        real(dp) :: y

        if (ieee_is_nan(n) .or. ieee_is_nan(k)) then
            y = 0.0_dp
            return
        end if
        if (k > n) then
            y = 0.0_dp
            return
        end if
        if (n < 0.0_dp) then
            y = 0.0_dp
            return
        end if
        if (k < 0.0_dp) then
            y = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(n) .and. .not. ieee_is_finite(k)) then
            y = quiet_nan(n)
            return
        end if

        y = poch(n - k + 1.0_dp, k)
    end function perm_real

    pure elemental function perm_integer(n, k) result(y)
        integer, intent(in) :: n !! integer number of available objects
        integer, intent(in) :: k !! integer number selected
        real(dp) :: y

        integer :: factor
        integer :: first_factor
        real(dp) :: factor_real

        if (k > n .or. n < 0 .or. k < 0) then
            y = 0.0_dp
            return
        end if
        if (k == 0) then
            y = 1.0_dp
            return
        end if

        first_factor = n - k + 1
        y = 1.0_dp
        do factor = first_factor, n
            factor_real = real(factor, dp)
            if (y > huge(1.0_dp) / factor_real) then
                y = positive_infinity(factor_real)
                return
            end if
            y = y * factor_real
        end do
    end function perm_integer

    pure elemental function add_minus_one(x, y) result(z)
        real(dp), intent(in) :: x !! first addend
        real(dp), intent(in) :: y !! second addend
        real(dp) :: z

        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            if (x > 0.0_dp .and. y < 0.0_dp) then
                z = quiet_nan(x)
            else if (x < 0.0_dp .and. y > 0.0_dp) then
                z = quiet_nan(x)
            else
                z = x + y - 1.0_dp
            end if
            return
        end if

        if (y > 0.0_dp) then
            if (x > huge(1.0_dp) - y) then
                z = positive_infinity(x)
                return
            end if
        else if (y < 0.0_dp) then
            if (x < -huge(1.0_dp) - y) then
                z = -positive_infinity(x)
                return
            end if
        end if
        z = x + y - 1.0_dp
    end function add_minus_one

end module scifort_factorials
