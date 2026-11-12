! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bernoulli numbers, Euler/secant numbers, and Stirling numbers of the second
! kind. The implementation is original SciFort code based on standard defining
! identities recorded in the NIST DLMF and Concrete Mathematics.

module scifort_integer_sequences
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity
    use scifort_special_elementary, only : gammaln
    use scifort_zeta, only : zeta
    implicit none
    private

    integer, parameter, public :: integer_sequences_status_success = 0
    integer, parameter, public :: integer_sequences_status_invalid_argument = 1

    public :: bernoulli
    public :: euler
    public :: stirling2

    real(dp), parameter :: log_huge_dp = log(huge(1.0_dp))
    real(dp), parameter :: log_two_pi = log(2.0_dp * scifort_pi)
    real(dp), parameter :: log_four = log(4.0_dp)
    real(dp), parameter :: log_pi = log(scifort_pi)

contains

    pure subroutine bernoulli(n, values, status)
        integer, intent(in) :: n !! highest Bernoulli index requested, >= 0
        real(dp), allocatable, intent(out) :: values(:) !! B_0 through B_n, indexed 0:n
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        integer :: k
        real(dp) :: log_abs
        real(dp) :: sign_value
        real(dp) :: zeta_value

        if (n < 0) then
            allocate(values(0))
            if (present(status)) status = integer_sequences_status_invalid_argument
            return
        end if

        allocate(values(0:n))
        values = 0.0_dp
        values(0) = 1.0_dp
        if (n >= 1) values(1) = -0.5_dp

        do k = 2, n, 2
            zeta_value = zeta(real(k, dp))
            log_abs = log(2.0_dp) + gammaln(real(k + 1, dp)) + log(zeta_value) - &
                real(k, dp) * log_two_pi
            sign_value = 1.0_dp
            if (mod(k, 4) == 0) sign_value = -1.0_dp
            if (log_abs > log_huge_dp) then
                values(k) = sign_value * positive_infinity(1.0_dp)
            else
                values(k) = sign_value * exp(log_abs)
            end if
        end do

        if (present(status)) status = integer_sequences_status_success
    end subroutine bernoulli

    pure subroutine euler(n, values, status)
        integer, intent(in) :: n !! highest Euler-number index requested, >= 0
        real(dp), allocatable, intent(out) :: values(:) !! E_0 through E_n, indexed 0:n
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        integer :: even_index
        integer :: j
        integer :: m
        real(dp) :: choose_value
        real(dp) :: log_abs
        real(dp) :: sign_value
        real(dp) :: total

        if (n < 0) then
            allocate(values(0))
            if (present(status)) status = integer_sequences_status_invalid_argument
            return
        end if

        allocate(values(0:n))
        values = 0.0_dp
        values(0) = 1.0_dp

        do even_index = 2, n, 2
            m = even_index / 2
            sign_value = 1.0_dp
            if (mod(m, 2) /= 0) sign_value = -1.0_dp

            if (even_index >= 188) then
                ! For s=even_index+1 >= 189, the Dirichlet beta factor in
                ! |E_n| = 4^(m+1) n! beta(n+1) / pi^(n+1) rounds to one in
                ! binary64. The log form avoids overflow in intermediate terms.
                log_abs = real(m + 1, dp) * log_four + gammaln(real(even_index + 1, dp)) - &
                    real(even_index + 1, dp) * log_pi
                if (log_abs > log_huge_dp) then
                    values(even_index) = sign_value * positive_infinity(1.0_dp)
                else
                    values(even_index) = sign_value * exp(log_abs)
                end if
                cycle
            end if

            total = values(0)
            choose_value = 1.0_dp
            do j = 1, m - 1
                choose_value = choose_value * real(even_index - 2 * j + 2, dp) * &
                    real(even_index - 2 * j + 1, dp) / &
                    (real(2 * j - 1, dp) * real(2 * j, dp))
                total = total + choose_value * values(2 * j)
            end do
            values(even_index) = -total
        end do

        if (present(status)) status = integer_sequences_status_success
    end subroutine euler

    pure elemental function stirling2(n, k) result(value)
        integer, intent(in) :: n !! number of set elements
        integer, intent(in) :: k !! number of nonempty subsets
        real(dp) :: value

        integer :: i
        integer :: j
        real(dp), allocatable :: row(:)
        real(dp) :: term

        if (n < 0 .or. k < 0 .or. k > n) then
            value = 0.0_dp
            return
        end if
        if (n == 0) then
            if (k == 0) then
                value = 1.0_dp
            else
                value = 0.0_dp
            end if
            return
        end if
        if (k == 0) then
            value = 0.0_dp
            return
        end if
        if (k == 1 .or. k == n) then
            value = 1.0_dp
            return
        end if

        allocate(row(0:k))
        row = 0.0_dp
        row(0) = 1.0_dp
        do i = 1, n
            do j = min(i, k), 1, -1
                if (.not. ieee_is_finite(row(j))) then
                    term = positive_infinity(1.0_dp)
                else if (row(j) > huge(1.0_dp) / real(j, dp)) then
                    term = positive_infinity(1.0_dp)
                else
                    term = real(j, dp) * row(j)
                end if

                if (.not. ieee_is_finite(term) .or. row(j - 1) > huge(1.0_dp) - term) then
                    row(j) = positive_infinity(1.0_dp)
                else
                    row(j) = row(j - 1) + term
                end if
            end do
            row(0) = 0.0_dp
        end do
        value = row(k)
    end function stirling2

end module scifort_integer_sequences
