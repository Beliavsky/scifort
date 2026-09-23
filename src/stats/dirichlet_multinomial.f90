! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_dirichlet_multinomial
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    public :: dirichlet_multinomial_cov
    public :: dirichlet_multinomial_logpmf
    public :: dirichlet_multinomial_mean
    public :: dirichlet_multinomial_pmf
    public :: dirichlet_multinomial_var

contains

    function dirichlet_multinomial_logpmf(x, alpha, n) result(y)
        real(dp), intent(in) :: x(:) !! category counts; each value must be a nonnegative integer
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp) :: y

        real(dp) :: alpha0

        if (.not. valid_alpha(alpha) .or. n < 0 .or. size(x) /= size(alpha)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (.not. valid_counts(x)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        if (sum(x) /= real(n, dp)) then
            y = negative_infinity(0.0_dp)
            return
        end if

        alpha0 = sum(alpha)
        y = gammaln(alpha0) + gammaln(real(n + 1, dp)) - &
            gammaln(real(n, dp) + alpha0) + &
            sum(gammaln(x + alpha) - gammaln(alpha) - gammaln(x + 1.0_dp))
    end function dirichlet_multinomial_logpmf

    function dirichlet_multinomial_pmf(x, alpha, n) result(y)
        real(dp), intent(in) :: x(:) !! category counts; each value must be a nonnegative integer
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp) :: y

        y = exp(dirichlet_multinomial_logpmf(x, alpha, n))
    end function dirichlet_multinomial_pmf

    function dirichlet_multinomial_mean(alpha, n) result(mu)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp) :: mu(size(alpha))

        if (.not. valid_alpha(alpha) .or. n < 0) then
            mu = quiet_nan(0.0_dp)
        else
            mu = real(n, dp) * alpha / sum(alpha)
        end if
    end function dirichlet_multinomial_mean

    function dirichlet_multinomial_var(alpha, n) result(v)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp) :: v(size(alpha))

        real(dp) :: alpha0
        real(dp) :: nn

        if (.not. valid_alpha(alpha) .or. n < 0) then
            v = quiet_nan(0.0_dp)
            return
        end if
        alpha0 = sum(alpha)
        nn = real(n, dp)
        v = nn * alpha / alpha0 * (1.0_dp - alpha / alpha0) * &
            (nn + alpha0) / (1.0_dp + alpha0)
    end function dirichlet_multinomial_var

    function dirichlet_multinomial_cov(alpha, n) result(cov)
        real(dp), intent(in) :: alpha(:) !! positive concentration parameters
        integer, intent(in) :: n !! total number of trials, >= 0
        real(dp) :: cov(size(alpha), size(alpha))

        integer :: i
        integer :: j
        real(dp) :: alpha0
        real(dp) :: nn
        real(dp) :: v(size(alpha))

        if (.not. valid_alpha(alpha) .or. n < 0) then
            cov = quiet_nan(0.0_dp)
            return
        end if
        alpha0 = sum(alpha)
        nn = real(n, dp)
        v = dirichlet_multinomial_var(alpha, n)
        do i = 1, size(alpha)
            do j = 1, size(alpha)
                cov(i, j) = -nn * alpha(i) * alpha(j) / (alpha0 * alpha0) * &
                    (nn + alpha0) / (1.0_dp + alpha0)
            end do
            cov(i, i) = v(i)
        end do
    end function dirichlet_multinomial_cov

    pure logical function valid_alpha(alpha) result(valid)
        real(dp), intent(in) :: alpha(:) !! candidate concentration vector

        valid = size(alpha) >= 1 .and. all(ieee_is_finite(alpha)) .and. all(alpha > 0.0_dp)
    end function valid_alpha

    pure logical function valid_counts(x) result(valid)
        real(dp), intent(in) :: x(:) !! candidate category counts

        valid = all(ieee_is_finite(x)) .and. all(x >= 0.0_dp) .and. all(x == floor(x))
    end function valid_counts

end module scifort_dirichlet_multinomial
