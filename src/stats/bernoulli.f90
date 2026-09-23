! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Bernoulli distribution matching scipy.stats.bernoulli. This is exactly the
! binomial distribution with one trial, so the implementation delegates to
! scifort_binomial and inherits its endpoint and tail behavior.

module scifort_bernoulli
    use scifort_binomial, only : binomial_cdf, binomial_isf, binomial_logcdf, &
        binomial_logpmf, binomial_logsf, binomial_pmf, binomial_ppf, binomial_sf
    use scifort_kinds, only : dp
    implicit none
    private

    public :: bernoulli_cdf
    public :: bernoulli_isf
    public :: bernoulli_logcdf
    public :: bernoulli_logpmf
    public :: bernoulli_logsf
    public :: bernoulli_pmf
    public :: bernoulli_ppf
    public :: bernoulli_sf

    interface bernoulli_pmf
        module procedure bernoulli_pmf_real
        module procedure bernoulli_pmf_int
    end interface bernoulli_pmf

    interface bernoulli_logpmf
        module procedure bernoulli_logpmf_real
        module procedure bernoulli_logpmf_int
    end interface bernoulli_logpmf

    interface bernoulli_cdf
        module procedure bernoulli_cdf_real
        module procedure bernoulli_cdf_int
    end interface bernoulli_cdf

    interface bernoulli_sf
        module procedure bernoulli_sf_real
        module procedure bernoulli_sf_int
    end interface bernoulli_sf

    interface bernoulli_logcdf
        module procedure bernoulli_logcdf_real
        module procedure bernoulli_logcdf_int
    end interface bernoulli_logcdf

    interface bernoulli_logsf
        module procedure bernoulli_logsf_real
        module procedure bernoulli_logsf_int
    end interface bernoulli_logsf

contains

    pure elemental function bernoulli_pmf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! success count; non-integer values have probability zero
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_pmf(k, 1.0_dp, p, loc)
    end function bernoulli_pmf_real

    pure elemental function bernoulli_pmf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! success count
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_pmf(k, 1, p, loc)
    end function bernoulli_pmf_int

    pure elemental function bernoulli_logpmf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! success count; non-integer values have probability zero
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_logpmf(k, 1.0_dp, p, loc)
    end function bernoulli_logpmf_real

    pure elemental function bernoulli_logpmf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! success count
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_logpmf(k, 1, p, loc)
    end function bernoulli_logpmf_int

    pure elemental function bernoulli_cdf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); the floor is used
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_cdf(k, 1.0_dp, p, loc)
    end function bernoulli_cdf_real

    pure elemental function bernoulli_cdf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_cdf(k, 1, p, loc)
    end function bernoulli_cdf_int

    pure elemental function bernoulli_sf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); the floor is used
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_sf(k, 1.0_dp, p, loc)
    end function bernoulli_sf_real

    pure elemental function bernoulli_sf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_sf(k, 1, p, loc)
    end function bernoulli_sf_int

    pure elemental function bernoulli_logcdf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k)); the floor is used
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_logcdf(k, 1.0_dp, p, loc)
    end function bernoulli_logcdf_real

    pure elemental function bernoulli_logcdf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_logcdf(k, 1, p, loc)
    end function bernoulli_logcdf_int

    pure elemental function bernoulli_logsf_real(k, p, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k)); the floor is used
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_logsf(k, 1.0_dp, p, loc)
    end function bernoulli_logsf_real

    pure elemental function bernoulli_logsf_int(k, p, loc) result(y)
        integer, intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_logsf(k, 1, p, loc)
    end function bernoulli_logsf_int

    pure elemental function bernoulli_ppf(probability, p, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_ppf(probability, 1.0_dp, p, loc)
    end function bernoulli_ppf

    pure elemental function bernoulli_isf(probability, p, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: p !! probability of success in [0, 1]
        real(dp), intent(in), optional :: loc !! shift of the support (default 0)
        real(dp) :: y

        y = binomial_isf(probability, 1.0_dp, p, loc)
    end function bernoulli_isf

end module scifort_bernoulli
