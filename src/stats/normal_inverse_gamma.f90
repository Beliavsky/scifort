! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_normal_inverse_gamma
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_invgamma, only : invgamma_ppf
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, quiet_nan
    use scifort_normal, only : normal_ppf
    use scifort_random, only : rng_state, rng_uniform
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    public :: normal_inverse_gamma_logpdf
    public :: normal_inverse_gamma_mean
    public :: normal_inverse_gamma_pdf
    public :: normal_inverse_gamma_rvs
    public :: normal_inverse_gamma_rvs_array
    public :: normal_inverse_gamma_var

contains

    pure elemental function normal_inverse_gamma_logpdf(x, s2, mu, lmbda, a, b) result(y)
        real(dp), intent(in) :: x !! normal component
        real(dp), intent(in) :: s2 !! variance component; must be positive
        real(dp), intent(in), optional :: mu !! location parameter; default 0
        real(dp), intent(in), optional :: lmbda !! positive normal precision multiplier; default 1
        real(dp), intent(in), optional :: a !! positive inverse-gamma shape; default 1
        real(dp), intent(in), optional :: b !! positive inverse-gamma scale; default 1
        real(dp) :: y

        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: lam
        real(dp) :: location

        call get_parameters(mu, lmbda, a, b, location, lam, aa, bb)
        if (.not. valid_parameters(location, lam, aa, bb)) then
            y = quiet_nan(x)
        else if (s2 <= 0.0_dp) then
            y = negative_infinity(x)
        else if (ieee_is_nan(x) .or. ieee_is_nan(s2)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(s2)) then
            y = negative_infinity(x)
        else
            y = 0.5_dp * (log(lam) - log(2.0_dp * acos(-1.0_dp) * s2)) + &
                aa * log(bb) - gammaln(aa) - (aa + 1.0_dp) * log(s2) - &
                (2.0_dp * bb + lam * (x - location)**2) / (2.0_dp * s2)
        end if
    end function normal_inverse_gamma_logpdf

    pure elemental function normal_inverse_gamma_pdf(x, s2, mu, lmbda, a, b) result(y)
        real(dp), intent(in) :: x !! normal component
        real(dp), intent(in) :: s2 !! variance component; must be positive
        real(dp), intent(in), optional :: mu !! location parameter; default 0
        real(dp), intent(in), optional :: lmbda !! positive normal precision multiplier; default 1
        real(dp), intent(in), optional :: a !! positive inverse-gamma shape; default 1
        real(dp), intent(in), optional :: b !! positive inverse-gamma scale; default 1
        real(dp) :: y

        y = exp(normal_inverse_gamma_logpdf(x, s2, mu, lmbda, a, b))
    end function normal_inverse_gamma_pdf

    pure subroutine normal_inverse_gamma_mean(mu, lmbda, a, b, mean_x, mean_s2)
        real(dp), intent(in), optional :: mu !! location parameter; default 0
        real(dp), intent(in), optional :: lmbda !! positive normal precision multiplier; default 1
        real(dp), intent(in), optional :: a !! inverse-gamma shape; mean requires a > 1
        real(dp), intent(in), optional :: b !! positive inverse-gamma scale; default 1
        real(dp), intent(out) :: mean_x !! mean of x
        real(dp), intent(out) :: mean_s2 !! mean of variance component

        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: lam
        real(dp) :: location

        call get_parameters(mu, lmbda, a, b, location, lam, aa, bb)
        if (.not. valid_parameters(location, lam, aa, bb) .or. aa <= 1.0_dp) then
            mean_x = quiet_nan(0.0_dp)
            mean_s2 = quiet_nan(0.0_dp)
        else
            mean_x = location
            mean_s2 = bb / (aa - 1.0_dp)
        end if
    end subroutine normal_inverse_gamma_mean

    pure subroutine normal_inverse_gamma_var(mu, lmbda, a, b, var_x, var_s2)
        real(dp), intent(in), optional :: mu !! location parameter; default 0
        real(dp), intent(in), optional :: lmbda !! positive normal precision multiplier; default 1
        real(dp), intent(in), optional :: a !! inverse-gamma shape
        real(dp), intent(in), optional :: b !! positive inverse-gamma scale; default 1
        real(dp), intent(out) :: var_x !! variance of x; finite for a > 1
        real(dp), intent(out) :: var_s2 !! variance of variance component; finite for a > 2

        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: lam
        real(dp) :: location

        call get_parameters(mu, lmbda, a, b, location, lam, aa, bb)
        if (.not. valid_parameters(location, lam, aa, bb)) then
            var_x = quiet_nan(0.0_dp)
            var_s2 = quiet_nan(0.0_dp)
            return
        end if
        if (aa > 1.0_dp) then
            var_x = bb / ((aa - 1.0_dp) * lam)
        else
            var_x = quiet_nan(0.0_dp)
        end if
        if (aa > 2.0_dp) then
            var_s2 = bb * bb / ((aa - 1.0_dp)**2 * (aa - 2.0_dp))
        else
            var_s2 = quiet_nan(0.0_dp)
        end if
    end subroutine normal_inverse_gamma_var

    subroutine normal_inverse_gamma_rvs(state, x, s2, mu, lmbda, a, b)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(out) :: x !! sampled normal component
        real(dp), intent(out) :: s2 !! sampled variance component
        real(dp), intent(in), optional :: mu !! location parameter; default 0
        real(dp), intent(in), optional :: lmbda !! positive normal precision multiplier; default 1
        real(dp), intent(in), optional :: a !! positive inverse-gamma shape; default 1
        real(dp), intent(in), optional :: b !! positive inverse-gamma scale; default 1

        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: lam
        real(dp) :: location

        call get_parameters(mu, lmbda, a, b, location, lam, aa, bb)
        if (.not. valid_parameters(location, lam, aa, bb)) then
            x = quiet_nan(0.0_dp)
            s2 = quiet_nan(0.0_dp)
            return
        end if
        s2 = invgamma_ppf(rng_uniform(state), aa, scale=bb)
        x = normal_ppf(rng_uniform(state), location, sqrt(s2 / lam))
    end subroutine normal_inverse_gamma_rvs

    subroutine normal_inverse_gamma_rvs_array(state, x, s2, mu, lmbda, a, b)
        type(rng_state), intent(inout) :: state !! explicit random-number generator state
        real(dp), intent(out) :: x(:) !! sampled normal components
        real(dp), intent(out) :: s2(:) !! sampled variance components, same length as x
        real(dp), intent(in), optional :: mu !! location parameter; default 0
        real(dp), intent(in), optional :: lmbda !! positive normal precision multiplier; default 1
        real(dp), intent(in), optional :: a !! positive inverse-gamma shape; default 1
        real(dp), intent(in), optional :: b !! positive inverse-gamma scale; default 1

        integer :: i

        if (size(s2) /= size(x)) then
            x = quiet_nan(0.0_dp)
            s2 = quiet_nan(0.0_dp)
            return
        end if
        do i = 1, size(x)
            call normal_inverse_gamma_rvs(state, x(i), s2(i), mu, lmbda, a, b)
        end do
    end subroutine normal_inverse_gamma_rvs_array

    pure subroutine get_parameters(mu, lmbda, a, b, location, lam, aa, bb)
        real(dp), intent(in), optional :: mu !! optional location
        real(dp), intent(in), optional :: lmbda !! optional precision multiplier
        real(dp), intent(in), optional :: a !! optional inverse-gamma shape
        real(dp), intent(in), optional :: b !! optional inverse-gamma scale
        real(dp), intent(out) :: location !! resolved location
        real(dp), intent(out) :: lam !! resolved precision multiplier
        real(dp), intent(out) :: aa !! resolved shape
        real(dp), intent(out) :: bb !! resolved scale

        location = 0.0_dp
        lam = 1.0_dp
        aa = 1.0_dp
        bb = 1.0_dp
        if (present(mu)) location = mu
        if (present(lmbda)) lam = lmbda
        if (present(a)) aa = a
        if (present(b)) bb = b
    end subroutine get_parameters

    pure logical function valid_parameters(mu, lmbda, a, b) result(valid)
        real(dp), intent(in) :: mu !! location parameter
        real(dp), intent(in) :: lmbda !! normal precision multiplier
        real(dp), intent(in) :: a !! inverse-gamma shape
        real(dp), intent(in) :: b !! inverse-gamma scale

        valid = ieee_is_finite(mu) .and. ieee_is_finite(lmbda) .and. &
            ieee_is_finite(a) .and. ieee_is_finite(b) .and. lmbda > 0.0_dp .and. &
            a > 0.0_dp .and. b > 0.0_dp
    end function valid_parameters

end module scifort_normal_inverse_gamma
