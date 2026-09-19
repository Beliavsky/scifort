! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Chi-square distribution with df > 0 degrees of freedom, location loc, and
! scale > 0. It is the gamma distribution with shape df / 2 and scale
! 2 * scale, matching scipy.stats.chi2.

module scifort_chi2
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_gamma, only : gamma_cdf, gamma_isf, gamma_logcdf, gamma_logpdf, &
        gamma_logsf, gamma_pdf, gamma_ppf, gamma_sf
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan, valid_loc_scale
    implicit none
    private

    public :: chi2_cdf
    public :: chi2_isf
    public :: chi2_logcdf
    public :: chi2_logpdf
    public :: chi2_logsf
    public :: chi2_pdf
    public :: chi2_ppf
    public :: chi2_sf

contains

    pure elemental function chi2_pdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else
            y = gamma_pdf(x, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_pdf

    pure elemental function chi2_logpdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else
            y = gamma_logpdf(x, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_logpdf

    pure elemental function chi2_cdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else
            y = gamma_cdf(x, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_cdf

    pure elemental function chi2_sf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else
            y = gamma_sf(x, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_sf

    pure elemental function chi2_logcdf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else
            y = gamma_logcdf(x, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_logcdf

    pure elemental function chi2_logsf(x, df, loc, scale) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(x)
        else
            y = gamma_logsf(x, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_logsf

    pure elemental function chi2_ppf(p, df, loc, scale) result(y)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(p)
        else
            y = gamma_ppf(p, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_ppf

    pure elemental function chi2_isf(p, df, loc, scale) result(y)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_df(df))) then
            y = quiet_nan(p)
        else
            y = gamma_isf(p, 0.5_dp * df, mu, 2.0_dp * sigma)
        end if
    end function chi2_isf

    pure elemental logical function valid_df(df) result(valid)
        real(dp), intent(in) :: df

        valid = ieee_is_finite(df) .and. df > 0.0_dp
    end function valid_df

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc
        real(dp), intent(in), optional :: scale
        real(dp), intent(out) :: mu
        real(dp), intent(out) :: sigma

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_chi2
