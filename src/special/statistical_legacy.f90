! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Historical scipy.special statistical wrappers.
!
! These entry points expose the same standard laws already implemented by the
! SciFort stats layer: finite-sample one-sided Kolmogorov-Smirnov, limiting
! two-sided Kolmogorov-Smirnov, and Tukey lambda.  No independent numerical
! approximation is duplicated here.
module scifort_special_statistical_legacy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_ksone, only : ksone_isf, ksone_sf
    use scifort_kstwobign, only : kstwobign_isf, kstwobign_sf
    use scifort_tukeylambda, only : tukeylambda_cdf
    implicit none
    private

    public :: kolmogi
    public :: kolmogorov
    public :: smirnov
    public :: smirnovi
    public :: tklmbda

contains

    pure elemental function smirnov(n, d) result(p)
        integer, intent(in) :: n !! positive sample size
        real(dp), intent(in) :: d !! one-sided KS deviation in [0,1]
        real(dp) :: p

        if (n <= 0 .or. ieee_is_nan(d) .or. d < 0.0_dp .or. d > 1.0_dp) then
            p = quiet_nan(d)
        else
            p = ksone_sf(d, real(n, dp))
        end if
    end function smirnov

    pure elemental function smirnovi(n, p) result(d)
        integer, intent(in) :: n !! positive sample size
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp) :: d

        d = ksone_isf(p, real(n, dp))
    end function smirnovi

    pure elemental function kolmogorov(y) result(p)
        real(dp), intent(in) :: y !! nonnegative limiting two-sided KS statistic
        real(dp) :: p

        p = kstwobign_sf(y)
    end function kolmogorov

    pure elemental function kolmogi(p) result(y)
        real(dp), intent(in) :: p !! upper-tail probability in [0,1]
        real(dp) :: y

        y = kstwobign_isf(p)
    end function kolmogi

    pure elemental function tklmbda(x, lmbda) result(p)
        real(dp), intent(in) :: x !! standardized Tukey-lambda variate
        real(dp), intent(in) :: lmbda !! finite Tukey-lambda shape parameter
        real(dp) :: p

        p = tukeylambda_cdf(x, lmbda)
    end function tklmbda

end module scifort_special_statistical_legacy
