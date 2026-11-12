! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Owen's T function for real arguments. This independent implementation uses
! the defining angular integral for |a| <= 1 and Owen's reciprocal identity
! to reduce larger slopes to that well-conditioned interval.

module scifort_owens_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi, scifort_sqrt_two
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    real(dp), parameter :: gl_x(32) = [ &
        -9.972638618494815699e-01_dp, -9.856115115452683817e-01_dp, &
        -9.647622555875063899e-01_dp, -9.349060759377396668e-01_dp, &
        -8.963211557660520912e-01_dp, -8.493676137325699704e-01_dp, &
        -7.944837959679423856e-01_dp, -7.321821187402897113e-01_dp, &
        -6.630442669302152314e-01_dp, -5.877157572407623043e-01_dp, &
        -5.068999089322293594e-01_dp, -4.213512761306353327e-01_dp, &
        -3.318686022821276671e-01_dp, -2.392873622521370647e-01_dp, &
        -1.444719615827964876e-01_dp, -4.830766568773832426e-02_dp, &
         4.830766568773832426e-02_dp,  1.444719615827964876e-01_dp, &
         2.392873622521370647e-01_dp,  3.318686022821276671e-01_dp, &
         4.213512761306353327e-01_dp,  5.068999089322293594e-01_dp, &
         5.877157572407623043e-01_dp,  6.630442669302152314e-01_dp, &
         7.321821187402897113e-01_dp,  7.944837959679423856e-01_dp, &
         8.493676137325699704e-01_dp,  8.963211557660520912e-01_dp, &
         9.349060759377396668e-01_dp,  9.647622555875063899e-01_dp, &
         9.856115115452683817e-01_dp,  9.972638618494815699e-01_dp]
    real(dp), parameter :: gl_w(32) = [ &
        7.018610009470505756e-03_dp, 1.627439473090574323e-02_dp, &
        2.539206530926202410e-02_dp, 3.427386291302176452e-02_dp, &
        4.283589802222683568e-02_dp, 5.099805926237609144e-02_dp, &
        5.868409347853556501e-02_dp, 6.582222277636168295e-02_dp, &
        7.234579410884833806e-02_dp, 7.819389578707022781e-02_dp, &
        8.331192422694670696e-02_dp, 8.765209300440378326e-02_dp, &
        9.117387869576377979e-02_dp, 9.384439908080451087e-02_dp, &
        9.563872007927470831e-02_dp, 9.654008851472765940e-02_dp, &
        9.654008851472765940e-02_dp, 9.563872007927470831e-02_dp, &
        9.384439908080451087e-02_dp, 9.117387869576377979e-02_dp, &
        8.765209300440378326e-02_dp, 8.331192422694670696e-02_dp, &
        7.819389578707022781e-02_dp, 7.234579410884833806e-02_dp, &
        6.582222277636168295e-02_dp, 5.868409347853556501e-02_dp, &
        5.099805926237609144e-02_dp, 4.283589802222683568e-02_dp, &
        3.427386291302176452e-02_dp, 2.539206530926202410e-02_dp, &
        1.627439473090574323e-02_dp, 7.018610009470505756e-03_dp]

    public :: owens_t

contains

    pure elemental function owens_t(h, a) result(t)
        real(dp), intent(in) :: h !! first Owen-T argument
        real(dp), intent(in) :: a !! slope parameter; any real value
        real(dp) :: t

        real(dp) :: aa
        real(dp) :: ah
        real(dp) :: hh
        real(dp) :: qah
        real(dp) :: qh
        real(dp) :: value

        if (ieee_is_nan(h) .or. ieee_is_nan(a)) then
            t = quiet_nan(h + a)
            return
        end if
        if (a == 0.0_dp) then
            t = a
            return
        end if

        hh = abs(h)
        aa = abs(a)
        if (.not. ieee_is_finite(aa)) then
            value = 0.25_dp * erfc(hh / scifort_sqrt_two)
            t = sign(value, a)
            return
        end if
        if (aa <= 1.0_dp) then
            value = owens_t_small(hh, aa)
        else
            if (hh > huge(1.0_dp) / aa) then
                ah = positive_infinity(hh)
            else
                ah = aa * hh
            end if
            qh = 0.5_dp * erfc(hh / scifort_sqrt_two)
            qah = 0.5_dp * erfc(ah / scifort_sqrt_two)
            value = 0.5_dp * (qh + qah) - qh * qah - owens_t_small(ah, 1.0_dp / aa)
        end if

        t = sign(value, a)
    end function owens_t

    pure function owens_t_small(h, a) result(t)
        real(dp), intent(in) :: h !! nonnegative first argument
        real(dp), intent(in) :: a !! slope in [0, 1]
        real(dp) :: t

        real(dp) :: center
        real(dp) :: ct
        real(dp) :: half
        real(dp) :: ratio
        real(dp) :: sum_value
        real(dp) :: theta
        integer :: i

        if (a == 0.0_dp) then
            t = 0.0_dp
            return
        end if

        half = 0.5_dp * atan(a)
        center = half
        sum_value = 0.0_dp
        do i = 1, size(gl_x)
            theta = center + half * gl_x(i)
            ct = cos(theta)
            ratio = h / ct
            if (abs(ratio) <= sqrt(2.0_dp) * sqrt(huge(1.0_dp))) then
                sum_value = sum_value + gl_w(i) * exp((-0.5_dp * ratio) * ratio)
            end if
        end do
        t = half * sum_value / (2.0_dp * scifort_pi)
    end function owens_t_small

end module scifort_owens_t
