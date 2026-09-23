! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Logarithmic gamma and beta function helpers shared by the special-function
! modules. Independent implementation from DLMF 5.7.3 (series for
! log(Gamma(1 + a))), DLMF 5.11.1 (Stirling series), and DLMF 4.6.1.
!
! This module is an implementation layer; it is not re-exported by
! scifort_special.

module scifort_log_gamma
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe
    implicit none
    private

    public :: log1p_minus_x
    public :: log_beta
    public :: log_gamma_difference
    public :: log_gamma_ratio_scaled
    public :: log_gamma_one_plus
    public :: stirling_remainder
    public :: stirling_threshold

    ! Arguments at or above this value use the Stirling series.
    real(dp), parameter :: stirling_threshold = 10.0_dp

    real(dp), parameter :: euler_gamma = &
        0.577215664901532860606512090082402431_dp

    ! zeta(k) - 1 for k = 2, ..., 30, generated with mpmath 1.3.0 at 40 digits
    ! and checked against zeta(2) = pi**2 / 6 and zeta(4) = pi**4 / 90.
    real(dp), parameter :: zeta_minus_one(2:30) = [ &
        0.644934066848226436472415166646_dp, &
        0.202056903159594285399738161511_dp, &
        0.0823232337111381915160036965412_dp, &
        0.036927755143369926331365486457_dp, &
        0.0173430619844491397145179297909_dp, &
        0.0083492773819228268397975498498_dp, &
        0.00407735619794433937868523850865_dp, &
        0.00200839282608221441785276923241_dp, &
        0.000994575127818085337145958900319_dp, &
        0.00049418860411946455870228252647_dp, &
        0.00024608655330804829863799804774_dp, &
        0.000122713347578489146751836526357_dp, &
        0.0000612481350587048292585451051353_dp, &
        0.0000305882363070204935517285106451_dp, &
        0.0000152822594086518717325714876367_dp, &
        0.00000763719763789976227360029356303_dp, &
        0.00000381729326499983985646164462194_dp, &
        0.0000019082127165539389256569577951_dp, &
        0.000000953962033872796113152038683449_dp, &
        0.000000476932986787806463116719604373_dp, &
        0.000000238450502727732990003648186753_dp, &
        0.000000119219925965311073067788718882_dp, &
        0.0000000596081890512594796124402079358_dp, &
        0.0000000298035035146522801860637050694_dp, &
        0.0000000149015548283650412346585066307_dp, &
        0.0000000074507117898354294919810041706_dp, &
        0.0000000037253340247884570548192040184_dp, &
        0.00000000186265972351304900640390994542_dp, &
        0.000000000931327432419668182871764735021_dp]

    ! B(2k) / (2k (2k - 1)) for k = 1, ..., 7 (DLMF 5.11.1).
    real(dp), parameter :: stirling_coefficients(7) = [ &
        1.0_dp / 12.0_dp, &
        -1.0_dp / 360.0_dp, &
        1.0_dp / 1260.0_dp, &
        -1.0_dp / 1680.0_dp, &
        1.0_dp / 1188.0_dp, &
        -691.0_dp / 360360.0_dp, &
        1.0_dp / 156.0_dp]

contains

    ! log(Gamma(1 + a)) for a >= 0. For a < 0.5 the series of DLMF 5.7.3 is
    ! used so that the result keeps full relative accuracy as a -> 0; the
    ! intrinsic would first round 1 + a.
    pure elemental function log_gamma_one_plus(a) result(y)
        real(dp), intent(in) :: a !! argument of log(Gamma(1 + a)), >= 0
        real(dp) :: y

        integer :: k
        real(dp) :: power
        real(dp) :: s

        if (a < 0.5_dp) then
            s = 0.0_dp
            power = -a
            do k = 2, 30
                power = -power * a
                s = s + zeta_minus_one(k) * power / real(k, dp)
            end do
            y = -log1p_safe(a) + a * (1.0_dp - euler_gamma) + s
        else
            y = log_gamma(1.0_dp + a)
        end if
    end function log_gamma_one_plus

    ! Stirling series remainder S(a) for a >= 10, from DLMF 5.11.1:
    ! log(Gamma(a)) = (a - 1/2) log(a) - a + log(sqrt(2 pi)) + S(a), and
    ! equivalently log(Gamma(a + 1)) = (a + 1/2) log(a) - a +
    ! log(sqrt(2 pi)) + S(a).
    pure elemental function stirling_remainder(a) result(y)
        real(dp), intent(in) :: a !! argument, >= 10
        real(dp) :: y

        integer :: k
        real(dp) :: inv_a2

        inv_a2 = 1.0_dp / (a * a)
        y = stirling_coefficients(size(stirling_coefficients))
        do k = size(stirling_coefficients) - 1, 1, -1
            y = stirling_coefficients(k) + inv_a2 * y
        end do
        y = y / a
    end function stirling_remainder

    ! log(1 + t) - t for t > -1. For |t| <= 0.25 the Maclaurin series of
    ! DLMF 4.6.1 is summed so that the result keeps full relative accuracy.
    pure elemental function log1p_minus_x(t) result(y)
        real(dp), intent(in) :: t !! argument of log(1 + t) - t, > -1
        real(dp) :: y

        integer :: k
        real(dp) :: power
        real(dp) :: term

        if (abs(t) > 0.25_dp) then
            y = log1p_safe(t) - t
            return
        end if
        power = t * t
        y = -0.5_dp * power
        do k = 3, 80
            power = -power * t
            term = power / real(k, dp)
            y = y - term
            if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(y)) exit
        end do
    end function log1p_minus_x

    ! log(B(a, b)) = log(Gamma(a)) + log(Gamma(b)) - log(Gamma(a + b)) for
    ! finite a, b > 0. When an argument is large, the Stirling series is used
    ! for the large arguments so that the terms of order n log(n) cancel
    ! analytically rather than numerically.
    pure elemental function log_beta(a, b) result(y)
        real(dp), intent(in) :: a !! first argument, finite and > 0
        real(dp), intent(in) :: b !! second argument, finite and > 0
        real(dp) :: y

        real(dp) :: big
        real(dp) :: n
        real(dp) :: small

        small = min(a, b)
        big = max(a, b)
        n = a + b
        if (small >= stirling_threshold) then
            ! Sum of the Stirling forms with n log(n) = a log(n) + b log(n).
            y = small * log(small / n) + big * log1p_safe(-small / n) + &
                0.5_dp * (log(n) - log(small) - log(big)) + &
                scifort_log_sqrt_two_pi + stirling_remainder(small) + &
                stirling_remainder(big) - stirling_remainder(n)
        else if (big >= stirling_threshold) then
            ! log(Gamma(big)) - log(Gamma(n)) from the Stirling forms.
            ! For a very small second argument, obtain log(Gamma(small))
            ! from log(Gamma(1 + small)) - log(small) so the intrinsic does
            ! not lose the small displacement from the pole at zero.
            if (small < 0.5_dp) then
                y = log_gamma_one_plus(small) - log(small)
            else
                y = log_gamma(small)
            end if
            y = y + (big - 0.5_dp) * log1p_safe(-small / n) - &
                small * log(n) + small + stirling_remainder(big) - &
                stirling_remainder(n)
        else
            y = log_gamma(a) + log_gamma(b) - log_gamma(n)
        end if
    end function log_beta

    ! log(Gamma(x + s)) - log(Gamma(x)) for x > 0 and s >= 0, with
    ! relative accuracy as s -> 0. For x >= 10 the Stirling forms are
    ! subtracted analytically:
    ! (x - 1/2) log1p(s/x) + s log(x + s) - s + S(x + s) - S(x),
    ! with each term of S(x + s) - S(x) formed through expm1. Smaller x is
    ! shifted upward with Gamma(x + 1) = x Gamma(x).
    pure elemental function log_gamma_difference(x, s) result(y)
        real(dp), intent(in) :: x !! base argument, > 0
        real(dp), intent(in) :: s !! increment in log(Gamma(x + s)) - log(Gamma(x)), >= 0
        real(dp) :: y

        integer :: j
        integer :: shift
        real(dp) :: log_ratio
        real(dp) :: z

        shift = 0
        if (x < stirling_threshold) shift = ceiling(stirling_threshold - x)
        z = x + real(shift, dp)

        log_ratio = log1p_safe(s / z)
        y = (z - 0.5_dp) * log_ratio + s * log(z + s) - s + &
            stirling_remainder_difference(z, log_ratio)

        do j = 0, shift - 1
            y = y - log1p_safe(s / (x + real(j, dp)))
        end do
    end function log_gamma_difference

    ! log(Gamma(x + s) / (Gamma(x) x**s)) for x >= 10 and s >= 0. From the
    ! Stirling forms, with t = s / x, it equals
    ! (x + s - 1/2) (log1p(t) - t) + s (s - 1/2) / x + S(x + s) - S(x),
    ! in which no terms of order s cancel.
    pure elemental function log_gamma_ratio_scaled(x, s) result(y)
        real(dp), intent(in) :: x !! base argument, >= 10
        real(dp), intent(in) :: s !! increment, >= 0
        real(dp) :: y

        real(dp) :: t

        t = s / x
        y = (x + s - 0.5_dp) * log1p_minus_x(t) + s * (s - 0.5_dp) / x + &
            stirling_remainder_difference(x, log1p_safe(t))
    end function log_gamma_ratio_scaled

    ! S(x + s) - S(x) for x >= 10 given log_ratio = log(1 + s/x), with each
    ! term c_k x**(1 - 2k) ((1 + s/x)**(1 - 2k) - 1) formed through expm1.
    pure elemental function stirling_remainder_difference(x, log_ratio) result(y)
        real(dp), intent(in) :: x !! base argument, >= 10
        real(dp), intent(in) :: log_ratio !! log(1 + s/x) for the increment s
        real(dp) :: y

        integer :: k

        y = 0.0_dp
        do k = 1, size(stirling_coefficients)
            y = y + stirling_coefficients(k) * x**(1 - 2 * k) * &
                expm1_safe(real(1 - 2 * k, dp) * log_ratio)
        end do
    end function stirling_remainder_difference

end module scifort_log_gamma
