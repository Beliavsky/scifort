! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Jahnke-Emden Lambda sequences and the Bessel-polynomial weighted integral.
! The implementation is original SciFort code based on the defining
! hypergeometric identity for Lambda_v and the defining Bessel integral.

module scifort_bessel_misc
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_gamma_family, only : rgamma
    use scifort_hypergeometric, only : hyp0f1
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_bessel_jy, only : jv
    implicit none
    private

    integer, parameter, public :: bessel_misc_status_success = 0
    integer, parameter, public :: bessel_misc_status_invalid_argument = 1
    integer, parameter :: sequence_limit = 100000

    real(dp), parameter :: gl16_nodes(8) = [ &
        0.09501250983763744051_dp, &
        0.2816035507792589154_dp, &
        0.4580167776572273697_dp, &
        0.6178762444026437706_dp, &
        0.7554044083550029987_dp, &
        0.8656312023878317552_dp, &
        0.9445750230732326003_dp, &
        0.9894009349916499385_dp]
    real(dp), parameter :: gl16_weights(8) = [ &
        0.1894506104550686409_dp, &
        0.1826034150449236393_dp, &
        0.1691565193950026469_dp, &
        0.1495959888165767082_dp, &
        0.1246289712555340712_dp, &
        0.09515851168249260528_dp, &
        0.06225352393864745648_dp, &
        0.02715245941175417621_dp]

    public :: besselpoly
    public :: lmbda

contains

    pure subroutine lmbda(v, x, values, derivatives, status)
        real(dp), intent(in) :: v !! highest Lambda order; finite and >= 0
        real(dp), intent(in) :: x !! finite real evaluation argument
        real(dp), allocatable, intent(out) :: values(:) !! Lambda values from frac(v) through v
        real(dp), allocatable, intent(out) :: derivatives(:) !! corresponding derivatives in x
        integer, intent(out), optional :: status !! zero on success; nonzero on invalid input

        integer :: k
        integer :: n
        real(dp) :: order
        real(dp) :: start_order
        real(dp) :: z

        if (.not. ieee_is_finite(v) .or. .not. ieee_is_finite(x) .or. v < 0.0_dp .or. &
                v > real(sequence_limit, dp)) then
            allocate(values(0), derivatives(0))
            if (present(status)) status = bessel_misc_status_invalid_argument
            return
        end if

        n = int(v)
        start_order = v - real(n, dp)
        allocate(values(n + 1), derivatives(n + 1))
        z = -0.25_dp * x * x
        do k = 0, n
            order = start_order + real(k, dp)
            values(k + 1) = hyp0f1(order + 1.0_dp, z)
            if (x == 0.0_dp) then
                derivatives(k + 1) = 0.0_dp
            else
                derivatives(k + 1) = -0.5_dp * x / (order + 1.0_dp) * &
                    hyp0f1(order + 2.0_dp, z)
            end if
        end do
        if (present(status)) status = bessel_misc_status_success
    end subroutine lmbda

    pure elemental function besselpoly(a, lmb, nu) result(value)
        real(dp), intent(in) :: a !! nonnegative Bessel scale factor
        real(dp), intent(in) :: lmb !! real power on the integration variable
        real(dp), intent(in) :: nu !! real Bessel order
        real(dp) :: value

        integer :: integer_order
        integer :: panel
        integer :: panels
        integer :: transform_power
        real(dp) :: effective_exponent
        real(dp) :: half
        real(dp) :: left
        real(dp) :: midpoint
        real(dp) :: right

        if (ieee_is_nan(a) .or. ieee_is_nan(lmb) .or. ieee_is_nan(nu) .or. &
                .not. ieee_is_finite(a) .or. .not. ieee_is_finite(lmb) .or. &
                .not. ieee_is_finite(nu) .or. a < 0.0_dp) then
            value = quiet_nan(a + lmb + nu)
            return
        end if

        if (a == 0.0_dp) then
            if (nu == 0.0_dp) then
                if (lmb > -1.0_dp) then
                    value = 1.0_dp / (lmb + 1.0_dp)
                else
                    value = positive_infinity(1.0_dp)
                end if
            else if (nu > 0.0_dp .or. is_negative_integer(nu)) then
                value = 0.0_dp
            else
                value = quiet_nan(a)
            end if
            return
        end if

        if (is_negative_integer(nu)) then
            integer_order = nint(-nu)
            effective_exponent = lmb + real(integer_order, dp)
        else
            effective_exponent = lmb + nu
        end if

        if (effective_exponent <= -1.0_dp) then
            value = divergent_besselpoly_sign(nu) * positive_infinity(1.0_dp)
            return
        end if

        if (effective_exponent < -0.5_dp) then
            value = besselpoly_endpoint_singular(a, nu, effective_exponent)
            return
        end if

        transform_power = max(1, min(16, ceiling(5.0_dp / (effective_exponent + 1.0_dp))))
        panels = max(1, ceiling(max(a, 1.0_dp) * real(transform_power, dp) / 2.0_dp))
        value = 0.0_dp
        do panel = 0, panels - 1
            left = real(panel, dp) / real(panels, dp)
            right = real(panel + 1, dp) / real(panels, dp)
            midpoint = 0.5_dp * (left + right)
            half = 0.5_dp * (right - left)
            call besselpoly_panel_pair(midpoint, half, transform_power, a, lmb, nu, value)
        end do
    end function besselpoly


    pure function besselpoly_endpoint_singular(a, nu, effective_exponent) result(value)
        real(dp), intent(in) :: a !! positive finite Bessel scale
        real(dp), intent(in) :: nu !! finite Bessel order
        real(dp), intent(in) :: effective_exponent !! leading endpoint power, > -1
        real(dp) :: value

        integer :: integer_order
        integer :: panel
        integer :: panels
        integer :: transform_power
        real(dp) :: coefficient
        real(dp) :: half
        real(dp) :: left
        real(dp) :: midpoint
        real(dp) :: normalized_order
        real(dp) :: residual_exponent
        real(dp) :: right
        real(dp) :: sign_value

        sign_value = 1.0_dp
        if (is_negative_integer(nu)) then
            integer_order = nint(-nu)
            normalized_order = real(integer_order, dp)
            if (mod(integer_order, 2) /= 0) sign_value = -1.0_dp
        else
            normalized_order = nu
        end if
        coefficient = sign_value * a**normalized_order * rgamma(normalized_order + 1.0_dp)
        value = coefficient / (effective_exponent + 1.0_dp)

        residual_exponent = effective_exponent + 2.0_dp
        transform_power = max(1, min(16, ceiling(5.0_dp / (residual_exponent + 1.0_dp))))
        panels = max(1, ceiling(max(a, 1.0_dp) * real(transform_power, dp) / 2.0_dp))
        do panel = 0, panels - 1
            left = real(panel, dp) / real(panels, dp)
            right = real(panel + 1, dp) / real(panels, dp)
            midpoint = 0.5_dp * (left + right)
            half = 0.5_dp * (right - left)
            call besselpoly_residual_panel(midpoint, half, transform_power, a, normalized_order, &
                effective_exponent, coefficient, value)
        end do
    end function besselpoly_endpoint_singular

    pure subroutine besselpoly_residual_panel(midpoint, half, transform_power, a, order, &
            endpoint_power, coefficient, total)
        real(dp), intent(in) :: midpoint !! midpoint of a transformed quadrature panel
        real(dp), intent(in) :: half !! half-width of a transformed quadrature panel
        integer, intent(in) :: transform_power !! x=t**transform_power smoothing power
        real(dp), intent(in) :: a !! positive finite Bessel scale
        real(dp), intent(in) :: order !! normalized non-pole Bessel order
        real(dp), intent(in) :: endpoint_power !! leading x exponent of the original integrand
        real(dp), intent(in) :: coefficient !! leading Bessel-series coefficient
        real(dp), intent(inout) :: total !! accumulated integral value

        integer :: i
        real(dp) :: t1
        real(dp) :: t2

        do i = 1, size(gl16_nodes)
            t1 = midpoint - half * gl16_nodes(i)
            t2 = midpoint + half * gl16_nodes(i)
            total = total + half * gl16_weights(i) * &
                (besselpoly_residual_integrand(t1, transform_power, a, order, endpoint_power, coefficient) + &
                 besselpoly_residual_integrand(t2, transform_power, a, order, endpoint_power, coefficient))
        end do
    end subroutine besselpoly_residual_panel

    pure function besselpoly_residual_integrand(t, transform_power, a, order, endpoint_power, &
            coefficient) result(value)
        real(dp), intent(in) :: t !! transformed integration coordinate in (0,1)
        integer, intent(in) :: transform_power !! x=t**transform_power smoothing power
        real(dp), intent(in) :: a !! positive finite Bessel scale
        real(dp), intent(in) :: order !! normalized non-pole Bessel order
        real(dp), intent(in) :: endpoint_power !! leading x exponent of the original integrand
        real(dp), intent(in) :: coefficient !! leading Bessel-series coefficient
        real(dp) :: value

        real(dp) :: lambda_value
        real(dp) :: x

        x = t**transform_power
        lambda_value = hyp0f1(order + 1.0_dp, -(a * x)**2)
        value = real(transform_power, dp) * t**(transform_power - 1) * coefficient * &
            x**endpoint_power * (lambda_value - 1.0_dp)
    end function besselpoly_residual_integrand

    pure subroutine besselpoly_panel_pair(midpoint, half, transform_power, a, lmb, nu, total)
        real(dp), intent(in) :: midpoint !! midpoint of a Gauss-Legendre panel in transformed t
        real(dp), intent(in) :: half !! half-width of a Gauss-Legendre panel in transformed t
        integer, intent(in) :: transform_power !! x=t**transform_power endpoint-smoothing power
        real(dp), intent(in) :: a !! nonnegative Bessel scale factor
        real(dp), intent(in) :: lmb !! power on the original integration variable
        real(dp), intent(in) :: nu !! Bessel order
        real(dp), intent(inout) :: total !! accumulated quadrature value

        integer :: i
        real(dp) :: t1
        real(dp) :: t2

        do i = 1, size(gl16_nodes)
            t1 = midpoint - half * gl16_nodes(i)
            t2 = midpoint + half * gl16_nodes(i)
            total = total + half * gl16_weights(i) * &
                (besselpoly_transformed_integrand(t1, transform_power, a, lmb, nu) + &
                 besselpoly_transformed_integrand(t2, transform_power, a, lmb, nu))
        end do
    end subroutine besselpoly_panel_pair

    pure function besselpoly_transformed_integrand(t, transform_power, a, lmb, nu) result(value)
        real(dp), intent(in) :: t !! transformed integration coordinate in (0,1)
        integer, intent(in) :: transform_power !! x=t**transform_power endpoint-smoothing power
        real(dp), intent(in) :: a !! nonnegative Bessel scale factor
        real(dp), intent(in) :: lmb !! power on the original integration variable
        real(dp), intent(in) :: nu !! Bessel order
        real(dp) :: value

        real(dp) :: x

        x = t**transform_power
        value = real(transform_power, dp) * t**(transform_power - 1) * x**lmb * &
            jv(nu, 2.0_dp * a * x)
    end function besselpoly_transformed_integrand

    pure elemental function is_negative_integer(x) result(is_integer)
        real(dp), intent(in) :: x !! finite real value to classify
        logical :: is_integer

        is_integer = x < 0.0_dp .and. x == real(nint(x), dp)
    end function is_negative_integer

    pure elemental function divergent_besselpoly_sign(nu) result(sign_value)
        real(dp), intent(in) :: nu !! finite Bessel order
        real(dp) :: sign_value

        integer :: n

        if (is_negative_integer(nu)) then
            n = nint(-nu)
            if (mod(n, 2) == 0) then
                sign_value = 1.0_dp
            else
                sign_value = -1.0_dp
            end if
        else
            ! The leading coefficient is a**nu/Gamma(nu+1). For the supported
            ! nonnegative a, only the reciprocal-Gamma sign is needed here.
            sign_value = reciprocal_gamma_sign(nu + 1.0_dp)
        end if
    end function divergent_besselpoly_sign

    pure elemental function reciprocal_gamma_sign(x) result(sign_value)
        real(dp), intent(in) :: x !! finite noninteger Gamma argument
        real(dp) :: sign_value

        integer :: interval

        if (x > 0.0_dp) then
            sign_value = 1.0_dp
        else
            interval = floor(-x)
            if (mod(interval, 2) == 0) then
                sign_value = -1.0_dp
            else
                sign_value = 1.0_dp
            end if
        end if
    end function reciprocal_gamma_sign

end module scifort_bessel_misc
