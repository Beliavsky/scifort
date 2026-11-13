! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2025, Alexey-Kuznetsov-math
! Copyright (c) 2026 SciFort contributors
!
! Euler-Maclaurin kernel and Bernoulli coefficients adapted from Alexey
! Kuznetsov's Riemann_zeta-and-Gamma, Fortran/gamma_zeta_module.f90,
! commit 1f1ad8934044cfb20d6d70137868ddd0c74652c4.
! Full BSD terms: LICENSES/kuznetsov-zeta-BSD-3-Clause.txt.
! Binary64 adaptation uses compensated summation, bounded work, and SciFort
! log-gamma. Continuation follows DLMF 25.4.2; summation follows 25.2.10.
module scifort_complex_zeta
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_math, only : quiet_nan
    use scifort_gamma_family, only : loggamma
    implicit none
    private
    public :: complex_zeta_value

    ! B_(2k)/(2k)!; checked against high-precision mpmath Bernoulli numbers.
    real(dp), parameter :: bernoulli_coeff(25) = [ &
        8.333333333333333333333333333333333333e-2_dp, &
        -1.388888888888888888888888888888888889e-3_dp, &
        3.306878306878306878306878306878306878e-5_dp, &
        -8.267195767195767195767195767195767196e-7_dp, &
        2.087675698786809897921009032120143231e-8_dp, &
        -5.284190138687493184847682202179556677e-10_dp, &
        1.338253653068467883282698097512912328e-11_dp, &
        -3.389680296322582866830195391249442500e-13_dp, &
        8.586062056277844564135905450425627134e-15_dp, &
        -2.174868698558061873041516423865917900e-16_dp, &
        5.509002828360229515202652608902254878e-18_dp, &
        -1.395446468581252334070768626406354976e-19_dp, &
        3.534707039629467471693229977803799215e-21_dp, &
        -8.953517427037546850402611318112741052e-23_dp, &
        2.267952452337683060310950738868166063e-24_dp, &
        -5.744790668872202445263881987607018400e-26_dp, &
        1.455172475614864901866264867271329336e-27_dp, &
        -3.685994940665310178181782479908660374e-29_dp, &
        9.336734257095044672032555152785623295e-31_dp, &
        -2.365022415700629934559635196369838240e-32_dp, &
        5.990671762482134304659912396819657826e-34_dp, &
        -1.517454884468290261710813135864718932e-35_dp, &
        3.843758125454188232229445290990232106e-37_dp, &
        -9.736353072646691035267621279250454181e-39_dp, &
        2.466247044200680957106400280288842886e-40_dp]

contains

    pure elemental function complex_zeta_value(s) result(value)
        implicit none
        complex(dp), intent(in) :: s !! Nonreal exponent; |Re(s)|<=1e4, |Im(s)|<=1e5.
        complex(dp) :: value, logarithm, reflected
        real(dp) :: nan_value, magnitude

        nan_value = quiet_nan(0.0_dp)
        value = cmplx(nan_value,nan_value,dp)
        if (.not. ieee_is_finite(real(s,dp)) .or. .not. ieee_is_finite(aimag(s))) return
        if (abs(real(s,dp)) > 1.0e4_dp .or. abs(aimag(s)) > 1.0e5_dp) return
        ! Guard the reciprocal pole before any complex division. This small
        ! excluded disk also keeps intermediate sums away from overflow.
        if (abs(s-1.0_dp) < 16.0_dp/huge(1.0_dp)) return
        if (real(s,dp) >= -0.25_dp) then
            value = euler_maclaurin(s)
        else
            reflected = euler_maclaurin(1.0_dp-s)
            if (abs(reflected) == 0.0_dp) return
            logarithm = log(2.0_dp)+(s-1.0_dp)*log(2.0_dp*scifort_pi) &
                + log_sin_half_pi(s)+loggamma(1.0_dp-s)+log(reflected)
            magnitude = real(logarithm,dp)
            if (magnitude > log(huge(1.0_dp))-1.0_dp) return
            if (magnitude < log(tiny(1.0_dp))+log(epsilon(1.0_dp))) then
                value = cmplx(0.0_dp,0.0_dp,dp)
            else
                value = exp(logarithm)
            end if
        end if
    end function complex_zeta_value

    pure function euler_maclaurin(s) result(value)
        implicit none
        complex(dp), intent(in) :: s !! Exponent with Re(s)>=-0.25, away from the pole.
        complex(dp) :: value, rising, correction, term, compensation, adjusted, next
        integer :: n, k
        real(dp) :: rn, logn

        n = int(abs(s))+32
        rn = real(n,dp)
        logn = log(rn)
        rising = s/rn
        correction = cmplx(0.0_dp,0.0_dp,dp)
        do k = 1, size(bernoulli_coeff)
            correction = correction+bernoulli_coeff(k)*rising
            rising = rising*((s+real(2*k-1,dp))/rn)*((s+real(2*k,dp))/rn)
        end do
        value = cmplx(0.0_dp,0.0_dp,dp)
        compensation = value
        do k = 1, n-1
            term = exp(-s*log(real(k,dp)))
            adjusted = term-compensation
            next = value+adjusted
            compensation = (next-value)-adjusted
            value = next
        end do
        value = value+exp(-s*logn)*(correction+0.5_dp) &
            +exp((1.0_dp-s)*logn)/(s-1.0_dp)
    end function euler_maclaurin

    pure function log_sin_half_pi(s) result(value)
        implicit none
        complex(dp), intent(in) :: s !! Nonreal bounded exponent in the reflected half-plane.
        complex(dp) :: value
        real(dp) :: r, x, y, sn, cs

        r = modulo(real(s,dp),4.0_dp)
        x = (0.5_dp*scifort_pi)*r
        y = (0.5_dp*scifort_pi)*aimag(s)
        if (abs(y) >= 20.0_dp) then
            ! The omitted log(1-exp(-2|y|+i phase)) is below 4.3e-18.
            value = cmplx(abs(y)-log(2.0_dp), &
                sign(1.0_dp,y)*(0.5_dp*scifort_pi-x),dp)
            return
        end if
        if (r == 0.0_dp .or. r == 2.0_dp) then
            sn = 0.0_dp
        else if (r <= 1.0_dp) then
            sn = sin(0.5_dp*scifort_pi*r)
        else if (r < 2.0_dp) then
            sn = sin(0.5_dp*scifort_pi*(2.0_dp-r))
        else if (r <= 3.0_dp) then
            sn = -sin(0.5_dp*scifort_pi*(r-2.0_dp))
        else
            sn = -sin(0.5_dp*scifort_pi*(4.0_dp-r))
        end if
        if (r == 1.0_dp .or. r == 3.0_dp) then
            cs = 0.0_dp
        else if (r == 0.0_dp) then
            cs = 1.0_dp
        else if (r == 2.0_dp) then
            cs = -1.0_dp
        else
            cs = cos(x)
        end if
        value = log(cmplx(sn*cosh(y),cs*sinh(y),dp))
    end function log_sin_half_pi
end module scifort_complex_zeta
