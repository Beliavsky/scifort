! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Lambert W and Wright omega functions.
!
! This is an independent implementation based on the defining equations,
! branch structure, branch-point expansion, and large-argument asymptotics in
! NIST DLMF 4.13. Iterative refinement uses Halley's method.

module scifort_lambert_w
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    public :: lambertw
    public :: wrightomega

    interface lambertw
        module procedure lambertw_real_argument
        module procedure lambertw_complex_argument
    end interface lambertw

    interface wrightomega
        module procedure wrightomega_real_argument
        module procedure wrightomega_complex_argument
    end interface wrightomega

    integer, parameter :: max_iterations = 100
    real(dp), parameter :: branch_point_radius = 0.5_dp
    real(dp), parameter :: direct_exp_limit = 680.0_dp

contains

    pure elemental function lambertw_real_argument(z, k, tol) result(value)
        real(dp), intent(in) :: z !! real argument
        integer, intent(in), optional :: k !! integer branch index; default 0
        real(dp), intent(in), optional :: tol !! positive convergence tolerance; default machine scale
        complex(dp) :: value

        value = lambertw_complex_argument(cmplx(z, 0.0_dp, dp), k, tol)
    end function lambertw_real_argument

    pure elemental function lambertw_complex_argument(z, k, tol) result(value)
        complex(dp), intent(in) :: z !! complex argument
        integer, intent(in), optional :: k !! integer branch index; default 0
        real(dp), intent(in), optional :: tol !! positive convergence tolerance; default machine scale
        complex(dp) :: value

        complex(dp) :: delta
        complex(dp) :: denominator
        complex(dp) :: ew
        complex(dp) :: f
        complex(dp) :: p
        complex(dp) :: q
        complex(dp) :: w
        complex(dp) :: wp1
        integer :: branch
        integer :: iteration
        real(dp) :: convergence_tolerance
        real(dp) :: residual_scale

        branch = 0
        if (present(k)) branch = k
        convergence_tolerance = 8.0_dp * epsilon(1.0_dp)
        if (present(tol)) then
            if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp) then
                value = complex_nan(z)
                return
            end if
            convergence_tolerance = max(convergence_tolerance, tol)
        end if

        if (complex_is_nan(z)) then
            value = complex_nan(z)
            return
        end if
        if (.not. complex_is_finite(z)) then
            value = lambertw_infinite_argument(z, branch)
            return
        end if
        if (z == cmplx(0.0_dp, 0.0_dp, dp)) then
            if (branch == 0) then
                value = cmplx(0.0_dp, 0.0_dp, dp)
            else
                value = cmplx(negative_infinity(real(z, dp)), 0.0_dp, dp)
            end if
            return
        end if

        q = exp(1.0_dp) * z + cmplx(1.0_dp, 0.0_dp, dp)
        if (branch == 0 .and. abs(q) < branch_point_radius) then
            p = sqrt(2.0_dp * q)
            w = branch_point_start(p)
        else if (branch == -1 .and. aimag(z) == 0.0_dp .and. &
                real(z, dp) >= -exp(-1.0_dp) .and. real(z, dp) < 0.0_dp) then
            p = cmplx(-sqrt(max(0.0_dp, 2.0_dp * real(q, dp))), 0.0_dp, dp)
            w = branch_point_start(p)
        else if (branch == 0 .and. abs(z) < 3.0_dp) then
            if (abs(cmplx(1.0_dp, 0.0_dp, dp) + z) > 0.2_dp) then
                w = log(cmplx(1.0_dp, 0.0_dp, dp) + z)
            else
                w = asymptotic_start(z, branch)
            end if
        else
            w = asymptotic_start(z, branch)
        end if

        residual_scale = max(abs(z), tiny(1.0_dp))
        do iteration = 1, max_iterations
            ew = exp(w)
            f = w * ew - z
            wp1 = w + cmplx(1.0_dp, 0.0_dp, dp)
            if (abs(wp1) < 1.0e-10_dp) then
                denominator = ew * wp1
            else
                denominator = ew * wp1 - (w + cmplx(2.0_dp, 0.0_dp, dp)) * f / (2.0_dp * wp1)
            end if
            if (denominator == cmplx(0.0_dp, 0.0_dp, dp)) exit

            delta = f / denominator
            w = w - delta
            if (abs(delta) <= convergence_tolerance * (1.0_dp + abs(w))) then
                value = w
                return
            end if
        end do

        ew = exp(w)
        if (abs(w * ew - z) <= 64.0_dp * convergence_tolerance * residual_scale) then
            value = w
        else
            value = complex_nan(z)
        end if
    end function lambertw_complex_argument

    pure elemental function wrightomega_real_argument(z) result(value)
        real(dp), intent(in) :: z !! real argument
        real(dp) :: value

        complex(dp) :: complex_value

        if (ieee_is_nan(z)) then
            value = quiet_nan(z)
            return
        end if
        if (.not. ieee_is_finite(z)) then
            if (z > 0.0_dp) then
                value = positive_infinity(z)
            else
                value = 0.0_dp
            end if
            return
        end if

        complex_value = wrightomega_complex_argument(cmplx(z, 0.0_dp, dp))
        value = real(complex_value, dp)
    end function wrightomega_real_argument

    pure elemental function wrightomega_complex_argument(z) result(value)
        complex(dp), intent(in) :: z !! complex argument
        complex(dp) :: value

        complex(dp) :: delta
        complex(dp) :: denominator
        complex(dp) :: f
        complex(dp) :: fp
        complex(dp) :: fpp
        complex(dp) :: w
        integer :: branch
        integer :: iteration
        real(dp) :: tolerance

        if (complex_is_nan(z)) then
            value = complex_nan(z)
            return
        end if
        if (.not. complex_is_finite(z)) then
            value = wrightomega_infinite_argument(z)
            return
        end if

        tolerance = 8.0_dp * epsilon(1.0_dp)
        if (is_wright_branch_point(z)) then
            value = cmplx(-1.0_dp, 0.0_dp, dp)
            return
        end if

        branch = wright_unwinding_number(aimag(z))
        if (real(z, dp) > -direct_exp_limit .and. real(z, dp) < direct_exp_limit) then
            value = lambertw_complex_argument(exp(z), branch)
            return
        end if

        if (real(z, dp) <= -direct_exp_limit .and. branch == 0) then
            value = exp(z)
            return
        end if

        w = z - log(z)
        do iteration = 1, max_iterations
            f = w + log(w) - z
            if (abs(f) <= tolerance * max(1.0_dp, abs(z))) then
                value = w
                return
            end if
            fp = cmplx(1.0_dp, 0.0_dp, dp) + 1.0_dp / w
            fpp = -1.0_dp / (w * w)
            denominator = 2.0_dp * fp * fp - f * fpp
            if (denominator == cmplx(0.0_dp, 0.0_dp, dp)) exit
            delta = 2.0_dp * f * fp / denominator
            w = w - delta
            if (abs(delta) <= tolerance * (1.0_dp + abs(w))) then
                value = w
                return
            end if
        end do

        if (abs(w + log(w) - z) <= 64.0_dp * tolerance * max(1.0_dp, abs(z))) then
            value = w
        else
            value = complex_nan(z)
        end if
    end function wrightomega_complex_argument

    pure function branch_point_start(p) result(w)
        complex(dp), intent(in) :: p !! square-root branch variable sqrt(2*(e*z+1))
        complex(dp) :: w

        w = -1.0_dp + p - p * p / 3.0_dp + 11.0_dp * p**3 / 72.0_dp - &
            43.0_dp * p**4 / 540.0_dp + 769.0_dp * p**5 / 17280.0_dp
    end function branch_point_start

    pure function asymptotic_start(z, branch) result(w)
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        integer, intent(in) :: branch !! integer Lambert-W branch index
        complex(dp) :: w

        complex(dp) :: l1
        complex(dp) :: l2

        l1 = log(z) + cmplx(0.0_dp, 2.0_dp * scifort_pi * real(branch, dp), dp)
        if (abs(l1) < 0.5_dp) then
            w = z
            return
        end if
        l2 = log(l1)
        w = l1 - l2 + l2 / l1 + l2 * (-2.0_dp + l2) / (2.0_dp * l1 * l1)
    end function asymptotic_start

    pure function lambertw_infinite_argument(z, branch) result(value)
        complex(dp), intent(in) :: z !! nonfinite complex argument
        integer, intent(in) :: branch !! integer branch index
        complex(dp) :: value

        real(dp) :: inf

        inf = positive_infinity(0.0_dp)
        if (ieee_is_finite(aimag(z)) .and. aimag(z) == 0.0_dp .and. &
                .not. ieee_is_finite(real(z, dp))) then
            if (real(z, dp) > 0.0_dp) then
                value = cmplx(inf, 2.0_dp * scifort_pi * real(branch, dp), dp)
            else
                value = cmplx(inf, (2.0_dp * real(branch, dp) + 1.0_dp) * scifort_pi, dp)
            end if
        else
            value = complex_nan(z)
        end if
    end function lambertw_infinite_argument

    pure function wrightomega_infinite_argument(z) result(value)
        complex(dp), intent(in) :: z !! nonfinite complex argument
        complex(dp) :: value

        if (ieee_is_finite(aimag(z)) .and. aimag(z) == 0.0_dp .and. &
                .not. ieee_is_finite(real(z, dp))) then
            if (real(z, dp) > 0.0_dp) then
                value = cmplx(positive_infinity(0.0_dp), 0.0_dp, dp)
            else
                value = cmplx(0.0_dp, 0.0_dp, dp)
            end if
        else
            value = complex_nan(z)
        end if
    end function wrightomega_infinite_argument

    pure function wright_unwinding_number(y) result(branch)
        real(dp), intent(in) :: y !! finite imaginary part of the Wright-omega argument
        integer :: branch

        real(dp) :: scaled

        scaled = (y - scifort_pi) / (2.0_dp * scifort_pi)
        if (abs(scaled) > real(huge(branch), dp) - 1.0_dp) then
            branch = 0
        else
            branch = ceiling(scaled)
        end if
    end function wright_unwinding_number

    pure function is_wright_branch_point(z) result(is_branch_point)
        complex(dp), intent(in) :: z !! finite complex argument
        logical :: is_branch_point

        real(dp) :: scale
        real(dp) :: tolerance

        scale = 1.0_dp + abs(z)
        tolerance = 32.0_dp * epsilon(1.0_dp) * scale
        is_branch_point = abs(real(z, dp) + 1.0_dp) <= tolerance .and. &
            abs(abs(aimag(z)) - scifort_pi) <= tolerance
    end function is_wright_branch_point

    pure function complex_is_nan(z) result(is_nan)
        complex(dp), intent(in) :: z !! complex value to inspect
        logical :: is_nan

        is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_is_nan

    pure function complex_is_finite(z) result(is_finite)
        complex(dp), intent(in) :: z !! complex value to inspect
        logical :: is_finite

        is_finite = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure function complex_nan(z) result(value)
        complex(dp), intent(in) :: z !! kind carrier for NaN result
        complex(dp) :: value

        real(dp) :: nan_value

        nan_value = quiet_nan(real(z, dp))
        value = cmplx(nan_value, nan_value, dp)
    end function complex_nan

end module scifort_lambert_w
