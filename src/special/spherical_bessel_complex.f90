! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Complex spherical Bessel functions from the half-integer cylindrical
! definitions in DLMF 10.47.3-10.47.9. The underlying complex J/Y/I/K
! functions provide principal-branch continuation away from the real axis.
! High-order spherical J near its turning region uses backward Miller
! recurrence from the DLMF 10.51 identities to avoid loss of accuracy.
module scifort_spherical_bessel_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_bessel_complex, only : iv_complex, jv_complex, kv_complex, yv_complex
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_spherical_bessel, only : spherical_in_real => spherical_in, &
        spherical_jn_real => spherical_jn, spherical_kn_real => spherical_kn, &
        spherical_yn_real => spherical_yn
    implicit none
    private

    public :: spherical_in_complex
    public :: spherical_jn_complex
    public :: spherical_kn_complex
    public :: spherical_yn_complex

contains

    pure elemental function spherical_jn_complex(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        logical :: want_derivative

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (invalid_inputs(n, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = spherical_j_or_i_at_zero(n, want_derivative)
        else if (aimag(z) == 0.0_dp) then
            value = cmplx(spherical_jn_real(n, real(z, dp), want_derivative), 0.0_dp, kind=dp)
        else if (.not. want_derivative) then
            value = spherical_jn_value(n, z)
        else if (n == 0) then
            value = -spherical_jn_value(1, z)
        else
            value = spherical_jn_value(n - 1, z) - &
                real(n + 1, dp) * spherical_jn_value(n, z) / z
        end if
    end function spherical_jn_complex

    pure elemental function spherical_yn_complex(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        logical :: want_derivative

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (invalid_inputs(n, z) .or. complex_is_zero(z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = cmplx(spherical_yn_real(n, real(z, dp), want_derivative), 0.0_dp, kind=dp)
        else if (.not. want_derivative) then
            value = spherical_yn_value(n, z)
        else if (n == 0) then
            value = -spherical_yn_value(1, z)
        else
            value = spherical_yn_value(n - 1, z) - &
                real(n + 1, dp) * spherical_yn_value(n, z) / z
        end if
    end function spherical_yn_complex

    pure elemental function spherical_in_complex(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        logical :: want_derivative

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (invalid_inputs(n, z)) then
            value = complex_nan()
        else if (complex_is_zero(z)) then
            value = spherical_j_or_i_at_zero(n, want_derivative)
        else if (aimag(z) == 0.0_dp) then
            value = cmplx(spherical_in_real(n, real(z, dp), want_derivative), 0.0_dp, kind=dp)
        else if (.not. want_derivative) then
            value = spherical_in_value(n, z)
        else if (n == 0) then
            value = spherical_in_value(1, z)
        else
            value = spherical_in_value(n - 1, z) - &
                real(n + 1, dp) * spherical_in_value(n, z) / z
        end if
    end function spherical_in_complex

    pure elemental function spherical_kn_complex(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        logical :: want_derivative

        want_derivative = .false.
        if (present(derivative)) want_derivative = derivative

        if (invalid_inputs(n, z) .or. complex_is_zero(z)) then
            value = complex_nan()
        else if (aimag(z) == 0.0_dp) then
            value = cmplx(spherical_kn_real(n, real(z, dp), want_derivative), 0.0_dp, kind=dp)
        else if (.not. want_derivative) then
            value = spherical_kn_value(n, z)
        else if (n == 0) then
            value = -spherical_kn_value(1, z)
        else
            value = -spherical_kn_value(n - 1, z) - &
                real(n + 1, dp) * spherical_kn_value(n, z) / z
        end if
    end function spherical_kn_complex

    pure elemental function spherical_jn_value(n, z) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero nonreal argument
        complex(dp) :: value

        if (n >= 8 .and. abs(z) <= 1.1_dp * real(n + 1, dp)) then
            value = spherical_jn_miller_complex(n, z)
        else
            value = spherical_prefactor(z) * jv_complex(real(n, dp) + 0.5_dp, z)
        end if
    end function spherical_jn_value

    pure function spherical_jn_miller_complex(n, z) result(value)
        integer, intent(in) :: n !! order >= 8 in the complex Miller regime
        complex(dp), intent(in) :: z !! finite nonzero nonreal argument
        complex(dp) :: value

        integer :: extra
        integer :: k
        integer :: m
        logical :: have_f1
        logical :: have_fn
        real(dp), parameter :: rescale = 1.0e-100_dp
        real(dp), parameter :: rescale_trigger = 1.0e100_dp
        complex(dp) :: b0
        complex(dp) :: b1
        complex(dp) :: f0
        complex(dp) :: f1
        complex(dp) :: fk
        complex(dp) :: fkm1
        complex(dp) :: fkp1
        complex(dp) :: fn

        extra = max(50, int(2.0_dp * sqrt(40.0_dp * real(n + 1, dp))) + 10)
        if (n > huge(n) - extra) then
            value = spherical_prefactor(z) * jv_complex(real(n, dp) + 0.5_dp, z)
            return
        end if

        m = n + extra
        fkp1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        fk = cmplx(1.0_dp, 0.0_dp, kind=dp)
        fn = cmplx(0.0_dp, 0.0_dp, kind=dp)
        f1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        have_fn = .false.
        have_f1 = .false.

        do k = m, 1, -1
            fkm1 = (real(2 * k + 1, dp) / z) * fk - fkp1
            if (k - 1 == n) then
                fn = fkm1
                have_fn = .true.
            end if
            if (k - 1 == 1) then
                f1 = fkm1
                have_f1 = .true.
            end if
            if (abs(fkm1) > rescale_trigger) then
                fkm1 = rescale * fkm1
                fk = rescale * fk
                if (have_fn) fn = rescale * fn
                if (have_f1) f1 = rescale * f1
            end if
            fkp1 = fk
            fk = fkm1
        end do

        f0 = fk
        b0 = spherical_prefactor(z) * jv_complex(0.5_dp, z)
        b1 = spherical_prefactor(z) * jv_complex(1.5_dp, z)
        if (abs(b0) >= abs(b1)) then
            value = fn * (b0 / f0)
        else
            value = fn * (b1 / f1)
        end if
    end function spherical_jn_miller_complex

    pure elemental function spherical_yn_value(n, z) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero nonreal argument
        complex(dp) :: value

        value = spherical_prefactor(z) * yv_complex(real(n, dp) + 0.5_dp, z)
    end function spherical_yn_value

    pure elemental function spherical_in_value(n, z) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero nonreal argument
        complex(dp) :: value

        value = spherical_prefactor(z) * iv_complex(real(n, dp) + 0.5_dp, z)
    end function spherical_in_value

    pure elemental function spherical_kn_value(n, z) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero nonreal argument
        complex(dp) :: value

        value = spherical_prefactor(z) * kv_complex(real(n, dp) + 0.5_dp, z)
    end function spherical_kn_value

    pure elemental function spherical_prefactor(z) result(value)
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        complex(dp) :: value

        value = sqrt(scifort_pi / (2.0_dp * z))
    end function spherical_prefactor

    pure elemental function spherical_j_or_i_at_zero(n, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        logical, intent(in) :: derivative !! whether the first derivative is requested
        complex(dp) :: value

        if (derivative) then
            if (n == 1) then
                value = cmplx(1.0_dp / 3.0_dp, 0.0_dp, kind=dp)
            else
                value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            end if
        else if (n == 0) then
            value = cmplx(1.0_dp, 0.0_dp, kind=dp)
        else
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end if
    end function spherical_j_or_i_at_zero

    pure elemental logical function invalid_inputs(n, z) result(invalid)
        integer, intent(in) :: n !! spherical-Bessel order to classify
        complex(dp), intent(in) :: z !! complex argument to classify

        invalid = n < 0 .or. .not. ieee_is_finite(real(z, dp)) .or. &
            .not. ieee_is_finite(aimag(z))
    end function invalid_inputs

    pure elemental logical function complex_is_zero(z) result(is_zero)
        complex(dp), intent(in) :: z !! complex value to compare with zero

        is_zero = real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp
    end function complex_is_zero

    pure elemental function complex_nan() result(value)
        complex(dp) :: value
        real(dp) :: nan_value

        nan_value = quiet_nan(0.0_dp)
        value = cmplx(nan_value, nan_value, kind=dp)
    end function complex_nan

end module scifort_spherical_bessel_complex
