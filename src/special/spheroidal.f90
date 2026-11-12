! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Prolate and oblate spheroidal characteristic values plus angular and radial
! functions. This is original SciFort code based on the finite
! tridiagonal eigenvalue construction, differential equations, and
! Ferrers/spherical-Bessel expansions in NIST DLMF sections 30.2, 30.8, 30.11,
! 30.14, and 30.16. No third-party implementation source is copied or adapted.

module scifort_spheroidal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_spherical_bessel, only : spherical_jn, spherical_yn
    implicit none
    private

    integer, parameter, public :: spheroidal_status_success = 0
    integer, parameter, public :: spheroidal_status_invalid_argument = 1
    integer, parameter, public :: spheroidal_status_numerical_failure = 2

    integer, parameter :: max_matrix_size = 256
    integer, parameter :: max_bisection_iterations = 256
    real(dp), parameter :: convergence_tolerance = 5.0e-14_dp
    real(dp), parameter :: angular_convergence_tolerance = 2.0e-13_dp
    real(dp), parameter :: radial_convergence_tolerance = 5.0e-13_dp
    real(dp), parameter :: radial2_convergence_tolerance = 2.0e-12_dp
    real(dp), parameter :: radial2_ode_tolerance = 5.0e-14_dp

    public :: obl_ang1
    public :: obl_ang1_cv
    public :: obl_cv
    public :: obl_cv_seq
    public :: obl_rad1
    public :: obl_rad1_cv
    public :: obl_rad2
    public :: obl_rad2_cv
    public :: pro_ang1
    public :: pro_ang1_cv
    public :: pro_cv
    public :: pro_cv_seq
    public :: pro_rad1
    public :: pro_rad1_cv
    public :: pro_rad2
    public :: pro_rad2_cv

contains

    pure elemental subroutine pro_ang1(m, n, c, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real prolate spheroidal parameter
        real(dp), intent(in) :: x !! angular argument, abs(x) < 1
        real(dp), intent(out) :: value !! prolate angular function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: characteristic

        characteristic = pro_cv(m, n, c)
        call pro_ang1_cv(m, n, c, characteristic, x, value, derivative)
    end subroutine pro_ang1

    pure elemental subroutine obl_ang1(m, n, c, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real oblate spheroidal parameter
        real(dp), intent(in) :: x !! angular argument, abs(x) < 1
        real(dp), intent(out) :: value !! oblate angular function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: characteristic

        characteristic = obl_cv(m, n, c)
        call obl_ang1_cv(m, n, c, characteristic, x, value, derivative)
    end subroutine obl_ang1

    pure elemental subroutine pro_ang1_cv(m, n, c, cv, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real prolate spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! angular argument, abs(x) < 1
        real(dp), intent(out) :: value !! prolate angular function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        call spheroidal_ang1(m, n, c, cv, x, 1.0_dp, value, derivative)
    end subroutine pro_ang1_cv

    pure elemental subroutine obl_ang1_cv(m, n, c, cv, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real oblate spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! angular argument, abs(x) < 1
        real(dp), intent(out) :: value !! oblate angular function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        call spheroidal_ang1(m, n, c, cv, x, -1.0_dp, value, derivative)
    end subroutine obl_ang1_cv

    pure elemental subroutine pro_rad1(m, n, c, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! nonnegative finite prolate spheroidal parameter
        real(dp), intent(in) :: x !! radial argument, x > 1
        real(dp), intent(out) :: value !! prolate radial function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: characteristic

        if (m < 0 .or. n < m .or. c < 0.0_dp .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(x) .or. x <= 1.0_dp) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        characteristic = pro_cv(m, n, c)
        call pro_rad1_cv(m, n, c, characteristic, x, value, derivative)
    end subroutine pro_rad1

    pure elemental subroutine obl_rad1(m, n, c, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! nonnegative finite oblate spheroidal parameter
        real(dp), intent(in) :: x !! radial argument, x >= 0
        real(dp), intent(out) :: value !! oblate radial function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: characteristic

        if (m < 0 .or. n < m .or. c < 0.0_dp .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        characteristic = obl_cv(m, n, c)
        call obl_rad1_cv(m, n, c, characteristic, x, value, derivative)
    end subroutine obl_rad1

    pure elemental subroutine pro_rad1_cv(m, n, c, cv, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! nonnegative finite prolate spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument, x > 1
        real(dp), intent(out) :: value !! prolate radial function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        call spheroidal_rad1(m, n, c, cv, x, 1.0_dp, .false., value, derivative)
    end subroutine pro_rad1_cv

    pure elemental subroutine obl_rad1_cv(m, n, c, cv, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! nonnegative finite oblate spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument, x >= 0
        real(dp), intent(out) :: value !! oblate radial function of the first kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        call spheroidal_rad1(m, n, c, cv, x, -1.0_dp, .true., value, derivative)
    end subroutine obl_rad1_cv

    pure elemental subroutine pro_rad2(m, n, c, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite prolate spheroidal parameter
        real(dp), intent(in) :: x !! radial argument, x > 1
        real(dp), intent(out) :: value !! prolate radial function of the second kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: characteristic

        if (m < 0 .or. n < m .or. c <= 0.0_dp .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(x) .or. x <= 1.0_dp) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        characteristic = pro_cv(m, n, c)
        call pro_rad2_cv(m, n, c, characteristic, x, value, derivative)
    end subroutine pro_rad2

    pure elemental subroutine obl_rad2(m, n, c, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite oblate spheroidal parameter
        real(dp), intent(in) :: x !! radial argument, x >= 0
        real(dp), intent(out) :: value !! oblate radial function of the second kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: characteristic

        if (m < 0 .or. n < m .or. c <= 0.0_dp .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(x) .or. x < 0.0_dp) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        characteristic = obl_cv(m, n, c)
        call obl_rad2_cv(m, n, c, characteristic, x, value, derivative)
    end subroutine obl_rad2

    pure elemental subroutine pro_rad2_cv(m, n, c, cv, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite prolate spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument, x > 1
        real(dp), intent(out) :: value !! prolate radial function of the second kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        call spheroidal_rad2(m, n, c, cv, x, 1.0_dp, .false., value, derivative)
    end subroutine pro_rad2_cv

    pure elemental subroutine obl_rad2_cv(m, n, c, cv, x, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite oblate spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument, x >= 0
        real(dp), intent(out) :: value !! oblate radial function of the second kind
        real(dp), intent(out) :: derivative !! derivative with respect to x

        call spheroidal_rad2(m, n, c, cv, x, -1.0_dp, .true., value, derivative)
    end subroutine obl_rad2_cv

    pure elemental function pro_cv(m, n, c) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real prolate spheroidal parameter
        real(dp) :: value

        value = spheroidal_cv(m, n, c, 1.0_dp)
    end function pro_cv

    pure elemental function obl_cv(m, n, c) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real oblate spheroidal parameter
        real(dp) :: value

        value = spheroidal_cv(m, n, c, -1.0_dp)
    end function obl_cv

    pure subroutine pro_cv_seq(m, n, c, values, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! largest spheroidal degree requested
        real(dp), intent(in) :: c !! finite real prolate spheroidal parameter
        real(dp), allocatable, intent(out) :: values(:) !! characteristic values for degrees m through n
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input or numerical failure

        call spheroidal_cv_sequence(m, n, c, 1.0_dp, values, status)
    end subroutine pro_cv_seq

    pure subroutine obl_cv_seq(m, n, c, values, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! largest spheroidal degree requested
        real(dp), intent(in) :: c !! finite real oblate spheroidal parameter
        real(dp), allocatable, intent(out) :: values(:) !! characteristic values for degrees m through n
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input or numerical failure

        call spheroidal_cv_sequence(m, n, c, -1.0_dp, values, status)
    end subroutine obl_cv_seq

    pure subroutine spheroidal_rad1(m, n, c, cv, x, gamma_sign, oblate, &
            value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! nonnegative finite spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        logical, intent(in) :: oblate !! true for the oblate radial coordinate
        real(dp), intent(out) :: value !! first-kind radial spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x

        integer :: dimension
        integer :: next_dimension
        integer :: status
        integer :: target
        real(dp) :: current_derivative
        real(dp) :: current_value
        real(dp) :: previous_derivative
        real(dp) :: previous_value
        real(dp) :: scale_derivative
        real(dp) :: scale_value
        logical :: have_previous

        if (m < 0 .or. n < m .or. c < 0.0_dp .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(cv) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        if ((oblate .and. x < 0.0_dp) .or. (.not. oblate .and. x <= 1.0_dp)) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        if (c == 0.0_dp) then
            if (m == 0 .and. n == 0) then
                value = 1.0_dp
            else
                value = 0.0_dp
            end if
            derivative = 0.0_dp
            return
        end if
        if (c > sqrt(huge(1.0_dp))) then
            value = quiet_nan(c)
            derivative = value
            return
        end if

        target = (n - m) / 2 + 1
        if (target > max_matrix_size - 8) then
            value = quiet_nan(c)
            derivative = value
            return
        end if

        dimension = max(24, target + 16)
        if (c < real(max_matrix_size - 16, dp)) then
            dimension = max(dimension, int(c) + 16)
        else
            dimension = max_matrix_size
        end if
        dimension = min(dimension, max_matrix_size)
        previous_value = 0.0_dp
        previous_derivative = 0.0_dp
        have_previous = .false.

        do
            call spheroidal_rad1_dimension(m, n, c, cv, x, gamma_sign, oblate, &
                dimension, current_value, current_derivative, status)
            if (status /= spheroidal_status_success) then
                value = quiet_nan(x)
                derivative = value
                return
            end if

            value = current_value
            derivative = current_derivative
            if (have_previous) then
                scale_value = max(1.0_dp, abs(value))
                scale_derivative = max(1.0_dp, abs(derivative))
                if (abs(value - previous_value) <= &
                        radial_convergence_tolerance * scale_value .and. &
                        abs(derivative - previous_derivative) <= &
                        radial_convergence_tolerance * scale_derivative) return
            end if

            if (dimension == max_matrix_size) then
                value = quiet_nan(x)
                derivative = value
                return
            end if

            previous_value = value
            previous_derivative = derivative
            have_previous = .true.
            next_dimension = max(dimension + 16, (3 * dimension + 1) / 2)
            dimension = min(max_matrix_size, next_dimension)
        end do
    end subroutine spheroidal_rad1

    pure subroutine spheroidal_rad1_dimension(m, n, c, cv, x, gamma_sign, oblate, &
            dimension, value, derivative, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! nonnegative finite spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        logical, intent(in) :: oblate !! true for the oblate radial coordinate
        integer, intent(in) :: dimension !! finite spherical-Bessel expansion dimension
        real(dp), intent(out) :: value !! first-kind radial spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: minus_coefficients(:)
        integer :: degree
        integer :: j
        integer :: target
        logical :: odd_family
        real(dp) :: bessel_derivative
        real(dp) :: bessel_value
        real(dp) :: correction_denominator
        real(dp) :: correction_derivative
        real(dp) :: correction_value
        real(dp) :: denominator
        real(dp) :: derivative_prefactor
        real(dp) :: prefactor
        real(dp) :: raw_derivative
        real(dp) :: raw_value
        real(dp) :: sign_factor
        real(dp) :: term
        real(dp) :: z

        target = (n - m) / 2 + 1
        odd_family = modulo(n - m, 2) == 1
        if (dimension < target + 1) then
            value = quiet_nan(x)
            derivative = value
            status = spheroidal_status_numerical_failure
            return
        end if

        call spheroidal_angular_coefficients(m, n, c, cv, gamma_sign, dimension, &
            coefficients, status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        call spheroidal_minus_order_coefficients(m, n, coefficients, minus_coefficients, &
            status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        denominator = 0.0_dp
        correction_denominator = 0.0_dp
        do j = 1, dimension
            if (modulo(j - target, 2) == 0) then
                sign_factor = 1.0_dp
            else
                sign_factor = -1.0_dp
            end if
            term = sign_factor * minus_coefficients(j)
            call compensated_add(denominator, correction_denominator, term)
        end do
        denominator = denominator + correction_denominator
        if (denominator == 0.0_dp .or. .not. ieee_is_finite(denominator)) then
            value = quiet_nan(x)
            derivative = value
            status = spheroidal_status_numerical_failure
            return
        end if

        if (oblate .and. x == 0.0_dp) then
            call oblate_rad1_origin(m, odd_family, c, minus_coefficients(1), &
                denominator, value, derivative, status)
            return
        end if

        z = c * x
        raw_value = 0.0_dp
        raw_derivative = 0.0_dp
        correction_value = 0.0_dp
        correction_derivative = 0.0_dp
        do j = 1, dimension
            if (odd_family) then
                degree = m + 2 * j - 1
            else
                degree = m + 2 * j - 2
            end if
            bessel_value = spherical_jn(degree, z)
            bessel_derivative = c * spherical_jn(degree, z, .true.)
            term = minus_coefficients(j) * bessel_value
            call compensated_add(raw_value, correction_value, term)
            term = minus_coefficients(j) * bessel_derivative
            call compensated_add(raw_derivative, correction_derivative, term)
        end do
        raw_value = (raw_value + correction_value) / denominator
        raw_derivative = (raw_derivative + correction_derivative) / denominator

        if (oblate) then
            prefactor = (sqrt(1.0_dp + x * x) / x) ** m
            if (m == 0) then
                derivative_prefactor = 0.0_dp
            else
                derivative_prefactor = -prefactor * real(m, dp) / &
                    (x * (1.0_dp + x * x))
            end if
        else
            prefactor = (sqrt((x - 1.0_dp) * (x + 1.0_dp)) / x) ** m
            if (m == 0) then
                derivative_prefactor = 0.0_dp
            else
                derivative_prefactor = prefactor * real(m, dp) / &
                    (x * (x - 1.0_dp) * (x + 1.0_dp))
            end if
        end if

        value = prefactor * raw_value
        derivative = derivative_prefactor * raw_value + prefactor * raw_derivative
        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(derivative)) then
            status = spheroidal_status_numerical_failure
            return
        end if
        status = spheroidal_status_success
    end subroutine spheroidal_rad1_dimension

    pure subroutine spheroidal_rad2(m, n, c, cv, x, gamma_sign, oblate, &
            value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        logical, intent(in) :: oblate !! true for the oblate radial coordinate
        real(dp), intent(out) :: value !! second-kind radial spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x

        integer :: target
        integer :: status
        real(dp) :: anchor
        real(dp) :: minimum_phase

        if (m < 0 .or. n < m .or. c <= 0.0_dp .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(cv) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        if ((oblate .and. x < 0.0_dp) .or. (.not. oblate .and. x <= 1.0_dp)) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        if (c > sqrt(huge(1.0_dp))) then
            value = quiet_nan(c)
            derivative = value
            return
        end if

        target = (n - m) / 2 + 1
        if (target > max_matrix_size - 8) then
            value = quiet_nan(c)
            derivative = value
            return
        end if

        minimum_phase = real(n, dp) + real(m, dp) + 12.0_dp
        if (c < minimum_phase / sqrt(huge(1.0_dp)) .or. &
                x > 0.5_dp * sqrt(huge(1.0_dp))) then
            value = quiet_nan(c)
            derivative = value
            return
        end if
        anchor = max(x, 2.0_dp, minimum_phase / c)

        call spheroidal_rad2_series(m, n, c, cv, anchor, gamma_sign, oblate, &
            value, derivative, status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        if (anchor > x) then
            call integrate_radial2_ode(m, c, cv, anchor, x, oblate, value, &
                derivative, status)
            if (status /= spheroidal_status_success) then
                value = quiet_nan(x)
                derivative = value
            end if
        end if
    end subroutine spheroidal_rad2

    pure subroutine spheroidal_rad2_series(m, n, c, cv, x, gamma_sign, oblate, &
            value, derivative, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! stable radial-series argument
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        logical, intent(in) :: oblate !! true for the oblate radial coordinate
        real(dp), intent(out) :: value !! second-kind radial spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer :: dimension
        integer :: next_dimension
        integer :: target
        real(dp) :: current_derivative
        real(dp) :: current_value
        real(dp) :: previous_derivative
        real(dp) :: previous_value
        real(dp) :: scale_derivative
        real(dp) :: scale_value
        logical :: have_previous

        target = (n - m) / 2 + 1
        dimension = max(24, target + 16)
        if (c < real(max_matrix_size - 16, dp)) then
            dimension = max(dimension, int(c) + 16)
        else
            dimension = max_matrix_size
        end if
        dimension = min(dimension, max_matrix_size)
        previous_value = 0.0_dp
        previous_derivative = 0.0_dp
        have_previous = .false.

        do
            call spheroidal_rad2_dimension(m, n, c, cv, x, gamma_sign, oblate, &
                dimension, current_value, current_derivative, status)
            if (status /= spheroidal_status_success) return

            value = current_value
            derivative = current_derivative
            if (have_previous) then
                scale_value = max(1.0_dp, abs(value))
                scale_derivative = max(1.0_dp, abs(derivative))
                if (abs(value - previous_value) <= &
                        radial2_convergence_tolerance * scale_value .and. &
                        abs(derivative - previous_derivative) <= &
                        radial2_convergence_tolerance * scale_derivative) then
                    status = spheroidal_status_success
                    return
                end if
            end if

            if (dimension == max_matrix_size) then
                status = spheroidal_status_numerical_failure
                return
            end if

            previous_value = value
            previous_derivative = derivative
            have_previous = .true.
            next_dimension = max(dimension + 12, (5 * dimension + 3) / 4)
            dimension = min(max_matrix_size, next_dimension)
        end do
    end subroutine spheroidal_rad2_series

    pure subroutine spheroidal_rad2_dimension(m, n, c, cv, x, gamma_sign, oblate, &
            dimension, value, derivative, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! positive finite spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial-series argument
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        logical, intent(in) :: oblate !! true for the oblate radial coordinate
        integer, intent(in) :: dimension !! finite spherical-Bessel expansion dimension
        real(dp), intent(out) :: value !! second-kind radial spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: minus_coefficients(:)
        integer :: degree
        integer :: j
        integer :: target
        logical :: odd_family
        real(dp) :: bessel_derivative
        real(dp) :: bessel_value
        real(dp) :: correction_denominator
        real(dp) :: correction_derivative
        real(dp) :: correction_value
        real(dp) :: denominator
        real(dp) :: derivative_prefactor
        real(dp) :: prefactor
        real(dp) :: raw_derivative
        real(dp) :: raw_value
        real(dp) :: sign_factor
        real(dp) :: term
        real(dp) :: z

        target = (n - m) / 2 + 1
        odd_family = modulo(n - m, 2) == 1
        if (dimension < target + 1) then
            value = quiet_nan(x)
            derivative = value
            status = spheroidal_status_numerical_failure
            return
        end if

        call spheroidal_angular_coefficients(m, n, c, cv, gamma_sign, dimension, &
            coefficients, status)
        if (status /= spheroidal_status_success) return
        call spheroidal_minus_order_coefficients(m, n, coefficients, minus_coefficients, &
            status)
        if (status /= spheroidal_status_success) return

        denominator = 0.0_dp
        correction_denominator = 0.0_dp
        do j = 1, dimension
            if (modulo(j - target, 2) == 0) then
                sign_factor = 1.0_dp
            else
                sign_factor = -1.0_dp
            end if
            term = sign_factor * minus_coefficients(j)
            call compensated_add(denominator, correction_denominator, term)
        end do
        denominator = denominator + correction_denominator
        if (denominator == 0.0_dp .or. .not. ieee_is_finite(denominator)) then
            status = spheroidal_status_numerical_failure
            return
        end if

        z = c * x
        raw_value = 0.0_dp
        raw_derivative = 0.0_dp
        correction_value = 0.0_dp
        correction_derivative = 0.0_dp
        do j = 1, dimension
            if (odd_family) then
                degree = m + 2 * j - 1
            else
                degree = m + 2 * j - 2
            end if
            bessel_value = spherical_yn(degree, z)
            bessel_derivative = c * spherical_yn(degree, z, .true.)
            if (.not. ieee_is_finite(bessel_value) .or. &
                    .not. ieee_is_finite(bessel_derivative)) then
                status = spheroidal_status_numerical_failure
                return
            end if
            term = minus_coefficients(j) * bessel_value
            call compensated_add(raw_value, correction_value, term)
            term = minus_coefficients(j) * bessel_derivative
            call compensated_add(raw_derivative, correction_derivative, term)
        end do
        raw_value = (raw_value + correction_value) / denominator
        raw_derivative = (raw_derivative + correction_derivative) / denominator

        if (oblate) then
            prefactor = (sqrt(1.0_dp + x * x) / x) ** m
            if (m == 0) then
                derivative_prefactor = 0.0_dp
            else
                derivative_prefactor = -prefactor * real(m, dp) / &
                    (x * (1.0_dp + x * x))
            end if
        else
            prefactor = (sqrt((x - 1.0_dp) * (x + 1.0_dp)) / x) ** m
            if (m == 0) then
                derivative_prefactor = 0.0_dp
            else
                derivative_prefactor = prefactor * real(m, dp) / &
                    (x * (x - 1.0_dp) * (x + 1.0_dp))
            end if
        end if

        value = prefactor * raw_value
        derivative = derivative_prefactor * raw_value + prefactor * raw_derivative
        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(derivative)) then
            status = spheroidal_status_numerical_failure
            return
        end if
        status = spheroidal_status_success
    end subroutine spheroidal_rad2_dimension

    pure subroutine integrate_radial2_ode(m, c, cv, x_start, x_end, oblate, &
            value, derivative, status)
        ! The embedded RK5(4) tableau is from Dormand and Prince,
        ! J. Comput. Appl. Math. 6 (1980), 19-26,
        ! doi:10.1016/0771-050X(80)90013-3.
        integer, intent(in) :: m !! nonnegative spheroidal order
        real(dp), intent(in) :: c !! positive finite spheroidal parameter
        real(dp), intent(in) :: cv !! SciPy-convention characteristic value
        real(dp), intent(in) :: x_start !! starting argument with stable series values
        real(dp), intent(in) :: x_end !! requested radial argument
        logical, intent(in) :: oblate !! true for the oblate radial equation
        real(dp), intent(inout) :: value !! value at start on entry, at end on return
        real(dp), intent(inout) :: derivative !! derivative at start on entry, at end on return
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer, parameter :: max_steps = 100000
        integer :: step_count
        real(dp) :: direction
        real(dp) :: error_norm
        real(dp) :: factor
        real(dp) :: h
        real(dp) :: remaining
        real(dp) :: scale1
        real(dp) :: scale2
        real(dp) :: x_current
        real(dp) :: y1
        real(dp) :: y1_fourth
        real(dp) :: y1_fifth
        real(dp) :: y2
        real(dp) :: y2_fourth
        real(dp) :: y2_fifth
        real(dp) :: k11, k12, k21, k22, k31, k32, k41, k42
        real(dp) :: k51, k52, k61, k62, k71, k72
        real(dp) :: t1, t2

        if (x_start == x_end) then
            status = spheroidal_status_success
            return
        end if

        direction = sign(1.0_dp, x_end - x_start)
        x_current = x_start
        y1 = value
        y2 = derivative
        h = direction * min(abs(x_end - x_start), &
            max(1.0e-4_dp, 0.1_dp / max(1.0_dp, c)))

        do step_count = 1, max_steps
            remaining = x_end - x_current
            if (direction * remaining <= 0.0_dp) exit
            if (abs(h) > abs(remaining)) h = remaining

            call radial_ode_rhs(m, c, cv, x_current, y1, y2, oblate, k11, k12)
            t1 = y1 + h * (1.0_dp / 5.0_dp) * k11
            t2 = y2 + h * (1.0_dp / 5.0_dp) * k12
            call radial_ode_rhs(m, c, cv, x_current + h / 5.0_dp, t1, t2, &
                oblate, k21, k22)

            t1 = y1 + h * (3.0_dp * k11 / 40.0_dp + 9.0_dp * k21 / 40.0_dp)
            t2 = y2 + h * (3.0_dp * k12 / 40.0_dp + 9.0_dp * k22 / 40.0_dp)
            call radial_ode_rhs(m, c, cv, x_current + 3.0_dp * h / 10.0_dp, &
                t1, t2, oblate, k31, k32)

            t1 = y1 + h * (44.0_dp * k11 / 45.0_dp - 56.0_dp * k21 / 15.0_dp + &
                32.0_dp * k31 / 9.0_dp)
            t2 = y2 + h * (44.0_dp * k12 / 45.0_dp - 56.0_dp * k22 / 15.0_dp + &
                32.0_dp * k32 / 9.0_dp)
            call radial_ode_rhs(m, c, cv, x_current + 4.0_dp * h / 5.0_dp, &
                t1, t2, oblate, k41, k42)

            t1 = y1 + h * (19372.0_dp * k11 / 6561.0_dp - &
                25360.0_dp * k21 / 2187.0_dp + 64448.0_dp * k31 / 6561.0_dp - &
                212.0_dp * k41 / 729.0_dp)
            t2 = y2 + h * (19372.0_dp * k12 / 6561.0_dp - &
                25360.0_dp * k22 / 2187.0_dp + 64448.0_dp * k32 / 6561.0_dp - &
                212.0_dp * k42 / 729.0_dp)
            call radial_ode_rhs(m, c, cv, x_current + 8.0_dp * h / 9.0_dp, &
                t1, t2, oblate, k51, k52)

            t1 = y1 + h * (9017.0_dp * k11 / 3168.0_dp - &
                355.0_dp * k21 / 33.0_dp + 46732.0_dp * k31 / 5247.0_dp + &
                49.0_dp * k41 / 176.0_dp - 5103.0_dp * k51 / 18656.0_dp)
            t2 = y2 + h * (9017.0_dp * k12 / 3168.0_dp - &
                355.0_dp * k22 / 33.0_dp + 46732.0_dp * k32 / 5247.0_dp + &
                49.0_dp * k42 / 176.0_dp - 5103.0_dp * k52 / 18656.0_dp)
            call radial_ode_rhs(m, c, cv, x_current + h, t1, t2, oblate, k61, k62)

            y1_fifth = y1 + h * (35.0_dp * k11 / 384.0_dp + &
                500.0_dp * k31 / 1113.0_dp + 125.0_dp * k41 / 192.0_dp - &
                2187.0_dp * k51 / 6784.0_dp + 11.0_dp * k61 / 84.0_dp)
            y2_fifth = y2 + h * (35.0_dp * k12 / 384.0_dp + &
                500.0_dp * k32 / 1113.0_dp + 125.0_dp * k42 / 192.0_dp - &
                2187.0_dp * k52 / 6784.0_dp + 11.0_dp * k62 / 84.0_dp)
            call radial_ode_rhs(m, c, cv, x_current + h, y1_fifth, y2_fifth, &
                oblate, k71, k72)

            y1_fourth = y1 + h * (5179.0_dp * k11 / 57600.0_dp + &
                7571.0_dp * k31 / 16695.0_dp + 393.0_dp * k41 / 640.0_dp - &
                92097.0_dp * k51 / 339200.0_dp + 187.0_dp * k61 / 2100.0_dp + &
                k71 / 40.0_dp)
            y2_fourth = y2 + h * (5179.0_dp * k12 / 57600.0_dp + &
                7571.0_dp * k32 / 16695.0_dp + 393.0_dp * k42 / 640.0_dp - &
                92097.0_dp * k52 / 339200.0_dp + 187.0_dp * k62 / 2100.0_dp + &
                k72 / 40.0_dp)

            if (.not. ieee_is_finite(y1_fifth) .or. &
                    .not. ieee_is_finite(y2_fifth) .or. &
                    .not. ieee_is_finite(y1_fourth) .or. &
                    .not. ieee_is_finite(y2_fourth)) then
                status = spheroidal_status_numerical_failure
                return
            end if

            scale1 = max(1.0_dp, abs(y1), abs(y1_fifth))
            scale2 = max(1.0_dp, abs(y2), abs(y2_fifth))
            error_norm = max(abs(y1_fifth - y1_fourth) / scale1, &
                abs(y2_fifth - y2_fourth) / scale2) / radial2_ode_tolerance

            if (error_norm <= 1.0_dp) then
                x_current = x_current + h
                y1 = y1_fifth
                y2 = y2_fifth
                if (x_current == x_end) exit
            end if

            if (error_norm == 0.0_dp) then
                factor = 5.0_dp
            else
                factor = 0.9_dp * error_norm ** (-0.2_dp)
                factor = min(5.0_dp, max(0.2_dp, factor))
            end if
            h = h * factor
            if (abs(h) < spacing(x_current)) then
                status = spheroidal_status_numerical_failure
                return
            end if
        end do

        if (direction * (x_end - x_current) > 8.0_dp * spacing(x_end)) then
            status = spheroidal_status_numerical_failure
            return
        end if
        value = y1
        derivative = y2
        status = spheroidal_status_success
    end subroutine integrate_radial2_ode

    pure subroutine radial_ode_rhs(m, c, cv, x, value, derivative, oblate, &
            value_rate, derivative_rate)
        integer, intent(in) :: m !! nonnegative spheroidal order
        real(dp), intent(in) :: c !! positive spheroidal parameter
        real(dp), intent(in) :: cv !! SciPy-convention characteristic value
        real(dp), intent(in) :: x !! radial argument
        real(dp), intent(in) :: value !! radial-function value
        real(dp), intent(in) :: derivative !! radial-function first derivative
        logical, intent(in) :: oblate !! true for the oblate radial equation
        real(dp), intent(out) :: value_rate !! derivative of value
        real(dp), intent(out) :: derivative_rate !! derivative of first derivative

        real(dp) :: denominator
        real(dp) :: lambda
        real(dp) :: potential

        value_rate = derivative
        if (oblate) then
            denominator = 1.0_dp + x * x
            lambda = cv + c * c
            potential = lambda / denominator - c * c - &
                real(m, dp) * real(m, dp) / (denominator * denominator)
        else
            denominator = (x - 1.0_dp) * (x + 1.0_dp)
            lambda = cv - c * c
            potential = lambda / denominator - c * c + &
                real(m, dp) * real(m, dp) / (denominator * denominator)
        end if
        derivative_rate = -2.0_dp * x * derivative / denominator + potential * value
    end subroutine radial_ode_rhs

    pure subroutine spheroidal_minus_order_coefficients(m, n, coefficients, &
            minus_coefficients, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: coefficients(:) !! relative coefficients a_n,k^m
        real(dp), allocatable, intent(out) :: minus_coefficients(:) !! scaled coefficients a_n,k^-m
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer :: degree
        integer :: j
        integer :: target
        real(dp), allocatable :: log_relation(:)
        real(dp) :: log_scale
        real(dp) :: log_weight
        real(dp) :: step

        allocate(minus_coefficients(size(coefficients)))
        allocate(log_relation(size(coefficients)))
        minus_coefficients = 0.0_dp
        log_relation = 0.0_dp
        target = (n - m) / 2 + 1
        if (target < 1 .or. target > size(coefficients)) then
            status = spheroidal_status_numerical_failure
            return
        end if

        do j = target + 1, size(coefficients)
            if (modulo(n - m, 2) == 1) then
                degree = m + 2 * (j - 1) - 1
            else
                degree = m + 2 * (j - 1) - 2
            end if
            step = log(real(degree + m + 1, dp)) + log(real(degree + m + 2, dp)) - &
                log(real(degree - m + 1, dp)) - log(real(degree - m + 2, dp))
            log_relation(j) = log_relation(j - 1) + step
        end do
        do j = target - 1, 1, -1
            if (modulo(n - m, 2) == 1) then
                degree = m + 2 * j - 1
            else
                degree = m + 2 * j - 2
            end if
            step = log(real(degree + m + 1, dp)) + log(real(degree + m + 2, dp)) - &
                log(real(degree - m + 1, dp)) - log(real(degree - m + 2, dp))
            log_relation(j) = log_relation(j + 1) - step
        end do

        log_scale = -huge(1.0_dp)
        do j = 1, size(coefficients)
            if (coefficients(j) /= 0.0_dp) then
                log_weight = log(abs(coefficients(j))) + log_relation(j)
                log_scale = max(log_scale, log_weight)
            end if
        end do
        if (.not. ieee_is_finite(log_scale)) then
            status = spheroidal_status_numerical_failure
            return
        end if

        do j = 1, size(coefficients)
            if (coefficients(j) /= 0.0_dp) then
                log_weight = log(abs(coefficients(j))) + log_relation(j) - log_scale
                minus_coefficients(j) = sign(exp(log_weight), coefficients(j))
            end if
        end do
        if (.not. all(ieee_is_finite(minus_coefficients))) then
            status = spheroidal_status_numerical_failure
            return
        end if
        status = spheroidal_status_success
    end subroutine spheroidal_minus_order_coefficients

    pure subroutine oblate_rad1_origin(m, odd_family, c, first_coefficient, denominator, &
            value, derivative, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        logical, intent(in) :: odd_family !! true when n-m is odd
        real(dp), intent(in) :: c !! positive oblate spheroidal parameter
        real(dp), intent(in) :: first_coefficient !! scaled a_n,k^-m coefficient at lowest degree
        real(dp), intent(in) :: denominator !! scaled A_n^-m normalization sum
        real(dp), intent(out) :: value !! oblate first-kind radial value at the origin
        real(dp), intent(out) :: derivative !! derivative with respect to x at the origin
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer :: k
        integer :: order
        real(dp) :: leading

        if (odd_family) then
            order = m + 1
        else
            order = m
        end if
        leading = 1.0_dp
        do k = 1, order
            leading = leading * c / real(2 * k + 1, dp)
            if (.not. ieee_is_finite(leading)) then
                value = quiet_nan(c)
                derivative = value
                status = spheroidal_status_numerical_failure
                return
            end if
        end do
        leading = first_coefficient * leading / denominator

        if (odd_family) then
            value = 0.0_dp
            derivative = leading
        else
            value = leading
            derivative = 0.0_dp
        end if
        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(derivative)) then
            status = spheroidal_status_numerical_failure
            return
        end if
        status = spheroidal_status_success
    end subroutine oblate_rad1_origin

    pure subroutine spheroidal_ang1(m, n, c, cv, x, gamma_sign, value, derivative)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! angular argument, abs(x) < 1
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        real(dp), intent(out) :: value !! first-kind angular spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x

        integer :: dimension
        integer :: next_dimension
        integer :: status
        integer :: target
        real(dp) :: current_derivative
        real(dp) :: current_value
        real(dp) :: previous_derivative
        real(dp) :: previous_value
        real(dp) :: scale_derivative
        real(dp) :: scale_value
        logical :: have_previous

        if (m < 0 .or. n < m .or. .not. ieee_is_finite(c) .or. &
                .not. ieee_is_finite(cv) .or. .not. ieee_is_finite(x) .or. &
                abs(x) >= 1.0_dp) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        if (abs(c) > sqrt(huge(1.0_dp))) then
            value = quiet_nan(c)
            derivative = value
            return
        end if

        target = (n - m) / 2 + 1
        if (target > max_matrix_size - 8) then
            value = quiet_nan(c)
            derivative = value
            return
        end if

        if (.not. (c < 0.0_dp .or. c > 0.0_dp) .or. c * c == 0.0_dp) then
            call phase_free_ferrers_value(m, n, x, value, derivative, status)
            if (status /= spheroidal_status_success) then
                value = quiet_nan(x)
                derivative = value
            end if
            return
        end if

        dimension = max(24, target + 16)
        if (abs(c) < real(max_matrix_size - 16, dp)) then
            dimension = max(dimension, int(abs(c)) + 16)
        else
            dimension = max_matrix_size
        end if
        dimension = min(dimension, max_matrix_size)
        have_previous = .false.
        previous_value = 0.0_dp
        previous_derivative = 0.0_dp

        do
            call spheroidal_ang1_dimension(m, n, c, cv, x, gamma_sign, dimension, &
                current_value, current_derivative, status)
            if (status /= spheroidal_status_success) then
                value = quiet_nan(x)
                derivative = value
                return
            end if

            value = current_value
            derivative = current_derivative
            if (have_previous) then
                scale_value = max(1.0_dp, abs(value))
                scale_derivative = max(1.0_dp, abs(derivative))
                if (abs(value - previous_value) <= &
                        angular_convergence_tolerance * scale_value .and. &
                        abs(derivative - previous_derivative) <= &
                        angular_convergence_tolerance * scale_derivative) return
            end if

            if (dimension == max_matrix_size) then
                value = quiet_nan(x)
                derivative = value
                return
            end if

            previous_value = value
            previous_derivative = derivative
            have_previous = .true.
            next_dimension = max(dimension + 16, (3 * dimension + 1) / 2)
            dimension = min(max_matrix_size, next_dimension)
        end do
    end subroutine spheroidal_ang1

    pure subroutine spheroidal_ang1_dimension(m, n, c, cv, x, gamma_sign, &
            dimension, value, derivative, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: x !! angular argument, abs(x) < 1
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        integer, intent(in) :: dimension !! finite Ferrers expansion dimension
        real(dp), intent(out) :: value !! first-kind angular spheroidal function
        real(dp), intent(out) :: derivative !! derivative with respect to x
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: derivatives(:)
        real(dp), allocatable :: derivatives_zero(:)
        real(dp), allocatable :: ferrers(:)
        real(dp), allocatable :: ferrers_zero(:)
        integer :: degree
        integer :: j
        integer :: max_degree
        integer :: target
        integer :: target_index
        logical :: odd_family
        real(dp) :: correction_derivative
        real(dp) :: correction_norm
        real(dp) :: correction_value
        real(dp) :: normalizer
        real(dp) :: raw_derivative
        real(dp) :: raw_norm
        real(dp) :: raw_value
        real(dp) :: sign_factor
        real(dp) :: target_norm
        real(dp) :: term

        odd_family = modulo(n - m, 2) == 1
        target = (n - m) / 2 + 1
        if (dimension < target + 1) then
            value = quiet_nan(x)
            derivative = value
            status = spheroidal_status_numerical_failure
            return
        end if

        call spheroidal_angular_coefficients(m, n, c, cv, gamma_sign, dimension, &
            coefficients, status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        if (odd_family) then
            max_degree = m + 2 * dimension - 1
        else
            max_degree = m + 2 * dimension - 2
        end if
        call ferrers_sequence(m, max_degree, x, ferrers, derivatives, status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        call ferrers_sequence(m, max_degree, 0.0_dp, ferrers_zero, derivatives_zero, &
            status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if

        raw_value = 0.0_dp
        raw_derivative = 0.0_dp
        raw_norm = 0.0_dp
        correction_value = 0.0_dp
        correction_derivative = 0.0_dp
        correction_norm = 0.0_dp
        do j = 1, dimension
            if (odd_family) then
                degree = m + 2 * j - 1
            else
                degree = m + 2 * j - 2
            end if
            if (modulo(j - target, 2) == 0) then
                sign_factor = 1.0_dp
            else
                sign_factor = -1.0_dp
            end if

            term = sign_factor * coefficients(j) * ferrers(degree - m + 1)
            call compensated_add(raw_value, correction_value, term)
            term = sign_factor * coefficients(j) * derivatives(degree - m + 1)
            call compensated_add(raw_derivative, correction_derivative, term)
            if (odd_family) then
                term = sign_factor * coefficients(j) * derivatives_zero(degree - m + 1)
            else
                term = sign_factor * coefficients(j) * ferrers_zero(degree - m + 1)
            end if
            call compensated_add(raw_norm, correction_norm, term)
        end do
        raw_value = raw_value + correction_value
        raw_derivative = raw_derivative + correction_derivative
        raw_norm = raw_norm + correction_norm

        target_index = n - m + 1
        if (odd_family) then
            target_norm = derivatives_zero(target_index)
        else
            target_norm = ferrers_zero(target_index)
        end if
        if (modulo(m, 2) == 1) target_norm = -target_norm

        if (raw_norm == 0.0_dp .or. .not. ieee_is_finite(raw_norm) .or. &
                .not. ieee_is_finite(target_norm)) then
            value = quiet_nan(x)
            derivative = value
            status = spheroidal_status_numerical_failure
            return
        end if

        normalizer = target_norm / raw_norm
        value = normalizer * raw_value
        derivative = normalizer * raw_derivative
        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(derivative)) then
            status = spheroidal_status_numerical_failure
            return
        end if
        status = spheroidal_status_success
    end subroutine spheroidal_ang1_dimension

    pure subroutine spheroidal_angular_coefficients(m, n, c, cv, gamma_sign, &
            dimension, coefficients, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real spheroidal parameter
        real(dp), intent(in) :: cv !! precomputed SciPy-convention characteristic value
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        integer, intent(in) :: dimension !! finite tridiagonal expansion dimension
        real(dp), allocatable, intent(out) :: coefficients(:) !! relative DLMF eigenvector
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        real(dp), allocatable :: ratios(:)
        integer :: j
        integer :: target
        logical :: odd_family
        real(dp) :: denominator
        real(dp) :: diagonal
        real(dp) :: gamma_squared
        real(dp) :: lambda
        real(dp) :: lower
        real(dp) :: next_ratio
        real(dp) :: scale
        real(dp) :: upper

        allocate(coefficients(dimension))
        allocate(ratios(dimension))
        coefficients = 0.0_dp
        ratios = 0.0_dp
        gamma_squared = gamma_sign * c * c
        lambda = cv - gamma_squared
        odd_family = modulo(n - m, 2) == 1
        target = (n - m) / 2 + 1

        if (target >= 2) then
            if (odd_family) then
                diagonal = odd_diagonal(m, 1, gamma_squared)
                upper = odd_upper(m, 1, gamma_squared)
            else
                diagonal = even_diagonal(m, 1, gamma_squared)
                upper = even_upper(m, 1, gamma_squared)
            end if
            if (upper == 0.0_dp) then
                status = spheroidal_status_numerical_failure
                return
            end if
            ratios(2) = -(diagonal - lambda) / upper

            do j = 3, target
                if (ratios(j - 1) == 0.0_dp) then
                    status = spheroidal_status_numerical_failure
                    return
                end if
                if (odd_family) then
                    lower = odd_lower_next(m, j - 2, gamma_squared)
                    diagonal = odd_diagonal(m, j - 1, gamma_squared)
                    upper = odd_upper(m, j - 1, gamma_squared)
                else
                    lower = even_lower_next(m, j - 2, gamma_squared)
                    diagonal = even_diagonal(m, j - 1, gamma_squared)
                    upper = even_upper(m, j - 1, gamma_squared)
                end if
                if (upper == 0.0_dp) then
                    status = spheroidal_status_numerical_failure
                    return
                end if
                ratios(j) = -(lower / ratios(j - 1) + diagonal - lambda) / upper
            end do
        end if

        next_ratio = 0.0_dp
        do j = dimension, target + 1, -1
            if (odd_family) then
                lower = odd_lower_next(m, j - 1, gamma_squared)
                diagonal = odd_diagonal(m, j, gamma_squared)
                if (j < dimension) then
                    upper = odd_upper(m, j, gamma_squared)
                else
                    upper = 0.0_dp
                end if
            else
                lower = even_lower_next(m, j - 1, gamma_squared)
                diagonal = even_diagonal(m, j, gamma_squared)
                if (j < dimension) then
                    upper = even_upper(m, j, gamma_squared)
                else
                    upper = 0.0_dp
                end if
            end if
            denominator = diagonal - lambda + upper * next_ratio
            if (denominator == 0.0_dp) then
                status = spheroidal_status_numerical_failure
                return
            end if
            ratios(j) = -lower / denominator
            next_ratio = ratios(j)
        end do

        coefficients(target) = 1.0_dp
        do j = target + 1, dimension
            coefficients(j) = coefficients(j - 1) * ratios(j)
        end do
        do j = target, 2, -1
            if (ratios(j) == 0.0_dp) then
                status = spheroidal_status_numerical_failure
                return
            end if
            coefficients(j - 1) = coefficients(j) / ratios(j)
        end do

        if (.not. all(ieee_is_finite(coefficients))) then
            status = spheroidal_status_numerical_failure
            return
        end if
        scale = maxval(abs(coefficients))
        if (scale == 0.0_dp .or. .not. ieee_is_finite(scale)) then
            status = spheroidal_status_numerical_failure
            return
        end if
        coefficients = coefficients / scale
        status = spheroidal_status_success
    end subroutine spheroidal_angular_coefficients

    pure subroutine phase_free_ferrers_value(m, n, x, value, derivative, status)
        integer, intent(in) :: m !! nonnegative associated-Legendre order
        integer, intent(in) :: n !! associated-Legendre degree, n >= m
        real(dp), intent(in) :: x !! Ferrers argument, abs(x) < 1
        real(dp), intent(out) :: value !! phase-free associated-Legendre value
        real(dp), intent(out) :: derivative !! derivative with respect to x
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        real(dp), allocatable :: derivatives(:)
        real(dp), allocatable :: ferrers(:)
        integer :: index

        call ferrers_sequence(m, n, x, ferrers, derivatives, status)
        if (status /= spheroidal_status_success) then
            value = quiet_nan(x)
            derivative = value
            return
        end if
        index = n - m + 1
        value = ferrers(index)
        derivative = derivatives(index)
        if (modulo(m, 2) == 1) then
            value = -value
            derivative = -derivative
        end if
    end subroutine phase_free_ferrers_value

    pure subroutine ferrers_sequence(m, max_degree, x, values, derivatives, status)
        integer, intent(in) :: m !! nonnegative associated-Legendre order
        integer, intent(in) :: max_degree !! largest degree requested, >= m
        real(dp), intent(in) :: x !! Ferrers argument, abs(x) < 1
        real(dp), allocatable, intent(out) :: values(:) !! P_l^m(x), l=m..max_degree
        real(dp), allocatable, intent(out) :: derivatives(:) !! derivatives for l=m..max_degree
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer :: count
        integer :: degree
        integer :: index
        integer :: k
        real(dp) :: current
        real(dp) :: denominator
        real(dp) :: previous
        real(dp) :: radial

        count = max_degree - m + 1
        if (m < 0 .or. count < 1 .or. abs(x) >= 1.0_dp .or. .not. ieee_is_finite(x)) then
            allocate(values(0))
            allocate(derivatives(0))
            status = spheroidal_status_invalid_argument
            return
        end if

        allocate(values(count))
        allocate(derivatives(count))
        values = 0.0_dp
        derivatives = 0.0_dp
        radial = sqrt(max(0.0_dp, 1.0_dp - x * x))
        current = 1.0_dp
        do k = 1, m
            current = -real(2 * k - 1, dp) * radial * current
            if (.not. ieee_is_finite(current)) then
                status = spheroidal_status_numerical_failure
                return
            end if
        end do
        values(1) = current

        if (count >= 2) then
            values(2) = real(2 * m + 1, dp) * x * values(1)
        end if
        do degree = m + 2, max_degree
            index = degree - m + 1
            values(index) = (real(2 * degree - 1, dp) * x * values(index - 1) - &
                real(degree + m - 1, dp) * values(index - 2)) / real(degree - m, dp)
        end do
        if (.not. all(ieee_is_finite(values))) then
            status = spheroidal_status_numerical_failure
            return
        end if

        denominator = x * x - 1.0_dp
        do degree = m, max_degree
            index = degree - m + 1
            if (degree == m) then
                previous = 0.0_dp
            else
                previous = values(index - 1)
            end if
            derivatives(index) = (real(degree, dp) * x * values(index) - &
                real(degree + m, dp) * previous) / denominator
        end do
        if (.not. all(ieee_is_finite(derivatives))) then
            status = spheroidal_status_numerical_failure
            return
        end if
        status = spheroidal_status_success
    end subroutine ferrers_sequence

    pure subroutine compensated_add(total, correction, term)
        real(dp), intent(inout) :: total !! running compensated sum
        real(dp), intent(inout) :: correction !! low-order correction to the running sum
        real(dp), intent(in) :: term !! next summand

        real(dp) :: updated

        updated = total + term
        if (abs(total) >= abs(term)) then
            correction = correction + (total - updated) + term
        else
            correction = correction + (term - updated) + total
        end if
        total = updated
    end subroutine compensated_add

    pure subroutine spheroidal_cv_sequence(m, n, c, gamma_sign, values, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! largest spheroidal degree requested
        real(dp), intent(in) :: c !! finite real spheroidal parameter
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        real(dp), allocatable, intent(out) :: values(:) !! characteristic values for degrees m through n
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input or numerical failure

        integer :: degree

        if (m < 0 .or. .not. ieee_is_finite(c)) then
            allocate(values(0))
            if (present(status)) status = spheroidal_status_invalid_argument
            return
        end if
        if (n < m) then
            allocate(values(0))
            if (present(status)) status = spheroidal_status_success
            return
        end if

        allocate(values(n - m + 1))
        do degree = m, n
            values(degree - m + 1) = spheroidal_cv(m, degree, c, gamma_sign)
            if (.not. ieee_is_finite(values(degree - m + 1))) then
                if (present(status)) status = spheroidal_status_numerical_failure
                return
            end if
        end do
        if (present(status)) status = spheroidal_status_success
    end subroutine spheroidal_cv_sequence

    pure function spheroidal_cv(m, n, c, gamma_sign) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: c !! finite real spheroidal parameter
        real(dp), intent(in) :: gamma_sign !! +1 prolate; -1 oblate
        real(dp) :: value

        integer :: dimension
        integer :: next_dimension
        integer :: status
        integer :: target
        real(dp) :: current
        real(dp) :: gamma_squared
        real(dp) :: previous
        real(dp) :: scale
        logical :: have_previous

        if (m < 0 .or. n < m .or. .not. ieee_is_finite(c)) then
            value = quiet_nan(c)
            return
        end if

        if (.not. (c < 0.0_dp .or. c > 0.0_dp)) then
            value = real(n, dp) * (real(n, dp) + 1.0_dp)
            return
        end if
        if (abs(c) > sqrt(sqrt(huge(1.0_dp)))) then
            value = quiet_nan(c)
            return
        end if

        target = (n - m) / 2 + 1
        if (target > max_matrix_size - 8) then
            value = quiet_nan(c)
            return
        end if

        gamma_squared = gamma_sign * c * c
        dimension = max(16, target + 12)
        if (abs(c) < real((max_matrix_size - 12) ** 2, dp)) then
            dimension = max(dimension, int(sqrt(abs(c))) + 12)
        else
            dimension = max_matrix_size
        end if
        dimension = min(dimension, max_matrix_size)
        previous = 0.0_dp
        have_previous = .false.

        do
            call spheroidal_eigenvalue(m, n, gamma_squared, dimension, current, status)
            if (status /= spheroidal_status_success) then
                value = quiet_nan(c)
                return
            end if

            value = current
            scale = max(1.0_dp, abs(value))
            if (have_previous) then
                if (abs(value - previous) <= convergence_tolerance * scale) return
            end if

            if (dimension == max_matrix_size) then
                value = quiet_nan(c)
                return
            end if

            previous = value
            have_previous = .true.
            next_dimension = max(dimension + 8, (3 * dimension + 1) / 2)
            dimension = min(max_matrix_size, next_dimension)
        end do
    end function spheroidal_cv

    pure subroutine spheroidal_eigenvalue(m, n, gamma_squared, dimension, value, status)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: n !! spheroidal degree, n >= m
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        integer, intent(in) :: dimension !! finite tridiagonal matrix dimension
        real(dp), intent(out) :: value !! SciPy-convention characteristic value
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: off_diagonal(:)
        integer :: j
        integer :: target
        logical :: odd_family
        real(dp) :: lower
        real(dp) :: product_root
        real(dp) :: upper

        allocate(diagonal(dimension))
        allocate(off_diagonal(max(0, dimension - 1)))
        odd_family = modulo(n - m, 2) == 1

        do j = 1, dimension
            if (odd_family) then
                diagonal(j) = odd_diagonal(m, j, gamma_squared) + gamma_squared
            else
                diagonal(j) = even_diagonal(m, j, gamma_squared) + gamma_squared
            end if

            if (j < dimension) then
                if (odd_family) then
                    upper = odd_upper(m, j, gamma_squared)
                    lower = odd_lower_next(m, j, gamma_squared)
                else
                    upper = even_upper(m, j, gamma_squared)
                    lower = even_lower_next(m, j, gamma_squared)
                end if
                if ((upper < 0.0_dp .and. lower > 0.0_dp) .or. &
                        (upper > 0.0_dp .and. lower < 0.0_dp)) then
                    value = quiet_nan(gamma_squared)
                    status = spheroidal_status_numerical_failure
                    return
                end if
                product_root = sqrt(abs(upper)) * sqrt(abs(lower))
                off_diagonal(j) = sign(product_root, upper)
            end if
        end do

        if (.not. all(ieee_is_finite(diagonal)) .or. &
                .not. all(ieee_is_finite(off_diagonal))) then
            value = quiet_nan(gamma_squared)
            status = spheroidal_status_numerical_failure
            return
        end if

        target = (n - m) / 2
        call symmetric_tridiagonal_eigenvalue(diagonal, off_diagonal, target, value, status)
    end subroutine spheroidal_eigenvalue

    pure subroutine symmetric_tridiagonal_eigenvalue(diagonal, off_diagonal, target, &
            value, status)
        real(dp), intent(in) :: diagonal(:) !! diagonal of a real symmetric tridiagonal matrix
        real(dp), intent(in) :: off_diagonal(:) !! first subdiagonal, size n-1
        integer, intent(in) :: target !! zero-based ordered eigenvalue index
        real(dp), intent(out) :: value !! requested eigenvalue
        integer, intent(out) :: status !! zero on success; nonzero on invalid/numerical failure

        integer :: count
        integer :: i
        integer :: iteration
        integer :: n
        real(dp) :: lower
        real(dp) :: margin
        real(dp) :: midpoint
        real(dp) :: off_left
        real(dp) :: off_right
        real(dp) :: pivmin
        real(dp) :: radius
        real(dp) :: scale
        real(dp) :: upper

        n = size(diagonal)
        if (n < 1 .or. size(off_diagonal) /= n - 1 .or. target < 0 .or. target >= n) then
            value = quiet_nan(0.0_dp)
            status = spheroidal_status_numerical_failure
            return
        end if

        if (n == 1) then
            lower = diagonal(1)
            upper = diagonal(1)
        else
            radius = abs(off_diagonal(1))
            lower = diagonal(1) - radius
            upper = diagonal(1) + radius
            do i = 2, n - 1
                off_left = abs(off_diagonal(i - 1))
                off_right = abs(off_diagonal(i))
                radius = off_left + off_right
                lower = min(lower, diagonal(i) - radius)
                upper = max(upper, diagonal(i) + radius)
            end do
            radius = abs(off_diagonal(n - 1))
            lower = min(lower, diagonal(n) - radius)
            upper = max(upper, diagonal(n) + radius)
        end if

        if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) then
            value = quiet_nan(0.0_dp)
            status = spheroidal_status_numerical_failure
            return
        end if

        scale = max(1.0_dp, abs(lower), abs(upper))
        pivmin = tridiagonal_pivmin(off_diagonal)
        margin = 8.0_dp * epsilon(1.0_dp) * scale + 8.0_dp * pivmin
        lower = lower - margin
        upper = upper + margin

        do iteration = 1, max_bisection_iterations
            midpoint = lower + 0.5_dp * (upper - lower)
            count = sturm_count(diagonal, off_diagonal, midpoint, pivmin)
            if (count <= target) then
                lower = midpoint
            else
                upper = midpoint
            end if

            if (upper - lower <= 8.0_dp * epsilon(1.0_dp) * &
                    max(1.0_dp, abs(lower), abs(upper))) then
                value = lower + 0.5_dp * (upper - lower)
                status = spheroidal_status_success
                return
            end if
        end do

        value = lower + 0.5_dp * (upper - lower)
        status = spheroidal_status_numerical_failure
    end subroutine symmetric_tridiagonal_eigenvalue

    pure integer function sturm_count(diagonal, off_diagonal, x, pivmin) result(count)
        real(dp), intent(in) :: diagonal(:) !! tridiagonal diagonal
        real(dp), intent(in) :: off_diagonal(:) !! tridiagonal first subdiagonal
        real(dp), intent(in) :: x !! trial eigenvalue
        real(dp), intent(in) :: pivmin !! minimum safe pivot magnitude

        integer :: i
        real(dp) :: pivot

        count = 0
        pivot = diagonal(1) - x
        if (abs(pivot) < pivmin) pivot = -pivmin
        if (pivot < 0.0_dp) count = count + 1

        do i = 2, size(diagonal)
            pivot = diagonal(i) - x - (off_diagonal(i - 1) / pivot) * off_diagonal(i - 1)
            if (abs(pivot) < pivmin) pivot = -pivmin
            if (pivot < 0.0_dp) count = count + 1
        end do
    end function sturm_count

    pure function tridiagonal_pivmin(off_diagonal) result(value)
        real(dp), intent(in) :: off_diagonal(:) !! tridiagonal first subdiagonal
        real(dp) :: value

        real(dp) :: largest
        real(dp) :: squared

        if (size(off_diagonal) == 0) then
            value = tiny(1.0_dp)
            return
        end if

        largest = maxval(abs(off_diagonal))
        if (largest > sqrt(huge(1.0_dp))) then
            value = sqrt(tiny(1.0_dp))
            return
        end if
        squared = largest * largest
        value = tiny(1.0_dp) * max(1.0_dp, squared)
    end function tridiagonal_pivmin

    pure function even_diagonal(m, j, gamma_squared) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: j !! one-based matrix row
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        real(dp) :: value

        real(dp) :: a
        real(dp) :: b
        real(dp) :: denominator
        real(dp) :: jj
        real(dp) :: mm

        mm = real(m, dp)
        jj = real(j, dp)
        a = mm + 2.0_dp * jj - 2.0_dp
        b = mm + 2.0_dp * jj - 1.0_dp
        denominator = (2.0_dp * mm + 4.0_dp * jj - 5.0_dp) * &
            (2.0_dp * mm + 4.0_dp * jj - 1.0_dp)
        value = a * b - 2.0_dp * gamma_squared * &
            (a * b - 1.0_dp + mm * mm) / denominator
    end function even_diagonal

    pure function even_upper(m, j, gamma_squared) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: j !! one-based matrix row before the superdiagonal
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        real(dp) :: value

        real(dp) :: jj
        real(dp) :: mm

        mm = real(m, dp)
        jj = real(j, dp)
        value = -gamma_squared * &
            (2.0_dp * mm + 2.0_dp * jj - 1.0_dp) * &
            (2.0_dp * mm + 2.0_dp * jj) / &
            ((2.0_dp * mm + 4.0_dp * jj - 1.0_dp) * &
            (2.0_dp * mm + 4.0_dp * jj + 1.0_dp))
    end function even_upper

    pure function even_lower_next(m, j, gamma_squared) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: j !! row whose superdiagonal partner is requested
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        real(dp) :: value

        real(dp) :: jj
        real(dp) :: mm

        mm = real(m, dp)
        jj = real(j, dp)
        value = -gamma_squared * (2.0_dp * jj - 1.0_dp) * (2.0_dp * jj) / &
            ((2.0_dp * mm + 4.0_dp * jj - 3.0_dp) * &
            (2.0_dp * mm + 4.0_dp * jj - 1.0_dp))
    end function even_lower_next

    pure function odd_diagonal(m, j, gamma_squared) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: j !! one-based matrix row
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        real(dp) :: value

        real(dp) :: a
        real(dp) :: b
        real(dp) :: denominator
        real(dp) :: jj
        real(dp) :: mm

        mm = real(m, dp)
        jj = real(j, dp)
        a = mm + 2.0_dp * jj - 1.0_dp
        b = mm + 2.0_dp * jj
        denominator = (2.0_dp * mm + 4.0_dp * jj - 3.0_dp) * &
            (2.0_dp * mm + 4.0_dp * jj + 1.0_dp)
        value = a * b - 2.0_dp * gamma_squared * &
            (a * b - 1.0_dp + mm * mm) / denominator
    end function odd_diagonal

    pure function odd_upper(m, j, gamma_squared) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: j !! one-based matrix row before the superdiagonal
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        real(dp) :: value

        real(dp) :: jj
        real(dp) :: mm

        mm = real(m, dp)
        jj = real(j, dp)
        value = -gamma_squared * &
            (2.0_dp * mm + 2.0_dp * jj) * &
            (2.0_dp * mm + 2.0_dp * jj + 1.0_dp) / &
            ((2.0_dp * mm + 4.0_dp * jj + 1.0_dp) * &
            (2.0_dp * mm + 4.0_dp * jj + 3.0_dp))
    end function odd_upper

    pure function odd_lower_next(m, j, gamma_squared) result(value)
        integer, intent(in) :: m !! nonnegative spheroidal order
        integer, intent(in) :: j !! row whose superdiagonal partner is requested
        real(dp), intent(in) :: gamma_squared !! signed DLMF spheroidal parameter squared
        real(dp) :: value

        real(dp) :: jj
        real(dp) :: mm

        mm = real(m, dp)
        jj = real(j, dp)
        value = -gamma_squared * (2.0_dp * jj) * (2.0_dp * jj + 1.0_dp) / &
            ((2.0_dp * mm + 4.0_dp * jj - 1.0_dp) * &
            (2.0_dp * mm + 4.0_dp * jj + 1.0_dp))
    end function odd_lower_next

end module scifort_spheroidal
