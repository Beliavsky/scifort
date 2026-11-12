! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Integer-order angular and modified/radial Mathieu functions. This is original
! SciFort code based on the Fourier recurrences and Bessel cross-product
! expansions in NIST DLMF chapter 28. No third-party implementation source is
! copied or adapted.

module scifort_mathieu
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_bessel_jy, only : jv, yv
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, symmetric_eigen_jacobi
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: mathieu_status_success = 0
    integer, parameter, public :: mathieu_status_invalid_argument = 1
    integer, parameter, public :: mathieu_status_numerical_failure = 2

    integer, parameter :: family_even = 1
    integer, parameter :: family_odd = 2
    integer, parameter :: max_matrix_size = 256
    real(dp), parameter :: sqrt_two = &
        1.414213562373095048801688724209698078570_dp

    public :: mathieu_a
    public :: mathieu_b
    public :: mathieu_cem
    public :: mathieu_modcem1
    public :: mathieu_modcem2
    public :: mathieu_even_coef
    public :: mathieu_odd_coef
    public :: mathieu_modsem1
    public :: mathieu_modsem2
    public :: mathieu_sem

contains

    pure elemental function mathieu_a(m, q) result(value)
        integer, intent(in) :: m !! nonnegative integer Mathieu order
        real(dp), intent(in) :: q !! finite real Mathieu parameter
        real(dp) :: value

        integer :: family
        integer :: status
        real(dp), allocatable :: coeff(:)

        if (m < 0 .or. .not. ieee_is_finite(q)) then
            value = quiet_nan(q)
            return
        end if

        family = family_even
        if (q < 0.0_dp .and. modulo(m, 2) == 1) family = family_odd
        call solve_positive_mode(m, abs(q), family, value, coeff, status)
        if (status /= mathieu_status_success) value = quiet_nan(q)
    end function mathieu_a

    pure elemental function mathieu_b(m, q) result(value)
        integer, intent(in) :: m !! positive integer Mathieu order
        real(dp), intent(in) :: q !! finite real Mathieu parameter
        real(dp) :: value

        integer :: family
        integer :: status
        real(dp), allocatable :: coeff(:)

        if (m < 1 .or. .not. ieee_is_finite(q)) then
            value = quiet_nan(q)
            return
        end if

        family = family_odd
        if (q < 0.0_dp .and. modulo(m, 2) == 1) family = family_even
        call solve_positive_mode(m, abs(q), family, value, coeff, status)
        if (status /= mathieu_status_success) value = quiet_nan(q)
    end function mathieu_b

    pure subroutine mathieu_even_coef(m, q, coeff, status)
        integer, intent(in) :: m !! nonnegative integer Mathieu order
        real(dp), intent(in) :: q !! nonnegative finite Mathieu parameter
        real(dp), allocatable, intent(out) :: coeff(:) !! normalized Fourier coefficients in ascending harmonic order
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input or eigensolver failure

        integer :: info
        real(dp) :: eigenvalue

        if (m < 0 .or. q < 0.0_dp .or. .not. ieee_is_finite(q)) then
            allocate(coeff(0))
            if (present(status)) status = mathieu_status_invalid_argument
            return
        end if

        call solve_positive_mode(m, q, family_even, eigenvalue, coeff, info)
        if (present(status)) status = info
    end subroutine mathieu_even_coef

    pure subroutine mathieu_odd_coef(m, q, coeff, status)
        integer, intent(in) :: m !! positive integer Mathieu order
        real(dp), intent(in) :: q !! nonnegative finite Mathieu parameter
        real(dp), allocatable, intent(out) :: coeff(:) !! normalized Fourier coefficients in ascending harmonic order
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input or eigensolver failure

        integer :: info
        real(dp) :: eigenvalue

        if (m < 1 .or. q < 0.0_dp .or. .not. ieee_is_finite(q)) then
            allocate(coeff(0))
            if (present(status)) status = mathieu_status_invalid_argument
            return
        end if

        call solve_positive_mode(m, q, family_odd, eigenvalue, coeff, info)
        if (present(status)) status = info
    end subroutine mathieu_odd_coef

    pure elemental subroutine mathieu_cem(m, q, x, value, derivative)
        integer, intent(in) :: m !! nonnegative integer Mathieu order
        real(dp), intent(in) :: q !! finite real Mathieu parameter
        real(dp), intent(in) :: x !! angular argument in degrees, matching SciPy
        real(dp), intent(out) :: value !! even Mathieu function ce_m(x,q)
        real(dp), intent(out) :: derivative !! derivative with respect to the radian argument

        integer :: k
        integer :: norder
        integer :: status
        real(dp), allocatable :: coeff(:)
        real(dp) :: angle
        real(dp) :: characteristic
        real(dp) :: frequency

        if (m < 0 .or. .not. ieee_is_finite(q) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        if (q >= 0.0_dp) then
            call solve_positive_mode(m, q, family_even, characteristic, coeff, status)
        else if (modulo(m, 2) == 0) then
            call solve_positive_mode(m, -q, family_even, characteristic, coeff, status)
            if (status == mathieu_status_success) then
                norder = m / 2
                do k = 1, size(coeff)
                    if (modulo(norder - (k - 1), 2) /= 0) coeff(k) = -coeff(k)
                end do
            end if
        else
            call solve_positive_mode(m, -q, family_odd, characteristic, coeff, status)
            if (status == mathieu_status_success) then
                norder = m / 2
                do k = 1, size(coeff)
                    if (modulo(norder - (k - 1), 2) /= 0) coeff(k) = -coeff(k)
                end do
            end if
        end if

        if (status /= mathieu_status_success) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        angle = x * scifort_pi / 180.0_dp
        value = 0.0_dp
        derivative = 0.0_dp
        do k = 1, size(coeff)
            if (modulo(m, 2) == 0) then
                frequency = 2.0_dp * real(k - 1, dp)
            else
                frequency = real(2 * k - 1, dp)
            end if
            value = value + coeff(k) * cos(frequency * angle)
            derivative = derivative - frequency * coeff(k) * sin(frequency * angle)
        end do

    end subroutine mathieu_cem

    pure elemental subroutine mathieu_sem(m, q, x, value, derivative)
        integer, intent(in) :: m !! positive integer Mathieu order
        real(dp), intent(in) :: q !! finite real Mathieu parameter
        real(dp), intent(in) :: x !! angular argument in degrees, matching SciPy
        real(dp), intent(out) :: value !! odd Mathieu function se_m(x,q)
        real(dp), intent(out) :: derivative !! derivative with respect to the radian argument

        integer :: k
        integer :: norder
        integer :: status
        real(dp), allocatable :: coeff(:)
        real(dp) :: angle
        real(dp) :: characteristic
        real(dp) :: frequency

        if (m < 1 .or. .not. ieee_is_finite(q) .or. .not. ieee_is_finite(x)) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        if (q >= 0.0_dp) then
            call solve_positive_mode(m, q, family_odd, characteristic, coeff, status)
        else if (modulo(m, 2) == 0) then
            call solve_positive_mode(m, -q, family_odd, characteristic, coeff, status)
            if (status == mathieu_status_success) then
                norder = m / 2 - 1
                do k = 1, size(coeff)
                    if (modulo(norder - (k - 1), 2) /= 0) coeff(k) = -coeff(k)
                end do
            end if
        else
            call solve_positive_mode(m, -q, family_even, characteristic, coeff, status)
            if (status == mathieu_status_success) then
                norder = m / 2
                do k = 1, size(coeff)
                    if (modulo(norder - (k - 1), 2) /= 0) coeff(k) = -coeff(k)
                end do
            end if
        end if

        if (status /= mathieu_status_success) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        angle = x * scifort_pi / 180.0_dp
        value = 0.0_dp
        derivative = 0.0_dp
        do k = 1, size(coeff)
            if (modulo(m, 2) == 1) then
                frequency = real(2 * k - 1, dp)
            else
                frequency = 2.0_dp * real(k, dp)
            end if
            value = value + coeff(k) * sin(frequency * angle)
            derivative = derivative + frequency * coeff(k) * cos(frequency * angle)
        end do

    end subroutine mathieu_sem

    pure elemental subroutine mathieu_modcem1(m, q, x, value, derivative)
        integer, intent(in) :: m !! nonnegative integer radial Mathieu order
        real(dp), intent(in) :: q !! positive finite Mathieu parameter
        real(dp), intent(in) :: x !! radial argument in radians, matching SciPy numerical behavior
        real(dp), intent(out) :: value !! even radial Mathieu function Mc1_m(x,q)
        real(dp), intent(out) :: derivative !! derivative with respect to the radian argument x

        call mathieu_radial(m, q, x, family_even, 1, value, derivative)
    end subroutine mathieu_modcem1

    pure elemental subroutine mathieu_modcem2(m, q, x, value, derivative)
        integer, intent(in) :: m !! nonnegative integer radial Mathieu order
        real(dp), intent(in) :: q !! positive finite Mathieu parameter
        real(dp), intent(in) :: x !! radial argument in radians, matching SciPy numerical behavior
        real(dp), intent(out) :: value !! even radial Mathieu function Mc2_m(x,q)
        real(dp), intent(out) :: derivative !! derivative with respect to the radian argument x

        call mathieu_radial(m, q, x, family_even, 2, value, derivative)
    end subroutine mathieu_modcem2

    pure elemental subroutine mathieu_modsem1(m, q, x, value, derivative)
        integer, intent(in) :: m !! positive integer radial Mathieu order
        real(dp), intent(in) :: q !! positive finite Mathieu parameter
        real(dp), intent(in) :: x !! radial argument in radians, matching SciPy numerical behavior
        real(dp), intent(out) :: value !! odd radial Mathieu function Ms1_m(x,q)
        real(dp), intent(out) :: derivative !! derivative with respect to the radian argument x

        call mathieu_radial(m, q, x, family_odd, 1, value, derivative)
    end subroutine mathieu_modsem1

    pure elemental subroutine mathieu_modsem2(m, q, x, value, derivative)
        integer, intent(in) :: m !! positive integer radial Mathieu order
        real(dp), intent(in) :: q !! positive finite Mathieu parameter
        real(dp), intent(in) :: x !! radial argument in radians, matching SciPy numerical behavior
        real(dp), intent(out) :: value !! odd radial Mathieu function Ms2_m(x,q)
        real(dp), intent(out) :: derivative !! derivative with respect to the radian argument x

        call mathieu_radial(m, q, x, family_odd, 2, value, derivative)
    end subroutine mathieu_modsem2

    pure subroutine mathieu_radial(m, q, x, family, cylinder_kind, value, derivative)
        integer, intent(in) :: m !! integer radial Mathieu order for the selected parity family
        real(dp), intent(in) :: q !! positive finite Mathieu parameter
        real(dp), intent(in) :: x !! finite radial argument in radians
        integer, intent(in) :: family !! family_even for Mc or family_odd for Ms
        integer, intent(in) :: cylinder_kind !! one for first-kind J expansion, two for second-kind Y
        real(dp), intent(out) :: value !! requested radial Mathieu value
        real(dp), intent(out) :: derivative !! derivative with respect to x

        real(dp) :: first_derivative
        real(dp) :: first_origin_derivative
        real(dp) :: first_origin_value
        real(dp) :: first_value
        real(dp) :: second_origin_derivative
        real(dp) :: second_origin_value
        real(dp) :: second_derivative
        real(dp) :: second_value
        real(dp) :: reflection_coefficient

        if (q <= 0.0_dp .or. .not. ieee_is_finite(q) .or. .not. ieee_is_finite(x) .or. &
                (family /= family_even .and. family /= family_odd) .or. &
                (cylinder_kind /= 1 .and. cylinder_kind /= 2) .or. m < 0 .or. &
                (family == family_odd .and. m < 1)) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        if (x >= 0.0_dp) then
            call mathieu_radial_nonnegative(m, q, x, family, cylinder_kind, value, derivative)
            return
        end if

        call mathieu_radial_nonnegative(m, q, -x, family, cylinder_kind, &
            second_value, second_derivative)
        if (cylinder_kind == 1) then
            if (family == family_even) then
                value = second_value
                derivative = -second_derivative
            else
                value = -second_value
                derivative = second_derivative
            end if
            return
        end if

        call mathieu_radial_nonnegative(m, q, -x, family, 1, first_value, first_derivative)
        call mathieu_radial_nonnegative(m, q, 0.0_dp, family, 1, &
            first_origin_value, first_origin_derivative)
        call mathieu_radial_nonnegative(m, q, 0.0_dp, family, 2, &
            second_origin_value, second_origin_derivative)

        if (family == family_even) then
            if (first_origin_value == 0.0_dp) then
                value = quiet_nan(q)
                derivative = value
                return
            end if
            reflection_coefficient = 2.0_dp * second_origin_value / first_origin_value
            value = -second_value + reflection_coefficient * first_value
            derivative = second_derivative - reflection_coefficient * first_derivative
        else
            if (first_origin_derivative == 0.0_dp) then
                value = quiet_nan(q)
                derivative = value
                return
            end if
            reflection_coefficient = -2.0_dp * second_origin_derivative / &
                first_origin_derivative
            value = second_value + reflection_coefficient * first_value
            derivative = -second_derivative - reflection_coefficient * first_derivative
        end if
    end subroutine mathieu_radial

    pure subroutine mathieu_radial_nonnegative(m, q, x, family, cylinder_kind, value, derivative)
        integer, intent(in) :: m !! integer radial Mathieu order for the selected parity family
        real(dp), intent(in) :: q !! positive finite Mathieu parameter
        real(dp), intent(in) :: x !! finite nonnegative radial argument in radians
        integer, intent(in) :: family !! family_even for Mc or family_odd for Ms
        integer, intent(in) :: cylinder_kind !! one for first-kind J expansion, two for second-kind Y
        real(dp), intent(out) :: value !! requested radial Mathieu value
        real(dp), intent(out) :: derivative !! derivative with respect to x

        integer :: anchor
        integer :: ell
        integer :: order_a1
        integer :: order_a2
        integer :: order_b1
        integer :: order_b2
        integer :: peak
        integer :: status
        real(dp), allocatable :: coeff(:)
        real(dp) :: a
        real(dp) :: b
        real(dp) :: characteristic
        real(dp) :: ratio
        real(dp) :: prefactor
        real(dp) :: sign_ell
        real(dp) :: term1
        real(dp) :: term2
        real(dp) :: dterm1
        real(dp) :: dterm2

        call solve_positive_mode(m, q, family, characteristic, coeff, status)
        if (status /= mathieu_status_success .or. size(coeff) == 0) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        peak = maxloc(abs(coeff), dim=1)
        anchor = peak
        do ell = 1, peak
            if (abs(coeff(ell)) >= 0.1_dp * abs(coeff(peak))) then
                anchor = ell
                exit
            end if
        end do
        if (coeff(anchor) == 0.0_dp) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        a = sqrt(q) * exp(-x)
        b = sqrt(q) * exp(x)
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = quiet_nan(q)
            derivative = value
            return
        end if

        value = 0.0_dp
        derivative = 0.0_dp
        do ell = 0, size(coeff) - 1
            ratio = coeff(ell + 1) / coeff(anchor)
            sign_ell = 1.0_dp
            if (modulo(ell, 2) /= 0) sign_ell = -1.0_dp

            if (family == family_even .and. modulo(m, 2) == 0) then
                order_a1 = ell - (anchor - 1)
                order_b1 = ell + (anchor - 1)
                order_a2 = order_b1
                order_b2 = order_a1
                prefactor = parity_sign(m / 2)
                if (anchor == 1) prefactor = 0.5_dp * prefactor
            else if (family == family_even) then
                order_a1 = ell - (anchor - 1)
                order_b1 = ell + (anchor - 1) + 1
                order_a2 = order_b1
                order_b2 = order_a1
                prefactor = parity_sign((m - 1) / 2)
            else if (modulo(m, 2) == 1) then
                order_a1 = ell - (anchor - 1)
                order_b1 = ell + (anchor - 1) + 1
                order_a2 = order_b1
                order_b2 = order_a1
                prefactor = parity_sign((m - 1) / 2)
            else
                order_a1 = ell - (anchor - 1)
                order_b1 = ell + (anchor - 1) + 2
                order_a2 = order_b1
                order_b2 = order_a1
                prefactor = parity_sign(m / 2 - 1)
            end if

            term1 = bessel_j_integer(order_a1, a) * &
                cylinder_integer(cylinder_kind, order_b1, b)
            term2 = bessel_j_integer(order_a2, a) * &
                cylinder_integer(cylinder_kind, order_b2, b)
            dterm1 = -a * bessel_j_integer_derivative(order_a1, a) * &
                cylinder_integer(cylinder_kind, order_b1, b) + &
                b * bessel_j_integer(order_a1, a) * &
                cylinder_integer_derivative(cylinder_kind, order_b1, b)
            dterm2 = -a * bessel_j_integer_derivative(order_a2, a) * &
                cylinder_integer(cylinder_kind, order_b2, b) + &
                b * bessel_j_integer(order_a2, a) * &
                cylinder_integer_derivative(cylinder_kind, order_b2, b)

            if (family == family_odd) then
                value = value + sign_ell * ratio * (term1 - term2)
                derivative = derivative + sign_ell * ratio * (dterm1 - dterm2)
            else
                value = value + sign_ell * ratio * (term1 + term2)
                derivative = derivative + sign_ell * ratio * (dterm1 + dterm2)
            end if
        end do

        value = prefactor * value
        derivative = prefactor * derivative
    end subroutine mathieu_radial_nonnegative

    pure elemental function bessel_j_integer(order, x) result(value)
        integer, intent(in) :: order !! integer order of J
        real(dp), intent(in) :: x !! nonnegative finite Bessel argument
        real(dp) :: value

        value = jv(real(order, dp), x)
    end function bessel_j_integer

    pure elemental function bessel_j_integer_derivative(order, x) result(value)
        integer, intent(in) :: order !! integer order of J derivative
        real(dp), intent(in) :: x !! positive finite Bessel argument
        real(dp) :: value

        value = 0.5_dp * (jv(real(order - 1, dp), x) - jv(real(order + 1, dp), x))
    end function bessel_j_integer_derivative

    pure elemental function cylinder_integer(kind, order, x) result(value)
        integer, intent(in) :: kind !! one for J, two for Y
        integer, intent(in) :: order !! integer cylinder-function order
        real(dp), intent(in) :: x !! positive finite cylinder-function argument
        real(dp) :: value

        if (kind == 1) then
            value = jv(real(order, dp), x)
        else
            value = yv(real(order, dp), x)
        end if
    end function cylinder_integer

    pure elemental function cylinder_integer_derivative(kind, order, x) result(value)
        integer, intent(in) :: kind !! one for J derivative, two for Y derivative
        integer, intent(in) :: order !! integer cylinder-function order
        real(dp), intent(in) :: x !! positive finite cylinder-function argument
        real(dp) :: value

        value = 0.5_dp * (cylinder_integer(kind, order - 1, x) - &
            cylinder_integer(kind, order + 1, x))
    end function cylinder_integer_derivative

    pure elemental function parity_sign(n) result(value)
        integer, intent(in) :: n !! integer exponent in (-1)^n
        real(dp) :: value

        value = 1.0_dp
        if (modulo(n, 2) /= 0) value = -1.0_dp
    end function parity_sign

    pure subroutine solve_positive_mode(m, q, family, eigenvalue, coeff, status)
        integer, intent(in) :: m !! nonnegative integer order; odd-family requires m >= 1
        real(dp), intent(in) :: q !! nonnegative finite Mathieu parameter
        integer, intent(in) :: family !! family_even for ce/a or family_odd for se/b
        real(dp), intent(out) :: eigenvalue !! characteristic value for the requested mode
        real(dp), allocatable, intent(out) :: coeff(:) !! normalized Fourier coefficients
        integer, intent(out) :: status !! mathieu status code

        integer :: eig_index
        integer :: i
        integer :: info
        integer :: matrix_size
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: matrix(:, :)

        eigenvalue = quiet_nan(q)
        if (m < 0 .or. q < 0.0_dp .or. .not. ieee_is_finite(q) .or. &
                (family /= family_even .and. family /= family_odd) .or. &
                (family == family_odd .and. m < 1)) then
            allocate(coeff(0))
            status = mathieu_status_invalid_argument
            return
        end if

        matrix_size = choose_matrix_size(m, q)
        if (matrix_size < 1 .or. matrix_size > max_matrix_size) then
            allocate(coeff(0))
            status = mathieu_status_numerical_failure
            return
        end if

        allocate(matrix(matrix_size, matrix_size), eigenvalues(matrix_size), &
            eigenvectors(matrix_size, matrix_size))
        matrix = 0.0_dp

        if (family == family_even) then
            eig_index = m / 2 + 1
            if (modulo(m, 2) == 0) then
                do i = 1, matrix_size
                    matrix(i, i) = real(2 * (i - 1), dp)**2
                end do
                if (matrix_size > 1) then
                    matrix(2, 1) = sqrt_two * q
                    do i = 3, matrix_size
                        matrix(i, i - 1) = q
                    end do
                end if
            else
                do i = 1, matrix_size
                    matrix(i, i) = real(2 * i - 1, dp)**2
                end do
                matrix(1, 1) = matrix(1, 1) + q
                do i = 2, matrix_size
                    matrix(i, i - 1) = q
                end do
            end if
        else
            if (modulo(m, 2) == 1) then
                eig_index = m / 2 + 1
                do i = 1, matrix_size
                    matrix(i, i) = real(2 * i - 1, dp)**2
                end do
                matrix(1, 1) = matrix(1, 1) - q
                do i = 2, matrix_size
                    matrix(i, i - 1) = q
                end do
            else
                eig_index = m / 2
                do i = 1, matrix_size
                    matrix(i, i) = real(2 * i, dp)**2
                end do
                do i = 2, matrix_size
                    matrix(i, i - 1) = q
                end do
            end if
        end if

        call symmetric_eigen_jacobi(matrix, eigenvalues, eigenvectors, info)
        if (info /= linalg_status_success .or. eig_index < 1 .or. eig_index > matrix_size) then
            allocate(coeff(0))
            status = mathieu_status_numerical_failure
            return
        end if

        eigenvalue = eigenvalues(eig_index)
        allocate(coeff(matrix_size))
        coeff = eigenvectors(:, eig_index)

        ! DLMF fixes the remaining sign ambiguity by continuity from q=0. For
        ! q >= 0 this is equivalent to taking the first Fourier coefficient
        ! positive for each of the four tridiagonal parity classes.
        if (coeff(1) < 0.0_dp) coeff = -coeff

        if (family == family_even .and. modulo(m, 2) == 0) then
            coeff(1) = coeff(1) / sqrt_two
        end if

        status = mathieu_status_success
    end subroutine solve_positive_mode

    pure integer function choose_matrix_size(m, q) result(n)
        integer, intent(in) :: m !! nonnegative integer Mathieu order
        real(dp), intent(in) :: q !! nonnegative finite Mathieu parameter

        integer :: target
        integer :: q_terms

        target = m / 2 + 1
        q_terms = ceiling(2.0_dp * sqrt(q)) + 24
        n = max(32, target + 24, q_terms)
    end function choose_matrix_size

end module scifort_mathieu
