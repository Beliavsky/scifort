! SPDX-License-Identifier: BSD-2-Clause
! Copyright (c) 2012, Jaydeep P. Bardhan
! Copyright (c) 2012, Matthew G. Knepley
! Copyright (c) 2014, Janani Padmanabhan
! Copyright (c) 2026 SciFort contributors
!
! The Lame-polynomial class matrices, coefficient scaling, and first-kind
! evaluation below are adapted from SciPy 1.17.0 scipy/special/_ellip_harm.pxd.
! SciFort replaces LAPACK DSTE(V)R with its pure symmetric Jacobi eigensolver.
! The second-kind and normalization quadrature uses the public defining
! integral formulas with independent endpoint-removing substitutions and
! adaptive Gauss-Legendre integration.

module scifort_ellipsoidal_harmonics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, symmetric_eigen_jacobi
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter, public :: ellipsoidal_status_success = 0
    integer, parameter, public :: ellipsoidal_status_invalid_argument = 1
    integer, parameter, public :: ellipsoidal_status_numerical_failure = 2

    integer, parameter :: class_k = 1
    integer, parameter :: class_l = 2
    integer, parameter :: class_m = 3
    integer, parameter :: class_n = 4
    integer, parameter :: max_quad_depth = 16
    real(dp), parameter :: quad_rel_tol = 3.0e-13_dp
    real(dp), parameter :: quad_abs_tol = 2.0e-15_dp

    ! Positive abscissas and weights of the 16-point Gauss-Legendre rule.
    real(dp), parameter :: gl_x(8) = [ &
        0.0950125098376374401853193354250_dp, 0.281603550779258913230460501460_dp, &
        0.458016777657227386342419442984_dp, 0.617876244402643748446671764049_dp, &
        0.755404408355003033895101194847_dp, 0.865631202387831743880467897712_dp, &
        0.944575023073232576077988415535_dp, 0.989400934991649932596154173450_dp]
    real(dp), parameter :: gl_w(8) = [ &
        0.189450610455068496285396723208_dp, 0.182603415044923588866763667969_dp, &
        0.169156519395002538189312079030_dp, 0.149595988816576732081501730547_dp, &
        0.124628971255533872052476282192_dp, 0.0951585116824927848099251076022_dp, &
        0.0622535239386478928628438369944_dp, 0.0271524594117540948517805724560_dp]

    public :: ellip_harm
    public :: ellip_harm_2
    public :: ellip_normal

contains

    pure elemental function ellip_harm(h2, k2, n, p, s, signm, signn) result(value)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis, finite and > 0
        real(dp), intent(in) :: k2 !! squared second focal semiaxis, finite and > h2
        integer, intent(in) :: n !! nonnegative polynomial degree
        integer, intent(in) :: p !! harmonic order in [1, 2*n+1]
        real(dp), intent(in) :: s !! real ellipsoidal coordinate
        integer, intent(in), optional :: signm !! L/N square-root sign, +1 or -1 (default +1)
        integer, intent(in), optional :: signn !! M/N square-root sign, +1 or -1 (default +1)
        real(dp) :: value

        integer :: cls
        integer :: info
        integer :: sm
        integer :: sn
        real(dp), allocatable :: coefficients(:)

        sm = 1
        sn = 1
        if (present(signm)) sm = signm
        if (present(signn)) sn = signn

        if (ieee_is_nan(s)) then
            value = quiet_nan(s)
            return
        end if

        call lame_coefficients(h2, k2, n, p, sm, sn, coefficients, cls, info)
        if (info /= ellipsoidal_status_success) then
            value = quiet_nan(h2)
            return
        end if

        value = evaluate_lame(h2, k2, n, cls, s, coefficients, sm, sn)
    end function ellip_harm

    pure elemental function ellip_harm_2(h2, k2, n, p, s) result(value)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis, finite and > 0
        real(dp), intent(in) :: k2 !! squared second focal semiaxis, finite and > h2
        integer, intent(in) :: n !! nonnegative polynomial degree
        integer, intent(in) :: p !! harmonic order in [1, 2*n+1]
        real(dp), intent(in) :: s !! real external ellipsoidal coordinate, |s| >= sqrt(k2)
        real(dp) :: value

        integer :: cls
        integer :: info
        real(dp) :: integral
        real(dp) :: magnitude
        real(dp) :: evalue
        real(dp), allocatable :: coefficients(:)

        if (ieee_is_nan(s)) then
            value = quiet_nan(s)
            return
        end if
        if (.not. ieee_is_finite(s)) then
            value = quiet_nan(s)
            return
        end if

        call lame_coefficients(h2, k2, n, p, 1, 1, coefficients, cls, info)
        if (info /= ellipsoidal_status_success) then
            value = quiet_nan(h2)
            return
        end if

        magnitude = abs(s)
        if (magnitude < sqrt(k2)) then
            value = quiet_nan(s)
            return
        end if

        integral = second_kind_integral(h2, k2, n, cls, magnitude, coefficients)
        if (.not. ieee_is_finite(integral)) then
            value = quiet_nan(s)
            return
        end if
        evalue = evaluate_lame(h2, k2, n, cls, s, coefficients, 1, 1)
        value = real(2 * n + 1, dp) * evalue * sign(integral, s)
    end function ellip_harm_2

    pure elemental function ellip_normal(h2, k2, n, p) result(value)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis, finite and > 0
        real(dp), intent(in) :: k2 !! squared second focal semiaxis, finite and > h2
        integer, intent(in) :: n !! nonnegative polynomial degree
        integer, intent(in) :: p !! harmonic order in [1, 2*n+1]
        real(dp) :: value

        integer :: cls
        integer :: info
        real(dp) :: a0
        real(dp) :: a2
        real(dp) :: b0
        real(dp) :: b2
        real(dp), allocatable :: coefficients(:)

        call lame_coefficients(h2, k2, n, p, 1, 1, coefficients, cls, info)
        if (info /= ellipsoidal_status_success) then
            value = quiet_nan(h2)
            return
        end if

        call normal_moments(h2, k2, n, cls, coefficients, a0, a2, b0, b2)
        if (.not. ieee_is_finite(a0) .or. .not. ieee_is_finite(a2) .or. &
                .not. ieee_is_finite(b0) .or. .not. ieee_is_finite(b2)) then
            value = quiet_nan(h2)
            return
        end if
        value = 8.0_dp * (a2 * b0 - a0 * b2)
    end function ellip_normal

    pure subroutine lame_coefficients(h2, k2, n, p, signm, signn, coefficients, cls, status)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! nonnegative polynomial degree
        integer, intent(in) :: p !! harmonic order in [1, 2*n+1]
        integer, intent(in) :: signm !! L/N square-root sign, +1 or -1
        integer, intent(in) :: signn !! M/N square-root sign, +1 or -1
        real(dp), allocatable, intent(out) :: coefficients(:) !! Lame polynomial coefficients
        integer, intent(out) :: cls !! internal K/L/M/N family code
        integer, intent(out) :: status !! ellipsoidal status code

        integer :: i
        integer :: info
        integer :: j
        integer :: r
        integer :: size_class
        integer :: target
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: gamma
        real(dp) :: normalizer
        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: f(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: scale(:)

        cls = 0
        if (.not. valid_parameters(h2, k2, n, p) .or. &
                abs(signm) /= 1 .or. abs(signn) /= 1) then
            allocate(coefficients(0))
            status = ellipsoidal_status_invalid_argument
            return
        end if

        r = n / 2
        if (p <= r + 1) then
            cls = class_k
            target = p
            size_class = r + 1
        else if (p <= r + 1 + n - r) then
            cls = class_l
            target = p - (r + 1)
            size_class = n - r
        else if (p <= r + 1 + 2 * (n - r)) then
            cls = class_m
            target = p - (r + 1) - (n - r)
            size_class = n - r
        else
            cls = class_n
            target = p - (r + 1) - 2 * (n - r)
            size_class = r
        end if

        if (size_class < 1 .or. target < 1 .or. target > size_class) then
            allocate(coefficients(0))
            status = ellipsoidal_status_invalid_argument
            return
        end if

        allocate(diagonal(size_class), f(size_class), g(size_class), scale(size_class))
        alpha = h2
        beta = k2 - h2
        gamma = alpha - beta

        do i = 1, size_class
            j = i - 1
            select case (cls)
            case (class_k)
                g(i) = -real((2 * j + 2) * (2 * j + 1), dp) * beta
                if (modulo(n, 2) == 1) then
                    f(i) = -alpha * real((2 * (r - (j + 1)) + 2) * &
                        (2 * ((j + 1) + r) + 1), dp)
                    diagonal(i) = real((2 * r + 1) * (2 * r + 2) - 4 * j * j, dp) * alpha + &
                        real((2 * j + 1) * (2 * j + 1), dp) * beta
                else
                    f(i) = -alpha * real((2 * (r - (j + 1)) + 2) * &
                        (2 * (r + (j + 1)) - 1), dp)
                    diagonal(i) = real(2 * r * (2 * r + 1), dp) * alpha - &
                        real(4 * j * j, dp) * gamma
                end if
            case (class_l)
                g(i) = -real((2 * j + 2) * (2 * j + 3), dp) * beta
                if (modulo(n, 2) == 1) then
                    f(i) = -alpha * real((2 * (r - (j + 1)) + 2) * &
                        (2 * ((j + 1) + r) + 1), dp)
                    diagonal(i) = real((2 * r + 1) * (2 * r + 2), dp) * alpha - &
                        real((2 * j + 1) * (2 * j + 1), dp) * gamma
                else
                    f(i) = -alpha * real((2 * (r - (j + 1))) * &
                        (2 * (r + (j + 1)) + 1), dp)
                    diagonal(i) = real(2 * r * (2 * r + 1) - &
                        (2 * j + 1) * (2 * j + 1), dp) * alpha + &
                        real((2 * j + 2) * (2 * j + 2), dp) * beta
                end if
            case (class_m)
                g(i) = -real((2 * j + 2) * (2 * j + 1), dp) * beta
                if (modulo(n, 2) == 1) then
                    f(i) = -alpha * real((2 * (r - (j + 1)) + 2) * &
                        (2 * ((j + 1) + r) + 1), dp)
                    diagonal(i) = real((2 * r + 1) * (2 * r + 2) - &
                        (2 * j + 1) * (2 * j + 1), dp) * alpha + &
                        real(4 * j * j, dp) * beta
                else
                    f(i) = -alpha * real((2 * (r - (j + 1))) * &
                        (2 * (r + (j + 1)) + 1), dp)
                    diagonal(i) = real(2 * r * (2 * r + 1), dp) * alpha - &
                        real((2 * j + 1) * (2 * j + 1), dp) * gamma
                end if
            case (class_n)
                g(i) = -real((2 * j + 2) * (2 * j + 3), dp) * beta
                if (modulo(n, 2) == 1) then
                    f(i) = -alpha * real((2 * (r - (j + 1))) * &
                        (2 * ((j + 1) + r) + 3), dp)
                    diagonal(i) = real((2 * r + 1) * (2 * r + 2), dp) * alpha - &
                        real((2 * j + 2) * (2 * j + 2), dp) * gamma
                else
                    f(i) = -alpha * real((2 * (r - (j + 1))) * &
                        (2 * (r + (j + 1)) + 1), dp)
                    diagonal(i) = real(2 * r * (2 * r + 1) - &
                        (2 * j + 2) * (2 * j + 2), dp) * alpha + &
                        real((2 * j + 1) * (2 * j + 1), dp) * beta
                end if
            end select
        end do

        scale(1) = 1.0_dp
        do i = 2, size_class
            if (abs(f(i - 1)) <= tiny(1.0_dp) .or. g(i - 1) / f(i - 1) <= 0.0_dp) then
                allocate(coefficients(0))
                status = ellipsoidal_status_numerical_failure
                return
            end if
            scale(i) = sqrt(g(i - 1) / f(i - 1)) * scale(i - 1)
        end do

        allocate(matrix(size_class, size_class), eigenvalues(size_class), &
            eigenvectors(size_class, size_class))
        matrix = 0.0_dp
        do i = 1, size_class
            matrix(i, i) = diagonal(i)
        end do
        do i = 1, size_class - 1
            matrix(i + 1, i) = g(i) * scale(i) / scale(i + 1)
        end do

        call symmetric_eigen_jacobi(matrix, eigenvalues, eigenvectors, info)
        if (info /= linalg_status_success) then
            allocate(coefficients(0))
            status = ellipsoidal_status_numerical_failure
            return
        end if

        allocate(coefficients(size_class))
        coefficients = eigenvectors(:, target) / scale
        normalizer = coefficients(size_class) / ((-h2)**(size_class - 1))
        if (.not. ieee_is_finite(normalizer) .or. abs(normalizer) <= tiny(1.0_dp)) then
            deallocate(coefficients)
            allocate(coefficients(0))
            status = ellipsoidal_status_numerical_failure
            return
        end if
        coefficients = coefficients / normalizer
        status = ellipsoidal_status_success
    end subroutine lame_coefficients

    pure function evaluate_lame(h2, k2, n, cls, s, coefficients, signm, signn) result(value)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! polynomial degree
        integer, intent(in) :: cls !! internal K/L/M/N family code
        real(dp), intent(in) :: s !! real coordinate
        real(dp), intent(in) :: coefficients(:) !! normalized polynomial coefficients
        integer, intent(in) :: signm !! L/N square-root sign
        integer, intent(in) :: signn !! M/N square-root sign
        real(dp) :: value

        integer :: i
        integer :: parity
        real(dp) :: lambda
        real(dp) :: poly
        real(dp) :: psi
        real(dp) :: s2
        real(dp) :: inf
        real(dp) :: sign_factor

        if (ieee_is_nan(s)) then
            value = quiet_nan(s)
            return
        end if
        if (.not. ieee_is_finite(s)) then
            if (n == 0) then
                value = 1.0_dp
                return
            end if
            parity = modulo(n, 2)
            sign_factor = 1.0_dp
            select case (cls)
            case (class_k)
                if (parity == 1) sign_factor = sign_factor * sign(1.0_dp, s)
            case (class_l)
                sign_factor = real(signm, dp)
                if (parity == 0) sign_factor = sign_factor * sign(1.0_dp, s)
            case (class_m)
                sign_factor = real(signn, dp)
                if (parity == 0) sign_factor = sign_factor * sign(1.0_dp, s)
            case (class_n)
                sign_factor = real(signm * signn, dp)
                if (parity == 1) sign_factor = sign_factor * sign(1.0_dp, s)
            end select
            inf = positive_infinity(1.0_dp)
            value = sign(inf, sign_factor)
            return
        end if

        s2 = s * s
        lambda = 1.0_dp - s2 / h2
        poly = coefficients(size(coefficients))
        do i = size(coefficients) - 1, 1, -1
            poly = poly * lambda + coefficients(i)
        end do

        parity = modulo(n, 2)
        select case (cls)
        case (class_k)
            if (parity == 0) then
                psi = 1.0_dp
            else
                psi = s
            end if
        case (class_l)
            psi = real(signm, dp) * sqrt(abs(s2 - h2))
            if (parity == 0) psi = psi * s
        case (class_m)
            psi = real(signn, dp) * sqrt(abs(s2 - k2))
            if (parity == 0) psi = psi * s
        case (class_n)
            psi = real(signm * signn, dp) * sqrt(abs(s2 - h2)) * sqrt(abs(s2 - k2))
            if (parity == 1) psi = psi * s
        case default
            value = quiet_nan(s)
            return
        end select
        value = poly * psi
    end function evaluate_lame

    pure function second_kind_integral(h2, k2, n, cls, s, coefficients) result(value)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! polynomial degree
        integer, intent(in) :: cls !! internal K/L/M/N family code
        real(dp), intent(in) :: s !! positive external coordinate, >= sqrt(k2)
        real(dp), intent(in) :: coefficients(:) !! normalized polynomial coefficients
        real(dp) :: value

        value = adaptive_integral(1, 0.0_dp, 0.5_dp * scifort_pi, h2, k2, n, cls, &
            s, coefficients, 0)
    end function second_kind_integral

    pure subroutine normal_moments(h2, k2, n, cls, coefficients, a0, a2, b0, b2)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! polynomial degree
        integer, intent(in) :: cls !! internal K/L/M/N family code
        real(dp), intent(in) :: coefficients(:) !! normalized polynomial coefficients
        real(dp), intent(out) :: a0 !! outer-coordinate zeroth moment
        real(dp), intent(out) :: a2 !! outer-coordinate second moment
        real(dp), intent(out) :: b0 !! inner-coordinate zeroth moment
        real(dp), intent(out) :: b2 !! inner-coordinate second moment

        real(dp), parameter :: half_pi = 0.5_dp * scifort_pi

        a0 = adaptive_integral(2, 0.0_dp, half_pi, h2, k2, n, cls, 0.0_dp, coefficients, 0)
        a2 = adaptive_integral(3, 0.0_dp, half_pi, h2, k2, n, cls, 0.0_dp, coefficients, 0)
        b0 = adaptive_integral(4, 0.0_dp, half_pi, h2, k2, n, cls, 0.0_dp, coefficients, 0)
        b2 = adaptive_integral(5, 0.0_dp, half_pi, h2, k2, n, cls, 0.0_dp, coefficients, 0)
    end subroutine normal_moments

    pure recursive function adaptive_integral(mode, left, right, h2, k2, n, cls, &
            s_external, coefficients, depth) result(value)
        integer, intent(in) :: mode !! integrand selector, 1 through 5
        real(dp), intent(in) :: left !! lower integration angle
        real(dp), intent(in) :: right !! upper integration angle
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! polynomial degree
        integer, intent(in) :: cls !! internal K/L/M/N family code
        real(dp), intent(in) :: s_external !! positive external coordinate for mode 1
        real(dp), intent(in) :: coefficients(:) !! normalized polynomial coefficients
        integer, intent(in) :: depth !! current adaptive subdivision depth
        real(dp) :: value

        real(dp) :: left_value
        real(dp) :: middle
        real(dp) :: right_value
        real(dp) :: split_value
        real(dp) :: tolerance
        real(dp) :: whole_value

        whole_value = gauss16_interval(mode, left, right, h2, k2, n, cls, &
            s_external, coefficients)
        middle = 0.5_dp * (left + right)
        left_value = gauss16_interval(mode, left, middle, h2, k2, n, cls, &
            s_external, coefficients)
        right_value = gauss16_interval(mode, middle, right, h2, k2, n, cls, &
            s_external, coefficients)
        split_value = left_value + right_value

        if (.not. ieee_is_finite(whole_value) .or. .not. ieee_is_finite(split_value)) then
            value = quiet_nan(left)
            return
        end if

        tolerance = quad_abs_tol + quad_rel_tol * abs(split_value)
        if (depth >= max_quad_depth .or. abs(split_value - whole_value) <= tolerance) then
            value = split_value
        else
            value = adaptive_integral(mode, left, middle, h2, k2, n, cls, &
                s_external, coefficients, depth + 1) + &
                adaptive_integral(mode, middle, right, h2, k2, n, cls, &
                s_external, coefficients, depth + 1)
        end if
    end function adaptive_integral

    pure function gauss16_interval(mode, left, right, h2, k2, n, cls, &
            s_external, coefficients) result(value)
        integer, intent(in) :: mode !! integrand selector, 1 through 5
        real(dp), intent(in) :: left !! lower integration angle
        real(dp), intent(in) :: right !! upper integration angle
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! polynomial degree
        integer, intent(in) :: cls !! internal K/L/M/N family code
        real(dp), intent(in) :: s_external !! positive external coordinate for mode 1
        real(dp), intent(in) :: coefficients(:) !! normalized polynomial coefficients
        real(dp) :: value

        integer :: i
        real(dp) :: center
        real(dp) :: half_width
        real(dp) :: theta1
        real(dp) :: theta2

        center = 0.5_dp * (left + right)
        half_width = 0.5_dp * (right - left)
        value = 0.0_dp
        do i = 1, 8
            theta1 = center - half_width * gl_x(i)
            theta2 = center + half_width * gl_x(i)
            value = value + gl_w(i) * (ellipsoidal_integrand(mode, theta1, h2, k2, n, cls, &
                s_external, coefficients) + ellipsoidal_integrand(mode, theta2, h2, k2, n, cls, &
                s_external, coefficients))
        end do
        value = half_width * value
    end function gauss16_interval

    pure function ellipsoidal_integrand(mode, theta, h2, k2, n, cls, &
            s_external, coefficients) result(value)
        integer, intent(in) :: mode !! integrand selector, 1 through 5
        real(dp), intent(in) :: theta !! transformed integration angle
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! polynomial degree
        integer, intent(in) :: cls !! internal K/L/M/N family code
        real(dp), intent(in) :: s_external !! positive external coordinate for mode 1
        real(dp), intent(in) :: coefficients(:) !! normalized polynomial coefficients
        real(dp) :: value

        real(dp) :: c
        real(dp) :: delta
        real(dp) :: denominator
        real(dp) :: evalue
        real(dp) :: h
        real(dp) :: k
        real(dp) :: s
        real(dp) :: st
        real(dp) :: u

        h = sqrt(h2)
        k = sqrt(k2)
        st = sin(theta)
        c = cos(theta)

        select case (mode)
        case (1)
            u = st / s_external
            if (abs(u) <= tiny(1.0_dp)) then
                value = 0.0_dp
                return
            end if
            s = 1.0_dp / u
            evalue = evaluate_lame(h2, k2, n, cls, s, coefficients, 1, 1)
            denominator = sqrt(max(0.0_dp, 1.0_dp - k2 * u * u)) * &
                sqrt(max(0.0_dp, 1.0_dp - h2 * u * u))
            if (denominator <= tiny(1.0_dp) .or. abs(evalue) <= tiny(1.0_dp)) then
                value = quiet_nan(theta)
            else
                value = (c / s_external) / (evalue * evalue * denominator)
            end if
        case (2, 3)
            delta = k2 - h2
            s = sqrt(h2 + delta * st * st)
            evalue = evaluate_lame(h2, k2, n, cls, s, coefficients, 1, 1)
            if (mode == 2) then
                value = evalue * evalue / s
            else
                value = s * evalue * evalue
            end if
        case (4, 5)
            s = h * st
            evalue = evaluate_lame(h2, k2, n, cls, s, coefficients, 1, 1)
            denominator = sqrt(max(0.0_dp, k2 - h2 * st * st))
            if (mode == 4) then
                value = evalue * evalue / denominator
            else
                value = s * s * evalue * evalue / denominator
            end if
        case default
            value = quiet_nan(theta)
        end select
    end function ellipsoidal_integrand

    pure logical function valid_parameters(h2, k2, n, p) result(valid)
        real(dp), intent(in) :: h2 !! squared first focal semiaxis
        real(dp), intent(in) :: k2 !! squared second focal semiaxis
        integer, intent(in) :: n !! proposed degree
        integer, intent(in) :: p !! proposed order

        valid = ieee_is_finite(h2) .and. ieee_is_finite(k2) .and. h2 > 0.0_dp .and. &
            k2 > h2 .and. n >= 0
        if (valid) valid = p >= 1 .and. p <= 2 * n + 1
    end function valid_parameters

end module scifort_ellipsoidal_harmonics
