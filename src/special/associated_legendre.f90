! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Associated Legendre and spherical-harmonic helpers.
!
! The real-degree Ferrers function lpmv is implemented independently from the
! hypergeometric representation and integer-order connection formulas in NIST
! DLMF 14.3 and 14.6. Integer-degree values use the classical three-term
! associated-Legendre recurrence to avoid gamma-function pole cancellation.
!
! sph_legendre_p and sph_harm_y use the fully normalized integer-degree
! recurrence corresponding to DLMF 14.30. The broader integer-degree
! assoc_legendre_p interface and derivative helpers use truncated Taylor jets,
! which differentiate the same recurrence without finite differences. The
! Legendre-Q sequences use DLMF 14.3.7 outside the unit disk and the DLMF 14.10
! degree/order/derivative recurrences elsewhere.

module scifort_associated_legendre
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_gamma_family, only : gammasgn
    use scifort_hypergeometric, only : hyp2f1
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, positive_infinity, quiet_nan
    use scifort_special_elementary, only : gammaln
    implicit none
    private

    public :: assoc_legendre_p
    public :: assoc_legendre_p_all
    public :: assoc_legendre_p_derivatives
    public :: associated_legendre_status_invalid_argument
    public :: associated_legendre_status_success
    public :: legendre_p
    public :: legendre_p_all
    public :: legendre_p_derivatives
    public :: lpmv
    public :: lqmn
    public :: lqn
    public :: sph_harm_y
    public :: sph_harm_y_all
    public :: sph_harm_y_derivatives
    public :: sph_legendre_p
    public :: sph_legendre_p_all
    public :: sph_legendre_p_derivatives

    interface assoc_legendre_p
        module procedure assoc_legendre_p_real
        module procedure assoc_legendre_p_complex
    end interface assoc_legendre_p

    interface assoc_legendre_p_derivatives
        module procedure assoc_legendre_p_derivatives_real
        module procedure assoc_legendre_p_derivatives_complex
    end interface assoc_legendre_p_derivatives

    interface assoc_legendre_p_all
        module procedure assoc_legendre_p_all_real
        module procedure assoc_legendre_p_all_complex
    end interface assoc_legendre_p_all

    interface legendre_p
        module procedure legendre_p_real
        module procedure legendre_p_complex
    end interface legendre_p

    interface legendre_p_derivatives
        module procedure legendre_p_derivatives_real
        module procedure legendre_p_derivatives_complex
    end interface legendre_p_derivatives

    interface legendre_p_all
        module procedure legendre_p_all_real
        module procedure legendre_p_all_complex
    end interface legendre_p_all

    interface lqmn
        module procedure lqmn_real
        module procedure lqmn_complex
    end interface lqmn

    interface lqn
        module procedure lqn_real
        module procedure lqn_complex
    end interface lqn

    integer, parameter :: associated_legendre_status_success = 0
    integer, parameter :: associated_legendre_status_invalid_argument = 1
    real(dp), parameter :: log_two = 0.693147180559945309417232121458176568_dp
    real(dp), parameter :: log_sqrt_pi = 0.572364942924700087071713675676529356_dp

contains

    pure elemental function lpmv(m, v, x) result(value)
        integer, intent(in) :: m !! integer order, including negative orders
        real(dp), intent(in) :: v !! real degree
        real(dp), intent(in) :: x !! real Ferrers argument, -1 <= x <= 1
        real(dp) :: value

        integer :: degree
        integer :: mm
        real(dp) :: hyper_value
        real(dp) :: log_abs
        real(dp) :: sign_value

        if (ieee_is_nan(v) .or. ieee_is_nan(x) .or. .not. ieee_is_finite(v) .or. &
                .not. ieee_is_finite(x) .or. abs(x) > 1.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (m == -huge(m) - 1) then
            value = quiet_nan(x)
            return
        end if

        mm = abs(m)

        if (is_representable_integer_degree(v)) then
            if (v >= 0.0_dp) then
                degree = int(v)
            else
                degree = int(-v - 1.0_dp)
            end if
            value = integer_degree_lpmv(degree, m, x)
            return
        end if

        if (x == -1.0_dp) then
            if (m == 0) then
                value = -positive_infinity(x)
            else
                value = positive_infinity(x)
            end if
            return
        end if

        if (m > 0) then
            ! Differentiate the Gauss representation of P_v(x) m times.
            ! This avoids the large gamma-ratio cancellation in the
            ! positive/negative-order connection formula near x = -1.
            hyper_value = hyp2f1(real(mm, dp) - v, real(mm, dp) + v + 1.0_dp, &
                real(mm + 1, dp), 0.5_dp * (1.0_dp - x))
            if (ieee_is_nan(hyper_value)) then
                value = quiet_nan(x)
                return
            end if
            if (hyper_value == 0.0_dp) then
                value = sign(0.0_dp, hyper_value)
                return
            end if
            if (x == 1.0_dp) then
                value = sign(0.0_dp, hyper_value)
                return
            end if

            log_abs = log(abs(hyper_value)) + 0.5_dp * real(mm, dp) * &
                (log1p_safe(-x) + log1p_safe(x)) - real(mm, dp) * log(2.0_dp) - &
                gammaln(real(mm + 1, dp))
            log_abs = log_abs + gammaln(-v + real(mm, dp)) - gammaln(-v) + &
                gammaln(v + real(mm, dp) + 1.0_dp) - gammaln(v + 1.0_dp)
            sign_value = sign(1.0_dp, hyper_value) * gammasgn(-v + real(mm, dp)) * &
                gammasgn(-v) * gammasgn(v + real(mm, dp) + 1.0_dp) * gammasgn(v + 1.0_dp)
            value = signed_exp(log_abs, sign_value)
            return
        end if

        hyper_value = hyp2f1(-v, v + 1.0_dp, real(mm + 1, dp), 0.5_dp * (1.0_dp - x))
        if (ieee_is_nan(hyper_value)) then
            value = quiet_nan(x)
            return
        end if
        if (hyper_value == 0.0_dp) then
            value = sign(0.0_dp, hyper_value)
            return
        end if

        log_abs = log(abs(hyper_value)) - gammaln(real(mm + 1, dp))
        if (mm > 0) then
            if (x == 1.0_dp) then
                log_abs = -positive_infinity(x)
            else
                log_abs = log_abs + 0.5_dp * real(mm, dp) * &
                    (log1p_safe(-x) - log1p_safe(x))
            end if
        end if
        sign_value = sign(1.0_dp, hyper_value)
        value = signed_exp(log_abs, sign_value)
    end function lpmv

    pure elemental function legendre_p_real(n, z) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: z !! real evaluation argument
        real(dp) :: value

        real(dp) :: values(0:0)

        call legendre_p_derivatives_real(n, z, 0, values)
        value = values(0)
    end function legendre_p_real

    pure elemental function legendre_p_complex(n, z) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: z !! complex evaluation argument
        complex(dp) :: value

        complex(dp) :: values(0:0)

        call legendre_p_derivatives_complex(n, z, 0, values)
        value = values(0)
    end function legendre_p_complex

    pure subroutine legendre_p_derivatives_real(n, z, diff_n, values)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: z !! real evaluation argument
        integer, intent(in) :: diff_n !! highest derivative order, >= 0
        real(dp), intent(out) :: values(0:diff_n) !! value and derivatives with respect to z

        integer :: derivative
        real(dp) :: derivative_factor
        real(dp) :: coefficients(0:diff_n)
        real(dp) :: zjet(0:diff_n)

        if (diff_n < 0) return
        if (n < 0) then
            values = quiet_nan(z)
            return
        end if

        if (ieee_is_nan(z) .or. .not. ieee_is_finite(z)) then
            call legendre_nonfinite_derivatives_real(n, z, values)
            return
        end if

        zjet = 0.0_dp
        zjet(0) = z
        if (diff_n >= 1) zjet(1) = 1.0_dp
        call legendre_jet_real(n, zjet, coefficients)

        derivative_factor = 1.0_dp
        do derivative = 0, diff_n
            if (derivative > 0) then
                derivative_factor = derivative_factor * real(derivative, dp)
            end if
            values(derivative) = derivative_factor * coefficients(derivative)
        end do
    end subroutine legendre_p_derivatives_real

    pure subroutine legendre_p_derivatives_complex(n, z, diff_n, values)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: z !! complex evaluation argument
        integer, intent(in) :: diff_n !! highest derivative order, >= 0
        complex(dp), intent(out) :: values(0:diff_n) !! value and derivatives with respect to z

        complex(dp) :: coefficients(0:diff_n)
        real(dp) :: derivative_factor
        integer :: derivative
        real(dp) :: nan_value
        complex(dp) :: zjet(0:diff_n)

        if (diff_n < 0) return
        if (n < 0) then
            values = complex_nan(z)
            return
        end if
        if (complex_is_nan(z) .or. .not. complex_is_finite(z)) then
            nan_value = quiet_nan(real(z, dp))
            values = cmplx(nan_value, nan_value, dp)
            if (diff_n >= n) then
                values(n) = cmplx(legendre_highest_derivative(n), 0.0_dp, dp)
            end if
            if (diff_n > n) values(n + 1:diff_n) = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if

        zjet = cmplx(0.0_dp, 0.0_dp, dp)
        zjet(0) = z
        if (diff_n >= 1) zjet(1) = cmplx(1.0_dp, 0.0_dp, dp)
        call legendre_jet_complex(n, zjet, coefficients)

        derivative_factor = 1.0_dp
        do derivative = 0, diff_n
            if (derivative > 0) then
                derivative_factor = derivative_factor * real(derivative, dp)
            end if
            values(derivative) = derivative_factor * coefficients(derivative)
        end do
    end subroutine legendre_p_derivatives_complex

    pure subroutine legendre_p_all_real(n, z, values, status, diff_n)
        integer, intent(in) :: n !! largest nonnegative degree requested
        real(dp), intent(in) :: z !! real evaluation argument
        real(dp), allocatable, intent(out) :: values(:, :) !! derivatives by degree, indexed 0:d,0:n
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid n or diff_n
        integer, intent(in), optional :: diff_n !! highest derivative order, >= 0 (default 0)

        integer :: degree
        integer :: derivative_count

        derivative_count = 0
        if (present(diff_n)) derivative_count = diff_n
        if (n < 0 .or. derivative_count < 0) then
            allocate(values(0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(values(0:derivative_count, 0:n))
        do degree = 0, n
            call legendre_p_derivatives_real(degree, z, derivative_count, values(:, degree))
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine legendre_p_all_real

    pure subroutine legendre_p_all_complex(n, z, values, status, diff_n)
        integer, intent(in) :: n !! largest nonnegative degree requested
        complex(dp), intent(in) :: z !! complex evaluation argument
        complex(dp), allocatable, intent(out) :: values(:, :) !! derivatives by degree, indexed 0:d,0:n
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid n or diff_n
        integer, intent(in), optional :: diff_n !! highest derivative order, >= 0 (default 0)

        integer :: degree
        integer :: derivative_count

        derivative_count = 0
        if (present(diff_n)) derivative_count = diff_n
        if (n < 0 .or. derivative_count < 0) then
            allocate(values(0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(values(0:derivative_count, 0:n))
        do degree = 0, n
            call legendre_p_derivatives_complex(degree, z, derivative_count, values(:, degree))
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine legendre_p_all_complex

    pure subroutine assoc_legendre_p_all_real(n, m, z, values, status, branch_cut, norm, diff_n)
        integer, intent(in) :: n !! largest nonnegative degree requested
        integer, intent(in) :: m !! largest nonnegative absolute order requested
        real(dp), intent(in) :: z !! real evaluation argument
        real(dp), allocatable, intent(out) :: values(:, :, :) !! derivatives, degree, order -m:m
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid shape options
        integer, intent(in), optional :: branch_cut !! branch convention, 2 (default) or 3
        logical, intent(in), optional :: norm !! use spherical normalization when true
        integer, intent(in), optional :: diff_n !! highest derivative order, >= 0 (default 0)

        integer :: bc
        integer :: degree
        integer :: derivative_count
        real(dp) :: factor
        integer :: order
        logical :: use_norm

        bc = 2
        if (present(branch_cut)) bc = branch_cut
        use_norm = .false.
        if (present(norm)) use_norm = norm
        derivative_count = 0
        if (present(diff_n)) derivative_count = diff_n

        if (n < 0 .or. m < 0 .or. derivative_count < 0 .or. &
                (bc /= 2 .and. bc /= 3)) then
            allocate(values(0, 0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(values(0:derivative_count, 0:n, -m:m))
        values = 0.0_dp
        do degree = 0, n
            do order = -min(m, degree), min(m, degree)
                if (order == 0) then
                    call legendre_p_derivatives_real(degree, z, derivative_count, &
                        values(:, degree, order))
                    if (use_norm) then
                        factor = sqrt(real(2 * degree + 1, dp) / 2.0_dp)
                        values(:, degree, order) = factor * values(:, degree, order)
                    end if
                else
                    call assoc_legendre_p_derivatives_real(degree, order, z, derivative_count, &
                        values(:, degree, order), bc, use_norm)
                end if
            end do
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine assoc_legendre_p_all_real

    pure subroutine assoc_legendre_p_all_complex(n, m, z, values, status, branch_cut, norm, diff_n)
        integer, intent(in) :: n !! largest nonnegative degree requested
        integer, intent(in) :: m !! largest nonnegative absolute order requested
        complex(dp), intent(in) :: z !! complex evaluation argument
        complex(dp), allocatable, intent(out) :: values(:, :, :) !! derivatives, degree, order -m:m
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid shape options
        integer, intent(in), optional :: branch_cut !! branch convention, 2 (default) or 3
        logical, intent(in), optional :: norm !! use spherical normalization when true
        integer, intent(in), optional :: diff_n !! highest derivative order, >= 0 (default 0)

        integer :: bc
        integer :: degree
        integer :: derivative_count
        real(dp) :: factor
        integer :: order
        logical :: use_norm

        bc = 2
        if (present(branch_cut)) bc = branch_cut
        use_norm = .false.
        if (present(norm)) use_norm = norm
        derivative_count = 0
        if (present(diff_n)) derivative_count = diff_n

        if (n < 0 .or. m < 0 .or. derivative_count < 0 .or. &
                (bc /= 2 .and. bc /= 3)) then
            allocate(values(0, 0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(values(0:derivative_count, 0:n, -m:m))
        values = cmplx(0.0_dp, 0.0_dp, dp)
        do degree = 0, n
            do order = -min(m, degree), min(m, degree)
                if (order == 0) then
                    call legendre_p_derivatives_complex(degree, z, derivative_count, &
                        values(:, degree, order))
                    if (use_norm) then
                        factor = sqrt(real(2 * degree + 1, dp) / 2.0_dp)
                        values(:, degree, order) = factor * values(:, degree, order)
                    end if
                else
                    call assoc_legendre_p_derivatives_complex(degree, order, z, &
                        derivative_count, values(:, degree, order), bc, use_norm)
                end if
            end do
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine assoc_legendre_p_all_complex

    pure elemental function assoc_legendre_p_real(n, m, z, branch_cut, norm) result(value)
        integer, intent(in) :: n !! nonnegative integer degree
        integer, intent(in) :: m !! integer order
        real(dp), intent(in) :: z !! real argument
        integer, intent(in), optional :: branch_cut !! branch convention, 2 (default) or 3
        logical, intent(in), optional :: norm !! use spherical normalization when true
        real(dp) :: value

        integer :: bc
        complex(dp) :: cvalue
        logical :: use_norm

        bc = 2
        if (present(branch_cut)) bc = branch_cut
        use_norm = .false.
        if (present(norm)) use_norm = norm

        if (ieee_is_nan(z) .or. .not. ieee_is_finite(z)) then
            value = quiet_nan(z)
            return
        end if
        if (n < 0) then
            value = 0.0_dp
            return
        end if
        if (m == -huge(m) - 1) then
            value = quiet_nan(z)
            return
        end if
        if (abs(m) > n) then
            value = 0.0_dp
            return
        end if
        if (bc /= 2 .and. bc /= 3) then
            value = quiet_nan(z)
            return
        end if
        if (m /= 0) then
            if (bc == 2 .and. abs(z) > 1.0_dp) then
                value = quiet_nan(z)
                return
            end if
            if (bc == 3 .and. abs(z) < 1.0_dp) then
                value = quiet_nan(z)
                return
            end if
        end if

        cvalue = assoc_legendre_p_complex(n, m, cmplx(z, 0.0_dp, dp), bc, use_norm)
        value = real(cvalue, dp)
    end function assoc_legendre_p_real

    pure elemental function assoc_legendre_p_complex(n, m, z, branch_cut, norm) result(value)
        integer, intent(in) :: n !! nonnegative integer degree
        integer, intent(in) :: m !! integer order
        complex(dp), intent(in) :: z !! complex argument
        integer, intent(in), optional :: branch_cut !! branch convention, 2 (default) or 3
        logical, intent(in), optional :: norm !! use spherical normalization when true
        complex(dp) :: value

        integer :: bc
        complex(dp) :: coefficients(0:0)
        complex(dp) :: zjet(0:0)
        logical :: use_norm

        bc = 2
        if (present(branch_cut)) bc = branch_cut
        use_norm = .false.
        if (present(norm)) use_norm = norm

        if (complex_is_nan(z) .or. .not. complex_is_finite(z)) then
            value = complex_nan(z)
            return
        end if
        if (n < 0) then
            value = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        if (m == -huge(m) - 1) then
            value = complex_nan(z)
            return
        end if
        if (abs(m) > n) then
            value = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        if (bc /= 2 .and. bc /= 3) then
            value = complex_nan(z)
            return
        end if

        zjet(0) = z
        call assoc_legendre_jet(n, m, zjet, bc, use_norm, coefficients)
        value = coefficients(0)
    end function assoc_legendre_p_complex

    pure subroutine assoc_legendre_p_derivatives_real(n, m, z, diff_n, values, branch_cut, norm)
        integer, intent(in) :: n !! nonnegative integer degree
        integer, intent(in) :: m !! integer order
        real(dp), intent(in) :: z !! real argument
        integer, intent(in) :: diff_n !! highest derivative order, >= 0
        real(dp), intent(out) :: values(0:diff_n) !! value and derivatives with respect to z
        integer, intent(in), optional :: branch_cut !! branch convention, 2 (default) or 3
        logical, intent(in), optional :: norm !! use spherical normalization when true

        integer :: bc
        integer :: k
        complex(dp) :: cvalues(0:diff_n)
        logical :: use_norm
        real(dp) :: nan_value

        if (diff_n < 0) return

        bc = 2
        if (present(branch_cut)) bc = branch_cut
        use_norm = .false.
        if (present(norm)) use_norm = norm

        nan_value = quiet_nan(z)
        if (ieee_is_nan(z) .or. .not. ieee_is_finite(z) .or. &
                m == -huge(m) - 1 .or. (bc /= 2 .and. bc /= 3)) then
            values = nan_value
            return
        end if
        if (n < 0 .or. abs(m) > n) then
            values = 0.0_dp
            return
        end if
        if (m /= 0) then
            if (bc == 2 .and. abs(z) > 1.0_dp) then
                values = nan_value
                return
            end if
            if (bc == 3 .and. abs(z) < 1.0_dp) then
                values = nan_value
                return
            end if
        end if

        call assoc_legendre_p_derivatives_complex(n, m, cmplx(z, 0.0_dp, dp), &
            diff_n, cvalues, bc, use_norm)
        do k = 0, diff_n
            values(k) = real(cvalues(k), dp)
        end do
    end subroutine assoc_legendre_p_derivatives_real

    pure subroutine assoc_legendre_p_derivatives_complex(n, m, z, diff_n, values, branch_cut, norm)
        integer, intent(in) :: n !! nonnegative integer degree
        integer, intent(in) :: m !! integer order
        complex(dp), intent(in) :: z !! complex argument
        integer, intent(in) :: diff_n !! highest derivative order, >= 0
        complex(dp), intent(out) :: values(0:diff_n) !! value and derivatives with respect to z
        integer, intent(in), optional :: branch_cut !! branch convention, 2 (default) or 3
        logical, intent(in), optional :: norm !! use spherical normalization when true

        integer :: bc
        complex(dp) :: coefficients(0:diff_n)
        real(dp) :: derivative_factor
        integer :: k
        logical :: use_norm
        complex(dp) :: zjet(0:diff_n)

        if (diff_n < 0) return

        bc = 2
        if (present(branch_cut)) bc = branch_cut
        use_norm = .false.
        if (present(norm)) use_norm = norm

        if (complex_is_nan(z) .or. .not. complex_is_finite(z) .or. &
                m == -huge(m) - 1 .or. (bc /= 2 .and. bc /= 3)) then
            values = complex_nan(z)
            return
        end if
        if (n < 0 .or. abs(m) > n) then
            values = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if

        zjet = cmplx(0.0_dp, 0.0_dp, dp)
        zjet(0) = z
        if (diff_n >= 1) zjet(1) = cmplx(1.0_dp, 0.0_dp, dp)
        call assoc_legendre_jet(n, m, zjet, bc, use_norm, coefficients)

        derivative_factor = 1.0_dp
        do k = 0, diff_n
            if (k > 0) derivative_factor = derivative_factor * real(k, dp)
            values(k) = derivative_factor * coefficients(k)
        end do
    end subroutine assoc_legendre_p_derivatives_complex

    pure subroutine sph_legendre_p_derivatives(n, m, theta, diff_n, values)
        integer, intent(in) :: n !! integer degree; negative values return zero
        integer, intent(in) :: m !! integer order
        real(dp), intent(in) :: theta !! polar argument in radians
        integer, intent(in) :: diff_n !! highest derivative order with respect to theta, >= 0
        real(dp), intent(out) :: values(0:diff_n) !! value and theta derivatives

        complex(dp) :: coefficients(0:diff_n)
        real(dp) :: derivative_factor
        integer :: k
        complex(dp) :: zjet(0:diff_n)

        if (diff_n < 0) return
        if (ieee_is_nan(theta) .or. .not. ieee_is_finite(theta)) then
            values = quiet_nan(theta)
            return
        end if
        if (m == -huge(m) - 1 .or. n < 0 .or. abs(m) > n) then
            values = 0.0_dp
            return
        end if

        call cosine_jet(theta, zjet)
        call assoc_legendre_jet(n, m, zjet, 2, .true., coefficients)

        derivative_factor = 1.0_dp / sqrt(2.0_dp * scifort_pi)
        do k = 0, diff_n
            if (k > 0) derivative_factor = derivative_factor * real(k, dp)
            values(k) = derivative_factor * real(coefficients(k), dp)
        end do
    end subroutine sph_legendre_p_derivatives

    pure subroutine sph_harm_y_derivatives(n, m, theta, phi, diff_n, values)
        integer, intent(in) :: n !! integer harmonic degree; negative values return zero
        integer, intent(in) :: m !! integer harmonic order
        real(dp), intent(in) :: theta !! polar coordinate in radians
        real(dp), intent(in) :: phi !! azimuthal coordinate in radians
        integer, intent(in) :: diff_n !! maximum total derivative order, >= 0
        complex(dp), intent(out) :: values(0:diff_n, 0:diff_n) !! mixed theta/phi derivatives

        complex(dp) :: phase
        complex(dp) :: phi_factor
        integer :: i
        integer :: j
        real(dp) :: nan_value
        real(dp) :: theta_values(0:diff_n)

        if (diff_n < 0) return
        values = cmplx(0.0_dp, 0.0_dp, dp)
        if (ieee_is_nan(theta) .or. ieee_is_nan(phi) .or. .not. ieee_is_finite(theta) .or. &
                .not. ieee_is_finite(phi)) then
            nan_value = quiet_nan(theta)
            values = cmplx(nan_value, nan_value, dp)
            return
        end if

        call sph_legendre_p_derivatives(n, m, theta, diff_n, theta_values)
        if (ieee_is_nan(theta_values(0))) then
            nan_value = quiet_nan(theta)
            values = cmplx(nan_value, nan_value, dp)
            return
        end if

        phase = cmplx(cos(real(m, dp) * phi), sin(real(m, dp) * phi), dp)
        do i = 0, diff_n
            phi_factor = cmplx(1.0_dp, 0.0_dp, dp)
            do j = 0, diff_n - i
                values(i, j) = theta_values(i) * phase * phi_factor
                phi_factor = phi_factor * cmplx(0.0_dp, real(m, dp), dp)
            end do
        end do
    end subroutine sph_harm_y_derivatives

    pure subroutine sph_legendre_p_all(n, m, theta, values, status, diff_n)
        integer, intent(in) :: n !! largest nonnegative spherical degree requested
        integer, intent(in) :: m !! largest nonnegative absolute order requested
        real(dp), intent(in) :: theta !! polar argument in radians
        real(dp), allocatable, intent(out) :: values(:, :, :) !! derivatives, degree, order -m:m
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid n, m, or diff_n
        integer, intent(in), optional :: diff_n !! highest theta-derivative order, >= 0 (default 0)

        integer :: degree
        integer :: derivative_count
        integer :: order

        derivative_count = 0
        if (present(diff_n)) derivative_count = diff_n
        if (n < 0 .or. m < 0 .or. derivative_count < 0) then
            allocate(values(0, 0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(values(0:derivative_count, 0:n, -m:m))
        values = 0.0_dp
        do degree = 0, n
            do order = -min(m, degree), min(m, degree)
                if (degree == 0 .and. order == 0 .and. &
                        (ieee_is_nan(theta) .or. .not. ieee_is_finite(theta))) then
                    values(:, degree, order) = 0.0_dp
                    values(0, degree, order) = 0.5_dp / sqrt(scifort_pi)
                else
                    call sph_legendre_p_derivatives(degree, order, theta, derivative_count, &
                        values(:, degree, order))
                end if
            end do
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine sph_legendre_p_all

    pure subroutine sph_harm_y_all(n, m, theta, phi, values, status, diff_n)
        integer, intent(in) :: n !! largest nonnegative spherical-harmonic degree requested
        integer, intent(in) :: m !! largest nonnegative absolute order requested
        real(dp), intent(in) :: theta !! polar coordinate in radians
        real(dp), intent(in) :: phi !! azimuthal coordinate in radians
        complex(dp), allocatable, intent(out) :: values(:, :, :, :) !! mixed derivatives by degree/order
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid n, m, or diff_n
        integer, intent(in), optional :: diff_n !! maximum total derivative order, >= 0 (default 0)

        integer :: degree
        integer :: derivative_count
        integer :: i
        integer :: order
        real(dp), allocatable :: theta_values(:)

        derivative_count = 0
        if (present(diff_n)) derivative_count = diff_n
        if (n < 0 .or. m < 0 .or. derivative_count < 0) then
            allocate(values(0, 0, 0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(values(0:derivative_count, 0:derivative_count, 0:n, -m:m))
        allocate(theta_values(0:derivative_count))
        values = cmplx(0.0_dp, 0.0_dp, dp)
        do degree = 0, n
            do order = -min(m, degree), min(m, degree)
                if (order == 0) then
                    call sph_legendre_p_derivatives(degree, order, theta, derivative_count, &
                        theta_values)
                    do i = 0, derivative_count
                        values(i, 0, degree, order) = cmplx(theta_values(i), 0.0_dp, dp)
                    end do
                else
                    call sph_harm_y_derivatives(degree, order, theta, phi, derivative_count, &
                        values(:, :, degree, order))
                end if
            end do
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine sph_harm_y_all

    pure subroutine lqn_real(n, z, qn, qnd, status)
        integer, intent(in) :: n !! highest nonnegative degree to return
        real(dp), intent(in) :: z !! real evaluation argument, excluding the singular points +/-1
        real(dp), allocatable, intent(out) :: qn(:) !! Q_k(z), k=0,...,n, with lower bound zero
        real(dp), allocatable, intent(out) :: qnd(:) !! derivatives dQ_k(z)/dz, k=0,...,n
        integer, intent(out), optional :: status !! zero on success; nonzero when n is negative

        real(dp), allocatable :: qtable(:, :)
        real(dp), allocatable :: dtable(:, :)
        integer :: local_status

        if (n < 0) then
            allocate(qn(0), qnd(0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        call lqmn_real(0, n, z, qtable, dtable, local_status)
        allocate(qn(0:n), qnd(0:n))
        qn = qtable(0, 0:n)
        qnd = dtable(0, 0:n)
        if (present(status)) status = local_status
    end subroutine lqn_real

    pure subroutine lqn_complex(n, z, qn, qnd, status)
        integer, intent(in) :: n !! highest nonnegative degree to return
        complex(dp), intent(in) :: z !! complex evaluation argument, excluding the singular points +/-1
        complex(dp), allocatable, intent(out) :: qn(:) !! Q_k(z), k=0,...,n, with lower bound zero
        complex(dp), allocatable, intent(out) :: qnd(:) !! derivatives dQ_k(z)/dz, k=0,...,n
        integer, intent(out), optional :: status !! zero on success; nonzero when n is negative

        complex(dp), allocatable :: qtable(:, :)
        complex(dp), allocatable :: dtable(:, :)
        integer :: local_status

        if (n < 0) then
            allocate(qn(0), qnd(0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        call lqmn_complex(0, n, z, qtable, dtable, local_status)
        allocate(qn(0:n), qnd(0:n))
        qn = qtable(0, 0:n)
        qnd = dtable(0, 0:n)
        if (present(status)) status = local_status
    end subroutine lqn_complex

    pure subroutine lqmn_real(m, n, z, qmn, qmnd, status)
        integer, intent(in) :: m !! highest nonnegative order to return
        integer, intent(in) :: n !! highest nonnegative degree to return
        real(dp), intent(in) :: z !! real evaluation argument, excluding the singular points +/-1
        real(dp), allocatable, intent(out) :: qmn(:, :) !! Q_degree^order(z), indexed 0:m,0:n
        real(dp), allocatable, intent(out) :: qmnd(:, :) !! derivatives with respect to z, indexed 0:m,0:n
        integer, intent(out), optional :: status !! zero on success; nonzero for negative m or n

        complex(dp), allocatable :: cq(:, :)
        complex(dp), allocatable :: cd(:, :)
        integer :: local_status

        if (m < 0 .or. n < 0) then
            allocate(qmn(0, 0), qmnd(0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        call lqmn_complex(m, n, cmplx(z, 0.0_dp, dp), cq, cd, local_status)
        allocate(qmn(0:m, 0:n), qmnd(0:m, 0:n))
        qmn = real(cq, dp)
        qmnd = real(cd, dp)
        if (present(status)) status = local_status
    end subroutine lqmn_real

    pure subroutine lqmn_complex(m, n, z, qmn, qmnd, status)
        integer, intent(in) :: m !! highest nonnegative order to return
        integer, intent(in) :: n !! highest nonnegative degree to return
        complex(dp), intent(in) :: z !! complex evaluation argument, excluding the singular points +/-1
        complex(dp), allocatable, intent(out) :: qmn(:, :) !! Q_degree^order(z), indexed 0:m,0:n
        complex(dp), allocatable, intent(out) :: qmnd(:, :) !! derivatives with respect to z, indexed 0:m,0:n
        integer, intent(out), optional :: status !! zero on success; nonzero for negative m or n

        complex(dp) :: denominator_inverse
        integer :: degree
        integer :: order
        real(dp) :: nan_value
        complex(dp), allocatable :: qwork(:, :)

        if (m < 0 .or. n < 0) then
            allocate(qmn(0, 0), qmnd(0, 0))
            if (present(status)) status = associated_legendre_status_invalid_argument
            return
        end if

        allocate(qmn(0:m, 0:n), qmnd(0:m, 0:n), qwork(0:m, 0:n + 1))
        if (complex_is_nan(z) .or. .not. complex_is_finite(z) .or. &
                z == cmplx(1.0_dp, 0.0_dp, dp) .or. z == cmplx(-1.0_dp, 0.0_dp, dp)) then
            nan_value = quiet_nan(real(z, dp))
            qmn = cmplx(nan_value, nan_value, dp)
            qmnd = cmplx(nan_value, nan_value, dp)
            if (present(status)) status = associated_legendre_status_success
            return
        end if

        call legendre_q_table(m, n + 1, z, qwork, denominator_inverse)
        qmn = qwork(:, 0:n)
        do order = 0, m
            do degree = 0, n
                qmnd(order, degree) = (real(order - degree - 1, dp) * &
                    qwork(order, degree + 1) + real(degree + 1, dp) * z * &
                    qwork(order, degree)) * denominator_inverse
            end do
        end do
        if (present(status)) status = associated_legendre_status_success
    end subroutine lqmn_complex

    pure subroutine legendre_q_table(m, n, z, q, denominator_inverse)
        integer, intent(in) :: m !! highest nonnegative order in q
        integer, intent(in) :: n !! highest nonnegative degree in q
        complex(dp), intent(in) :: z !! finite nonsingular complex argument
        complex(dp), intent(out) :: q(0:m, 0:n) !! associated Legendre-Q table
        complex(dp), intent(out) :: denominator_inverse !! reciprocal of 1-z**2

        complex(dp) :: derivative
        integer :: degree
        complex(dp) :: inverse_z
        logical :: outside
        integer :: order
        complex(dp) :: radial

        outside = abs(z) >= 1.0_dp
        if (outside) then
            inverse_z = cmplx(1.0_dp, 0.0_dp, dp) / z
            denominator_inverse = -inverse_z * inverse_z / &
                (cmplx(1.0_dp, 0.0_dp, dp) - inverse_z * inverse_z)
            if (abs(z) > 1.0_dp) then
                do degree = 0, n
                    q(0, degree) = legendre_q_zero_outside(degree, z)
                end do
            else
                q(0, 0) = 0.5_dp * log((cmplx(1.0_dp, 0.0_dp, dp) + inverse_z) / &
                    (cmplx(1.0_dp, 0.0_dp, dp) - inverse_z))
                if (n >= 1) q(0, 1) = z * q(0, 0) - cmplx(1.0_dp, 0.0_dp, dp)
                do degree = 1, n - 1
                    q(0, degree + 1) = (real(2 * degree + 1, dp) * z * q(0, degree) - &
                        real(degree, dp) * q(0, degree - 1)) / real(degree + 1, dp)
                end do
            end if
            radial = outside_radial(z)
        else
            denominator_inverse = cmplx(1.0_dp, 0.0_dp, dp) / &
                (cmplx(1.0_dp, 0.0_dp, dp) - z * z)
            q(0, 0) = 0.5_dp * log((cmplx(1.0_dp, 0.0_dp, dp) + z) / &
                (cmplx(1.0_dp, 0.0_dp, dp) - z))
            if (n >= 1) q(0, 1) = z * q(0, 0) - cmplx(1.0_dp, 0.0_dp, dp)
            do degree = 1, n - 1
                q(0, degree + 1) = (real(2 * degree + 1, dp) * z * q(0, degree) - &
                    real(degree, dp) * q(0, degree - 1)) / real(degree + 1, dp)
            end do
            radial = sqrt(cmplx(1.0_dp, 0.0_dp, dp) - z * z)
        end if

        if (m == 0) return

        derivative = denominator_inverse
        if (outside) then
            q(1, 0) = radial * derivative
        else
            q(1, 0) = -radial * derivative
        end if
        do degree = 1, n
            derivative = real(degree, dp) * (q(0, degree - 1) - z * q(0, degree)) * &
                denominator_inverse
            if (outside) then
                q(1, degree) = radial * derivative
            else
                q(1, degree) = -radial * derivative
            end if
        end do

        do order = 0, m - 2
            do degree = 0, n
                if (outside) then
                    q(order + 2, degree) = -2.0_dp * real(order + 1, dp) * z / radial * &
                        q(order + 1, degree) + real((degree - order) * &
                        (degree + order + 1), dp) * q(order, degree)
                else
                    q(order + 2, degree) = -2.0_dp * real(order + 1, dp) * z / radial * &
                        q(order + 1, degree) - real((degree - order) * &
                        (degree + order + 1), dp) * q(order, degree)
                end if
            end do
        end do
    end subroutine legendre_q_table

    pure function legendre_q_zero_outside(n, z) result(value)
        integer, intent(in) :: n !! nonnegative integer degree
        complex(dp), intent(in) :: z !! argument with magnitude greater than one
        complex(dp) :: value

        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        complex(dp) :: hyper_value
        complex(dp) :: inverse_z
        real(dp) :: log_coefficient
        real(dp) :: real_hyper
        real(dp) :: sign_factor

        a = 0.5_dp * real(n, dp) + 1.0_dp
        b = 0.5_dp * real(n, dp) + 0.5_dp
        c = real(n, dp) + 1.5_dp
        inverse_z = cmplx(1.0_dp, 0.0_dp, dp) / z
        log_coefficient = log_sqrt_pi + gammaln(real(n + 1, dp)) - &
            real(n + 1, dp) * log_two - gammaln(real(n, dp) + 1.5_dp)

        if (aimag(z) == 0.0_dp) then
            real_hyper = hyp2f1(a, b, c, real(inverse_z * inverse_z, dp))
            sign_factor = 1.0_dp
            if (real(z, dp) < 0.0_dp .and. mod(n + 1, 2) == 1) sign_factor = -1.0_dp
            value = cmplx(sign_factor * exp(log_coefficient - real(n + 1, dp) * &
                log(abs(real(z, dp)))) * real_hyper, 0.0_dp, dp)
        else
            hyper_value = hyp2f1(a, b, c, inverse_z * inverse_z)
            value = exp(cmplx(log_coefficient, 0.0_dp, dp) - &
                real(n + 1, dp) * log(z)) * hyper_value
        end if
    end function legendre_q_zero_outside

    pure function outside_radial(z) result(value)
        complex(dp), intent(in) :: z !! complex argument with magnitude at least one and z /= +/-1
        complex(dp) :: value

        complex(dp) :: inverse_z
        complex(dp) :: root_z_squared

        inverse_z = cmplx(1.0_dp, 0.0_dp, dp) / z
        if (real(z, dp) > 0.0_dp .or. &
                (real(z, dp) == 0.0_dp .and. aimag(z) >= 0.0_dp)) then
            root_z_squared = z
        else
            root_z_squared = -z
        end if
        value = root_z_squared * sqrt(cmplx(1.0_dp, 0.0_dp, dp) - inverse_z * inverse_z)
    end function outside_radial

    pure elemental function sph_legendre_p(n, m, theta) result(value)
        integer, intent(in) :: n !! integer degree; negative values return zero
        integer, intent(in) :: m !! integer order
        real(dp), intent(in) :: theta !! real polar argument in radians
        real(dp) :: value

        integer :: k
        integer :: mm
        real(dp) :: a
        real(dp) :: b
        real(dp) :: current
        real(dp) :: next_value
        real(dp) :: previous
        real(dp) :: radial
        real(dp) :: x

        if (ieee_is_nan(theta) .or. .not. ieee_is_finite(theta)) then
            value = quiet_nan(theta)
            return
        end if
        if (m == -huge(m) - 1) then
            value = 0.0_dp
            return
        end if
        if (n < 0) then
            value = 0.0_dp
            return
        end if

        mm = abs(m)
        if (mm > n) then
            value = 0.0_dp
            return
        end if

        x = cos(theta)
        radial = sqrt(max(0.0_dp, 1.0_dp - x * x))

        current = 0.5_dp / sqrt(scifort_pi)
        if (mm > 0) then
            do k = 1, mm
                current = -sqrt((2.0_dp * real(k, dp) + 1.0_dp) / &
                    (2.0_dp * real(k, dp))) * radial * current
            end do
        end if

        if (n == mm) then
            value = current
        else
            previous = current
            current = sqrt(2.0_dp * real(mm, dp) + 3.0_dp) * x * previous
            if (n == mm + 1) then
                value = current
            else
                do k = mm + 2, n
                    a = sqrt((4.0_dp * real(k, dp) * real(k, dp) - 1.0_dp) / &
                        (real(k, dp) * real(k, dp) - real(mm, dp) * real(mm, dp)))
                    b = sqrt(((real(k - 1, dp) * real(k - 1, dp)) - &
                        real(mm, dp) * real(mm, dp)) / &
                        (4.0_dp * real(k - 1, dp) * real(k - 1, dp) - 1.0_dp))
                    next_value = a * (x * current - b * previous)
                    previous = current
                    current = next_value
                end do
                value = current
            end if
        end if

        if (m < 0 .and. mod(mm, 2) == 1) value = -value
    end function sph_legendre_p

    pure elemental function sph_harm_y(n, m, theta, phi) result(value)
        integer, intent(in) :: n !! integer harmonic degree; negative values return zero
        integer, intent(in) :: m !! integer harmonic order
        real(dp), intent(in) :: theta !! polar coordinate in radians
        real(dp), intent(in) :: phi !! azimuthal coordinate in radians
        complex(dp) :: value

        real(dp) :: amplitude
        real(dp) :: angle
        real(dp) :: nan_value

        if (ieee_is_nan(theta) .or. ieee_is_nan(phi) .or. .not. ieee_is_finite(theta) .or. &
                .not. ieee_is_finite(phi)) then
            nan_value = quiet_nan(theta)
            value = cmplx(nan_value, nan_value, dp)
            return
        end if

        amplitude = sph_legendre_p(n, m, theta)
        if (ieee_is_nan(amplitude)) then
            nan_value = quiet_nan(theta)
            value = cmplx(nan_value, nan_value, dp)
            return
        end if
        if (amplitude == 0.0_dp) then
            value = cmplx(amplitude, 0.0_dp, dp)
            return
        end if

        angle = real(m, dp) * phi
        value = amplitude * cmplx(cos(angle), sin(angle), dp)
    end function sph_harm_y

    pure subroutine legendre_nonfinite_derivatives_real(n, z, values)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: z !! NaN or infinite real argument
        real(dp), intent(out) :: values(0:) !! value and derivatives with respect to z

        integer :: derivative
        real(dp) :: highest_derivative
        real(dp) :: infinity_value
        real(dp) :: nan_value

        highest_derivative = legendre_highest_derivative(n)
        nan_value = quiet_nan(z)
        infinity_value = positive_infinity(z)

        do derivative = 0, ubound(values, 1)
            if (derivative > n) then
                values(derivative) = 0.0_dp
            else if (derivative == n) then
                values(derivative) = highest_derivative
            else if (ieee_is_nan(z)) then
                values(derivative) = nan_value
            else if (z > 0.0_dp .or. mod(n - derivative, 2) == 0) then
                values(derivative) = infinity_value
            else
                values(derivative) = -infinity_value
            end if
        end do
    end subroutine legendre_nonfinite_derivatives_real

    pure function legendre_highest_derivative(n) result(value)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp) :: value

        real(dp) :: factor
        integer :: k

        value = 1.0_dp
        do k = 1, n
            factor = 0.5_dp * (real(n, dp) + real(k, dp))
            if (value > huge(value) / factor) then
                value = positive_infinity(value)
                return
            end if
            value = value * factor
        end do
    end function legendre_highest_derivative

    pure subroutine legendre_jet_real(n, zjet, coefficients)
        integer, intent(in) :: n !! nonnegative polynomial degree
        real(dp), intent(in) :: zjet(0:) !! Taylor coefficients of z+t
        real(dp), intent(out) :: coefficients(0:) !! Taylor coefficients of P_n(z+t)

        real(dp) :: current(0:ubound(coefficients, 1))
        integer :: degree
        real(dp) :: next_value(0:ubound(coefficients, 1))
        real(dp) :: previous(0:ubound(coefficients, 1))
        real(dp) :: product(0:ubound(coefficients, 1))

        coefficients = 0.0_dp
        previous = 0.0_dp
        previous(0) = 1.0_dp
        if (n == 0) then
            coefficients = previous
            return
        end if

        current = zjet
        if (n == 1) then
            coefficients = current
            return
        end if

        do degree = 2, n
            call legendre_argument_product_real(zjet(0), current, product)
            next_value = (real(2 * degree - 1, dp) * product - &
                real(degree - 1, dp) * previous) / real(degree, dp)
            previous = current
            current = next_value
        end do
        coefficients = current
    end subroutine legendre_jet_real

    pure subroutine legendre_jet_complex(n, zjet, coefficients)
        integer, intent(in) :: n !! nonnegative polynomial degree
        complex(dp), intent(in) :: zjet(0:) !! Taylor coefficients of z+t
        complex(dp), intent(out) :: coefficients(0:) !! Taylor coefficients of P_n(z+t)

        complex(dp) :: current(0:ubound(coefficients, 1))
        integer :: degree
        complex(dp) :: next_value(0:ubound(coefficients, 1))
        complex(dp) :: previous(0:ubound(coefficients, 1))
        complex(dp) :: product(0:ubound(coefficients, 1))

        coefficients = cmplx(0.0_dp, 0.0_dp, dp)
        previous = cmplx(0.0_dp, 0.0_dp, dp)
        previous(0) = cmplx(1.0_dp, 0.0_dp, dp)
        if (n == 0) then
            coefficients = previous
            return
        end if

        current = zjet
        if (n == 1) then
            coefficients = current
            return
        end if

        do degree = 2, n
            call legendre_argument_product_complex(zjet(0), current, product)
            next_value = (real(2 * degree - 1, dp) * product - &
                real(degree - 1, dp) * previous) / real(degree, dp)
            previous = current
            current = next_value
        end do
        coefficients = current
    end subroutine legendre_jet_complex

    pure subroutine legendre_argument_product_real(z, input, product)
        real(dp), intent(in) :: z !! finite expansion point
        real(dp), intent(in) :: input(0:) !! Taylor coefficients of the multiplicand
        real(dp), intent(out) :: product(0:) !! Taylor coefficients of (z+t)*input

        integer :: k

        product = 0.0_dp
        product(0) = z * input(0)
        do k = 1, ubound(product, 1)
            product(k) = z * input(k) + input(k - 1)
        end do
    end subroutine legendre_argument_product_real

    pure subroutine legendre_argument_product_complex(z, input, product)
        complex(dp), intent(in) :: z !! finite expansion point
        complex(dp), intent(in) :: input(0:) !! Taylor coefficients of the multiplicand
        complex(dp), intent(out) :: product(0:) !! Taylor coefficients of (z+t)*input

        integer :: k

        product = cmplx(0.0_dp, 0.0_dp, dp)
        product(0) = z * input(0)
        do k = 1, ubound(product, 1)
            product(k) = z * input(k) + input(k - 1)
        end do
    end subroutine legendre_argument_product_complex

    pure subroutine assoc_legendre_jet(n, m, zjet, branch_cut, norm, coefficients)
        integer, intent(in) :: n !! nonnegative integer degree
        integer, intent(in) :: m !! integer order
        complex(dp), intent(in) :: zjet(0:) !! Taylor coefficients of the argument
        integer, intent(in) :: branch_cut !! branch convention, 2 or 3
        logical, intent(in) :: norm !! use spherical normalization when true
        complex(dp), intent(out) :: coefficients(0:) !! Taylor coefficients of the result

        complex(dp) :: current(0:ubound(zjet, 1))
        integer :: degree
        real(dp) :: factor
        real(dp) :: log_factor
        integer :: mm
        complex(dp) :: next_value(0:ubound(zjet, 1))
        complex(dp) :: previous(0:ubound(zjet, 1))
        complex(dp) :: q1(0:ubound(zjet, 1))
        complex(dp) :: q2(0:ubound(zjet, 1))
        complex(dp) :: r1(0:ubound(zjet, 1))
        complex(dp) :: r2(0:ubound(zjet, 1))
        complex(dp) :: radial(0:ubound(zjet, 1))
        complex(dp) :: product(0:ubound(zjet, 1))

        coefficients = cmplx(0.0_dp, 0.0_dp, dp)
        if (n < 0 .or. m == -huge(m) - 1) return
        mm = abs(m)
        if (mm > n) return

        q1 = -zjet
        q2 = zjet
        if (branch_cut == 2) then
            q1(0) = q1(0) + cmplx(1.0_dp, 0.0_dp, dp)
            q2(0) = q2(0) + cmplx(1.0_dp, 0.0_dp, dp)
        else
            q1 = zjet
            q2 = zjet
            q1(0) = q1(0) - cmplx(1.0_dp, 0.0_dp, dp)
            q2(0) = q2(0) + cmplx(1.0_dp, 0.0_dp, dp)
        end if
        call jet_sqrt(q1, r1)
        call jet_sqrt(q2, r2)
        call jet_multiply(r1, r2, radial)

        current = cmplx(0.0_dp, 0.0_dp, dp)
        current(0) = cmplx(1.0_dp, 0.0_dp, dp)
        do degree = 1, mm
            call jet_multiply(radial, current, product)
            factor = real(2 * degree - 1, dp)
            if (branch_cut == 2) factor = -factor
            current = factor * product
        end do

        if (n == mm) then
            coefficients = current
        else
            previous = current
            call jet_multiply(zjet, previous, product)
            current = real(2 * mm + 1, dp) * product
            if (n == mm + 1) then
                coefficients = current
            else
                do degree = mm + 2, n
                    call jet_multiply(zjet, current, product)
                    next_value = (real(2 * degree - 1, dp) * product - &
                        real(degree + mm - 1, dp) * previous) / real(degree - mm, dp)
                    previous = current
                    current = next_value
                end do
                coefficients = current
            end if
        end if

        if (norm) then
            log_factor = 0.5_dp * (log(real(2 * n + 1, dp)) - log(2.0_dp) + &
                gammaln(real(n - mm + 1, dp)) - gammaln(real(n + mm + 1, dp)))
            factor = exp(log_factor)
            if (m < 0 .and. branch_cut == 2 .and. mod(mm, 2) == 1) factor = -factor
            coefficients = factor * coefficients
        else if (m < 0) then
            log_factor = gammaln(real(n - mm + 1, dp)) - gammaln(real(n + mm + 1, dp))
            factor = exp(log_factor)
            if (branch_cut == 2 .and. mod(mm, 2) == 1) factor = -factor
            coefficients = factor * coefficients
        end if
    end subroutine assoc_legendre_jet

    pure subroutine jet_multiply(a, b, product)
        complex(dp), intent(in) :: a(0:) !! first Taylor-coefficient sequence
        complex(dp), intent(in) :: b(0:) !! second Taylor-coefficient sequence
        complex(dp), intent(out) :: product(0:) !! Taylor coefficients of a*b

        integer :: j
        integer :: k

        product = cmplx(0.0_dp, 0.0_dp, dp)
        do k = 0, ubound(product, 1)
            do j = 0, k
                product(k) = product(k) + a(j) * b(k - j)
            end do
        end do
    end subroutine jet_multiply

    pure subroutine jet_sqrt(q, root)
        complex(dp), intent(in) :: q(0:) !! Taylor coefficients of the radicand
        complex(dp), intent(out) :: root(0:) !! principal-square-root Taylor coefficients

        complex(dp) :: accumulated
        integer :: j
        integer :: k
        real(dp) :: nan_value

        root = cmplx(0.0_dp, 0.0_dp, dp)
        root(0) = sqrt(q(0))
        if (ubound(root, 1) == 0) return
        if (abs(root(0)) == 0.0_dp) then
            nan_value = quiet_nan(real(q(0), dp))
            root(1:) = cmplx(nan_value, nan_value, dp)
            return
        end if

        do k = 1, ubound(root, 1)
            accumulated = cmplx(0.0_dp, 0.0_dp, dp)
            do j = 1, k - 1
                accumulated = accumulated + root(j) * root(k - j)
            end do
            root(k) = (q(k) - accumulated) / (2.0_dp * root(0))
        end do
    end subroutine jet_sqrt

    pure subroutine cosine_jet(theta, coefficients)
        real(dp), intent(in) :: theta !! expansion point in radians
        complex(dp), intent(out) :: coefficients(0:) !! Taylor coefficients of cos(theta+t)

        real(dp) :: factorial_value
        integer :: k
        real(dp) :: derivative_value

        factorial_value = 1.0_dp
        do k = 0, ubound(coefficients, 1)
            if (k > 0) factorial_value = factorial_value * real(k, dp)
            select case (mod(k, 4))
            case (0)
                derivative_value = cos(theta)
            case (1)
                derivative_value = -sin(theta)
            case (2)
                derivative_value = -cos(theta)
            case default
                derivative_value = sin(theta)
            end select
            coefficients(k) = cmplx(derivative_value / factorial_value, 0.0_dp, dp)
        end do
    end subroutine cosine_jet

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
        complex(dp), intent(in) :: z !! reference value used to create a quiet NaN
        complex(dp) :: value

        real(dp) :: nan_value

        nan_value = quiet_nan(real(z, dp))
        value = cmplx(nan_value, nan_value, dp)
    end function complex_nan

    pure function integer_degree_lpmv(n, m, x) result(value)
        integer, intent(in) :: n !! nonnegative integer degree
        integer, intent(in) :: m !! integer order
        real(dp), intent(in) :: x !! real Ferrers argument, -1 <= x <= 1
        real(dp) :: value

        integer :: k
        integer :: mm
        real(dp) :: current
        real(dp) :: log_factor
        real(dp) :: next_value
        real(dp) :: previous
        real(dp) :: radial

        mm = abs(m)
        if (m < 0 .and. mm > n) then
            value = quiet_nan(x)
            return
        end if
        if (m >= 0 .and. mm > n) then
            value = 0.0_dp
            return
        end if

        radial = sqrt(max(0.0_dp, 1.0_dp - x * x))
        current = 1.0_dp
        if (mm > 0) then
            do k = 1, mm
                current = -real(2 * k - 1, dp) * radial * current
            end do
        end if

        if (n == mm) then
            value = current
        else
            previous = current
            current = real(2 * mm + 1, dp) * x * previous
            if (n == mm + 1) then
                value = current
            else
                do k = mm + 2, n
                    next_value = (real(2 * k - 1, dp) * x * current - &
                        real(k + mm - 1, dp) * previous) / real(k - mm, dp)
                    previous = current
                    current = next_value
                end do
                value = current
            end if
        end if

        if (m < 0) then
            log_factor = gammaln(real(n - mm + 1, dp)) - gammaln(real(n + mm + 1, dp))
            if (mod(mm, 2) == 1) then
                value = -value * exp(log_factor)
            else
                value = value * exp(log_factor)
            end if
        end if
    end function integer_degree_lpmv

    pure logical function is_representable_integer_degree(v) result(is_integer)
        real(dp), intent(in) :: v !! finite real degree to classify

        real(dp), parameter :: max_integer_degree = real(huge(0) - 1, dp)

        is_integer = abs(v) <= max_integer_degree .and. v == anint(v)
    end function is_representable_integer_degree

    pure function signed_exp(log_abs, sign_value) result(value)
        real(dp), intent(in) :: log_abs !! logarithm of the desired magnitude
        real(dp), intent(in) :: sign_value !! sign of the desired result
        real(dp) :: value

        if (ieee_is_nan(log_abs) .or. ieee_is_nan(sign_value)) then
            value = quiet_nan(log_abs)
        else if (log_abs > log(huge(1.0_dp))) then
            value = sign(positive_infinity(0.0_dp), sign_value)
        else if (log_abs == -positive_infinity(0.0_dp)) then
            value = sign(0.0_dp, sign_value)
        else
            value = sign(exp(log_abs), sign_value)
        end if
    end function signed_exp

end module scifort_associated_legendre
