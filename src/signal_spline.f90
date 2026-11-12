! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Cubic and quadratic cardinal B-spline coefficients and evaluation.
! The forward/backward recurrences and mirror-end conventions follow
! scipy.signal._spline_filters (SciPy BSD-3-Clause). See CODE_PROVENANCE.md.
module scifort_signal_spline
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: spline_status_success = 0
    integer, parameter, public :: spline_status_invalid_input = 1

    interface cspline1d
        module procedure cubic_coeff_real
        module procedure cubic_coeff_complex
    end interface
    interface qspline1d
        module procedure quadratic_coeff_real
        module procedure quadratic_coeff_complex
    end interface
    interface cspline1d_eval
        module procedure cubic_eval_real
        module procedure cubic_eval_complex
        module procedure cubic_eval_real_scalar
        module procedure cubic_eval_complex_scalar
    end interface
    interface qspline1d_eval
        module procedure quadratic_eval_real
        module procedure quadratic_eval_complex
        module procedure quadratic_eval_real_scalar
        module procedure quadratic_eval_complex_scalar
    end interface
    public :: cspline1d, qspline1d, cspline1d_eval, qspline1d_eval

contains

    function cubic_coeff_real(signal, lamb, status) result(coeff)
        real(dp), intent(in) :: signal(:) !! sampled real signal
        real(dp), intent(in), optional :: lamb !! nonnegative cubic smoothing parameter, default zero
        integer, intent(out), optional :: status !! status, zero on success
        real(dp), allocatable :: coeff(:)
        complex(dp), allocatable :: tmp(:)
        integer :: code

        tmp = spline_coeff_core(cmplx(signal, 0.0_dp, dp), 3, lamb, code)
        coeff = real(tmp, dp)
        if (present(status)) status = code
    end function cubic_coeff_real

    function quadratic_coeff_real(signal, lamb, status) result(coeff)
        real(dp), intent(in) :: signal(:) !! sampled real signal
        real(dp), intent(in), optional :: lamb !! must be zero, as in SciPy
        integer, intent(out), optional :: status !! status, zero on success
        real(dp), allocatable :: coeff(:)
        complex(dp), allocatable :: tmp(:)
        integer :: code

        tmp = spline_coeff_core(cmplx(signal, 0.0_dp, dp), 2, lamb, code)
        coeff = real(tmp, dp)
        if (present(status)) status = code
    end function quadratic_coeff_real

    function cubic_coeff_complex(signal, lamb, status) result(coeff)
        complex(dp), intent(in) :: signal(:) !! sampled complex signal
        real(dp), intent(in), optional :: lamb !! nonnegative cubic smoothing parameter, default zero
        integer, intent(out), optional :: status !! status, zero on success
        complex(dp), allocatable :: coeff(:)
        integer :: code

        coeff = spline_coeff_core(signal, 3, lamb, code)
        if (present(status)) status = code
    end function cubic_coeff_complex

    function quadratic_coeff_complex(signal, lamb, status) result(coeff)
        complex(dp), intent(in) :: signal(:) !! sampled complex signal
        real(dp), intent(in), optional :: lamb !! must be zero, as in SciPy
        integer, intent(out), optional :: status !! status, zero on success
        complex(dp), allocatable :: coeff(:)
        integer :: code

        coeff = spline_coeff_core(signal, 2, lamb, code)
        if (present(status)) status = code
    end function quadratic_coeff_complex

    function spline_coeff_core(signal, degree, lamb, status) result(coeff)
        complex(dp), intent(in) :: signal(:) !! sampled values
        integer, intent(in) :: degree !! spline degree, two or three
        real(dp), intent(in), optional :: lamb !! cubic smoothing, or zero for interpolation
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: coeff(:), forward(:)
        complex(dp) :: total, power
        real(dp) :: lam, pole, scale
        integer :: n, i

        status = spline_status_invalid_input
        allocate(coeff(0))
        n = size(signal)
        if (n == 0) return
        if (any(.not. ieee_is_finite(real(signal, dp))) .or. &
            any(.not. ieee_is_finite(aimag(signal)))) return
        lam = 0.0_dp
        if (present(lamb)) lam = lamb
        if (.not. ieee_is_finite(lam) .or. lam < 0.0_dp) return
        if (degree == 2 .and. lam /= 0.0_dp) return
        if (degree /= 2 .and. degree /= 3) return

        if (lam > 0.0_dp) then
            ! SciPy's smoothed cubic recursion requires real omega and at
            ! least three knots. The limiting value at 1/144 is singular.
            if (degree /= 3 .or. n < 3 .or. lam <= 1.0_dp / 144.0_dp) return
            deallocate(coeff)
            coeff = smooth_cubic(signal, lam)
            if (size(coeff) == 0) return
        else
            pole = -2.0_dp + sqrt(3.0_dp)
            scale = 6.0_dp
            if (degree == 2) then
                pole = -3.0_dp + 2.0_dp * sqrt(2.0_dp)
                scale = 8.0_dp
            end if
            allocate(forward(n))
            power = cmplx(1.0_dp, 0.0_dp, dp)
            total = cmplx(0.0_dp, 0.0_dp, dp)
            do i = 1, n
                total = total + power * signal(i)
                power = power * pole
            end do
            forward(1) = signal(1) + pole * total
            if (n == 1) then
                ! SciPy's special K=1 path omits the final spline scale.
                deallocate(coeff)
                allocate(coeff(1))
                coeff(1) = pole / (pole - 1.0_dp) * forward(1)
            else
                do i = 2, n
                    forward(i) = signal(i) + pole * forward(i - 1)
                end do
                deallocate(coeff)
                allocate(coeff(n))
                coeff(n) = pole / (pole - 1.0_dp) * forward(n)
                do i = n - 1, 1, -1
                    coeff(i) = pole * (coeff(i + 1) - forward(i))
                end do
                coeff = scale * coeff
            end if
        end if
        if (any(.not. ieee_is_finite(real(coeff, dp))) .or. &
            any(.not. ieee_is_finite(aimag(coeff)))) then
            deallocate(coeff)
            allocate(coeff(0))
            return
        end if
        status = spline_status_success
    end function spline_coeff_core

    function smooth_cubic(signal, lam) result(coeff)
        complex(dp), intent(in) :: signal(:) !! sampled signal
        real(dp), intent(in) :: lam !! cubic smoothing coefficient above 1/144
        complex(dp), allocatable :: coeff(:)
        complex(dp), allocatable :: forward(:)
        complex(dp) :: edge2, edge1
        real(dp) :: xi, rho, omega, gain, a, b, c0, gamma
        real(dp) :: factor, rho2, den, sn, co
        integer :: n, i, k

        n = size(signal)
        xi = 1.0_dp - 96.0_dp * lam + 24.0_dp * lam * sqrt(3.0_dp + 144.0_dp * lam)
        omega = atan2(sqrt(144.0_dp * lam - 1.0_dp), sqrt(xi))
        rho = (24.0_dp * lam - 1.0_dp - sqrt(xi)) / (24.0_dp * lam)
        rho = rho * sqrt((48.0_dp * lam + 24.0_dp * lam * sqrt(3.0_dp + 144.0_dp * lam)) / xi)
        ! A few small smoothing values give |rho| > 1 in SciPy's closed
        ! form, making the recurrence explosive for long signals. Reject
        ! the unstable branch instead of producing overflow or Inf values.
        if (.not. ieee_is_finite(rho) .or. abs(rho) >= 1.0_dp) then
            allocate(coeff(0))
            return
        end if
        rho2 = rho*rho
        gain = 1.0_dp - 2.0_dp * rho * cos(omega) + rho2
        a = 2.0_dp * rho * cos(omega)
        b = rho2
        sn = sin(omega)
        co = cos(omega)
        allocate(forward(n), coeff(n))
        edge2 = gain / sn * sn * signal(1)
        edge1 = edge2 + gain / sn * rho * sin(2.0_dp * omega) * signal(2)
        factor = rho
        do k = 0, n - 1
            edge2 = edge2 + gain / sn * factor * sin(omega * real(k + 2, dp)) * signal(k + 1)
            edge1 = edge1 + gain / sn * (factor*rho) * sin(omega * real(k + 3, dp)) * signal(k + 1)
            factor = factor * rho
        end do
        forward(1) = edge2
        forward(2) = edge1
        do i = 3, n
            forward(i) = gain * signal(i) + a * forward(i - 1) - b * forward(i - 2)
        end do

        den = 1.0_dp - 2.0_dp*rho2*cos(2.0_dp*omega) + rho2*rho2
        c0 = gain*gain*(1.0_dp+rho2)/(1.0_dp-rho2)/den
        gamma = (1.0_dp-rho2)/(1.0_dp+rho2)/tan(omega)
        edge2 = cmplx(0.0_dp, 0.0_dp, dp)
        edge1 = cmplx(0.0_dp, 0.0_dp, dp)
        do k = 0, n - 1
            edge2 = edge2 + (hs(k) + hs(k+1)) * signal(n-k)
            edge1 = edge1 + (hs(k-1) + hs(k+2)) * signal(n-k)
        end do
        coeff(n-1) = edge1
        coeff(n) = edge2
        do i = n-2, 1, -1
            coeff(i) = gain * forward(i) + a * coeff(i+1) - b * coeff(i+2)
        end do
    contains
        real(dp) function hs(k) result(value)
            integer, intent(in) :: k !! offset from spline knot
            integer :: m
            m = abs(k)
            value = c0 * rho**m * (cos(omega*real(m,dp)) + gamma*sin(omega*real(m,dp)))
        end function hs
    end function smooth_cubic

    function cubic_eval_real(cj, newx, dx, x0, status) result(res)
        real(dp), intent(in) :: cj(:) !! cubic spline coefficients
        real(dp), intent(in) :: newx(:) !! positions at which to evaluate
        real(dp), intent(in), optional :: dx !! original knot spacing, default one
        real(dp), intent(in), optional :: x0 !! original knot origin, default zero
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: res(:)
        complex(dp), allocatable :: tmp(:)
        integer :: code
        tmp = spline_eval_core(cmplx(cj,0.0_dp,dp), newx, 3, dx, x0, code)
        res = real(tmp, dp)
        if (present(status)) status = code
    end function cubic_eval_real

    function quadratic_eval_real(cj, newx, dx, x0, status) result(res)
        real(dp), intent(in) :: cj(:) !! quadratic spline coefficients
        real(dp), intent(in) :: newx(:) !! positions at which to evaluate
        real(dp), intent(in), optional :: dx !! original knot spacing, default one
        real(dp), intent(in), optional :: x0 !! original knot origin, default zero
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: res(:)
        complex(dp), allocatable :: tmp(:)
        integer :: code
        tmp = spline_eval_core(cmplx(cj,0.0_dp,dp), newx, 2, dx, x0, code)
        res = real(tmp, dp)
        if (present(status)) status = code
    end function quadratic_eval_real

    function cubic_eval_complex(cj, newx, dx, x0, status) result(res)
        complex(dp), intent(in) :: cj(:) !! cubic spline coefficients
        real(dp), intent(in) :: newx(:) !! positions at which to evaluate
        real(dp), intent(in), optional :: dx !! original knot spacing, default one
        real(dp), intent(in), optional :: x0 !! original knot origin, default zero
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: res(:)
        integer :: code
        res = spline_eval_core(cj, newx, 3, dx, x0, code)
        if (present(status)) status = code
    end function cubic_eval_complex

    function quadratic_eval_complex(cj, newx, dx, x0, status) result(res)
        complex(dp), intent(in) :: cj(:) !! quadratic spline coefficients
        real(dp), intent(in) :: newx(:) !! positions at which to evaluate
        real(dp), intent(in), optional :: dx !! original knot spacing, default one
        real(dp), intent(in), optional :: x0 !! original knot origin, default zero
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: res(:)
        integer :: code
        res = spline_eval_core(cj, newx, 2, dx, x0, code)
        if (present(status)) status = code
    end function quadratic_eval_complex

    function spline_eval_core(cj, newx, degree, dx, x0, status) result(res)
        complex(dp), intent(in) :: cj(:) !! spline coefficients
        real(dp), intent(in) :: newx(:) !! physical knot locations for evaluation
        integer, intent(in) :: degree !! spline degree
        real(dp), intent(in), optional :: dx !! knot spacing, default one
        real(dp), intent(in), optional :: x0 !! origin, default zero
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: res(:)
        real(dp) :: spacing, origin, x, period
        integer :: i, n

        status = spline_status_invalid_input
        allocate(res(0))
        n = size(cj)
        if (n < 1) return
        if (any(.not. ieee_is_finite(real(cj, dp))) .or. any(.not. ieee_is_finite(aimag(cj)))) return
        if (any(.not. ieee_is_finite(newx))) return
        spacing = 1.0_dp
        origin = 0.0_dp
        if (present(dx)) spacing = dx
        if (present(x0)) origin = x0
        if (.not. ieee_is_finite(spacing) .or. spacing == 0.0_dp) return
        if (.not. ieee_is_finite(origin)) return
        deallocate(res)
        allocate(res(size(newx)))
        period = real(2 * (n - 1), dp)
        do i = 1, size(newx)
            x = (newx(i) - origin) / spacing
            if (.not. ieee_is_finite(x)) then
                deallocate(res)
                allocate(res(0))
                return
            end if
            if (n == 1) then
                x = 0.0_dp
            else
                x = modulo(x, period)
                if (x > real(n - 1,dp)) x = period - x
            end if
            res(i) = eval_at(cj, x, degree)
        end do
        status = spline_status_success
    end function spline_eval_core

    pure complex(dp) function eval_at(cj, x, degree) result(value)
        complex(dp), intent(in) :: cj(:) !! spline coefficients
        real(dp), intent(in) :: x !! reflected coordinate in [0, n-1]
        integer, intent(in) :: degree !! spline order
        real(dp) :: t
        integer :: k, j, count, lower, m, n

        n = size(cj)
        count = degree + 1
        lower = floor(x - real(degree + 1,dp)/2.0_dp) + 1
        value = cmplx(0.0_dp,0.0_dp,dp)
        do k = 0, count - 1
            j = lower + k
            m = max(1, min(n, j + 1))
            t = abs(x - real(j,dp))
            if (degree == 3) then
                if (t < 1.0_dp) then
                    value = value + cj(m) * (4.0_dp - 6.0_dp*t*t + 3.0_dp*t**3) / 6.0_dp
                else if (t < 2.0_dp) then
                    value = value + cj(m) * (2.0_dp-t)**3 / 6.0_dp
                end if
            else
                if (t < 0.5_dp) then
                    value = value + cj(m) * (0.75_dp-t*t)
                else if (t < 1.5_dp) then
                    value = value + cj(m) * 0.5_dp * (1.5_dp-t)**2
                end if
            end if
        end do
    end function eval_at

    function cubic_eval_real_scalar(cj, newx, dx, x0, status) result(res)
        real(dp), intent(in) :: cj(:) !! cubic spline coefficients
        real(dp), intent(in) :: newx !! evaluation coordinate
        real(dp), intent(in), optional :: dx !! knot spacing
        real(dp), intent(in), optional :: x0 !! knot origin
        integer, intent(out), optional :: status !! status of evaluation
        real(dp) :: res
        real(dp), allocatable :: values(:)
        integer :: code
        values = cubic_eval_real(cj, [newx], dx, x0, code)
        res = ieee_value(0.0_dp, ieee_quiet_nan)
        if (code == spline_status_success) res = values(1)
        if (present(status)) status = code
    end function cubic_eval_real_scalar

    function quadratic_eval_real_scalar(cj, newx, dx, x0, status) result(res)
        real(dp), intent(in) :: cj(:) !! quadratic spline coefficients
        real(dp), intent(in) :: newx !! evaluation coordinate
        real(dp), intent(in), optional :: dx !! knot spacing
        real(dp), intent(in), optional :: x0 !! knot origin
        integer, intent(out), optional :: status !! status of evaluation
        real(dp) :: res
        real(dp), allocatable :: values(:)
        integer :: code
        values = quadratic_eval_real(cj, [newx], dx, x0, code)
        res = ieee_value(0.0_dp, ieee_quiet_nan)
        if (code == spline_status_success) res = values(1)
        if (present(status)) status = code
    end function quadratic_eval_real_scalar

    function cubic_eval_complex_scalar(cj, newx, dx, x0, status) result(res)
        complex(dp), intent(in) :: cj(:) !! cubic spline coefficients
        real(dp), intent(in) :: newx !! evaluation coordinate
        real(dp), intent(in), optional :: dx !! knot spacing
        real(dp), intent(in), optional :: x0 !! knot origin
        integer, intent(out), optional :: status !! status of evaluation
        complex(dp) :: res
        complex(dp), allocatable :: values(:)
        integer :: code
        values = cubic_eval_complex(cj, [newx], dx, x0, code)
        res = cmplx(ieee_value(0.0_dp, ieee_quiet_nan),0.0_dp,dp)
        if (code == spline_status_success) res = values(1)
        if (present(status)) status = code
    end function cubic_eval_complex_scalar

    function quadratic_eval_complex_scalar(cj, newx, dx, x0, status) result(res)
        complex(dp), intent(in) :: cj(:) !! quadratic spline coefficients
        real(dp), intent(in) :: newx !! evaluation coordinate
        real(dp), intent(in), optional :: dx !! knot spacing
        real(dp), intent(in), optional :: x0 !! knot origin
        integer, intent(out), optional :: status !! status of evaluation
        complex(dp) :: res
        complex(dp), allocatable :: values(:)
        integer :: code
        values = quadratic_eval_complex(cj, [newx], dx, x0, code)
        res = cmplx(ieee_value(0.0_dp, ieee_quiet_nan),0.0_dp,dp)
        if (code == spline_status_success) res = values(1)
        if (present(status)) status = code
    end function quadratic_eval_complex_scalar
end module scifort_signal_spline
