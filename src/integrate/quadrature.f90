! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Adapted portions: SciPy 1.17.0 scipy/integrate/_quadrature.py numerical formulas
! for Simpson end correction and Newton-Cotes built-in coefficients.
! Retained BSD terms: THIRD_PARTY_LICENSES.md; details: CODE_PROVENANCE.md.

module scifort_integrate_quadrature
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_orthogonal_roots, only : orthogonal_roots_status_success, roots_legendre
    implicit none
    private

    integer, parameter, public :: integrate_status_success = 0
    integer, parameter, public :: integrate_status_invalid_shape = 1
    integer, parameter, public :: integrate_status_invalid_argument = 2
    integer, parameter, public :: integrate_status_numerical_failure = 3

    public :: cumulative_simpson
    public :: cumulative_trapezoid
    public :: fixed_quad
    public :: newton_cotes
    public :: romb
    public :: simpson
    public :: trapezoid

    abstract interface
        function scalar_integrand(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x !! integration point
            real(dp) :: value
        end function scalar_integrand
    end interface

    interface newton_cotes
        module procedure newton_cotes_order
        module procedure newton_cotes_nodes
    end interface newton_cotes

contains

    function trapezoid(y, x, dx, status) result(value)
        real(dp), intent(in) :: y(:) !! sampled function values in integration order
        real(dp), intent(in), optional :: x(:) !! sample coordinates; same size as y
        real(dp), intent(in), optional :: dx !! uniform spacing when x is absent; default 1
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input
        real(dp) :: value

        integer :: i
        real(dp) :: spacing

        if (present(status)) status = integrate_status_success
        if (present(x)) then
            if (size(x) /= size(y)) then
                value = quiet_nan(0.0_dp)
                if (present(status)) status = integrate_status_invalid_shape
                return
            end if
            value = 0.0_dp
            do i = 1, size(y) - 1
                value = value + 0.5_dp * (x(i + 1) - x(i)) * (y(i) + y(i + 1))
            end do
        else
            spacing = 1.0_dp
            if (present(dx)) spacing = dx
            if (size(y) <= 1) then
                value = 0.0_dp
            else
                value = spacing * (0.5_dp * y(1) + sum(y(2:size(y) - 1)) + 0.5_dp * y(size(y)))
            end if
        end if
    end function trapezoid

    function cumulative_trapezoid(y, x, dx, initial, status) result(values)
        real(dp), intent(in) :: y(:) !! sampled function values in integration order
        real(dp), intent(in), optional :: x(:) !! sample coordinates; same size as y
        real(dp), intent(in), optional :: dx !! uniform spacing when x is absent; default 1
        real(dp), intent(in), optional :: initial !! optional initial value; must be zero
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input
        real(dp), allocatable :: values(:)

        integer :: i
        integer :: offset
        integer :: n
        real(dp) :: spacing
        real(dp) :: term

        if (present(status)) status = integrate_status_success
        n = size(y)
        if (n == 0) then
            allocate(values(0))
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if
        if (present(initial)) then
            if (initial /= 0.0_dp) then
                allocate(values(0))
                if (present(status)) status = integrate_status_invalid_argument
                return
            end if
        end if
        if (present(x)) then
            if (size(x) /= n) then
                allocate(values(0))
                if (present(status)) status = integrate_status_invalid_shape
                return
            end if
        end if

        offset = 0
        if (present(initial)) offset = 1
        allocate(values(n - 1 + offset))
        if (offset == 1) values(1) = 0.0_dp
        if (n == 1) return

        spacing = 1.0_dp
        if (present(dx)) spacing = dx
        do i = 1, n - 1
            if (present(x)) then
                term = 0.5_dp * (x(i + 1) - x(i)) * (y(i) + y(i + 1))
            else
                term = 0.5_dp * spacing * (y(i) + y(i + 1))
            end if
            if (i == 1) then
                values(i + offset) = term
            else
                values(i + offset) = values(i - 1 + offset) + term
            end if
        end do
    end function cumulative_trapezoid

    function simpson(y, x, dx, status) result(value)
        real(dp), intent(in) :: y(:) !! sampled function values in integration order
        real(dp), intent(in), optional :: x(:) !! sample coordinates; same size as y
        real(dp), intent(in), optional :: dx !! uniform spacing when x is absent; default 1
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input
        real(dp) :: value

        integer :: i
        integer :: n
        real(dp) :: h0
        real(dp) :: h1
        real(dp) :: spacing
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: eta

        if (present(status)) status = integrate_status_success
        n = size(y)
        if (n == 0) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if
        if (present(x)) then
            if (size(x) /= n) then
                value = quiet_nan(0.0_dp)
                if (present(status)) status = integrate_status_invalid_shape
                return
            end if
        end if
        if (n == 1) then
            value = 0.0_dp
            return
        end if
        if (n == 2) then
            value = trapezoid(y, x, dx)
            return
        end if

        spacing = 1.0_dp
        if (present(dx)) spacing = dx
        value = 0.0_dp
        if (mod(n, 2) == 1) then
            do i = 1, n - 2, 2
                call add_simpson_panel(i, y, x, spacing, value)
            end do
        else
            do i = 1, n - 3, 2
                call add_simpson_panel(i, y, x, spacing, value)
            end do
            if (present(x)) then
                h0 = x(n - 1) - x(n - 2)
                h1 = x(n) - x(n - 1)
            else
                h0 = spacing
                h1 = spacing
            end if
            alpha = 0.0_dp
            if (h0 + h1 /= 0.0_dp) then
                alpha = (2.0_dp * h1 * h1 + 3.0_dp * h0 * h1) / (6.0_dp * (h0 + h1))
            end if
            beta = 0.0_dp
            if (h0 /= 0.0_dp) then
                beta = (h1 * h1 + 3.0_dp * h0 * h1) / (6.0_dp * h0)
            end if
            eta = 0.0_dp
            if (h0 /= 0.0_dp .and. h0 + h1 /= 0.0_dp) then
                eta = h1 * h1 * h1 / (6.0_dp * h0 * (h0 + h1))
            end if
            value = value + alpha * y(n) + beta * y(n - 1) - eta * y(n - 2)
        end if
    end function simpson

    function cumulative_simpson(y, x, dx, initial, status) result(values)
        real(dp), intent(in) :: y(:) !! sampled function values in integration order
        real(dp), intent(in), optional :: x(:) !! strictly increasing coordinates; same size as y
        real(dp), intent(in), optional :: dx !! uniform spacing when x is absent; default 1
        real(dp), intent(in), optional :: initial !! optional value prepended and added to all integrals
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input
        real(dp), allocatable :: values(:)

        integer :: i
        integer :: n
        integer :: offset
        real(dp), allocatable :: interval_integrals(:)
        real(dp), allocatable :: widths(:)
        real(dp) :: base
        real(dp) :: spacing

        if (present(status)) status = integrate_status_success
        n = size(y)
        if (n == 0) then
            allocate(values(0))
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if
        if (present(x)) then
            if (size(x) /= n) then
                allocate(values(0))
                if (present(status)) status = integrate_status_invalid_shape
                return
            end if
            if (n > 1 .and. any(x(2:n) <= x(1:n - 1))) then
                allocate(values(0))
                if (present(status)) status = integrate_status_invalid_argument
                return
            end if
        end if

        base = 0.0_dp
        offset = 0
        if (present(initial)) then
            base = initial
            offset = 1
        end if
        allocate(values(n - 1 + offset))
        if (offset == 1) values(1) = base
        if (n == 1) return

        spacing = 1.0_dp
        if (present(dx)) spacing = dx
        if (n == 2) then
            if (present(x)) then
                values(1 + offset) = base + 0.5_dp * (x(2) - x(1)) * (y(1) + y(2))
            else
                values(1 + offset) = base + 0.5_dp * spacing * (y(1) + y(2))
            end if
            return
        end if

        allocate(widths(n - 1), interval_integrals(n - 1))
        if (present(x)) then
            widths = x(2:n) - x(1:n - 1)
        else
            widths = spacing
        end if

        if (.not. present(x) .and. spacing == 0.0_dp) then
            interval_integrals = 0.0_dp
        else
            do i = 1, n - 2
                interval_integrals(i) = simpson_subinterval_left(y(i), y(i + 1), y(i + 2), &
                    widths(i), widths(i + 1))
            end do
            do i = 2, n - 1, 2
                interval_integrals(i) = simpson_subinterval_right(y(i - 1), y(i), y(i + 1), &
                    widths(i - 1), widths(i))
            end do
            if (mod(n - 1, 2) == 1) then
                interval_integrals(n - 1) = simpson_subinterval_right(y(n - 2), y(n - 1), y(n), &
                    widths(n - 2), widths(n - 1))
            end if
        end if

        do i = 1, n - 1
            if (i == 1) then
                values(i + offset) = base + interval_integrals(i)
            else
                values(i + offset) = values(i - 1 + offset) + interval_integrals(i)
            end if
        end do
    end function cumulative_simpson

    function romb(y, dx, status) result(value)
        real(dp), intent(in) :: y(:) !! 2**k+1 equally spaced samples
        real(dp), intent(in), optional :: dx !! sample spacing; default 1
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid sample count
        real(dp) :: value

        integer :: i
        integer :: j
        integer :: k
        integer :: ninterval
        integer :: step
        integer :: start
        real(dp) :: h
        real(dp) :: spacing
        real(dp), allocatable :: previous(:)
        real(dp), allocatable :: current(:)

        if (present(status)) status = integrate_status_success
        ninterval = size(y) - 1
        if (ninterval < 1 .or. iand(ninterval, ninterval - 1) /= 0) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if

        k = 0
        step = ninterval
        do while (step > 1)
            step = step / 2
            k = k + 1
        end do
        allocate(previous(k + 1), current(k + 1))
        previous = 0.0_dp
        current = 0.0_dp
        spacing = 1.0_dp
        if (present(dx)) spacing = dx
        h = real(ninterval, dp) * spacing
        previous(1) = 0.5_dp * h * (y(1) + y(size(y)))
        step = ninterval
        start = ninterval + 1

        do i = 1, k
            start = start / 2
            current = 0.0_dp
            current(1) = 0.5_dp * previous(1)
            j = start + 1
            do while (j < size(y))
                current(1) = current(1) + 0.5_dp * h * y(j)
                j = j + step
            end do
            do j = 1, i
                current(j + 1) = current(j) + (current(j) - previous(j)) / &
                    (real(2 ** (2 * j), dp) - 1.0_dp)
            end do
            previous = current
            h = 0.5_dp * h
            step = step / 2
        end do
        value = previous(k + 1)
    end function romb

    function fixed_quad(f, a, b, n, status) result(value)
        procedure(scalar_integrand) :: f !! scalar integrand function
        real(dp), intent(in) :: a !! finite lower integration limit
        real(dp), intent(in) :: b !! finite upper integration limit
        integer, intent(in), optional :: n !! Gauss-Legendre order; default 5
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid input
        real(dp) :: value

        integer :: i
        integer :: info
        integer :: order
        real(dp) :: midpoint
        real(dp) :: halfwidth
        real(dp), allocatable :: nodes(:)
        real(dp), allocatable :: weights(:)

        if (present(status)) status = integrate_status_success
        order = 5
        if (present(n)) order = n
        if (order < 1 .or. .not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_invalid_argument
            return
        end if

        call roots_legendre(order, nodes, weights, status=info)
        if (info /= orthogonal_roots_status_success) then
            value = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_numerical_failure
            return
        end if
        midpoint = 0.5_dp * (a + b)
        halfwidth = 0.5_dp * (b - a)
        value = 0.0_dp
        do i = 1, order
            value = value + weights(i) * f(midpoint + halfwidth * nodes(i))
        end do
        value = halfwidth * value
    end function fixed_quad

    subroutine newton_cotes_order(n, weights, error_coeff, status)
        integer, intent(in) :: n !! rule order; returns n+1 equally spaced weights
        real(dp), allocatable, intent(out) :: weights(:) !! Newton-Cotes weights on relative nodes 0:n
        real(dp), intent(out) :: error_coeff !! leading derivative error coefficient
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid order

        real(dp), allocatable :: nodes(:)
        integer :: i

        if (n < 1) then
            allocate(weights(0))
            error_coeff = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_invalid_argument
            return
        end if
        allocate(nodes(n + 1))
        do i = 0, n
            nodes(i + 1) = real(i, dp)
        end do
        call newton_cotes_impl(nodes, .true., weights, error_coeff, status)
    end subroutine newton_cotes_order

    subroutine newton_cotes_nodes(rn, weights, error_coeff, equal, status)
        real(dp), intent(in) :: rn(:) !! relative sample positions, starting at 0 and ending at n
        real(dp), allocatable, intent(out) :: weights(:) !! Newton-Cotes weights for rn
        real(dp), intent(out) :: error_coeff !! leading derivative error coefficient
        logical, intent(in), optional :: equal !! treat nodes as equally spaced; default false
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid nodes

        logical :: equally_spaced

        equally_spaced = .false.
        if (present(equal)) equally_spaced = equal
        if (.not. equally_spaced .and. size(rn) > 1) then
            equally_spaced = all(abs((rn(2:size(rn)) - rn(1:size(rn) - 1)) - 1.0_dp) &
                <= 8.0_dp * epsilon(1.0_dp))
        end if
        call newton_cotes_impl(rn, equally_spaced, weights, error_coeff, status)
    end subroutine newton_cotes_nodes

    subroutine newton_cotes_impl(rn, equal, weights, error_coeff, status)
        real(dp), intent(in) :: rn(:) !! relative sample positions
        logical, intent(in) :: equal !! true when nodes are equally spaced
        real(dp), allocatable, intent(out) :: weights(:) !! Newton-Cotes weights
        real(dp), intent(out) :: error_coeff !! leading derivative error coefficient
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        integer :: i
        integer :: info
        integer :: iteration
        integer :: n
        integer :: power
        real(dp) :: coefficient
        real(dp) :: order_real
        real(dp), allocatable :: correction(:)
        real(dp), allocatable :: matrix(:, :)
        real(dp), allocatable :: moments(:)
        real(dp), allocatable :: nodes(:)
        real(dp), allocatable :: residual(:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: work_matrix(:, :)

        if (present(status)) status = integrate_status_success
        n = size(rn) - 1
        if (n < 1 .or. rn(1) /= 0.0_dp .or. rn(size(rn)) /= real(n, dp)) then
            allocate(weights(0))
            error_coeff = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_invalid_argument
            return
        end if

        if (equal .and. n <= 14) then
            call builtin_newton_cotes(n, weights, error_coeff)
            return
        end if

        order_real = real(n, dp)
        allocate(matrix(n + 1, n + 1), moments(n + 1), nodes(n + 1), &
            residual(n + 1), rhs(n + 1), work_matrix(n + 1, n + 1), &
            weights(n + 1), correction(n + 1))
        nodes = 2.0_dp * rn / order_real - 1.0_dp
        do i = 0, n
            matrix(i + 1, :) = nodes ** i
            if (mod(i, 2) == 0) then
                moments(i + 1) = 2.0_dp / real(i + 1, dp)
            else
                moments(i + 1) = 0.0_dp
            end if
        end do

        work_matrix = matrix
        rhs = moments
        call solve_linear_system(work_matrix, rhs, weights, info)
        if (info /= 0) then
            weights = quiet_nan(0.0_dp)
            error_coeff = quiet_nan(0.0_dp)
            if (present(status)) status = integrate_status_numerical_failure
            return
        end if
        do iteration = 1, 2
            residual = moments - matmul(matrix, weights)
            work_matrix = matrix
            rhs = residual
            call solve_linear_system(work_matrix, rhs, correction, info)
            if (info /= 0) then
                weights = quiet_nan(0.0_dp)
                error_coeff = quiet_nan(0.0_dp)
                if (present(status)) status = integrate_status_numerical_failure
                return
            end if
            weights = weights + correction
        end do
        weights = 0.5_dp * order_real * weights

        if (equal .and. mod(n, 2) == 0) then
            power = n + 2
            coefficient = order_real / real(n + 3, dp)
        else
            power = n + 1
            coefficient = order_real / real(n + 2, dp)
        end if
        coefficient = coefficient - sum((rn / order_real) ** power * weights)
        error_coeff = coefficient * exp(real(power, dp) * log(order_real) - &
            log_gamma(real(power + 1, dp)))
    end subroutine newton_cotes_impl

    subroutine add_simpson_panel(i, y, x, spacing, value)
        integer, intent(in) :: i !! first index of the three-point panel
        real(dp), intent(in) :: y(:) !! sampled values
        real(dp), intent(in), optional :: x(:) !! optional sample coordinates
        real(dp), intent(in) :: spacing !! uniform spacing when x is absent
        real(dp), intent(inout) :: value !! accumulated integral

        real(dp) :: h0
        real(dp) :: h0_div_h1
        real(dp) :: h1
        real(dp) :: hsum
        real(dp) :: hsum_over_hprod
        real(dp) :: inv_h0_div_h1

        if (present(x)) then
            h0 = x(i + 1) - x(i)
            h1 = x(i + 2) - x(i + 1)
            hsum = h0 + h1
            h0_div_h1 = 0.0_dp
            if (h1 /= 0.0_dp) h0_div_h1 = h0 / h1
            inv_h0_div_h1 = 0.0_dp
            if (h0_div_h1 /= 0.0_dp) inv_h0_div_h1 = 1.0_dp / h0_div_h1
            hsum_over_hprod = 0.0_dp
            if (h0 * h1 /= 0.0_dp) hsum_over_hprod = hsum / (h0 * h1)
            value = value + hsum / 6.0_dp * ( &
                y(i) * (2.0_dp - inv_h0_div_h1) + &
                y(i + 1) * hsum * hsum_over_hprod + &
                y(i + 2) * (2.0_dp - h0_div_h1))
        else
            value = value + spacing / 3.0_dp * (y(i) + 4.0_dp * y(i + 1) + y(i + 2))
        end if
    end subroutine add_simpson_panel

    pure function simpson_subinterval_left(f1, f2, f3, h1, h2) result(value)
        real(dp), intent(in) :: f1 !! value at left point
        real(dp), intent(in) :: f2 !! value at middle point
        real(dp), intent(in) :: f3 !! value at right point
        real(dp), intent(in) :: h1 !! left interval width
        real(dp), intent(in) :: h2 !! right interval width
        real(dp) :: value

        real(dp) :: ratio
        real(dp) :: total

        total = h1 + h2
        ratio = h1 / total
        value = h1 / 6.0_dp * ((3.0_dp - ratio) * f1 + &
            (3.0_dp + h1 * h1 / (h2 * total) + ratio) * f2 - &
            h1 * h1 / (h2 * total) * f3)
    end function simpson_subinterval_left

    pure function simpson_subinterval_right(f1, f2, f3, h1, h2) result(value)
        real(dp), intent(in) :: f1 !! value at left point
        real(dp), intent(in) :: f2 !! value at middle point
        real(dp), intent(in) :: f3 !! value at right point
        real(dp), intent(in) :: h1 !! left interval width
        real(dp), intent(in) :: h2 !! right interval width
        real(dp) :: value

        value = simpson_subinterval_left(f3, f2, f1, h2, h1)
    end function simpson_subinterval_right

    subroutine solve_linear_system(a, b, x, info)
        real(dp), intent(inout) :: a(:, :) !! square coefficient matrix, overwritten by elimination
        real(dp), intent(inout) :: b(:) !! right-hand side, overwritten by elimination
        real(dp), intent(out) :: x(:) !! solution vector
        integer, intent(out) :: info !! zero on success; nonzero for singular system

        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot
        real(dp) :: factor
        real(dp) :: scale
        real(dp) :: temp

        info = 0
        n = size(b)
        do k = 1, n - 1
            pivot = k
            do i = k + 1, n
                if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
            end do
            scale = max(1.0_dp, maxval(abs(a(:, k))))
            if (abs(a(pivot, k)) <= epsilon(1.0_dp) * scale) then
                info = 1
                x = quiet_nan(0.0_dp)
                return
            end if
            if (pivot /= k) then
                do j = k, n
                    temp = a(k, j)
                    a(k, j) = a(pivot, j)
                    a(pivot, j) = temp
                end do
                temp = b(k)
                b(k) = b(pivot)
                b(pivot) = temp
            end if
            do i = k + 1, n
                factor = a(i, k) / a(k, k)
                a(i, k) = 0.0_dp
                a(i, k + 1:n) = a(i, k + 1:n) - factor * a(k, k + 1:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do
        if (abs(a(n, n)) <= epsilon(1.0_dp) * max(1.0_dp, abs(a(n, n)))) then
            info = 1
            x = quiet_nan(0.0_dp)
            return
        end if
        x(n) = b(n) / a(n, n)
        do i = n - 1, 1, -1
            x(i) = (b(i) - sum(a(i, i + 1:n) * x(i + 1:n))) / a(i, i)
        end do
    end subroutine solve_linear_system

    subroutine builtin_newton_cotes(n, weights, error_coeff)
        integer, intent(in) :: n !! equally spaced rule order, 1 through 14
        real(dp), allocatable, intent(out) :: weights(:) !! exact-tabulated Newton-Cotes weights
        real(dp), intent(out) :: error_coeff !! exact-tabulated leading error coefficient

        allocate(weights(n + 1))
        select case (n)
        case (1)
            weights = [1.0_dp, 1.0_dp] / 2.0_dp
            error_coeff = -1.0_dp / 12.0_dp
        case (2)
            weights = [1.0_dp, 4.0_dp, 1.0_dp] / 3.0_dp
            error_coeff = -1.0_dp / 90.0_dp
        case (3)
            weights = 3.0_dp * [1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp] / 8.0_dp
            error_coeff = -3.0_dp / 80.0_dp
        case (4)
            weights = 2.0_dp * [7.0_dp, 32.0_dp, 12.0_dp, 32.0_dp, 7.0_dp] / 45.0_dp
            error_coeff = -8.0_dp / 945.0_dp
        case (5)
            weights = 5.0_dp * [19.0_dp, 75.0_dp, 50.0_dp, 50.0_dp, 75.0_dp, 19.0_dp] / 288.0_dp
            error_coeff = -275.0_dp / 12096.0_dp
        case (6)
            weights = [41.0_dp, 216.0_dp, 27.0_dp, 272.0_dp, 27.0_dp, 216.0_dp, 41.0_dp] / 140.0_dp
            error_coeff = -9.0_dp / 1400.0_dp
        case (7)
            weights = 7.0_dp * [751.0_dp, 3577.0_dp, 1323.0_dp, 2989.0_dp, &
                2989.0_dp, 1323.0_dp, 3577.0_dp, 751.0_dp] / 17280.0_dp
            error_coeff = -8183.0_dp / 518400.0_dp
        case (8)
            weights = 4.0_dp * [989.0_dp, 5888.0_dp, -928.0_dp, 10496.0_dp, -4540.0_dp, &
                10496.0_dp, -928.0_dp, 5888.0_dp, 989.0_dp] / 14175.0_dp
            error_coeff = -2368.0_dp / 467775.0_dp
        case (9)
            weights = 9.0_dp * [2857.0_dp, 15741.0_dp, 1080.0_dp, 19344.0_dp, 5778.0_dp, &
                5778.0_dp, 19344.0_dp, 1080.0_dp, 15741.0_dp, 2857.0_dp] / 89600.0_dp
            error_coeff = -4671.0_dp / 394240.0_dp
        case (10)
            weights = 5.0_dp * [16067.0_dp, 106300.0_dp, -48525.0_dp, 272400.0_dp, &
                -260550.0_dp, 427368.0_dp, -260550.0_dp, 272400.0_dp, -48525.0_dp, &
                106300.0_dp, 16067.0_dp] / 299376.0_dp
            error_coeff = -673175.0_dp / 163459296.0_dp
        case (11)
            weights = 11.0_dp * [2171465.0_dp, 13486539.0_dp, -3237113.0_dp, &
                25226685.0_dp, -9595542.0_dp, 15493566.0_dp, 15493566.0_dp, &
                -9595542.0_dp, 25226685.0_dp, -3237113.0_dp, 13486539.0_dp, &
                2171465.0_dp] / 87091200.0_dp
            error_coeff = -2224234463.0_dp / 237758976000.0_dp
        case (12)
            weights = [1364651.0_dp, 9903168.0_dp, -7587864.0_dp, 35725120.0_dp, &
                -51491295.0_dp, 87516288.0_dp, -87797136.0_dp, 87516288.0_dp, &
                -51491295.0_dp, 35725120.0_dp, -7587864.0_dp, 9903168.0_dp, &
                1364651.0_dp] / 5255250.0_dp
            error_coeff = -3012.0_dp / 875875.0_dp
        case (13)
            weights = 13.0_dp * [8181904909.0_dp, 56280729661.0_dp, -31268252574.0_dp, &
                156074417954.0_dp, -151659573325.0_dp, 206683437987.0_dp, &
                -43111992612.0_dp, -43111992612.0_dp, 206683437987.0_dp, &
                -151659573325.0_dp, 156074417954.0_dp, -31268252574.0_dp, &
                56280729661.0_dp, 8181904909.0_dp] / 402361344000.0_dp
            error_coeff = -2639651053.0_dp / 344881152000.0_dp
        case (14)
            weights = 7.0_dp * [90241897.0_dp, 710986864.0_dp, -770720657.0_dp, &
                3501442784.0_dp, -6625093363.0_dp, 12630121616.0_dp, &
                -16802270373.0_dp, 19534438464.0_dp, -16802270373.0_dp, &
                12630121616.0_dp, -6625093363.0_dp, 3501442784.0_dp, &
                -770720657.0_dp, 710986864.0_dp, 90241897.0_dp] / 2501928000.0_dp
            error_coeff = -3740727473.0_dp / 1275983280000.0_dp
        end select
    end subroutine builtin_newton_cotes

end module scifort_integrate_quadrature
