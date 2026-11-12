! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Principal logarithms and real scalar powers of small dense matrices.
! Independently implemented inverse scaling-and-squaring with matrix square
! roots and the atanh series; this is not the Schur-Pade method in SciPy.
module scifort_linalg_logarithm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_matrix_functions, only : sqrtm, matrix_function_max_order
    use scifort_linalg_expm_frechet, only : expm_frechet
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter, public :: logm_max_order = matrix_function_max_order
    integer, parameter :: default_max_sqrt = 32
    integer, parameter :: default_max_terms = 80
    real(dp), parameter :: default_tol = 8.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: series_threshold = 0.25_dp

    public :: logm, logm_real, fractional_matrix_power, fractional_matrix_power_real

    interface logm
        module procedure logm_real_input
        module procedure logm_complex_input
    end interface

    interface fractional_matrix_power
        module procedure fractional_power_real_input
        module procedure fractional_power_complex_input
    end interface

contains

    subroutine logm_real_input(a, x, status, residual, square_roots, terms, &
            tol, max_sqrt, max_terms)
        real(dp), intent(in) :: a(:, :) !! square real matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! principal matrix logarithm (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(out), optional :: residual !! relative exp(log(A))-A reconstruction error
        integer, intent(out), optional :: square_roots !! number of matrix square roots
        integer, intent(out), optional :: terms !! odd-power series terms used
        real(dp), intent(in), optional :: tol !! relative series tolerance, default 8*epsilon
        integer, intent(in), optional :: max_sqrt !! maximum square roots, default 32
        integer, intent(in), optional :: max_terms !! maximum series terms, default 80

        call logm_core(cmplx(a, kind=dp), x, status, residual, square_roots, &
            terms, tol, max_sqrt, max_terms, real_input=.true.)
    end subroutine logm_real_input

    subroutine logm_complex_input(a, x, status, residual, square_roots, terms, &
            tol, max_sqrt, max_terms)
        complex(dp), intent(in) :: a(:, :) !! square complex matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! principal matrix logarithm (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(out), optional :: residual !! relative exp(log(A))-A reconstruction error
        integer, intent(out), optional :: square_roots !! number of matrix square roots
        integer, intent(out), optional :: terms !! odd-power series terms used
        real(dp), intent(in), optional :: tol !! relative series tolerance, default 8*epsilon
        integer, intent(in), optional :: max_sqrt !! maximum square roots, default 32
        integer, intent(in), optional :: max_terms !! maximum series terms, default 80

        call logm_core(a, x, status, residual, square_roots, terms, tol, max_sqrt, max_terms)
    end subroutine logm_complex_input

    subroutine logm_real(a, x, status, residual)
        real(dp), intent(in) :: a(:, :) !! square real matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! real principal logarithm, if it exists
        integer, intent(out) :: status !! nonzero when logarithm is non-real or undefined
        real(dp), intent(out), optional :: residual !! relative reconstruction error
        complex(dp), allocatable :: z(:, :)
        real(dp) :: magnitude

        x = 0.0_dp
        allocate(z(size(x, 1), size(x, 2)))
        call logm(a, z, status, residual)
        if (status /= linalg_status_success) return
        magnitude = max(1.0_dp, maxval(abs(z)))
        if (maxval(abs(aimag(z))) > 1.0e-9_dp*magnitude) then
            status = linalg_status_no_convergence
            return
        end if
        x = real(z, dp)
    end subroutine logm_real

    subroutine logm_core(a, x, status, residual, square_roots, terms, &
            tol, max_sqrt, max_terms, real_input)
        complex(dp), intent(in) :: a(:, :) !! square complex input (n,n)
        complex(dp), intent(out) :: x(:, :) !! principal logarithm (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(out), optional :: residual !! relative reconstruction error on scaled input
        integer, intent(out), optional :: square_roots !! successive square roots taken
        integer, intent(out), optional :: terms !! series terms evaluated
        real(dp), intent(in), optional :: tol !! relative series tolerance
        integer, intent(in), optional :: max_sqrt !! maximum repeated square roots
        integer, intent(in), optional :: max_terms !! maximum odd-power terms
        logical, intent(in), optional :: real_input !! prefer +imaginary on real negative spectrum
        complex(dp), allocatable :: b(:, :), root(:, :), eye(:, :), y(:, :), y2(:, :)
        complex(dp), allocatable :: power(:, :), sum_series(:, :), exp_b(:, :), zero(:, :)
        complex(dp), allocatable :: inverse_check(:, :), log_b(:, :)
        real(dp) :: alpha, tolerance, term_norm, err
        integer :: n, i, k, roots, root_limit, term_limit, code
        logical :: converged, snap_real

        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(residual)) residual = huge(1.0_dp)
        if (present(square_roots)) square_roots = 0
        if (present(terms)) terms = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (n < 1 .or. n > logm_max_order .or. size(a, 2) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= n) return
        tolerance = default_tol
        if (present(tol)) tolerance = tol
        root_limit = default_max_sqrt
        if (present(max_sqrt)) root_limit = max_sqrt
        term_limit = default_max_terms
        if (present(max_terms)) term_limit = max_terms
        if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp .or. &
                root_limit < 0 .or. root_limit > 50 .or. term_limit < 1) return
        if (.not. finite_matrix(a)) return
        alpha = maxval(abs(a))
        if (alpha == 0.0_dp) then
            status = linalg_status_singular
            return
        end if
        if (.not. ieee_is_finite(alpha)) then
            status = linalg_status_no_convergence
            return
        end if
        allocate(b(n, n), root(n, n), eye(n, n), y(n, n), y2(n, n), &
            power(n, n), sum_series(n, n), exp_b(n, n), zero(n, n), &
            inverse_check(n, n), log_b(n, n))
        eye = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            eye(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        ! Positive scalar scaling keeps subnormal/huge diagonal magnitudes
        ! out of the matrix-square-root and exponential recurrences.
        b = a/alpha
        call solve_dense(b, eye, inverse_check, code)
        if (code /= linalg_status_success) then
            status = code
            return
        end if
        snap_real = .false.
        if (present(real_input)) snap_real = real_input
        roots = 0
        do while (matrix_onenorm(b-eye) > series_threshold)
            if (roots >= root_limit) then
                status = linalg_status_no_convergence
                return
            end if
            if (roots == 0 .and. snap_real) then
                call sqrtm(real(b, dp), root, code)
            else
                call sqrtm(b, root, code)
            end if
            if (code /= linalg_status_success) then
                status = code
                return
            end if
            b = root
            roots = roots+1
        end do
        if (present(square_roots)) square_roots = roots
        ! Y = (B+I)^(-1) (B-I), so log(B) = 2*atanh(Y).
        ! The inverse is applied via partial-pivoted linear solves.
        call solve_dense(b+eye, b-eye, y, code)
        if (code /= linalg_status_success) then
            status = code
            return
        end if
        if (matrix_onenorm(y) >= 0.5_dp) then
            status = linalg_status_no_convergence
            return
        end if
        y2 = matmul(y, y)
        power = y
        sum_series = y
        converged = all(y == cmplx(0.0_dp, 0.0_dp, dp))
        if (present(terms)) terms = 1
        do k = 1, term_limit-1
            if (converged) exit
            power = matmul(power, y2)
            term_norm = matrix_onenorm(power)/real(2*k+1, dp)
            if (.not. ieee_is_finite(term_norm)) exit
            sum_series = sum_series+power/real(2*k+1, dp)
            if (.not. finite_matrix(sum_series)) exit
            if (present(terms)) terms = k+1
            if (term_norm <= tolerance*max(tiny(1.0_dp), matrix_onenorm(sum_series))) &
                converged = .true.
        end do
        if (.not. converged) then
            status = linalg_status_no_convergence
            return
        end if
        log_b = scale(2.0_dp, roots)*sum_series
        if (.not. finite_matrix(log_b)) then
            status = linalg_status_no_convergence
            return
        end if
        ! Verify a backward residual before accepting a result. Working with
        ! A/alpha prevents overflow from exp(log(alpha)).
        zero = cmplx(0.0_dp, 0.0_dp, dp)
        call expm_frechet(log_b, zero, exp_b, inverse_check, code)
        if (code /= linalg_status_success) then
            status = code
            return
        end if
        err = matrix_onenorm(exp_b-a/alpha)/max(tiny(1.0_dp), matrix_onenorm(a/alpha))
        if (present(residual)) residual = err
        if (.not. ieee_is_finite(err) .or. err > max(2.0e-9_dp, 100.0_dp*tolerance)) then
            status = linalg_status_no_convergence
            return
        end if
        x = log_b+log(alpha)*eye
        if (.not. finite_matrix(x)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        status = linalg_status_success
    end subroutine logm_core

    subroutine fractional_power_real_input(a, t, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(in) :: t !! real exponent; finite
        complex(dp), intent(out) :: x(:, :) !! principal A**t (n,n), possibly complex
        integer, intent(out) :: status !! linalg_status_* result

        call fractional_power_complex_input(cmplx(a, kind=dp), t, x, status, real_input=.true.)
    end subroutine fractional_power_real_input

    subroutine fractional_matrix_power_real(a, t, x, status)
        real(dp), intent(in) :: a(:, :) !! square real matrix (n,n)
        real(dp), intent(in) :: t !! real finite exponent
        real(dp), intent(out) :: x(:, :) !! real principal A**t, if real-valued
        integer, intent(out) :: status !! nonzero if result is not real or undefined
        complex(dp), allocatable :: z(:, :)
        real(dp) :: magnitude

        x = 0.0_dp
        allocate(z(size(x, 1), size(x, 2)))
        call fractional_matrix_power(a, t, z, status)
        if (status /= linalg_status_success) return
        magnitude = max(1.0_dp, maxval(abs(z)))
        if (maxval(abs(aimag(z))) > 1.0e-9_dp*magnitude) then
            status = linalg_status_no_convergence
            return
        end if
        x = real(z, dp)
    end subroutine fractional_matrix_power_real

    subroutine fractional_power_complex_input(a, t, x, status, real_input)
        complex(dp), intent(in) :: a(:, :) !! square complex matrix (n,n)
        real(dp), intent(in) :: t !! finite real exponent
        complex(dp), intent(out) :: x(:, :) !! principal matrix power A**t (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: real_input !! real input, for square-root branch
        complex(dp), allocatable :: logarithm(:, :), unused(:, :), zero(:, :), eye(:, :)
        real(dp) :: alpha
        integer :: n, k, code
        logical :: snap_real

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (n < 1 .or. n > logm_max_order .or. size(a, 2) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= n) return
        if (.not. ieee_is_finite(t) .or. .not. finite_matrix(a)) return
        allocate(logarithm(n, n), unused(n, n), zero(n, n), eye(n, n))
        eye = cmplx(0.0_dp, 0.0_dp, dp)
        do k = 1, n
            eye(k, k) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        if (t == 0.0_dp) then
            x = eye
            status = linalg_status_success
            return
        end if
        if (t == 1.0_dp) then
            x = a
            status = linalg_status_success
            return
        end if
        snap_real = .false.
        if (present(real_input)) snap_real = real_input
        if (t == 0.5_dp) then
            ! The existing sqrtm checks an absolute triangular pivot size;
            ! normalize here so extremely small nonzero matrices work.
            alpha = maxval(abs(a))
            if (alpha == 0.0_dp) then
                status = linalg_status_success
                return
            end if
            if (snap_real) then
                call sqrtm(real(a/alpha, dp), x, status)
            else
                call sqrtm(a/alpha, x, status)
            end if
            if (status /= linalg_status_success) return
            x = sqrt(alpha)*x
            if (.not. finite_matrix(x)) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                status = linalg_status_no_convergence
            end if
            return
        end if
        if (abs(t) <= real(huge(0), dp) .and. t == anint(t)) then
            call integral_matrix_power(a, int(t), x, status)
            return
        end if
        alpha = maxval(abs(a))
        if (alpha == 0.0_dp) then
            status = linalg_status_singular
            return
        end if
        call logm_core(a, logarithm, code, real_input=snap_real)
        if (code /= linalg_status_success) then
            status = code
            return
        end if
        if (.not. finite_matrix(t*logarithm)) then
            status = linalg_status_no_convergence
            return
        end if
        zero = cmplx(0.0_dp, 0.0_dp, dp)
        call expm_frechet(t*logarithm, zero, x, unused, status)
        if (status /= linalg_status_success) x = cmplx(0.0_dp, 0.0_dp, dp)
    end subroutine fractional_power_complex_input

    subroutine integral_matrix_power(a, exponent, x, status)
        complex(dp), intent(in) :: a(:, :) !! square matrix (n,n)
        integer, intent(in) :: exponent !! nonzero integer exponent
        complex(dp), intent(out) :: x(:, :) !! A raised to integer exponent
        integer, intent(out) :: status !! singular for negative power if not invertible
        complex(dp), allocatable :: b(:, :), eye(:, :), work(:, :)
        integer :: n, i, k, code

        n = size(a, 1)
        allocate(b(n, n), eye(n, n), work(n, n))
        eye = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            eye(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        x = eye
        if (exponent < 0) then
            call solve_dense(a, eye, b, code)
            if (code /= linalg_status_success) then
                status = code
                x = cmplx(0.0_dp, 0.0_dp, dp)
                return
            end if
            ! -huge(0) is safe; the most negative integer is excluded by int(t)
            k = -exponent
        else
            b = a
            k = exponent
        end if
        do while (k > 0)
            if (btest(k, 0)) then
                work = matmul(x, b)
                x = work
                if (.not. finite_matrix(x)) exit
            end if
            k = ishft(k, -1)
            if (k == 0) exit
            work = matmul(b, b)
            b = work
            if (.not. finite_matrix(b)) exit
        end do
        status = linalg_status_success
        if (.not. finite_matrix(x) .or. k /= 0) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
        end if
    end subroutine integral_matrix_power

    subroutine solve_dense(a, rhs, x, status)
        complex(dp), intent(in) :: a(:, :) !! square coefficient matrix (n,n)
        complex(dp), intent(in) :: rhs(:, :) !! right-hand sides (n,m)
        complex(dp), intent(out) :: x(:, :) !! solutions (n,m)
        integer, intent(out) :: status !! success, singular, or invalid shape
        complex(dp), allocatable :: b(:, :), r(:, :), row(:)
        complex(dp) :: factor
        real(dp) :: threshold, scale_a
        integer :: n, m, k, i, pivot

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 1)
        m = size(rhs, 2)
        if (n < 1 .or. size(a, 2) /= n .or. size(rhs, 1) /= n .or. &
                size(x, 1) /= n .or. size(x, 2) /= m) return
        if (.not. finite_matrix(a) .or. .not. finite_matrix(rhs)) return
        allocate(b(n, n), r(n, m), row(max(n,m)))
        b = a
        r = rhs
        scale_a = maxval(abs(a))
        threshold = max(32.0_dp*tiny(1.0_dp), 8.0_dp*epsilon(1.0_dp)*scale_a)
        status = linalg_status_singular
        do k = 1, n
            pivot = k-1+maxloc(abs(b(k:n, k)), dim=1)
            if (abs(b(pivot, k)) <= threshold) return
            if (pivot /= k) then
                row(1:n) = b(k, :)
                b(k, :) = b(pivot, :)
                b(pivot, :) = row(1:n)
                row(1:m) = r(k, :)
                r(k, :) = r(pivot, :)
                r(pivot, :) = row(1:m)
            end if
            do i = k+1, n
                factor = b(i, k)/b(k, k)
                b(i, k) = cmplx(0.0_dp, 0.0_dp, dp)
                b(i, k+1:n) = b(i, k+1:n)-factor*b(k, k+1:n)
                r(i, :) = r(i, :)-factor*r(k, :)
            end do
        end do
        x = r
        do k = n, 1, -1
            x(k, :) = (x(k, :)-matmul(b(k, k+1:n), x(k+1:n, :)))/b(k, k)
        end do
        if (.not. finite_matrix(x)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        status = linalg_status_success
    end subroutine solve_dense

    pure real(dp) function matrix_onenorm(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! matrix with finite entries
        value = maxval(sum(abs(a), dim=1))
    end function matrix_onenorm

    pure logical function finite_matrix(a) result(valid)
        complex(dp), intent(in) :: a(:, :) !! matrix to check for finite real and imaginary parts
        valid = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix
end module scifort_linalg_logarithm
