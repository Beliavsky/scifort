! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense matrix trigonometric and hyperbolic functions.
! Independently implemented paired Taylor series with argument halving and
! repeated double-angle identities. No LAPACK/BLAS or SciPy source is used.
!
module scifort_linalg_trigonometric
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_matrix_functions, only : matrix_function_max_order, dense_inverse
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter :: max_terms = 40
    integer, parameter :: max_scalings = 60
    integer, parameter :: operation_sine = 1
    integer, parameter :: operation_cosine = 2
    integer, parameter :: operation_tangent = 3
    real(dp), parameter :: taylor_tolerance = 4.0_dp*epsilon(1.0_dp)

    public :: sinm, cosm, tanm, sinhm, coshm, tanhm

    interface sinm
        module procedure sinm_real
        module procedure sinm_complex
    end interface sinm
    interface cosm
        module procedure cosm_real
        module procedure cosm_complex
    end interface cosm
    interface tanm
        module procedure tanm_real
        module procedure tanm_complex
    end interface tanm
    interface sinhm
        module procedure sinhm_real
        module procedure sinhm_complex
    end interface sinhm
    interface coshm
        module procedure coshm_real
        module procedure coshm_complex
    end interface coshm
    interface tanhm
        module procedure tanhm_real
        module procedure tanhm_complex
    end interface tanhm

contains

    subroutine sinm_real(a, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! sine of the input matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_real(a, x, status, .false., operation_sine)
    end subroutine sinm_real

    subroutine sinm_complex(a, x, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! sine of the input matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_complex(a, x, status, .false., operation_sine)
    end subroutine sinm_complex

    subroutine cosm_real(a, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! cosine of the input matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_real(a, x, status, .false., operation_cosine)
    end subroutine cosm_real

    subroutine cosm_complex(a, x, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! cosine of the input matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_complex(a, x, status, .false., operation_cosine)
    end subroutine cosm_complex

    subroutine tanm_real(a, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! tangent of the input matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result; singular at cosine poles
        call trig_real(a, x, status, .false., operation_tangent)
    end subroutine tanm_real

    subroutine tanm_complex(a, x, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! tangent of the input matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result; singular at cosine poles
        call trig_complex(a, x, status, .false., operation_tangent)
    end subroutine tanm_complex

    subroutine sinhm_real(a, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! hyperbolic sine of the matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_real(a, x, status, .true., operation_sine)
    end subroutine sinhm_real

    subroutine sinhm_complex(a, x, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! hyperbolic sine of the matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_complex(a, x, status, .true., operation_sine)
    end subroutine sinhm_complex

    subroutine coshm_real(a, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! hyperbolic cosine of the matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_real(a, x, status, .true., operation_cosine)
    end subroutine coshm_real

    subroutine coshm_complex(a, x, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! hyperbolic cosine of the matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        call trig_complex(a, x, status, .true., operation_cosine)
    end subroutine coshm_complex

    subroutine tanhm_real(a, x, status)
        real(dp), intent(in) :: a(:, :) !! real square matrix (n,n)
        real(dp), intent(out) :: x(:, :) !! hyperbolic tangent of the matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result; singular at cosh poles
        call trig_real(a, x, status, .true., operation_tangent)
    end subroutine tanhm_real

    subroutine tanhm_complex(a, x, status)
        complex(dp), intent(in) :: a(:, :) !! complex square matrix (n,n)
        complex(dp), intent(out) :: x(:, :) !! hyperbolic tangent of the matrix (n,n)
        integer, intent(out) :: status !! linalg_status_* result; singular at cosh poles
        call trig_complex(a, x, status, .true., operation_tangent)
    end subroutine tanhm_complex

    subroutine trig_real(a, x, status, hyperbolic, operation)
        real(dp), intent(in) :: a(:, :) !! real square argument (n,n)
        real(dp), intent(out) :: x(:, :) !! real matrix function (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in) :: hyperbolic !! true for hyperbolic rather than circular functions
        integer, intent(in) :: operation !! sine=1, cosine=2, tangent=3
        complex(dp), allocatable :: xc(:, :)
        real(dp) :: magnitude

        x = 0.0_dp
        allocate(xc(size(x, 1), size(x, 2)))
        call trig_complex(cmplx(a, kind=dp), xc, status, hyperbolic, operation)
        if (status /= linalg_status_success) return
        magnitude = max(1.0_dp, maxval(abs(xc)))
        if (maxval(abs(aimag(xc))) > 1.0e-12_dp*magnitude) then
            status = linalg_status_no_convergence
            return
        end if
        x = real(xc, dp)
    end subroutine trig_real

    subroutine trig_complex(a, x, status, hyperbolic, operation)
        complex(dp), intent(in) :: a(:, :) !! complex square argument (n,n)
        complex(dp), intent(out) :: x(:, :) !! matrix function (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in) :: hyperbolic !! true for hyperbolic rather than circular functions
        integer, intent(in) :: operation !! sine=1, cosine=2, tangent=3
        complex(dp), allocatable :: s(:, :), c(:, :), inverse(:, :)
        integer :: n, code
        real(dp) :: error

        x = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (n < 1 .or. n > matrix_function_max_order) return
        if (size(a, 2) /= n .or. size(x, 1) /= n .or. size(x, 2) /= n) return
        if (operation < operation_sine .or. operation > operation_tangent) return
        if (.not. finite_matrix(a)) then
            status = linalg_status_no_convergence
            return
        end if
        allocate(s(n, n), c(n, n))
        call scaled_pair(a, s, c, status, hyperbolic)
        if (status /= linalg_status_success) return
        select case (operation)
        case (operation_sine)
            x = s
        case (operation_cosine)
            x = c
        case (operation_tangent)
            allocate(inverse(n, n))
            call dense_inverse(c, inverse, code)
            if (code /= linalg_status_success) then
                status = code
                return
            end if
            x = matmul(s, inverse)
            if (.not. finite_matrix(x)) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                status = linalg_status_no_convergence
                return
            end if
            ! Reject unreliable tangent solves near poles or with large cancellation.
            error = maxval(abs(matmul(x, c)-s))/ &
                max(tiny(1.0_dp), maxval(abs(s)), maxval(abs(x))*maxval(abs(c)))
            if (.not. ieee_is_finite(error) .or. error > 1.0e-8_dp) then
                x = cmplx(0.0_dp, 0.0_dp, dp)
                status = linalg_status_no_convergence
                return
            end if
        end select
    end subroutine trig_complex

    subroutine scaled_pair(a, s, c, status, hyperbolic)
        complex(dp), intent(in) :: a(:, :) !! finite square matrix, order <= 40
        complex(dp), intent(out) :: s(:, :) !! sine or hyperbolic sine
        complex(dp), intent(out) :: c(:, :) !! cosine or hyperbolic cosine
        integer, intent(out) :: status !! success, overflow or nonconvergence
        logical, intent(in) :: hyperbolic !! selects sign in Taylor and doubling formulas
        complex(dp), allocatable :: b(:, :), b2(:, :), st(:, :), ct(:, :), sn(:, :), cn(:, :)
        real(dp) :: anorm, factor, sign_term, denominator_s, denominator_c
        integer :: n, steps, k, j, i
        logical :: converged

        n = size(a, 1)
        s = cmplx(0.0_dp, 0.0_dp, dp)
        c = cmplx(0.0_dp, 0.0_dp, dp)
        status = linalg_status_no_convergence
        anorm = maxval(sum(abs(a), dim=1))
        if (.not. ieee_is_finite(anorm)) return
        steps = 0
        if (anorm > 0.5_dp) then
            steps = ceiling(log(anorm/0.5_dp)/log(2.0_dp))
        end if
        if (steps > max_scalings) return
        factor = scale(1.0_dp, -steps)
        allocate(b(n, n), b2(n, n), st(n, n), ct(n, n), sn(n, n), cn(n, n))
        b = factor*a
        b2 = matmul(b, b)
        st = b
        ct = cmplx(0.0_dp, 0.0_dp, dp)
        s = b
        do i = 1, n
            ct(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
            c(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        sign_term = -1.0_dp
        if (hyperbolic) sign_term = 1.0_dp
        converged = .false.
        do k = 1, max_terms
            denominator_c = real(2*k*(2*k-1), dp)
            denominator_s = real(2*k*(2*k+1), dp)
            ct = sign_term*matmul(ct, b2)/denominator_c
            st = sign_term*matmul(st, b2)/denominator_s
            c = c+ct
            s = s+st
            if (.not. finite_matrix(c) .or. .not. finite_matrix(s)) return
            if (maxval(abs(ct)) <= taylor_tolerance*max(1.0_dp, maxval(abs(c))) .and. &
                    maxval(abs(st)) <= taylor_tolerance*max(1.0_dp, maxval(abs(s)))) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) return
        do j = 1, steps
            sn = 2.0_dp*matmul(s, c)
            if (hyperbolic) then
                cn = matmul(c, c)+matmul(s, s)
            else
                cn = matmul(c, c)-matmul(s, s)
            end if
            s = sn
            c = cn
            if (.not. finite_matrix(c) .or. .not. finite_matrix(s)) return
        end do
        status = linalg_status_success
    end subroutine scaled_pair

    pure logical function finite_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! complex array being checked
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix

end module scifort_linalg_trigonometric
