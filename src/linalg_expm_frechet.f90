! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Frechet derivative of exp(A) along E, using differentiated Taylor
! polynomials with scaling-and-squaring (no LAPACK/BLAS dependencies).
! For very large ||A|| the exponential can overflow in binary64.
module scifort_linalg_expm_frechet
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    integer, parameter, public :: expm_frechet_max_order = 48
    integer, parameter :: default_terms = 120
    real(dp), parameter :: default_tol = 8.0_dp*epsilon(1.0_dp)
    public :: expm_frechet

    interface expm_frechet
        module procedure expm_frechet_real
        module procedure expm_frechet_complex
    end interface expm_frechet

contains

    subroutine expm_frechet_real(a, e, exp_a, frechet, status, tol, max_terms, terms, scalings)
        real(dp), intent(in) :: a(:, :) !! real square matrix A (n,n)
        real(dp), intent(in) :: e(:, :) !! real perturbation E (n,n)
        real(dp), intent(out) :: exp_a(:, :) !! exp(A) (n,n)
        real(dp), intent(out) :: frechet(:, :) !! derivative L_exp(A,E) (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! relative Taylor tolerance, default 8 eps
        integer, intent(in), optional :: max_terms !! maximum Taylor terms, default 120
        integer, intent(out), optional :: terms !! Taylor iterations used
        integer, intent(out), optional :: scalings !! number of squarings
        complex(dp), allocatable :: cexp(:, :), cder(:, :)

        exp_a = 0.0_dp
        frechet = 0.0_dp
        allocate(cexp(size(exp_a, 1), size(exp_a, 2)), cder(size(frechet, 1), size(frechet, 2)))
        call expm_frechet_complex(cmplx(a, kind=dp), cmplx(e, kind=dp), cexp, cder, &
            status, tol, max_terms, terms, scalings)
        if (status /= linalg_status_success) return
        exp_a = real(cexp, dp)
        frechet = real(cder, dp)
    end subroutine expm_frechet_real

    subroutine expm_frechet_complex(a, e, exp_a, frechet, status, tol, max_terms, terms, scalings)
        complex(dp), intent(in) :: a(:, :) !! square matrix A (n,n)
        complex(dp), intent(in) :: e(:, :) !! direction E (n,n)
        complex(dp), intent(out) :: exp_a(:, :) !! exp(A) (n,n)
        complex(dp), intent(out) :: frechet(:, :) !! derivative L_exp(A,E) (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        real(dp), intent(in), optional :: tol !! relative Taylor tolerance, default 8 eps
        integer, intent(in), optional :: max_terms !! maximum Taylor terms, default 120
        integer, intent(out), optional :: terms !! Taylor iterations used
        integer, intent(out), optional :: scalings !! number of squarings
        complex(dp), allocatable :: as(:, :), es(:, :), t(:, :), d(:, :), tn(:, :), dn(:, :), &
            next_x(:, :), next_d(:, :)
        real(dp) :: tolerance, anorm, enorm, factor, xnorm, dnorm
        integer :: n, i, k, steps, term_limit
        logical :: converged

        exp_a = cmplx(0.0_dp, 0.0_dp, dp)
        frechet = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(terms)) terms = 0
        if (present(scalings)) scalings = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (n < 1 .or. n > expm_frechet_max_order) return
        if (size(a, 2) /= n .or. size(e, 1) /= n .or. size(e, 2) /= n) return
        if (size(exp_a, 1) /= n .or. size(exp_a, 2) /= n .or. &
                size(frechet, 1) /= n .or. size(frechet, 2) /= n) return
        tolerance = default_tol
        if (present(tol)) tolerance = tol
        term_limit = default_terms
        if (present(max_terms)) term_limit = max_terms
        if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp .or. &
                term_limit < 1) return
        if (.not. finite_complex(a) .or. .not. finite_complex(e)) return

        allocate(as(n, n), es(n, n), t(n, n), d(n, n), tn(n, n), dn(n, n), &
            next_x(n, n), next_d(n, n))
        anorm = maxval(sum(abs(a), dim=1))
        if (.not. ieee_is_finite(anorm)) then
            status = linalg_status_no_convergence
            return
        end if
        steps = 0
        if (anorm > 0.5_dp) then
            steps = max(0, ceiling(log(anorm/0.5_dp)/log(2.0_dp)))
        end if
        ! The squaring count is bounded for work and overflow safety.
        if (steps > 60) then
            status = linalg_status_no_convergence
            return
        end if
        if (present(scalings)) scalings = steps
        factor = scale(1.0_dp, -steps)
        as = a*factor
        es = e*factor
        enorm = fro_norm(es)
        t = cmplx(0.0_dp, 0.0_dp, dp)
        d = cmplx(0.0_dp, 0.0_dp, dp)
        exp_a = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            exp_a(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
            t(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        converged = .false.
        do k = 1, term_limit
            ! Differentiate T_k=A^k/k! in direction E.
            dn = (matmul(es, t)+matmul(as, d))/real(k, dp)
            tn = matmul(as, t)/real(k, dp)
            exp_a = exp_a+tn
            frechet = frechet+dn
            if (present(terms)) terms = k
            if (.not. finite_complex(exp_a) .or. .not. finite_complex(frechet)) exit
            xnorm = fro_norm(exp_a)
            dnorm = fro_norm(frechet)
            if (fro_norm(tn) <= tolerance*max(1.0_dp, xnorm) .and. &
                    fro_norm(dn) <= tolerance*max(tiny(1.0_dp), enorm, dnorm)) then
                converged = .true.
                exit
            end if
            t = tn
            d = dn
        end do
        if (.not. converged) then
            exp_a = cmplx(0.0_dp, 0.0_dp, dp)
            frechet = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
            return
        end if
        do i = 1, steps
            next_d = matmul(exp_a, frechet)+matmul(frechet, exp_a)
            next_x = matmul(exp_a, exp_a)
            exp_a = next_x
            frechet = next_d
            if (.not. finite_complex(exp_a) .or. .not. finite_complex(frechet)) then
                exp_a = cmplx(0.0_dp, 0.0_dp, dp)
                frechet = cmplx(0.0_dp, 0.0_dp, dp)
                status = linalg_status_no_convergence
                return
            end if
        end do
        status = linalg_status_success
    end subroutine expm_frechet_complex

    pure real(dp) function fro_norm(a) result(value)
        complex(dp), intent(in) :: a(:, :) !! complex matrix
        value = sqrt(sum(abs(a)**2))
    end function fro_norm

    pure logical function finite_complex(a)
        complex(dp), intent(in) :: a(:, :) !! complex matrix
        finite_complex = all(ieee_is_finite(real(a, dp))) .and. &
            all(ieee_is_finite(aimag(a)))
    end function finite_complex
end module scifort_linalg_expm_frechet
