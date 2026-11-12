! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! Matrix-free real Krylov methods for the nonlinear-root Jacobian solves.
! BiCGSTAB and CGS recurrences follow SciPy 1.17.0's BSD-3-Clause
! scipy.sparse.linalg._isolve.iterative implementation; TFQMR follows
! scipy.sparse.linalg._isolve.tfqmr; MINRES follows the Paige-Saunders
! recurrence in scipy.sparse.linalg._isolve.minres (also BSD-3-Clause).
! MINRES requires a symmetric Jacobian and a positive-definite preconditioner.
!
module scifort_optimize_krylov_iterative
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    abstract interface
        subroutine krylov_operator(vec, product, ok)
            import dp
            implicit none
            real(dp), intent(in) :: vec(:) !! input vector
            real(dp), intent(out) :: product(:) !! operator product
            logical, intent(out) :: ok !! whether application succeeded
        end subroutine krylov_operator
    end interface

    public :: krylov_operator
    public :: krylov_iterative_solve

contains

    subroutine krylov_iterative_solve(method, matvec, preconditioner, rhs, &
            rtol, maxiter, solution, ok, verify_residual, right_preconditioned)
        character(len=*), intent(in) :: method !! bicgstab, cgs, minres, or tfqmr
        procedure(krylov_operator) :: matvec !! matrix-free Jacobian action
        procedure(krylov_operator) :: preconditioner !! approximate inverse, identity if absent
        real(dp), intent(in) :: rhs(:) !! linear-system right-hand side
        real(dp), intent(in) :: rtol !! relative residual stopping tolerance
        integer, intent(in) :: maxiter !! maximum inner iterations
        real(dp), intent(out) :: solution(:) !! approximate solution
        logical, intent(out) :: ok !! false on invalid input, operator failure, or breakdown
        logical, intent(in), optional :: verify_residual !! require true TFQMR residual convergence
        logical, intent(in), optional :: right_preconditioned !! use the consistent A*M TFQMR recurrence

        logical :: verify, right

        solution = 0.0_dp
        ok = .false.
        if (size(rhs) <= 0 .or. size(solution) /= size(rhs)) return
        if (maxiter <= 0 .or. .not. ieee_is_finite(rtol)) return
        if (rtol < 0.0_dp .or. .not. all(ieee_is_finite(rhs))) return
        verify = .false.
        if (present(verify_residual)) verify = verify_residual
        ! The opt-in true-residual criterion applies only to TFQMR. Reject
        ! unsupported methods instead of silently ignoring the request.
        if (verify .and. trim(method) /= 'tfqmr') return
        right = .false.
        if (present(right_preconditioned)) right = right_preconditioned
        if (right .and. trim(method) /= 'tfqmr') return
        select case (trim(method))
        case ('bicgstab')
            call solve_bicgstab(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        case ('cgs')
            call solve_cgs(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        case ('tfqmr')
            call solve_tfqmr(matvec, preconditioner, rhs, rtol, maxiter, solution, ok, verify, right)
        case ('minres')
            call solve_minres(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        end select
    end subroutine krylov_iterative_solve

    subroutine apply_checked(operator, vec, product, ok)
        procedure(krylov_operator) :: operator !! linear operator to apply
        real(dp), intent(in) :: vec(:) !! input vector
        real(dp), intent(out) :: product(:) !! output vector
        logical, intent(out) :: ok !! true only for a finite successful product

        ok = .false.
        if (.not. all(ieee_is_finite(vec))) return
        call operator(vec, product, ok)
        if (ok) ok = all(ieee_is_finite(product))
    end subroutine apply_checked

    subroutine solve_bicgstab(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(krylov_operator) :: matvec !! linear operator
        procedure(krylov_operator) :: preconditioner !! inverse preconditioner
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! relative residual tolerance
        integer, intent(in) :: maxiter !! maximum BiCGSTAB iterations
        real(dp), intent(out) :: solution(:) !! approximate solution
        logical, intent(out) :: ok !! false on breakdown or operator failure

        real(dp) :: r(size(rhs)), rtilde(size(rhs)), p(size(rhs)), v(size(rhs))
        real(dp) :: phat(size(rhs)), s(size(rhs)), shat(size(rhs)), t(size(rhs))
        real(dp) :: alpha, beta, omega, rho, rho_prev, rv, tt, atol
        integer :: iteration
        logical :: valid

        solution = 0.0_dp
        ok = .false.
        r = rhs
        rtilde = r
        p = 0.0_dp
        v = 0.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp
        rho_prev = 1.0_dp
        atol = min(1.0_dp, rtol) * norm2(rhs)
        do iteration = 1, maxiter
            if (norm2(r) <= atol ) exit
            rho = dot_product(rtilde, r)
            if (.not. ieee_is_finite(rho) .or. abs(rho) < epsilon(1.0_dp)**2) return
            if (iteration > 1) then
                if (abs(omega) <= tiny(1.0_dp)) return
                beta = (rho / rho_prev) * (alpha / omega)
                if (.not. ieee_is_finite(beta)) return
                p = r + beta * (p - omega * v)
            else
                p = r
            end if
            call apply_checked(preconditioner, p, phat, valid)
            if (.not. valid) return
            call apply_checked(matvec, phat, v, valid)
            if (.not. valid) return
            rv = dot_product(rtilde, v)
            if (abs(rv) <= tiny(1.0_dp) .or. .not. ieee_is_finite(rv)) then
                ok = .false.
                return
            end if
            alpha = rho / rv
            if (.not. ieee_is_finite(alpha)) then
                ok = .false.
                return
            end if
            s = r - alpha * v
            if (norm2(s) <= atol ) then
                solution = solution + alpha * phat
                ok = all(ieee_is_finite(solution))
                return
            end if
            call apply_checked(preconditioner, s, shat, valid)
            if (.not. valid) return
            call apply_checked(matvec, shat, t, valid)
            if (.not. valid) return
            tt = dot_product(t, t)
            if (tt <= tiny(1.0_dp) .or. .not. ieee_is_finite(tt)) then
                ok = .false.
                return
            end if
            omega = dot_product(t, s) / tt
            if (.not. ieee_is_finite(omega)) then
                ok = .false.
                return
            end if
            solution = solution + alpha * phat + omega * shat
            r = s - omega * t
            if (.not. all(ieee_is_finite(solution)) .or. .not. all(ieee_is_finite(r))) then
                ok = .false.
                return
            end if
            rho_prev = rho
        end do
        ! Like scipy.optimize.KrylovJacobian, accept finite incomplete inner
        ! iterates: the outer Newton line search controls global convergence.
        ok = all(ieee_is_finite(solution))
    end subroutine solve_bicgstab

    subroutine solve_cgs(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(krylov_operator) :: matvec !! linear operator
        procedure(krylov_operator) :: preconditioner !! inverse preconditioner
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! relative residual tolerance
        integer, intent(in) :: maxiter !! maximum CGS iterations
        real(dp), intent(out) :: solution(:) !! approximate solution
        logical, intent(out) :: ok !! false on breakdown or operator failure

        real(dp) :: r(size(rhs)), rtilde(size(rhs)), p(size(rhs)), u(size(rhs))
        real(dp) :: q(size(rhs)), phat(size(rhs)), v(size(rhs)), uhat(size(rhs))
        real(dp) :: trial(size(rhs)), ax(size(rhs))
        real(dp) :: alpha, beta, rho, rho_prev, rv, atol
        integer :: iteration
        logical :: valid

        solution = 0.0_dp
        ok = .false.
        r = rhs
        rtilde = r
        p = 0.0_dp
        u = 0.0_dp
        q = 0.0_dp
        rho_prev = 1.0_dp
        atol = min(1.0_dp, rtol) * norm2(rhs)
        do iteration = 1, maxiter
            if (norm2(r) <= atol ) exit
            rho = dot_product(rtilde, r)
            if (.not. ieee_is_finite(rho) .or. abs(rho) < epsilon(1.0_dp)**2) return
            if (iteration > 1) then
                beta = rho / rho_prev
                if (.not. ieee_is_finite(beta)) return
                u = r + beta * q
                p = u + beta * (q + beta * p)
            else
                u = r
                p = r
            end if
            call apply_checked(preconditioner, p, phat, valid)
            if (.not. valid) return
            call apply_checked(matvec, phat, v, valid)
            if (.not. valid) return
            rv = dot_product(rtilde, v)
            if (abs(rv) <= tiny(1.0_dp) .or. .not. ieee_is_finite(rv)) then
                ok = .false.
                return
            end if
            alpha = rho / rv
            if (.not. ieee_is_finite(alpha)) then
                ok = .false.
                return
            end if
            q = u - alpha * v
            trial = u + q
            call apply_checked(preconditioner, trial, uhat, valid)
            if (.not. valid) return
            solution = solution + alpha * uhat
            if (.not. all(ieee_is_finite(solution))) then
                ok = .false.
                return
            end if
            ! SciPy recomputes b - A*x rather than accumulating r - alpha*A*uhat.
            call apply_checked(matvec, solution, ax, valid)
            if (.not. valid) return
            r = rhs - ax
            rho_prev = rho
        end do
        ok = all(ieee_is_finite(solution))
    end subroutine solve_cgs

    subroutine apply_tfqmr_operator(matvec, preconditioner, vec, product, work, ok, right_preconditioned)
        procedure(krylov_operator) :: matvec !! matrix-free Jacobian action
        procedure(krylov_operator) :: preconditioner !! inverse preconditioner
        real(dp), intent(in) :: vec(:) !! Krylov recurrence operand
        real(dp), intent(out) :: product(:) !! preconditioned operator product
        real(dp), intent(out) :: work(:) !! scratch vector
        logical, intent(out) :: ok !! both operator applications succeeded
        logical, intent(in) :: right_preconditioned !! use A*M rather than SciPy's M*A

        ! The TFQMR iterate is always updated as x = x + eta*M*d. Thus
        ! right preconditioning must use the composed operator A*M on d.
        ! SciPy 1.17 uses M*A in its recurrence, retained as the default
        ! for compatibility, but M*A generally does not solve A*x = b with
        ! the M*d update for a nontrivial preconditioner.
        if (right_preconditioned) then
            call apply_checked(preconditioner, vec, work, ok)
            if (.not. ok) return
            call apply_checked(matvec, work, product, ok)
        else
            call apply_checked(matvec, vec, work, ok)
            if (.not. ok) return
            call apply_checked(preconditioner, work, product, ok)
        end if
    end subroutine apply_tfqmr_operator

    subroutine solve_tfqmr(matvec, preconditioner, rhs, rtol, maxiter, solution, ok, &
            verify_residual, right_preconditioned)
        procedure(krylov_operator) :: matvec !! linear operator
        procedure(krylov_operator) :: preconditioner !! inverse preconditioner
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! relative residual tolerance
        integer, intent(in) :: maxiter !! maximum TFQMR half-iterations
        real(dp), intent(out) :: solution(:) !! approximate solution
        logical, intent(out) :: ok !! false on breakdown or operator failure
        logical, intent(in) :: verify_residual !! check original residual before acceptance
        logical, intent(in) :: right_preconditioned !! use consistent right preconditioning

        real(dp) :: rstar(size(rhs)), u(size(rhs)), w(size(rhs)), v(size(rhs))
        real(dp) :: uhat(size(rhs)), unext(size(rhs)), d(size(rhs))
        real(dp) :: z(size(rhs)), av(size(rhs))
        real(dp) :: alpha, beta, c, eta, rho, rho_last, tau, theta, vtrstar, atol
        integer :: iteration
        logical :: valid
        logical :: checked_last_residual

        solution = 0.0_dp
        ok = .false.
        if (maxval(abs(rhs)) <= 0.0_dp) then
            ok = .true.
            return
        end if
        rstar = rhs
        u = rhs
        w = rhs
        d = 0.0_dp
        call apply_tfqmr_operator(matvec, preconditioner, rhs, v, av, valid, &
            right_preconditioned)
        if (.not. valid) return
        uhat = v
        rho = dot_product(rhs, rhs)
        rho_last = rho
        tau = sqrt(rho)
        theta = 0.0_dp
        eta = 0.0_dp
        alpha = 0.0_dp
        atol = min(1.0_dp, rtol) * tau
        checked_last_residual = .false.
        do iteration = 1, maxiter
            checked_last_residual = .false.
            if (mod(iteration, 2) == 1) then
                vtrstar = dot_product(rstar, v)
                if (abs(vtrstar) <= tiny(1.0_dp) .or. .not. ieee_is_finite(vtrstar)) then
                    ok = .false.
                    return
                end if
                alpha = rho / vtrstar
                if (abs(alpha) <= tiny(1.0_dp) .or. .not. ieee_is_finite(alpha)) then
                    ok = .false.
                    return
                end if
                unext = u - alpha * v
            end if
            w = w - alpha * uhat
            d = u + (theta * theta * eta / alpha) * d
            if (tau <= tiny(1.0_dp)) exit
            theta = norm2(w) / tau
            c = 1.0_dp / hypot(1.0_dp, theta)
            tau = tau * theta * c
            eta = c * c * alpha
            call apply_checked(preconditioner, d, z, valid)
            if (.not. valid) return
            solution = solution + eta * z
            if (.not. all(ieee_is_finite(solution))) then
                ok = .false.
                return
            end if
            ! Match SciPy's estimator-based stopping by default. With
            ! verification, accept only when the original (unpreconditioned)
            ! residual norm ||rhs - A*solution|| meets the requested tolerance.
            ! This adds one matrix-free action per candidate convergence test.
            if (tau * sqrt(real(iteration, dp)) <= atol .or. tau <= tiny(1.0_dp)) then
                if (.not. verify_residual) exit
                call apply_checked(matvec, solution, av, valid)
                if (.not. valid) return
                checked_last_residual = .true.
                if (norm2(rhs - av) <= atol) then
                    ok = .true.
                    return
                end if
                ! A zero estimator cannot support further reliable recurrences.
                if (tau <= tiny(1.0_dp)) return
            end if
            if (mod(iteration, 2) == 0) then
                rho = dot_product(rstar, w)
                if (abs(rho_last) <= tiny(1.0_dp) .or. .not. ieee_is_finite(rho)) then
                    ok = .false.
                    return
                end if
                beta = rho / rho_last
                if (.not. ieee_is_finite(beta)) then
                    ok = .false.
                    return
                end if
                u = w + beta * u
                v = beta * uhat + beta * beta * v
                call apply_tfqmr_operator(matvec, preconditioner, u, uhat, av, &
                    valid, right_preconditioned)
                if (.not. valid) return
                v = v + uhat
            else
                call apply_tfqmr_operator(matvec, preconditioner, unext, uhat, av, &
                    valid, right_preconditioned)
                if (.not. valid) return
                u = unext
                rho_last = rho
            end if
        end do
        if (verify_residual) then
            ! A finite truncated iterate is acceptable in SciPy-compatible
            ! mode, but never claim verified convergence without checking it.
            if (checked_last_residual) return
            call apply_checked(matvec, solution, av, valid)
            if (.not. valid) return
            ok = norm2(rhs - av) <= atol
        else
            ok = all(ieee_is_finite(solution))
        end if
    end subroutine solve_tfqmr

    subroutine solve_minres(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(krylov_operator) :: matvec !! symmetric linear operator
        procedure(krylov_operator) :: preconditioner !! positive-definite inverse
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(in) :: rtol !! relative MINRES tolerance
        integer, intent(in) :: maxiter !! maximum Lanczos iterations
        real(dp), intent(out) :: solution(:) !! approximate solution
        logical, intent(out) :: ok !! false if preconditioner indefinite or operator fails

        real(dp) :: r1(size(rhs)), r2(size(rhs)), y(size(rhs)), v(size(rhs))
        real(dp) :: av(size(rhs)), w(size(rhs)), w1(size(rhs)), w2(size(rhs))
        real(dp) :: alpha, anorm, beta, beta1, beta2, cs, dbar, delta, epsln
        real(dp) :: gamma, gbar, oldb, oldeps, phi, phibar, root, sn
        real(dp) :: test1, test2, tnorm2, ynorm
        integer :: iteration
        logical :: valid

        solution = 0.0_dp
        ok = .false.
        r1 = rhs
        r2 = rhs
        if (maxval(abs(rhs)) <= 0.0_dp) then
            ok = .true.
            return
        end if
        call apply_checked(preconditioner, r1, y, valid)
        if (.not. valid) return
        beta1 = dot_product(r1, y)
        if (.not. ieee_is_finite(beta1) .or. beta1 <= 0.0_dp) then
            ok = .false.
            return
        end if
        beta1 = sqrt(beta1)
        beta = beta1
        oldb = 0.0_dp
        dbar = 0.0_dp
        epsln = 0.0_dp
        tnorm2 = 0.0_dp
        phibar = beta1
        cs = -1.0_dp
        sn = 0.0_dp
        w = 0.0_dp
        w2 = 0.0_dp

        do iteration = 1, maxiter
            if (beta <= 0.0_dp) exit
            v = y / beta
            call apply_checked(matvec, v, av, valid)
            if (.not. valid) return
            if (iteration >= 2) av = av - (beta / oldb) * r1
            alpha = dot_product(v, av)
            av = av - (alpha / beta) * r2
            r1 = r2
            r2 = av
            call apply_checked(preconditioner, r2, y, valid)
            if (.not. valid) return
            oldb = beta
            beta2 = dot_product(r2, y)
            if (.not. ieee_is_finite(beta2) .or. beta2 < 0.0_dp) then
                ok = .false.
                return
            end if
            beta = sqrt(beta2)
            tnorm2 = tnorm2 + alpha * alpha + oldb * oldb + beta * beta
            oldeps = epsln
            delta = cs * dbar + sn * alpha
            gbar = sn * dbar - cs * alpha
            epsln = sn * beta
            dbar = -cs * beta
            root = hypot(gbar, dbar)
            gamma = max(hypot(gbar, beta), epsilon(1.0_dp))
            cs = gbar / gamma
            sn = beta / gamma
            phi = cs * phibar
            phibar = sn * phibar
            w1 = w2
            w2 = w
            w = (v - oldeps * w1 - delta * w2) / gamma
            solution = solution + phi * w
            if (.not. all(ieee_is_finite(solution))) then
                ok = .false.
                return
            end if
            anorm = sqrt(tnorm2)
            ynorm = norm2(solution)
            test1 = huge(1.0_dp)
            test2 = huge(1.0_dp)
            if (anorm > 0.0_dp .and. ynorm > 0.0_dp) &
                test1 = abs(phibar) / (anorm * ynorm)
            if (anorm > 0.0_dp) test2 = root / anorm
            if (test1 <= rtol .or. test2 <= rtol) exit
            if (abs(phibar) <= epsilon(1.0_dp) * beta1) exit
        end do
        ok = all(ieee_is_finite(solution))
    end subroutine solve_minres

end module scifort_optimize_krylov_iterative
