! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Copyright (c) 2010-2014 David Fong, Michael Saunders, Stanford University
!
! Matrix-free LSMR iteration used by scipy.optimize least-squares adapters.
! The algorithm follows scipy.sparse.linalg.lsmr (SciPy 1.17) and the
! original Fong-Saunders method. The callback interface is deliberately
! minimal so sparse/LinearOperator Jacobian adapters can reuse it later.

module scifort_optimize_lsmr
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    implicit none
    private

    abstract interface
        subroutine lsmr_matvec_function(x, y)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! input vector
            real(dp), intent(out) :: y(:) !! matrix-vector product
        end subroutine lsmr_matvec_function
    end interface

    public :: lsmr_solve

contains

    subroutine lsmr_solve(matvec, rmatvec, m, n, b, x, damp, atol, btol, conlim, &
            maxiter, istop, itn, normr, normar, norma, conda, normx)
        procedure(lsmr_matvec_function) :: matvec !! compute A*x
        procedure(lsmr_matvec_function) :: rmatvec !! compute transpose(A)*x
        integer, intent(in) :: m !! row count
        integer, intent(in) :: n !! column count
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! least-squares solution
        real(dp), intent(in), optional :: damp !! Tikhonov damping; default zero
        real(dp), intent(in), optional :: atol !! matrix backward-error tolerance
        real(dp), intent(in), optional :: btol !! rhs backward-error tolerance
        real(dp), intent(in), optional :: conlim !! condition-number limit
        integer, intent(in), optional :: maxiter !! iteration cap; default min(m,n)
        integer, intent(out), optional :: istop !! SciPy LSMR stopping code
        integer, intent(out), optional :: itn !! completed iterations
        real(dp), intent(out), optional :: normr !! estimated residual norm
        real(dp), intent(out), optional :: normar !! estimated normal-equation residual
        real(dp), intent(out), optional :: norma !! estimated Frobenius norm of Abar
        real(dp), intent(out), optional :: conda !! estimated condition number
        real(dp), intent(out), optional :: normx !! solution norm

        real(dp), allocatable :: h(:), hbar(:), u(:), v(:), work_m(:), work_n(:)
        real(dp) :: alpha, alphabar, alphahat, atol_use, beta, betaacute, betacheck
        real(dp) :: betad, betadd, betahat, btol_use, c, cbar, chat, conda_use
        real(dp) :: conlim_use, ctildeold, ctol, d, damp_use, inf, maxrbar, minrbar
        real(dp) :: norma_use, norma2, normar_use, normb, normr_use, normx_use
        real(dp) :: rho, rhobar, rhobarold, rhodold, rhoold, rhotemp, rhotildeold
        real(dp) :: rtol, s, sbar, shat, stildeold, t1, taud, tautildeold
        real(dp) :: test1, test2, test3, thetabar, thetanew, thetatilde
        real(dp) :: thetatildeold, zeta, zetabar, zetaold
        integer :: istop_use, itn_use, maxiter_use

        x = 0.0_dp
        istop_use = 0
        itn_use = 0
        normr_use = 0.0_dp
        normar_use = 0.0_dp
        norma_use = 0.0_dp
        conda_use = 1.0_dp
        normx_use = 0.0_dp
        if (m <= 0 .or. n <= 0 .or. size(b) /= m .or. size(x) /= n) then
            call copy_outputs()
            return
        end if

        damp_use = 0.0_dp
        atol_use = 1.0e-6_dp
        btol_use = 1.0e-6_dp
        conlim_use = 1.0e8_dp
        maxiter_use = min(m, n)
        if (present(damp)) damp_use = damp
        if (present(atol)) atol_use = atol
        if (present(btol)) btol_use = btol
        if (present(conlim)) conlim_use = conlim
        if (present(maxiter)) maxiter_use = maxiter
        if (damp_use < 0.0_dp .or. atol_use < 0.0_dp .or. btol_use < 0.0_dp .or. &
                conlim_use < 0.0_dp .or. maxiter_use < 0) then
            istop_use = 7
            call copy_outputs()
            return
        end if

        allocate(h(n), hbar(n), u(m), v(n), work_m(m), work_n(n))
        u = b
        normb = vector_norm2(u)
        beta = normb
        if (beta > 0.0_dp) then
            u = u / beta
            call rmatvec(u, v)
            alpha = vector_norm2(v)
        else
            v = 0.0_dp
            alpha = 0.0_dp
        end if
        if (alpha > 0.0_dp) v = v / alpha

        zetabar = alpha * beta
        alphabar = alpha
        rho = 1.0_dp
        rhobar = 1.0_dp
        cbar = 1.0_dp
        sbar = 0.0_dp
        h = v
        hbar = 0.0_dp

        betadd = beta
        betad = 0.0_dp
        rhodold = 1.0_dp
        tautildeold = 0.0_dp
        thetatilde = 0.0_dp
        zeta = 0.0_dp
        d = 0.0_dp

        norma2 = alpha * alpha
        maxrbar = 0.0_dp
        minrbar = huge(1.0_dp)
        norma_use = sqrt(norma2)
        conda_use = 1.0_dp
        normx_use = 0.0_dp
        ctol = 0.0_dp
        if (conlim_use > 0.0_dp) ctol = 1.0_dp / conlim_use
        normr_use = beta
        normar_use = alpha * beta
        if (same_real(normar_use, 0.0_dp) .or. same_real(normb, 0.0_dp) .or. maxiter_use == 0) then
            if (maxiter_use == 0 .and. .not. same_real(normar_use, 0.0_dp)) istop_use = 7
            call copy_outputs()
            return
        end if

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        do while (itn_use < maxiter_use)
            itn_use = itn_use + 1

            call matvec(v, work_m)
            u = work_m - alpha * u
            beta = vector_norm2(u)
            if (beta > 0.0_dp) then
                u = u / beta
                call rmatvec(u, work_n)
                v = work_n - beta * v
                alpha = vector_norm2(v)
                if (alpha > 0.0_dp) v = v / alpha
            end if

            call sym_ortho(alphabar, damp_use, chat, shat, alphahat)
            rhoold = rho
            call sym_ortho(alphahat, beta, c, s, rho)
            thetanew = s * alpha
            alphabar = c * alpha

            rhobarold = rhobar
            zetaold = zeta
            thetabar = sbar * rho
            rhotemp = cbar * rho
            call sym_ortho(cbar * rho, thetanew, cbar, sbar, rhobar)
            zeta = cbar * zetabar
            zetabar = -sbar * zetabar

            hbar = hbar * (-(thetabar * rho / (rhoold * rhobarold))) + h
            x = x + (zeta / (rho * rhobar)) * hbar
            h = v - (thetanew / rho) * h

            betaacute = chat * betadd
            betacheck = -shat * betadd
            betahat = c * betaacute
            betadd = -s * betaacute

            thetatildeold = thetatilde
            call sym_ortho(rhodold, thetabar, ctildeold, stildeold, rhotildeold)
            thetatilde = stildeold * rhobar
            rhodold = ctildeold * rhobar
            betad = -stildeold * betad + ctildeold * betahat

            tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold
            taud = (zeta - thetatilde * tautildeold) / rhodold
            d = d + betacheck * betacheck
            normr_use = sqrt(d + (betad - taud) ** 2 + betadd * betadd)

            norma2 = norma2 + beta * beta
            norma_use = sqrt(norma2)
            norma2 = norma2 + alpha * alpha

            maxrbar = max(maxrbar, rhobarold)
            if (itn_use > 1) minrbar = min(minrbar, rhobarold)
            conda_use = max(maxrbar, rhotemp) / min(minrbar, rhotemp)

            normar_use = abs(zetabar)
            normx_use = vector_norm2(x)
            test1 = normr_use / normb
            if (.not. same_real(norma_use * normr_use, 0.0_dp)) then
                test2 = normar_use / (norma_use * normr_use)
            else
                test2 = inf
            end if
            test3 = 1.0_dp / conda_use
            t1 = test1 / (1.0_dp + norma_use * normx_use / normb)
            rtol = btol_use + atol_use * norma_use * normx_use / normb

            istop_use = 0
            if (itn_use >= maxiter_use) istop_use = 7
            if (1.0_dp + test3 <= 1.0_dp) istop_use = 6
            if (1.0_dp + test2 <= 1.0_dp) istop_use = 5
            if (1.0_dp + t1 <= 1.0_dp) istop_use = 4
            if (test3 <= ctol) istop_use = 3
            if (test2 <= atol_use) istop_use = 2
            if (test1 <= rtol) istop_use = 1
            if (istop_use > 0) exit
        end do

        call copy_outputs()

    contains

        subroutine copy_outputs()
            if (present(istop)) istop = istop_use
            if (present(itn)) itn = itn_use
            if (present(normr)) normr = normr_use
            if (present(normar)) normar = normar_use
            if (present(norma)) norma = norma_use
            if (present(conda)) conda = conda_use
            if (present(normx)) normx = normx_use
        end subroutine copy_outputs

    end subroutine lsmr_solve

    subroutine sym_ortho(a, b, c, s, r)
        real(dp), intent(in) :: a !! first rotation element
        real(dp), intent(in) :: b !! second rotation element
        real(dp), intent(out) :: c !! cosine
        real(dp), intent(out) :: s !! sine
        real(dp), intent(out) :: r !! rotated magnitude

        real(dp) :: tau

        if (same_real(b, 0.0_dp)) then
            c = sign(1.0_dp, a)
            s = 0.0_dp
            r = abs(a)
        else if (same_real(a, 0.0_dp)) then
            c = 0.0_dp
            s = sign(1.0_dp, b)
            r = abs(b)
        else if (abs(b) > abs(a)) then
            tau = a / b
            s = sign(1.0_dp, b) / sqrt(1.0_dp + tau * tau)
            c = s * tau
            r = b / s
        else
            tau = b / a
            c = sign(1.0_dp, a) / sqrt(1.0_dp + tau * tau)
            s = c * tau
            r = a / c
        end if
    end subroutine sym_ortho

    pure logical function same_real(a, b) result(equal)
        real(dp), intent(in) :: a, b
        equal = a <= b .and. a >= b
    end function same_real

    pure real(dp) function vector_norm2(x) result(value)
        real(dp), intent(in) :: x(:) !! vector
        value = sqrt(max(0.0_dp, dot_product(x, x)))
    end function vector_norm2

end module scifort_optimize_lsmr
