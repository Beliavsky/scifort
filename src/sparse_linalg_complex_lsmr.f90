! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. and 2003 SciPy Developers
! Copyright (c) 2010-2014 David Fong and Michael Saunders, Stanford University
!
! Complex binary64 LSMR for canonical CSR/CSC matrices. The Golub-Kahan/
! orthogonal-rotation recurrences are adapted from SciFort's existing
! scifort_optimize_lsmr, itself derived from the SciPy LSMR implementation.
! The Hermitian adjoint A^H is required for complex least-squares problems.
! This module retains the BSD-3-Clause terms: see LICENSES/scipy-lsmr-BSD-3-Clause.txt.
module scifort_sparse_linalg_complex_lsmr
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr, csr_complex_hermitian
    implicit none
    private

    integer, parameter, public :: complex_lsmr_invalid_input = -1
    integer, parameter, public :: complex_lsmr_breakdown = -2
    public :: lsmr_complex

    interface lsmr_complex
        module procedure lsmr_complex_csr
        module procedure lsmr_complex_csc
    end interface lsmr_complex

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate complex value.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure real(dp) function complex_norm(z) result(value)
        complex(dp), intent(in) :: z(:) !! Complex input vector.
        real(dp) :: scale, ssq, component
        integer :: i, part
        scale = 0.0_dp
        ssq = 1.0_dp
        do i = 1,size(z)
            do part = 1,2
                if (part == 1) then
                    component = abs(real(z(i),dp))
                else
                    component = abs(aimag(z(i)))
                end if
                if (component == 0.0_dp) cycle
                if (scale < component) then
                    ssq = 1.0_dp+ssq*(scale/component)**2
                    scale = component
                else
                    ssq = ssq+(component/scale)**2
                end if
            end do
        end do
        value = scale*sqrt(ssq)
    end function complex_norm

    pure subroutine product(a,x,y)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        complex(dp), intent(in) :: x(:) !! Input, length a%n_cols.
        complex(dp), intent(out) :: y(:) !! Output, length a%n_rows.
        integer :: i, p
        y = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1,a%n_rows
            do p = a%indptr(i),a%indptr(i+1)-1
                y(i) = y(i)+a%data(p)*x(a%indices(p))
            end do
        end do
    end subroutine product

    subroutine sym_ortho(a,b,c,s,r)
        real(dp), intent(in) :: a !! First Givens entry.
        real(dp), intent(in) :: b !! Second Givens entry.
        real(dp), intent(out) :: c !! Rotation cosine.
        real(dp), intent(out) :: s !! Rotation sine.
        real(dp), intent(out) :: r !! Rotated length.
        real(dp) :: tau
        if (b == 0.0_dp) then
            c = sign(1.0_dp,a)
            s = 0.0_dp
            r = abs(a)
        else if (a == 0.0_dp) then
            c = 0.0_dp
            s = sign(1.0_dp,b)
            r = abs(b)
        else if (abs(b) > abs(a)) then
            tau = a/b
            s = sign(1.0_dp,b)/sqrt(1.0_dp+tau*tau)
            c = s*tau
            r = b/s
        else
            tau = b/a
            c = sign(1.0_dp,a)/sqrt(1.0_dp+tau*tau)
            s = c*tau
            r = a/c
        end if
    end subroutine sym_ortho

    subroutine lsmr_complex_csc(a,b,x,istop,itn,normr,normar,norma,conda,normx, &
            damp,atol,btol,conlim,maxiter,x0)
        type(csc_complex_matrix), intent(in) :: a !! Canonical complex CSC matrix (m,n).
        complex(dp), intent(in) :: b(:) !! Right-hand side, length m.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution, length n; empty on failure.
        integer, intent(out) :: istop !! SciPy-style 0..7, or negative failure.
        integer, intent(out), optional :: itn !! Completed iterations.
        real(dp), intent(out), optional :: normr !! True augmented residual norm.
        real(dp), intent(out), optional :: normar !! True normal-equation residual norm.
        real(dp), intent(out), optional :: norma !! Estimate of ||[A;damp*I]||.
        real(dp), intent(out), optional :: conda !! Estimated condition number.
        real(dp), intent(out), optional :: normx !! Norm of computed solution.
        real(dp), intent(in), optional :: damp !! Nonnegative regularization, default zero.
        real(dp), intent(in), optional :: atol !! Matrix backward-error tolerance, default 1e-6.
        real(dp), intent(in), optional :: btol !! Right-hand-side tolerance, default 1e-6.
        real(dp), intent(in), optional :: conlim !! Condition threshold, default 1e8.
        integer, intent(in), optional :: maxiter !! Iteration cap, default min(m,n).
        complex(dp), intent(in), optional :: x0(:) !! Initial iterate/regularization center.
        type(csr_complex_matrix) :: csr
        integer :: status

        istop = complex_lsmr_invalid_input
        allocate(x(0))
        if (present(itn)) itn = 0
        if (present(normr)) normr = huge(1.0_dp)
        if (present(normar)) normar = huge(1.0_dp)
        if (present(norma)) norma = 0.0_dp
        if (present(conda)) conda = 0.0_dp
        if (present(normx)) normx = 0.0_dp
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(x)
        call lsmr_complex_csr(csr,b,x,istop,itn,normr,normar,norma,conda,normx, &
            damp,atol,btol,conlim,maxiter,x0)
    end subroutine lsmr_complex_csc

    subroutine lsmr_complex_csr(a,b,x,istop,itn,normr,normar,norma,conda,normx, &
            damp,atol,btol,conlim,maxiter,x0)
        type(csr_complex_matrix), intent(in) :: a !! Canonical complex CSR matrix (m,n).
        complex(dp), intent(in) :: b(:) !! Right-hand side, length m.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution, length n; empty on failure.
        integer, intent(out) :: istop !! SciPy-style 0..7 or negative error.
        integer, intent(out), optional :: itn !! Completed bidiagonalization steps.
        real(dp), intent(out), optional :: normr !! True augmented residual norm.
        real(dp), intent(out), optional :: normar !! True normal-equation residual norm.
        real(dp), intent(out), optional :: norma !! Estimated norm of augmented A.
        real(dp), intent(out), optional :: conda !! Estimated condition number.
        real(dp), intent(out), optional :: normx !! True norm of returned x.
        real(dp), intent(in), optional :: damp !! Tikhonov damping, default zero.
        real(dp), intent(in), optional :: atol !! Backward-error tolerance, default 1e-6.
        real(dp), intent(in), optional :: btol !! Backward-error tolerance, default 1e-6.
        real(dp), intent(in), optional :: conlim !! Condition threshold, default 1e8 (zero disables).
        integer, intent(in), optional :: maxiter !! Iteration limit, default min(m,n).
        complex(dp), intent(in), optional :: x0(:) !! Initial iterate and Tikhonov center.

        type(csr_complex_matrix) :: ah
        complex(dp), allocatable :: h(:), hbar(:), u(:), v(:), wm(:), wn(:), base(:), r(:), grad(:)
        real(dp) :: alpha, alphabar, alphahat, beta, betaacute, betacheck
        real(dp) :: betad, betadd, betahat, c, cbar, chat, ctildeold, ctol, d
        real(dp) :: maxrbar, minrbar, normb, norma_use, norma2, normr_est, normar_est
        real(dp) :: normx_est, conda_use, rho, rhobar, rhobarold, rhodold, rhoold
        real(dp) :: rhotemp, rhotildeold, rtol, s, sbar, shat, stildeold, t1
        real(dp) :: taud, tautildeold, test1, test2, test3, thetabar, thetanew
        real(dp) :: thetatilde, thetatildeold, zeta, zetabar, zetaold, damp_use
        real(dp) :: atol_use, btol_use, conlim_use
        integer :: m, n, cap, steps, status

        istop = complex_lsmr_invalid_input
        if (present(itn)) itn = 0
        if (present(normr)) normr = huge(1.0_dp)
        if (present(normar)) normar = huge(1.0_dp)
        if (present(norma)) norma = 0.0_dp
        if (present(conda)) conda = 0.0_dp
        if (present(normx)) normx = 0.0_dp
        allocate(x(0))
        if (.not. valid_csr_complex(a)) return
        m = a%n_rows
        n = a%n_cols
        if (size(b) /= m .or. .not. all(finite_complex(b))) return
        damp_use = 0.0_dp
        atol_use = 1.0e-6_dp
        btol_use = 1.0e-6_dp
        conlim_use = 1.0e8_dp
        cap = min(m,n)
        if (present(damp)) damp_use = damp
        if (present(atol)) atol_use = atol
        if (present(btol)) btol_use = btol
        if (present(conlim)) conlim_use = conlim
        if (present(maxiter)) cap = maxiter
        if (.not. all(ieee_is_finite([damp_use,atol_use,btol_use,conlim_use]))) return
        if (minval([damp_use,atol_use,btol_use,conlim_use]) < 0.0_dp .or. cap < 0) return
        allocate(base(n))
        base = cmplx(0.0_dp,0.0_dp,dp)
        if (present(x0)) then
            if (size(x0) /= n) return
            if (.not. all(finite_complex(x0))) return
            base = x0
        end if
        ah = csr_complex_hermitian(a,status)
        if (status /= 0) return
        deallocate(x)
        allocate(x(n),h(n),hbar(n),u(m),v(n),wm(m),wn(n),r(m),grad(n))
        x = cmplx(0.0_dp,0.0_dp,dp)
        steps = 0
        call product(a,base,wm)
        if (.not. all(finite_complex(wm))) goto 800
        u = b-wm
        normb = complex_norm(u)
        beta = normb
        if (.not. ieee_is_finite(beta)) goto 800
        steps = 0
        normr_est = beta
        normar_est = 0.0_dp
        norma_use = 0.0_dp
        conda_use = 1.0_dp
        istop = 0
        if (beta == 0.0_dp .or. m == 0 .or. n == 0) goto 900
        u = u/beta
        call product(ah,u,v)
        if (.not. all(finite_complex(v))) goto 800
        alpha = complex_norm(v)
        if (.not. ieee_is_finite(alpha)) goto 800
        normar_est = alpha*beta
        if (alpha == 0.0_dp) goto 900
        v = v/alpha
        if (cap == 0) then
            istop = 7
            goto 900
        end if

        zetabar = alpha*beta
        alphabar = alpha
        rho = 1.0_dp
        rhobar = 1.0_dp
        cbar = 1.0_dp
        sbar = 0.0_dp
        h = v
        hbar = cmplx(0.0_dp,0.0_dp,dp)
        betadd = beta
        betad = 0.0_dp
        rhodold = 1.0_dp
        tautildeold = 0.0_dp
        thetatilde = 0.0_dp
        zeta = 0.0_dp
        d = 0.0_dp
        norma2 = alpha*alpha
        maxrbar = 0.0_dp
        minrbar = huge(1.0_dp)
        norma_use = sqrt(norma2)
        ctol = 0.0_dp
        if (conlim_use > 0.0_dp) ctol = 1.0_dp/conlim_use
        istop = 7
        do steps = 1,cap
            call product(a,v,wm)
            u = wm-alpha*u
            if (.not. all(finite_complex(u))) goto 800
            beta = complex_norm(u)
            if (.not. ieee_is_finite(beta)) goto 800
            if (beta > 0.0_dp) then
                u = u/beta
                call product(ah,u,wn)
                v = wn-beta*v
                if (.not. all(finite_complex(v))) goto 800
                alpha = complex_norm(v)
                if (.not. ieee_is_finite(alpha)) goto 800
                if (alpha > 0.0_dp) v = v/alpha
            else
                alpha = 0.0_dp
                v = cmplx(0.0_dp,0.0_dp,dp)
            end if

            call sym_ortho(alphabar,damp_use,chat,shat,alphahat)
            rhoold = rho
            call sym_ortho(alphahat,beta,c,s,rho)
            thetanew = s*alpha
            alphabar = c*alpha
            rhobarold = rhobar
            zetaold = zeta
            thetabar = sbar*rho
            rhotemp = cbar*rho
            call sym_ortho(cbar*rho,thetanew,cbar,sbar,rhobar)
            zeta = cbar*zetabar
            zetabar = -sbar*zetabar
            if (rhoold == 0.0_dp .or. rhobarold == 0.0_dp .or. &
                rho == 0.0_dp .or. rhobar == 0.0_dp) goto 800
            hbar = hbar*(-thetabar*rho/(rhoold*rhobarold))+h
            x = x+(zeta/(rho*rhobar))*hbar
            h = v-(thetanew/rho)*h
            if (.not. all(finite_complex(x)) .or. .not. all(finite_complex(h))) goto 800

            betaacute = chat*betadd
            betacheck = -shat*betadd
            betahat = c*betaacute
            betadd = -s*betaacute
            thetatildeold = thetatilde
            call sym_ortho(rhodold,thetabar,ctildeold,stildeold,rhotildeold)
            thetatilde = stildeold*rhobar
            rhodold = ctildeold*rhobar
            betad = -stildeold*betad+ctildeold*betahat
            if (rhotildeold == 0.0_dp .or. rhodold == 0.0_dp) goto 800
            tautildeold = (zetaold-thetatildeold*tautildeold)/rhotildeold
            taud = (zeta-thetatilde*tautildeold)/rhodold
            d = d+betacheck*betacheck
            normr_est = sqrt(d+(betad-taud)**2+betadd*betadd)
            norma2 = norma2+beta*beta
            norma_use = sqrt(norma2)
            norma2 = norma2+alpha*alpha
            maxrbar = max(maxrbar,rhobarold)
            if (steps > 1) minrbar = min(minrbar,rhobarold)
            if (min(minrbar,abs(rhotemp)) > 0.0_dp) then
                conda_use = max(maxrbar,abs(rhotemp))/min(minrbar,abs(rhotemp))
            else
                conda_use = huge(1.0_dp)
            end if
            normar_est = abs(zetabar)
            normx_est = complex_norm(x)
            if (.not. all(ieee_is_finite([normr_est,norma_use,conda_use,normar_est,normx_est]))) goto 800
            test1 = normr_est/normb
            test2 = huge(1.0_dp)
            if (norma_use*normr_est > 0.0_dp) test2 = normar_est/(norma_use*normr_est)
            test3 = 1.0_dp/conda_use
            t1 = test1/(1.0_dp+norma_use*normx_est/normb)
            rtol = btol_use+atol_use*norma_use*normx_est/normb
            istop = 0
            if (steps >= cap) istop = 7
            if (1.0_dp+test3 <= 1.0_dp) istop = 6
            if (1.0_dp+test2 <= 1.0_dp) istop = 5
            if (1.0_dp+t1 <= 1.0_dp) istop = 4
            if (test3 <= ctol) istop = 3
            if (test2 <= atol_use) istop = 2
            if (test1 <= rtol) istop = 1
            if (istop > 0) exit
        end do
        if (steps > cap) steps = cap
900     continue
        x = base+x
        if (.not. all(finite_complex(x))) goto 800
        call product(a,x,wm)
        r = b-wm
        if (.not. all(finite_complex(r))) goto 800
        if (present(normr)) normr = hypot(complex_norm(r),damp_use*complex_norm(x-base))
        if (present(normar)) then
            call product(ah,r,grad)
            grad = grad-damp_use**2*(x-base)
            if (.not. all(finite_complex(grad))) goto 800
            normar = complex_norm(grad)
        end if
        if (present(itn)) itn = steps
        if (present(norma)) norma = norma_use
        if (present(conda)) conda = conda_use
        if (present(normx)) normx = complex_norm(x)
        return
800     continue
        istop = complex_lsmr_breakdown
        if (allocated(x)) deallocate(x)
        allocate(x(0))
        if (present(itn)) itn = steps
    end subroutine lsmr_complex_csr

end module scifort_sparse_linalg_complex_lsmr
