! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dependency-free complex binary64 LSQR for rectangular CSR/CSC matrices.
! Uses Golub-Kahan bidiagonalization with the conjugate-transpose operator.
! The algorithm is written independently from the standard published LSQR
! recurrences (Paige & Saunders, 1982); no third-party code is incorporated.
module scifort_sparse_linalg_complex_lsqr
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr, csr_complex_hermitian
    implicit none
    private

    integer, parameter, public :: complex_lsqr_invalid_input = -1
    integer, parameter, public :: complex_lsqr_breakdown = -2
    public :: lsqr_complex

    interface lsqr_complex
        module procedure lsqr_complex_csr, lsqr_complex_csc
    end interface lsqr_complex

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate scalar.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    ! Scale the real and imaginary parts separately to avoid squaring overflow.
    pure real(dp) function complex_norm(z) result(value)
        complex(dp), intent(in) :: z(:) !! Vector to measure.
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

    ! Internal product for already validated CSR input, without allocation.
    pure subroutine product(a, x, y)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR operator.
        complex(dp), intent(in) :: x(:) !! Input vector (n_cols).
        complex(dp), intent(out) :: y(:) !! Output vector (n_rows).
        integer :: i, p

        y = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1,a%n_rows
            do p = a%indptr(i),a%indptr(i+1)-1
                y(i) = y(i)+a%data(p)*x(a%indices(p))
            end do
        end do
    end subroutine product

    subroutine lsqr_complex_csc(a, b, x, istop, itn, r1norm, r2norm, &
            damp, atol, btol, conlim, iter_lim, x0, anorm, acond, arnorm, xnorm, var, calc_var)
        type(csc_complex_matrix), intent(in) :: a !! Complex CSC matrix (m,n).
        complex(dp), intent(in) :: b(:) !! Complex right-hand side (m).
        complex(dp), allocatable, intent(out) :: x(:) !! Approximate solution (n), empty on invalid input.
        integer, intent(out) :: istop !! LSQR stop code 0,1,2,3,7, or negative failure.
        integer, intent(out), optional :: itn !! Number of bidiagonalization iterations.
        real(dp), intent(out), optional :: r1norm !! True norm of b-A*x.
        real(dp), intent(out), optional :: r2norm !! Norm of augmented residual.
        real(dp), intent(in), optional :: damp !! Nonnegative damping strength, default 0.
        real(dp), intent(in), optional :: atol !! Relative matrix precision, default 1e-6.
        real(dp), intent(in), optional :: btol !! Relative right-hand-side precision, default 1e-6.
        real(dp), intent(in), optional :: conlim !! Condition threshold, default 1e8; 0 disables.
        integer, intent(in), optional :: iter_lim !! Maximum steps, default 2*min(m,n).
        complex(dp), intent(in), optional :: x0(:) !! Initial solution and Tikhonov center (n).
        real(dp), intent(out), optional :: anorm !! Estimated Frobenius norm of augmented operator.
        real(dp), intent(out), optional :: acond !! Estimated condition number.
        real(dp), intent(out), optional :: arnorm !! True augmented normal-equation residual.
        real(dp), intent(out), optional :: xnorm !! Norm of returned solution.
        real(dp), allocatable, intent(out), optional :: var(:) !! Real variance-diagonal estimate.
        logical, intent(in), optional :: calc_var !! Accumulate variance estimate, default false.
        type(csr_complex_matrix) :: csr
        integer :: status

        istop = complex_lsqr_invalid_input
        allocate(x(0))
        if (present(itn)) itn = 0
        if (present(r1norm)) r1norm = huge(1.0_dp)
        if (present(r2norm)) r2norm = huge(1.0_dp)
        if (present(anorm)) anorm = 0.0_dp
        if (present(acond)) acond = 0.0_dp
        if (present(arnorm)) arnorm = huge(1.0_dp)
        if (present(xnorm)) xnorm = 0.0_dp
        if (present(var)) allocate(var(0))
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(x)
        if (present(var)) deallocate(var)
        call lsqr_complex_csr(csr,b,x,istop,itn,r1norm,r2norm, &
            damp,atol,btol,conlim,iter_lim,x0,anorm,acond,arnorm,xnorm,var,calc_var)
    end subroutine lsqr_complex_csc

    subroutine lsqr_complex_csr(a, b, x, istop, itn, r1norm, r2norm, &
            damp, atol, btol, conlim, iter_lim, x0, anorm, acond, arnorm, xnorm, var, calc_var)
        type(csr_complex_matrix), intent(in) :: a !! Complex CSR matrix (m,n).
        complex(dp), intent(in) :: b(:) !! Complex right-hand side (m).
        complex(dp), allocatable, intent(out) :: x(:) !! Approximate solution (n), empty on invalid input.
        integer, intent(out) :: istop !! LSQR stop code 0,1,2,3,7, or negative failure.
        integer, intent(out), optional :: itn !! Bidiagonalization iteration count.
        real(dp), intent(out), optional :: r1norm !! True norm of b-A*x.
        real(dp), intent(out), optional :: r2norm !! Norm of augmented residual.
        real(dp), intent(in), optional :: damp !! Nonnegative damping, default 0.
        real(dp), intent(in), optional :: atol !! Relative matrix precision, default 1e-6.
        real(dp), intent(in), optional :: btol !! Relative RHS precision, default 1e-6.
        real(dp), intent(in), optional :: conlim !! Condition limit, default 1e8; 0 disables.
        integer, intent(in), optional :: iter_lim !! Maximum steps, default 2*min(m,n).
        complex(dp), intent(in), optional :: x0(:) !! Initial guess/Tikhonov center.
        real(dp), intent(out), optional :: anorm !! Estimated augmented matrix norm.
        real(dp), intent(out), optional :: acond !! Estimated augmented condition number.
        real(dp), intent(out), optional :: arnorm !! True normal-equation residual norm.
        real(dp), intent(out), optional :: xnorm !! Euclidean norm of solution.
        real(dp), allocatable, intent(out), optional :: var(:) !! Real variance estimate diagonal.
        logical, intent(in), optional :: calc_var !! Enable variance accumulation, default false.
        type(csr_complex_matrix) :: ah
        complex(dp), allocatable :: base(:), dx(:), u(:), v(:), w(:), work_m(:), &
            work_n(:), gradient(:)
        real(dp) :: damping, tol_a, tol_b, cond_limit, alpha, beta, beta0
        real(dp) :: rhobar, phibar, rhobar1, cs1, sn1, psi, rho, cs, sn
        real(dp) :: theta, phi, tau, ddnorm, residual2, est_rnorm, est_arnorm
        real(dp) :: norm_a, condition, norm_x, test_bound
        integer :: m, n, cap, steps, status
        logical :: estimate_var

        istop = complex_lsqr_invalid_input
        if (present(itn)) itn = 0
        if (present(r1norm)) r1norm = huge(1.0_dp)
        if (present(r2norm)) r2norm = huge(1.0_dp)
        if (present(anorm)) anorm = 0.0_dp
        if (present(acond)) acond = 0.0_dp
        if (present(arnorm)) arnorm = huge(1.0_dp)
        if (present(xnorm)) xnorm = 0.0_dp
        if (present(var)) allocate(var(0))
        allocate(x(0))
        if (.not. valid_csr_complex(a)) return
        m = a%n_rows
        n = a%n_cols
        if (size(b) /= m .or. .not. all(finite_complex(b))) return
        damping = 0.0_dp
        tol_a = 1.0e-6_dp
        tol_b = 1.0e-6_dp
        cond_limit = 1.0e8_dp
        if (present(damp)) damping = damp
        if (present(atol)) tol_a = atol
        if (present(btol)) tol_b = btol
        if (present(conlim)) cond_limit = conlim
        if (.not. all(ieee_is_finite([damping,tol_a,tol_b,cond_limit]))) return
        if (minval([damping,tol_a,tol_b,cond_limit]) < 0.0_dp) return
        if (min(m,n) > huge(cap)/2) then
            cap = huge(cap)
        else
            cap = 2*min(m,n)
        end if
        if (present(iter_lim)) cap = iter_lim
        if (cap < 0) return
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
        allocate(x(n),dx(n),u(m),v(n),w(n),work_m(m),work_n(n),gradient(n))
        x = base
        dx = cmplx(0.0_dp,0.0_dp,dp)
        estimate_var = .false.
        if (present(calc_var)) estimate_var = calc_var
        if (present(var)) then
            deallocate(var)
            allocate(var(n))
            var = 0.0_dp
        end if
        steps = 0
        call product(a,base,work_m)
        if (.not. all(finite_complex(work_m))) goto 800
        u = b-work_m
        beta = complex_norm(u)
        beta0 = beta
        norm_a = 0.0_dp
        condition = 0.0_dp
        norm_x = 0.0_dp
        steps = 0
        istop = 0
        if (.not. ieee_is_finite(beta)) goto 800
        if (beta == 0.0_dp) goto 900
        if (m == 0 .or. n == 0) then
            istop = 2
            goto 900
        end if
        u = u/beta
        call product(ah,u,work_n)
        if (.not. all(finite_complex(work_n))) goto 800
        alpha = complex_norm(work_n)
        if (.not. ieee_is_finite(alpha)) goto 800
        if (alpha == 0.0_dp) then
            istop = 2
            goto 900
        end if
        v = work_n/alpha
        w = v
        rhobar = alpha
        phibar = beta
        ddnorm = 0.0_dp
        residual2 = 0.0_dp
        istop = 7
        do steps = 1,cap
            call product(a,v,work_m)
            u = work_m-alpha*u
            if (.not. all(finite_complex(u))) goto 800
            beta = complex_norm(u)
            if (.not. ieee_is_finite(beta)) goto 800
            if (beta /= 0.0_dp) u = u/beta
            if (beta /= 0.0_dp) then
                call product(ah,u,work_n)
            else
                work_n = cmplx(0.0_dp,0.0_dp,dp)
            end if
            v = work_n-beta*v
            if (.not. all(finite_complex(v))) goto 800
            alpha = complex_norm(v)
            if (.not. ieee_is_finite(alpha)) goto 800
            if (alpha /= 0.0_dp) v = v/alpha
            norm_a = hypot(norm_a,hypot(alpha,hypot(beta,damping)))

            rhobar1 = hypot(rhobar,damping)
            if (rhobar1 == 0.0_dp) goto 800
            cs1 = rhobar/rhobar1
            sn1 = damping/rhobar1
            psi = sn1*phibar
            phibar = cs1*phibar
            rho = hypot(rhobar1,beta)
            if (rho == 0.0_dp) goto 800
            cs = rhobar1/rho
            sn = beta/rho
            theta = sn*alpha
            rhobar = -cs*alpha
            phi = cs*phibar
            phibar = sn*phibar
            tau = sn*phi
            dx = dx+(phi/rho)*w
            if (present(var)) then
                if (estimate_var) var = var+abs(w/rho)**2
            end if
            ddnorm = ddnorm+complex_norm(w/rho)**2
            w = v-(theta/rho)*w
            residual2 = residual2+psi**2
            est_rnorm = hypot(phibar,sqrt(residual2))
            est_arnorm = abs(alpha*tau)
            norm_x = complex_norm(dx)
            condition = norm_a*sqrt(ddnorm)
            if (.not. all(finite_complex(dx)) .or. &
                .not. all(ieee_is_finite([norm_x,norm_a,condition,est_rnorm,est_arnorm]))) goto 800
            test_bound = tol_b*beta0+tol_a*norm_a*norm_x
            if (est_rnorm <= test_bound) then
                istop = 1
                exit
            end if
            if (est_arnorm <= tol_a*norm_a*est_rnorm) then
                istop = 2
                exit
            end if
            if (cond_limit > 0.0_dp .and. condition >= cond_limit) then
                istop = 3
                exit
            end if
        end do
        if (steps > cap) steps = cap
        x = base+dx
        if (.not. all(finite_complex(x))) goto 800
900     continue
        call product(a,x,work_m)
        if (.not. all(finite_complex(work_m))) goto 800
        work_m = b-work_m
        if (.not. all(finite_complex(work_m))) goto 800
        if (present(r1norm)) r1norm = complex_norm(work_m)
        if (present(r2norm)) r2norm = hypot(complex_norm(work_m),damping*complex_norm(x-base))
        if (present(itn)) itn = steps
        if (present(anorm)) anorm = norm_a
        if (present(acond)) acond = condition
        if (present(xnorm)) xnorm = complex_norm(x)
        if (present(arnorm)) then
            call product(ah,work_m,gradient)
            gradient = gradient-damping**2*(x-base)
            if (.not. all(finite_complex(gradient))) goto 800
            arnorm = complex_norm(gradient)
        end if
        return
800     continue
        istop = complex_lsqr_breakdown
        if (allocated(x)) deallocate(x)
        allocate(x(0))
        if (present(var)) then
            if (allocated(var)) deallocate(var)
            allocate(var(0))
        end if
        if (present(itn)) itn = steps
    end subroutine lsqr_complex_csr

end module scifort_sparse_linalg_complex_lsqr
