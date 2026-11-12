! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense continuous/discrete stabilizing algebraic Riccati equations.
! CARE: matrix-sign iteration on the Hamiltonian, with Newton correction.
! DARE: structure-preserving doubling, with Newton correction.
! Self-contained pivoted solves; no BLAS/LAPACK requirement.  Intended for
! small and moderate, well-separated stabilizable control problems.
!
module scifort_linalg_riccati
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence, linalg_status_not_positive_definite
    use scifort_linalg_matrix_equations, only : solve_sylvester, solve_discrete_lyapunov, &
        matrix_equation_max_unknowns
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: solve_continuous_are
    public :: solve_discrete_are
    integer, parameter, public :: riccati_max_states = 20
    integer, parameter :: default_max_iter = 100
    real(dp), parameter :: default_tol = 1.0e-11_dp

    interface solve_continuous_are
        module procedure continuous_are_real
        module procedure continuous_are_complex
    end interface
    interface solve_discrete_are
        module procedure discrete_are_real
        module procedure discrete_are_complex
    end interface

contains

    subroutine continuous_are_real(a, b, q, r, x, status, e, s, max_iter, tol, iterations, residual)
        real(dp), intent(in) :: a(:, :) !! system dynamics (states,states)
        real(dp), intent(in) :: b(:, :) !! control matrix (states,controls)
        real(dp), intent(in) :: q(:, :) !! Hermitian positive-semidefinite state cost (states,states)
        real(dp), intent(in) :: r(:, :) !! positive-definite control cost (controls,controls)
        real(dp), intent(out) :: x(:, :) !! stabilizing symmetric solution (states,states)
        integer, intent(out) :: status !! zero success, otherwise linalg status code
        real(dp), intent(in), optional :: e(:, :) !! optional nonsingular descriptor matrix
        real(dp), intent(in), optional :: s(:, :) !! optional state-control cross cost
        integer, intent(in), optional :: max_iter !! maximum sign iterations, default 100
        real(dp), intent(in), optional :: tol !! relative tolerance, default 1e-11
        integer, intent(out), optional :: iterations !! matrix-sign iterations performed
        real(dp), intent(out), optional :: residual !! normalized residual of original equation
        complex(dp), allocatable :: ec(:, :), sc(:, :), xc(:, :)
        integer :: n, m

        n = size(a, 1)
        m = size(b, 2)
        allocate(ec(n, n), sc(n, m), xc(size(x, 1), size(x, 2)))
        ec = identity(n)
        sc = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(e)) then
            if (size(e, 1) /= n .or. size(e, 2) /= n) then
                x = 0.0_dp
                status = linalg_status_invalid_shape
                if (present(iterations)) iterations = 0
                if (present(residual)) residual = huge(1.0_dp)
                return
            end if
            ec = cmplx(e, kind=dp)
        end if
        if (present(s)) then
            if (size(s, 1) /= n .or. size(s, 2) /= m) then
                x = 0.0_dp
                status = linalg_status_invalid_shape
                if (present(iterations)) iterations = 0
                if (present(residual)) residual = huge(1.0_dp)
                return
            end if
            sc = cmplx(s, kind=dp)
        end if
        call continuous_are_complex(cmplx(a, kind=dp), cmplx(b, kind=dp), &
            cmplx(q, kind=dp), cmplx(r, kind=dp), xc, status, e=ec, s=sc, &
            max_iter=max_iter, tol=tol, iterations=iterations, residual=residual)
        x = real(xc, dp)
    end subroutine continuous_are_real

    subroutine discrete_are_real(a, b, q, r, x, status, e, s, max_iter, tol, iterations, residual)
        real(dp), intent(in) :: a(:, :) !! system dynamics (states,states)
        real(dp), intent(in) :: b(:, :) !! control matrix (states,controls)
        real(dp), intent(in) :: q(:, :) !! Hermitian positive-semidefinite state cost (states,states)
        real(dp), intent(in) :: r(:, :) !! positive-definite control cost (controls,controls)
        real(dp), intent(out) :: x(:, :) !! stabilizing symmetric solution (states,states)
        integer, intent(out) :: status !! zero success, otherwise linalg status code
        real(dp), intent(in), optional :: e(:, :) !! optional nonsingular descriptor matrix
        real(dp), intent(in), optional :: s(:, :) !! optional state-control cross cost
        integer, intent(in), optional :: max_iter !! maximum doubling iterations, default 100
        real(dp), intent(in), optional :: tol !! relative tolerance, default 1e-11
        integer, intent(out), optional :: iterations !! doubling iterations performed
        real(dp), intent(out), optional :: residual !! normalized residual of original equation
        complex(dp), allocatable :: ec(:, :), sc(:, :), xc(:, :)
        integer :: n, m

        n = size(a, 1)
        m = size(b, 2)
        allocate(ec(n, n), sc(n, m), xc(size(x, 1), size(x, 2)))
        ec = identity(n)
        sc = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(e)) then
            if (size(e, 1) /= n .or. size(e, 2) /= n) then
                x = 0.0_dp
                status = linalg_status_invalid_shape
                if (present(iterations)) iterations = 0
                if (present(residual)) residual = huge(1.0_dp)
                return
            end if
            ec = cmplx(e, kind=dp)
        end if
        if (present(s)) then
            if (size(s, 1) /= n .or. size(s, 2) /= m) then
                x = 0.0_dp
                status = linalg_status_invalid_shape
                if (present(iterations)) iterations = 0
                if (present(residual)) residual = huge(1.0_dp)
                return
            end if
            sc = cmplx(s, kind=dp)
        end if
        call discrete_are_complex(cmplx(a, kind=dp), cmplx(b, kind=dp), &
            cmplx(q, kind=dp), cmplx(r, kind=dp), xc, status, e=ec, s=sc, &
            max_iter=max_iter, tol=tol, iterations=iterations, residual=residual)
        x = real(xc, dp)
    end subroutine discrete_are_real

    subroutine continuous_are_complex(a, b, q, r, x, status, e, s, max_iter, tol, iterations, residual)
        complex(dp), intent(in) :: a(:, :) !! complex system matrix (n,n)
        complex(dp), intent(in) :: b(:, :) !! control matrix (n,m)
        complex(dp), intent(in) :: q(:, :) !! Hermitian state cost (n,n)
        complex(dp), intent(in) :: r(:, :) !! Hermitian positive-definite control cost (m,m)
        complex(dp), intent(out) :: x(:, :) !! stabilizing Hermitian Riccati solution (n,n)
        integer, intent(out) :: status !! zero success, otherwise linalg status
        complex(dp), intent(in), optional :: e(:, :) !! descriptor matrix (n,n), identity by default
        complex(dp), intent(in), optional :: s(:, :) !! cross cost (n,m), zero by default
        integer, intent(in), optional :: max_iter !! maximum sign iterations, default 100
        real(dp), intent(in), optional :: tol !! relative tolerance, default 1e-11
        integer, intent(out), optional :: iterations !! sign iterations performed
        real(dp), intent(out), optional :: residual !! normalized CARE residual
        complex(dp), allocatable :: ab(:, :), bb(:, :), qb(:, :), g(:, :), einv(:, :)
        complex(dp), allocatable :: z(:, :), zi(:, :), zn(:, :), rhs(:, :), y(:, :)
        complex(dp), allocatable :: f(:, :), acl(:, :), delta(:, :)
        real(dp) :: epsilon_tol, mu, err, stopping
        integer :: n, iter, maxit, code, refine
        logical :: converged

        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(iterations)) iterations = 0
        if (present(residual)) residual = huge(1.0_dp)
        call prepare_equation(a, b, q, r, x, e, s, ab, bb, qb, g, einv, status)
        if (status /= linalg_status_success) return
        n = size(a, 1)
        epsilon_tol = default_tol
        if (present(tol)) epsilon_tol = tol
        maxit = default_max_iter
        if (present(max_iter)) maxit = max_iter
        if (.not. ieee_is_finite(epsilon_tol) .or. epsilon_tol <= 0.0_dp .or. maxit < 1) then
            status = linalg_status_invalid_shape
            return
        end if
        allocate(z(2*n, 2*n), zi(2*n, 2*n), zn(2*n, 2*n), rhs(n,n), y(n,n))
        z = cmplx(0.0_dp, 0.0_dp, dp)
        z(1:n,1:n) = ab
        z(1:n,n+1:2*n) = -g
        z(n+1:2*n,1:n) = -qb
        z(n+1:2*n,n+1:2*n) = -adjoint(ab)
        converged = .false.
        do iter = 1, maxit
            call solve_dense(z, identity(2*n), zi, code)
            if (code /= linalg_status_success) then
                status = linalg_status_no_convergence
                return
            end if
            mu = sqrt(fnorm(zi) / max(fnorm(z), tiny(1.0_dp)))
            zn = 0.5_dp * (mu*z + zi/mu)
            if (.not. finite_matrix(zn)) exit
            err = fnorm(zn-z) / max(1.0_dp, fnorm(zn))
            z = zn
            if (err <= max(epsilon_tol, 128.0_dp*epsilon(1.0_dp)) .and. &
                    fnorm(matmul(z,z)-identity(2*n)) <= 32.0_dp * epsilon_tol * real(2*n,dp)) then
                converged = .true.
                exit
            end if
        end do
        if (present(iterations)) iterations = min(iter, maxit)
        if (.not. converged) then
            status = linalg_status_no_convergence
            return
        end if
        rhs = -(z(n+1:2*n,1:n))
        call solve_dense(z(n+1:2*n,n+1:2*n)+identity(n), rhs, y, code)
        if (code /= linalg_status_success) then
            ! A stabilizing graph can exist even when (I+S22) is singular,
            ! e.g. scalar Q=0 with unstable dynamics.  Use the top block.
            call solve_dense(z(1:n,n+1:2*n), -(identity(n)+z(1:n,1:n)), y, code)
        end if
        if (code /= linalg_status_success) then
            ! As a last resort, use both subspace equations together.
            block
                complex(dp), allocatable :: top(:, :), bottom(:, :), first(:, :)
                top = z(1:n,n+1:2*n)
                bottom = z(n+1:2*n,n+1:2*n)+identity(n)
                first = -(identity(n)+z(1:n,1:n))
                call solve_dense(matmul(adjoint(top),top)+matmul(adjoint(bottom),bottom), &
                    matmul(adjoint(top),first)+matmul(adjoint(bottom),rhs), y, code)
            end block
        end if
        if (code /= linalg_status_success) then
            status = linalg_status_no_convergence
            return
        end if
        y = hermitian_part(y)
        ! Newton corrections: F(Y)=A^H Y+YA-Y G Y+Q.
        allocate(f(n,n), acl(n,n), delta(n,n))
        do refine = 1, 6
            f = matmul(adjoint(ab),y)+matmul(y,ab)-matmul(matmul(y,g),y)+qb
            stopping = fnorm(f) / max(1.0_dp, fnorm(qb), fnorm(matmul(y,ab)), &
                fnorm(matmul(matmul(y,g),y)))
            if (stopping <= max(2.0e-12_dp, 4.0_dp*epsilon_tol)) exit
            acl = ab-matmul(g,y)
            call solve_sylvester(adjoint(acl), acl, -f, delta, code)
            if (code /= linalg_status_success) exit
            if (.not. finite_matrix(delta)) exit
            y = hermitian_part(y+delta)
        end do
        x = matmul(matmul(adjoint(einv),y),einv)
        x = hermitian_part(x)
        call verify_residual(a,b,q,r,x,.true.,e,s,err,code)
        if (present(residual)) residual = err
        acl = ab-matmul(g,y)
        ! The converged Hamiltonian sign selects the stable invariant graph.
        ! A power-norm stability screen is not reliable for nonnormal feedback.
        if (code /= linalg_status_success .or. err > max(2.0e-8_dp, 128.0_dp*epsilon_tol)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
        else
            status = linalg_status_success
        end if
    end subroutine continuous_are_complex

    subroutine discrete_are_complex(a, b, q, r, x, status, e, s, max_iter, tol, iterations, residual)
        complex(dp), intent(in) :: a(:, :) !! complex system matrix (n,n)
        complex(dp), intent(in) :: b(:, :) !! control matrix (n,m)
        complex(dp), intent(in) :: q(:, :) !! Hermitian state cost (n,n)
        complex(dp), intent(in) :: r(:, :) !! Hermitian positive-definite control cost (m,m)
        complex(dp), intent(out) :: x(:, :) !! stabilizing Hermitian Riccati solution (n,n)
        integer, intent(out) :: status !! zero success, otherwise linalg status
        complex(dp), intent(in), optional :: e(:, :) !! descriptor matrix (n,n), identity by default
        complex(dp), intent(in), optional :: s(:, :) !! cross cost (n,m), zero by default
        integer, intent(in), optional :: max_iter !! maximum doubling iterations, default 100
        real(dp), intent(in), optional :: tol !! relative tolerance, default 1e-11
        integer, intent(out), optional :: iterations !! doubling iterations performed
        real(dp), intent(out), optional :: residual !! normalized DARE residual
        complex(dp), allocatable :: ab(:, :), bb(:, :), qb(:, :), g(:, :), einv(:, :)
        complex(dp), allocatable :: h(:, :), ak(:, :), an(:, :), gn(:, :), hn(:, :)
        complex(dp), allocatable :: t(:, :), u(:, :), y(:, :), delta(:, :), f(:, :)
        complex(dp), allocatable :: rb(:, :), k(:, :), acl(:, :)
        real(dp) :: epsilon_tol, err
        integer :: n, m, iter, maxit, code, refine
        logical :: converged

        x = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(iterations)) iterations = 0
        if (present(residual)) residual = huge(1.0_dp)
        call prepare_equation(a,b,q,r,x,e,s,ab,bb,qb,g,einv,status)
        if (status /= linalg_status_success) return
        n = size(a,1)
        m = size(b,2)
        epsilon_tol = default_tol
        if (present(tol)) epsilon_tol = tol
        maxit = default_max_iter
        if (present(max_iter)) maxit = max_iter
        if (.not. ieee_is_finite(epsilon_tol) .or. epsilon_tol <= 0.0_dp .or. maxit < 1) then
            status = linalg_status_invalid_shape
            return
        end if
        allocate(h(n,n), ak(n,n), an(n,n), gn(n,n), hn(n,n), t(n,n), u(n,n), y(n,n))
        h = qb
        ak = ab
        converged = .false.
        do iter = 1, maxit
            ! A_next = A (I+GH)^-1 A
            call solve_dense(identity(n)+matmul(g,h),ak,t,code)
            if (code /= linalg_status_success) exit
            an = matmul(ak,t)
            ! H_next = H + A^H (I+HG)^-1 H A
            call solve_dense(identity(n)+matmul(h,g),matmul(h,ak),u,code)
            if (code /= linalg_status_success) exit
            hn = hermitian_part(h+matmul(adjoint(ak),u))
            ! G_next = G + A G (I+HG)^-1 A^H
            call solve_dense(identity(n)+matmul(h,g),adjoint(ak),u,code)
            if (code /= linalg_status_success) exit
            gn = hermitian_part(g+matmul(matmul(ak,g),u))
            if (.not. (finite_matrix(an) .and. finite_matrix(hn) .and. finite_matrix(gn))) exit
            err = fnorm(hn-h)/max(1.0_dp,fnorm(hn))
            ak = an
            h = hn
            g = gn
            if (err <= epsilon_tol .and. fnorm(ak) <= 1.0_dp) then
                converged = .true.
                exit
            end if
        end do
        if (present(iterations)) iterations = min(iter,maxit)
        if (.not. converged) then
            status = linalg_status_no_convergence
            return
        end if
        y = h
        allocate(delta(n,n), f(n,n), rb(m,m), k(m,n), acl(n,n))
        ! Newton refinement by discrete Lyapunov solution around the feedback.
        do refine = 1, 6
            rb = r+matmul(matmul(adjoint(bb),y),bb)
            call solve_dense(rb,matmul(matmul(adjoint(bb),y),ab),k,code)
            if (code /= linalg_status_success) exit
            f = matmul(matmul(adjoint(ab),y),ab)-y+qb - &
                matmul(matmul(matmul(adjoint(ab),y),bb),k)
            if (fnorm(f)/max(1.0_dp,fnorm(qb),fnorm(y),fnorm(matmul(adjoint(ab),matmul(y,ab)))) &
                    <= max(2.0e-12_dp,4.0_dp*epsilon_tol)) exit
            acl = ab-matmul(bb,k)
            call solve_discrete_lyapunov(adjoint(acl), f, delta, code)
            if (code /= linalg_status_success) exit
            if (.not. finite_matrix(delta)) exit
            y = hermitian_part(y+delta)
        end do
        x = matmul(matmul(adjoint(einv),y),einv)
        x = hermitian_part(x)
        call verify_residual(a,b,q,r,x,.false.,e,s,err,code)
        if (present(residual)) residual = err
        rb = r+matmul(matmul(adjoint(bb),y),bb)
        call solve_dense(rb,matmul(matmul(adjoint(bb),y),ab),k,code)
        if (code == linalg_status_success) acl = ab-matmul(bb,k)
        if (code /= linalg_status_success .or. err > max(2.0e-8_dp,128.0_dp*epsilon_tol)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
        else if (fnorm(qb) <= 1.0e-14_dp*max(1.0_dp,fnorm(q)) .and. &
                .not. is_discrete_stable(acl)) then
            x = cmplx(0.0_dp, 0.0_dp, dp)
            status = linalg_status_no_convergence
        else
            status = linalg_status_success
        end if
    end subroutine discrete_are_complex

    subroutine prepare_equation(a,b,q,r,x,e,s,ab,bb,qb,g,einv,status)
        complex(dp), intent(in) :: a(:, :) !! system matrix (n,n)
        complex(dp), intent(in) :: b(:, :) !! input matrix (n,m)
        complex(dp), intent(in) :: q(:, :) !! state cost (n,n)
        complex(dp), intent(in) :: r(:, :) !! control cost (m,m)
        complex(dp), intent(in) :: x(:, :) !! output shape template (n,n)
        complex(dp), intent(in), optional :: e(:, :) !! descriptor matrix (n,n)
        complex(dp), intent(in), optional :: s(:, :) !! cross cost (n,m)
        complex(dp), allocatable, intent(out) :: ab(:, :) !! transformed dynamics (n,n)
        complex(dp), allocatable, intent(out) :: bb(:, :) !! transformed control matrix (n,m)
        complex(dp), allocatable, intent(out) :: qb(:, :) !! transformed state cost (n,n)
        complex(dp), allocatable, intent(out) :: g(:, :) !! B R^-1 B^H (n,n)
        complex(dp), allocatable, intent(out) :: einv(:, :) !! inverse descriptor matrix (n,n)
        integer, intent(out) :: status !! status of checks and factorizations
        complex(dp), allocatable :: rinv(:, :), smat(:, :)
        integer :: n, m
        n = size(a,1)
        m = size(b,2)
        status = linalg_status_invalid_shape
        if (n < 1 .or. m < 1 .or. n > riccati_max_states) return
        if (n > matrix_equation_max_unknowns/n) return
        if (size(a,2) /= n .or. size(b,1) /= n) return
        if (size(q,1) /= n .or. size(q,2) /= n) return
        if (size(r,1) /= m .or. size(r,2) /= m) return
        if (size(x,1) /= n .or. size(x,2) /= n) return
        if (present(e)) then
            if (size(e,1) /= n .or. size(e,2) /= n) return
        end if
        if (present(s)) then
            if (size(s,1) /= n .or. size(s,2) /= m) return
        end if
        status = linalg_status_no_convergence
        if (.not. (finite_matrix(a) .and. finite_matrix(b) .and. &
            finite_matrix(q) .and. finite_matrix(r))) return
        if (present(e)) then
            if (.not. finite_matrix(e)) return
        end if
        if (present(s)) then
            if (.not. finite_matrix(s)) return
        end if
        status = linalg_status_invalid_shape
        if (fnorm(q-adjoint(q)) > 1.0e-10_dp*max(1.0_dp,fnorm(q))) return
        if (fnorm(r-adjoint(r)) > 1.0e-10_dp*max(1.0_dp,fnorm(r))) return
        if (.not. positive_definite(r)) then
            status = linalg_status_not_positive_definite
            return
        end if
        allocate(einv(n,n),rinv(m,m), smat(n,m),ab(n,n),bb(n,m),qb(n,n),g(n,n))
        if (present(e)) then
            call solve_dense(e,identity(n),einv,status)
            if (status /= linalg_status_success) return
        else
            einv = identity(n)
        end if
        call solve_dense(r,identity(m),rinv,status)
        if (status /= linalg_status_success) return
        smat = cmplx(0.0_dp,0.0_dp,dp)
        if (present(s)) smat = s
        bb = matmul(einv,b)
        ab = matmul(einv,a)-matmul(matmul(bb,rinv),adjoint(smat))
        qb = hermitian_part(q-matmul(matmul(smat,rinv),adjoint(smat)))
        g = hermitian_part(matmul(matmul(bb,rinv),adjoint(bb)))
        status = linalg_status_success
    end subroutine prepare_equation

    subroutine verify_residual(a,b,q,r,x,continuous,e,s,residual,status)
        complex(dp), intent(in) :: a(:, :) !! original dynamics (n,n)
        complex(dp), intent(in) :: b(:, :) !! original controls (n,m)
        complex(dp), intent(in) :: q(:, :) !! original state costs (n,n)
        complex(dp), intent(in) :: r(:, :) !! original control costs (m,m)
        complex(dp), intent(in) :: x(:, :) !! candidate solution (n,n)
        logical, intent(in) :: continuous !! continuous or discrete equation
        complex(dp), intent(in), optional :: e(:, :) !! descriptor (n,n)
        complex(dp), intent(in), optional :: s(:, :) !! cross costs (n,m)
        real(dp), intent(out) :: residual !! normalized equation residual
        integer, intent(out) :: status !! zero if solve completed
        complex(dp), allocatable :: ee(:, :), ss(:, :), left(:, :), right(:, :), rr(:, :)
        complex(dp), allocatable :: k(:, :), f(:, :), aa(:, :)
        integer :: n, m
        n=size(a,1)
        m=size(b,2)
        allocate(ee(n,n),ss(n,m),left(n,m),right(m,n),rr(m,m),k(m,n),f(n,n),aa(n,n))
        ee=identity(n)
        ss=cmplx(0.0_dp,0.0_dp,dp)
        if (present(e)) ee=e
        if (present(s)) ss=s
        if (continuous) then
            left=matmul(matmul(adjoint(ee),x),b)+ss
            right=matmul(matmul(adjoint(b),x),ee)+adjoint(ss)
            rr=r
            aa=matmul(matmul(adjoint(ee),x),a)+matmul(matmul(adjoint(a),x),ee)+q
        else
            left=matmul(matmul(adjoint(a),x),b)+ss
            right=matmul(matmul(adjoint(b),x),a)+adjoint(ss)
            rr=r+matmul(matmul(adjoint(b),x),b)
            aa=matmul(matmul(adjoint(a),x),a)-matmul(matmul(adjoint(ee),x),ee)+q
        end if
        call solve_dense(rr,right,k,status)
        residual=huge(1.0_dp)
        if (status/=linalg_status_success) return
        f=aa-matmul(left,k)
        residual=fnorm(f)/max(1.0_dp,fnorm(q),fnorm(aa),fnorm(matmul(left,k)))
        if (.not. ieee_is_finite(residual)) status=linalg_status_no_convergence
    end subroutine verify_residual

    subroutine solve_dense(a,b,x,status)
        complex(dp), intent(in) :: a(:, :) !! square nonsingular coefficient (n,n)
        complex(dp), intent(in) :: b(:, :) !! matrix right-hand side (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solution (n,nrhs)
        integer, intent(out) :: status !! success, invalid shape, singular, or nonfinite
        complex(dp), allocatable :: u(:, :), rhs(:, :)
        complex(dp) :: pivot_row, factor
        integer :: n, nrhs, i, j, k, pivot
        real(dp) :: scale, pivot_tol

        x=cmplx(0.0_dp,0.0_dp,dp)
        n=size(a,1)
        nrhs=size(b,2)
        status=linalg_status_invalid_shape
        if (n < 1 .or. size(a,2)/=n .or. size(b,1)/=n .or. &
                size(x,1)/=n .or. size(x,2)/=nrhs) return
        status=linalg_status_no_convergence
        if (.not. (finite_matrix(a) .and. finite_matrix(b))) return
        allocate(u(n,n),rhs(n,nrhs))
        u=a
        rhs=b
        scale=maxval(abs(u))
        pivot_tol=max(64.0_dp*tiny(1.0_dp),8.0_dp*epsilon(1.0_dp)*scale)
        do k=1,n
            pivot=k-1+maxloc(abs(u(k:n,k)),dim=1)
            if (abs(u(pivot,k))<=pivot_tol) then
                status=linalg_status_singular
                return
            end if
            if (pivot/=k) then
                do j=k,n
                    pivot_row=u(k,j)
                    u(k,j)=u(pivot,j)
                    u(pivot,j)=pivot_row
                end do
                do j=1,nrhs
                    pivot_row=rhs(k,j)
                    rhs(k,j)=rhs(pivot,j)
                    rhs(pivot,j)=pivot_row
                end do
            end if
            do i=k+1,n
                factor=u(i,k)/u(k,k)
                u(i,k)=cmplx(0.0_dp,0.0_dp,dp)
                u(i,k+1:n)=u(i,k+1:n)-factor*u(k,k+1:n)
                rhs(i,:)=rhs(i,:)-factor*rhs(k,:)
            end do
        end do
        do k=n,1,-1
            x(k,:)=(rhs(k,:)-matmul(u(k,k+1:n),x(k+1:n,:)))/u(k,k)
        end do
        status=linalg_status_success
        if (.not. finite_matrix(x)) status=linalg_status_no_convergence
    end subroutine solve_dense

    logical function is_discrete_stable(a)
        complex(dp), intent(in) :: a(:, :) !! closed-loop discrete-time transition
        complex(dp), allocatable :: power(:, :)
        integer :: i
        power = a
        is_discrete_stable = .false.
        do i = 1, 12
            if (.not. finite_matrix(power)) return
            if (fnorm(power) < 0.8_dp) then
                is_discrete_stable = .true.
                return
            end if
            if (fnorm(power) > 1.0e100_dp) return
            power = matmul(power,power)
        end do
        if (finite_matrix(power)) is_discrete_stable = fnorm(power) < 0.8_dp
    end function is_discrete_stable

    pure logical function positive_definite(a)
        complex(dp), intent(in) :: a(:, :) !! Hermitian matrix (n,n)
        complex(dp), allocatable :: l(:, :)
        complex(dp) :: val
        integer :: n, i, j
        real(dp) :: bound
        n=size(a,1)
        allocate(l(n,n))
        l=cmplx(0.0_dp,0.0_dp,dp)
        positive_definite=.false.
        bound=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
        do i=1,n
            val=a(i,i)
            if (i>1) val=val-sum(l(i,1:i-1)*conjg(l(i,1:i-1)))
            if (abs(aimag(val))>bound .or. real(val,dp)<=bound) return
            l(i,i)=cmplx(sqrt(real(val,dp)),0.0_dp,dp)
            do j=i+1,n
                val=a(j,i)
                if (i>1) val=val-sum(l(j,1:i-1)*conjg(l(i,1:i-1)))
                l(j,i)=val/l(i,i)
            end do
        end do
        positive_definite=.true.
    end function positive_definite

    pure function identity(n) result(a)
        integer, intent(in) :: n !! matrix dimension
        complex(dp) :: a(n,n) !! identity matrix
        integer :: i
        a=cmplx(0.0_dp,0.0_dp,dp)
        do i=1,n
            a(i,i)=cmplx(1.0_dp,0.0_dp,dp)
        end do
    end function identity

    pure function adjoint(a) result(b)
        complex(dp), intent(in) :: a(:, :) !! complex matrix
        complex(dp) :: b(size(a,2),size(a,1)) !! conjugate transpose
        b=transpose(conjg(a))
    end function adjoint

    pure function hermitian_part(a) result(b)
        complex(dp), intent(in) :: a(:, :) !! nearly Hermitian matrix
        complex(dp) :: b(size(a,1),size(a,2)) !! symmetrized matrix
        b=0.5_dp*(a+adjoint(a))
    end function hermitian_part

    pure real(dp) function fnorm(a)
        complex(dp), intent(in) :: a(:, :) !! input matrix
        fnorm=sqrt(sum(abs(a)**2))
    end function fnorm

    pure logical function finite_matrix(a)
        complex(dp), intent(in) :: a(:, :) !! input matrix
        finite_matrix=all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix

end module scifort_linalg_riccati
