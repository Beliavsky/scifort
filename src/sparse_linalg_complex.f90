! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! Complex binary64 CG, BiCG, BiCGSTAB, QMR and restarted GMRES for CSR/CSC.
! All solvers use true residual checks and conjugating DOT_PRODUCT.
! Complex QMR extends the existing SciFort real QMR recurrence, which adapts
! BSD-3-Clause SciPy QMR. Preserve LICENSES/scipy-qmr-BSD-3-Clause.txt and
! CODE_PROVENANCE.md. Other new complex solver/preconditioner code follows
! standard published methods independently.
module scifort_sparse_linalg_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr
    use scifort_sparse_complex_ilu, only : complex_ilu0_factor, &
        valid_ilu0_complex, ilu0_complex_apply
    implicit none
    private

    integer, parameter, public :: sparse_complex_linalg_success = 0
    integer, parameter, public :: sparse_complex_linalg_invalid_input = -1
    integer, parameter, public :: sparse_complex_linalg_breakdown = -2
    integer, parameter, public :: sparse_complex_linalg_bad_diagonal = -3
    public :: cg_complex, bicgstab_complex, gmres_complex, bicg_complex, qmr_complex

    interface bicg_complex
        module procedure bicg_complex_csr, bicg_complex_csc
    end interface bicg_complex
    interface qmr_complex
        module procedure qmr_complex_csr, qmr_complex_csc
    end interface qmr_complex
    interface cg_complex
        module procedure cg_complex_csr, cg_complex_csc
    end interface cg_complex
    interface bicgstab_complex
        module procedure bicgstab_complex_csr, bicgstab_complex_csc
    end interface bicgstab_complex
    interface gmres_complex
        module procedure gmres_complex_csr, gmres_complex_csc
    end interface gmres_complex

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate complex scalar.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    ! Scaled norm avoids overflow when forming abs(z)**2 for large components.
    pure real(dp) function complex_norm(x) result(value)
        complex(dp), intent(in) :: x(:) !! Complex vector.
        real(dp) :: scale, ssq, a
        integer :: i, part

        scale = 0.0_dp
        ssq = 1.0_dp
        do i = 1, size(x)
            do part = 1, 2
                if (part == 1) then
                    a = abs(real(x(i),dp))
                else
                    a = abs(aimag(x(i)))
                end if
                if (a == 0.0_dp) cycle
                if (scale < a) then
                    ssq = 1.0_dp+ssq*(scale/a)**2
                    scale = a
                else
                    ssq = ssq+(a/scale)**2
                end if
            end do
        end do
        value = scale*sqrt(ssq)
    end function complex_norm

    ! Internal product assumes validated matrix/vector and avoids allocating.
    subroutine product(a, x, y)
        type(csr_complex_matrix), intent(in) :: a !! Prevalidated CSR matrix.
        complex(dp), intent(in) :: x(:) !! Input vector of length n_cols.
        complex(dp), intent(out) :: y(:) !! Result vector of length n_rows.
        integer :: i, k

        y = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                y(i) = y(i) + a%data(k)*x(a%indices(k))
            end do
        end do
    end subroutine product

    subroutine prepare(a, b, x, info, threshold, cap, invdiag, rtol, atol, maxiter, x0, jacobi, ilu)
        type(csr_complex_matrix), intent(in) :: a !! Canonical square CSR matrix.
        complex(dp), intent(in) :: b(:) !! Right hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Initial guess or zero vector.
        integer, intent(out) :: info !! Status code.
        real(dp), intent(out) :: threshold !! Absolute true residual threshold.
        integer, intent(out) :: cap !! Maximum iterations or restart cycles.
        complex(dp), allocatable, intent(out) :: invdiag(:) !! Inverse Jacobi diagonal, or ones.
        real(dp), intent(in), optional :: rtol !! Relative tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Iteration/cycle cap, default max(1,10*n).
        complex(dp), intent(in), optional :: x0(:) !! Starting guess.
        logical, intent(in), optional :: jacobi !! Use diagonal preconditioning.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0), exclusive with Jacobi.
        real(dp) :: reltol, abstol
        integer :: n, i, k
        logical :: enabled, found

        info = sparse_complex_linalg_invalid_input
        threshold = huge(1.0_dp)
        cap = 0
        allocate(x(0),invdiag(0))
        if (.not. valid_csr_complex(a)) return
        if (a%n_rows /= a%n_cols .or. size(b) /= a%n_rows) return
        if (.not. all(finite_complex(b))) return
        n = size(b)
        if (present(x0)) then
            if (size(x0) /= n) return
            if (.not. all(finite_complex(x0))) return
        end if
        reltol = 1.0e-5_dp
        abstol = 0.0_dp
        if (present(rtol)) reltol = rtol
        if (present(atol)) abstol = atol
        if (.not. ieee_is_finite(reltol) .or. .not. ieee_is_finite(abstol)) return
        if (reltol < 0.0_dp .or. abstol < 0.0_dp) return
        if (n > huge(0)/10) then
            cap = huge(0)
        else
            cap = max(1,10*n)
        end if
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        enabled = .false.
        if (present(jacobi)) enabled = jacobi
        if (present(ilu)) then
            if (enabled) return
            if (.not. valid_ilu0_complex(ilu)) return
            if (ilu%lu%n_rows /= n) return
        end if
        deallocate(x,invdiag)
        allocate(x(n),invdiag(n))
        x = cmplx(0.0_dp,0.0_dp,dp)
        if (present(x0)) x = x0
        invdiag = cmplx(1.0_dp,0.0_dp,dp)
        threshold = max(abstol,reltol*complex_norm(b))
        if (.not. ieee_is_finite(threshold)) then
            info = sparse_complex_linalg_invalid_input
            return
        end if
        if (enabled) then
            do i = 1, n
                found = .false.
                do k = a%indptr(i),a%indptr(i+1)-1
                    if (a%indices(k) /= i) cycle
                    if (abs(a%data(k)) == 0.0_dp) exit
                    invdiag(i) = 1.0_dp/a%data(k)
                    found = finite_complex(invdiag(i))
                    exit
                end do
                if (.not. found) then
                    info = sparse_complex_linalg_bad_diagonal
                    return
                end if
            end do
        end if
        info = sparse_complex_linalg_success
    end subroutine prepare

    ! Apply M^(-1), or M^(-H) for the shadow system in BiCG/QMR.
    pure subroutine precondition_complex(invdiag, v, z, ilu, adjoint)
        complex(dp), intent(in) :: invdiag(:) !! Jacobi inverse diagonal, or all ones.
        complex(dp), intent(in) :: v(:) !! Input vector.
        complex(dp), intent(out) :: z(:) !! Preconditioned vector.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Validated complex ILU(0).
        logical, intent(in), optional :: adjoint !! Use adjoint of inverse preconditioner.
        logical :: transpose_apply

        transpose_apply = .false.
        if (present(adjoint)) transpose_apply = adjoint
        if (present(ilu)) then
            call ilu0_complex_apply(ilu,v,z,transpose_apply)
        else if (transpose_apply) then
            z = conjg(invdiag)*v
        else
            z = invdiag*v
        end if
    end subroutine precondition_complex

    subroutine product_adjoint(a, x, y)
        type(csr_complex_matrix), intent(in) :: a !! Canonical CSR matrix.
        complex(dp), intent(in) :: x(:) !! Vector of length n_rows.
        complex(dp), intent(out) :: y(:) !! Conjugate-transposed product, length n_cols.
        integer :: i, k

        y = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, a%n_rows
            do k = a%indptr(i),a%indptr(i+1)-1
                y(a%indices(k)) = y(a%indices(k))+conjg(a%data(k))*x(i)
            end do
        end do
    end subroutine product_adjoint

    subroutine cg_complex_csr(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm)
        type(csr_complex_matrix), intent(in) :: a !! Hermitian positive-definite matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Computed solution.
        integer, intent(out) :: info !! 0 success, positive no convergence, negative error.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance.
        integer, intent(in), optional :: maxiter !! Maximum CG iterations.
        complex(dp), intent(in), optional :: x0(:) !! Starting guess.
        logical, intent(in), optional :: jacobi !! Use diagonal preconditioner.
        integer, intent(out), optional :: iterations !! Number of CG steps.
        real(dp), intent(out), optional :: residual_norm !! Final true-residual norm.
        complex(dp), allocatable :: invdiag(:), r(:), z(:), p(:), ap(:), ax(:)
        complex(dp) :: rho, rho_new, denom, beta
        real(dp) :: threshold, err, alpha, denom_real, rho_real
        integer :: cap, n, it, steps

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,threshold,cap,invdiag,rtol,atol,maxiter,x0,jacobi)
        if (info /= 0) return
        n = size(b)
        allocate(r(n),z(n),p(n),ap(n),ax(n))
        call product(a,x,ax)
        r = b-ax
        err = complex_norm(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_complex_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            z = invdiag*r
            rho = dot_product(r,z)
            p = z
            do it = 1, cap
                steps = it
                rho_real = real(rho,dp)
                if (.not. finite_complex(rho) .or. rho_real <= 0.0_dp) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (abs(aimag(rho)) > 1.0e-8_dp*abs(rho_real)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                call product(a,p,ap)
                denom = dot_product(p,ap)
                denom_real = real(denom,dp)
                if (.not. finite_complex(denom) .or. denom_real <= 0.0_dp) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (abs(aimag(denom)) > 1.0e-8_dp*abs(denom_real)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                alpha = rho_real/denom_real
                x = x+alpha*p
                r = r-alpha*ap
                err = complex_norm(r)
                if (.not. all(finite_complex(x)) .or. .not. ieee_is_finite(err)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    call product(a,x,ax)
                    err = complex_norm(b-ax)
                    if (err <= threshold) then
                        info = 0
                        exit
                    end if
                    r = b-ax
                end if
                z = invdiag*r
                rho_new = dot_product(r,z)
                if (.not. finite_complex(rho_new) .or. rho_new == 0.0_dp) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                beta = rho_new/rho
                if (.not. finite_complex(beta)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                p = z+beta*p
                rho = rho_new
            end do
        end if
        call product(a,x,ax)
        err = complex_norm(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine cg_complex_csr

    subroutine bicgstab_complex_csr(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm, ilu)
        type(csr_complex_matrix), intent(in) :: a !! Nonsymmetric complex CSR matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Computed solution.
        integer, intent(out) :: info !! 0 success, positive no convergence, negative error.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance.
        integer, intent(in), optional :: maxiter !! Maximum BiCGSTAB steps.
        complex(dp), intent(in), optional :: x0(:) !! Initial estimate.
        logical, intent(in), optional :: jacobi !! Right Jacobi preconditioner.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Validated ILU(0), exclusive with Jacobi.
        integer, intent(out), optional :: iterations !! Number of steps.
        real(dp), intent(out), optional :: residual_norm !! Final true-residual norm.
        complex(dp), allocatable :: invdiag(:), r(:), rhat(:), p(:), v(:), s(:), t(:), &
            phat(:), shat(:), ax(:)
        complex(dp) :: alpha, omega, rho, rho_prev, beta, denom
        real(dp) :: threshold, err
        integer :: cap, n, it, steps

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,threshold,cap,invdiag,rtol,atol,maxiter,x0,jacobi,ilu)
        if (info /= 0) return
        n = size(b)
        allocate(r(n),rhat(n),p(n),v(n),s(n),t(n),phat(n),shat(n),ax(n))
        call product(a,x,ax)
        r = b-ax
        rhat = r
        p = cmplx(0.0_dp,0.0_dp,dp)
        v = p
        alpha = cmplx(1.0_dp,0.0_dp,dp)
        omega = alpha
        rho_prev = alpha
        err = complex_norm(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_complex_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            do it = 1, cap
                steps = it
                rho = dot_product(rhat,r)
                if (rho == 0.0_dp .or. .not. finite_complex(rho)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (it == 1) then
                    p = r
                else
                    if (omega == 0.0_dp .or. .not. finite_complex(omega)) then
                        info = sparse_complex_linalg_breakdown
                        exit
                    end if
                    beta = (rho/rho_prev)*(alpha/omega)
                    p = r+beta*(p-omega*v)
                end if
                call precondition_complex(invdiag,p,phat,ilu)
                call product(a,phat,v)
                denom = dot_product(rhat,v)
                if (denom == 0.0_dp .or. .not. finite_complex(denom)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                alpha = rho/denom
                s = r-alpha*v
                if (complex_norm(s) <= threshold) then
                    x = x+alpha*phat
                    call product(a,x,ax)
                    err = complex_norm(b-ax)
                    if (err <= threshold) then
                        info = 0
                    else
                        info = sparse_complex_linalg_breakdown
                    end if
                    exit
                end if
                call precondition_complex(invdiag,s,shat,ilu)
                call product(a,shat,t)
                denom = dot_product(t,t)
                if (denom == 0.0_dp .or. .not. finite_complex(denom)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                omega = dot_product(t,s)/denom
                if (omega == 0.0_dp .or. .not. finite_complex(omega)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                x = x+alpha*phat+omega*shat
                r = s-omega*t
                err = complex_norm(r)
                if (.not. ieee_is_finite(err) .or. .not. all(finite_complex(x))) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    call product(a,x,ax)
                    err = complex_norm(b-ax)
                    if (err <= threshold) then
                        info = 0
                        exit
                    end if
                    r = b-ax
                end if
                rho_prev = rho
            end do
        end if
        call product(a,x,ax)
        err = complex_norm(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine bicgstab_complex_csr

    ! Compute complex Givens rotation G*[a,b]^T=[r,0]^T where
    ! G = [ c, s; -conjg(s), c ], with real c >= 0.
    subroutine complex_givens(a, b, c, s, r)
        complex(dp), intent(in) :: a !! Upper entry.
        complex(dp), intent(in) :: b !! Lower entry.
        real(dp), intent(out) :: c !! Real cosine.
        complex(dp), intent(out) :: s !! Complex sine.
        complex(dp), intent(out) :: r !! Rotated upper entry.
        real(dp) :: aa, bb, scale, normab
        complex(dp) :: phase

        aa = abs(a)
        bb = abs(b)
        if (bb == 0.0_dp) then
            c = 1.0_dp
            s = cmplx(0.0_dp,0.0_dp,dp)
            r = a
        else if (aa == 0.0_dp) then
            c = 0.0_dp
            s = conjg(b)/bb
            r = cmplx(bb,0.0_dp,dp)
        else
            scale = max(aa,bb)
            normab = scale*sqrt((aa/scale)**2+(bb/scale)**2)
            phase = a/aa
            c = aa/normab
            s = phase*conjg(b)/normab
            r = phase*normab
        end if
    end subroutine complex_givens

    subroutine gmres_complex_csr(a, b, x, info, rtol, atol, maxiter, restart, x0, jacobi, iterations, residual_norm, ilu)
        type(csr_complex_matrix), intent(in) :: a !! General complex CSR matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Computed solution.
        integer, intent(out) :: info !! 0 success, positive no convergence, negative error.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance.
        integer, intent(in), optional :: maxiter !! Maximum restart cycles.
        integer, intent(in), optional :: restart !! Arnoldi dimension, default min(20,n).
        complex(dp), intent(in), optional :: x0(:) !! Starting guess.
        logical, intent(in), optional :: jacobi !! Apply right Jacobi preconditioner.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Validated ILU(0), exclusive with Jacobi.
        integer, intent(out), optional :: iterations !! Number of Arnoldi steps.
        real(dp), intent(out), optional :: residual_norm !! Final true-residual norm.
        complex(dp), allocatable :: invdiag(:), v(:,:), z(:,:), h(:,:), g(:), y(:), &
            w(:), r(:), ax(:), xtry(:), s(:)
        real(dp), allocatable :: c(:)
        complex(dp) :: temp, diag
        real(dp) :: threshold, err, beta, hnext
        integer :: cap, n, m, it, j, i, pass, steps, completed
        logical :: finish, bad

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,threshold,cap,invdiag,rtol,atol,maxiter,x0,jacobi,ilu)
        if (info /= 0) return
        n = size(b)
        m = min(20,n)
        if (present(restart)) m = restart
        if (m <= 0 .and. n > 0) then
            info = sparse_complex_linalg_invalid_input
            return
        end if
        m = min(m,n)
        allocate(v(n,m+1),z(n,m),h(m+1,m),g(m+1),y(m),w(n),r(n),ax(n),xtry(n),c(m),s(m))
        call product(a,x,ax)
        r = b-ax
        err = complex_norm(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_complex_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            do it = 1, cap
                call product(a,x,ax)
                r = b-ax
                beta = complex_norm(r)
                if (beta == 0.0_dp .or. .not. ieee_is_finite(beta)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                v(:,1) = r/beta
                h = cmplx(0.0_dp,0.0_dp,dp)
                g = cmplx(0.0_dp,0.0_dp,dp)
                g(1) = cmplx(beta,0.0_dp,dp)
                completed = 0
                finish = .false.
                do j = 1, m
                    steps = steps+1
                    completed = j
                    call precondition_complex(invdiag,v(:,j),z(:,j),ilu)
                    call product(a,z(:,j),w)
                    if (.not. all(finite_complex(w))) then
                        info = sparse_complex_linalg_breakdown
                        finish = .true.
                        exit
                    end if
                    ! Twice modified Gram-Schmidt for complex Arnoldi.
                    do pass = 1, 2
                        do i = 1, j
                            temp = dot_product(v(:,i),w)
                            h(i,j) = h(i,j)+temp
                            w = w-temp*v(:,i)
                        end do
                    end do
                    hnext = complex_norm(w)
                    h(j+1,j) = cmplx(hnext,0.0_dp,dp)
                    if (hnext > 0.0_dp) v(:,j+1) = w/hnext
                    do i = 1, j-1
                        temp = c(i)*h(i,j)+s(i)*h(i+1,j)
                        h(i+1,j) = -conjg(s(i))*h(i,j)+c(i)*h(i+1,j)
                        h(i,j) = temp
                    end do
                    call complex_givens(h(j,j),h(j+1,j),c(j),s(j),diag)
                    if (.not. finite_complex(diag) .or. diag == 0.0_dp) then
                        info = sparse_complex_linalg_breakdown
                        finish = .true.
                        exit
                    end if
                    h(j,j) = diag
                    h(j+1,j) = cmplx(0.0_dp,0.0_dp,dp)
                    temp = c(j)*g(j)+s(j)*g(j+1)
                    g(j+1) = -conjg(s(j))*g(j)+c(j)*g(j+1)
                    g(j) = temp
                    y(:j) = g(:j)
                    bad = .false.
                    do i = j, 1, -1
                        if (h(i,i) == 0.0_dp) then
                            bad = .true.
                            exit
                        end if
                        ! SUM instead of DOT_PRODUCT: no conjugation in back substitution.
                        y(i) = (y(i)-sum(h(i,i+1:j)*y(i+1:j)))/h(i,i)
                    end do
                    if (bad .or. .not. all(finite_complex(y(:j)))) then
                        info = sparse_complex_linalg_breakdown
                        finish = .true.
                        exit
                    end if
                    xtry = x+matmul(z(:,:j),y(:j))
                    call product(a,xtry,ax)
                    err = complex_norm(b-ax)
                    if (.not. ieee_is_finite(err) .or. .not. all(finite_complex(xtry))) then
                        info = sparse_complex_linalg_breakdown
                        finish = .true.
                        exit
                    end if
                    if (err <= threshold) then
                        x = xtry
                        info = 0
                        finish = .true.
                        exit
                    end if
                    if (hnext <= epsilon(1.0_dp)*beta) then
                        ! Arnoldi happy breakdown without true convergence.
                        info = sparse_complex_linalg_breakdown
                        finish = .true.
                        exit
                    end if
                end do
                if (finish) exit
                if (completed > 0) x = xtry
            end do
        end if
        call product(a,x,ax)
        err = complex_norm(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine gmres_complex_csr

    ! Preconditioned complex BiCG. For A and A^H the shadow recurrence must
    ! conjugate alpha and beta; DOT_PRODUCT conjugates its first argument.
    subroutine bicg_complex_csr(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
        type(csr_complex_matrix), intent(in) :: a !! General nonsymmetric CSR matrix.
        complex(dp), intent(in) :: b(:) !! RHS of Ax=b.
        complex(dp), allocatable, intent(out) :: x(:) !! Final approximation.
        integer, intent(out) :: info !! 0 success; positive cap; negative failure.
        real(dp), intent(in), optional :: rtol !! Relative tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute tolerance, default zero.
        integer, intent(in), optional :: maxiter !! Maximum BiCG iterations.
        complex(dp), intent(in), optional :: x0(:) !! Initial approximation.
        logical, intent(in), optional :: jacobi !! Use Jacobi instead of identity.
        integer, intent(out), optional :: iterations !! Executed iterations.
        real(dp), intent(out), optional :: residual_norm !! True Euclidean residual norm.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0), exclusive with Jacobi.
        complex(dp), allocatable :: invdiag(:),r(:),rh(:),z(:),zh(:),p(:),ph(:),&
            ap(:),aph(:),ax(:)
        complex(dp) :: rho, rho_next, alpha, beta, denominator
        real(dp) :: threshold, err
        integer :: cap, n, it, steps

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,threshold,cap,invdiag,rtol,atol,maxiter,x0,jacobi,ilu)
        if (info /= sparse_complex_linalg_success) return
        n = size(b)
        allocate(r(n),rh(n),z(n),zh(n),p(n),ph(n),ap(n),aph(n),ax(n))
        call product(a,x,ax)
        r = b-ax
        err = complex_norm(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_complex_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_complex_linalg_success
        else
            info = max(1,cap)
            rh = r
            call precondition_complex(invdiag,r,z,ilu)
            call precondition_complex(invdiag,rh,zh,ilu,adjoint=.true.)
            rho = dot_product(rh,z)
            p = z
            ph = zh
            do it = 1, cap
                steps = it
                if (.not. finite_complex(rho) .or. abs(rho) <= tiny(1.0_dp)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (.not. all(finite_complex(p)) .or. .not. all(finite_complex(ph))) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                call product(a,p,ap)
                call product_adjoint(a,ph,aph)
                denominator = dot_product(ph,ap)
                if (.not. finite_complex(denominator) .or. abs(denominator) <= tiny(1.0_dp)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                alpha = rho/denominator
                x = x+alpha*p
                r = r-alpha*ap
                rh = rh-conjg(alpha)*aph
                err = complex_norm(r)
                if (.not. all(finite_complex(x)) .or. .not. ieee_is_finite(err)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    call product(a,x,ax)
                    err = complex_norm(b-ax)
                    if (err <= threshold) then
                        info = sparse_complex_linalg_success
                        exit
                    end if
                    ! Restarting the residual without restarting the dual
                    ! recurrence would violate biorthogonality: do not do it.
                end if
                call precondition_complex(invdiag,r,z,ilu)
                call precondition_complex(invdiag,rh,zh,ilu,adjoint=.true.)
                rho_next = dot_product(rh,z)
                if (.not. finite_complex(rho_next) .or. abs(rho_next) <= tiny(1.0_dp)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                beta = rho_next/rho
                p = z+beta*p
                ph = zh+conjg(beta)*ph
                rho = rho_next
            end do
        end if
        call product(a,x,ax)
        err = complex_norm(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_complex_linalg_success
        end if
    end subroutine bicg_complex_csr

    ! Complex two-sided Lanczos QMR: adaptation of SciFort's real QMR
    ! (BSD-3-Clause SciPy attribution retained in LICENSES/).
    ! For stability, the algorithm
    ! works on B=M^(-1)A, whose adjoint is B^H=A^H*M^(-H). All reported
    ! stopping criteria use the true residual of the original A*x=b.
    subroutine qmr_complex_csr(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
        type(csr_complex_matrix), intent(in) :: a !! General nonsymmetric CSR matrix.
        complex(dp), intent(in) :: b(:) !! RHS of Ax=b.
        complex(dp), allocatable, intent(out) :: x(:) !! Final approximation.
        integer, intent(out) :: info !! 0 success; positive cap; negative failure.
        real(dp), intent(in), optional :: rtol !! Relative tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute tolerance, default zero.
        integer, intent(in), optional :: maxiter !! Maximum QMR iterations.
        complex(dp), intent(in), optional :: x0(:) !! Starting guess.
        logical, intent(in), optional :: jacobi !! Use Jacobi preconditioning.
        integer, intent(out), optional :: iterations !! Number of QMR steps.
        real(dp), intent(out), optional :: residual_norm !! Final true residual norm.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0), exclusive with Jacobi.
        complex(dp), allocatable :: invdiag(:), r(:),vt(:),wt(:),y(:),z(:), &
            v(:),w(:),p(:),q(:),pt(:),yt(:),zt(:),d(:),ax(:),temp(:)
        complex(dp) :: delta, eps_prev, eps_now, beta, eta
        real(dp) :: threshold, err, rho, rho_prev, xi, theta, theta_prev
        real(dp) :: gamma, gamma_prev, weight
        integer :: n, cap, it, steps

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,threshold,cap,invdiag,rtol,atol,maxiter,x0,jacobi,ilu)
        if (info /= sparse_complex_linalg_success) return
        n = size(b)
        allocate(r(n),vt(n),wt(n),y(n),z(n),v(n),w(n),p(n),q(n), &
            pt(n),yt(n),zt(n),d(n),ax(n),temp(n))
        call product(a,x,ax)
        err = complex_norm(b-ax)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_complex_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_complex_linalg_success
        else
            info = max(1,cap)
            call precondition_complex(invdiag,b-ax,r,ilu)
            vt = r
            wt = r
            y = vt
            z = wt
            rho = complex_norm(y)
            xi = complex_norm(z)
            gamma = 1.0_dp
            eta = cmplx(-1.0_dp,0.0_dp,dp)
            theta = 0.0_dp
            eps_prev = cmplx(0.0_dp,0.0_dp,dp)
            d = cmplx(0.0_dp,0.0_dp,dp)
            do it = 1, cap
                steps = it
                if (rho <= tiny(1.0_dp) .or. xi <= tiny(1.0_dp) .or. &
                    .not. ieee_is_finite(rho) .or. .not. ieee_is_finite(xi)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                v = vt/rho
                y = y/rho
                w = wt/xi
                z = z/xi
                delta = dot_product(z,y)
                if (.not. finite_complex(delta) .or. abs(delta) <= tiny(1.0_dp)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                yt = y
                zt = z
                if (it > 1) then
                    if (abs(eps_prev) <= tiny(1.0_dp)) then
                        info = sparse_complex_linalg_breakdown
                        exit
                    end if
                    yt = yt-(xi*delta/eps_prev)*p
                    zt = zt-(rho*conjg(delta)/conjg(eps_prev))*q
                end if
                p = yt
                q = zt
                call product(a,p,ax)
                call precondition_complex(invdiag,ax,pt,ilu)
                eps_now = dot_product(q,pt)
                if (.not. finite_complex(eps_now) .or. abs(eps_now) <= tiny(1.0_dp)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                beta = eps_now/delta
                if (.not. finite_complex(beta) .or. abs(beta) <= tiny(1.0_dp)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                vt = pt-beta*v
                y = vt
                rho_prev = rho
                rho = complex_norm(y)
                call precondition_complex(invdiag,q,temp,ilu,adjoint=.true.)
                call product_adjoint(a,temp,wt)
                wt = wt-conjg(beta)*w
                z = wt
                xi = complex_norm(z)
                gamma_prev = gamma
                theta_prev = theta
                theta = rho/(gamma_prev*abs(beta))
                gamma = 1.0_dp/hypot(1.0_dp,theta)
                if (gamma == 0.0_dp .or. .not. ieee_is_finite(gamma)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                eta = -eta*(rho_prev/beta)*(gamma/gamma_prev)**2
                weight = (theta_prev*gamma)**2
                if (.not. finite_complex(eta) .or. .not. ieee_is_finite(weight)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                d = eta*p+weight*d
                x = x+d
                if (.not. all(finite_complex(x))) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                call product(a,x,ax)
                err = complex_norm(b-ax)
                if (.not. ieee_is_finite(err)) then
                    info = sparse_complex_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    info = sparse_complex_linalg_success
                    exit
                end if
                eps_prev = eps_now
            end do
        end if
        call product(a,x,ax)
        err = complex_norm(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_complex_linalg_success
        end if
    end subroutine qmr_complex_csr

    ! CSC overloads use a sparse CSR conversion once per call.
    subroutine bicg_complex_csc(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
        type(csc_complex_matrix), intent(in) :: a !! Canonical complex CSC matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution.
        integer, intent(out) :: info !! Convergence status.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: maxiter !! Iteration limit.
        complex(dp), intent(in), optional :: x0(:) !! Starting guess.
        logical, intent(in), optional :: jacobi !! Jacobi preconditioner.
        integer, intent(out), optional :: iterations !! Number of iterations.
        real(dp), intent(out), optional :: residual_norm !! Norm of original residual.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0) preconditioner.
        type(csr_complex_matrix) :: csr
        integer :: code

        info = sparse_complex_linalg_invalid_input
        allocate(x(0))
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,code)
        if (code /= sparse_success) return
        deallocate(x)
        call bicg_complex_csr(csr,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
    end subroutine bicg_complex_csc

    subroutine qmr_complex_csc(a,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
        type(csc_complex_matrix), intent(in) :: a !! Canonical complex CSC matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution.
        integer, intent(out) :: info !! Convergence status.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: maxiter !! Iteration limit.
        complex(dp), intent(in), optional :: x0(:) !! Starting guess.
        logical, intent(in), optional :: jacobi !! Jacobi preconditioner.
        integer, intent(out), optional :: iterations !! Number of iterations.
        real(dp), intent(out), optional :: residual_norm !! Norm of original residual.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0) preconditioner.
        type(csr_complex_matrix) :: csr
        integer :: code

        info = sparse_complex_linalg_invalid_input
        allocate(x(0))
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,code)
        if (code /= sparse_success) return
        deallocate(x)
        call qmr_complex_csr(csr,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
    end subroutine qmr_complex_csc

    ! CSC wrappers convert once to CSR; all arithmetic remains sparse.
    subroutine cg_complex_csc(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution.
        integer, intent(out) :: info !! Solver status.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: maxiter !! Maximum steps.
        complex(dp), intent(in), optional :: x0(:) !! Initial estimate.
        logical, intent(in), optional :: jacobi !! Diagonal preconditioner.
        integer, intent(out), optional :: iterations !! Steps.
        real(dp), intent(out), optional :: residual_norm !! True residual norm.
        type(csr_complex_matrix) :: csr
        integer :: code

        info = sparse_complex_linalg_invalid_input
        allocate(x(0))
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,code)
        if (code /= sparse_success) return
        deallocate(x)
        call cg_complex_csr(csr,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm)
    end subroutine cg_complex_csc

    subroutine bicgstab_complex_csc(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm, ilu)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution.
        integer, intent(out) :: info !! Solver status.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: maxiter !! Maximum steps.
        complex(dp), intent(in), optional :: x0(:) !! Initial estimate.
        logical, intent(in), optional :: jacobi !! Diagonal preconditioner.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Validated ILU(0), exclusive with Jacobi.
        integer, intent(out), optional :: iterations !! Steps.
        real(dp), intent(out), optional :: residual_norm !! True residual norm.
        type(csr_complex_matrix) :: csr
        integer :: code

        info = sparse_complex_linalg_invalid_input
        allocate(x(0))
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,code)
        if (code /= sparse_success) return
        deallocate(x)
        call bicgstab_complex_csr(csr,b,x,info,rtol,atol,maxiter,x0,jacobi,iterations,residual_norm,ilu)
    end subroutine bicgstab_complex_csc

    subroutine gmres_complex_csc(a, b, x, info, rtol, atol, maxiter, restart, x0, jacobi, iterations, residual_norm, ilu)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC matrix.
        complex(dp), intent(in) :: b(:) !! Right-hand side.
        complex(dp), allocatable, intent(out) :: x(:) !! Solution.
        integer, intent(out) :: info !! Solver status.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: maxiter !! Maximum restart cycles.
        integer, intent(in), optional :: restart !! Arnoldi dimension.
        complex(dp), intent(in), optional :: x0(:) !! Initial estimate.
        logical, intent(in), optional :: jacobi !! Diagonal preconditioner.
        type(complex_ilu0_factor), intent(in), optional :: ilu !! Validated ILU(0), exclusive with Jacobi.
        integer, intent(out), optional :: iterations !! Total Arnoldi steps.
        real(dp), intent(out), optional :: residual_norm !! True residual norm.
        type(csr_complex_matrix) :: csr
        integer :: code

        info = sparse_complex_linalg_invalid_input
        allocate(x(0))
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,code)
        if (code /= sparse_success) return
        deallocate(x)
        call gmres_complex_csr(csr,b,x,info,rtol,atol,maxiter,restart,x0,jacobi,iterations,residual_norm,ilu)
    end subroutine gmres_complex_csc

end module scifort_sparse_linalg_complex
