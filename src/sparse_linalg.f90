! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
!
! Real binary64 iterative sparse solvers for canonical one-based CSR matrices.
! This adapter is independent of the candidate kiss implementation.  The
! LSMR entry point reuses SciFort's existing BSD-3-Clause LSMR implementation.
! General numerical methods: Hestenes--Stiefel CG (1952), van der Vorst
! BiCGSTAB (1992), Saad--Schultz GMRES (1986).
! The QMR and TFQMR recurrences adapt BSD-3-Clause SciPy algorithms;
! see LICENSES/scipy-qmr-BSD-3-Clause.txt,
! LICENSES/scipy-tfqmr-BSD-3-Clause.txt, and CODE_PROVENANCE.md.
! The MINRES Lanczos/QR recurrences follow Paige--Saunders (1975) and are
! adapted with attribution from the BSD-3-Clause SciPy 1.17 MINRES kernel;
! see LICENSES/scipy-minres-BSD-3-Clause.txt and CODE_PROVENANCE.md.
module scifort_sparse_linalg
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, valid_csr, csr_transpose, &
        csc_matrix, valid_csc, csc_to_csr, csr_to_csc, csr_from_coo, &
        csr_eye, csr_matmul, sparse_success
    use scifort_optimize_lsmr, only : lsmr_solve
    implicit none
    private

    integer, parameter, public :: sparse_linalg_success = 0
    integer, parameter, public :: sparse_linalg_invalid_input = -1
    integer, parameter, public :: sparse_linalg_breakdown = -2
    integer, parameter, public :: sparse_linalg_bad_diagonal = -3
    integer, parameter, public :: sparse_linalg_singular = -4
    public :: cg, bicg, cgs, qmr, tfqmr, bicgstab, gmres, lgmres, gcrotmk, minres, eigsh, eigsh_lanczos, eigs, eigs_arnoldi
    public :: lsmr, lsqr, lobpcg
    public :: ilu0_factor, valid_ilu0, ilu0, ilu0_solve, spilu, spilu_solve
    public :: svds, expm_multiply, expm_multiply_interval, onenormest, expm
    public :: funm_multiply_krylov, dense_matrix_function, real_matvec_operator
    public :: spsolve_triangular, spsolve, splu, splu_solve, splu_factor, valid_splu
    public :: matrix_power, inv, norm, norm_axis
    public :: spbandwidth, is_sptriangular

    interface spbandwidth
        module procedure spbandwidth_csr, spbandwidth_csc
    end interface spbandwidth

    interface is_sptriangular
        module procedure is_sptriangular_csr, is_sptriangular_csc
    end interface is_sptriangular

    interface norm
        module procedure sparse_norm_csr
        module procedure sparse_norm_csc
        module procedure sparse_norm_csr_int
        module procedure sparse_norm_csc_int
        module procedure sparse_norm_csr_named
        module procedure sparse_norm_csc_named
    end interface norm

    interface norm_axis
        module procedure sparse_norm_axis_csr
        module procedure sparse_norm_axis_csc
        module procedure sparse_norm_axis_csr_int
        module procedure sparse_norm_axis_csc_int
    end interface norm_axis

    interface matrix_power
        module procedure matrix_power_csr
        module procedure matrix_power_csc
    end interface matrix_power

    interface inv
        module procedure inv_csr
        module procedure inv_csc
    end interface inv

    interface spilu
        module procedure spilu_csr
        module procedure spilu_csc
    end interface spilu

    interface spilu_solve
        module procedure spilu_solve_vector
        module procedure spilu_solve_matrix
    end interface spilu_solve

    interface splu
        module procedure splu_csr
        module procedure splu_csc
    end interface splu

    interface spsolve_triangular
        module procedure spsolve_triangular_vector
        module procedure spsolve_triangular_matrix
        module procedure spsolve_triangular_csc_vector
        module procedure spsolve_triangular_csc_matrix
    end interface spsolve_triangular

    interface splu_solve
        module procedure splu_solve_vector
        module procedure splu_solve_matrix
    end interface splu_solve

    interface spsolve
        module procedure spsolve_vector
        module procedure spsolve_matrix
        module procedure spsolve_csc_vector
        module procedure spsolve_csc_matrix
    end interface spsolve

    ! A dense matrix-function callback receives a small projected matrix and
    ! must fill its same-shaped result.  status=0 indicates success.
    abstract interface
        subroutine dense_matrix_function(matrix, value, status)
            import dp
            real(dp), intent(in) :: matrix(:,:) !! Dense square argument to f.
            real(dp), intent(out) :: value(:,:) !! f(matrix), same square shape.
            integer, intent(out) :: status !! Zero on success, nonzero on failure.
        end subroutine dense_matrix_function

        subroutine real_matvec_operator(x, y)
            import dp
            real(dp), intent(in) :: x(:) !! Vector to which A is applied, length n.
            real(dp), intent(out) :: y(:) !! A*x, length n.
        end subroutine real_matvec_operator
    end interface

    interface funm_multiply_krylov
        module procedure funm_multiply_krylov_csr
        module procedure funm_multiply_krylov_operator
    end interface funm_multiply_krylov

    interface expm_multiply
        module procedure expm_multiply_vector
        module procedure expm_multiply_matrix
    end interface expm_multiply

    interface expm_multiply_interval
        module procedure expm_multiply_interval_vector
        module procedure expm_multiply_interval_matrix
    end interface expm_multiply_interval

    type :: sparse_lu_row
        integer, allocatable :: columns(:) !! Sorted, one-based indices; capacity may exceed used.
        real(dp), allocatable :: values(:) !! Sparse row entries, sorted by column.
        integer :: used = 0 !! Number of occupied entries.
    end type sparse_lu_row

    type :: splu_factor
        type(csr_matrix) :: lu !! Combined unit-lower and upper factors in CSR.
        integer, allocatable :: row_perm(:) !! Permuted RHS: rhs(i) = b(row_perm(i)).
    end type splu_factor

    type :: ilu0_factor
        type(csr_matrix) :: lu  !! Shared sorted sparsity pattern, unit-diagonal L and U.
        integer, allocatable :: diag_ptr(:)  !! Stored one-based U pivots in each row.
    end type ilu0_factor

contains

    ! All solver kernels call this only on a previously validated CSR input.
    ! Do not allocate a temporary vector or rescan the matrix on each iteration.
    pure subroutine csr_product(a, x, y)
        type(csr_matrix), intent(in) :: a !! Canonical CSR operator, prevalidated.
        real(dp), intent(in) :: x(:) !! Vector with a%n_cols entries.
        real(dp), intent(out) :: y(:) !! Output vector with a%n_rows entries.
        integer :: i, k

        y = 0.0_dp
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                y(i) = y(i) + a%data(k)*x(a%indices(k))
            end do
        end do
    end subroutine csr_product

    ! Shared validation for square systems.  An optional boolean requests
    ! Jacobi left preconditioning; a missing/zero diagonal is a distinct error.
    subroutine prepare(a, b, x, info, rtol, atol, maxiter, x0, jacobi, invdiag, threshold, cap, ilu)
        type(csr_matrix), intent(in) :: a !! Candidate square CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! RHS of length n with finite entries.
        real(dp), allocatable, intent(out) :: x(:) !! Initial iterate; zero if no x0.
        integer, intent(out) :: info !! Zero when valid; negative on invalid input.
        real(dp), intent(in), optional :: rtol !! Relative tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Iteration or restart-cycle cap.
        real(dp), intent(in), optional :: x0(:) !! Optional finite starting guess.
        logical, intent(in), optional :: jacobi !! Apply inverse diagonal, default false.
        real(dp), allocatable, intent(out) :: invdiag(:) !! Multiplicative preconditioner.
        real(dp), intent(out) :: threshold !! Max of absolute and relative residual targets.
        integer, intent(out) :: cap !! Iteration cap, default max(1,10*n).
        type(ilu0_factor), intent(in), optional :: ilu !! Validated ILU(0) factor, exclusive with Jacobi.
        integer :: n, i, k
        real(dp) :: relative, absolute, diagonal
        logical :: precond

        n = size(b)
        allocate(x(n), invdiag(n))
        x = 0.0_dp
        invdiag = 1.0_dp
        threshold = huge(1.0_dp)
        cap = 0
        info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        if (a%n_rows /= a%n_cols .or. a%n_rows /= n) return
        if (.not. all(ieee_is_finite(b))) return
        relative = 1.0e-5_dp
        absolute = 0.0_dp
        if (present(rtol)) relative = rtol
        if (present(atol)) absolute = atol
        if (.not. ieee_is_finite(relative) .or. .not. ieee_is_finite(absolute)) return
        if (relative < 0.0_dp .or. absolute < 0.0_dp) return
        if (n > huge(n)/10) then
            cap = huge(n)
        else
            cap = max(1, 10*n)
        end if
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        if (present(x0)) then
            if (size(x0) /= n) return
            if (.not. all(ieee_is_finite(x0))) return
            x = x0
        end if
        threshold = max(absolute, relative*norm2(b))
        if (.not. ieee_is_finite(threshold)) return
        precond = .false.
        if (present(jacobi)) precond = jacobi
        if (present(ilu)) then
            if (precond .or. .not. valid_ilu0(ilu)) return
            if (ilu%lu%n_rows /= n) return
        end if
        if (precond) then
            do i = 1, n
                diagonal = 0.0_dp
                do k = a%indptr(i), a%indptr(i+1)-1
                    if (a%indices(k) == i) then
                        diagonal = a%data(k)
                        exit
                    end if
                end do
                if (diagonal == 0.0_dp) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
                invdiag(i) = 1.0_dp/diagonal
                if (.not. ieee_is_finite(invdiag(i))) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
            end do
        end if
        info = sparse_linalg_success
    end subroutine prepare

    subroutine cg(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm)
        type(csr_matrix), intent(in) :: a !! Real symmetric positive definite CSR matrix.
        real(dp), intent(in) :: b(:) !! RHS vector of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Approximate solution, length n.
        integer, intent(out) :: info !! 0 success, positive nonconvergence, negative invalid/breakdown.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Iteration cap, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Starting guess, default zero.
        logical, intent(in), optional :: jacobi !! Enable diagonal inverse preconditioner.
        integer, intent(out), optional :: iterations !! Number of performed CG steps.
        real(dp), intent(out), optional :: residual_norm !! Norm of actual final residual.
        real(dp), allocatable :: invdiag(:), r(:), z(:), p(:), ap(:)
        real(dp) :: threshold, rho, rho_new, denom, alpha, beta, err
        integer :: n, cap, i, steps

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap)
        if (info /= sparse_linalg_success) return
        n = size(b)
        allocate(r(n),z(n),p(n),ap(n))
        call csr_product(a,x,ap)
        r = b-ap
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            z = invdiag*r
            p = z
            rho = dot_product(r,z)
            if (rho <= 0.0_dp .or. .not. ieee_is_finite(rho)) info = sparse_linalg_breakdown
            if (info /= sparse_linalg_breakdown) then
                do i = 1, cap
                    steps = i
                    call csr_product(a,p,ap)
                    denom = dot_product(p,ap)
                    if (denom <= 0.0_dp .or. .not. ieee_is_finite(denom)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    alpha = rho/denom
                    if (.not. ieee_is_finite(alpha)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    x = x+alpha*p
                    r = r-alpha*ap
                    err = norm2(r)
                    if (.not. ieee_is_finite(err) .or. .not. all(ieee_is_finite(x))) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    if (err <= threshold) then
                        call csr_product(a,x,ap)
                        r = b-ap
                        err = norm2(r)
                        if (err <= threshold) then
                            info = 0
                            exit
                        end if
                    end if
                    z = invdiag*r
                    rho_new = dot_product(r,z)
                    if (rho_new <= 0.0_dp .or. .not. ieee_is_finite(rho_new)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    beta = rho_new/rho
                    if (.not. ieee_is_finite(beta)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    p = z+beta*p
                    rho = rho_new
                end do
            end if
        end if
        call csr_product(a,x,ap)
        err = norm2(b-ap)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine cg

    ! Real BiConjugate Gradient for nonsymmetric square systems.  This
    ! method explicitly uses both A and A**T; the shadow preconditioner is
    ! the transpose inverse of the diagonal Jacobi preconditioner.
    ! Reference: Barrett et al., Templates for the Solution of Linear Systems,
    ! SIAM (1994), sections on BiCG and left preconditioning.
    subroutine bicg(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm)
        type(csr_matrix), intent(in) :: a !! Real, canonical square CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Last iterate, length n.
        integer, intent(out) :: info !! 0 converged, positive nonconvergence, negative failure.
        real(dp), intent(in), optional :: rtol !! Relative true residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true residual tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Maximum iterations, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Finite starting guess, default zero.
        logical, intent(in), optional :: jacobi !! Enable diagonal preconditioning, default false.
        integer, intent(out), optional :: iterations !! Number of performed BiCG iterations.
        real(dp), intent(out), optional :: residual_norm !! Euclidean norm of the true residual.
        type(csr_matrix) :: transpose_a
        real(dp), allocatable :: invdiag(:), r(:), shadow(:), z(:), zshadow(:)
        real(dp), allocatable :: p(:), pshadow(:), ap(:), atp(:), ax(:)
        real(dp) :: threshold, err, rho, rho_old, alpha, beta, denominator
        integer :: n, cap, step, steps, ierr
        logical :: start_recurrence

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap)
        if (info /= sparse_linalg_success) return
        transpose_a = csr_transpose(a,ierr)
        if (ierr /= sparse_success) then
            info = sparse_linalg_invalid_input
            return
        end if
        n = size(b)
        allocate(r(n),shadow(n),z(n),zshadow(n),p(n),pshadow(n),ap(n),atp(n),ax(n))
        call csr_product(a,x,ax)
        r = b-ax
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_linalg_success
        else
            info = max(1,cap)
            shadow = r
            start_recurrence = .true.
            rho_old = 0.0_dp
            do step = 1, cap
                steps = step
                z = invdiag*r
                zshadow = invdiag*shadow
                rho = dot_product(shadow,z)
                if (rho == 0.0_dp .or. .not. ieee_is_finite(rho)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (start_recurrence) then
                    p = z
                    pshadow = zshadow
                    start_recurrence = .false.
                else
                    beta = rho/rho_old
                    if (.not. ieee_is_finite(beta)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    p = z+beta*p
                    pshadow = zshadow+beta*pshadow
                end if
                if (.not. all(ieee_is_finite(p)) .or. &
                    .not. all(ieee_is_finite(pshadow))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call csr_product(a,p,ap)
                call csr_product(transpose_a,pshadow,atp)
                denominator = dot_product(pshadow,ap)
                if (denominator == 0.0_dp .or. .not. ieee_is_finite(denominator)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                alpha = rho/denominator
                if (.not. ieee_is_finite(alpha)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                x = x+alpha*p
                r = r-alpha*ap
                shadow = shadow-alpha*atp
                if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(r)) .or. &
                    .not. all(ieee_is_finite(shadow))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                ! Verify convergence against the actual (unpreconditioned)
                ! residual; recursive residuals can severely underestimate it.
                call csr_product(a,x,ax)
                err = norm2(b-ax)
                if (.not. ieee_is_finite(err)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    info = sparse_linalg_success
                    exit
                end if
                if (norm2(r) <= threshold) then
                    ! Reliable restart after an optimistic recursive residual.
                    r = b-ax
                    shadow = r
                    start_recurrence = .true.
                end if
                rho_old = rho
            end do
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_linalg_success
        end if
    end subroutine bicg

    ! Conjugate Gradient Squared uses two forward A products and no transpose.
    ! It can have irregular residual behavior; only a checked true residual
    ! (not its squared Krylov polynomial recurrence) establishes convergence.
    ! Reference: P. Sonneveld (1989); Barrett et al., Templates (1994).
    subroutine cgs(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm, ilu)
        type(csr_matrix), intent(in) :: a !! Real, canonical square CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Last iterate, length n.
        integer, intent(out) :: info !! 0 converged, positive nonconvergence, negative failure.
        real(dp), intent(in), optional :: rtol !! Relative true residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true residual tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Maximum iterations, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Finite starting guess, default zero.
        logical, intent(in), optional :: jacobi !! Enable diagonal preconditioning, default false.
        integer, intent(out), optional :: iterations !! Number of performed CGS iterations.
        real(dp), intent(out), optional :: residual_norm !! Euclidean norm of the true residual.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0), exclusive with jacobi.
        real(dp), allocatable :: invdiag(:), r(:), shadow(:), p(:), q(:), u(:)
        real(dp), allocatable :: phat(:), v(:), uhat(:), ax(:)
        real(dp) :: threshold, err, rho, rho_old, beta, alpha, denominator
        integer :: n, cap, step, steps
        logical :: start_recurrence

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap,ilu)
        if (info /= sparse_linalg_success) return
        n = size(b)
        allocate(r(n),shadow(n),p(n),q(n),u(n),phat(n),v(n),uhat(n),ax(n))
        call csr_product(a,x,ax)
        r = b-ax
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_linalg_success
        else
            info = max(1,cap)
            shadow = r
            start_recurrence = .true.
            rho_old = 0.0_dp
            do step = 1, cap
                steps = step
                rho = dot_product(shadow,r)
                if (rho == 0.0_dp .or. .not. ieee_is_finite(rho)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (start_recurrence) then
                    u = r
                    p = u
                    start_recurrence = .false.
                else
                    beta = rho/rho_old
                    if (.not. ieee_is_finite(beta)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    u = r+beta*q
                    p = u+beta*(q+beta*p)
                end if
                if (.not. all(ieee_is_finite(p))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call apply_square_preconditioner(invdiag,p,phat,ilu)
                if (.not. all(ieee_is_finite(phat))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call csr_product(a,phat,v)
                denominator = dot_product(shadow,v)
                if (denominator == 0.0_dp .or. .not. ieee_is_finite(denominator)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                alpha = rho/denominator
                if (.not. ieee_is_finite(alpha)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                q = u-alpha*v
                call apply_square_preconditioner(invdiag,u+q,uhat,ilu)
                if (.not. all(ieee_is_finite(uhat))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                x = x+alpha*uhat
                call csr_product(a,uhat,v)
                r = r-alpha*v
                if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(r))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call csr_product(a,x,ax)
                err = norm2(b-ax)
                if (.not. ieee_is_finite(err)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    info = sparse_linalg_success
                    exit
                end if
                if (norm2(r) <= threshold) then
                    r = b-ax
                    shadow = r
                    start_recurrence = .true.
                end if
                rho_old = rho
            end do
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_linalg_success
        end if
    end subroutine cgs

    ! Quasi-Minimal Residual iteration for real nonsymmetric sparse systems.
    ! Two-sided Lanczos recurrence: A and A**T products are both required.
    ! The Jacobi option is applied on the right, A*D**(-1), with its transpose
    ! acting on the shadow space. Algorithm adapted from SciPy 1.17 qmr,
    ! scipy/sparse/linalg/_isolve/iterative.py (BSD-3-Clause): see
    ! LICENSES/scipy-qmr-BSD-3-Clause.txt and CODE_PROVENANCE.md.
    subroutine qmr(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm)
        type(csr_matrix), intent(in) :: a !! Real, canonical square CSR matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side vector, length n.
        real(dp), allocatable, intent(out) :: x(:) !! Last approximation, length n.
        integer, intent(out) :: info !! 0 converged; positive cap reached; negative failure.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default zero.
        integer, intent(in), optional :: maxiter !! Maximum QMR iterations, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Initial approximation, default zero.
        logical, intent(in), optional :: jacobi !! Right Jacobi inverse, default false.
        integer, intent(out), optional :: iterations !! Number of executed QMR steps.
        real(dp), intent(out), optional :: residual_norm !! Norm of true final residual ||b-A*x||.
        type(csr_matrix) :: at
        real(dp), allocatable :: invdiag(:), r(:), vt(:), wt(:), y(:), z(:)
        real(dp), allocatable :: v(:), w(:), p(:), q(:), pt(:), yt(:), zt(:)
        real(dp), allocatable :: d(:), s(:), ax(:)
        real(dp) :: threshold, err, rho, rho_prev, xi, delta, epsilon_prev
        real(dp) :: epsilon_now, beta, theta, theta_prev, gamma, gamma_prev, eta, weight
        integer :: n, cap, i, steps, ierr

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap)
        if (info /= sparse_linalg_success) return
        n = size(b)
        allocate(r(n),vt(n),wt(n),y(n),z(n),v(n),w(n),p(n),q(n))
        allocate(pt(n),yt(n),zt(n),d(n),s(n),ax(n))
        call csr_product(a,x,ax)
        r = b-ax
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_linalg_success
        else
            at = csr_transpose(a,ierr)
            if (ierr /= sparse_success) then
                info = sparse_linalg_invalid_input
                return
            end if
            info = max(1,cap)
            vt = r
            wt = r
            y = vt
            z = invdiag*wt
            rho = norm2(y)
            xi = norm2(z)
            gamma = 1.0_dp
            eta = -1.0_dp
            theta = 0.0_dp
            epsilon_prev = 0.0_dp
            d = 0.0_dp
            s = 0.0_dp
            do i = 1, cap
                steps = i
                if (rho == 0.0_dp .or. xi == 0.0_dp .or. &
                    .not. ieee_is_finite(rho) .or. .not. ieee_is_finite(xi)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                v = vt/rho
                y = y/rho
                w = wt/xi
                z = z/xi
                delta = dot_product(z,y)
                if (delta == 0.0_dp .or. .not. ieee_is_finite(delta)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                yt = invdiag*y
                zt = z
                if (i > 1) then
                    if (epsilon_prev == 0.0_dp) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    yt = yt-(xi*delta/epsilon_prev)*p
                    zt = zt-(rho*delta/epsilon_prev)*q
                end if
                p = yt
                q = zt
                if (.not. all(ieee_is_finite(p)) .or. .not. all(ieee_is_finite(q))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call csr_product(a,p,pt)
                epsilon_now = dot_product(q,pt)
                if (epsilon_now == 0.0_dp .or. .not. ieee_is_finite(epsilon_now)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                beta = epsilon_now/delta
                if (beta == 0.0_dp .or. .not. ieee_is_finite(beta)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                vt = pt-beta*v
                y = vt
                rho_prev = rho
                rho = norm2(y)
                call csr_product(at,q,wt)
                wt = wt-beta*w
                z = invdiag*wt
                xi = norm2(z)
                gamma_prev = gamma
                theta_prev = theta
                theta = rho/(gamma_prev*abs(beta))
                gamma = 1.0_dp/hypot(1.0_dp,theta)
                if (gamma == 0.0_dp .or. .not. ieee_is_finite(gamma)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                eta = -eta*(rho_prev/beta)*(gamma/gamma_prev)**2
                weight = (theta_prev*gamma)**2
                if (.not. ieee_is_finite(eta) .or. .not. ieee_is_finite(weight)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                d = eta*p+weight*d
                s = eta*pt+weight*s
                x = x+d
                if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(s))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                ! Always use the actual (not quasi-minimal-estimated) residual.
                call csr_product(a,x,ax)
                err = norm2(b-ax)
                if (.not. ieee_is_finite(err)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    info = sparse_linalg_success
                    exit
                end if
                epsilon_prev = epsilon_now
            end do
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_linalg_success
        end if
    end subroutine qmr

    ! Transpose-Free Quasi-Minimal Residual iteration (Freund 1993).
    ! This adapts the BSD-3-Clause SciPy TFQMR recurrence via SciFort's
    ! existing BSD-3-Clause matrix-free implementation in
    ! src/optimize/krylov_iterative.f90. The A*M right-preconditioned variant
    ! is mathematically consistent. One iteration means a half-step.
    ! See LICENSES/scipy-tfqmr-BSD-3-Clause.txt.
    subroutine tfqmr(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm, ilu)
        type(csr_matrix), intent(in) :: a !! Real, canonical square CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side, length n.
        real(dp), allocatable, intent(out) :: x(:) !! Last approximation, length n.
        integer, intent(out) :: info !! 0 converged; positive cap reached; negative failure.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default zero.
        integer, intent(in), optional :: maxiter !! Maximum half-iterations, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Initial approximation, default zero.
        logical, intent(in), optional :: jacobi !! Right Jacobi inverse, default false.
        integer, intent(out), optional :: iterations !! Number of executed TFQMR half-steps.
        real(dp), intent(out), optional :: residual_norm !! Norm of true final residual ||b-A*x||.
        type(ilu0_factor), intent(in), optional :: ilu !! Right ILU(0), exclusive with Jacobi.
        real(dp), allocatable :: invdiag(:), r(:), shadow(:), u(:), w(:), v(:)
        real(dp), allocatable :: uhat(:), unext(:), d(:), z(:), ax(:)
        real(dp) :: threshold, err, rho, rho_prev, tau, theta, eta
        real(dp) :: alpha, beta, c, denominator, correction
        integer :: n, cap, i, steps
        logical :: odd

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap,ilu)
        if (info /= sparse_linalg_success) return
        n = size(b)
        allocate(r(n),shadow(n),u(n),w(n),v(n),uhat(n),unext(n),d(n),z(n),ax(n))
        call csr_product(a,x,ax)
        r = b-ax
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_linalg_success
        else
            info = max(1,cap)
            shadow = r
            u = r
            w = r
            d = 0.0_dp
            call apply_square_preconditioner(invdiag,u,z,ilu)
            call csr_product(a,z,v)
            uhat = v
            rho = dot_product(shadow,r)
            rho_prev = rho
            tau = err
            theta = 0.0_dp
            eta = 0.0_dp
            alpha = 0.0_dp
            odd = .true.
            if (.not. all(ieee_is_finite(v))) then
                info = sparse_linalg_breakdown
            else
                do i = 1, cap
                    steps = i
                    if (odd) then
                        denominator = dot_product(shadow,v)
                        if (denominator == 0.0_dp .or. .not. ieee_is_finite(denominator)) then
                            info = sparse_linalg_breakdown
                            exit
                        end if
                        alpha = rho/denominator
                        if (alpha == 0.0_dp .or. .not. ieee_is_finite(alpha)) then
                            info = sparse_linalg_breakdown
                            exit
                        end if
                        unext = u-alpha*v
                    end if
                    w = w-alpha*uhat
                    correction = (theta*theta*eta)/alpha
                    d = u+correction*d
                    if (.not. ieee_is_finite(correction) .or. .not. all(ieee_is_finite(d))) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    if (tau <= 0.0_dp .or. .not. ieee_is_finite(tau)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    theta = norm2(w)/tau
                    c = 1.0_dp/hypot(1.0_dp,theta)
                    tau = tau*theta*c
                    eta = c*c*alpha
                    if (.not. ieee_is_finite(eta) .or. .not. ieee_is_finite(tau)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    call apply_square_preconditioner(invdiag,d,z,ilu)
                    x = x+eta*z
                    if (.not. all(ieee_is_finite(x))) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    call csr_product(a,x,ax)
                    err = norm2(b-ax)
                    if (.not. ieee_is_finite(err)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    if (err <= threshold) then
                        info = sparse_linalg_success
                        exit
                    end if
                    if (tau <= tiny(1.0_dp)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    if (odd) then
                        call apply_square_preconditioner(invdiag,unext,z,ilu)
                        call csr_product(a,z,uhat)
                        u = unext
                        rho_prev = rho
                    else
                        if (rho_prev == 0.0_dp) then
                            info = sparse_linalg_breakdown
                            exit
                        end if
                        rho = dot_product(shadow,w)
                        beta = rho/rho_prev
                        if (.not. ieee_is_finite(beta)) then
                            info = sparse_linalg_breakdown
                            exit
                        end if
                        u = w+beta*u
                        v = beta*uhat+beta*beta*v
                        call apply_square_preconditioner(invdiag,u,z,ilu)
                        call csr_product(a,z,uhat)
                        v = v+uhat
                    end if
                    if (.not. all(ieee_is_finite(v)) .or. &
                        .not. all(ieee_is_finite(uhat)) .or. &
                        .not. all(ieee_is_finite(u))) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    odd = .not. odd
                end do
            end if
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_linalg_success
        end if
    end subroutine tfqmr

    subroutine bicgstab(a, b, x, info, rtol, atol, maxiter, x0, jacobi, iterations, residual_norm, ilu)
        type(csr_matrix), intent(in) :: a !! Real square CSR matrix (not necessarily symmetric).
        real(dp), intent(in) :: b(:) !! RHS vector of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Approximate solution, length n.
        integer, intent(out) :: info !! 0 success, positive nonconvergence, negative invalid/breakdown.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Iteration cap, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Starting guess, default zero.
        logical, intent(in), optional :: jacobi !! Enable diagonal inverse preconditioner.
        integer, intent(out), optional :: iterations !! Number of BiCGSTAB steps.
        real(dp), intent(out), optional :: residual_norm !! Norm of actual final residual.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0), exclusive with jacobi.
        real(dp), allocatable :: invdiag(:), r(:), rhat(:), p(:), v(:), ph(:), s(:), sh(:), t(:)
        real(dp) :: threshold, err, rho, rho_old, alpha, beta, omega, rv, tt
        integer :: n, cap, i, steps
        logical :: restart_recurrence

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap,ilu)
        if (info /= 0) return
        n = size(b)
        allocate(r(n),rhat(n),p(n),v(n),ph(n),s(n),sh(n),t(n))
        call csr_product(a,x,t)
        r = b-t
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            rhat = r
            p = 0.0_dp
            v = 0.0_dp
            rho_old = 1.0_dp
            alpha = 1.0_dp
            omega = 1.0_dp
            restart_recurrence = .true.
            do i = 1, cap
                steps = i
                rho = dot_product(rhat,r)
                if (rho == 0.0_dp .or. .not. ieee_is_finite(rho)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (restart_recurrence) then
                    p = r
                    restart_recurrence = .false.
                else
                    if (omega == 0.0_dp) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    beta = (rho/rho_old)*(alpha/omega)
                    p = r+beta*(p-omega*v)
                end if
                call apply_square_preconditioner(invdiag,p,ph,ilu)
                if (.not. all(ieee_is_finite(ph))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call csr_product(a,ph,v)
                rv = dot_product(rhat,v)
                if (rv == 0.0_dp .or. .not. ieee_is_finite(rv)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                alpha = rho/rv
                s = r-alpha*v
                if (.not. ieee_is_finite(alpha) .or. .not. all(ieee_is_finite(s))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (norm2(s) <= threshold) then
                    x = x+alpha*ph
                    call csr_product(a,x,t)
                    r = b-t
                    err = norm2(r)
                    if (err <= threshold) then
                        info = 0
                        exit
                    end if
                    ! Recursive residual was optimistic; restart from true residual.
                    rhat = r
                    p = 0.0_dp
                    v = 0.0_dp
                    rho_old = 1.0_dp
                    alpha = 1.0_dp
                    omega = 1.0_dp
                    restart_recurrence = .true.
                    cycle
                end if
                call apply_square_preconditioner(invdiag,s,sh,ilu)
                if (.not. all(ieee_is_finite(sh))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                call csr_product(a,sh,t)
                tt = dot_product(t,t)
                if (tt == 0.0_dp .or. .not. ieee_is_finite(tt)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                omega = dot_product(t,s)/tt
                if (omega == 0.0_dp .or. .not. ieee_is_finite(omega)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                x = x+alpha*ph+omega*sh
                r = s-omega*t
                err = norm2(r)
                if (.not. ieee_is_finite(err) .or. .not. all(ieee_is_finite(x))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (err <= threshold) then
                    call csr_product(a,x,t)
                    r = b-t
                    err = norm2(r)
                    if (err <= threshold) then
                        info = 0
                        exit
                    end if
                end if
                rho_old = rho
            end do
        end if
        call csr_product(a,x,t)
        err = norm2(b-t)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine bicgstab

    subroutine gmres(a, b, x, info, rtol, atol, maxiter, restart, x0, jacobi, &
            iterations, residual_norm, ilu)
        type(csr_matrix), intent(in) :: a !! Real square CSR matrix.
        real(dp), intent(in) :: b(:) !! RHS vector of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Approximate solution, length n.
        integer, intent(out) :: info !! 0 success, positive cycles without convergence, negative error.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Restart-cycle cap, default max(1,10*n).
        integer, intent(in), optional :: restart !! Arnoldi dimension, default min(20,n).
        real(dp), intent(in), optional :: x0(:) !! Starting guess, default zero.
        logical, intent(in), optional :: jacobi !! Apply left diagonal preconditioning.
        integer, intent(out), optional :: iterations !! Total Arnoldi iterations.
        real(dp), intent(out), optional :: residual_norm !! Norm of actual final residual.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0), exclusive with jacobi.
        real(dp), allocatable :: invdiag(:), r(:), w(:), ax(:), z(:), xtry(:), v(:,:), h(:,:)
        real(dp), allocatable :: c(:), s(:), g(:), y(:)
        real(dp) :: threshold, err, beta, denom, temp, hij
        integer :: cap, n, m, i, j, k, cycle, completed, steps
        logical :: ended, singular

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap,ilu)
        if (info /= 0) return
        n = size(b)
        m = min(20,n)
        if (present(restart)) m = restart
        if (m <= 0 .and. n > 0) then
            info = sparse_linalg_invalid_input
            return
        end if
        m = min(m,n)
        allocate(r(n),w(n),ax(n),z(n),xtry(n),v(n,m+1),h(m+1,m),c(m),s(m),g(m+1),y(m))
        call csr_product(a,x,ax)
        r = b-ax
        err = norm2(r)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            do cycle = 1, cap
                r = b-ax
                call apply_square_preconditioner(invdiag,r,z,ilu)
                beta = norm2(z)
                if (beta == 0.0_dp .or. .not. ieee_is_finite(beta)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                v(:,1) = z/beta
                h = 0.0_dp
                c = 0.0_dp
                s = 0.0_dp
                g = 0.0_dp
                g(1) = beta
                completed = 0
                ended = .false.
                do j = 1, m
                    steps = steps+1
                    completed = j
                    call csr_product(a,v(:,j),w)
                    call apply_square_preconditioner(invdiag,w,z,ilu)
                    w = z
                    if (.not. all(ieee_is_finite(w))) then
                        info = sparse_linalg_breakdown
                        ended = .true.
                        exit
                    end if
                    ! Twice modified Gram-Schmidt for improved orthogonality.
                    do k = 1, 2
                        do i = 1, j
                            hij = dot_product(v(:,i),w)
                            h(i,j) = h(i,j)+hij
                            w = w-hij*v(:,i)
                        end do
                    end do
                    h(j+1,j) = norm2(w)
                    if (h(j+1,j) > 0.0_dp) v(:,j+1) = w/h(j+1,j)
                    do i = 1, j-1
                        temp = c(i)*h(i,j)+s(i)*h(i+1,j)
                        h(i+1,j) = -s(i)*h(i,j)+c(i)*h(i+1,j)
                        h(i,j) = temp
                    end do
                    denom = norm2(h(j:j+1,j))
                    if (denom == 0.0_dp .or. .not. ieee_is_finite(denom)) then
                        info = sparse_linalg_breakdown
                        ended = .true.
                        exit
                    end if
                    c(j) = h(j,j)/denom
                    s(j) = h(j+1,j)/denom
                    h(j,j) = denom
                    h(j+1,j) = 0.0_dp
                    g(j+1) = -s(j)*g(j)
                    g(j) = c(j)*g(j)
                    if (.not. all(ieee_is_finite(g(:j+1)))) then
                        info = sparse_linalg_breakdown
                        ended = .true.
                        exit
                    end if
                    ! Form and verify the candidate with the *true* residual;
                    ! the rotated preconditioned norm alone is not sufficient.
                    y(:j) = g(:j)
                    singular = .false.
                    do i = j, 1, -1
                        if (h(i,i) == 0.0_dp) then
                            singular = .true.
                            exit
                        end if
                        y(i) = (y(i)-dot_product(h(i,i+1:j),y(i+1:j)))/h(i,i)
                    end do
                    if (singular .or. .not. all(ieee_is_finite(y(:j)))) then
                        info = sparse_linalg_breakdown
                        ended = .true.
                        exit
                    end if
                    xtry = x+matmul(v(:,:j),y(:j))
                    call csr_product(a,xtry,ax)
                    err = norm2(b-ax)
                    if (.not. ieee_is_finite(err) .or. .not. all(ieee_is_finite(xtry))) then
                        info = sparse_linalg_breakdown
                        ended = .true.
                        exit
                    end if
                    if (err <= threshold) then
                        x = xtry
                        info = 0
                        ended = .true.
                        exit
                    end if
                    if (norm2(w) == 0.0_dp) then
                        ! Arnoldi invariant subspace: no new direction exists.
                        info = sparse_linalg_breakdown
                        ended = .true.
                        exit
                    end if
                end do
                if (ended) exit
                if (completed > 0) x = xtry
                call csr_product(a,x,ax)
                err = norm2(b-ax)
                if (err <= threshold) then
                    info = 0
                    exit
                end if
            end do
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine gmres

    ! Solve the unpreconditioned, right-flexible Arnoldi least-squares problem
    ! in the projected space orthogonal to recycle_c. The first inner_m
    ! directions are M^{-1} v_j; the remaining directions are retained LGMRES
    ! corrections. Unlike left-preconditioned GMRES, all estimates here refer
    ! to the original Euclidean residual. Two MGS passes stabilize Arnoldi.
    subroutine recycled_arnoldi_cycle(a, r, invdiag, inner_m, outer_v, recycle_c, &
            target, dx, projections, steps, info, ilu)
        type(csr_matrix), intent(in) :: a !! Prevalidated canonical square CSR operator.
        real(dp), intent(in) :: r(:) !! Starting residual orthogonal to recycle_c, length n.
        real(dp), intent(in) :: invdiag(:) !! Jacobi inverse diagonal or identity, length n.
        integer, intent(in) :: inner_m !! Maximum newly generated Krylov directions, at least one.
        real(dp), intent(in) :: outer_v(:,:) !! Optional normalized LGMRES corrections, n by p.
        real(dp), intent(in) :: recycle_c(:,:) !! Orthonormal GCROT image basis, n by k.
        real(dp), intent(in) :: target !! Target norm of the original residual, nonnegative.
        real(dp), intent(out) :: dx(:) !! Right-preconditioned correction, length n.
        real(dp), intent(out) :: projections(:) !! recycle_c**T A dx coefficients, length k.
        integer, intent(out) :: steps !! Number of newly computed Arnoldi columns.
        integer, intent(out) :: info !! Zero for a valid correction, negative on breakdown.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional validated right ILU(0).
        real(dp), allocatable :: v(:,:), z(:,:), h(:,:), bproj(:,:), &
            cs(:), sn(:), g(:), y(:), w(:)
        real(dp) :: beta, hij, wnorm, hnorm, denom, temp, scale
        integer :: n, nrec, nouter, ninner, ncols, first_outer, j, i, pass, done
        logical :: dependent

        n = size(r)
        nrec = size(recycle_c, 2)
        dx = 0.0_dp
        projections = 0.0_dp
        steps = 0
        info = sparse_linalg_breakdown
        beta = norm2(r)
        if (beta == 0.0_dp .or. .not. ieee_is_finite(beta)) return
        ! Do not produce more columns than the dimension of C's complement.
        nouter = min(size(outer_v, 2), max(0, n-nrec-1))
        ninner = min(inner_m, max(0, n-nrec-nouter))
        if (ninner <= 0) return
        ncols = ninner+nouter
        first_outer = size(outer_v, 2)-nouter
        allocate(v(n,ncols+1), z(n,ncols), h(ncols+1,ncols), &
            bproj(nrec,ncols), cs(ncols), sn(ncols), g(ncols+1), &
            y(ncols), w(n))
        v = 0.0_dp
        v(:,1) = r/beta
        h = 0.0_dp
        bproj = 0.0_dp
        cs = 0.0_dp
        sn = 0.0_dp
        g = 0.0_dp
        g(1) = beta
        done = 0
        do j = 1, ncols
            if (j <= ninner) then
                call apply_square_preconditioner(invdiag,v(:,j),z(:,j),ilu)
            else
                z(:,j) = outer_v(:,first_outer+j-ninner)
            end if
            if (.not. all(ieee_is_finite(z(:,j)))) return
            call csr_product(a,z(:,j),w)
            if (.not. all(ieee_is_finite(w))) return
            wnorm = norm2(w)
            ! Enforce C^T w = 0. Sum both passes into the projection matrix.
            do pass = 1, 2
                do i = 1, nrec
                    hij = dot_product(recycle_c(:,i),w)
                    bproj(i,j) = bproj(i,j)+hij
                    w = w-hij*recycle_c(:,i)
                end do
            end do
            do pass = 1, 2
                do i = 1, j
                    hij = dot_product(v(:,i),w)
                    h(i,j) = h(i,j)+hij
                    w = w-hij*v(:,i)
                end do
            end do
            hnorm = norm2(w)
            if (.not. ieee_is_finite(hnorm)) return
            scale = max(wnorm, tiny(1.0_dp))
            dependent = hnorm <= 32.0_dp*epsilon(1.0_dp)*scale
            if (.not. dependent) then
                h(j+1,j) = hnorm
                v(:,j+1) = w/hnorm
            else
                h(j+1,j) = 0.0_dp
            end if
            do i = 1, j-1
                temp = cs(i)*h(i,j)+sn(i)*h(i+1,j)
                h(i+1,j) = -sn(i)*h(i,j)+cs(i)*h(i+1,j)
                h(i,j) = temp
            end do
            denom = norm2(h(j:j+1,j))
            steps = steps+1
            if (denom == 0.0_dp .or. .not. ieee_is_finite(denom)) return
            cs(j) = h(j,j)/denom
            sn(j) = h(j+1,j)/denom
            h(j,j) = denom
            h(j+1,j) = 0.0_dp
            g(j+1) = -sn(j)*g(j)
            g(j) = cs(j)*g(j)
            if (.not. all(ieee_is_finite(g(:j+1)))) return
            done = j
            if (abs(g(j+1)) <= target .or. dependent) exit
        end do
        if (done == 0) return
        y(:done) = g(:done)
        do i = done, 1, -1
            if (abs(h(i,i)) <= tiny(1.0_dp)) return
            y(i) = (y(i)-dot_product(h(i,i+1:done),y(i+1:done)))/h(i,i)
        end do
        if (.not. all(ieee_is_finite(y(:done)))) return
        dx = matmul(z(:,:done),y(:done))
        projections = matmul(bproj(:,:done),y(:done))
        if (.not. all(ieee_is_finite(dx)) .or. &
                .not. all(ieee_is_finite(projections))) return
        if (norm2(dx) == 0.0_dp) return
        info = sparse_linalg_success
    end subroutine recycled_arnoldi_cycle

    ! Baker--Jessup--Manteuffel (2005) LGMRES: augment each restarted
    ! flexible-GMRES cycle with normalized corrections from previous cycles.
    ! If supplied, outer_v transports corrections across related systems;
    ! image products are always recomputed, so changing A is supported.
    subroutine lgmres(a, b, x, info, rtol, atol, maxiter, inner_m, outer_k, &
            x0, jacobi, iterations, residual_norm, ilu, outer_v)
        type(csr_matrix), intent(in) :: a !! Square canonical CSR matrix, n by n.
        real(dp), intent(in) :: b(:) !! Finite RHS vector of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Approximate solution, length n.
        integer, intent(out) :: info !! 0 converged, >0 cycle cap, <0 invalid/breakdown.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default zero.
        integer, intent(in), optional :: maxiter !! Maximum outer cycles, default max(1,10*n).
        integer, intent(in), optional :: inner_m !! Krylov columns per cycle, default min(30,n).
        integer, intent(in), optional :: outer_k !! Maximum retained corrections, default three.
        real(dp), intent(in), optional :: x0(:) !! Optional starting solution, length n.
        logical, intent(in), optional :: jacobi !! Right diagonal preconditioner, default false.
        integer, intent(out), optional :: iterations !! Total Arnoldi column evaluations.
        real(dp), intent(out), optional :: residual_norm !! Norm2(b-A*x), unpreconditioned.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional right ILU(0), exclusive with Jacobi.
        real(dp), allocatable, intent(inout), optional :: outer_v(:,:) !! Reused corrections, n by <=outer_k.
        real(dp), allocatable :: invdiag(:), previous(:,:), empty(:,:), empty_coeff(:)
        real(dp), allocatable :: r(:), ax(:), dx(:), test_ax(:)
        real(dp) :: threshold, err, step_norm, corr
        integer :: n, m, kmax, cap, count, cyc, j, start, used, status, total_steps
        logical :: duplicate

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap,ilu)
        if (info /= 0) return
        n = size(b)
        m = min(30,n)
        if (present(inner_m)) m = inner_m
        kmax = min(3,max(0,n-1))
        if (present(outer_k)) kmax = outer_k
        if ((m <= 0 .and. n > 0) .or. kmax < 0) then
            info = sparse_linalg_invalid_input
            return
        end if
        m = min(m,n)
        kmax = min(kmax,max(0,n-1))
        if (present(outer_v)) then
            if (allocated(outer_v)) then
                if (size(outer_v,1) /= n) then
                    info = sparse_linalg_invalid_input
                    return
                end if
                if (.not. all(ieee_is_finite(outer_v))) then
                    info = sparse_linalg_invalid_input
                    return
                end if
            end if
        end if
        allocate(previous(n,kmax),empty(n,0),empty_coeff(0),r(n),ax(n),dx(n),test_ax(n))
        count = 0
        if (present(outer_v)) then
            if (allocated(outer_v)) then
                start = max(1,size(outer_v,2)-kmax+1)
                do j = start, size(outer_v,2)
                    step_norm = norm2(outer_v(:,j))
                    if (step_norm <= tiny(1.0_dp) .or. &
                            .not. ieee_is_finite(step_norm)) cycle
                    count = count+1
                    previous(:,count) = outer_v(:,j)/step_norm
                end do
            end if
        end if
        total_steps = 0
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            do cyc = 1, cap
                r = b-ax
                call recycled_arnoldi_cycle(a,r,invdiag,m,previous(:,:count), &
                    empty,threshold,dx,empty_coeff,used,status,ilu)
                total_steps = total_steps+used
                if (status /= 0) then
                    info = status
                    exit
                end if
                call csr_product(a,x+dx,test_ax)
                err = norm2(b-test_ax)
                if (.not. ieee_is_finite(err)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                ! Save the full-cycle correction (not individual Arnoldi vectors).
                step_norm = norm2(dx)
                if (step_norm > tiny(1.0_dp) .and. ieee_is_finite(step_norm) &
                        .and. kmax > 0) then
                    duplicate = .false.
                    do j = 1, count
                        corr = dot_product(dx/step_norm,previous(:,j))
                        if (abs(corr) >= 1.0_dp-1.0e-10_dp) duplicate = .true.
                    end do
                    if (.not. duplicate) then
                        if (count == kmax) then
                            if (kmax > 1) previous(:,1:kmax-1) = previous(:,2:kmax)
                            count = count-1
                        end if
                        count = count+1
                        previous(:,count) = dx/step_norm
                    end if
                end if
                x = x+dx
                ax = test_ax
                if (err <= threshold) then
                    info = 0
                    exit
                end if
            end do
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(outer_v)) outer_v = previous(:,:count)
        if (present(iterations)) iterations = total_steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine lgmres

    ! GCROT(m,k), with oldest-first truncation of an A-image-orthonormal
    ! recycled (C,U) space: C = A*U, C**T C = I. We expose U rather than
    ! cached C because its images must be rebuilt for a changed operator.
    ! GCROT uses a projected right-preconditioned flexible Arnoldi cycle.
    subroutine gcrotmk(a, b, x, info, rtol, atol, maxiter, m, k, x0, jacobi, &
            iterations, residual_norm, ilu, recycle_u)
        type(csr_matrix), intent(in) :: a !! Square canonical CSR matrix, n by n.
        real(dp), intent(in) :: b(:) !! Finite RHS vector of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Approximate solution, length n.
        integer, intent(out) :: info !! 0 converged, >0 cycle cap, <0 invalid/breakdown.
        real(dp), intent(in), optional :: rtol !! Relative true-residual tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! Absolute true-residual tolerance, default zero.
        integer, intent(in), optional :: maxiter !! Maximum outer cycles, default max(1,10*n).
        integer, intent(in), optional :: m !! GMRES columns per cycle, default min(20,n).
        integer, intent(in), optional :: k !! Maximum recycled directions, default min(m,n-1).
        real(dp), intent(in), optional :: x0(:) !! Optional starting solution, length n.
        logical, intent(in), optional :: jacobi !! Right diagonal preconditioner, default false.
        integer, intent(out), optional :: iterations !! Total Arnoldi column evaluations.
        real(dp), intent(out), optional :: residual_norm !! Norm2(b-A*x), unpreconditioned.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional right ILU(0), exclusive with Jacobi.
        real(dp), allocatable, intent(inout), optional :: recycle_u(:,:) !! Recycled U basis, n by <=k.
        real(dp), allocatable :: invdiag(:), u(:,:), c(:,:), empty(:,:), r(:), ax(:)
        real(dp), allocatable :: dx(:), du(:), dc(:), projection(:)
        real(dp) :: threshold, err, gamma, delta, original_norm, image_norm
        integer :: n, mmax, kmax, cap, count, cyc, j, i, pass, start, used, status, total_steps

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a,b,x,info,rtol,atol,maxiter,x0,jacobi,invdiag,threshold,cap,ilu)
        if (info /= 0) return
        n = size(b)
        mmax = min(20,n)
        if (present(m)) mmax = m
        kmax = min(mmax,max(0,n-1))
        if (present(k)) kmax = k
        if ((mmax <= 0 .and. n > 0) .or. kmax < 0) then
            info = sparse_linalg_invalid_input
            return
        end if
        mmax = min(mmax,n)
        kmax = min(kmax,max(0,n-1))
        if (present(recycle_u)) then
            if (allocated(recycle_u)) then
                if (size(recycle_u,1) /= n) then
                    info = sparse_linalg_invalid_input
                    return
                end if
                if (.not. all(ieee_is_finite(recycle_u))) then
                    info = sparse_linalg_invalid_input
                    return
                end if
            end if
        end if
        allocate(u(n,kmax),c(n,kmax),empty(n,0),r(n),ax(n), &
            dx(n),du(n),dc(n),projection(kmax))
        count = 0
        ! Recompute C for the current A, orthogonalizing U simultaneously.
        if (present(recycle_u)) then
            if (allocated(recycle_u)) then
                start = max(1,size(recycle_u,2)-kmax+1)
                do j = start, size(recycle_u,2)
                    du = recycle_u(:,j)
                    call csr_product(a,du,dc)
                    original_norm = norm2(dc)
                    if (.not. ieee_is_finite(original_norm)) cycle
                    do pass = 1, 2
                        do i = 1, count
                            delta = dot_product(c(:,i),dc)
                            dc = dc-delta*c(:,i)
                            du = du-delta*u(:,i)
                        end do
                    end do
                    image_norm = norm2(dc)
                    if (image_norm <= 32.0_dp*epsilon(1.0_dp)*original_norm) cycle
                    if (.not. ieee_is_finite(image_norm)) cycle
                    count = count+1
                    c(:,count) = dc/image_norm
                    u(:,count) = du/image_norm
                end do
            end if
        end if
        call csr_product(a,x,ax)
        r = b-ax
        ! A recycled pair improves the initial guess by orthogonal projection.
        do j = 1, count
            gamma = dot_product(c(:,j),r)
            x = x+gamma*u(:,j)
            r = r-gamma*c(:,j)
        end do
        call csr_product(a,x,ax)
        r = b-ax
        err = norm2(r)
        total_steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = 0
        else
            info = max(1,cap)
            do cyc = 1, cap
                ! Project away roundoff in directions already recycled.
                do j = 1, count
                    gamma = dot_product(c(:,j),r)
                    x = x+gamma*u(:,j)
                    r = r-gamma*c(:,j)
                end do
                call recycled_arnoldi_cycle(a,r,invdiag, &
                    min(n-count,mmax+max(0,kmax-count)),empty,c(:,:count), &
                    threshold,dx,projection(:count),used,status,ilu)
                total_steps = total_steps+used
                if (status /= 0) then
                    info = status
                    exit
                end if
                du = dx
                do j = 1, count
                    du = du-projection(j)*u(:,j)
                end do
                call csr_product(a,du,dc)
                original_norm = norm2(dc)
                if (.not. ieee_is_finite(original_norm)) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                ! Restore C-orthogonality while preserving A*du = dc.
                do pass = 1, 2
                    do j = 1, count
                        delta = dot_product(c(:,j),dc)
                        dc = dc-delta*c(:,j)
                        du = du-delta*u(:,j)
                    end do
                end do
                image_norm = norm2(dc)
                if (.not. ieee_is_finite(image_norm) .or. &
                        image_norm <= 32.0_dp*epsilon(1.0_dp)*original_norm) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                dc = dc/image_norm
                du = du/image_norm
                gamma = dot_product(dc,r)
                x = x+gamma*du
                call csr_product(a,x,ax)
                r = b-ax
                err = norm2(r)
                if (.not. ieee_is_finite(err) .or. .not. all(ieee_is_finite(x))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                if (kmax > 0) then
                    if (count == kmax) then
                        if (kmax > 1) then
                            c(:,1:kmax-1) = c(:,2:kmax)
                            u(:,1:kmax-1) = u(:,2:kmax)
                        end if
                        count = count-1
                    end if
                    count = count+1
                    c(:,count) = dc
                    u(:,count) = du
                end if
                if (err <= threshold) then
                    info = 0
                    exit
                end if
            end do
        end if
        call csr_product(a,x,ax)
        err = norm2(b-ax)
        if (present(recycle_u)) recycle_u = u(:,:count)
        if (present(iterations)) iterations = total_steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. err <= threshold) info = 0
    end subroutine gcrotmk

    ! MINRES for real symmetric (possibly indefinite or singular) systems.
    ! Lanczos three-term recurrence and a streamed QR factorization of its
    ! tridiagonal projection.  No Krylov basis or dense matrix is stored.
    ! Reference: C. C. Paige and M. A. Saunders (1975), SIAM J. Numer. Anal.
    ! SciPy's BSD-3-Clause minres.py serves as the procedural reference;
    ! this version adds true shifted-residual checks and explicit CSR storage.
    subroutine minres(a, b, x, info, rtol, atol, maxiter, x0, shift, jacobi, &
            check_symmetry, iterations, residual_norm)
        type(csr_matrix), intent(in) :: a !! Real symmetric one-based CSR matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side, length n.
        real(dp), allocatable, intent(out) :: x(:) !! Solution or last iterate, length n.
        integer, intent(out) :: info !! 0 converged; >0 iteration cap; <0 invalid/breakdown.
        real(dp), intent(in), optional :: rtol !! True-residual relative tolerance, default 1e-5.
        real(dp), intent(in), optional :: atol !! True-residual absolute tolerance, default 0.
        integer, intent(in), optional :: maxiter !! Iteration cap, default max(1,10*n).
        real(dp), intent(in), optional :: x0(:) !! Optional finite starting guess.
        real(dp), intent(in), optional :: shift !! Solve (A-shift*I)x=b; default 0.
        logical, intent(in), optional :: jacobi !! Inverse-diagonal SPD preconditioner, default false.
        logical, intent(in), optional :: check_symmetry !! Validate A=A**T, default false.
        integer, intent(out), optional :: iterations !! Performed Lanczos steps.
        real(dp), intent(out), optional :: residual_norm !! True norm of b-(A-shift*I)x.
        real(dp), allocatable :: invdiag(:), r_previous(:), r_current(:), y(:), v(:), &
            work(:), w_old(:), w_current(:), w_next(:)
        real(dp) :: threshold, sigma, beta_initial, beta, beta_old, beta_squared
        real(dp) :: alpha, cs, sn, dbar, epsln, oldeps, delta, gbar, gamma
        real(dp) :: phi, phibar, err, eps_scale, mirror
        integer :: n, cap, step, steps, i, j, k, position
        logical :: use_jacobi, verify_symmetry

        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = huge(1.0_dp)
        call prepare(a, b, x, info, rtol, atol, maxiter, x0, &
            jacobi, invdiag, threshold, cap)
        if (info /= sparse_linalg_success) return
        n = size(b)
        sigma = 0.0_dp
        if (present(shift)) sigma = shift
        if (.not. ieee_is_finite(sigma)) then
            info = sparse_linalg_invalid_input
            return
        end if
        use_jacobi = .false.
        if (present(jacobi)) use_jacobi = jacobi
        if (use_jacobi) then
            ! MINRES requires a symmetric *positive-definite* preconditioner,
            ! unlike GMRES/BiCGSTAB, where signed diagonal inverses suffice.
            if (any(invdiag <= 0.0_dp)) then
                info = sparse_linalg_bad_diagonal
                return
            end if
        end if
        verify_symmetry = .false.
        if (present(check_symmetry)) verify_symmetry = check_symmetry
        if (verify_symmetry) then
            do i = 1, n
                do k = a%indptr(i), a%indptr(i+1)-1
                    j = a%indices(k)
                    position = row_position(a, j, i)
                    mirror = 0.0_dp
                    if (position > 0) mirror = a%data(position)
                    eps_scale = max(1.0_dp, abs(a%data(k)), abs(mirror))
                    if (abs(a%data(k)-mirror) > 1.0e-12_dp*eps_scale) then
                        info = sparse_linalg_invalid_input
                        return
                    end if
                end do
            end do
        end if
        allocate(r_previous(n), r_current(n), y(n), v(n), work(n), &
            w_old(n), w_current(n), w_next(n))
        call csr_product(a,x,work)
        r_current = b-work+sigma*x
        err = norm2(r_current)
        steps = 0
        if (.not. ieee_is_finite(err)) then
            info = sparse_linalg_breakdown
        else if (err <= threshold) then
            info = sparse_linalg_success
        else
            y = invdiag*r_current
            beta_squared = dot_product(r_current,y)
            if (beta_squared <= 0.0_dp .or. .not. ieee_is_finite(beta_squared)) then
                info = sparse_linalg_breakdown
            else
                beta_initial = sqrt(beta_squared)
                beta = beta_initial
                beta_old = 1.0_dp
                phibar = beta_initial
                cs = -1.0_dp
                sn = 0.0_dp
                dbar = 0.0_dp
                epsln = 0.0_dp
                r_previous = 0.0_dp
                w_old = 0.0_dp
                w_current = 0.0_dp
                info = max(1,cap)
                do step = 1, cap
                    steps = step
                    v = y/beta
                    call csr_product(a, v, work)
                    work = work-sigma*v
                    if (step >= 2) work = work-(beta/beta_old)*r_previous
                    alpha = dot_product(v,work)
                    work = work-(alpha/beta)*r_current
                    r_previous = r_current
                    r_current = work
                    y = invdiag*r_current
                    beta_squared = dot_product(r_current,y)
                    if (.not. ieee_is_finite(beta_squared) .or. &
                        beta_squared < -100.0_dp*epsilon(1.0_dp)*dot_product(r_current,r_current)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    beta_old = beta
                    beta = sqrt(max(0.0_dp,beta_squared))
                    oldeps = epsln
                    delta = cs*dbar+sn*alpha
                    gbar = sn*dbar-cs*alpha
                    epsln = sn*beta
                    dbar = -cs*beta
                    gamma = hypot(gbar,beta)
                    if (.not. ieee_is_finite(gamma) .or. gamma <= tiny(1.0_dp)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    cs = gbar/gamma
                    sn = beta/gamma
                    phi = cs*phibar
                    phibar = sn*phibar
                    w_next = (v-oldeps*w_old-delta*w_current)/gamma
                    x = x+phi*w_next
                    w_old = w_current
                    w_current = w_next
                    if (.not. all(ieee_is_finite(x))) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    call csr_product(a,x,work)
                    err = norm2(b-work+sigma*x)
                    if (.not. ieee_is_finite(err)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    if (err <= threshold) then
                        info = sparse_linalg_success
                        exit
                    end if
                    ! Exact Lanczos termination without achieving the required
                    ! residual denotes an inconsistent/singular solve.
                    if (beta <= tiny(1.0_dp)) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                end do
            end if
        end if
        call csr_product(a,x,work)
        err = norm2(b-work+sigma*x)
        if (present(iterations)) iterations = steps
        if (present(residual_norm)) residual_norm = err
        if (info >= 0 .and. ieee_is_finite(err)) then
            if (err <= threshold) info = sparse_linalg_success
        end if
    end subroutine minres

    ! Symmetric sparse Ritz pairs via block orthogonal subspace iteration.
    ! A dependency-free real eigsh subset; not a substitute for implicitly
    ! restarted ARPACK on difficult clustered/large-scale spectral problems.
    subroutine eigsh(a, k, eigenvalues, eigenvectors, info, which, tol, &
            maxiter, ncv, v0, iterations, residual_norms, check_symmetry, &
            sigma, inner_rtol, inner_maxiter, method)
        type(csr_matrix), intent(in) :: a !! Symmetric, finite real n-by-n CSR matrix.
        integer, intent(in) :: k !! Number of eigenpairs, 1 <= k < n.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! Ritz values, length k.
        real(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Unit eigenvectors, n by k.
        integer, intent(out) :: info !! 0 if every true residual converged; >0 cap; -1 invalid.
        character(len=*), intent(in), optional :: which !! LA, SA, LM, or SM; default LM; shifted modes use inverse spectrum.
        real(dp), intent(in), optional :: tol !! Relative residual bound, default 1e-8.
        integer, intent(in), optional :: maxiter !! Maximum block updates, default max(300,20*n).
        integer, intent(in), optional :: ncv !! Block subspace dimension, k+1 through n.
        real(dp), intent(in), optional :: v0(:) !! Optional initial direction of length n.
        integer, intent(out), optional :: iterations !! Number of block iterations.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True ||A*v-lambda*v||.
        logical, intent(in), optional :: check_symmetry !! Validate symmetry; default true.
        real(dp), intent(in), optional :: sigma !! Shift for (A-sigma*I)**(-1); shifted SM needs ncv=n.
        real(dp), intent(in), optional :: inner_rtol !! MINRES shifted-solve tolerance, default 1e-12.
        integer, intent(in), optional :: inner_maxiter !! MINRES inner iteration cap, default 10*n.
        character(len=*), intent(in), optional :: method !! 'subspace' (default) or 'lanczos' (supports shift).
        real(dp), allocatable :: q(:,:), aq(:,:), bq(:,:), small(:,:), &
            vecs(:,:), values(:), chosen(:,:), tmp(:,:), residuals(:), inverse(:)
        real(dp) :: tolerance, bound, row_sum, factor, scale, mirror, solve_tol
        integer, allocatable :: order(:)
        integer :: n, p, cap, iter, i, j, t, ii, jj, status, position, inner_cap
        character(len=2) :: target
        logical :: validate_symmetry, converged

        ! Keep the original subspace implementation as the default;
        ! shifted restarted Lanczos uses repeated shifted MINRES solves.
        if (present(method)) then
            select case (method)
            case ('lanczos', 'lanczos_locked')
                call eigsh_lanczos(a,k,eigenvalues,eigenvectors,info, &
                    which=which,tol=tol,maxiter=maxiter,ncv=ncv,v0=v0, &
                    iterations=iterations,residual_norms=residual_norms, &
                    check_symmetry=check_symmetry,sigma=sigma, &
                    inner_rtol=inner_rtol,inner_maxiter=inner_maxiter, &
                    hard_locking=(method == 'lanczos_locked'))
                return
            case ('subspace')
                continue
            case default
                info = sparse_linalg_invalid_input
                if (present(iterations)) iterations = 0
                allocate(eigenvalues(0),eigenvectors(0,0))
                if (present(residual_norms)) allocate(residual_norms(0))
                return
            end select
        end if
        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (n /= a%n_cols .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'SM' .and. target /= 'LA' .and. target /= 'SA') return
        solve_tol = 1.0e-12_dp
        if (present(inner_rtol)) solve_tol = inner_rtol
        if (.not. ieee_is_finite(solve_tol) .or. solve_tol < 0.0_dp) return
        inner_cap = max(1,10*n)
        if (present(inner_maxiter)) inner_cap = inner_maxiter
        if (inner_cap < 0) return
        if (present(sigma)) then
            if (.not. ieee_is_finite(sigma)) return
        end if
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        ! Accept tol=0 as machine-precision request, like SciPy.
        if (tolerance == 0.0_dp) tolerance = 100.0_dp*epsilon(1.0_dp)
        cap = max(300,20*n)
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        p = min(n,max(2*k+6,16))
        if (present(ncv)) p = ncv
        if (p <= k .or. p > n) return
        ! In normal shift-invert, SM selects the *smallest inverse*
        ! eigenvalues (farthest from the shift).  Reduced inverse-subspace
        ! iteration does not reliably retain those modes.
        if (present(sigma)) then
            if (target == 'SM' .and. p /= n) return
        end if
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(ieee_is_finite(v0))) return
            if (norm2(v0) == 0.0_dp) return
        end if
        validate_symmetry = .true.
        if (present(check_symmetry)) validate_symmetry = check_symmetry
        if (validate_symmetry) then
            do i = 1,n
                do t = a%indptr(i),a%indptr(i+1)-1
                    j = a%indices(t)
                    position = row_position(a,j,i)
                    mirror = 0.0_dp
                    if (position > 0) mirror = a%data(position)
                    if (abs(a%data(t)-mirror) > 1.0e-12_dp* &
                        max(1.0_dp,abs(a%data(t)),abs(mirror))) return
                end do
            end do
        end if
        bound = 0.0_dp
        do i = 1,n
            row_sum = sum(abs(a%data(a%indptr(i):a%indptr(i+1)-1)))
            if (.not. ieee_is_finite(row_sum)) return
            bound = max(bound,row_sum)
        end do
        if (bound > sqrt(huge(bound))/2.0_dp) return
        scale = max(1.0_dp,bound)
        allocate(q(n,p),aq(n,p),bq(n,p),small(p,p),vecs(p,p), &
            values(p),order(p),chosen(p,p),tmp(n,p),residuals(k))
        do j = 1,p
            do i = 1,n
                q(i,j) = sin(real(i+1,dp)*real(j+3,dp)*0.6180339887498948_dp) &
                    + cos(real(i+3,dp)*real(j+1,dp)*0.337123_dp)
            end do
        end do
        if (present(v0)) q(:,1) = v0
        call orthonormal_columns(q)
        ! Move allocations only after all arguments have been checked.
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = 0.0_dp
        eigenvectors = 0.0_dp
        residuals = huge(1.0_dp)
        info = max(1,cap)
        do iter = 1,cap
            do j = 1,p
                call csr_product(a,q(:,j),aq(:,j))
            end do
            if (present(sigma)) then
                ! Shift-invert: MINRES solves each real RHS without assembling
                ! a dense inverse or changing the original sparse matrix.
                do j = 1,p
                    call minres(a,q(:,j),inverse,status,rtol=solve_tol, &
                        maxiter=inner_cap,shift=sigma)
                    if (status /= sparse_linalg_success) exit
                    bq(:,j) = inverse
                end do
                if (status /= sparse_linalg_success) then
                    info = sparse_linalg_breakdown
                    exit
                end if
            else
            select case(target)
            case('LA')
                bq = aq+bound*q
            case('SA')
                bq = bound*q-aq
            case('LM')
                do j = 1,p
                    call csr_product(a,aq(:,j),bq(:,j))
                end do
            case('SM')
                do j = 1,p
                    call csr_product(a,aq(:,j),bq(:,j))
                end do
                bq = (bound*bound)*q-bq
            end select
            end if
            if (.not. all(ieee_is_finite(bq))) then
                info = sparse_linalg_breakdown
                exit
            end if
            q = bq
            call orthonormal_columns(q)
            do j = 1,p
                call csr_product(a,q(:,j),aq(:,j))
            end do
            small = matmul(transpose(q),aq)
            small = 0.5_dp*(small+transpose(small))
            call symmetric_jacobi_small(small,values,vecs,status)
            if (status /= 0) then
                info = sparse_linalg_breakdown
                exit
            end if
            order = [(j,j=1,p)]
            ! A block iteration selects the leading eigenpairs of its
            ! transformed operator.  Rotate the orthonormal Ritz subspace
            ! into descending transformed-eigenvalue order.
            do i = 1,p-1
                ii = i
                do j = i+1,p
                    if (present(sigma)) then
                        if (symmetric_shift_priority(values(order(j)),sigma,target) > &
                                symmetric_shift_priority(values(order(ii)),sigma,target)) ii = j
                    else
                        if (spectral_priority(values(order(j)),target) > &
                            spectral_priority(values(order(ii)),target)) ii = j
                    end if
                end do
                jj = order(i)
                order(i) = order(ii)
                order(ii) = jj
            end do
            do j = 1,p
                chosen(:,j) = vecs(:,order(j))
            end do
            tmp = matmul(q,chosen)
            q = tmp
            tmp = matmul(aq,chosen)
            converged = .true.
            do j = 1,k
                eigenvalues(j) = values(order(j))
                eigenvectors(:,j) = q(:,j)
                residuals(j) = norm2(tmp(:,j)-eigenvalues(j)*q(:,j))
                if (.not. ieee_is_finite(residuals(j))) then
                    info = sparse_linalg_breakdown
                    converged = .false.
                    exit
                end if
                if (residuals(j) > tolerance*max(scale,abs(eigenvalues(j)))) converged = .false.
            end do
            if (info == sparse_linalg_breakdown) exit
            if (converged) then
                info = sparse_linalg_success
                exit
            end if
        end do
        if (present(iterations)) iterations = min(iter,cap)
        if (present(residual_norms)) residual_norms = residuals
        ! SciPy returns eigenvalues in ascending algebraic order when
        ! eigenvectors are requested, even if 'LM' was used to select them.
        do i = 1,k-1
            ii = i
            do j = i+1,k
                if (eigenvalues(j) < eigenvalues(ii)) ii = j
            end do
            if (ii /= i) then
                factor = eigenvalues(i)
                eigenvalues(i) = eigenvalues(ii)
                eigenvalues(ii) = factor
                tmp(:,1) = eigenvectors(:,i)
                eigenvectors(:,i) = eigenvectors(:,ii)
                eigenvectors(:,ii) = tmp(:,1)
                factor = residuals(i)
                residuals(i) = residuals(ii)
                residuals(ii) = factor
            end if
        end do
        if (present(residual_norms)) residual_norms = residuals
    end subroutine eigsh

    ! Small complex dense projected eigensolver and nonsymmetric real CSR
    ! eigensolver.  The matrix remains sparse; only the p-by-p Rayleigh--Ritz
    ! projection is dense.  The QR/Schur routines are new MIT-licensed code.
    ! This is a subspace-iteration implementation, not implicitly restarted
    ! ARPACK: strongly nonnormal or clustered problems may not converge.
    ! Symmetric thick-restarted Lanczos with Rayleigh--Ritz extraction.
    !
    ! After the first Krylov cycle, select and retain leading (and already
    ! converged) Ritz vectors.  Full two-pass reorthogonalization is used
    ! when extending the next Krylov chain; this stabilizes clustered and
    ! nearly multiple eigenvalues and costs O(n*ncv) auxiliary storage.
    ! Convergence uses *true* original-operator residuals.  The selected
    ! vectors are soft-locked, not exactly deflated as in ARPACK.
    subroutine eigsh_lanczos(a, k, eigenvalues, eigenvectors, info, which, tol, &
            maxiter, ncv, v0, iterations, residual_norms, locked_count, &
            matvecs, check_symmetry, sigma, inner_rtol, inner_maxiter, &
            hard_locking, block_size)
        type(csr_matrix), intent(in) :: a !! Real symmetric canonical CSR matrix.
        integer, intent(in) :: k !! Requested eigenpairs, 1 <= k < n.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! Eigenvalues in ascending order.
        real(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Normalized vectors, n by k.
        integer, intent(out) :: info !! 0 converged; >0 cycles exhausted; <0 invalid/breakdown.
        character(len=*), intent(in), optional :: which !! LM, LA, SA; SM requires ncv=n; shifted ranks inverse eigenvalues.
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance, default 1e-8.
        integer, intent(in), optional :: maxiter !! Maximum restarts, default 150.
        integer, intent(in), optional :: ncv !! Krylov subspace width k+1 .. n.
        real(dp), intent(in), optional :: v0(:) !! Optional nonzero real starting vector.
        integer, intent(out), optional :: iterations !! Completed Krylov cycles.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True residuals.
        integer, intent(out), optional :: locked_count !! Converged selected Ritz pairs.
        integer, intent(out), optional :: matvecs !! Outer products or inverse solves; excludes inner Krylov steps.
        logical, intent(in), optional :: check_symmetry !! Default true.
        real(dp), intent(in), optional :: sigma !! Shift for normal-mode inverse (A-sigma*I)**(-1).
        real(dp), intent(in), optional :: inner_rtol !! Shifted MINRES relative tolerance, default 1e-12.
        integer, intent(in), optional :: inner_maxiter !! Shifted MINRES iteration limit, default 10*n.
        logical, intent(in), optional :: hard_locking !! Deflate converged symmetric Ritz vectors, default false.
        integer, intent(in), optional :: block_size !! Number of deterministic Krylov seeds for hard locking, default 1.
        real(dp), allocatable :: q(:,:), aq(:,:), h(:,:), vecs(:,:), values(:), &
            retained(:,:), work(:), candidates(:,:), residual(:,:), errors(:), &
            swap(:), inverse(:), original_work(:)
        integer, allocatable :: order(:)
        integer :: n, p, cap, i, j, t, mirror_pos, cycle, nr, r, chosen, &
            status, products, count, tmp, best, inner_cap
        real(dp) :: tolerance, norm_a, norm_v, diagonal_mirror, threshold, &
            old_value, solve_tol, transformed_value
        character(len=2) :: target
        logical :: validate_symmetry, converged, shifted

        if (present(hard_locking)) then
            if (hard_locking) then
                call eigsh_lanczos_deflated(a,k,eigenvalues,eigenvectors,info, &
                    which=which,tol=tol,maxiter=maxiter,ncv=ncv,v0=v0, &
                    iterations=iterations,residual_norms=residual_norms, &
                    locked_count=locked_count,matvecs=matvecs, &
                    check_symmetry=check_symmetry,sigma=sigma, &
                    inner_rtol=inner_rtol,inner_maxiter=inner_maxiter, &
                    block_size=block_size)
                return
            end if
        end if
        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        if (present(locked_count)) locked_count = 0
        if (present(matvecs)) matvecs = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (present(block_size)) return
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (n /= a%n_cols .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'LA' .and. &
            target /= 'SA' .and. target /= 'SM') return
        p = min(n,max(2*k+8,20))
        if (present(ncv)) p = ncv
        if (p <= k .or. p > n) return
        ! A small Ritz space cannot reliably recover interior |lambda| minima
        ! without an inverse spectral transform; do not promise that it can.
        if (target == 'SM' .and. p /= n) return
        shifted = present(sigma)
        if (shifted) then
            if (.not. ieee_is_finite(sigma)) return
        end if
        solve_tol = 1.0e-12_dp
        if (present(inner_rtol)) solve_tol = inner_rtol
        if (.not. ieee_is_finite(solve_tol) .or. solve_tol <= 0.0_dp) return
        inner_cap = max(1,10*n)
        if (present(inner_maxiter)) inner_cap = inner_maxiter
        if (inner_cap < 1) return
        cap = 150
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (tolerance == 0.0_dp) tolerance = 100.0_dp*epsilon(1.0_dp)
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(ieee_is_finite(v0))) return
            if (norm2(v0) <= tiny(1.0_dp)) return
        end if
        validate_symmetry = .true.
        if (present(check_symmetry)) validate_symmetry = check_symmetry
        if (validate_symmetry) then
            do i = 1,n
                do t = a%indptr(i),a%indptr(i+1)-1
                    j = a%indices(t)
                    mirror_pos = row_position(a,j,i)
                    diagonal_mirror = 0.0_dp
                    if (mirror_pos > 0) diagonal_mirror = a%data(mirror_pos)
                    if (abs(a%data(t)-diagonal_mirror) > 1.0e-12_dp* &
                            max(1.0_dp,abs(a%data(t)),abs(diagonal_mirror))) return
                end do
            end do
        end if
        norm_a = 0.0_dp
        do i = 1,n
            norm_a = max(norm_a,sum(abs(a%data(a%indptr(i):a%indptr(i+1)-1))))
        end do
        if (.not. ieee_is_finite(norm_a)) return
        if (norm_a > sqrt(huge(norm_a))/4.0_dp) return
        allocate(q(n,p),aq(n,p),h(p,p),vecs(p,p),values(p), &
            retained(n,p),work(n),original_work(n),candidates(n,k), &
            residual(n,k),errors(k),swap(n),order(p))
        do i = 1,n
            q(i,1) = sin(real(i,dp)*0.716_dp)+cos(real(i,dp)*0.319_dp)
        end do
        if (present(v0)) q(:,1) = v0
        q(:,1) = q(:,1)/norm2(q(:,1))
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = 0.0_dp
        eigenvectors = 0.0_dp
        errors = huge(1.0_dp)
        count = 0
        products = 0
        nr = 1
        info = max(1,cap)
        do cycle = 1,cap
            ! Grow one sparse Krylov chain from the start vector (first
            ! cycle) or the residual appended to retained Ritz vectors.
            do j = nr,p-1
                if (shifted) then
                    call minres(a,q(:,j),inverse,status,rtol=solve_tol, &
                        maxiter=inner_cap,shift=sigma)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    work = inverse
                else
                    call csr_product(a,q(:,j),work)
                end if
                products = products+1
                q(:,j+1) = work
                call lanczos_orthogonalize(q,j+1)
            end do
            if (info == sparse_linalg_breakdown) exit
            do j = 1,p
                if (shifted) then
                    call minres(a,q(:,j),inverse,status,rtol=solve_tol, &
                        maxiter=inner_cap,shift=sigma)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    aq(:,j) = inverse
                else
                    call csr_product(a,q(:,j),aq(:,j))
                end if
                products = products+1
            end do
            if (info == sparse_linalg_breakdown) exit
            h = matmul(transpose(q),aq)
            h = 0.5_dp*(h+transpose(h))
            call symmetric_jacobi_small(h,values,vecs,status)
            if (status /= sparse_linalg_success) then
                info = sparse_linalg_breakdown
                exit
            end if
            order = [(j,j=1,p)]
            do i = 1,p-1
                best = i
                do j = i+1,p
                    if (spectral_priority(values(order(j)),target) > &
                        spectral_priority(values(order(best)),target)) best = j
                end do
                tmp = order(i)
                order(i) = order(best)
                order(best) = tmp
            end do
            converged = .true.
            count = 0
            do j = 1,k
                transformed_value = values(order(j))
                if (shifted) then
                    if (abs(transformed_value) <= 1.0_dp/sqrt(huge(1.0_dp))) exit
                    eigenvalues(j) = sigma + 1.0_dp/transformed_value
                    if (.not. ieee_is_finite(eigenvalues(j))) exit
                else
                    eigenvalues(j) = transformed_value
                end if
                candidates(:,j) = matmul(q,vecs(:,order(j)))
                norm_v = norm2(candidates(:,j))
                if (.not. ieee_is_finite(norm_v) .or. norm_v <= tiny(1.0_dp)) exit
                candidates(:,j) = candidates(:,j)/norm_v
                eigenvectors(:,j) = candidates(:,j)
                ! Retain the transformed residual for the next restart,
                ! but test convergence against the original operator.
                residual(:,j) = matmul(aq,vecs(:,order(j)))/norm_v - &
                    transformed_value*candidates(:,j)
                if (shifted) then
                    call csr_product(a,candidates(:,j),original_work)
                    errors(j) = norm2(original_work-eigenvalues(j)*candidates(:,j))
                else
                    errors(j) = norm2(residual(:,j))
                end if
                if (.not. ieee_is_finite(errors(j))) exit
                threshold = tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))
                if (errors(j) <= threshold) then
                    count = count+1
                else
                    converged = .false.
                end if
            end do
            if (j <= k) then
                info = sparse_linalg_breakdown
                exit
            end if
            if (converged) then
                info = sparse_linalg_success
                exit
            end if
            ! Retain converged first, followed by the remaining best Ritz
            ! vectors and guards. This is soft locking, not hard deflation.
            r = min(p-1,k+1)
            chosen = 0
            do j = 1,k
                if (errors(j) <= tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) then
                    chosen = chosen+1
                    retained(:,chosen) = candidates(:,j)
                end if
            end do
            do j = 1,p
                if (chosen >= r) exit
                if (j <= k) then
                    if (errors(j) <= tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) cycle
                end if
                chosen = chosen+1
                retained(:,chosen) = matmul(q,vecs(:,order(j)))
            end do
            call orthonormal_columns(retained(:,:r))
            q(:,1:r) = retained(:,1:r)
            work = 0.0_dp
            do j = 1,k
                if (errors(j) > tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) &
                    work = work + residual(:,j)/max(errors(j),tiny(1.0_dp))
            end do
            q(:,r+1) = work
            call lanczos_orthogonalize(q,r+1)
            nr = r+1
        end do
        if (present(iterations)) iterations = min(cycle,cap)
        if (present(locked_count)) locked_count = count
        if (present(matvecs)) matvecs = products
        ! SciPy returns real symmetric Ritz pairs in algebraically ascending
        ! order when vectors are requested.  Keep residuals paired on swap.
        do i = 1,k-1
            best = i
            do j = i+1,k
                if (eigenvalues(j) < eigenvalues(best)) best = j
            end do
            if (best == i) cycle
            old_value = eigenvalues(i)
            eigenvalues(i) = eigenvalues(best)
            eigenvalues(best) = old_value
            old_value = errors(i)
            errors(i) = errors(best)
            errors(best) = old_value
            swap = eigenvectors(:,i)
            eigenvectors(:,i) = eigenvectors(:,best)
            eigenvectors(:,best) = swap
        end do
        if (present(residual_norms)) residual_norms = errors
    end subroutine eigsh_lanczos

    ! Hard-locked, thick-restarted symmetric Rayleigh--Ritz iteration.
    !
    ! Converged selected vectors are frozen in an orthonormal locked basis L.
    ! New Krylov directions lie in the orthogonal complement of span(L).
    ! The active projection is Q^T A Q, with Q^T L = 0.  When L is an
    ! approximate invariant subspace (certified by the original residual),
    ! this is an approximate block deflation of A.  Locked vectors never
    ! undergo further power/inverse iteration, unlike soft locking.
    !
    ! This method is valid for *symmetric* A.  Orthogonally projecting out
    ! right eigenvectors alone does not give safe deflation for a general
    ! nonsymmetric operator; the Arnoldi routine keeps its existing behavior.
    subroutine eigsh_lanczos_deflated(a, k, eigenvalues, eigenvectors, info, which, tol, &
            maxiter, ncv, v0, iterations, residual_norms, locked_count, &
            matvecs, check_symmetry, sigma, inner_rtol, inner_maxiter, block_size)
        type(csr_matrix), intent(in) :: a !! Canonical real symmetric sparse operator.
        integer, intent(in) :: k !! Number of eigenpairs requested, between 1 and n-1.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! Selected eigenvalues, ascending.
        real(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Unit vectors, n by k.
        integer, intent(out) :: info !! 0 converged, >0 exhausted, <0 invalid/breakdown.
        character(len=*), intent(in), optional :: which !! LM, LA, SA, or SM; default LM.
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance, default 1e-8.
        integer, intent(in), optional :: maxiter !! Restart-cycle limit, default 150.
        integer, intent(in), optional :: ncv !! Active subspace width, k+1 through n.
        real(dp), intent(in), optional :: v0(:) !! Optional nonzero starting direction.
        integer, intent(out), optional :: iterations !! Completed search cycles.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True residual norms.
        integer, intent(out), optional :: locked_count !! Persistently converged eigenpairs.
        integer, intent(out), optional :: matvecs !! Outer products and shifted solves.
        logical, intent(in), optional :: check_symmetry !! Verify symmetry, default true.
        real(dp), intent(in), optional :: sigma !! Real shift for inverse spectral transformation.
        real(dp), intent(in), optional :: inner_rtol !! MINRES solve tolerance, default 1e-12.
        integer, intent(in), optional :: inner_maxiter !! MINRES maximum iterations, default 10*n.
        integer, intent(in), optional :: block_size !! Independent seeds, 1..k; default 1.
        real(dp), allocatable :: q(:,:), aq(:,:), projected(:,:), small_vectors(:,:), &
            small_values(:), frozen(:,:), frozen_values(:), frozen_errors(:), &
            trial_vectors(:,:), trial_values(:), trial_errors(:), &
            next_basis(:,:), work(:), original_work(:), inverse(:), &
            transformed_residual(:,:), errors(:), swap(:)
        integer, allocatable :: order(:)
        logical, allocatable :: accepted(:)
        integer :: n, p, cap, inner_cap, active, needed, first_new, nlocked, &
            i, j, t, pos, status, cycle, steps, products, keep, nretained, &
            trial_count, old_active, best, old_index, block_width, source
        real(dp) :: tolerance, solve_tol, bound, diagonal_mirror, threshold, &
            norm_v, mu, old_value
        logical :: shifted, verify_symmetry, ok
        character(len=2) :: target

        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        if (present(locked_count)) locked_count = 0
        if (present(matvecs)) matvecs = 0
        allocate(eigenvalues(0),eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (n /= a%n_cols .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'LA' .and. &
                target /= 'SA' .and. target /= 'SM') return
        p = min(n,max(2*k+8,20))
        if (present(ncv)) p = ncv
        if (p <= k .or. p > n) return
        if (present(block_size)) then
            if (block_size < 1 .or. block_size > k) return
        end if
        shifted = present(sigma)
        if (shifted) then
            if (.not. ieee_is_finite(sigma)) return
        end if
        if (target == 'SM' .and. p /= n) return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (tolerance == 0.0_dp) tolerance = 100.0_dp*epsilon(1.0_dp)
        solve_tol = 1.0e-12_dp
        if (present(inner_rtol)) solve_tol = inner_rtol
        if (.not. ieee_is_finite(solve_tol) .or. solve_tol <= 0.0_dp) return
        inner_cap = max(1,10*n)
        if (present(inner_maxiter)) inner_cap = inner_maxiter
        if (inner_cap < 1) return
        cap = 150
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(ieee_is_finite(v0))) return
            if (norm2(v0) <= tiny(1.0_dp)) return
        end if
        verify_symmetry = .true.
        if (present(check_symmetry)) verify_symmetry = check_symmetry
        if (verify_symmetry) then
            do i = 1,n
                do t = a%indptr(i),a%indptr(i+1)-1
                    j = a%indices(t)
                    pos = row_position(a,j,i)
                    diagonal_mirror = 0.0_dp
                    if (pos > 0) diagonal_mirror = a%data(pos)
                    if (abs(a%data(t)-diagonal_mirror) > 1.0e-12_dp* &
                            max(1.0_dp,abs(a%data(t)),abs(diagonal_mirror))) return
                end do
            end do
        end if
        bound = 0.0_dp
        do i = 1,n
            bound = max(bound,sum(abs(a%data(a%indptr(i):a%indptr(i+1)-1))))
        end do
        if (.not. ieee_is_finite(bound)) return
        if (bound > sqrt(huge(bound))/4.0_dp) return
        allocate(q(n,p),aq(n,p),projected(p,p),small_vectors(p,p),small_values(p), &
            frozen(n,k),frozen_values(k),frozen_errors(k), &
            trial_vectors(n,k),trial_values(k),trial_errors(k), &
            next_basis(n,p),work(n),original_work(n),transformed_residual(n,k), &
            errors(k),swap(n),order(p),accepted(k))
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = 0.0_dp
        eigenvectors = 0.0_dp
        errors = huge(1.0_dp)
        frozen = 0.0_dp
        frozen_values = 0.0_dp
        frozen_errors = huge(1.0_dp)
        nlocked = 0
        products = 0
        steps = 0
        ! Optionally use several deterministic seed vectors.  A single
        ! Krylov seed cannot span an exactly repeated eigenspace of
        ! dimension greater than one, even in exact arithmetic.  Multiple
        ! seeds build a small block-Krylov space in the first cycle.
        block_width = 1
        if (present(block_size)) block_width = block_size
        first_new = block_width
        active = p
        q = 0.0_dp
        do i = 1,n
            q(i,1) = sin(real(i,dp)*0.716_dp)+cos(real(i,dp)*0.319_dp)
        end do
        if (present(v0)) q(:,1) = v0
        do j = 1,block_width
            if (j > 1) then
                do i = 1,n
                    q(i,j) = sin((real(i,dp)+3.0_dp)*real(j,dp)*0.2718281828459045_dp) &
                        +cos(real(i,dp)*(1.0_dp+2.0_dp*real(j,dp))*0.395928_dp)
                end do
            end if
            call orthogonalize_deflated(q,j,frozen,nlocked,ok)
            if (.not. ok) then
                info = sparse_linalg_breakdown
                return
            end if
        end do
        info = max(1,cap)
        do cycle = 1,cap
            steps = cycle
            ! Complete the active Krylov basis, orthogonal to locked modes.
            do j = first_new,active-1
                source = j
                if (cycle == 1) source = j+1-block_width
                if (shifted) then
                    call minres(a,q(:,source),inverse,status,rtol=solve_tol, &
                        maxiter=inner_cap,shift=sigma)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    work = inverse
                else
                    call csr_product(a,q(:,source),work)
                end if
                products = products+1
                q(:,j+1) = work
                call orthogonalize_deflated(q,j+1,frozen,nlocked,ok)
                if (.not. ok) then
                    info = sparse_linalg_breakdown
                    exit
                end if
            end do
            if (info == sparse_linalg_breakdown) exit
            do j = 1,active
                if (shifted) then
                    call minres(a,q(:,j),inverse,status,rtol=solve_tol, &
                        maxiter=inner_cap,shift=sigma)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    aq(:,j) = inverse
                else
                    call csr_product(a,q(:,j),aq(:,j))
                end if
                products = products+1
            end do
            if (info == sparse_linalg_breakdown) exit
            projected(:active,:active) = matmul(transpose(q(:,:active)),aq(:,:active))
            projected(:active,:active) = 0.5_dp* &
                (projected(:active,:active)+transpose(projected(:active,:active)))
            call symmetric_jacobi_small(projected(:active,:active), &
                small_values(:active),small_vectors(:active,:active),status)
            if (status /= sparse_linalg_success) then
                info = sparse_linalg_breakdown
                exit
            end if
            order(:active) = [(j,j=1,active)]
            do i = 1,active-1
                best = i
                do j = i+1,active
                    if (spectral_priority(small_values(order(j)),target) > &
                            spectral_priority(small_values(order(best)),target)) best = j
                end do
                old_index = order(i)
                order(i) = order(best)
                order(best) = old_index
            end do
            needed = k-nlocked
            trial_count = 0
            accepted = .false.
            do j = 1,needed
                mu = small_values(order(j))
                if (shifted) then
                    if (abs(mu) <= 1.0_dp/sqrt(huge(1.0_dp))) exit
                    trial_values(j) = sigma+1.0_dp/mu
                    if (.not. ieee_is_finite(trial_values(j))) exit
                else
                    trial_values(j) = mu
                end if
                trial_vectors(:,j) = matmul(q(:,:active), &
                    small_vectors(:active,order(j)))
                norm_v = norm2(trial_vectors(:,j))
                if (norm_v <= tiny(norm_v) .or. .not. ieee_is_finite(norm_v)) exit
                trial_vectors(:,j) = trial_vectors(:,j)/norm_v
                call csr_product(a,trial_vectors(:,j),original_work)
                trial_errors(j) = norm2(original_work-trial_values(j)*trial_vectors(:,j))
                if (.not. ieee_is_finite(trial_errors(j))) exit
                transformed_residual(:,j) = matmul(aq(:,:active), &
                    small_vectors(:active,order(j)))/norm_v - mu*trial_vectors(:,j)
                threshold = tolerance*max(1.0_dp,bound,abs(trial_values(j)))
                accepted(j) = trial_errors(j) <= threshold
                trial_count = trial_count+1
            end do
            if (trial_count /= needed) then
                info = sparse_linalg_breakdown
                exit
            end if
            ! Snapshot the last iteration's selected modes, even on exhaustion.
            do j = 1,nlocked
                eigenvalues(j) = frozen_values(j)
                eigenvectors(:,j) = frozen(:,j)
                errors(j) = frozen_errors(j)
            end do
            do j = 1,needed
                eigenvalues(nlocked+j) = trial_values(j)
                eigenvectors(:,nlocked+j) = trial_vectors(:,j)
                errors(nlocked+j) = trial_errors(j)
            end do
            ! Hard lock newly converged Ritz pairs and exclude them from the
            ! next Krylov space.  Only selected Ritz pairs may be frozen.
            do j = 1,needed
                if (.not. accepted(j)) cycle
                nlocked = nlocked+1
                frozen(:,nlocked) = trial_vectors(:,j)
                frozen_values(nlocked) = trial_values(j)
                frozen_errors(nlocked) = trial_errors(j)
            end do
            if (nlocked == k) then
                info = sparse_linalg_success
                eigenvalues = frozen_values
                eigenvectors = frozen
                errors = frozen_errors
                exit
            end if
            ! Retain unconverged selected Ritz vectors plus one guard Ritz
            ! vector.  All new trial vectors are explicitly deflated against
            ! the enlarged locked subspace, and a Ritz residual starts the
            ! next Krylov expansion.
            old_active = active
            active = min(p,n-nlocked)
            keep = min(active-1,k-nlocked+1)
            next_basis = 0.0_dp
            nretained = 0
            do j = 1,needed
                if (accepted(j)) cycle
                if (nretained >= keep) exit
                nretained = nretained+1
                next_basis(:,nretained) = trial_vectors(:,j)
            end do
            do j = needed+1,old_active
                if (nretained >= keep) exit
                nretained = nretained+1
                next_basis(:,nretained) = matmul(q(:,:old_active), &
                    small_vectors(:old_active,order(j)))
            end do
            ! In the rare case where many vectors were newly frozen, fill
            ! missing guards with deterministic orthogonal complements.
            do j = nretained+1,keep
                next_basis(:,j) = 0.0_dp
            end do
            work = 0.0_dp
            do j = 1,needed
                if (.not. accepted(j)) then
                    work = work + transformed_residual(:,j)/ &
                        max(norm2(transformed_residual(:,j)),tiny(1.0_dp))
                end if
            end do
            next_basis(:,keep+1) = work
            q(:,1:keep+1) = next_basis(:,1:keep+1)
            do j = 1,keep+1
                call orthogonalize_deflated(q,j,frozen,nlocked,ok)
                if (.not. ok) then
                    info = sparse_linalg_breakdown
                    exit
                end if
            end do
            if (info == sparse_linalg_breakdown) exit
            first_new = keep+1
        end do
        if (present(iterations)) iterations = steps
        if (present(locked_count)) locked_count = nlocked
        if (present(matvecs)) matvecs = products
        ! Report selected modes in the same ascending order as eigsh.
        do i = 1,k-1
            best = i
            do j = i+1,k
                if (eigenvalues(j) < eigenvalues(best)) best = j
            end do
            if (best == i) cycle
            old_value = eigenvalues(i)
            eigenvalues(i) = eigenvalues(best)
            eigenvalues(best) = old_value
            old_value = errors(i)
            errors(i) = errors(best)
            errors(best) = old_value
            swap = eigenvectors(:,i)
            eigenvectors(:,i) = eigenvectors(:,best)
            eigenvectors(:,best) = swap
        end do
        if (present(residual_norms)) residual_norms = errors
    end subroutine eigsh_lanczos_deflated

    ! Fully reorthogonalize one active basis direction against all locked
    ! symmetric modes and the preceding active directions.  Deterministic
    ! coordinate complements handle invariant-subspace breakdown without RNG.
    subroutine orthogonalize_deflated(q, j, frozen, nlocked, ok)
        real(dp), intent(inout) :: q(:,:) !! Active basis, with candidate in column j.
        integer, intent(in) :: j !! One-based candidate column index.
        real(dp), intent(in) :: frozen(:,:) !! Orthogonal hard-locked eigenvectors.
        integer, intent(in) :: nlocked !! Number of frozen basis columns.
        logical, intent(out) :: ok !! False if no independent vector exists.
        integer :: n, trial, pass, i
        real(dp) :: magnitude
        n = size(q,1)
        ok = .false.
        do trial = 0,n
            do pass = 1,2
                do i = 1,nlocked
                    q(:,j) = q(:,j)-dot_product(frozen(:,i),q(:,j))*frozen(:,i)
                end do
                do i = 1,j-1
                    q(:,j) = q(:,j)-dot_product(q(:,i),q(:,j))*q(:,i)
                end do
            end do
            magnitude = norm2(q(:,j))
            if (ieee_is_finite(magnitude)) then
                if (magnitude > 100.0_dp*epsilon(1.0_dp)) then
                    q(:,j) = q(:,j)/magnitude
                    ok = .true.
                    return
                end if
            end if
            q(:,j) = 0.0_dp
            if (trial < n) q(mod(j+trial-1,n)+1,j) = 1.0_dp
        end do
    end subroutine orthogonalize_deflated

    ! Two-pass modified Gram--Schmidt for a Lanczos expansion direction.
    ! Fallback to a deterministic orthogonal complement on exact breakdown,
    ! allowing a basis to span invariant subspaces and zero matrices.
    subroutine lanczos_orthogonalize(q,j)
        real(dp), intent(inout) :: q(:,:) !! Existing orthonormal columns; normalize column j.
        integer, intent(in) :: j !! One-based candidate column, j >= 2.
        integer :: pass, i, trial, n
        real(dp) :: magnitude
        n = size(q,1)
        do trial = 0,n
            do pass = 1,2
                do i = 1,j-1
                    q(:,j) = q(:,j)-dot_product(q(:,i),q(:,j))*q(:,i)
                end do
            end do
            magnitude = norm2(q(:,j))
            if (magnitude > 100.0_dp*epsilon(1.0_dp)) then
                q(:,j) = q(:,j)/magnitude
                return
            end if
            q(:,j) = 0.0_dp
            q(mod(j+trial-1,n)+1,j) = 1.0_dp
        end do
        q(:,j) = 0.0_dp
    end subroutine lanczos_orthogonalize

    subroutine eigs(a, k, eigenvalues, eigenvectors, info, which, tol, &
            maxiter, ncv, v0, iterations, residual_norms, sigma, &
            inner_rtol, inner_maxiter, inner_restart, method, block_size)
        type(csr_matrix), intent(in) :: a !! Finite real n-by-n canonical CSR matrix.
        integer, intent(in) :: k !! Number of eigenpairs; 1 <= k < n.
        complex(dp), allocatable, intent(out) :: eigenvalues(:) !! Complex Ritz values, length k.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Complex unit vectors, n by k.
        integer, intent(out) :: info !! 0 success, >0 nonconvergence, -1 invalid, -2 breakdown.
        character(len=*), intent(in), optional :: which !! LM, SM, LR, SR, LI, SI; reduced SM uses zero-shift inverse.
        real(dp), intent(in), optional :: tol !! Relative true-residual bound; default 1e-8.
        integer, intent(in), optional :: maxiter !! Maximum subspace updates, default 600.
        integer, intent(in), optional :: ncv !! Subspace size, k+1 .. n; default max(20,2*k+8).
        complex(dp), intent(in), optional :: v0(:) !! Nonzero finite initial search direction, length n.
        integer, intent(out), optional :: iterations !! Number of subspace updates.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True eigenpair residuals.
        complex(dp), intent(in), optional :: sigma !! Optional complex shift; use which='LM' for nearest eigenvalues.
        real(dp), intent(in), optional :: inner_rtol !! Relative shifted-GMRES residual target, default 1e-11.
        integer, intent(in), optional :: inner_maxiter !! Total inner Arnoldi iterations, default max(200,10*n).
        integer, intent(in), optional :: inner_restart !! Restart dimension of inner complex GMRES, default min(40,n).
        character(len=*), intent(in), optional :: method !! 'subspace', 'arnoldi', or 'arnoldi_schur'.
        integer, intent(in), optional :: block_size !! Independent Arnoldi search directions; default 1.
        complex(dp), allocatable :: q(:,:), aq(:,:), work(:,:), small(:,:), &
            vecs(:,:), projected_vals(:), residuals_vec(:), cand(:)
        real(dp), allocatable :: errors(:)
        integer, allocatable :: order(:)
        real(dp) :: bound, tolerance, scale, r, solve_tol
        complex(dp), allocatable :: shifted_rhs(:)
        complex(dp) :: effective_sigma
        integer :: solve_cap, solve_restart
        integer :: n, p, cap, i, j, t, ii, it, code
        character(len=2) :: target
        logical :: converged, shifted

        ! Restarted Arnoldi uses the existing complex shifted-GMRES kernel
        ! when requested; the original subspace method remains default.
        if (present(method)) then
            select case (method)
            case ('arnoldi','arnoldi_schur')
                call eigs_arnoldi(a,k,eigenvalues,eigenvectors,info, &
                    which=which,tol=tol,maxiter=maxiter,ncv=ncv,v0=v0, &
                    iterations=iterations,residual_norms=residual_norms, &
                    sigma=sigma,inner_rtol=inner_rtol, &
                    inner_maxiter=inner_maxiter,inner_restart=inner_restart, &
                    schur_locking=(method == 'arnoldi_schur'),block_size=block_size)
                return
            case ('subspace')
                continue
            case default
                info = sparse_linalg_invalid_input
                if (present(iterations)) iterations = 0
                allocate(eigenvalues(0),eigenvectors(0,0))
                if (present(residual_norms)) allocate(residual_norms(0))
                return
            end select
        end if
        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        ! Block search belongs to restarted Arnoldi, not block subspace iteration.
        if (present(block_size)) return
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (n /= a%n_cols .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'LR' .and. target /= 'SR' .and. &
                target /= 'LI' .and. target /= 'SI' .and. target /= 'SM') return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (tolerance == 0.0_dp) tolerance = 100.0_dp*epsilon(1.0_dp)
        cap = 600
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        p = min(n,max(2*k+8,20))
        if (present(ncv)) p = ncv
        if (p <= k .or. p > n) return
        if (present(sigma)) then
            if (.not. ieee_is_finite(real(sigma,dp)) .or. &
                    .not. ieee_is_finite(aimag(sigma))) return
            ! A real operator with a complex shift uses a complex inverse;
            ! SciPy's default real-OPpart transformed operator differs.
            if (target /= 'LM') return
        end if
        shifted = present(sigma) .or. (target == 'SM' .and. p < n)
        effective_sigma = cmplx(0.0_dp,0.0_dp,kind=dp)
        if (present(sigma)) effective_sigma = sigma
        solve_tol = 1.0e-11_dp
        if (present(inner_rtol)) solve_tol = inner_rtol
        if (.not. ieee_is_finite(solve_tol) .or. solve_tol <= 0.0_dp) return
        solve_cap = max(200,10*n)
        if (present(inner_maxiter)) solve_cap = inner_maxiter
        if (solve_cap < 1) return
        solve_restart = min(n,40)
        if (present(inner_restart)) solve_restart = inner_restart
        if (solve_restart < 1) return
        solve_restart = min(n,solve_restart)
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(ieee_is_finite(real(v0,dp))) .or. &
                    .not. all(ieee_is_finite(aimag(v0)))) return
            if (complex_norm(v0) == 0.0_dp) return
        end if
        bound = 0.0_dp
        do i = 1,n
            r = sum(abs(a%data(a%indptr(i):a%indptr(i+1)-1)))
            if (.not. ieee_is_finite(r)) return
            bound = max(bound,r)
        end do
        if (bound > sqrt(huge(bound))/4.0_dp) return
        scale = max(1.0_dp,bound)
        allocate(q(n,p),aq(n,p),work(n,p),small(p,p), &
            vecs(p,p),projected_vals(p),residuals_vec(n), &
            cand(n),errors(k),order(p))
        do j = 1,p
            do i = 1,n
                q(i,j) = cmplx(sin(real(i+1,dp)*real(j+3,dp)*0.6180339887498948_dp), &
                    cos(real(i+3,dp)*real(j+1,dp)*0.337123_dp),kind=dp)
            end do
        end do
        if (present(v0)) q(:,1) = v0
        call orthonormal_complex(q)
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = cmplx(0.0_dp,0.0_dp,kind=dp)
        eigenvectors = cmplx(0.0_dp,0.0_dp,kind=dp)
        errors = huge(1.0_dp)
        info = max(1,cap)
        do it = 1,cap
            do j = 1,p
                call csr_complex_product(a,q(:,j),aq(:,j))
            end do
            if (shifted) then
                do j = 1,p
                    call complex_shift_gmres(a,effective_sigma,q(:,j),shifted_rhs,code, &
                        solve_tol,solve_cap,solve_restart)
                    if (code /= sparse_linalg_success) exit
                    work(:,j) = shifted_rhs
                end do
                if (code /= sparse_linalg_success) then
                    info = sparse_linalg_breakdown
                    exit
                end if
            else
            select case(target)
            case('LM','SM')
                work = aq
            case('LR')
                work = aq+cmplx(bound,0.0_dp,kind=dp)*q
            case('SR')
                work = cmplx(bound,0.0_dp,kind=dp)*q-aq
            case('LI')
                work = aq+cmplx(0.0_dp,bound,kind=dp)*q
            case('SI')
                work = aq-cmplx(0.0_dp,bound,kind=dp)*q
            end select
            end if
            if (.not. all(ieee_is_finite(real(work,dp))) .or. &
                    .not. all(ieee_is_finite(aimag(work)))) then
                info = sparse_linalg_breakdown
                exit
            end if
            q = work
            call orthonormal_complex(q)
            do j = 1,p
                call csr_complex_product(a,q(:,j),aq(:,j))
            end do
            small = matmul(conjg(transpose(q)),aq)
            call complex_schur_eigensystem(small,projected_vals,vecs,code)
            if (code /= sparse_linalg_success) then
                info = sparse_linalg_breakdown
                exit
            end if
            order = [(j,j=1,p)]
            do i = 1,p-1
                ii = i
                do j = i+1,p
                    if (shifted) then
                        if (abs(projected_vals(order(j))-effective_sigma) < &
                                abs(projected_vals(order(ii))-effective_sigma)) ii = j
                    else
                        if (nonsymmetric_priority(projected_vals(order(j)),target) > &
                                nonsymmetric_priority(projected_vals(order(ii)),target)) ii = j
                    end if
                end do
                t = order(i)
                order(i) = order(ii)
                order(ii) = t
            end do
            converged = .true.
            do j = 1,k
                eigenvalues(j) = projected_vals(order(j))
                cand = matmul(q,vecs(:,order(j)))
                r = complex_norm(cand)
                if (r <= tiny(r) .or. .not. ieee_is_finite(r)) then
                    info = sparse_linalg_breakdown
                    converged = .false.
                    exit
                end if
                eigenvectors(:,j) = cand/r
                call csr_complex_product(a,eigenvectors(:,j),residuals_vec)
                errors(j) = complex_norm(residuals_vec-eigenvalues(j)*eigenvectors(:,j))
                if (.not. ieee_is_finite(errors(j))) then
                    info = sparse_linalg_breakdown
                    converged = .false.
                    exit
                end if
                if (errors(j) > tolerance*max(scale,abs(eigenvalues(j)))) converged = .false.
            end do
            if (info == sparse_linalg_breakdown) exit
            if (converged) then
                info = sparse_linalg_success
                exit
            end if
        end do
        if (present(iterations)) iterations = min(it,cap)
        if (present(residual_norms)) residual_norms = errors
    end subroutine eigs


    ! Thick-restarted Arnoldi/Rayleigh--Ritz solver for a real CSR operator.
    !
    ! A full Arnoldi basis is grown initially.  After Rayleigh--Ritz, the
    ! leading Ritz vectors are retained explicitly and reorthogonalized, and
    ! the next Krylov expansion is seeded by a Ritz residual.  By default,
    ! converged Ritz vectors are kept first (soft locking).  With the opt-in
    ! schur_locking mode, an approximately invariant Schur subspace is frozen
    ! and new search vectors are projected into its orthogonal complement;
    ! full Rayleigh--Ritz extraction retains the nonsymmetric coupling from
    ! the Schur subspace to the complement.  This is NOT implicit QR restarting.
    ! Convergence is always checked with ||A*x-lambda*x||, not a projected
    ! residual estimate, including in complex shift-invert mode. Optional
    ! block_size > 1 advances independent Krylov search directions and
    ! restarts from separate Ritz residuals to capture repeated eigenspaces.
    ! The original block-iteration eigs remains available.
    subroutine eigs_arnoldi(a, k, eigenvalues, eigenvectors, info, which, tol, &
            maxiter, ncv, v0, iterations, residual_norms, locked_count, &
            matvecs, sigma, inner_rtol, inner_maxiter, inner_restart, &
            schur_locking, schur_locked_count, block_size)
        type(csr_matrix), intent(in) :: a !! Square real canonical CSR matrix.
        integer, intent(in) :: k !! Number of eigenpairs, 1 <= k < n.
        complex(dp), allocatable, intent(out) :: eigenvalues(:) !! Complex Ritz values, length k.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Complex Ritz vectors, n by k.
        integer, intent(out) :: info !! 0 true residuals pass; positive exhaustion; negative error.
        character(len=*), intent(in), optional :: which !! LM, LR, SR, LI, SI; SM needs ncv=n; shifted ranks inverse Ritz values.
        real(dp), intent(in), optional :: tol !! Relative true residual target; default 1e-8.
        integer, intent(in), optional :: maxiter !! Max restart cycles, default 150.
        integer, intent(in), optional :: ncv !! Subspace dimension, k+1 through n.
        complex(dp), intent(in), optional :: v0(:) !! Optional finite nonzero initial vector.
        integer, intent(out), optional :: iterations !! Number of restart cycles performed.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True residual norms.
        integer, intent(out), optional :: locked_count !! Wanted Ritz pairs passing tolerance.
        integer, intent(out), optional :: matvecs !! Outer products or inverse solves, excluding inner Krylov steps.
        complex(dp), intent(in), optional :: sigma !! Complex normal-mode shift for inverse (A-sigma*I)**(-1).
        real(dp), intent(in), optional :: inner_rtol !! Complex GMRES relative tolerance, default 1e-11.
        integer, intent(in), optional :: inner_maxiter !! Complex GMRES iteration limit, default max(200,10*n).
        integer, intent(in), optional :: inner_restart !! Complex GMRES restart dimension, default min(40,n).
        logical, intent(in), optional :: schur_locking !! Freeze certified invariant Schur subspace, default false.
        integer, intent(out), optional :: schur_locked_count !! Dimension of frozen Schur subspace.
        integer, intent(in), optional :: block_size !! Number of independent Arnoldi seeds, 1..k; default 1.
        complex(dp), allocatable :: q(:,:), aq(:,:), h(:,:), y(:,:), values(:), &
            retained(:,:), work(:), cand(:,:), residual(:,:), initial(:), &
            inverse(:), original_work(:), locked(:,:), trial(:,:), op_trial(:,:), &
            lock_matrix(:,:), lock_basis(:,:), lock_vectors(:,:), &
            lock_values(:), new_locked(:,:), candidate(:)
        real(dp), allocatable :: errors(:)
        integer, allocatable :: order(:)
        real(dp) :: threshold, norm_a, t, norm_v, tolerance, solve_tol, leakage, op_scale
        integer :: n, p, cap, cycle, j, i, ii, tmp, nr, count, products, &
            status, inner_cap, inner_restart_size, nlocked, previous_locks, pass, m, l, &
            block_width, source_column, retained_count, seed_count
        complex(dp) :: transformed_value
        character(len=2) :: target
        logical :: converged, shifted, use_schur, acceptable

        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        if (present(locked_count)) locked_count = 0
        if (present(matvecs)) matvecs = 0
        if (present(schur_locked_count)) schur_locked_count = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (a%n_cols /= n .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'LR' .and. target /= 'SR' .and. &
            target /= 'LI' .and. target /= 'SI' .and. target /= 'SM') return
        p = min(n,max(2*k+8,20))
        if (present(ncv)) p = ncv
        if (p <= k .or. p > n) return
        block_width = 1
        if (present(block_size)) block_width = block_size
        if (block_width < 1 .or. block_width > k) return
        ! A restarted block needs k wanted vectors, b residual seeds and
        ! one expansion slot; a full-size Krylov space needs no restart.
        if (block_width > 1 .and. p < k+block_width+1 .and. p /= n) return
        if (target == 'SM' .and. p /= n) return
        use_schur = .false.
        if (present(schur_locking)) use_schur = schur_locking
        shifted = present(sigma)
        if (shifted) then
            if (.not. ieee_is_finite(real(sigma,dp)) .or. &
                .not. ieee_is_finite(aimag(sigma))) return
        end if
        solve_tol = 1.0e-11_dp
        if (present(inner_rtol)) solve_tol = inner_rtol
        if (.not. ieee_is_finite(solve_tol) .or. solve_tol <= 0.0_dp) return
        inner_cap = max(200,10*n)
        if (present(inner_maxiter)) inner_cap = inner_maxiter
        if (inner_cap < 1) return
        inner_restart_size = min(n,40)
        if (present(inner_restart)) inner_restart_size = inner_restart
        if (inner_restart_size < 1) return
        inner_restart_size = min(n,inner_restart_size)
        cap = 150
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (tolerance == 0.0_dp) tolerance = 100.0_dp*epsilon(1.0_dp)
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(ieee_is_finite(real(v0,dp))) .or. &
                .not. all(ieee_is_finite(aimag(v0)))) return
            if (complex_norm(v0) == 0.0_dp) return
        end if
        norm_a = 0.0_dp
        do i = 1,n
            t = sum(abs(a%data(a%indptr(i):a%indptr(i+1)-1)))
            if (.not. ieee_is_finite(t)) return
            norm_a = max(norm_a,t)
        end do
        if (norm_a > sqrt(huge(norm_a))/4.0_dp) return
        allocate(q(n,p), aq(n,p), h(p,p), y(p,p), values(p), &
            retained(n,p), work(n), original_work(n), cand(n,k), residual(n,k), &
            initial(n), errors(k), order(p))
        if (use_schur) allocate(locked(n,k), trial(n,k), op_trial(n,k), candidate(n))
        do i = 1,n
            initial(i) = cmplx(sin(real(i,dp)*0.716_dp), &
                cos(real(i,dp)*0.319_dp),kind=dp)
        end do
        if (present(v0)) initial = v0
        initial = initial/complex_norm(initial)
        q = cmplx(0.0_dp,0.0_dp,kind=dp)
        q(:,1) = initial
        ! Multiple independent starting vectors preserve an eigenspace even
        ! when its eigenvalues are exactly repeated.  The basis is generated
        ! deterministically, without changing the global random-number state.
        do j = 2,block_width
            do i = 1,n
                q(i,j) = cmplx(sin(real(i+2*j,dp)*real(j+3,dp)*0.6180339887498948_dp), &
                    cos(real(i+3*j,dp)*real(j+1,dp)*0.337123_dp),kind=dp)
            end do
            call arnoldi_orthogonalize(q,j)
        end do
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = cmplx(0.0_dp,0.0_dp,kind=dp)
        eigenvectors = cmplx(0.0_dp,0.0_dp,kind=dp)
        errors = huge(1.0_dp)
        count = 0
        products = 0
        nr = block_width
        nlocked = 0
        info = max(1,cap)
        do cycle = 1,cap
            ! Block Arnoldi with a two-pass orthogonalization.  In each sweep
            ! A*q(:,j-b+1) generates the next column, so b independent
            ! directions evolve in parallel rather than collapsing onto a
            ! single Krylov sequence.  b=1 is precisely the original path.
            do j = nr,p-1
                source_column = j-block_width+1
                if (shifted) then
                    call complex_shift_gmres(a,sigma,q(:,source_column),inverse,status, &
                        solve_tol,inner_cap,inner_restart_size)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    work = inverse
                else
                    call csr_complex_product(a,q(:,source_column),work)
                end if
                products = products+1
                q(:,j+1) = work
                call arnoldi_orthogonalize(q,j+1)
            end do
            if (info == sparse_linalg_breakdown) exit
            do j = 1,p
                if (shifted) then
                    call complex_shift_gmres(a,sigma,q(:,j),inverse,status, &
                        solve_tol,inner_cap,inner_restart_size)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    aq(:,j) = inverse
                else
                    call csr_complex_product(a,q(:,j),aq(:,j))
                end if
                products = products+1
            end do
            if (info == sparse_linalg_breakdown) exit
            h = matmul(conjg(transpose(q)),aq)
            call complex_schur_eigensystem(h,values,y,status)
            if (status /= sparse_linalg_success) then
                info = sparse_linalg_breakdown
                exit
            end if
            order = [(j,j=1,p)]
            do i = 1,p-1
                ii = i
                do j = i+1,p
                    if (nonsymmetric_priority(values(order(j)),target) > &
                        nonsymmetric_priority(values(order(ii)),target)) ii = j
                end do
                tmp = order(i)
                order(i) = order(ii)
                order(ii) = tmp
            end do
            converged = .true.
            count = 0
            do j = 1,k
                transformed_value = values(order(j))
                if (shifted) then
                    if (abs(transformed_value) <= 1.0_dp/sqrt(huge(1.0_dp))) exit
                    eigenvalues(j) = sigma + 1.0_dp/transformed_value
                    if (.not. ieee_is_finite(real(eigenvalues(j),dp)) .or. &
                        .not. ieee_is_finite(aimag(eigenvalues(j)))) exit
                else
                    eigenvalues(j) = transformed_value
                end if
                cand(:,j) = matmul(q,y(:,order(j)))
                norm_v = complex_norm(cand(:,j))
                if (norm_v <= tiny(norm_v) .or. .not. ieee_is_finite(norm_v)) exit
                cand(:,j) = cand(:,j)/norm_v
                eigenvectors(:,j) = cand(:,j)
                ! Reuse A*Q to form the true residual of each Ritz vector.
                ! Restart with the inverse-operator Ritz residual, but
                ! certify eigenpairs with the original unshifted CSR A.
                residual(:,j) = matmul(aq,y(:,order(j)))/norm_v &
                    -transformed_value*cand(:,j)
                if (shifted) then
                    call csr_complex_product(a,cand(:,j),original_work)
                    errors(j) = complex_norm(original_work-eigenvalues(j)*cand(:,j))
                else
                    errors(j) = complex_norm(residual(:,j))
                end if
                if (.not. ieee_is_finite(errors(j))) exit
                threshold = tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))
                if (errors(j) <= threshold) then
                    count = count+1
                else
                    converged = .false.
                end if
            end do
            if (j <= k) then
                info = sparse_linalg_breakdown
                exit
            end if
            if (converged) then
                info = sparse_linalg_success
                exit
            end if
            ! A full search space cannot be restarted with additional block
            ! directions; return the nonconvergence status without writing
            ! beyond the allocated basis if its Ritz solve was inaccurate.
            if (block_width > 1 .and. p < k+block_width+1) exit
            if (use_schur) then
                ! Lock only Ritz vectors that satisfy true residual bounds.
                ! Right eigenvectors need not be mutually orthogonal for a
                ! nonnormal matrix; their ORTHONORMAL SPAN, not their individual
                ! directions, is the invariant Schur subspace.  The projection
                ! is used only for searching: locked vectors remain explicitly
                ! in the full Rayleigh--Ritz problem to retain the upper-right
                ! coupling needed to reconstruct correct right eigenvectors.
                previous_locks = nlocked
                do j = 1,k
                    threshold = tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))
                    if (errors(j) > threshold) cycle
                    if (nlocked >= k) exit
                    candidate = cand(:,j)
                    do pass = 1,2
                        do l = 1,nlocked
                            candidate = candidate-dot_product(locked(:,l),candidate)*locked(:,l)
                        end do
                    end do
                    norm_v = complex_norm(candidate)
                    ! Do not magnify a small near-dependent component of a
                    ! nonnormal eigenvector into an inaccurate locked vector.
                    if (norm_v <= 1.0e-5_dp) cycle
                    trial(:,1:nlocked) = locked(:,1:nlocked)
                    trial(:,nlocked+1) = candidate/norm_v
                    acceptable = .true.
                    do l = 1,nlocked+1
                        if (shifted) then
                            call complex_shift_gmres(a,sigma,trial(:,l),inverse,status, &
                                solve_tol,inner_cap,inner_restart_size)
                            if (status /= sparse_linalg_success) then
                                acceptable = .false.
                                exit
                            end if
                            op_trial(:,l) = inverse
                        else
                            call csr_complex_product(a,trial(:,l),op_trial(:,l))
                        end if
                        products = products+1
                    end do
                    if (.not. acceptable) cycle
                    ! Check the invariant SUBSPACE residual, not just the
                    ! separate eigenpair residuals.  The threshold scales with
                    ! the active operator (A or its complex shifted inverse).
                    m = nlocked+1
                    leakage = sqrt(sum(abs(op_trial(:,:m)-matmul(trial(:,:m), &
                        matmul(conjg(transpose(trial(:,:m))),op_trial(:,:m))))**2))
                    op_scale = sqrt(sum(abs(op_trial(:,:m))**2))
                    if (.not. ieee_is_finite(leakage)) cycle
                    if (leakage > 5.0_dp*tolerance*max(1.0_dp,op_scale)) cycle
                    nlocked = m
                    locked(:,nlocked) = trial(:,nlocked)
                end do
                if (nlocked > previous_locks) then
                    ! Unitary Schur rotation of the locked invariant subspace.
                    ! This produces an orthonormal basis with a (numerically)
                    ! upper-triangular projected operator; no left-eigenvector
                    ! deflation or invalid right-eigenvector projection is used.
                    allocate(lock_matrix(nlocked,nlocked),lock_basis(nlocked,nlocked), &
                        lock_vectors(nlocked,nlocked),lock_values(nlocked),new_locked(n,nlocked))
                    do l = 1,nlocked
                        if (shifted) then
                            call complex_shift_gmres(a,sigma,locked(:,l),inverse,status, &
                                solve_tol,inner_cap,inner_restart_size)
                            if (status /= sparse_linalg_success) then
                                info = sparse_linalg_breakdown
                                exit
                            end if
                            op_trial(:,l) = inverse
                        else
                            call csr_complex_product(a,locked(:,l),op_trial(:,l))
                        end if
                        products = products+1
                    end do
                    if (info == sparse_linalg_breakdown) exit
                    lock_matrix = matmul(conjg(transpose(locked(:,:nlocked))),op_trial(:,:nlocked))
                    call complex_schur_eigensystem(lock_matrix,lock_values,lock_vectors,status, &
                        schur_vectors=lock_basis)
                    if (status /= sparse_linalg_success) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    new_locked = matmul(locked(:,:nlocked),lock_basis)
                    locked(:,:nlocked) = new_locked
                    deallocate(lock_matrix,lock_basis,lock_vectors,lock_values,new_locked)
                end if
                ! Retain the entire invariant block and additional unlocked
                ! wanted/guard Ritz vectors.  Modified Gram--Schmidt deflates
                ! each new search direction against the frozen Schur vectors.
                if (block_width == 1) then
                    nr = min(p-1,max(k+1,nlocked+1))
                else
                    nr = min(p-block_width-1,max(k,nlocked))
                end if
                retained(:,1:nlocked) = locked(:,1:nlocked)
                tmp = nlocked
                do j = 1,p
                    if (tmp >= nr) exit
                    candidate = matmul(q,y(:,order(j)))
                    do pass = 1,2
                        do l = 1,tmp
                            candidate = candidate-dot_product(retained(:,l),candidate)*retained(:,l)
                        end do
                    end do
                    norm_v = complex_norm(candidate)
                    if (norm_v <= 1.0e-8_dp) cycle
                    tmp = tmp+1
                    retained(:,tmp) = candidate/norm_v
                end do
                do while (tmp < nr)
                    tmp = tmp+1
                    retained(:,tmp) = cmplx(0.0_dp,0.0_dp,kind=dp)
                    retained(mod(tmp-1,n)+1,tmp) = cmplx(1.0_dp,0.0_dp,kind=dp)
                    call arnoldi_orthogonalize(retained(:,:tmp),tmp)
                end do
                q(:,1:nr) = retained(:,1:nr)
            else
            ! Soft lock: keep already converged desired Ritz vectors ahead
            ! of other desired and guard Ritz vectors.  Preserve more than k
            ! vectors so a cluster of eigenvalues does not get discarded.
            if (block_width == 1) then
                nr = min(p-1,k+1)
            else
                nr = min(p-block_width-1,k)
            end if
            tmp = 0
            do j = 1,k
                if (errors(j) <= tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) then
                    tmp = tmp+1
                    retained(:,tmp) = cand(:,j)
                end if
            end do
            do j = 1,p
                if (tmp >= nr) exit
                if (j <= k) then
                    if (errors(j) <= tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) cycle
                end if
                tmp = tmp+1
                retained(:,tmp) = matmul(q,y(:,order(j)))
            end do
            call orthonormal_complex(retained(:,:nr))
            q(:,1:nr) = retained(:,:nr)
            end if
            ! Restart with independent Ritz residuals, not their sum: the
            ! latter is a single Krylov direction and can miss degeneracies.
            ! Complete a short residual block with deterministic guard seeds
            ! when fewer than b wanted eigenpairs are unconverged.
            if (block_width == 1) then
                work = cmplx(0.0_dp,0.0_dp,kind=dp)
                do j = 1,k
                    if (errors(j) > tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) &
                        work = work+residual(:,j)/max(errors(j),tiny(1.0_dp))
                end do
                q(:,nr+1) = work
                call arnoldi_orthogonalize(q,nr+1)
                nr = nr+1
            else
                if (nr < k) exit ! A full-size search space has no valid restart.
                retained_count = nr
                seed_count = 0
                do j = 1,k
                    if (seed_count >= block_width) exit
                    if (errors(j) <= tolerance*max(1.0_dp,norm_a,abs(eigenvalues(j)))) cycle
                    seed_count = seed_count+1
                    q(:,retained_count+seed_count) = residual(:,j)
                    call arnoldi_orthogonalize(q,retained_count+seed_count)
                end do
                do while (seed_count < block_width)
                    seed_count = seed_count+1
                    do i = 1,n
                        q(i,retained_count+seed_count) = cmplx( &
                            sin(real(i+cycle+seed_count,dp)*0.337123_dp*real(seed_count+2,dp)), &
                            cos(real(2*i+cycle+seed_count,dp)*0.618034_dp),kind=dp)
                    end do
                    call arnoldi_orthogonalize(q,retained_count+seed_count)
                end do
                nr = retained_count+block_width
            end if
        end do
        if (present(iterations)) iterations = min(cycle,cap)
        if (present(locked_count)) locked_count = count
        if (present(schur_locked_count)) schur_locked_count = nlocked
        if (present(matvecs)) matvecs = products
        if (present(residual_norms)) residual_norms = errors
    end subroutine eigs_arnoldi

    ! Two-pass modified Gram--Schmidt with deterministic complement fallback
    ! when an Arnoldi direction breaks down.  The fallback never relies on a
    ! global random-number generator.
    subroutine arnoldi_orthogonalize(q, j)
        complex(dp), intent(inout) :: q(:,:) !! Basis; normalize column j against 1..j-1.
        integer, intent(in) :: j !! One-based column to orthogonalize.
        integer :: i, pass, trial, n
        real(dp) :: magnitude
        n = size(q,1)
        do trial = 0,n
            do pass = 1,2
                do i = 1,j-1
                    q(:,j) = q(:,j)-dot_product(q(:,i),q(:,j))*q(:,i)
                end do
            end do
            magnitude = complex_norm(q(:,j))
            if (magnitude > 100.0_dp*epsilon(1.0_dp)) then
                q(:,j) = q(:,j)/magnitude
                return
            end if
            q(:,j) = cmplx(0.0_dp,0.0_dp,kind=dp)
            q(mod(j+trial-1,n)+1,j) = cmplx(1.0_dp,0.0_dp,kind=dp)
        end do
        q(:,j) = cmplx(0.0_dp,0.0_dp,kind=dp)
    end subroutine arnoldi_orthogonalize

    ! Restarted, complex-valued GMRES for (A-sigma*I)x=b.  It operates on
    ! a real CSR input but accepts complex shifts, right-hand sides and
    ! solutions, without constructing a shifted CSR matrix or dense inverse.
    ! The inner iteration count is the total Arnoldi-step budget.
    subroutine complex_shift_gmres(a, sigma, b, x, info, rtol, maxiter, restart)
        type(csr_matrix), intent(in) :: a !! Real square canonical CSR operator.
        complex(dp), intent(in) :: sigma !! Complex shift for the linear system.
        complex(dp), intent(in) :: b(:) !! Complex right-hand side of length n.
        complex(dp), allocatable, intent(out) :: x(:) !! Complex solution of length n.
        integer, intent(out) :: info !! 0 if residual passed, -2 on failure.
        real(dp), intent(in) :: rtol !! Positive relative true-residual tolerance.
        integer, intent(in) :: maxiter !! Total number of Arnoldi iterations available.
        integer, intent(in) :: restart !! Maximum Arnoldi steps per restart cycle.
        complex(dp), allocatable :: v(:,:), h(:,:), s(:), g(:), y(:), &
            w(:), r(:), candidate(:)
        complex(dp) :: alpha, old, diagonal
        real(dp), allocatable :: c(:)
        real(dp) :: beta, threshold, a_abs, b_abs, hyp, err
        integer :: n, m, i, j, pass, cycle_steps, used, steps
        logical :: terminated

        n = size(b)
        allocate(x(n))
        x = cmplx(0.0_dp,0.0_dp,kind=dp)
        info = sparse_linalg_breakdown
        beta = complex_norm(b)
        if (.not. ieee_is_finite(beta)) return
        if (beta == 0.0_dp) then
            info = sparse_linalg_success
            return
        end if
        threshold = rtol*beta
        m = min(n,restart)
        allocate(v(n,m+1),h(m+1,m),s(m),c(m),g(m+1),y(m), &
            w(n),r(n),candidate(n))
        r = b
        used = 0
        do while (used < maxiter)
            beta = complex_norm(r)
            if (beta <= threshold) then
                info = sparse_linalg_success
                exit
            end if
            if (beta <= tiny(beta) .or. .not. ieee_is_finite(beta)) exit
            v(:,1) = r/beta
            h = cmplx(0.0_dp,0.0_dp,kind=dp)
            g = cmplx(0.0_dp,0.0_dp,kind=dp)
            g(1) = cmplx(beta,0.0_dp,kind=dp)
            steps = min(m,maxiter-used)
            cycle_steps = 0
            terminated = .false.
            do j = 1,steps
                cycle_steps = j
                call csr_complex_product(a,v(:,j),w)
                w = w-sigma*v(:,j)
                do pass = 1,2
                    do i = 1,j
                        old = dot_product(v(:,i),w)
                        h(i,j) = h(i,j)+old
                        w = w-old*v(:,i)
                    end do
                end do
                h(j+1,j) = cmplx(complex_norm(w),0.0_dp,kind=dp)
                if (.not. ieee_is_finite(real(h(j+1,j),dp))) return
                if (real(h(j+1,j),dp) > 0.0_dp) v(:,j+1) = w/real(h(j+1,j),dp)
                do i = 1,j-1
                    old = c(i)*h(i,j)+s(i)*h(i+1,j)
                    h(i+1,j) = -conjg(s(i))*h(i,j)+c(i)*h(i+1,j)
                    h(i,j) = old
                end do
                ! Stable complex Givens: [c,s;-conjg(s),c] zeroes h(j+1,j).
                a_abs = abs(h(j,j))
                b_abs = abs(h(j+1,j))
                hyp = hypot(a_abs,b_abs)
                if (hyp == 0.0_dp .or. .not. ieee_is_finite(hyp)) return
                c(j) = a_abs/hyp
                alpha = cmplx(1.0_dp,0.0_dp,kind=dp)
                if (a_abs /= 0.0_dp) alpha = h(j,j)/a_abs
                s(j) = alpha*conjg(h(j+1,j))/hyp
                h(j,j) = c(j)*h(j,j)+s(j)*h(j+1,j)
                h(j+1,j) = cmplx(0.0_dp,0.0_dp,kind=dp)
                g(j+1) = -conjg(s(j))*g(j)
                g(j) = c(j)*g(j)
                if (abs(g(j+1)) <= threshold .or. b_abs <= tiny(b_abs)) then
                    terminated = .true.
                    exit
                end if
            end do
            y(:cycle_steps) = g(:cycle_steps)
            do i = cycle_steps,1,-1
                diagonal = h(i,i)
                if (abs(diagonal) <= tiny(1.0_dp)) return
                y(i) = y(i)/diagonal
                if (i > 1) y(:i-1) = y(:i-1)-h(:i-1,i)*y(i)
            end do
            candidate = x+matmul(v(:,:cycle_steps),y(:cycle_steps))
            call csr_complex_product(a,candidate,w)
            r = b-w+sigma*candidate
            err = complex_norm(r)
            if (.not. ieee_is_finite(err)) return
            x = candidate
            used = used+cycle_steps
            if (err <= threshold) then
                info = sparse_linalg_success
                exit
            end if
            if (terminated .and. cycle_steps < steps) then
                ! Early stopping was based on the projected residual, but
                ! its true residual did not meet the user's tolerance.
                cycle
            end if
        end do
    end subroutine complex_shift_gmres

    pure subroutine csr_complex_product(a, x, y)
        type(csr_matrix), intent(in) :: a !! Validated real CSR matrix.
        complex(dp), intent(in) :: x(:) !! Complex input vector.
        complex(dp), intent(out) :: y(:) !! Complex product vector.
        integer :: i, t
        y = cmplx(0.0_dp,0.0_dp,kind=dp)
        do i = 1,a%n_rows
            do t = a%indptr(i),a%indptr(i+1)-1
                y(i) = y(i)+a%data(t)*x(a%indices(t))
            end do
        end do
    end subroutine csr_complex_product

    subroutine orthonormal_complex(q)
        complex(dp), intent(inout) :: q(:,:) !! Complex n-by-p basis, orthonormalized in place.
        real(dp) :: magnitude
        integer :: i,j,attempt,pass,n
        n = size(q,1)
        do j = 1,size(q,2)
            do attempt = 0,n
                do pass = 1,2
                    do i = 1,j-1
                        q(:,j) = q(:,j)-dot_product(q(:,i),q(:,j))*q(:,i)
                    end do
                end do
                magnitude = complex_norm(q(:,j))
                if (magnitude > 100.0_dp*epsilon(1.0_dp)) exit
                q(:,j) = cmplx(0.0_dp,0.0_dp,kind=dp)
                q(mod(j+attempt-1,n)+1,j) = cmplx(1.0_dp,0.0_dp,kind=dp)
            end do
            if (magnitude > 100.0_dp*epsilon(1.0_dp)) then
                q(:,j) = q(:,j)/magnitude
            else
                q(:,j) = cmplx(0.0_dp,0.0_dp,kind=dp)
            end if
        end do
    end subroutine orthonormal_complex

    pure real(dp) function nonsymmetric_priority(value, target) result(priority)
        complex(dp), intent(in) :: value !! Complex eigenvalue.
        character(len=2), intent(in) :: target !! Mode LM, SM, LR, SR, LI, SI.
        select case(target)
        case('LM')
            priority = abs(value)
        case('SM')
            priority = -abs(value)
        case('LR')
            priority = real(value,dp)
        case('SR')
            priority = -real(value,dp)
        case('LI')
            priority = aimag(value)
        case('SI')
            priority = -aimag(value)
        end select
    end function nonsymmetric_priority

    ! Unit Householder reflector I-2*v*v^H that maps x to an e_1 multiple.
    subroutine complex_reflector(x, v)
        complex(dp), intent(in) :: x(:) !! Column segment to be eliminated.
        complex(dp), intent(out) :: v(:) !! Unit Householder direction, zero when no operation.
        complex(dp) :: phase
        real(dp) :: magnitude, denominator
        v = cmplx(0.0_dp,0.0_dp,kind=dp)
        magnitude = complex_norm(x)
        if (magnitude <= tiny(magnitude)) return
        phase = cmplx(1.0_dp,0.0_dp,kind=dp)
        if (abs(x(1)) > 0.0_dp) phase = x(1)/abs(x(1))
        v = x
        v(1) = v(1)+phase*magnitude
        denominator = complex_norm(v)
        if (denominator > 0.0_dp) v = v/denominator
    end subroutine complex_reflector

    ! Form U (unitary), R (upper triangular) with B=U*R, from Householder QR.
    subroutine complex_qr(b, u, r)
        complex(dp), intent(in) :: b(:,:) !! Square matrix to factorize.
        complex(dp), intent(out) :: u(:,:) !! Unitary Q factor.
        complex(dp), intent(out) :: r(:,:) !! Upper triangular R factor.
        complex(dp) :: v(size(b,1)), scalar
        integer :: n, j, col, row
        n = size(b,1)
        u = cmplx(0.0_dp,0.0_dp,kind=dp)
        do j = 1,n
            u(j,j) = cmplx(1.0_dp,0.0_dp,kind=dp)
        end do
        r = b
        do j = 1,n-1
            call complex_reflector(r(j:n,j),v(j:n))
            if (complex_norm(v(j:n)) == 0.0_dp) cycle
            do col = j,n
                scalar = 2.0_dp*dot_product(v(j:n),r(j:n,col))
                r(j:n,col) = r(j:n,col)-v(j:n)*scalar
            end do
            do row = 1,n
                scalar = 2.0_dp*sum(u(row,j:n)*v(j:n))
                u(row,j:n) = u(row,j:n)-scalar*conjg(v(j:n))
            end do
        end do
    end subroutine complex_qr

    ! Complex Schur triangularization with Wilkinson-shift QR; eigenvectors
    ! are recovered from triangular back-substitution.  p is intentionally
    ! small (the subspace dimension), and no external LAPACK is required.
    subroutine complex_schur_eigensystem(h, values, vectors, info, schur_vectors)
        complex(dp), intent(in) :: h(:,:) !! Small complex square Rayleigh--Ritz matrix.
        complex(dp), intent(out) :: values(:) !! Eigenvalues (unsorted).
        complex(dp), intent(out) :: vectors(:,:) !! Right complex unit eigenvectors.
        integer, intent(out) :: info !! Zero success; -2 on QR failure.
        complex(dp), intent(out), optional :: schur_vectors(:,:) !! Unitary Schur basis of projected matrix.
        complex(dp), allocatable :: t(:,:), z(:,:), u(:,:), r(:,:), &
            block(:,:), updated(:,:), col(:), v(:), y(:), tmp(:,:)
        complex(dp) :: alpha, delta, beta, gamma, mean, disc, mu, scalar
        real(dp) :: threshold, magnitude
        integer :: n, i,j,l,m,iter
        n = size(h,1)
        info = sparse_linalg_breakdown
        allocate(t(n,n),z(n,n),u(n,n),r(n,n),block(n,n),updated(n,n), &
            col(n),v(n),y(n),tmp(n,n))
        t = h
        z = cmplx(0.0_dp,0.0_dp,kind=dp)
        do i=1,n
            z(i,i)=cmplx(1.0_dp,0.0_dp,kind=dp)
        end do
        ! Reduce to upper Hessenberg form, preserving the unitary similarity.
        do j=1,n-2
            call complex_reflector(t(j+1:n,j),v(j+1:n))
            if (complex_norm(v(j+1:n)) == 0.0_dp) cycle
            do l=j,n
                scalar=2.0_dp*dot_product(v(j+1:n),t(j+1:n,l))
                t(j+1:n,l)=t(j+1:n,l)-v(j+1:n)*scalar
            end do
            do i=1,n
                scalar=2.0_dp*sum(t(i,j+1:n)*v(j+1:n))
                t(i,j+1:n)=t(i,j+1:n)-scalar*conjg(v(j+1:n))
                scalar=2.0_dp*sum(z(i,j+1:n)*v(j+1:n))
                z(i,j+1:n)=z(i,j+1:n)-scalar*conjg(v(j+1:n))
            end do
            t(j+2:n,j)=cmplx(0.0_dp,0.0_dp,kind=dp)
        end do
        m=n
        iter=0
        do while(m>1)
            threshold=32.0_dp*epsilon(1.0_dp)* &
                (abs(t(m-1,m-1))+abs(t(m,m))+1.0_dp)
            if(abs(t(m,m-1))<=threshold) then
                t(m,m-1)=cmplx(0.0_dp,0.0_dp,kind=dp)
                m=m-1
                cycle
            end if
            iter=iter+1
            if(iter>max(200,200*n)) return
            alpha=t(m-1,m-1)
            beta=t(m-1,m)
            gamma=t(m,m-1)
            delta=t(m,m)
            mean=0.5_dp*(alpha+delta)
            disc=sqrt((0.5_dp*(alpha-delta))**2+beta*gamma)
            mu=mean+disc
            if(abs(mean-disc-delta)<abs(mu-delta)) mu=mean-disc
            block(:m,:m)=t(:m,:m)
            do i=1,m
                block(i,i)=block(i,i)-mu
            end do
            call complex_qr(block(:m,:m),u(:m,:m),r(:m,:m))
            t(:m,:m)=matmul(r(:m,:m),u(:m,:m))
            do i=1,m
                t(i,i)=t(i,i)+mu
            end do
            if(m<n) then
                updated(:m,m+1:n)=matmul(conjg(transpose(u(:m,:m))),t(:m,m+1:n))
                t(:m,m+1:n)=updated(:m,m+1:n)
            end if
            tmp(:,1:m)=matmul(z(:,1:m),u(:m,:m))
            z(:,1:m)=tmp(:,1:m)
            ! Remove harmless roundoff beneath the first subdiagonal.
            do j=1,m-2
                t(j+2:m,j)=cmplx(0.0_dp,0.0_dp,kind=dp)
            end do
            if (.not. all(ieee_is_finite(real(t,dp))) .or. &
                    .not. all(ieee_is_finite(aimag(t)))) return
        end do
        do i=1,n
            values(i)=t(i,i)
            y=cmplx(0.0_dp,0.0_dp,kind=dp)
            y(i)=cmplx(1.0_dp,0.0_dp,kind=dp)
            do j=i-1,1,-1
                scalar=t(j,j)-values(i)
                threshold=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(values(i)),abs(t(j,j)))
                if(abs(scalar)<threshold) scalar=cmplx(threshold,0.0_dp,kind=dp)
                y(j)=-sum(t(j,j+1:i)*y(j+1:i))/scalar
                if(abs(y(j))>1.0e100_dp) y(1:i)=y(1:i)*1.0e-100_dp
            end do
            col=matmul(z,y)
            magnitude=complex_norm(col)
            if(magnitude<=tiny(magnitude) .or. .not. ieee_is_finite(magnitude)) return
            vectors(:,i)=col/magnitude
        end do
        if (present(schur_vectors)) schur_vectors = z
        info=sparse_linalg_success
    end subroutine complex_schur_eigensystem

    pure real(dp) function complex_norm(x) result(magnitude)
        complex(dp), intent(in) :: x(:) !! Complex vector whose Euclidean norm is needed.
        real(dp) :: scale
        scale = maxval(abs(x))
        magnitude = 0.0_dp
        if (scale > 0.0_dp) magnitude = scale*sqrt(sum((abs(x)/scale)**2))
    end function complex_norm

    pure real(dp) function symmetric_shift_priority(value, sigma, target) result(priority)
        real(dp), intent(in) :: value !! Rayleigh--Ritz value of original A.
        real(dp), intent(in) :: sigma !! Real shift of inverse operator.
        character(len=2), intent(in) :: target !! LM, SM, LA, or SA, applied to inverse eigenvalues.
        real(dp) :: delta
        delta = value-sigma
        select case(target)
        case('LM')
            priority = -abs(delta)
        case('SM')
            priority = abs(delta)
        case('LA')
            if (delta == 0.0_dp) then
                priority = huge(1.0_dp)
            else
                priority = 1.0_dp/delta
            end if
        case('SA')
            if (delta == 0.0_dp) then
                priority = -huge(1.0_dp)
            else
                priority = -1.0_dp/delta
            end if
        end select
    end function symmetric_shift_priority

    pure real(dp) function spectral_priority(value, target) result(priority)
        real(dp), intent(in) :: value !! Eigenvalue of original symmetric operator.
        character(len=2), intent(in) :: target !! One of LA, SA, LM, SM.
        select case(target)
        case('LA')
            priority = value
        case('SA')
            priority = -value
        case('LM')
            priority = abs(value)
        case('SM')
            priority = -abs(value)
        end select
    end function spectral_priority

    subroutine orthonormal_columns(q)
        real(dp), intent(inout) :: q(:,:) !! n-by-p columns replaced by orthonormal vectors.
        real(dp) :: mag
        integer :: i, j, attempt, pass, n
        n = size(q,1)
        do j = 1,size(q,2)
            do attempt = 0,n
                do pass = 1,2
                    do i = 1,j-1
                        q(:,j) = q(:,j)-dot_product(q(:,i),q(:,j))*q(:,i)
                    end do
                end do
                mag = norm2(q(:,j))
                if (mag > 100.0_dp*epsilon(1.0_dp)) exit
                q(:,j) = 0.0_dp
                q(mod(j+attempt-1,n)+1,j) = 1.0_dp
            end do
            if (mag > 100.0_dp*epsilon(1.0_dp)) then
                q(:,j) = q(:,j)/mag
            else
                q(:,j) = 0.0_dp
            end if
        end do
    end subroutine orthonormal_columns

    subroutine symmetric_jacobi_small(a, eigenvalues, eigenvectors, info)
        real(dp), intent(inout) :: a(:,:) !! Small dense symmetric matrix, overwritten.
        real(dp), intent(out) :: eigenvalues(:) !! Eigenvalues of original matrix.
        real(dp), intent(out) :: eigenvectors(:,:) !! Orthonormal eigenvector columns.
        integer, intent(out) :: info !! Zero on convergence, -2 on failure.
        real(dp) :: off, max_diag, tau, t, cs, sn, app, aqq, apq, r, s
        integer :: n, i, j, p, q, sweep
        n = size(a,1)
        eigenvectors = 0.0_dp
        do i = 1,n
            eigenvectors(i,i) = 1.0_dp
        end do
        info = sparse_linalg_breakdown
        do sweep = 1,40*n*n
            off = 0.0_dp
            p = 1
            q = 1
            do i = 1,n-1
                do j = i+1,n
                    if (abs(a(i,j)) > off) then
                        off = abs(a(i,j))
                        p = i
                        q = j
                    end if
                end do
            end do
            max_diag = 1.0_dp
            do i = 1,n
                max_diag = max(max_diag,abs(a(i,i)))
            end do
            if (off <= 10.0_dp*epsilon(1.0_dp)*max_diag) then
                info = sparse_linalg_success
                exit
            end if
            app = a(p,p)
            aqq = a(q,q)
            apq = a(p,q)
            tau = (aqq-app)/(2.0_dp*apq)
            t = sign(1.0_dp,tau)/(abs(tau)+hypot(1.0_dp,tau))
            cs = 1.0_dp/sqrt(1.0_dp+t*t)
            sn = t*cs
            a(p,p) = app-t*apq
            a(q,q) = aqq+t*apq
            a(p,q) = 0.0_dp
            a(q,p) = 0.0_dp
            do j = 1,n
                if (j == p .or. j == q) cycle
                r = a(j,p)
                s = a(j,q)
                a(j,p) = cs*r-sn*s
                a(p,j) = a(j,p)
                a(j,q) = sn*r+cs*s
                a(q,j) = a(j,q)
            end do
            do j = 1,n
                r = eigenvectors(j,p)
                s = eigenvectors(j,q)
                eigenvectors(j,p) = cs*r-sn*s
                eigenvectors(j,q) = sn*r+cs*s
            end do
        end do
        do i = 1,n
            eigenvalues(i) = a(i,i)
        end do
    end subroutine symmetric_jacobi_small

    ! Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG).
    ! The Rayleigh--Ritz search space is span{X, M*R, P}, where P stores
    ! conjugate search directions.  Columns are explicitly B-orthogonalized
    ! (B=I when omitted); only the at-most-3*k projected eigenproblem is dense.
    ! Independent implementation of the Knyazev (2001) method. No BLAS.
    subroutine lobpcg(a, x0, eigenvalues, eigenvectors, info, largest, tol, &
            maxiter, b, m, y, jacobi, iterations, residual_norms, &
            check_symmetry)
        type(csr_matrix), intent(in) :: a !! Real symmetric n-by-n CSR operator.
        real(dp), intent(in) :: x0(:,:) !! Initial n-by-k linearly independent vectors.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! k Ritz values, descending if largest.
        real(dp), allocatable, intent(out) :: eigenvectors(:,:) !! n-by-k B-orthonormal eigenvectors.
        integer, intent(out) :: info !! 0 converged; positive iteration cap; negative invalid/breakdown.
        logical, intent(in), optional :: largest !! Largest eigenvalues by default; false selects smallest.
        real(dp), intent(in), optional :: tol !! Absolute eigen-residual tolerance; default n*sqrt(epsilon).
        integer, intent(in), optional :: maxiter !! Block iterations, default 20.
        type(csr_matrix), intent(in), optional :: b !! Optional symmetric positive-definite B (A*x=lambda*B*x).
        type(csr_matrix), intent(in), optional :: m !! Optional symmetric positive-definite preconditioner action.
        real(dp), intent(in), optional :: y(:,:) !! Constraint directions, B-orthogonal complement.
        logical, intent(in), optional :: jacobi !! Positive diagonal preconditioner; exclusive with m.
        integer, intent(out), optional :: iterations !! Number of completed block iterations.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! ||A*x-lambda*B*x||; projected when Y is supplied.
        logical, intent(in), optional :: check_symmetry !! Validate A and B symmetry, default true.
        real(dp), allocatable :: x(:,:), p(:,:), w(:,:), r(:,:), ax(:,:), bx(:,:), &
            s(:,:), bs(:,:), as(:,:), qy(:,:), bqy(:,:), empty(:,:), &
            h(:,:), small_vec(:,:), small_val(:), residuals(:), best_res(:), &
            best_x(:,:), best_val(:), invdiag(:)
        real(dp) :: threshold, best_error, current_error, diag
        integer, allocatable :: permutation(:)
        integer :: n, k, ny, i, j, pos, count, nconstraints, cap, step, ierr
        logical :: pick_largest, use_jacobi, verify_symmetry, accepted

        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        allocate(eigenvalues(0), eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr(a)) return
        n = a%n_rows
        k = size(x0,2)
        ny = 0
        if (present(y)) ny = size(y,2)
        if (n /= a%n_cols .or. n <= 0 .or. k < 1 .or. k > n-ny) return
        if (k > huge(k)/3) return
        if (size(x0,1) /= n .or. .not. all(ieee_is_finite(x0))) return
        if (present(y)) then
            if (size(y,1) /= n .or. .not. all(ieee_is_finite(y))) return
        end if
        if (present(b)) then
            if (.not. valid_csr(b)) return
            if (b%n_rows /= n .or. b%n_cols /= n) return
        end if
        if (present(m)) then
            if (.not. valid_csr(m)) return
            if (m%n_rows /= n .or. m%n_cols /= n) return
        end if
        use_jacobi = .false.
        if (present(jacobi)) use_jacobi = jacobi
        if (use_jacobi .and. present(m)) return
        verify_symmetry = .true.
        if (present(check_symmetry)) verify_symmetry = check_symmetry
        if (verify_symmetry) then
            if (.not. symmetric_csr(a)) return
            if (present(b)) then
                if (.not. symmetric_csr(b)) return
            end if
        end if
        threshold = real(n,dp)*sqrt(epsilon(1.0_dp))
        if (present(tol)) threshold = tol
        if (.not. ieee_is_finite(threshold) .or. threshold < 0.0_dp) return
        if (threshold == 0.0_dp) threshold = 10.0_dp*epsilon(1.0_dp)
        cap = 20
        if (present(maxiter)) cap = maxiter
        if (cap < 0) return
        pick_largest = .true.
        if (present(largest)) pick_largest = largest
        allocate(invdiag(n))
        invdiag = 1.0_dp
        if (use_jacobi) then
            do i=1,n
                pos = row_position(a,i,i)
                if (pos == 0) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
                diag = abs(a%data(pos))
                if (diag <= tiny(diag)) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
                invdiag(i) = 1.0_dp/diag
                if (.not. ieee_is_finite(invdiag(i))) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
            end do
        end if
        allocate(x(n,k),p(n,k),w(n,k),r(n,k),ax(n,k),bx(n,k), &
            s(n,3*k),bs(n,3*k),as(n,3*k),qy(n,ny),bqy(n,ny), &
            empty(n,0),h(3*k,3*k),small_vec(3*k,3*k),small_val(3*k), &
            permutation(3*k),residuals(k),best_res(k),best_x(n,k),best_val(k))
        qy = 0.0_dp
        bqy = 0.0_dp
        s = 0.0_dp
        bs = 0.0_dp
        p = 0.0_dp
        count = 0
        nconstraints = 0
        if (present(y)) then
            do j=1,ny
                call append_direction(y(:,j),qy,bqy,nconstraints,empty,empty,0,accepted)
                if (.not. accepted) return
            end do
        end if
        do j=1,k
            call append_direction(x0(:,j),s,bs,count,qy,bqy,nconstraints,accepted)
            if (.not. accepted) return
        end do
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues=0.0_dp
        eigenvectors=0.0_dp
        ! Diagonalize the initial k-dimensional Rayleigh quotient.  The
        ! iteration starts with genuine Ritz vectors and residuals.
        call update_ritz(count,.false.,ierr)
        if (ierr /= 0) then
            info = sparse_linalg_breakdown
            return
        end if
        call compute_residuals()
        if (.not. all(ieee_is_finite(residuals))) then
            info = sparse_linalg_breakdown
            return
        end if
        best_error = maxval(residuals)
        best_val = eigenvalues
        best_x = x
        best_res = residuals
        info = max(1,cap)
        if (best_error <= threshold) info = sparse_linalg_success
        do step=1,cap
            if (info == sparse_linalg_success) exit
            if (present(m)) then
                do j=1,k
                    call csr_product(m,r(:,j),w(:,j))
                end do
            else
                w = r
                if (use_jacobi) then
                    do j=1,k
                        w(:,j) = w(:,j)*invdiag
                    end do
                end if
            end if
            if (.not. all(ieee_is_finite(w))) then
                info = sparse_linalg_breakdown
                exit
            end if
            count = 0
            do j=1,k
                call append_direction(x(:,j),s,bs,count,qy,bqy,nconstraints,accepted)
                if (.not. accepted) exit
            end do
            if (count /= k) then
                info = sparse_linalg_breakdown
                exit
            end if
            do j=1,k
                call append_direction(w(:,j),s,bs,count,qy,bqy,nconstraints,accepted)
            end do
            do j=1,k
                call append_direction(p(:,j),s,bs,count,qy,bqy,nconstraints,accepted)
            end do
            if (count == k) then
                ! No independent correction vector: report honest stagnation.
                exit
            end if
            call update_ritz(count,.true.,ierr)
            if (ierr /= 0) then
                info = sparse_linalg_breakdown
                exit
            end if
            call compute_residuals()
            current_error = maxval(residuals)
            if (.not. ieee_is_finite(current_error)) then
                info = sparse_linalg_breakdown
                exit
            end if
            if (current_error < best_error) then
                best_error = current_error
                best_val = eigenvalues
                best_x = x
                best_res = residuals
            end if
            if (present(iterations)) iterations = step
            if (current_error <= threshold) then
                info = sparse_linalg_success
                exit
            end if
        end do
        eigenvalues = best_val
        eigenvectors = best_x
        if (present(residual_norms)) then
            deallocate(residual_norms)
            allocate(residual_norms(k),source=best_res)
        end if

    contains

        subroutine apply_b(v,bv)
            real(dp), intent(in) :: v(:) !! Input length n.
            real(dp), intent(out) :: bv(:) !! Output B*v, or v if B absent.
            if (present(b)) then
                call csr_product(b,v,bv)
            else
                bv = v
            end if
        end subroutine apply_b

        logical function symmetric_csr(mat) result(ok)
            type(csr_matrix), intent(in) :: mat !! Canonical square CSR matrix.
            integer :: row, entry, column, other
            real(dp) :: opposite
            ok = .false.
            do row=1,n
                do entry=mat%indptr(row),mat%indptr(row+1)-1
                    column=mat%indices(entry)
                    other=row_position(mat,column,row)
                    opposite=0.0_dp
                    if (other > 0) opposite=mat%data(other)
                    if (abs(mat%data(entry)-opposite) > 1.0e-12_dp* &
                            max(1.0_dp,abs(mat%data(entry)),abs(opposite))) return
                end do
            end do
            ok = .true.
        end function symmetric_csr

        subroutine append_direction(vin, basis, bbasis, used, guard, bguard, ng, ok)
            real(dp), intent(in) :: vin(:) !! Candidate length n.
            real(dp), intent(inout) :: basis(:,:) !! Accumulating B-orthonormal directions.
            real(dp), intent(inout) :: bbasis(:,:) !! Matching B*basis products.
            integer, intent(inout) :: used !! Current basis size.
            real(dp), intent(in) :: guard(:,:) !! B-orthonormal constraints.
            real(dp), intent(in) :: bguard(:,:) !! B*guard.
            integer, intent(in) :: ng !! Number of constraints.
            logical, intent(out) :: ok !! True if direction is independent.
            real(dp) :: v(size(vin)), bv(size(vin)), initial, magnitude, coeff
            integer :: pass, j
            ok = .false.
            if (used >= size(basis,2)) return
            v = vin
            call apply_b(v,bv)
            if (.not. all(ieee_is_finite(bv))) return
            initial = dot_product(v,bv)
            if (initial <= 0.0_dp .or. .not. ieee_is_finite(initial)) return
            do pass=1,2
                do j=1,ng
                    coeff=dot_product(bguard(:,j),v)
                    v=v-coeff*guard(:,j)
                    bv=bv-coeff*bguard(:,j)
                end do
                do j=1,used
                    coeff=dot_product(bbasis(:,j),v)
                    v=v-coeff*basis(:,j)
                    bv=bv-coeff*bbasis(:,j)
                end do
            end do
            magnitude = dot_product(v,bv)
            if (.not. ieee_is_finite(magnitude)) return
            if (magnitude <= max(100.0_dp*epsilon(1.0_dp)*initial, &
                    tiny(magnitude))) return
            magnitude = sqrt(magnitude)
            used=used+1
            basis(:,used)=v/magnitude
            bbasis(:,used)=bv/magnitude
            ok = .true.
        end subroutine append_direction

        subroutine update_ritz(ns, retain_p, error_code)
            integer, intent(in) :: ns !! Current B-orthonormal subspace dimension.
            logical, intent(in) :: retain_p !! Preserve W/P search corrections.
            integer, intent(out) :: error_code !! Zero on success.
            integer :: aidx, bidx, choice, swap, iidx, l
            error_code=sparse_linalg_breakdown
            do aidx=1,ns
                call csr_product(a,s(:,aidx),as(:,aidx))
            end do
            h(1:ns,1:ns)=matmul(transpose(s(:,1:ns)),as(:,1:ns))
            h(1:ns,1:ns)=0.5_dp*(h(1:ns,1:ns)+transpose(h(1:ns,1:ns)))
            if (.not. all(ieee_is_finite(h(1:ns,1:ns)))) return
            call symmetric_jacobi_small(h(1:ns,1:ns),small_val(1:ns), &
                small_vec(1:ns,1:ns),error_code)
            if (error_code /= 0) return
            do aidx=1,ns
                permutation(aidx)=aidx
            end do
            ! Largest algebraic eigenvalues first, or smallest first.
            do aidx=1,k
                choice=aidx
                do bidx=aidx+1,ns
                    iidx=permutation(bidx)
                    l=permutation(choice)
                    if (pick_largest) then
                        if (small_val(iidx)>small_val(l)) choice=bidx
                    else
                        if (small_val(iidx)<small_val(l)) choice=bidx
                    end if
                end do
                swap=permutation(aidx)
                permutation(aidx)=permutation(choice)
                permutation(choice)=swap
            end do
            do aidx=1,k
                iidx=permutation(aidx)
                eigenvalues(aidx)=small_val(iidx)
                x(:,aidx)=matmul(s(:,1:ns),small_vec(1:ns,iidx))
                p(:,aidx)=0.0_dp
                if (retain_p .and. ns>k) then
                    p(:,aidx)=matmul(s(:,k+1:ns),small_vec(k+1:ns,iidx))
                end if
            end do
            error_code=sparse_linalg_success
        end subroutine update_ritz

        subroutine compute_residuals()
            integer :: j
            do j=1,k
                call csr_product(a,x(:,j),ax(:,j))
                call apply_b(x(:,j),bx(:,j))
                r(:,j)=ax(:,j)-eigenvalues(j)*bx(:,j)
                ! With constraints, the eigen-equation includes a term B*Y*mu.
                ! Remove that Lagrange-multiplier component from the residual.
                ! QY is B-orthonormal, hence Y^T*(B*Y)=I.
                if (nconstraints > 0) then
                    r(:,j)=r(:,j)-matmul(bqy(:,1:nconstraints), &
                        matmul(transpose(qy(:,1:nconstraints)),r(:,j)))
                end if
                residuals(j)=norm2(r(:,j))
            end do
        end subroutine compute_residuals

    end subroutine lobpcg

    ! Partial singular-value decomposition by Golub--Kahan bidiagonalization
    ! with double reorthogonalization and a one-sided Jacobi Rayleigh--Ritz
    ! extraction.  No sparse Gram matrix or dense m-by-n copy is constructed.
    ! The current method grows a Krylov basis without implicit restarting:
    ! ncv is a strict upper bound on the basis size.  If the requested triplets
    ! do not converge within this space, return a positive info code.
    ! Smallest singular values require a full right Krylov basis; see docs.
    subroutine svds(a, k, u, s, vt, info, which, tol, ncv, maxiter, v0, &
            iterations, residual_norms)
        type(csr_matrix), intent(in) :: a !! Canonical finite m-by-n CSR matrix.
        integer, intent(in) :: k !! Requested singular triplets, 1 <= k <= min(m,n).
        real(dp), allocatable, intent(out) :: u(:,:) !! Left vectors, shape (m,k).
        real(dp), allocatable, intent(out) :: s(:) !! Singular values, length k, ascending.
        real(dp), allocatable, intent(out) :: vt(:,:) !! Transposed right vectors, (k,n).
        integer, intent(out) :: info !! 0 success; positive steps on nonconvergence; negative on invalid input/breakdown.
        character(len=*), intent(in), optional :: which !! 'LM' (default) or 'SM'.
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance, default 1e-8; 0 means machine precision.
        integer, intent(in), optional :: ncv !! Max Golub--Kahan basis width, default min(q,max(40,8*k+16)).
        integer, intent(in), optional :: maxiter !! Maximum bidiagonalization steps; default ncv.
        real(dp), intent(in), optional :: v0(:) !! Initial vector on shorter side, length q=min(m,n).
        integer, intent(out), optional :: iterations !! Number of Golub--Kahan steps used.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! Combined true left/right singular residuals.
        type(csr_matrix) :: op, opt
        real(dp), allocatable :: vbasis(:,:), ubasis(:,:), work(:), nextv(:), &
            left(:,:), right(:,:), singular(:), residuals(:)
        real(dp) :: alpha, beta, prev_beta, threshold, tolerance, scale
        integer :: m, n, q, p, width, steps, j, i, pass, stat
        character(len=2) :: target
        logical :: is_transposed, converged, valid, seeded

        info = sparse_linalg_invalid_input
        if (present(iterations)) iterations = 0
        allocate(u(0,0), s(0), vt(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr(a)) return
        m = a%n_rows
        n = a%n_cols
        q = min(m,n)
        p = max(m,n)
        if (q < 1 .or. k < 1 .or. k > q) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'SM') return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        tolerance = max(tolerance,100.0_dp*epsilon(1.0_dp))
        width = min(q,max(40,8*k+16))
        if (target == 'SM') width = q
        if (present(ncv)) width = ncv
        if (width < k .or. width > q) return
        if (target == 'SM' .and. width /= q) return
        steps = width
        if (present(maxiter)) steps = min(width,maxiter)
        if (steps < k) return
        if (target == 'SM' .and. steps /= q) return
        if (present(v0)) then
            if (size(v0) /= q) return
            if (.not. all(ieee_is_finite(v0))) return
            if (norm2(v0) == 0.0_dp) return
        end if
        scale = 0.0_dp
        if (size(a%data) > 0) scale = maxval(abs(a%data))
        if (scale == 0.0_dp) then
            deallocate(u,s,vt)
            allocate(u(m,k),s(k),vt(k,n))
            u = 0.0_dp
            s = 0.0_dp
            vt = 0.0_dp
            do j = 1,k
                u(j,j) = 1.0_dp
                vt(j,j) = 1.0_dp
            end do
            if (present(residual_norms)) then
                deallocate(residual_norms)
                allocate(residual_norms(k))
                residual_norms = 0.0_dp
            end if
            info = sparse_linalg_success
            return
        end if
        threshold = 100.0_dp*epsilon(1.0_dp)*scale
        is_transposed = m < n
        if (is_transposed) then
            op = csr_transpose(a)
            opt = a
        else
            op = a
            opt = csr_transpose(a)
        end if
        allocate(vbasis(q,steps),ubasis(p,steps),work(p),nextv(q))
        vbasis = 0.0_dp
        ubasis = 0.0_dp
        if (present(v0)) then
            vbasis(:,1) = v0/norm2(v0)
        else
            do i = 1,q
                vbasis(i,1) = sin(0.6180339887498948_dp*real(i+1,dp)) &
                    + cos(0.3316624790355399_dp*real(i+3,dp))
            end do
            vbasis(:,1) = vbasis(:,1)/norm2(vbasis(:,1))
        end if
        prev_beta = 0.0_dp
        info = max(1,steps)
        do j = 1,steps
            call csr_product(op,vbasis(:,j),work)
            if (j > 1) work = work-prev_beta*ubasis(:,j-1)
            do pass = 1,2
                do i = 1,j-1
                    work = work-dot_product(ubasis(:,i),work)*ubasis(:,i)
                end do
            end do
            alpha = norm2(work)
            beta = 0.0_dp
            if (alpha > threshold) then
                ubasis(:,j) = work/alpha
                call csr_product(opt,ubasis(:,j),nextv)
                nextv = nextv-alpha*vbasis(:,j)
                do pass = 1,2
                    do i = 1,j
                        nextv = nextv-dot_product(vbasis(:,i),nextv)*vbasis(:,i)
                    end do
                end do
                beta = norm2(nextv)
            end if
            if (.not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(beta)) then
                info = sparse_linalg_breakdown
                exit
            end if
            ! Compute true two-sided singular residuals; do not infer
            ! convergence from projected Ritz values alone.
            if (j >= k .and. (j == steps .or. mod(j,8) == 0)) then
                if (target == 'LM' .or. j == q) then
                    call projected_svds(op,opt,vbasis(:,:j),k,target,tolerance, &
                        left,singular,right,residuals,converged,valid)
                    if (.not. valid) then
                        info = sparse_linalg_breakdown
                        exit
                    end if
                    deallocate(u,s,vt)
                    allocate(u(m,k),s(k),vt(k,n))
                    s = singular
                    if (is_transposed) then
                        u = right
                        vt = transpose(left)
                    else
                        u = left
                        vt = transpose(right)
                    end if
                    if (present(residual_norms)) then
                        deallocate(residual_norms)
                        residual_norms = residuals
                    end if
                    if (converged) then
                        info = sparse_linalg_success
                        exit
                    end if
                end if
            end if
            if (j == steps) exit
            if (beta > threshold) then
                vbasis(:,j+1) = nextv/beta
                prev_beta = beta
            else
                ! Lucky breakdown can occur before k distinct singular
                ! directions have been explored (e.g. repeated values).
                ! Start another orthogonal Krylov chain, rather than
                ! returning fewer than the requested k triplets.
                call independent_seed(vbasis(:,:j),nextv,seeded)
                if (.not. seeded) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                vbasis(:,j+1) = nextv
                prev_beta = 0.0_dp
            end if
        end do
        if (present(iterations)) iterations = j
    end subroutine svds

    ! Find a deterministic coordinate direction outside the current right
    ! Krylov space.  This permits new chains after an exact beta breakdown.
    subroutine independent_seed(q, v, ok)
        real(dp), intent(in) :: q(:,:) !! Previously orthonormal directions, n by j.
        real(dp), intent(out) :: v(:) !! New unit direction, length n.
        logical, intent(out) :: ok !! True if an orthogonal complement exists.
        integer :: i, best, pass, j
        real(dp) :: available, largest, length

        best = 0
        largest = 0.0_dp
        do i = 1,size(q,1)
            available = max(0.0_dp,1.0_dp-sum(q(i,:)**2))
            if (available > largest) then
                largest = available
                best = i
            end if
        end do
        ok = .false.
        v = 0.0_dp
        if (best == 0 .or. largest < 1000.0_dp*epsilon(1.0_dp)) return
        v(best) = 1.0_dp
        do pass = 1,2
            do j = 1,size(q,2)
                v = v-dot_product(q(:,j),v)*q(:,j)
            end do
        end do
        length = norm2(v)
        if (length < 100.0_dp*epsilon(1.0_dp)) return
        v = v/length
        ok = .true.
    end subroutine independent_seed

    ! One-sided Jacobi SVD of the small projected operator A*V.  Although
    ! this is a tall p-by-j array, j is at most the user-selected ncv; we
    ! never form A**T*A, which would square the conditioning of A.
    subroutine projected_svds(a,at,basis,k,target,tolerance,left,singular,right,residuals,converged,valid)
        type(csr_matrix), intent(in) :: a !! Tall p-by-q sparse operator.
        type(csr_matrix), intent(in) :: at !! Transposed q-by-p CSR operator.
        real(dp), intent(in) :: basis(:,:) !! Orthonormal q-by-t trial right basis.
        integer, intent(in) :: k !! Number of singular triplets requested.
        character(len=2), intent(in) :: target !! 'LM' or 'SM'.
        real(dp), intent(in) :: tolerance !! Relative true-residual threshold.
        real(dp), allocatable, intent(out) :: left(:,:) !! p-by-k left vectors.
        real(dp), allocatable, intent(out) :: singular(:) !! Selected singular values.
        real(dp), allocatable, intent(out) :: right(:,:) !! q-by-k right vectors.
        real(dp), allocatable, intent(out) :: residuals(:) !! True two-sided residual norms.
        logical, intent(out) :: converged !! Every selected true residual meets tolerance.
        logical, intent(out) :: valid !! Jacobi projected SVD converged and is finite.
        real(dp), allocatable :: w(:,:), rotation(:,:), sl(:), ul(:,:), &
            left_candidate(:), probe(:), v(:), av(:), atu(:)
        integer, allocatable :: order(:)
        real(dp) :: pp, qq, cross, tau, tangent, cs, sn, &
            sigmax, cutoff, bound, rleft, rright, probe_norm, best_norm
        integer :: p, q, t, i, j, h, sweep, pass, choice, rank, best, old
        logical :: changed, small_converged

        valid = .false.
        converged = .false.
        p = a%n_rows
        q = a%n_cols
        t = size(basis,2)
        allocate(w(p,t),rotation(t,t),sl(t),ul(p,t),order(t), &
            left_candidate(p),probe(p),v(q),av(p),atu(q))
        rotation = 0.0_dp
        do j = 1,t
            rotation(j,j) = 1.0_dp
            call csr_product(a,basis(:,j),w(:,j))
        end do
        small_converged = .false.
        do sweep = 1,60
            changed = .false.
            do i = 1,t-1
                do j = i+1,t
                    pp = dot_product(w(:,i),w(:,i))
                    qq = dot_product(w(:,j),w(:,j))
                    cross = dot_product(w(:,i),w(:,j))
                    if (.not. ieee_is_finite(cross) .or. .not. ieee_is_finite(pp) &
                            .or. .not. ieee_is_finite(qq)) return
                    if (pp == 0.0_dp .or. qq == 0.0_dp) cycle
                    if (abs(cross) <= 8.0_dp*epsilon(1.0_dp)*sqrt(pp)*sqrt(qq)) cycle
                    changed = .true.
                    tau = (qq-pp)/(2.0_dp*cross)
                    tangent = sign(1.0_dp,tau)/(abs(tau)+hypot(1.0_dp,tau))
                    cs = 1.0_dp/sqrt(1.0_dp+tangent*tangent)
                    sn = cs*tangent
                    left_candidate = cs*w(:,i)-sn*w(:,j)
                    w(:,j) = sn*w(:,i)+cs*w(:,j)
                    w(:,i) = left_candidate
                    do h = 1,t
                        pp = rotation(h,i)
                        qq = rotation(h,j)
                        rotation(h,i) = cs*pp-sn*qq
                        rotation(h,j) = sn*pp+cs*qq
                    end do
                end do
            end do
            if (.not. changed) then
                small_converged = .true.
                exit
            end if
        end do
        if (.not. small_converged) return
        do j = 1,t
            sl(j) = norm2(w(:,j))
        end do
        if (.not. all(ieee_is_finite(sl))) return
        sigmax = maxval(sl)
        cutoff = 100.0_dp*epsilon(1.0_dp)*sigmax
        ul = 0.0_dp
        rank = 0
        do j = 1,t
            if (sl(j) > cutoff) then
                ul(:,j) = w(:,j)/sl(j)
                rank = rank+1
            else
                sl(j) = 0.0_dp
            end if
        end do
        ! Fill missing zero-singular-value left vectors using the orthogonal
        ! complement of the nonzero left singular subspace.  When t=q this
        ! is an actual left nullspace of A**T, including rank-deficient cases.
        do j = 1,t
            if (sl(j) /= 0.0_dp) cycle
            best = 0
            best_norm = 0.0_dp
            do i = 1,p
                probe = 0.0_dp
                probe(i) = 1.0_dp
                do pass = 1,2
                    do h = 1,t
                        if (h == j) cycle
                        if (any(ul(:,h) /= 0.0_dp)) then
                            probe = probe-dot_product(ul(:,h),probe)*ul(:,h)
                        end if
                    end do
                end do
                probe_norm = norm2(probe)
                if (probe_norm > best_norm) then
                    best_norm = probe_norm
                    best = i
                end if
            end do
            if (best == 0 .or. best_norm < 100.0_dp*epsilon(1.0_dp)) return
            ul(best,j) = 1.0_dp
            do pass = 1,2
                do h = 1,t
                    if (h == j) cycle
                    ul(:,j) = ul(:,j)-dot_product(ul(:,h),ul(:,j))*ul(:,h)
                end do
            end do
            ul(:,j) = ul(:,j)/norm2(ul(:,j))
        end do
        order = [(i,i=1,t)]
        ! Retain largest/smallest targets, but report selected values in
        ! ascending order as SciPy typically does with vectors requested.
        do i = 1,t-1
            choice = i
            do j = i+1,t
                if (target == 'LM') then
                    if (sl(order(j)) > sl(order(choice))) choice = j
                else
                    if (sl(order(j)) < sl(order(choice))) choice = j
                end if
            end do
            old = order(i)
            order(i) = order(choice)
            order(choice) = old
        end do
        allocate(left(p,k),singular(k),right(q,k),residuals(k))
        do j = 1,k
            h = order(j)
            singular(j) = sl(h)
            left(:,j) = ul(:,h)
            right(:,j) = matmul(basis,rotation(:,h))
        end do
        ! Ascending singular values regardless of target selection.
        do i = 1,k-1
            do j = i+1,k
                if (singular(j) < singular(i)) then
                    pp = singular(i)
                    singular(i) = singular(j)
                    singular(j) = pp
                    left_candidate = left(:,i)
                    left(:,i) = left(:,j)
                    left(:,j) = left_candidate
                    v = right(:,i)
                    right(:,i) = right(:,j)
                    right(:,j) = v
                end if
            end do
        end do
        converged = .true.
        bound = tolerance*max(sigmax,tiny(1.0_dp))
        do j = 1,k
            call csr_product(a,right(:,j),av)
            call csr_product(at,left(:,j),atu)
            rleft = norm2(av-singular(j)*left(:,j))
            rright = norm2(atu-singular(j)*right(:,j))
            residuals(j) = hypot(rleft,rright)
            if (.not. ieee_is_finite(residuals(j))) return
            if (residuals(j) > bound) converged = .false.
        end do
        valid = .true.
    end subroutine projected_svds

    subroutine lsmr(a, b, x, istop, itn, normr, damp, atol, btol, conlim, maxiter)
        type(csr_matrix), intent(in) :: a !! Real m-by-n canonical CSR matrix.
        real(dp), intent(in) :: b(:) !! RHS vector of length m.
        real(dp), allocatable, intent(out) :: x(:) !! Least-squares vector of length n.
        integer, intent(out) :: istop !! SciPy LSMR stop code; -1 invalid input.
        integer, intent(out), optional :: itn !! Number of Golub-Kahan steps.
        real(dp), intent(out), optional :: normr !! Estimated residual norm.
        real(dp), intent(in), optional :: damp !! Tikhonov damping, default zero.
        real(dp), intent(in), optional :: atol !! Matrix backward-error tolerance, default 1e-6.
        real(dp), intent(in), optional :: btol !! RHS backward-error tolerance, default 1e-6.
        real(dp), intent(in), optional :: conlim !! Condition-number stopping limit, default 1e8.
        integer, intent(in), optional :: maxiter !! Iteration cap, default min(m,n).
        type(csr_matrix) :: at
        real(dp) :: damping, abs_tol, rel_tol, limit, norm_out
        integer :: max_steps, steps, status

        allocate(x(0))
        istop = sparse_linalg_invalid_input
        if (present(itn)) itn = 0
        if (present(normr)) normr = huge(1.0_dp)
        if (.not. valid_csr(a)) return
        if (size(b) /= a%n_rows .or. .not. all(ieee_is_finite(b))) return
        damping = 0.0_dp
        abs_tol = 1.0e-6_dp
        rel_tol = 1.0e-6_dp
        limit = 1.0e8_dp
        max_steps = min(a%n_rows,a%n_cols)
        if (present(damp)) damping = damp
        if (present(atol)) abs_tol = atol
        if (present(btol)) rel_tol = btol
        if (present(conlim)) limit = conlim
        if (present(maxiter)) max_steps = maxiter
        if (.not. all(ieee_is_finite([damping,abs_tol,rel_tol,limit]))) return
        if (minval([damping,abs_tol,rel_tol,limit]) < 0.0_dp .or. max_steps < 0) return
        deallocate(x)
        allocate(x(a%n_cols))
        x = 0.0_dp
        if (a%n_cols == 0 .or. a%n_rows == 0) then
            istop = 0
            if (present(normr)) normr = norm2(b)
            return
        end if
        at = csr_transpose(a,status)
        if (status /= 0) then
            istop = sparse_linalg_invalid_input
            return
        end if
        call lsmr_solve(matvec, rmatvec, a%n_rows, a%n_cols, b, x, &
                damping, abs_tol, rel_tol, limit, max_steps, istop, steps, normr=norm_out)
        if (present(itn)) itn = steps
        if (present(normr)) normr = normr_true()
    contains
        subroutine matvec(v, product)
            real(dp), intent(in) :: v(:) !! Vector for A multiplication.
            real(dp), intent(out) :: product(:) !! Resulting product A*v.
            call csr_product(a,v,product)
        end subroutine matvec

        subroutine rmatvec(v, product)
            real(dp), intent(in) :: v(:) !! Vector for transpose(A) multiplication.
            real(dp), intent(out) :: product(:) !! Resulting A^T*v.
            call csr_product(at,v,product)
        end subroutine rmatvec

        real(dp) function normr_true() result(rnorm)
            real(dp) :: product(a%n_rows)
            call csr_product(a,x,product)
            rnorm = norm2(b-product)
        end function normr_true
    end subroutine lsmr
    ! Zero-fill incomplete LU factorization.  All entries are retained on the
    ! canonical CSR sparsity pattern, with strict lower/upper factors stored
    ! together.  This is ILU(0), not a wrapper around SuperLU's spilu (ILUTP).
    pure logical function valid_ilu0(factor) result(ok)
        type(ilu0_factor), intent(in) :: factor !! Factorization to inspect.
        integer :: i, p

        ok = .false.
        if (.not. valid_csr(factor%lu)) return
        if (factor%lu%n_rows /= factor%lu%n_cols) return
        if (.not. allocated(factor%diag_ptr)) return
        if (size(factor%diag_ptr) /= factor%lu%n_rows) return
        do i = 1, factor%lu%n_rows
            p = factor%diag_ptr(i)
            if (p < factor%lu%indptr(i) .or. p >= factor%lu%indptr(i+1)) return
            if (factor%lu%indices(p) /= i) return
            if (abs(factor%lu%data(p)) <= tiny(1.0_dp)) return
        end do
        ok = .true.
    end function valid_ilu0

    ! Find a column's stored index in an already-sorted CSR row.
    pure integer function row_position(a, row, column) result(position)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: row !! One-based row number.
        integer, intent(in) :: column !! One-based column number.
        integer :: lo, hi, mid

        position = 0
        lo = a%indptr(row)
        hi = a%indptr(row+1)-1
        do while (lo <= hi)
            mid = lo+(hi-lo)/2
            if (a%indices(mid) == column) then
                position = mid
                return
            else if (a%indices(mid) < column) then
                lo = mid+1
            else
                hi = mid-1
            end if
        end do
    end function row_position

    subroutine ilu0(a, factor, info, diagonal_shift)
        type(csr_matrix), intent(in) :: a !! Square canonical CSR input, finite entries.
        type(ilu0_factor), intent(out) :: factor !! L (unit diagonal) and U in CSR pattern.
        integer, intent(out) :: info !! Zero success; -1 invalid, -3 zero or unstable pivot.
        real(dp), intent(in), optional :: diagonal_shift !! Nonnegative shift added to the diagonal.
        integer :: n, i, p, j, q, target
        real(dp) :: shift, multiplier, pivot, updated

        info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (n /= a%n_cols) return
        shift = 0.0_dp
        if (present(diagonal_shift)) shift = diagonal_shift
        if (.not. ieee_is_finite(shift) .or. shift < 0.0_dp) return
        allocate(factor%diag_ptr(n))
        do i = 1, n
            factor%diag_ptr(i) = row_position(a, i, i)
            if (factor%diag_ptr(i) == 0) then
                info = sparse_linalg_bad_diagonal
                return
            end if
        end do
        factor%lu = a
        do i = 1, n
            if (shift > 0.0_dp) then
                p = factor%diag_ptr(i)
                factor%lu%data(p) = factor%lu%data(p)+shift
                if (.not. ieee_is_finite(factor%lu%data(p))) then
                    info = sparse_linalg_breakdown
                    return
                end if
            end if
        end do
        do i = 1, n
            do p = factor%lu%indptr(i), factor%diag_ptr(i)-1
                j = factor%lu%indices(p)
                pivot = factor%lu%data(factor%diag_ptr(j))
                if (abs(pivot) <= tiny(1.0_dp)) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
                multiplier = factor%lu%data(p)/pivot
                if (.not. ieee_is_finite(multiplier)) then
                    info = sparse_linalg_breakdown
                    return
                end if
                factor%lu%data(p) = multiplier
                do q = factor%diag_ptr(j)+1, factor%lu%indptr(j+1)-1
                    target = row_position(factor%lu,i,factor%lu%indices(q))
                    if (target == 0) cycle  ! No fill: ILU(0).
                    updated = factor%lu%data(target)-multiplier*factor%lu%data(q)
                    if (.not. ieee_is_finite(updated)) then
                        info = sparse_linalg_breakdown
                        return
                    end if
                    factor%lu%data(target) = updated
                end do
            end do
            if (abs(factor%lu%data(factor%diag_ptr(i))) <= tiny(1.0_dp)) then
                info = sparse_linalg_bad_diagonal
                return
            end if
        end do
        info = sparse_linalg_success
    end subroutine ilu0

    ! Internal solve for a previously validated ILU(0) factor; avoids an
    ! allocation and a full validation scan on each Krylov iteration.
    pure subroutine ilu0_apply(factor, b, x)
        type(ilu0_factor), intent(in) :: factor !! Previously validated factors.
        real(dp), intent(in) :: b(:) !! Right-hand side, length n.
        real(dp), intent(out) :: x(:) !! M^{-1} b, length n.
        integer :: i, k, n

        n = size(b)
        x = b
        do i = 1, n
            do k = factor%lu%indptr(i), factor%diag_ptr(i)-1
                x(i) = x(i)-factor%lu%data(k)*x(factor%lu%indices(k))
            end do
        end do
        do i = n, 1, -1
            do k = factor%diag_ptr(i)+1, factor%lu%indptr(i+1)-1
                x(i) = x(i)-factor%lu%data(k)*x(factor%lu%indices(k))
            end do
            x(i) = x(i)/factor%lu%data(factor%diag_ptr(i))
        end do
    end subroutine ilu0_apply

    subroutine ilu0_solve(factor, b, x, info)
        type(ilu0_factor), intent(in) :: factor !! Previously factorized square matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side, length n.
        real(dp), allocatable, intent(out) :: x(:) !! Preconditioned result.
        integer, intent(out) :: info !! Zero success; negative validation/breakdown.

        allocate(x(0))
        info = sparse_linalg_invalid_input
        if (.not. valid_ilu0(factor)) return
        if (size(b) /= factor%lu%n_rows) return
        if (.not. all(ieee_is_finite(b))) return
        deallocate(x)
        allocate(x(size(b)))
        call ilu0_apply(factor,b,x)
        if (.not. all(ieee_is_finite(x))) then
            info = sparse_linalg_breakdown
        else
            info = sparse_linalg_success
        end if
    end subroutine ilu0_solve

    pure subroutine apply_square_preconditioner(invdiag, r, z, ilu)
        real(dp), intent(in) :: invdiag(:) !! Jacobi inverse (all ones if disabled).
        real(dp), intent(in) :: r(:) !! Vector to precondition.
        real(dp), intent(out) :: z(:) !! Preconditioned vector.
        type(ilu0_factor), intent(in), optional :: ilu !! Optional ILU(0) factor.

        if (present(ilu)) then
            call ilu0_apply(ilu,r,z)
        else
            z = invdiag*r
        end if
    end subroutine apply_square_preconditioner

    ! Incomplete LU with threshold dropping and bounded fill (ILUT).
    !
    ! This is a dependency-free, natural-order, no-pivot ILUT implementation.
    ! Unlike scipy.sparse.linalg.spilu (SuperLU's ILUTP), it does not reorder,
    ! equilibrate, or pivot.  A caller may add a nonnegative diagonal shift
    ! to avoid a small or missing pivot.  The factor type is deliberately the
    ! same as ILU(0), allowing the current GMRES/BiCGSTAB preconditioner API
    ! and triangular application to be reused without conversion.
    !
    ! Threshold dropping compares off-diagonal entries to drop_tol times the
    ! largest absolute entry of the ORIGINAL row.  Each row retains at most
    ! int(fill_factor*max(1,original_count)) strict-lower and strict-upper
    ! entries separately (capped by n), chosen by magnitude.  The diagonal
    ! is never dropped.  These are per-row budgets, not SuperLU's global
    ! storage accounting, so the sparsity patterns need not match SciPy.
    subroutine spilu_csr(a, factor, info, drop_tol, fill_factor, diagonal_shift)
        type(csr_matrix), intent(in) :: a !! Canonical square one-based real CSR matrix.
        type(ilu0_factor), intent(out) :: factor !! Combined ILUT factors compatible with ilu0_solve.
        integer, intent(out) :: info !! 0 success; -1 invalid; -3 failed pivot; -2 numerical failure.
        real(dp), intent(in), optional :: drop_tol !! Relative dropping threshold [0,1], default 1e-4.
        real(dp), intent(in), optional :: fill_factor !! Per-row fill multiplier >=1, default 10.
        real(dp), intent(in), optional :: diagonal_shift !! Nonnegative diagonal shift, default 0.
        type(sparse_lu_row), allocatable :: rows(:)
        integer :: n, i, j, k, p, nnz, ofs, nlo, nhi, limit_lo, limit_hi, used
        real(dp) :: tol, fill, shift, row_scale, pivot, multiplier, value, updated

        info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (n /= a%n_cols) return
        tol = 1.0e-4_dp
        if (present(drop_tol)) tol = drop_tol
        fill = 10.0_dp
        if (present(fill_factor)) fill = fill_factor
        shift = 0.0_dp
        if (present(diagonal_shift)) shift = diagonal_shift
        if (.not. ieee_is_finite(tol) .or. .not. ieee_is_finite(fill) .or. &
            .not. ieee_is_finite(shift)) return
        if (tol < 0.0_dp .or. tol > 1.0_dp .or. fill < 1.0_dp .or. shift < 0.0_dp) return

        allocate(rows(n),factor%diag_ptr(n))
        do i = 1,n
            nnz = a%indptr(i+1)-a%indptr(i)
            allocate(rows(i)%columns(nnz),rows(i)%values(nnz))
            rows(i)%used = nnz
            rows(i)%columns = a%indices(a%indptr(i):a%indptr(i+1)-1)
            rows(i)%values = a%data(a%indptr(i):a%indptr(i+1)-1)
            row_scale = 0.0_dp
            if (nnz > 0) row_scale = maxval(abs(rows(i)%values))
            if (shift > 0.0_dp) then
                value = sparse_row_value(rows(i),i)+shift
                if (.not. ieee_is_finite(value)) then
                    info = sparse_linalg_breakdown
                    return
                end if
                call sparse_row_set(rows(i),i,value)
            end if
            nlo = count(rows(i)%columns(:rows(i)%used) < i)
            nhi = count(rows(i)%columns(:rows(i)%used) > i)
            limit_lo = spilu_budget(fill,nlo,n)
            limit_hi = spilu_budget(fill,nhi,n)

            p = 1
            do while (p <= rows(i)%used)
                j = rows(i)%columns(p)
                if (j >= i) exit
                pivot = sparse_row_value(rows(j),j)
                if (abs(pivot) <= tiny(1.0_dp)) then
                    info = sparse_linalg_bad_diagonal
                    return
                end if
                multiplier = rows(i)%values(p)/pivot
                if (.not. ieee_is_finite(multiplier)) then
                    info = sparse_linalg_breakdown
                    return
                end if
                if (abs(multiplier) <= tol*row_scale) then
                    call sparse_row_set(rows(i),j,0.0_dp)
                    cycle
                end if
                rows(i)%values(p) = multiplier
                p = p+1
                do k = 1,rows(j)%used
                    if (rows(j)%columns(k) <= j) cycle
                    value = sparse_row_value(rows(i),rows(j)%columns(k))
                    updated = value-multiplier*rows(j)%values(k)
                    if (.not. ieee_is_finite(updated)) then
                        info = sparse_linalg_breakdown
                        return
                    end if
                    call sparse_row_set(rows(i),rows(j)%columns(k),updated)
                end do
            end do
            pivot = sparse_row_value(rows(i),i)
            if (abs(pivot) <= tiny(1.0_dp)) then
                info = sparse_linalg_bad_diagonal
                return
            end if
            call spilu_prune_row(rows(i),i,tol*row_scale,limit_lo,limit_hi)
        end do

        nnz = 0
        do i = 1,n
            if (rows(i)%used > huge(nnz)-nnz) then
                info = sparse_linalg_breakdown
                return
            end if
            nnz = nnz+rows(i)%used
        end do
        factor%lu%n_rows = n
        factor%lu%n_cols = n
        allocate(factor%lu%indptr(n+1),factor%lu%indices(nnz),factor%lu%data(nnz))
        ofs = 1
        do i = 1,n
            factor%lu%indptr(i) = ofs
            used = rows(i)%used
            factor%lu%indices(ofs:ofs+used-1) = rows(i)%columns(:used)
            factor%lu%data(ofs:ofs+used-1) = rows(i)%values(:used)
            factor%diag_ptr(i) = 0
            do j = 1,used
                if (rows(i)%columns(j) == i) factor%diag_ptr(i) = ofs+j-1
            end do
            ofs = ofs+used
        end do
        factor%lu%indptr(n+1) = ofs
        if (.not. valid_ilu0(factor)) then
            info = sparse_linalg_breakdown
        else
            info = sparse_linalg_success
        end if
    end subroutine spilu_csr

    ! Conversion adapter: preserves the one-based canonical CSC interface.
    subroutine spilu_csc(a, factor, info, drop_tol, fill_factor, diagonal_shift)
        type(csc_matrix), intent(in) :: a !! Square canonical one-based CSC matrix.
        type(ilu0_factor), intent(out) :: factor !! Resulting ILUT factors.
        integer, intent(out) :: info !! Factorization status.
        real(dp), intent(in), optional :: drop_tol !! Relative drop tolerance.
        real(dp), intent(in), optional :: fill_factor !! Per-row fill multiplier.
        real(dp), intent(in), optional :: diagonal_shift !! Nonnegative diagonal perturbation.
        type(csr_matrix) :: a_csr

        a_csr = csc_to_csr(a)
        call spilu_csr(a_csr,factor,info,drop_tol,fill_factor,diagonal_shift)
    end subroutine spilu_csc

    pure integer function spilu_budget(fill, originals, n) result(budget)
        real(dp), intent(in) :: fill !! Finite fill multiplier >= 1.
        integer, intent(in) :: originals !! Number of off-diagonal entries in the original row.
        integer, intent(in) :: n !! Matrix dimension.

        if (fill >= real(n,dp)) then
            budget = n
        else
            budget = int(min(real(n,dp), fill*real(max(1,originals),dp)))
        end if
    end function spilu_budget

    subroutine spilu_prune_row(row, diagonal, threshold, max_lower, max_upper)
        type(sparse_lu_row), intent(inout) :: row !! Completed ILU row, sorted by column.
        integer, intent(in) :: diagonal !! Diagonal column, never dropped.
        real(dp), intent(in) :: threshold !! Absolute drop threshold for off-diagonal entries.
        integer, intent(in) :: max_lower !! Maximum number of strict-lower entries.
        integer, intent(in) :: max_upper !! Maximum number of strict-upper entries.
        logical, allocatable :: keep(:)
        integer :: k, j, count_keep, limit, smallest, selected, used
        real(dp) :: magnitude

        used = row%used
        allocate(keep(used))
        keep = .false.
        do k = 1,used
            keep(k) = row%columns(k) == diagonal .or. abs(row%values(k)) > threshold
        end do
        do j = 1,2
            limit = max_lower
            if (j == 2) limit = max_upper
            count_keep = 0
            do k = 1,used
                if (.not. keep(k)) cycle
                if ((j == 1 .and. row%columns(k) < diagonal) .or. &
                    (j == 2 .and. row%columns(k) > diagonal)) count_keep = count_keep+1
            end do
            do while (count_keep > limit)
                smallest = 0
                magnitude = huge(1.0_dp)
                do k = 1,used
                    if (.not. keep(k)) cycle
                    if (j == 1 .and. row%columns(k) >= diagonal) cycle
                    if (j == 2 .and. row%columns(k) <= diagonal) cycle
                    if (abs(row%values(k)) < magnitude) then
                        magnitude = abs(row%values(k))
                        smallest = k
                    end if
                end do
                if (smallest == 0) exit
                keep(smallest) = .false.
                count_keep = count_keep-1
            end do
        end do
        selected = 0
        do k = 1,used
            if (.not. keep(k)) cycle
            selected = selected+1
            row%columns(selected) = row%columns(k)
            row%values(selected) = row%values(k)
        end do
        row%used = selected
    end subroutine spilu_prune_row

    subroutine spilu_solve_vector(factor, b, x, info)
        type(ilu0_factor), intent(in) :: factor !! Valid ILUT or ILU(0) factors.
        real(dp), intent(in) :: b(:) !! Finite RHS vector of length n.
        real(dp), allocatable, intent(out) :: x(:) !! Preconditioned vector.
        integer, intent(out) :: info !! 0 success, -1 invalid, -2 numerical failure.

        call ilu0_solve(factor,b,x,info)
    end subroutine spilu_solve_vector

    subroutine spilu_solve_matrix(factor, b, x, info)
        type(ilu0_factor), intent(in) :: factor !! Valid ILUT or ILU(0) factors.
        real(dp), intent(in) :: b(:,:) !! Finite RHS array, shape (n,nrhs).
        real(dp), allocatable, intent(out) :: x(:,:) !! Preconditioned block.
        integer, intent(out) :: info !! 0 success, -1 invalid, -2 numerical failure.
        integer :: n, j

        info = sparse_linalg_invalid_input
        allocate(x(0,0))
        if (.not. valid_ilu0(factor)) return
        n = factor%lu%n_rows
        if (size(b,1) /= n) return
        if (.not. all(ieee_is_finite(b))) return
        deallocate(x)
        allocate(x(n,size(b,2)))
        do j = 1,size(b,2)
            call ilu0_apply(factor,b(:,j),x(:,j))
        end do
        info = sparse_linalg_success
        if (.not. all(ieee_is_finite(x))) info = sparse_linalg_breakdown
    end subroutine spilu_solve_matrix

    ! Golub--Kahan bidiagonalization with Paige--Saunders LSQR recurrences.
    ! Method implemented directly, without copying external LSQR source code.
    ! Optional x0 is the center of the Tikhonov penalty, as in SciPy.
    subroutine lsqr(a, b, x, istop, itn, r1norm, r2norm, damp, atol, btol, &
            conlim, iter_lim, x0, anorm, acond, arnorm, xnorm, var, calc_var)
        type(csr_matrix), intent(in) :: a !! Canonical rectangular real CSR operator.
        real(dp), intent(in) :: b(:) !! Right-hand side of length m.
        real(dp), allocatable, intent(out) :: x(:) !! Approximate least-squares solution, length n.
        integer, intent(out) :: istop !! SciPy LSQR codes 0..3, 7; negative on invalid/breakdown.
        integer, intent(out), optional :: itn !! Number of bidiagonalization steps.
        real(dp), intent(out), optional :: r1norm !! True norm of b-A*x.
        real(dp), intent(out), optional :: r2norm !! Norm of augmented residual.
        real(dp), intent(in), optional :: damp !! Nonnegative Tikhonov penalty, default 0.
        real(dp), intent(in), optional :: atol !! Matrix-error stopping tolerance, default 1e-6.
        real(dp), intent(in), optional :: btol !! RHS-error stopping tolerance, default 1e-6.
        real(dp), intent(in), optional :: conlim !! Condition-number bound, default 1e8 (0 disables).
        integer, intent(in), optional :: iter_lim !! Iteration cap, default 2*min(m,n).
        real(dp), intent(in), optional :: x0(:) !! Initial estimate and damping center, default zero.
        real(dp), intent(out), optional :: anorm !! Estimated Frobenius norm of [A; damp*I].
        real(dp), intent(out), optional :: acond !! Estimated condition number of augmented operator.
        real(dp), intent(out), optional :: arnorm !! True normal-equation residual norm.
        real(dp), intent(out), optional :: xnorm !! Norm of final solution.
        real(dp), allocatable, intent(out), optional :: var(:) !! Optional variance-diagonal estimate.
        logical, intent(in), optional :: calc_var !! Calculate variance-diagonal estimate.
        type(csr_matrix) :: at
        real(dp), allocatable :: base(:), dx(:), u(:), v(:), w(:), work_m(:), work_n(:), gradient(:)
        real(dp) :: damping, tol_a, tol_b, cond_limit, alpha, beta, beta0
        real(dp) :: rhobar, phibar, rhobar1, cs1, sn1, psi, rho, cs, sn
        real(dp) :: theta, phi, tau, ddnorm, residual2, est_rnorm, est_arnorm
        real(dp) :: norm_a, condition, norm_x, test_bound
        integer :: m, n, cap, steps, status
        logical :: estimate_var

        istop = sparse_linalg_invalid_input
        if (present(itn)) itn = 0
        if (present(r1norm)) r1norm = huge(1.0_dp)
        if (present(r2norm)) r2norm = huge(1.0_dp)
        if (present(anorm)) anorm = 0.0_dp
        if (present(acond)) acond = 0.0_dp
        if (present(arnorm)) arnorm = huge(1.0_dp)
        if (present(xnorm)) xnorm = 0.0_dp
        if (present(var)) allocate(var(0))
        allocate(x(0))
        if (.not. valid_csr(a)) return
        m = a%n_rows
        n = a%n_cols
        if (size(b) /= m .or. .not. all(ieee_is_finite(b))) return
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
        base = 0.0_dp
        if (present(x0)) then
            if (size(x0) /= n) return
            if (.not. all(ieee_is_finite(x0))) return
            base = x0
        end if
        deallocate(x)
        allocate(x(n),dx(n),u(m),v(n),w(n),work_m(m),work_n(n),gradient(n))
        x = base
        dx = 0.0_dp
        estimate_var = .false.
        if (present(calc_var)) estimate_var = calc_var
        if (present(var)) then
            deallocate(var)
            allocate(var(n))
            var = 0.0_dp
        end if
        at = csr_transpose(a,status)
        if (status /= 0) return
        call csr_product(a,base,work_m)
        u = b-work_m
        beta = norm2(u)
        beta0 = beta
        norm_a = 0.0_dp
        condition = 0.0_dp
        norm_x = 0.0_dp
        steps = 0
        istop = 0
        if (beta == 0.0_dp .or. n == 0 .or. m == 0) goto 900
        u = u/beta
        call csr_product(at,u,work_n)
        alpha = norm2(work_n)
        if (alpha == 0.0_dp) goto 900
        v = work_n/alpha
        w = v
        rhobar = alpha
        phibar = beta
        ddnorm = 0.0_dp
        residual2 = 0.0_dp
        istop = 7
        do steps = 1, cap
            call csr_product(a,v,work_m)
            u = work_m-alpha*u
            beta = norm2(u)
            if (beta /= 0.0_dp) u = u/beta
            if (beta /= 0.0_dp) then
                call csr_product(at,u,work_n)
            else
                work_n = 0.0_dp
            end if
            v = work_n-beta*v
            alpha = norm2(v)
            if (alpha /= 0.0_dp) v = v/alpha
            norm_a = hypot(norm_a,hypot(alpha,hypot(beta,damping)))

            rhobar1 = hypot(rhobar,damping)
            if (rhobar1 == 0.0_dp) then
                istop = sparse_linalg_breakdown
                exit
            end if
            cs1 = rhobar/rhobar1
            sn1 = damping/rhobar1
            psi = sn1*phibar
            phibar = cs1*phibar
            rho = hypot(rhobar1,beta)
            if (rho == 0.0_dp) then
                istop = sparse_linalg_breakdown
                exit
            end if
            cs = rhobar1/rho
            sn = beta/rho
            theta = sn*alpha
            rhobar = -cs*alpha
            phi = cs*phibar
            phibar = sn*phibar
            tau = sn*phi
            dx = dx+(phi/rho)*w
            if (present(var)) then
                if (estimate_var) var = var+(w/rho)**2
            end if
            ddnorm = ddnorm+dot_product(w/rho,w/rho)
            w = v-(theta/rho)*w
            residual2 = residual2+psi**2
            est_rnorm = hypot(phibar,sqrt(residual2))
            est_arnorm = abs(alpha*tau)
            norm_x = norm2(dx)
            condition = norm_a*sqrt(ddnorm)
            if (.not. all(ieee_is_finite([norm_x,norm_a,condition,est_rnorm,est_arnorm]))) then
                istop = sparse_linalg_breakdown
                exit
            end if
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
900     continue
        call csr_product(a,x,work_m)
        work_m = b-work_m
        if (present(r1norm)) r1norm = norm2(work_m)
        if (present(r2norm)) r2norm = hypot(norm2(work_m),damping*norm2(x-base))
        if (present(itn)) itn = steps
        if (present(anorm)) anorm = norm_a
        if (present(acond)) acond = condition
        if (present(xnorm)) xnorm = norm2(x)
        if (present(arnorm)) then
            call csr_product(at,work_m,gradient)
            gradient = gradient-damping**2*(x-base)
            arnorm = norm2(gradient)
        end if
    end subroutine lsqr

    ! Act with exp(t*A) on one or more vectors without forming exp(t*A).
    ! For each scaled step, use the Taylor series on the trace-shifted
    ! operator; an induced infinity-norm bound controls the omitted tail.
    ! This is an independently written, conservative variant of the
    ! scaling/Taylor action approach of Al-Mohy and Higham (2011).
    pure subroutine csr_block_product(a, x, y)
        type(csr_matrix), intent(in) :: a !! Prevalidated n-by-n CSR matrix.
        real(dp), intent(in) :: x(:,:) !! Dense block of RHS vectors, shape (n,nrhs).
        real(dp), intent(out) :: y(:,:) !! Matrix-block product of shape (n,nrhs).
        integer :: i, k

        y = 0.0_dp
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                y(i,:) = y(i,:) + a%data(k)*x(a%indices(k),:)
            end do
        end do
    end subroutine csr_block_product

    pure subroutine csr_shift_infinity_norm(a, mu, operator_norm, info)
        type(csr_matrix), intent(in) :: a !! Prevalidated real CSR input.
        real(dp), intent(in) :: mu !! Scalar diagonal shift.
        real(dp), intent(out) :: operator_norm !! Infinity norm of A - mu*I.
        integer, intent(out) :: info !! Zero if all shifted row sums are finite.
        integer :: i, k
        real(dp) :: diagonal, row_sum

        info = sparse_linalg_success
        operator_norm = 0.0_dp
        do i = 1, a%n_rows
            row_sum = 0.0_dp
            diagonal = -mu
            do k = a%indptr(i), a%indptr(i+1)-1
                if (a%indices(k) == i) then
                    diagonal = a%data(k) - mu
                else
                    row_sum = row_sum + abs(a%data(k))
                end if
            end do
            row_sum = row_sum + abs(diagonal)
            if (.not. ieee_is_finite(row_sum)) then
                info = sparse_linalg_breakdown
                return
            end if
            operator_norm = max(operator_norm,row_sum)
        end do
    end subroutine csr_shift_infinity_norm

    !> Approximate f(t*A)*b using short Arnoldi/Lanczos projections and the
    !! restarted matrix-function Krylov construction of Eiermann--Ernst (2006).
    !! This overload accepts a real, canonical, one-based CSR matrix.  It
    !! delegates to the matrix-free version, with the same options and status.
    !! A full cycle may be approximate; info>0 signals an unverified truncation.
    subroutine funm_multiply_krylov_csr(a, b, f, y, info, t, rtol, atol, &
                                       restart_every_m, max_restarts, assume_hermitian, &
                                       cycles, products, update_norm)
        type(csr_matrix), intent(in) :: a !! Finite square, canonical CSR matrix.
        real(dp), intent(in) :: b(:) !! Vector of length a%n_rows.
        procedure(dense_matrix_function) :: f !! Callable f for a small dense square matrix.
        real(dp), allocatable, intent(out) :: y(:) !! Approximation to f(t*A)*b.
        integer, intent(out) :: info !! 0 success, positive unconverged cycles, negative error.
        real(dp), intent(in), optional :: t !! Matrix scaling t, default 1.
        real(dp), intent(in), optional :: rtol !! Relative update tolerance, default 1e-6.
        real(dp), intent(in), optional :: atol !! Absolute update tolerance, default 0.
        integer, intent(in), optional :: restart_every_m !! Krylov vectors per cycle, default min(20,n).
        integer, intent(in), optional :: max_restarts !! Maximum Krylov cycles, default 20.
        logical, intent(in), optional :: assume_hermitian !! Use Lanczos for symmetric A, default false.
        integer, intent(out), optional :: cycles !! Number of performed Krylov cycles.
        integer, intent(out), optional :: products !! Number of sparse matvec operations.
        real(dp), intent(out), optional :: update_norm !! Norm of final Krylov update.

        if (.not. valid_csr(a)) then
            allocate(y(0))
            info = sparse_linalg_invalid_input
            if (present(cycles)) cycles = 0
            if (present(products)) products = 0
            if (present(update_norm)) update_norm = 0.0_dp
            return
        end if
        if (a%n_rows /= a%n_cols) then
            allocate(y(0))
            info = sparse_linalg_invalid_input
            if (present(cycles)) cycles = 0
            if (present(products)) products = 0
            if (present(update_norm)) update_norm = 0.0_dp
            return
        end if
        call funm_multiply_krylov_operator(a%n_rows, b, csr_apply, f, y, info, &
            t, rtol, atol, restart_every_m, max_restarts, assume_hermitian, &
            cycles, products, update_norm)
    contains
        subroutine csr_apply(x, result)
            real(dp), intent(in) :: x(:) !! Candidate Krylov vector.
            real(dp), intent(out) :: result(:) !! Sparse matrix-vector result.
            call csr_product(a, x, result)
        end subroutine csr_apply
    end subroutine funm_multiply_krylov_csr

    !> Generalized restarted matrix-function Krylov action with user matvec.
    !! Full (modified-Gram-Schmidt/reorthogonalized) Arnoldi handles general
    !! matrices; three-term Lanczos handles user-declared symmetric matrices.
    !! A callback f evaluates the function only on projected dense matrices.
    !! Following SciPy's convention, a small correction between successive
    !! cycles ends the iteration, but that is an *update test*, not a rigorous
    !! a-posteriori error estimate for arbitrary functions.
    !! Positive info means the cycle budget was exhausted without satisfying
    !! that test.  Negative info denotes invalid data or failed callbacks.
    subroutine funm_multiply_krylov_operator(n, b, matvec, f, y, info, t, rtol, atol, &
                                            restart_every_m, max_restarts, assume_hermitian, &
                                            cycles, products, update_norm)
        integer, intent(in) :: n !! Square operator order, >=0.
        real(dp), intent(in) :: b(:) !! Finite right-hand vector of length n.
        procedure(real_matvec_operator) :: matvec !! Real matrix-free x -> A*x.
        procedure(dense_matrix_function) :: f !! Dense f(H) evaluator, status=0 when valid.
        real(dp), allocatable, intent(out) :: y(:) !! Krylov approximation, shape (n).
        integer, intent(out) :: info !! 0 converged; >0 budget exceeded; <0 invalid/breakdown.
        real(dp), intent(in), optional :: t !! Real time/argument scale, default 1.
        real(dp), intent(in), optional :: rtol !! Relative update tolerance, default 1e-6.
        real(dp), intent(in), optional :: atol !! Absolute update tolerance, default 0.
        integer, intent(in), optional :: restart_every_m !! Local basis order, default min(20,n).
        integer, intent(in), optional :: max_restarts !! Maximum number of cycles, default 20.
        logical, intent(in), optional :: assume_hermitian !! Use symmetric Lanczos, default false.
        integer, intent(out), optional :: cycles !! Number of finished cycles, 0 on input failure.
        integer, intent(out), optional :: products !! Total calls to matvec.
        real(dp), intent(out), optional :: update_norm !! Norm of most recent correction.

        integer :: m, cap, capacity, next_capacity, ic, base, j, global_j, &
                   nruns, nproducts, function_status
        real(dp) :: time_scale, reltol, abstol, bnorm, threshold, change
        real(dp), allocatable :: v(:,:), h(:,:), fh(:,:), seed(:), candidate(:), grown_h(:,:)
        logical :: symmetric, happy

        allocate(y(0))
        info = sparse_linalg_invalid_input
        nruns = 0
        nproducts = 0
        change = 0.0_dp
        if (present(cycles)) cycles = 0
        if (present(products)) products = 0
        if (present(update_norm)) update_norm = 0.0_dp
        if (n < 0 .or. n /= size(b)) return
        if (.not. all(ieee_is_finite(b))) return
        time_scale = 1.0_dp
        reltol = 1.0e-6_dp
        abstol = 0.0_dp
        cap = 20
        symmetric = .false.
        if (present(t)) time_scale = t
        if (present(rtol)) reltol = rtol
        if (present(atol)) abstol = atol
        if (present(max_restarts)) cap = max_restarts
        if (present(assume_hermitian)) symmetric = assume_hermitian
        if (.not. ieee_is_finite(time_scale) .or. .not. ieee_is_finite(reltol) .or. &
            .not. ieee_is_finite(abstol)) return
        if (reltol < 0.0_dp .or. abstol < 0.0_dp .or. cap <= 0) return
        if (present(restart_every_m)) then
            if (restart_every_m <= 0) return
        end if
        if (n == 0) then
            deallocate(y)
            allocate(y(0))
            info = sparse_linalg_success
            return
        end if
        m = min(20, n)
        if (present(restart_every_m)) m = min(restart_every_m,n)
        ! A complete n-dimensional basis needs no additional restart.  For
        ! shorter bases, respect the user's cycle budget: later corrections
        ! can be useful even beyond n/m+1 cycles (notably for sin and exp).
        if (m >= n) cap = 1
        bnorm = norm2(b)
        if (.not. ieee_is_finite(bnorm)) then
            info = sparse_linalg_breakdown
            return
        end if
        deallocate(y)
        allocate(y(n))
        y = 0.0_dp
        if (bnorm == 0.0_dp) then
            info = sparse_linalg_success
            return
        end if
        threshold = max(abstol, reltol*bnorm)
        if (.not. ieee_is_finite(threshold)) then
            info = sparse_linalg_invalid_input
            return
        end if
        ! A large requested restart budget should not allocate all future
        ! dense projected matrices in advance.  Grow only if more cycles are
        ! actually needed; most calls use just one or two cycles.
        capacity = min(cap,2)
        if (m > (huge(m)-1)/capacity) then
            info = sparse_linalg_invalid_input
            return
        end if
        allocate(v(n,m+1), h(m*capacity+1,m*capacity), &
                 fh(m*capacity,m*capacity), seed(n), candidate(n))
        h = 0.0_dp
        seed = b/bnorm
        info = cap
        do ic = 1, cap
            if (ic > capacity) then
                next_capacity = int(min(int(cap,int64),2_int64*int(capacity,int64)))
                if (m > (huge(m)-1)/next_capacity) then
                    info = sparse_linalg_invalid_input
                    exit
                end if
                allocate(grown_h(m*next_capacity+1,m*next_capacity))
                grown_h = 0.0_dp
                grown_h(1:size(h,1),1:size(h,2)) = h
                call move_alloc(grown_h,h)
                deallocate(fh)
                allocate(fh(m*next_capacity,m*next_capacity))
                capacity = next_capacity
            end if
            base = (ic-1)*m
            call funm_krylov_cycle(n, matvec, seed, v, h, base, m, symmetric, &
                                  j, happy, nproducts, function_status)
            if (function_status /= sparse_linalg_success) then
                info = function_status
                exit
            end if
            global_j = base+j
            if (.not. all(ieee_is_finite(time_scale*h(1:global_j,1:global_j)))) then
                info = sparse_linalg_breakdown
                exit
            end if
            call f(time_scale*h(1:global_j,1:global_j), fh(1:global_j,1:global_j), function_status)
            if (function_status /= 0) then
                info = sparse_linalg_breakdown
                exit
            end if
            if (.not. all(ieee_is_finite(fh(1:global_j,1:global_j)))) then
                info = sparse_linalg_breakdown
                exit
            end if
            candidate = bnorm*matmul(v(:,1:j),fh(base+1:global_j,1))
            if (.not. all(ieee_is_finite(candidate))) then
                info = sparse_linalg_breakdown
                exit
            end if
            if (.not. all(ieee_is_finite(y+candidate))) then
                info = sparse_linalg_breakdown
                exit
            end if
            y = y+candidate
            change = norm2(candidate)
            nruns = ic
            if (happy .or. (ic > 1 .and. change <= threshold)) then
                info = sparse_linalg_success
                exit
            end if
            if (ic == cap) exit
            seed = v(:,m+1)
        end do
        if (present(cycles)) cycles = nruns
        if (present(products)) products = nproducts
        if (present(update_norm)) update_norm = change
    end subroutine funm_multiply_krylov_operator

    ! Construct one block of the extended block-lower-Hessenberg matrix.
    ! Each new basis starts from the previous block's last Arnoldi vector;
    ! the preceding coupling h(base+1,base) is retained in the global H.
    subroutine funm_krylov_cycle(n, matvec, seed, v, h, base, m, symmetric, &
                                used, happy, nproducts, info)
        integer, intent(in) :: n !! Number of rows of the operator.
        procedure(real_matvec_operator) :: matvec !! Sparse or matrix-free operator.
        real(dp), intent(in) :: seed(:) !! Start vector, usually already unit-norm.
        real(dp), intent(out) :: v(:,:) !! Local n-by-(m+1) Krylov basis.
        real(dp), intent(inout) :: h(:,:) !! Global extended Hessenberg matrix.
        integer, intent(in) :: base !! Number of columns in preceding blocks.
        integer, intent(in) :: m !! Maximum local Krylov order.
        logical, intent(in) :: symmetric !! Whether Lanczos is valid for A.
        integer, intent(out) :: used !! Number of columns actually constructed.
        logical, intent(out) :: happy !! Exact or numerically safe Arnoldi breakdown.
        integer, intent(inout) :: nproducts !! Accumulated matvec count.
        integer, intent(out) :: info !! Zero on success, negative on nonfinite products.

        integer :: k, j, pass
        real(dp) :: beta, alpha, projection, before
        real(dp) :: w(n)

        v = 0.0_dp
        v(:,1) = seed / norm2(seed)
        used = 0
        happy = .false.
        info = sparse_linalg_success
        do k = 1, m
            call matvec(v(:,k),w)
            nproducts = nproducts+1
            if (.not. all(ieee_is_finite(w))) then
                info = sparse_linalg_breakdown
                return
            end if
            before = norm2(w)
            if (symmetric) then
                if (k > 1) w = w-h(base+k,base+k-1)*v(:,k-1)
                alpha = dot_product(v(:,k),w)
                h(base+k,base+k) = alpha
                w = w-alpha*v(:,k)
            else
                ! Modified Gram--Schmidt twice reduces loss of orthogonality.
                do pass = 1,2
                    do j = 1,k
                        projection = dot_product(v(:,j),w)
                        h(base+j,base+k) = h(base+j,base+k)+projection
                        w = w-projection*v(:,j)
                    end do
                end do
            end if
            beta = norm2(w)
            if (.not. ieee_is_finite(beta)) then
                info = sparse_linalg_breakdown
                return
            end if
            used = k
            h(base+k+1,base+k) = beta
            happy = (beta <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,before)) .or. k == n
            if (happy) return
            v(:,k+1) = w/beta
            if (symmetric .and. k < m) h(base+k,base+k+1) = beta
        end do
    end subroutine funm_krylov_cycle

    subroutine expm_shift_norm(a, tracea, mu, operator_norm, info)
        type(csr_matrix), intent(in) :: a !! Prevalidated square CSR input.
        real(dp), intent(in), optional :: tracea !! Optional candidate trace shift.
        real(dp), intent(out) :: mu !! Selected shift, usually trace(A)/n.
        real(dp), intent(out) :: operator_norm !! Infinity norm of A - mu*I.
        integer, intent(out) :: info !! Zero if the shifted operator is valid.
        integer :: i, k, optional_info
        real(dp) :: exact_trace, compensation, addend, scratch
        real(dp) :: proposed_mu, proposed_norm

        info = sparse_linalg_success
        mu = 0.0_dp
        operator_norm = 0.0_dp
        if (present(tracea)) then
            if (.not. ieee_is_finite(tracea)) then
                info = sparse_linalg_invalid_input
                return
            end if
        end if
        if (a%n_rows == 0) return
        exact_trace = 0.0_dp
        compensation = 0.0_dp
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                if (a%indices(k) == i) then
                    addend = a%data(k)/real(a%n_rows,dp) - compensation
                    scratch = exact_trace + addend
                    compensation = (scratch - exact_trace) - addend
                    exact_trace = scratch
                end if
            end do
        end do
        if (.not. ieee_is_finite(exact_trace)) then
            info = sparse_linalg_breakdown
            return
        end if
        mu = exact_trace
        call csr_shift_infinity_norm(a,mu,operator_norm,info)
        if (info /= sparse_linalg_success) return
        if (present(tracea)) then
            proposed_mu = tracea/real(a%n_rows,dp)
            if (.not. ieee_is_finite(proposed_mu)) return
            call csr_shift_infinity_norm(a,proposed_mu,proposed_norm,optional_info)
            ! For an explicit sparse matrix the exact trace is inexpensive.
            ! Treat tracea as a *hint*, avoiding inaccurate user-supplied
            ! shifts that would increase cancellation and Taylor work.
            if (optional_info == sparse_linalg_success .and. &
                proposed_norm < operator_norm) then
                mu = proposed_mu
                operator_norm = proposed_norm
            end if
        end if
    end subroutine expm_shift_norm

    subroutine expm_multiply_matrix(a, b, y, info, t, tracea, tol, &
                                    scaling_steps, products)
        type(csr_matrix), intent(in) :: a !! Square, real, canonical CSR matrix.
        real(dp), intent(in) :: b(:,:) !! Finite RHS block, shape (n,nrhs).
        real(dp), allocatable, intent(out) :: y(:,:) !! exp(t*A)*b, shape (n,nrhs).
        integer, intent(out) :: info !! 0 success; -1 invalid; -2 overflow; 1 no Taylor convergence.
        real(dp), intent(in), optional :: t !! Scalar time, default 1.
        real(dp), intent(in), optional :: tracea !! Trace for shifting, default exact trace.
        real(dp), intent(in), optional :: tol !! Relative Taylor-tail target, default 4e-15.
        integer, intent(out), optional :: scaling_steps !! Number of scaling substeps.
        integer, intent(out), optional :: products !! Number of CSR times dense-block products.
        real(dp), allocatable :: term(:,:), sum_terms(:,:), product(:,:)
        real(dp) :: time, threshold, mu, norm_a, scaled_norm, h, shift_factor
        real(dp) :: beta, tail, ratio, denominator, term_norm, sum_norm
        integer :: n, n_rhs, s, step, j, nprod
        logical :: converged
        integer, parameter :: max_taylor_terms = 80
        real(dp), parameter :: max_step_norm = 2.0_dp
        integer, parameter :: max_scaling_steps = 200000

        info = sparse_linalg_invalid_input
        if (present(scaling_steps)) scaling_steps = 0
        if (present(products)) products = 0
        allocate(y(0,0))
        if (.not. valid_csr(a)) return
        n = a%n_rows
        n_rhs = size(b,2)
        if (a%n_cols /= n .or. size(b,1) /= n) return
        if (.not. all(ieee_is_finite(b))) return
        if (present(tracea)) then
            if (.not. ieee_is_finite(tracea)) return
        end if
        time = 1.0_dp
        if (present(t)) time = t
        threshold = 4.0e-15_dp
        if (present(tol)) threshold = tol
        if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(threshold)) return
        if (threshold < 0.0_dp .or. threshold >= 1.0_dp) return
        threshold = max(threshold,4.0_dp*epsilon(1.0_dp))
        deallocate(y)
        allocate(y(n,n_rhs))
        y = b
        info = sparse_linalg_success
        if (n == 0 .or. n_rhs == 0 .or. time == 0.0_dp) return
        if (all(b == 0.0_dp)) return

        call expm_shift_norm(a,tracea,mu,norm_a,info)
        if (info /= sparse_linalg_success) return
        if (norm_a == 0.0_dp) then
            shift_factor = exp(time*mu)
            if (.not. ieee_is_finite(shift_factor)) then
                info = sparse_linalg_breakdown
                return
            end if
            y = shift_factor*b
            if (.not. all(ieee_is_finite(y))) info = sparse_linalg_breakdown
            return
        end if
        if (abs(time) > huge(1.0_dp)/norm_a) then
            info = sparse_linalg_breakdown
            return
        end if
        scaled_norm = abs(time)*norm_a
        if (scaled_norm > real(max_scaling_steps,dp)*max_step_norm) then
            info = sparse_linalg_invalid_input
            return
        end if
        s = max(1,ceiling(scaled_norm/max_step_norm))
        h = time/real(s,dp)
        beta = abs(h)*norm_a
        shift_factor = exp(h*mu)
        if (.not. ieee_is_finite(shift_factor)) then
            info = sparse_linalg_breakdown
            return
        end if
        if (present(scaling_steps)) scaling_steps = s
        allocate(term(n,n_rhs),sum_terms(n,n_rhs),product(n,n_rhs))
        nprod = 0
        do step = 1,s
            term = y
            sum_terms = y
            converged = .false.
            do j = 1,max_taylor_terms
                call csr_block_product(a,term,product)
                nprod = nprod + 1
                term = (h/real(j,dp))*(product-mu*term)
                sum_terms = sum_terms+term
                if (.not. all(ieee_is_finite(sum_terms)) .or. &
                    .not. all(ieee_is_finite(term))) then
                    info = sparse_linalg_breakdown
                    exit
                end if
                ! A bound on the entire omitted remainder, using
                ! ||M^k||_infty <= ||M||_infty^k, M=h*(A-mu*I).
                ratio = beta/real(j+1,dp)
                denominator = 1.0_dp-beta/real(j+2,dp)
                if (denominator > 0.0_dp .and. ratio < 1.0_dp) then
                    term_norm = maxval(abs(term))
                    sum_norm = maxval(abs(sum_terms))
                    tail = term_norm*ratio/denominator
                    if (tail <= threshold*sum_norm) then
                        converged = .true.
                        exit
                    end if
                end if
            end do
            if (info /= sparse_linalg_success) exit
            if (.not. converged) then
                info = 1
                exit
            end if
            y = shift_factor*sum_terms
            if (.not. all(ieee_is_finite(y))) then
                info = sparse_linalg_breakdown
                exit
            end if
        end do
        if (present(products)) products = nprod
    end subroutine expm_multiply_matrix

    subroutine expm_multiply_vector(a, b, y, info, t, tracea, tol, &
                                    scaling_steps, products)
        type(csr_matrix), intent(in) :: a !! Square, real, canonical CSR matrix.
        real(dp), intent(in) :: b(:) !! Finite input vector of length n.
        real(dp), allocatable, intent(out) :: y(:) !! exp(t*A)*b, length n.
        integer, intent(out) :: info !! 0 success; -1 invalid; -2 overflow; 1 nonconvergence.
        real(dp), intent(in), optional :: t !! Scalar time, default 1.
        real(dp), intent(in), optional :: tracea !! Optional exact or approximate trace.
        real(dp), intent(in), optional :: tol !! Relative Taylor-tail target.
        integer, intent(out), optional :: scaling_steps !! Number of scaling substeps.
        integer, intent(out), optional :: products !! Number of sparse products.
        real(dp), allocatable :: block(:,:), result(:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call expm_multiply_matrix(a,block,result,info,t,tracea,tol,scaling_steps,products)
        if (size(result,2) == 1) then
            y = result(:,1)
        else
            allocate(y(0))
        end if
    end subroutine expm_multiply_vector

    subroutine expm_multiply_interval_matrix(a, b, start, stop, num, y, info, &
                                              endpoint, tracea, tol, products)
        type(csr_matrix), intent(in) :: a !! Square real CSR coefficient matrix.
        real(dp), intent(in) :: b(:,:) !! RHS block (n,nrhs).
        real(dp), intent(in) :: start !! First time in the grid.
        real(dp), intent(in) :: stop !! Last bound of the grid.
        integer, intent(in) :: num !! Number of grid samples, at least one.
        real(dp), allocatable, intent(out) :: y(:,:,:) !! Shape (n,nrhs,num).
        integer, intent(out) :: info !! Zero if all requested actions succeed.
        logical, intent(in), optional :: endpoint !! Include stop in the grid, default true.
        real(dp), intent(in), optional :: tracea !! Optional trace for scalar shifting.
        real(dp), intent(in), optional :: tol !! Relative Taylor-tail target.
        integer, intent(out), optional :: products !! Total CSR times block products.
        real(dp), allocatable :: current(:,:), next(:,:)
        integer :: i, anchor, count_products, nproducts, denominator
        real(dp) :: delta, grid_time, anchor_time
        logical :: include_endpoint

        info = sparse_linalg_invalid_input
        if (present(products)) products = 0
        allocate(y(0,0,0))
        if (.not. valid_csr(a)) return
        if (a%n_rows /= a%n_cols .or. size(b,1) /= a%n_rows) return
        if (.not. all(ieee_is_finite(b))) return
        if (num < 1 .or. .not. ieee_is_finite(start) .or. &
            .not. ieee_is_finite(stop)) return
        include_endpoint = .true.
        if (present(endpoint)) include_endpoint = endpoint
        delta = 0.0_dp
        if (num > 1) then
            denominator = num
            if (include_endpoint) denominator = num-1
            delta = (stop-start)/real(denominator,dp)
            if (.not. ieee_is_finite(delta)) return
        end if
        ! Propagate from the point closest to t=0, in both directions.
        ! Otherwise, an exp(t*A)*b that underflows at a remote first time
        ! cannot recover when the interval later approaches zero.
        anchor = 1
        anchor_time = start
        do i = 2,num
            grid_time = start+real(i-1,dp)*delta
            if (.not. ieee_is_finite(grid_time)) return
            if (abs(grid_time) < abs(anchor_time)) then
                anchor = i
                anchor_time = grid_time
            end if
        end do
        call expm_multiply_matrix(a,b,current,info,t=anchor_time,tracea=tracea, &
            tol=tol,products=nproducts)
        if (info /= sparse_linalg_success) return
        count_products = nproducts
        deallocate(y)
        allocate(y(a%n_rows,size(b,2),num))
        y = 0.0_dp
        y(:,:,anchor) = current
        do i = anchor+1,num
            call expm_multiply_matrix(a,current,next,info,t=delta,tracea=tracea, &
                tol=tol,products=nproducts)
            count_products = count_products+nproducts
            if (info /= sparse_linalg_success) exit
            y(:,:,i) = next
            call move_alloc(next,current)
        end do
        if (info == sparse_linalg_success) then
            current = y(:,:,anchor)
            do i = anchor-1,1,-1
                call expm_multiply_matrix(a,current,next,info,t=-delta,tracea=tracea, &
                    tol=tol,products=nproducts)
                count_products = count_products+nproducts
                if (info /= sparse_linalg_success) exit
                y(:,:,i) = next
                call move_alloc(next,current)
            end do
        end if
        if (present(products)) products = count_products
    end subroutine expm_multiply_interval_matrix

    subroutine expm_multiply_interval_vector(a, b, start, stop, num, y, info, &
                                              endpoint, tracea, tol, products)
        type(csr_matrix), intent(in) :: a !! Square real CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! RHS vector of length n.
        real(dp), intent(in) :: start !! First requested time.
        real(dp), intent(in) :: stop !! Last requested bound.
        integer, intent(in) :: num !! Requested number of output times.
        real(dp), allocatable, intent(out) :: y(:,:) !! Shape (n,num).
        integer, intent(out) :: info !! Zero if all time points succeed.
        logical, intent(in), optional :: endpoint !! Include stop, default true.
        real(dp), intent(in), optional :: tracea !! Optional trace for scalar shifting.
        real(dp), intent(in), optional :: tol !! Relative Taylor-tail target.
        integer, intent(out), optional :: products !! Total sparse products.
        real(dp), allocatable :: block(:,:), result(:,:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call expm_multiply_interval_matrix(a,block,start,stop,num,result,info, &
            endpoint,tracea,tol,products)
        if (size(result,2) == 1) then
            y = result(:,1,:)
        else
            allocate(y(0,0))
        end if
    end subroutine expm_multiply_interval_vector


    ! Deterministic block power/adjoint estimator for the induced matrix 1-norm.
    ! All probes have unit 1-norm; thus every recorded value is a lower bound.
    ! This uses the sign/adjoint-direction principle of Higham--Tisseur (2000).
    ! For a materialized CSR input the exact norm can also be found in O(nnz),
    ! but this estimator is groundwork for future matrix-free operators.
    function onenormest(a, t, itmax, v, w, info) result(estimate)
        type(csr_matrix), intent(in) :: a !! Square, real, canonical one-based CSR matrix.
        integer, intent(in), optional :: t !! Number of probe vectors, default 2, >= 1.
        integer, intent(in), optional :: itmax !! Sign/adjoint iterations, default 5, >= 1.
        real(dp), allocatable, intent(out), optional :: v(:) !! Witness with ||v||_1=1.
        real(dp), allocatable, intent(out), optional :: w(:) !! A*v; ||w||_1=estimate.
        integer, intent(out), optional :: info !! 0 success, negative invalid input/breakdown.
        real(dp) :: estimate
        real(dp), allocatable :: x(:,:), y(:,:), signs(:,:), z(:,:), h(:), best_v(:), best_w(:)
        real(dp) :: candidate, highest
        integer, allocatable :: selected(:)
        logical, allocatable :: visited(:)
        integer :: n, width, cycles, i, j, k, idx, step, status

        estimate = 0.0_dp
        status = sparse_linalg_invalid_input
        if (present(info)) info = status
        if (present(v)) allocate(v(0))
        if (present(w)) allocate(w(0))
        if (.not. valid_csr(a)) return
        n = a%n_cols
        if (a%n_rows /= n) return
        width = 2
        if (present(t)) width = t
        cycles = 5
        if (present(itmax)) cycles = itmax
        if (width <= 0 .or. cycles <= 0) return
        if (n == 0) then
            status = sparse_linalg_success
            if (present(info)) info = status
            return
        end if
        width = min(width,n)
        allocate(h(n),best_v(n),best_w(n))
        best_v = 0.0_dp
        best_v(1) = 1.0_dp
        best_w = 0.0_dp
        if (width == n) then
            ! A full set of coordinate probes returns the exact 1-norm in
            ! O(nnz) time; there is no reason to allocate an n-by-n block.
            h = 0.0_dp
            do i = 1,n
                do k = a%indptr(i),a%indptr(i+1)-1
                    j = a%indices(k)
                    h(j) = h(j)+abs(a%data(k))
                end do
            end do
            if (.not. all(ieee_is_finite(h))) then
                status = sparse_linalg_breakdown
            else
                idx = maxloc(h,dim=1)
                estimate = h(idx)
                best_v = 0.0_dp
                best_v(idx) = 1.0_dp
                call csr_product(a,best_v,best_w)
                status = sparse_linalg_success
            end if
        else
            allocate(x(n,width),y(n,width),signs(n,width),z(n,width), &
                selected(width),visited(n))
            visited = .false.
            ! Start with a normalized all-ones probe and deterministic signed
            ! companions. Subsequent probes are coordinate vectors selected
            ! from the adjoint action of the matrix on the output signs.
            x(:,1) = 1.0_dp/real(n,dp)
            do j = 2,width
                do i = 1,n
                    k = mod(7*mod(i,17)+11*mod(j,19)+ &
                        3*mod(i,23)*mod(j,13),19)
                    x(i,j) = merge(1.0_dp,-1.0_dp,k < 9)/real(n,dp)
                end do
            end do
            status = sparse_linalg_success
            do step = 1,cycles
                call csr_block_product(a,x,y)
                if (.not. all(ieee_is_finite(y))) then
                    status = sparse_linalg_breakdown
                    exit
                end if
                do j = 1,width
                    candidate = sum(abs(y(:,j)))
                    if (.not. ieee_is_finite(candidate)) then
                        status = sparse_linalg_breakdown
                        exit
                    end if
                    if (candidate > estimate .or. (step == 1 .and. j == 1)) then
                        estimate = candidate
                        best_v = x(:,j)
                        best_w = y(:,j)
                    end if
                end do
                if (status /= sparse_linalg_success) exit
                if (step == cycles) exit
                signs = merge(1.0_dp,-1.0_dp,y >= 0.0_dp)
                z = 0.0_dp
                do i = 1,n
                    do k = a%indptr(i),a%indptr(i+1)-1
                        z(a%indices(k),:) = z(a%indices(k),:) + a%data(k)*signs(i,:)
                    end do
                end do
                if (.not. all(ieee_is_finite(z))) then
                    status = sparse_linalg_breakdown
                    exit
                end if
                h = maxval(abs(z),dim=2)
                ! This gives a safe stopping test: no adjoint direction can
                ! improve the current estimate significantly.
                if (maxval(h) <= estimate*(1.0_dp+32.0_dp*epsilon(estimate))) exit
                selected = 0
                do j = 1,width
                    idx = 0
                    highest = -1.0_dp
                    do i = 1,n
                        if (visited(i)) cycle
                        if (h(i) > highest) then
                            highest = h(i)
                            idx = i
                        end if
                    end do
                    if (idx == 0) exit
                    selected(j) = idx
                    visited(idx) = .true.
                end do
                if (all(selected == 0)) exit
                x = 0.0_dp
                do j = 1,width
                    if (selected(j) /= 0) then
                        x(selected(j),j) = 1.0_dp
                    else
                        x(1,j) = 1.0_dp
                    end if
                end do
            end do
        end if
        if (status == sparse_linalg_success) then
            if (present(v)) v = best_v
            if (present(w)) w = best_w
        end if
        if (present(info)) info = status
    end function onenormest

    ! Form the sparse exponential by applying the already-tested Taylor
    ! exponential action to blocks of identity columns.  Unlike SciPy's
    ! sparse expm this is not a rational Pade method; the result may have
    ! significant fill-in and cost O(n) exponential actions.
    function expm(a, info, drop_tol, block_size, tol) result(result)
        type(csr_matrix), intent(in) :: a !! Square one-based real CSR operator.
        integer, intent(out), optional :: info !! 0 success; -1 invalid; -2 overflow; 1 no convergence.
        real(dp), intent(in), optional :: drop_tol !! Absolute sparsification cutoff, default 0.
        integer, intent(in), optional :: block_size !! Identity columns per action call, default 8.
        real(dp), intent(in), optional :: tol !! Relative exponential-action tolerance.
        type(csr_matrix) :: result
        integer, allocatable :: rows(:), cols(:), row_tmp(:), col_tmp(:)
        real(dp), allocatable :: values(:), val_tmp(:), rhs(:,:), output(:,:), diag(:)
        real(dp) :: cutoff, e
        integer :: status, n, i, j, k, col, first, nb, width, used, capacity, newcap
        integer :: sparsestatus
        integer(int64) :: max_entries, next_capacity
        logical :: diagonal

        status = sparse_linalg_invalid_input
        if (present(info)) info = status
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (a%n_cols /= n) return
        cutoff = 0.0_dp
        if (present(drop_tol)) cutoff = drop_tol
        width = 8
        if (present(block_size)) width = block_size
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        if (width <= 0) return
        if (present(tol)) then
            if (.not. ieee_is_finite(tol)) return
            if (tol < 0.0_dp .or. tol >= 1.0_dp) return
        end if
        if (n == 0) then
            result = csr_eye(0)
            if (present(info)) info = sparse_linalg_success
            return
        end if

        ! Diagonal inputs are handled in O(n+nnz) rather than using n RHS
        ! columns; stored off-diagonal zeros are ignored here.
        diagonal = .true.
        allocate(diag(n))
        diag = 0.0_dp
        do i = 1,n
            do k = a%indptr(i),a%indptr(i+1)-1
                if (a%indices(k) == i) then
                    diag(i) = a%data(k)
                else if (a%data(k) /= 0.0_dp) then
                    diagonal = .false.
                end if
            end do
        end do
        if (diagonal) then
            allocate(rows(n),cols(n),values(n))
            used = 0
            do i = 1,n
                e = exp(diag(i))
                if (.not. ieee_is_finite(e)) then
                    if (present(info)) info = sparse_linalg_breakdown
                    return
                end if
                if (abs(e) <= cutoff) cycle
                used = used+1
                rows(used) = i
                cols(used) = i
                values(used) = e
            end do
            result = csr_from_coo(n,n,rows(:used),cols(:used),values(:used),sparsestatus)
            if (sparsestatus /= sparse_success) then
                if (present(info)) info = sparse_linalg_breakdown
            else
                if (present(info)) info = sparse_linalg_success
            end if
            return
        end if

        max_entries = min(int(n,int64)*int(n,int64),int(huge(n)-1,int64))
        capacity = int(min(max_entries,max(int(n,int64), &
            2_int64*int(size(a%data),int64))))
        allocate(rows(capacity),cols(capacity),values(capacity))
        width = min(width,n)
        allocate(rhs(n,width))
        used = 0
        status = sparse_linalg_success
        do first = 1,n,width
            nb = min(width,n-first+1)
            rhs = 0.0_dp
            do j = 1,nb
                rhs(first+j-1,j) = 1.0_dp
            end do
            call expm_multiply(a,rhs(:,:nb),output,status,tol=tol)
            if (status /= sparse_linalg_success) exit
            do j = 1,nb
                col = first+j-1
                do i = 1,n
                    if (abs(output(i,j)) <= cutoff) cycle
                    if (used == capacity) then
                        next_capacity = min(max_entries,max(int(capacity,int64)+1_int64, &
                            2_int64*int(capacity,int64)))
                        if (next_capacity <= int(capacity,int64)) then
                            status = sparse_linalg_breakdown
                            exit
                        end if
                        newcap = int(next_capacity)
                        allocate(row_tmp(newcap),col_tmp(newcap),val_tmp(newcap))
                        row_tmp(:used) = rows(:used)
                        col_tmp(:used) = cols(:used)
                        val_tmp(:used) = values(:used)
                        call move_alloc(row_tmp,rows)
                        call move_alloc(col_tmp,cols)
                        call move_alloc(val_tmp,values)
                        capacity = newcap
                    end if
                    used = used+1
                    rows(used) = i
                    cols(used) = col
                    values(used) = output(i,j)
                end do
                if (status /= sparse_linalg_success) exit
            end do
            if (status /= sparse_linalg_success) exit
        end do
        if (status == sparse_linalg_success) then
            result = csr_from_coo(n,n,rows(:used),cols(:used),values(:used),sparsestatus)
            if (sparsestatus /= sparse_success) status = sparse_linalg_breakdown
        end if
        if (present(info)) info = status
    end function expm

    ! A sparse LU factor stores P*A=L*U in a single canonical CSR pattern:
    ! entries below the diagonal belong to unit-diagonal L; entries on or
    ! above the diagonal belong to U.  perm(i) is the original RHS row at i.
    ! The row permutation is kept explicitly, unlike ILU(0) (no pivoting).
    
    ! CSC inputs are converted to canonical CSR once, so the same triangular
    ! and pivoted direct kernels serve both sparse storage formats.
    subroutine spsolve_triangular_csc_vector(a, b, x, info, lower, unit_diagonal)
        type(csc_matrix), intent(in) :: a !! Canonical CSC triangular operator.
        real(dp), intent(in) :: b(:) !! Finite right-hand side vector.
        real(dp), allocatable, intent(out) :: x(:) !! Computed solution.
        integer, intent(out) :: info !! Zero success; negative error status.
        logical, intent(in), optional :: lower !! Select lower vs upper triangle.
        logical, intent(in), optional :: unit_diagonal !! Treat diagonal as all ones.
        type(csr_matrix) :: csr
        integer :: status

        csr = csc_to_csr(a,status)
        if (status /= sparse_success) then
            info = sparse_linalg_invalid_input
            allocate(x(0))
            return
        end if
        call spsolve_triangular_vector(csr,b,x,info,lower,unit_diagonal)
    end subroutine spsolve_triangular_csc_vector

    subroutine spsolve_triangular_csc_matrix(a, b, x, info, lower, unit_diagonal)
        type(csc_matrix), intent(in) :: a !! Canonical CSC triangular operator.
        real(dp), intent(in) :: b(:,:) !! Finite block right-hand side.
        real(dp), allocatable, intent(out) :: x(:,:) !! Computed solution block.
        integer, intent(out) :: info !! Zero success; negative error status.
        logical, intent(in), optional :: lower !! Select lower vs upper triangle.
        logical, intent(in), optional :: unit_diagonal !! Treat diagonal as all ones.
        type(csr_matrix) :: csr
        integer :: status

        csr = csc_to_csr(a,status)
        if (status /= sparse_success) then
            info = sparse_linalg_invalid_input
            allocate(x(0,0))
            return
        end if
        call spsolve_triangular_matrix(csr,b,x,info,lower,unit_diagonal)
    end subroutine spsolve_triangular_csc_matrix

    subroutine splu_csc(a, factor, info, pivot_tol)
        type(csc_matrix), intent(in) :: a !! Square canonical CSC matrix.
        type(splu_factor), intent(out) :: factor !! Sparse pivoted factors.
        integer, intent(out) :: info !! Zero success; negative error status.
        real(dp), intent(in), optional :: pivot_tol !! Minimum absolute pivot.
        type(csr_matrix) :: csr
        integer :: status

        csr = csc_to_csr(a,status)
        if (status /= sparse_success) then
            info = sparse_linalg_invalid_input
            return
        end if
        call splu_csr(csr,factor,info,pivot_tol)
    end subroutine splu_csc

    subroutine spsolve_csc_vector(a, b, x, info, pivot_tol)
        type(csc_matrix), intent(in) :: a !! Square canonical CSC matrix.
        real(dp), intent(in) :: b(:) !! Finite right-hand side vector.
        real(dp), allocatable, intent(out) :: x(:) !! Computed solution.
        integer, intent(out) :: info !! Zero success; negative error status.
        real(dp), intent(in), optional :: pivot_tol !! Minimum absolute pivot.
        type(splu_factor) :: factor

        call splu_csc(a,factor,info,pivot_tol)
        if (info /= sparse_linalg_success) then
            allocate(x(0))
            return
        end if
        call splu_solve_vector(factor,b,x,info)
    end subroutine spsolve_csc_vector

    subroutine spsolve_csc_matrix(a, b, x, info, pivot_tol)
        type(csc_matrix), intent(in) :: a !! Square canonical CSC matrix.
        real(dp), intent(in) :: b(:,:) !! Finite right-hand side matrix.
        real(dp), allocatable, intent(out) :: x(:,:) !! Computed solution block.
        integer, intent(out) :: info !! Zero success; negative error status.
        real(dp), intent(in), optional :: pivot_tol !! Minimum absolute pivot.
        type(splu_factor) :: factor

        call splu_csc(a,factor,info,pivot_tol)
        if (info /= sparse_linalg_success) then
            allocate(x(0,0))
            return
        end if
        call splu_solve_matrix(factor,b,x,info)
    end subroutine spsolve_csc_matrix

    ! Solve triangular CSR matrices, vector and multi-RHS overloads.
    subroutine spsolve_triangular_vector(a, b, x, info, lower, unit_diagonal)
        type(csr_matrix), intent(in) :: a !! Canonical square CSR matrix.
        real(dp), intent(in) :: b(:) !! Finite RHS vector, length n.
        real(dp), allocatable, intent(out) :: x(:) !! Solution, length n on success.
        integer, intent(out) :: info !! 0 success, -1 invalid, -3 zero/missing pivot, -2 overflow.
        logical, intent(in), optional :: lower !! Solve lower triangle (default true).
        logical, intent(in), optional :: unit_diagonal !! Assume all diagonal entries are 1.
        real(dp), allocatable :: block(:,:)
        real(dp), allocatable :: rhs(:,:)

        allocate(rhs(size(b),1))
        rhs(:,1) = b
        call spsolve_triangular_matrix(a,rhs,block,info,lower,unit_diagonal)
        if (info == sparse_linalg_success) then
            x = block(:,1)
        else
            allocate(x(0))
        end if
    end subroutine spsolve_triangular_vector

    subroutine spsolve_triangular_matrix(a, b, x, info, lower, unit_diagonal)
        type(csr_matrix), intent(in) :: a !! Canonical square CSR matrix.
        real(dp), intent(in) :: b(:,:) !! Finite RHS matrix with n rows.
        real(dp), allocatable, intent(out) :: x(:,:) !! Solution with shape b on success.
        integer, intent(out) :: info !! 0 success, -1 invalid, -3 zero/missing pivot, -2 overflow.
        logical, intent(in), optional :: lower !! Use lower triangle, default true.
        logical, intent(in), optional :: unit_diagonal !! Assume unit diagonal, default false.
        integer :: i, j, k, n, first, last, increment
        real(dp) :: diagonal
        logical :: is_lower, is_unit, has_diagonal

        info = sparse_linalg_invalid_input
        allocate(x(0,0))
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (a%n_cols /= n .or. size(b,1) /= n) return
        if (.not. all(ieee_is_finite(b))) return
        is_lower = .true.
        if (present(lower)) is_lower = lower
        is_unit = .false.
        if (present(unit_diagonal)) is_unit = unit_diagonal
        deallocate(x)
        allocate(x(n,size(b,2)))
        x = b
        info = sparse_linalg_success
        if (is_lower) then
            first = 1
            last = n
            increment = 1
        else
            first = n
            last = 1
            increment = -1
        end if
        do i = first, last, increment
            diagonal = 1.0_dp
            has_diagonal = is_unit
            do k = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(k)
                if (j == i .and. .not. is_unit) then
                    diagonal = a%data(k)
                    has_diagonal = .true.
                else if ((is_lower .and. j < i) .or. &
                    (.not. is_lower .and. j > i)) then
                    x(i,:) = x(i,:) - a%data(k)*x(j,:)
                end if
            end do
            if (.not. has_diagonal .or. diagonal == 0.0_dp) then
                info = sparse_linalg_bad_diagonal
                return
            end if
            x(i,:) = x(i,:)/diagonal
            if (.not. all(ieee_is_finite(x(i,:)))) then
                info = sparse_linalg_breakdown
                return
            end if
        end do
    end subroutine spsolve_triangular_matrix

    ! Binary search within a mutable row.  The negative return value is the
    ! insertion position, with 1-based columns sorted strictly increasing.
    pure integer function sparse_row_position(row, column) result(pos)
        type(sparse_lu_row), intent(in) :: row !! Sorted mutable sparse row.
        integer, intent(in) :: column !! One-based desired column.
        integer :: left, right, mid

        left = 1
        right = row%used
        do while (left <= right)
            mid = left + (right-left)/2
            if (row%columns(mid) == column) then
                pos = mid
                return
            else if (row%columns(mid) < column) then
                left = mid+1
            else
                right = mid-1
            end if
        end do
        pos = -left
    end function sparse_row_position

    pure real(dp) function sparse_row_value(row, column) result(value)
        type(sparse_lu_row), intent(in) :: row !! Row to search.
        integer, intent(in) :: column !! One-based column to read.
        integer :: pos

        pos = sparse_row_position(row,column)
        value = 0.0_dp
        if (pos > 0) value = row%values(pos)
    end function sparse_row_value

    subroutine sparse_row_set(row, column, value)
        type(sparse_lu_row), intent(inout) :: row !! Row updated in place.
        integer, intent(in) :: column !! Column to replace or insert.
        real(dp), intent(in) :: value !! Finite value to store; exact zero deletes.
        integer :: pos, capacity, next_capacity
        integer, allocatable :: temp_columns(:)
        real(dp), allocatable :: temp_values(:)

        pos = sparse_row_position(row,column)
        if (pos > 0) then
            if (value == 0.0_dp) then
                row%columns(pos:row%used-1) = row%columns(pos+1:row%used)
                row%values(pos:row%used-1) = row%values(pos+1:row%used)
                row%used = row%used-1
            else
                row%values(pos) = value
            end if
            return
        end if
        if (value == 0.0_dp) return
        pos = -pos
        capacity = size(row%columns)
        if (row%used == capacity) then
            if (capacity > huge(capacity)/2) error stop 'sparse LU row capacity overflow'
            next_capacity = max(4,2*capacity)
            allocate(temp_columns(next_capacity),temp_values(next_capacity))
            temp_columns(:row%used) = row%columns(:row%used)
            temp_values(:row%used) = row%values(:row%used)
            call move_alloc(temp_columns,row%columns)
            call move_alloc(temp_values,row%values)
        end if
        row%columns(pos+1:row%used+1) = row%columns(pos:row%used)
        row%values(pos+1:row%used+1) = row%values(pos:row%used)
        row%columns(pos) = column
        row%values(pos) = value
        row%used = row%used+1
    end subroutine sparse_row_set

    ! Gaussian elimination with partial row pivoting and sparse, growable
    ! rows.  This is a direct factorization, not ILU(0): necessary fill-in is
    ! retained.  Complexity depends on pivot fill; pathological cases can be
    ! expensive, just as for other non-multifrontal sparse LU algorithms.
    subroutine splu_csr(a, factor, info, pivot_tol)
        type(csr_matrix), intent(in) :: a !! Square canonical real CSR matrix.
        type(splu_factor), intent(out) :: factor !! LU factors and row permutation.
        integer, intent(out) :: info !! 0 success; -1 invalid; -4 singular; -2 overflow.
        real(dp), intent(in), optional :: pivot_tol !! Nonnegative minimum absolute pivot; default 0.
        type(sparse_lu_row), allocatable :: rows(:)
        type(sparse_lu_row) :: swaprow
        integer :: n, i, j, k, p, idx, nnz, ofs, tmp_perm
        real(dp) :: tol, pivot, candidate, multiplier, value, best

        info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (a%n_cols /= n) return
        tol = 0.0_dp
        if (present(pivot_tol)) tol = pivot_tol
        if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) return
        allocate(rows(n),factor%row_perm(n))
        factor%row_perm = [(i,i=1,n)]
        do i = 1,n
            nnz = a%indptr(i+1)-a%indptr(i)
            allocate(rows(i)%columns(nnz),rows(i)%values(nnz))
            rows(i)%used = nnz
            rows(i)%columns = a%indices(a%indptr(i):a%indptr(i+1)-1)
            rows(i)%values = a%data(a%indptr(i):a%indptr(i+1)-1)
        end do

        do k = 1,n
            p = k
            best = abs(sparse_row_value(rows(k),k))
            do i = k+1,n
                candidate = abs(sparse_row_value(rows(i),k))
                if (candidate > best) then
                    p = i
                    best = candidate
                end if
            end do
            if (best <= tol) then
                info = sparse_linalg_singular
                return
            end if
            if (p /= k) then
                swaprow = rows(k)
                rows(k) = rows(p)
                rows(p) = swaprow
                tmp_perm = factor%row_perm(k)
                factor%row_perm(k) = factor%row_perm(p)
                factor%row_perm(p) = tmp_perm
            end if
            pivot = sparse_row_value(rows(k),k)
            do i = k+1,n
                value = sparse_row_value(rows(i),k)
                if (value == 0.0_dp) cycle
                multiplier = value/pivot
                if (.not. ieee_is_finite(multiplier)) then
                    info = sparse_linalg_breakdown
                    return
                end if
                call sparse_row_set(rows(i),k,multiplier)
                do idx = 1,rows(k)%used
                    j = rows(k)%columns(idx)
                    if (j <= k) cycle
                    value = sparse_row_value(rows(i),j) - multiplier*rows(k)%values(idx)
                    if (.not. ieee_is_finite(value)) then
                        info = sparse_linalg_breakdown
                        return
                    end if
                    call sparse_row_set(rows(i),j,value)
                end do
            end do
        end do

        nnz = 0
        do i = 1,n
            if (rows(i)%used > huge(nnz)-nnz) then
                info = sparse_linalg_breakdown
                return
            end if
            nnz = nnz + rows(i)%used
        end do
        factor%lu%n_rows = n
        factor%lu%n_cols = n
        allocate(factor%lu%indptr(n+1),factor%lu%indices(nnz),factor%lu%data(nnz))
        ofs = 1
        do i = 1,n
            factor%lu%indptr(i) = ofs
            j = rows(i)%used
            factor%lu%indices(ofs:ofs+j-1) = rows(i)%columns(:j)
            factor%lu%data(ofs:ofs+j-1) = rows(i)%values(:j)
            ofs = ofs+j
        end do
        factor%lu%indptr(n+1) = ofs
        info = sparse_linalg_success
    end subroutine splu_csr

    pure logical function valid_splu(factor) result(ok)
        type(splu_factor), intent(in) :: factor !! Candidate pivoted LU factor.
        integer :: i, k, n
        logical, allocatable :: seen(:)
        logical :: diagonal_present

        ok = .false.
        if (.not. valid_csr(factor%lu)) return
        n = factor%lu%n_rows
        if (factor%lu%n_cols /= n) return
        if (.not. allocated(factor%row_perm)) return
        if (size(factor%row_perm) /= n) return
        if (any(factor%row_perm < 1) .or. any(factor%row_perm > n)) return
        allocate(seen(n))
        seen = .false.
        do i = 1,n
            if (seen(factor%row_perm(i))) return
            seen(factor%row_perm(i)) = .true.
            diagonal_present = .false.
            do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                if (factor%lu%indices(k) == i .and. factor%lu%data(k) /= 0.0_dp) diagonal_present = .true.
            end do
            if (.not. diagonal_present) return
        end do
        ok = .true.
    end function valid_splu

    subroutine splu_solve_vector(factor, b, x, info)
        type(splu_factor), intent(in) :: factor !! Valid pivoted sparse LU.
        real(dp), intent(in) :: b(:) !! Finite RHS vector.
        real(dp), allocatable, intent(out) :: x(:) !! Solution on success.
        integer, intent(out) :: info !! 0 success; -1 invalid; -2 overflow.
        real(dp), allocatable :: rhs(:,:), block(:,:)

        allocate(rhs(size(b),1))
        rhs(:,1) = b
        call splu_solve_matrix(factor,rhs,block,info)
        if (info == sparse_linalg_success) then
            x = block(:,1)
        else
            allocate(x(0))
        end if
    end subroutine splu_solve_vector

    subroutine splu_solve_matrix(factor, b, x, info)
        type(splu_factor), intent(in) :: factor !! Valid pivoted sparse LU.
        real(dp), intent(in) :: b(:,:) !! Finite RHS block (n,nrhs).
        real(dp), allocatable, intent(out) :: x(:,:) !! Solution block.
        integer, intent(out) :: info !! 0 success; -1 invalid; -2 overflow.
        integer :: n, i, j, k
        real(dp) :: diagonal

        info = sparse_linalg_invalid_input
        allocate(x(0,0))
        if (.not. valid_splu(factor)) return
        n = factor%lu%n_rows
        if (size(b,1) /= n) return
        if (.not. all(ieee_is_finite(b))) return
        deallocate(x)
        allocate(x(n,size(b,2)))
        do i = 1,n
            x(i,:) = b(factor%row_perm(i),:)
        end do
        ! Forward substitution with implicit unit diagonal of L.
        do i = 1,n
            do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                j = factor%lu%indices(k)
                if (j >= i) exit
                x(i,:) = x(i,:) - factor%lu%data(k)*x(j,:)
            end do
        end do
        ! Back substitution with the nonunit upper factor U.
        do i = n,1,-1
            diagonal = 0.0_dp
            do k = factor%lu%indptr(i),factor%lu%indptr(i+1)-1
                j = factor%lu%indices(k)
                if (j == i) diagonal = factor%lu%data(k)
                if (j > i) x(i,:) = x(i,:) - factor%lu%data(k)*x(j,:)
            end do
            x(i,:) = x(i,:)/diagonal
            if (.not. all(ieee_is_finite(x(i,:)))) then
                info = sparse_linalg_breakdown
                return
            end if
        end do
        info = sparse_linalg_success
    end subroutine splu_solve_matrix

    subroutine spsolve_vector(a, b, x, info, pivot_tol)
        type(csr_matrix), intent(in) :: a !! Square CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! RHS vector.
        real(dp), allocatable, intent(out) :: x(:) !! Pivoted direct solution.
        integer, intent(out) :: info !! 0 success; -1 invalid; -4 singular; -2 overflow.
        real(dp), intent(in), optional :: pivot_tol !! Optional minimum absolute pivot, default 0.
        type(splu_factor) :: factor

        call splu(a,factor,info,pivot_tol)
        if (info /= sparse_linalg_success) then
            allocate(x(0))
            return
        end if
        call splu_solve(factor,b,x,info)
    end subroutine spsolve_vector

    subroutine spsolve_matrix(a, b, x, info, pivot_tol)
        type(csr_matrix), intent(in) :: a !! Square CSR coefficient matrix.
        real(dp), intent(in) :: b(:,:) !! Multiple RHS, shape (n,nrhs).
        real(dp), allocatable, intent(out) :: x(:,:) !! Pivoted direct solution.
        integer, intent(out) :: info !! 0 success; -1 invalid; -4 singular; -2 overflow.
        real(dp), intent(in), optional :: pivot_tol !! Optional minimum absolute pivot, default 0.
        type(splu_factor) :: factor

        call splu(a,factor,info,pivot_tol)
        if (info /= sparse_linalg_success) then
            allocate(x(0,0))
            return
        end if
        call splu_solve(factor,b,x,info)
    end subroutine spsolve_matrix


    ! Repeated squaring uses O(log(power)) CSR products.  Negative powers are
    ! not defined for sparse matrices in SciPy's matrix_power interface.
    function matrix_power_csr(a, power, info) result(out)
        type(csr_matrix), intent(in) :: a !! Square, canonical one-based CSR input.
        integer, intent(in) :: power !! Integer exponent, >= 0.
        integer, intent(out), optional :: info !! 0 success; -1 invalid; -2 overflow.
        type(csr_matrix) :: out
        type(csr_matrix) :: base, next
        integer :: k, ierr

        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        if (a%n_rows /= a%n_cols .or. power < 0) return
        if (power == 0) then
            out = csr_eye(a%n_rows)
            if (present(info)) info = sparse_linalg_success
            return
        end if
        if (power == 1) then
            out = a
            if (present(info)) info = sparse_linalg_success
            return
        end if

        out = csr_eye(a%n_rows)
        base = a
        k = power
        do while (k > 0)
            if (btest(k,0)) then
                next = csr_matmul(out,base,ierr)
                if (ierr /= sparse_success) then
                    out = csr_matrix()
                    if (present(info)) info = sparse_linalg_breakdown
                    return
                end if
                out = next
            end if
            k = ishft(k,-1)
            if (k == 0) exit
            next = csr_matmul(base,base,ierr)
            if (ierr /= sparse_success) then
                out = csr_matrix()
                if (present(info)) info = sparse_linalg_breakdown
                return
            end if
            base = next
        end do
        if (present(info)) info = sparse_linalg_success
    end function matrix_power_csr

    function matrix_power_csc(a, power, info) result(out)
        type(csc_matrix), intent(in) :: a !! Square, canonical one-based CSC input.
        integer, intent(in) :: power !! Integer exponent, >= 0.
        integer, intent(out), optional :: info !! 0 success; -1 invalid; -2 overflow.
        type(csc_matrix) :: out
        type(csr_matrix) :: csr, powered
        integer :: ierr

        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csc(a)) return
        if (a%n_rows /= a%n_cols .or. power < 0) return
        csr = csc_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        powered = matrix_power_csr(csr,power,ierr)
        if (ierr /= sparse_linalg_success) then
            if (present(info)) info = ierr
            return
        end if
        out = csr_to_csc(powered,ierr)
        if (ierr == sparse_success) then
            if (present(info)) info = sparse_linalg_success
        else
            if (present(info)) info = sparse_linalg_breakdown
        end if
    end function matrix_power_csc

    ! Sparse inverse uses one sparse LU factorization with multiple identity
    ! right-hand sides.  The inverse itself may be mathematically dense.
    function inv_csr(a, info, drop_tol, block_size, pivot_tol) result(out)
        type(csr_matrix), intent(in) :: a !! Nonsingular square CSR matrix.
        integer, intent(out), optional :: info !! 0 success; -1 invalid; -4 singular; -2 breakdown.
        real(dp), intent(in), optional :: drop_tol !! Absolute cutoff for inverse entries, default 0.
        integer, intent(in), optional :: block_size !! Number of identity columns per solve, default 8.
        real(dp), intent(in), optional :: pivot_tol !! Minimum pivot magnitude, default 0.
        type(csr_matrix) :: out
        type(splu_factor) :: factor
        integer, allocatable :: rows(:), cols(:), new_rows(:), new_cols(:)
        real(dp), allocatable :: values(:), new_values(:), diagonal(:)
        real(dp), allocatable :: rhs(:,:), solution(:,:)
        real(dp) :: cutoff, pivot, reciprocal
        integer(int64) :: cap64, max_entries, grow
        integer :: n, width, first, j, i, idx, used, capacity, ierr
        logical :: is_diagonal

        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        n = a%n_rows
        if (a%n_cols /= n) return
        cutoff = 0.0_dp
        if (present(drop_tol)) cutoff = drop_tol
        pivot = 0.0_dp
        if (present(pivot_tol)) pivot = pivot_tol
        width = 8
        if (present(block_size)) width = block_size
        if (.not. ieee_is_finite(cutoff) .or. cutoff < 0.0_dp) return
        if (.not. ieee_is_finite(pivot) .or. pivot < 0.0_dp) return
        if (width <= 0) return
        if (n == 0) then
            out = csr_eye(0)
            if (present(info)) info = sparse_linalg_success
            return
        end if

        ! For diagonal matrices, avoid O(n) triangular solves and factor fill.
        is_diagonal = .true.
        allocate(diagonal(n))
        diagonal = 0.0_dp
        do i = 1,n
            do idx = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(idx)
                if (j == i) then
                    diagonal(i) = a%data(idx)
                else if (a%data(idx) /= 0.0_dp) then
                    is_diagonal = .false.
                end if
            end do
        end do
        if (is_diagonal) then
            if (any(abs(diagonal) <= pivot)) then
                if (present(info)) info = sparse_linalg_singular
                return
            end if
            allocate(rows(n),cols(n),values(n))
            used = 0
            do i = 1,n
                reciprocal = 1.0_dp/diagonal(i)
                if (.not. ieee_is_finite(reciprocal)) then
                    if (present(info)) info = sparse_linalg_breakdown
                    return
                end if
                if (abs(reciprocal) <= cutoff) cycle
                used = used + 1
                rows(used) = i
                cols(used) = i
                values(used) = reciprocal
            end do
            out = csr_from_coo(n,n,rows(:used),cols(:used),values(:used),ierr)
            if (ierr /= sparse_success) then
                if (present(info)) info = sparse_linalg_breakdown
            else
                if (present(info)) info = sparse_linalg_success
            end if
            return
        end if

        call splu(a,factor,ierr,pivot_tol=pivot)
        if (ierr /= sparse_linalg_success) then
            if (present(info)) info = ierr
            return
        end if
        width = min(width,n)
        allocate(rhs(n,width))
        max_entries = min(int(n,int64)*int(n,int64),int(huge(0)-1,int64))
        cap64 = min(max_entries,max(int(n,int64),2_int64*int(size(a%data),int64)))
        capacity = int(cap64)
        allocate(rows(capacity),cols(capacity),values(capacity))
        used = 0
        do first = 1,n,width
            rhs = 0.0_dp
            do j = 1,min(width,n-first+1)
                rhs(first+j-1,j) = 1.0_dp
            end do
            call splu_solve(factor,rhs(:,1:min(width,n-first+1)),solution,ierr)
            if (ierr /= sparse_linalg_success) then
                if (present(info)) info = ierr
                return
            end if
            do j = 1,size(solution,2)
                do i = 1,n
                    if (abs(solution(i,j)) <= cutoff) cycle
                    if (used == capacity) then
                        grow = min(max_entries,max(2_int64*int(capacity,int64), &
                            int(capacity,int64)+1_int64))
                        if (grow <= int(capacity,int64)) then
                            if (present(info)) info = sparse_linalg_breakdown
                            return
                        end if
                        allocate(new_rows(int(grow)),new_cols(int(grow)),new_values(int(grow)))
                        new_rows(:used) = rows(:used)
                        new_cols(:used) = cols(:used)
                        new_values(:used) = values(:used)
                        call move_alloc(new_rows,rows)
                        call move_alloc(new_cols,cols)
                        call move_alloc(new_values,values)
                        capacity = int(grow)
                    end if
                    used = used+1
                    rows(used) = i
                    cols(used) = first+j-1
                    values(used) = solution(i,j)
                end do
            end do
        end do
        out = csr_from_coo(n,n,rows(:used),cols(:used),values(:used),ierr)
        if (ierr /= sparse_success) then
            if (present(info)) info = sparse_linalg_breakdown
        else
            if (present(info)) info = sparse_linalg_success
        end if
    end function inv_csr

    function inv_csc(a, info, drop_tol, block_size, pivot_tol) result(out)
        type(csc_matrix), intent(in) :: a !! Nonsingular square CSC matrix.
        integer, intent(out), optional :: info !! 0 success; -1 invalid; -4 singular; -2 breakdown.
        real(dp), intent(in), optional :: drop_tol !! Absolute cutoff for inverse entries, default 0.
        integer, intent(in), optional :: block_size !! Identity columns per solve, default 8.
        real(dp), intent(in), optional :: pivot_tol !! Minimum pivot magnitude, default 0.
        type(csc_matrix) :: out
        type(csr_matrix) :: csr, inverse
        integer :: ierr

        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csc(a)) return
        if (a%n_rows /= a%n_cols) return
        csr = csc_to_csr(a,ierr)
        if (ierr /= sparse_success) return
        inverse = inv_csr(csr,ierr,drop_tol,block_size,pivot_tol)
        if (ierr /= sparse_linalg_success) then
            if (present(info)) info = ierr
            return
        end if
        out = csr_to_csc(inverse,ierr)
        if (ierr == sparse_success) then
            if (present(info)) info = sparse_linalg_success
        else
            if (present(info)) info = sparse_linalg_breakdown
        end if
    end function inv_csc

    ! Scalar matrix norms for scipy.sparse.linalg.norm compatibility.
    ! Numeric ord supports +/-1, +/-infinity and spectral 2; omission means
    ! Frobenius.  A character ord='fro' is accepted by a separate overload.
    function sparse_norm_csr(a, ord, info) result(value)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        real(dp), intent(in), optional :: ord !! Numeric order, +/-1, +/-infinity or 2.
        integer, intent(out), optional :: info !! 0 on success; negative on invalid/nonconvergent calculation.
        real(dp) :: value
        real(dp), allocatable :: sums(:), u(:,:), s(:), vt(:,:)
        real(dp) :: order, scale
        integer :: i, k, ierr

        value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        if (.not. present(ord)) then
            ! Scaling avoids intermediate overflow in the sum of squares.
            scale = 0.0_dp
            if (size(a%data) > 0) scale = maxval(abs(a%data))
            value = 0.0_dp
            if (scale > 0.0_dp) value = scale*sqrt(sum((a%data/scale)**2))
        else
            order = ord
            if (order == 1.0_dp .or. order == -1.0_dp) then
                if (a%n_cols == 0) return
                allocate(sums(a%n_cols))
                sums = 0.0_dp
                do k = 1,size(a%data)
                    sums(a%indices(k)) = sums(a%indices(k))+abs(a%data(k))
                end do
                if (order == 1.0_dp) then
                    value = maxval(sums)
                else
                    value = minval(sums)
                end if
            else if (order == ieee_value(0.0_dp,ieee_positive_inf) .or. &
                    order == ieee_value(0.0_dp,ieee_negative_inf)) then
                if (a%n_rows == 0) return
                allocate(sums(a%n_rows))
                sums = 0.0_dp
                do i = 1,a%n_rows
                    do k = a%indptr(i),a%indptr(i+1)-1
                        sums(i) = sums(i)+abs(a%data(k))
                    end do
                end do
                if (order > 0.0_dp) then
                    value = maxval(sums)
                else
                    value = minval(sums)
                end if
            else if (order == 2.0_dp) then
                if (min(a%n_rows,a%n_cols) == 0) return
                ! The existing Golub--Kahan solver supplies the approximate
                ! largest singular value, including true residual checks.
                call svds(a,1,u,s,vt,ierr,which='LM')
                if (ierr /= sparse_linalg_success) then
                    if (present(info)) info = ierr
                    return
                end if
                value = s(1)
            else
                return
            end if
        end if
        if (.not. ieee_is_finite(value)) then
            value = ieee_value(0.0_dp,ieee_quiet_nan)
            if (present(info)) info = sparse_linalg_breakdown
            return
        end if
        if (present(info)) info = sparse_linalg_success
    end function sparse_norm_csr

    function sparse_norm_csc(a, ord, info) result(value)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        real(dp), intent(in), optional :: ord !! Numeric order, see CSR overload.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp) :: value
        type(csr_matrix) :: csr

        value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csc(a)) return
        csr = csc_to_csr(a)
        value = sparse_norm_csr(csr,ord,info)
    end function sparse_norm_csc

    function sparse_norm_csr_int(a,ord,info) result(value)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: ord !! One of -1,1,2.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp) :: value

        value = sparse_norm_csr(a,real(ord,dp),info)
    end function sparse_norm_csr_int

    function sparse_norm_csc_int(a,ord,info) result(value)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(in) :: ord !! One of -1,1,2.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp) :: value

        value = sparse_norm_csc(a,real(ord,dp),info)
    end function sparse_norm_csc_int

    function sparse_norm_csr_named(a,ord,info) result(value)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        character(len=*), intent(in) :: ord !! Frobenius norm: 'fro'.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp) :: value

        value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(info)) info = sparse_linalg_invalid_input
        if (ord /= 'fro') return
        value = sparse_norm_csr(a,info=info)
    end function sparse_norm_csr_named

    function sparse_norm_csc_named(a,ord,info) result(value)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        character(len=*), intent(in) :: ord !! Frobenius norm: 'fro'.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp) :: value

        value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(info)) info = sparse_linalg_invalid_input
        if (ord /= 'fro') return
        value = sparse_norm_csc(a,info=info)
    end function sparse_norm_csc_named

    ! Fortran axes 1=columns (reduce rows), 2=rows (reduce columns).
    ! Unlike the scalar matrix-norm interface, ord=0 counts nonzero entries.
    function sparse_norm_axis_csr(a,axis,ord,info) result(values)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: axis !! One-based axis to reduce, 1 (rows) or 2 (columns).
        real(dp), intent(in), optional :: ord !! Vector norm order; default 2.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp), allocatable :: values(:)
        integer, allocatable :: counts(:)
        real(dp) :: order, weight
        integer :: i, k, j, count_dims, reduced_dims
        logical :: posinf, neginf

        if (present(info)) info = sparse_linalg_invalid_input
        allocate(values(0))
        if (.not. valid_csr(a) .or. (axis /= 1 .and. axis /= 2)) return
        order = 2.0_dp
        if (present(ord)) order = ord
        posinf = order == ieee_value(0.0_dp,ieee_positive_inf)
        neginf = order == ieee_value(0.0_dp,ieee_negative_inf)
        if (.not. ieee_is_finite(order) .and. .not. posinf .and. .not. neginf) return
        count_dims = a%n_cols
        reduced_dims = a%n_rows
        if (axis == 2) then
            count_dims = a%n_rows
            reduced_dims = a%n_cols
        end if
        if (reduced_dims == 0 .and. (posinf .or. neginf .or. order < 0.0_dp)) return
        deallocate(values)
        allocate(values(count_dims),counts(count_dims))
        values = 0.0_dp
        counts = 0
        if (neginf) values = huge(0.0_dp)
        do i = 1,a%n_rows
            do k = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(k)
                if (axis == 2) j = i
                weight = abs(a%data(k))
                if (weight > 0.0_dp) counts(j) = counts(j)+1
                if (order == 0.0_dp) then
                    if (weight /= 0.0_dp) values(j) = values(j)+1.0_dp
                else if (posinf) then
                    values(j) = max(values(j),weight)
                else if (neginf) then
                    values(j) = min(values(j),weight)
                else if (order == 1.0_dp) then
                    values(j) = values(j)+weight
                else if (order > 0.0_dp) then
                    values(j) = values(j)+weight**order
                else
                    if (weight > 0.0_dp) values(j) = values(j)+weight**order
                end if
            end do
        end do
        if (neginf) then
            where (counts < reduced_dims) values = 0.0_dp
        else if (order < 0.0_dp) then
            where (counts < reduced_dims) values = huge(0.0_dp)
            where (values > 0.0_dp) values = values**(1.0_dp/order)
            where (counts < reduced_dims) values = 0.0_dp
        else if (order > 0.0_dp .and. ieee_is_finite(order) .and. order /= 1.0_dp) then
            values = values**(1.0_dp/order)
        end if
        if (.not. all(ieee_is_finite(values))) then
            deallocate(values)
            allocate(values(0))
            if (present(info)) info = sparse_linalg_breakdown
            return
        end if
        if (present(info)) info = sparse_linalg_success
    end function sparse_norm_axis_csr

    function sparse_norm_axis_csc(a,axis,ord,info) result(values)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(in) :: axis !! One-based reduction axis, 1 or 2.
        real(dp), intent(in), optional :: ord !! Vector norm order.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp), allocatable :: values(:)
        type(csr_matrix) :: csr

        if (present(info)) info = sparse_linalg_invalid_input
        allocate(values(0))
        if (.not. valid_csc(a)) return
        csr = csc_to_csr(a)
        values = sparse_norm_axis_csr(csr,axis,ord,info)
    end function sparse_norm_axis_csc

    function sparse_norm_axis_csr_int(a,axis,ord,info) result(values)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(in) :: axis !! One-based reduction axis 1 or 2.
        integer, intent(in) :: ord !! Integer vector order.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp), allocatable :: values(:)

        values = sparse_norm_axis_csr(a,axis,real(ord,dp),info)
    end function sparse_norm_axis_csr_int

    function sparse_norm_axis_csc_int(a,axis,ord,info) result(values)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(in) :: axis !! One-based reduction axis 1 or 2.
        integer, intent(in) :: ord !! Integer vector order.
        integer, intent(out), optional :: info !! 0 on success.
        real(dp), allocatable :: values(:)

        values = sparse_norm_axis_csc(a,axis,real(ord,dp),info)
    end function sparse_norm_axis_csc_int

    ! Structural bandwidth; explicitly stored zeros count as entries, as
    ! in scipy.sparse.linalg.spbandwidth. Rectangular matrices are accepted.
    subroutine spbandwidth_csr(a, below, above, info)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        integer, intent(out) :: below !! Maximum offset below the main diagonal.
        integer, intent(out) :: above !! Maximum offset above the main diagonal.
        integer, intent(out), optional :: info !! 0 on success.
        integer :: i, k, j

        below = 0
        above = 0
        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csr(a)) return
        do i = 1, a%n_rows
            do k = a%indptr(i), a%indptr(i+1)-1
                j = a%indices(k)
                below = max(below,i-j)
                above = max(above,j-i)
            end do
        end do
        if (present(info)) info = sparse_linalg_success
    end subroutine spbandwidth_csr

    subroutine spbandwidth_csc(a, below, above, info)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        integer, intent(out) :: below !! Maximum offset below the main diagonal.
        integer, intent(out) :: above !! Maximum offset above the main diagonal.
        integer, intent(out), optional :: info !! 0 on success.
        integer :: i, k, j

        below = 0
        above = 0
        if (present(info)) info = sparse_linalg_invalid_input
        if (.not. valid_csc(a)) return
        do j = 1, a%n_cols
            do k = a%indptr(j), a%indptr(j+1)-1
                i = a%indices(k)
                below = max(below,i-j)
                above = max(above,j-i)
            end do
        end do
        if (present(info)) info = sparse_linalg_success
    end subroutine spbandwidth_csc

    ! Structural triangularity (not numerical triangularity). SciPy ignores
    ! explicit values and tests only whether entries lie above or below.
    subroutine is_sptriangular_csr(a, lower, upper, info)
        type(csr_matrix), intent(in) :: a !! Canonical CSR matrix.
        logical, intent(out) :: lower !! No structurally stored entry above diagonal.
        logical, intent(out) :: upper !! No structurally stored entry below diagonal.
        integer, intent(out), optional :: info !! 0 on success.
        integer :: below, above, istat

        lower = .false.
        upper = .false.
        call spbandwidth_csr(a,below,above,istat)
        if (present(info)) info = istat
        if (istat /= sparse_linalg_success) return
        lower = above == 0
        upper = below == 0
    end subroutine is_sptriangular_csr

    subroutine is_sptriangular_csc(a, lower, upper, info)
        type(csc_matrix), intent(in) :: a !! Canonical CSC matrix.
        logical, intent(out) :: lower !! No structurally stored entry above diagonal.
        logical, intent(out) :: upper !! No structurally stored entry below diagonal.
        integer, intent(out), optional :: info !! 0 on success.
        integer :: below, above, istat

        lower = .false.
        upper = .false.
        call spbandwidth_csc(a,below,above,istat)
        if (present(info)) info = istat
        if (istat /= sparse_linalg_success) return
        lower = above == 0
        upper = below == 0
    end subroutine is_sptriangular_csc

end module scifort_sparse_linalg
