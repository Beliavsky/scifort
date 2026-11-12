! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Partial nonsymmetric complex sparse eigensolver.  An explicit, bounded
! Arnoldi/Rayleigh--Ritz search is restarted with retained Ritz directions.
! The small projected eigenproblem uses the dependency-free complex Schur
! machinery adapted from SciFort's real-matrix Arnoldi implementation.
! This is NOT ARPACK implicit restarting; shift-invert uses complex sparse LU.
module scifort_sparse_linalg_complex_eigs
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csr_complex_matvec, &
        csc_complex_to_csr, csr_complex_from_coo
    use scifort_sparse_linalg_complex_direct, only : splu_complex_factor, &
        splu_complex, splu_complex_solve, valid_splu_complex, &
        complex_direct_singular, complex_direct_success
    implicit none
    private

    integer, parameter, public :: complex_eigs_success = 0
    integer, parameter, public :: complex_eigs_no_convergence = 1
    integer, parameter, public :: complex_eigs_invalid_input = -1
    integer, parameter, public :: complex_eigs_breakdown = -2
    integer, parameter, public :: complex_eigs_singular_shift = -3
    public :: eigs_complex, eigs_complex_operator, eigs_complex_prepare_shift

    abstract interface
        subroutine complex_eigs_apply(x,y,info)
            import dp
            complex(dp), intent(in) :: x(:) !! Complex input vector, length n.
            complex(dp), intent(out) :: y(:) !! Output of A*x, length n.
            integer, intent(out) :: info !! Zero on success, nonzero for callback failure.
        end subroutine complex_eigs_apply
    end interface

    type, public :: complex_eigs_shift_cache
        private
        type(splu_complex_factor) :: factor !! Cached pivoted LU of A-sigma*I.
        type(csr_complex_matrix) :: source !! CSR snapshot to reject stale factors.
        complex(dp) :: sigma = cmplx(0.0_dp,0.0_dp,dp) !! Shift used in LU.
        logical :: ready = .false. !! True if preparation succeeded.
    end type complex_eigs_shift_cache

    interface eigs_complex_prepare_shift
        module procedure prepare_shift_csr, prepare_shift_csc
    end interface eigs_complex_prepare_shift

    interface eigs_complex
        module procedure eigs_complex_csr, eigs_complex_csc, eigs_complex_operator
    end interface eigs_complex
contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate scalar.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure real(dp) function cnorm(x) result(value)
        complex(dp), intent(in) :: x(:) !! Input vector.
        value = norm2(abs(x))
    end function cnorm

    pure real(dp) function eigen_priority(z, target) result(priority)
        complex(dp), intent(in) :: z !! Eigenvalue to rank.
        character(len=2), intent(in) :: target !! LM, SM, LR, SR, LI, or SI.
        priority = -huge(1.0_dp)
        select case (target)
        case ('LM')
            priority = abs(z)
        case ('SM')
            priority = -abs(z)
        case ('LR')
            priority = real(z,dp)
        case ('SR')
            priority = -real(z,dp)
        case ('LI')
            priority = aimag(z)
        case ('SI')
            priority = -aimag(z)
        end select
    end function eigen_priority

    subroutine reorthogonalize(x, basis, count, magnitude)
        complex(dp), intent(inout) :: x(:) !! Vector, replaced by its normalized orthogonal part.
        complex(dp), intent(in) :: basis(:,:) !! Orthonormal basis vectors (n,count).
        integer, intent(in) :: count !! Number of preceding basis vectors.
        real(dp), intent(out) :: magnitude !! Norm before normalization, zero for failure.
        integer :: i, pass
        do pass = 1,2
            do i = 1,count
                x = x-dot_product(basis(:,i),x)*basis(:,i)
            end do
        end do
        magnitude = cnorm(x)
        if (magnitude > 100.0_dp*epsilon(1.0_dp)) then
            x = x/magnitude
        else
            magnitude = 0.0_dp
        end if
    end subroutine reorthogonalize

    subroutine complete_direction(x,basis,count,seed,status)
        complex(dp), intent(inout) :: x(:) !! Candidate direction; normalized on return.
        complex(dp), intent(in) :: basis(:,:) !! Existing orthonormal columns.
        integer, intent(in) :: count !! Number of columns already accepted.
        integer, intent(in) :: seed !! Deterministic starting index for complement search.
        integer, intent(out) :: status !! Zero on success, -2 if no complement exists.
        integer :: i, n, trial
        real(dp) :: magnitude
        n = size(x)
        call reorthogonalize(x,basis,count,magnitude)
        if (magnitude /= 0.0_dp) then
            status = complex_eigs_success
            return
        end if
        ! Exact Arnoldi breakdown can indicate an invariant subspace, not an
        ! error. Expand from a deterministic orthogonal complement to allow
        ! independent eigenvectors and other spectral components to be found.
        do trial = 0,n-1
            i = mod(seed+trial-1,n)+1
            x = cmplx(0.0_dp,0.0_dp,dp)
            x(i) = cmplx(1.0_dp,0.0_dp,dp)
            call reorthogonalize(x,basis,count,magnitude)
            if (magnitude /= 0.0_dp) then
                status = complex_eigs_success
                return
            end if
        end do
        status = complex_eigs_breakdown
    end subroutine complete_direction

    subroutine factor_shifted_matrix(a,sigma,factor,info)
        type(csr_complex_matrix), intent(in) :: a !! Square, validated complex CSR matrix.
        complex(dp), intent(in) :: sigma !! Complex shift in (A-sigma*I).
        type(splu_complex_factor), intent(out) :: factor !! Reusable pivoted LU factors.
        integer, intent(out) :: info !! 0 success, -2 breakdown, -3 singular shift.
        type(csr_complex_matrix) :: shifted
        integer, allocatable :: rows(:), columns(:)
        complex(dp), allocatable :: values(:)
        integer :: n, nnz, row, j, p, status

        info = complex_eigs_breakdown
        n = a%n_rows
        nnz = size(a%data)
        if (n >= huge(0)-nnz) return
        allocate(rows(nnz+n),columns(nnz+n),values(nnz+n))
        columns(:nnz) = a%indices
        values(:nnz) = a%data
        do row = 1,n
            do p = a%indptr(row),a%indptr(row+1)-1
                rows(p) = row
            end do
        end do
        do j = 1,n
            rows(nnz+j) = j
            columns(nnz+j) = j
            values(nnz+j) = -sigma
        end do
        shifted = csr_complex_from_coo(n,n,rows,columns,values,status)
        if (status /= 0) return
        call splu_complex(shifted,factor,status)
        if (status == complex_direct_singular) then
            info = complex_eigs_singular_shift
        else if (status == complex_direct_success) then
            info = complex_eigs_success
        end if
    end subroutine factor_shifted_matrix

    subroutine eigs_complex_csc(a,k,eigenvalues,eigenvectors,info,which,tol,ncv, &
            maxiter,v0,iterations,matvecs,residual_norms,sigma,shift_cache)
        type(csc_complex_matrix), intent(in) :: a !! Canonical, square complex CSC matrix.
        integer, intent(in) :: k !! Number of requested eigenpairs, 1 <= k < n.
        complex(dp), allocatable, intent(out) :: eigenvalues(:) !! Complex eigenvalues, length k.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Right eigenvectors (n,k).
        integer, intent(out) :: info !! 0 converged, 1 exhausted, negative invalid/breakdown.
        character(len=*), intent(in), optional :: which !! LM(default),SM,LR,SR,LI,SI.
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance, >=0.
        integer, intent(in), optional :: ncv !! Search-space width, k+1 .. n.
        integer, intent(in), optional :: maxiter !! Maximum explicit restart cycles, >=1.
        complex(dp), intent(in), optional :: v0(:) !! Finite nonzero initial direction (n).
        integer, intent(out), optional :: iterations !! Number of restart cycles.
        integer, intent(out), optional :: matvecs !! Number of sparse matrix-vector products.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True residuals, length k.
        complex(dp), intent(in), optional :: sigma !! Shift; eigenvalues near sigma via LU.
        type(complex_eigs_shift_cache), intent(in), optional :: shift_cache !! Cached shifted LU.
        type(csr_complex_matrix) :: csr
        integer :: status

        info = complex_eigs_invalid_input
        if (present(iterations)) iterations = 0
        if (present(matvecs)) matvecs = 0
        allocate(eigenvalues(0),eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(eigenvalues,eigenvectors)
        if (present(residual_norms)) deallocate(residual_norms)
        call eigs_complex_csr(csr,k,eigenvalues,eigenvectors,info,which,tol, &
            ncv,maxiter,v0,iterations,matvecs,residual_norms,sigma,shift_cache)
    end subroutine eigs_complex_csc

    pure logical function matches_shift_cache(a,sigma,cache) result(ok)
        type(csr_complex_matrix), intent(in) :: a !! Matrix to match.
        complex(dp), intent(in) :: sigma !! Requested shift.
        type(complex_eigs_shift_cache), intent(in) :: cache !! Cached matrix and LU factors.
        ok = .false.
        if (.not. cache%ready .or. cache%sigma /= sigma) return
        if (.not. valid_csr_complex(cache%source)) return
        if (.not. valid_splu_complex(cache%factor)) return
        if (a%n_rows /= cache%source%n_rows .or. a%n_cols /= cache%source%n_cols) return
        if (size(a%indptr) /= size(cache%source%indptr) .or. &
            size(a%indices) /= size(cache%source%indices)) return
        if (any(a%indptr /= cache%source%indptr)) return
        if (any(a%indices /= cache%source%indices)) return
        if (any(a%data /= cache%source%data)) return
        ok = .true.
    end function matches_shift_cache

    subroutine prepare_shift_csr(a,sigma,cache,info)
        type(csr_complex_matrix), intent(in) :: a !! Canonical square complex CSR matrix.
        complex(dp), intent(in) :: sigma !! Complex spectral shift.
        type(complex_eigs_shift_cache), intent(out) :: cache !! Reusable shifted LU and matrix snapshot.
        integer, intent(out) :: info !! 0 success, -1 invalid matrix/shift, -3 singular shift.
        info = complex_eigs_invalid_input
        if (.not. valid_csr_complex(a)) return
        if (a%n_rows /= a%n_cols .or. a%n_rows < 2) return
        if (.not. finite_complex(sigma)) return
        call factor_shifted_matrix(a,sigma,cache%factor,info)
        if (info /= complex_eigs_success) return
        cache%source = a
        cache%sigma = sigma
        cache%ready = .true.
    end subroutine prepare_shift_csr

    subroutine prepare_shift_csc(a,sigma,cache,info)
        type(csc_complex_matrix), intent(in) :: a !! Canonical square complex CSC matrix.
        complex(dp), intent(in) :: sigma !! Complex spectral shift.
        type(complex_eigs_shift_cache), intent(out) :: cache !! Reusable shifted LU and matrix snapshot.
        integer, intent(out) :: info !! Preparation status.
        type(csr_complex_matrix) :: csr
        integer :: status
        info = complex_eigs_invalid_input
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        call prepare_shift_csr(csr,sigma,cache,info)
    end subroutine prepare_shift_csc

    subroutine eigs_complex_csr(a,k,eigenvalues,eigenvectors,info,which,tol,ncv, &
            maxiter,v0,iterations,matvecs,residual_norms,sigma,shift_cache)
        type(csr_complex_matrix), intent(in) :: a !! Canonical square complex CSR matrix.
        integer, intent(in) :: k !! Number of eigenpairs, 1 <= k < n.
        complex(dp), allocatable, intent(out) :: eigenvalues(:) !! Eigenvalues.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Right eigenvectors.
        integer, intent(out) :: info !! 0 success, 1 nonconvergence, negative error.
        character(len=*), intent(in), optional :: which !! LM,SM,LR,SR,LI,SI.
        real(dp), intent(in), optional :: tol !! Relative true-residual tolerance.
        integer, intent(in), optional :: ncv !! Search-space width.
        integer, intent(in), optional :: maxiter !! Maximum explicit restart cycles.
        complex(dp), intent(in), optional :: v0(:) !! Initial Arnoldi vector.
        integer, intent(out), optional :: iterations !! Number of restart cycles.
        integer, intent(out), optional :: matvecs !! Arnoldi operator applications.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! True original-system residuals.
        complex(dp), intent(in), optional :: sigma !! Shift in (A-sigma*I)^(-1).
        type(complex_eigs_shift_cache), intent(in), optional, target :: shift_cache !! Prepared reusable LU.
        type(splu_complex_factor), target :: local_factor
        type(splu_complex_factor), pointer :: factor => null()
        complex(dp), allocatable :: result(:)
        real(dp) :: anorm, score
        integer :: n, row, ierr

        info = complex_eigs_invalid_input
        if (present(iterations)) iterations = 0
        if (present(matvecs)) matvecs = 0
        allocate(eigenvalues(0),eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (.not. valid_csr_complex(a)) return
        n = a%n_rows
        if (n < 2 .or. n /= a%n_cols .or. k < 1 .or. k >= n) return
        anorm = 0.0_dp
        do row = 1,n
            score = sum(abs(a%data(a%indptr(row):a%indptr(row+1)-1)))
            if (.not. ieee_is_finite(score)) return
            anorm = max(anorm,score)
        end do
        if (.not. ieee_is_finite(anorm) .or. anorm > sqrt(huge(anorm))/4.0_dp) return
        if (present(shift_cache) .and. .not. present(sigma)) return
        if (present(sigma)) then
            if (.not. finite_complex(sigma)) return
            if (present(shift_cache)) then
                if (.not. matches_shift_cache(a,sigma,shift_cache)) return
                factor => shift_cache%factor
            else
                call factor_shifted_matrix(a,sigma,local_factor,ierr)
                if (ierr /= complex_eigs_success) then
                    info = ierr
                    return
                end if
                factor => local_factor
            end if
        end if
        deallocate(eigenvalues,eigenvectors)
        if (present(residual_norms)) deallocate(residual_norms)
        if (present(sigma)) then
            call eigs_complex_operator(n,k,apply_original,eigenvalues,eigenvectors,info, &
                which,tol,ncv,maxiter,v0,iterations,matvecs,residual_norms, &
                sigma=sigma,inverse_matvec=apply_shift,norm_estimate=anorm)
        else
            call eigs_complex_operator(n,k,apply_original,eigenvalues,eigenvectors,info, &
                which,tol,ncv,maxiter,v0,iterations,matvecs,residual_norms, &
                norm_estimate=anorm)
        end if
    contains
        subroutine apply_original(x,y,status)
            complex(dp), intent(in) :: x(:) !! Input vector.
            complex(dp), intent(out) :: y(:) !! A times input.
            integer, intent(out) :: status !! Zero on success.
            result = csr_complex_matvec(a,x,status)
            if (status == 0) then
                if (size(result) == size(y)) then
                    y = result
                else
                    status = complex_eigs_breakdown
                end if
            end if
        end subroutine apply_original

        subroutine apply_shift(x,y,status)
            complex(dp), intent(in) :: x(:) !! Input vector.
            complex(dp), intent(out) :: y(:) !! (A-sigma*I)^(-1) times input.
            integer, intent(out) :: status !! Zero on success.
            call splu_complex_solve(factor,x,result,status)
            if (status == 0) then
                if (size(result) == size(y)) then
                    y = result
                else
                    status = complex_eigs_breakdown
                end if
            end if
        end subroutine apply_shift
    end subroutine eigs_complex_csr

    ! Matrix-free complex Arnoldi.  Callback outputs must have the same size
    ! as their inputs, be finite, and represent a fixed linear operator.
    subroutine eigs_complex_operator(n,k,matvec,eigenvalues,eigenvectors,info, &
            which,tol,ncv,maxiter,v0,iterations,matvecs,residual_norms, &
            sigma,inverse_matvec,norm_estimate)
        integer, intent(in) :: n !! Dimension of the square operator, at least 2.
        integer, intent(in) :: k !! Requested eigenpairs, 1 <= k < n.
        procedure(complex_eigs_apply) :: matvec !! Applies original A (x,y,info).
        complex(dp), allocatable, intent(out) :: eigenvalues(:) !! Eigenvalues, length k on valid input.
        complex(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Right eigenvectors (n,k).
        integer, intent(out) :: info !! Convergence status; negative means invalid/breakdown.
        character(len=*), intent(in), optional :: which !! LM,SM,LR,SR,LI,SI.
        real(dp), intent(in), optional :: tol !! Scaled true-residual tolerance, >=0.
        integer, intent(in), optional :: ncv !! Arnoldi basis width, k+1..n.
        integer, intent(in), optional :: maxiter !! Maximum restart cycles.
        complex(dp), intent(in), optional :: v0(:) !! Finite nonzero initial direction.
        integer, intent(out), optional :: iterations !! Number of Arnoldi cycles.
        integer, intent(out), optional :: matvecs !! Count of Arnoldi operator calls.
        real(dp), allocatable, intent(out), optional :: residual_norms(:) !! Original-system residual norms.
        complex(dp), intent(in), optional :: sigma !! Original eigenvalues are sigma + 1/mu.
        procedure(complex_eigs_apply), optional :: inverse_matvec !! Applies (A-sigma*I)^(-1).
        real(dp), intent(in), optional :: norm_estimate !! Known estimate/bound on ||A||, default 1.
        complex(dp), allocatable :: basis(:,:), abasis(:,:), projected(:,:), &
            values(:), vectors(:,:), candidate(:,:), residual(:,:), work(:), &
            initial(:), old(:,:), start(:)
        real(dp), allocatable :: errors(:), priorities(:)
        integer, allocatable :: order(:)
        integer :: width, cap, cycle, col, row, j, i, keep, nstart, &
            status, work_status, products, tmp, best_unconverged
        real(dp) :: tolerance, anorm, magnitude
        complex(dp) :: transformed_value
        character(len=2) :: target
        logical :: converged

        info = complex_eigs_invalid_input
        if (present(iterations)) iterations = 0
        if (present(matvecs)) matvecs = 0
        allocate(eigenvalues(0),eigenvectors(0,0))
        if (present(residual_norms)) allocate(residual_norms(0))
        if (n < 2 .or. k < 1 .or. k >= n) return
        target = 'LM'
        if (present(which)) then
            if (len_trim(which) /= 2) return
            target = which
        end if
        if (target /= 'LM' .and. target /= 'SM' .and. target /= 'LR' .and. &
            target /= 'SR' .and. target /= 'LI' .and. target /= 'SI') return
        width = min(n,max(k+min(k+8,n-k),24))
        if (present(ncv)) width = ncv
        if (width <= k .or. width > n) return
        ! Unshifted small-eigenvalue Arnoldi is not generally reliable.
        if (target == 'SM' .and. width /= n) return
        if (present(sigma)) then
            if (.not. finite_complex(sigma)) return
            if (.not. present(inverse_matvec)) return
        else
            if (present(inverse_matvec)) return
        end if
        cap = 30
        if (present(maxiter)) cap = maxiter
        if (cap < 1) return
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        tolerance = max(tolerance,100.0_dp*epsilon(1.0_dp))
        if (present(v0)) then
            if (size(v0) /= n) return
            if (.not. all(finite_complex(v0))) return
            if (cnorm(v0) <= tiny(1.0_dp)) return
        end if
        ! The caller may provide a valid bound on the operator norm.
        ! This is a convergence scaling hint, not a substitute for true residuals.
        anorm = 1.0_dp
        if (present(norm_estimate)) anorm = norm_estimate
        if (.not. ieee_is_finite(anorm) .or. anorm < 0.0_dp .or. &
            anorm > sqrt(huge(anorm))/4.0_dp) return
        deallocate(eigenvalues,eigenvectors)
        allocate(eigenvalues(k),eigenvectors(n,k))
        eigenvalues = cmplx(0.0_dp,0.0_dp,dp)
        eigenvectors = cmplx(0.0_dp,0.0_dp,dp)
        if (present(residual_norms)) then
            deallocate(residual_norms)
            allocate(residual_norms(k))
            residual_norms = huge(1.0_dp)
        end if
        allocate(basis(n,width),abasis(n,width),projected(width,width), &
            values(width),vectors(width,width),candidate(n,k),residual(n,k), &
            work(n),initial(n),errors(k),priorities(width),order(width))
        basis = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1,n
            initial(i) = cmplx(sin(0.632_dp*real(i,dp))+cos(0.211_dp*real(i+3,dp)), &
                cos(0.437_dp*real(i+2,dp)),dp)
        end do
        if (present(v0)) initial = v0
        basis(:,1) = initial/cnorm(initial)
        nstart = 1
        products = 0
        info = complex_eigs_no_convergence
        do cycle = 1,cap
            ! Fresh Arnoldi expansion and an explicit Rayleigh--Ritz matrix.
            ! On later cycles the retained vectors need not form a Krylov
            ! sequence, so we evaluate every projected matrix column.
            do col = nstart+1,width
                if (present(sigma)) then
                    call inverse_matvec(basis(:,col-1),work,work_status)
                else
                    call matvec(basis(:,col-1),work,work_status)
                end if
                if (work_status /= 0) then
                    info = complex_eigs_breakdown
                    exit
                end if
                if (.not. all(finite_complex(work))) then
                    info = complex_eigs_breakdown
                    exit
                end if
                products = products+1
                call complete_direction(work,basis,col-1,col,status)
                if (status /= complex_eigs_success) then
                    info = complex_eigs_breakdown
                    exit
                end if
                basis(:,col) = work
            end do
            if (info == complex_eigs_breakdown) exit
            do col = 1,width
                if (present(sigma)) then
                    call inverse_matvec(basis(:,col),work,work_status)
                else
                    call matvec(basis(:,col),work,work_status)
                end if
                if (work_status /= 0) then
                    info = complex_eigs_breakdown
                    exit
                end if
                if (.not. all(finite_complex(work))) then
                    info = complex_eigs_breakdown
                    exit
                end if
                abasis(:,col) = work
                products = products+1
            end do
            if (info == complex_eigs_breakdown) exit
            projected = matmul(conjg(transpose(basis)),abasis)
            call schur_eigensystem(projected,values,vectors,status)
            if (status /= complex_eigs_success) then
                info = complex_eigs_breakdown
                exit
            end if
            do j = 1,width
                order(j) = j
                priorities(j) = eigen_priority(values(j),target)
            end do
            do i = 1,width-1
                row = i
                do j = i+1,width
                    if (priorities(order(j)) > priorities(order(row))) row = j
                end do
                tmp = order(i)
                order(i) = order(row)
                order(row) = tmp
            end do
            converged = .true.
            do j = 1,k
                col = order(j)
                candidate(:,j) = matmul(basis,vectors(:,col))
                magnitude = cnorm(candidate(:,j))
                if (magnitude <= tiny(magnitude) .or. .not. ieee_is_finite(magnitude)) then
                    info = complex_eigs_breakdown
                    exit
                end if
                candidate(:,j) = candidate(:,j)/magnitude
                transformed_value = values(col)
                residual(:,j) = matmul(abasis,vectors(:,col))/magnitude &
                    -transformed_value*candidate(:,j)
                if (present(sigma)) then
                    if (abs(transformed_value) <= tiny(1.0_dp)) then
                        info = complex_eigs_breakdown
                        exit
                    end if
                    eigenvalues(j) = sigma + 1.0_dp/transformed_value
                    if (.not. finite_complex(eigenvalues(j))) then
                        info = complex_eigs_breakdown
                        exit
                    end if
                    call matvec(candidate(:,j),work,work_status)
                    if (work_status /= 0) then
                        info = complex_eigs_breakdown
                        exit
                    end if
                    if (.not. all(finite_complex(work))) then
                        info = complex_eigs_breakdown
                        exit
                    end if
                    errors(j) = cnorm(work-eigenvalues(j)*candidate(:,j))
                else
                    eigenvalues(j) = transformed_value
                    errors(j) = cnorm(residual(:,j))
                end if
                if (.not. ieee_is_finite(errors(j))) then
                    info = complex_eigs_breakdown
                    exit
                end if
                eigenvectors(:,j) = candidate(:,j)
                if (errors(j) > tolerance*max(1.0_dp,anorm,abs(eigenvalues(j)))) converged = .false.
            end do
            if (info == complex_eigs_breakdown) exit
            if (present(residual_norms)) residual_norms = errors
            if (converged) then
                info = complex_eigs_success
                exit
            end if
            if (width == n .or. cycle == cap) exit
            ! Explicit thick restart: keep k desired Ritz directions and
            ! a few guards, then seed expansion with a true Ritz residual.
            keep = min(width-2,max(k+1,width/3))
            allocate(old(n,keep),start(n))
            old = matmul(basis,vectors(:,order(:keep)))
            best_unconverged = 1
            do j = 1,k
                if (errors(j) > tolerance*max(1.0_dp,anorm,abs(eigenvalues(j)))) then
                    best_unconverged = j
                    exit
                end if
            end do
            start = residual(:,best_unconverged)
            basis(:,1:keep) = cmplx(0.0_dp,0.0_dp,dp)
            do col = 1,keep
                work = old(:,col)
                call complete_direction(work,basis,col-1,col,status)
                basis(:,col) = work
            end do
            call complete_direction(start,basis,keep,keep+1,status)
            basis(:,keep+1) = start
            nstart = keep+1
            deallocate(old,start)
        end do
        if (present(iterations)) iterations = min(cycle,cap)
        if (present(matvecs)) matvecs = products
    end subroutine eigs_complex_operator

    subroutine reflector(x,v)
        complex(dp), intent(in) :: x(:) !! Column segment to eliminate.
        complex(dp), intent(out) :: v(:) !! Unit Householder vector.
        complex(dp) :: phase
        real(dp) :: magnitude, denominator
        v = cmplx(0.0_dp,0.0_dp,dp)
        magnitude = cnorm(x)
        if (magnitude <= tiny(magnitude)) return
        phase = cmplx(1.0_dp,0.0_dp,dp)
        if (abs(x(1)) > 0.0_dp) phase = x(1)/abs(x(1))
        v = x
        v(1) = v(1)+phase*magnitude
        denominator = cnorm(v)
        if (denominator > 0.0_dp) v = v/denominator
    end subroutine reflector

    subroutine complex_qr(b,u,r)
        complex(dp), intent(in) :: b(:,:) !! Square input matrix.
        complex(dp), intent(out) :: u(:,:) !! Unitary Q factor.
        complex(dp), intent(out) :: r(:,:) !! Upper triangular R factor.
        complex(dp) :: v(size(b,1)), scalar
        integer :: n, j, col, row
        n = size(b,1)
        u = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1,n
            u(j,j) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        r = b
        do j = 1,n-1
            call reflector(r(j:n,j),v(j:n))
            if (cnorm(v(j:n)) == 0.0_dp) cycle
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

    ! Reused/adapted from SciFort's own MIT-licensed dense projected Schur code.
    subroutine schur_eigensystem(h, values, vectors, info)
        complex(dp), intent(in) :: h(:,:) !! Small projected Rayleigh--Ritz matrix.
        complex(dp), intent(out) :: values(:) !! Eigenvalues, length n.
        complex(dp), intent(out) :: vectors(:,:) !! Right eigenvectors (n,n).
        integer, intent(out) :: info !! Zero for success, -2 for failure.
        complex(dp), allocatable :: t(:,:),z(:,:),u(:,:),r(:,:),block(:,:), &
            updated(:,:),col(:),v(:),y(:),tmp(:,:)
        complex(dp) :: alpha,beta,gamma,delta,mean,disc,mu,scalar
        real(dp) :: threshold,magnitude
        integer :: n,i,j,l,m,iter
        n = size(h,1)
        info = complex_eigs_breakdown
        allocate(t(n,n),z(n,n),u(n,n),r(n,n),block(n,n),updated(n,n), &
            col(n),v(n),y(n),tmp(n,n))
        t = h
        z = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1,n
            z(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        do j = 1,n-2
            call reflector(t(j+1:n,j),v(j+1:n))
            if (cnorm(v(j+1:n)) == 0.0_dp) cycle
            do l = j,n
                scalar = 2.0_dp*dot_product(v(j+1:n),t(j+1:n,l))
                t(j+1:n,l) = t(j+1:n,l)-v(j+1:n)*scalar
            end do
            do i = 1,n
                scalar = 2.0_dp*sum(t(i,j+1:n)*v(j+1:n))
                t(i,j+1:n) = t(i,j+1:n)-scalar*conjg(v(j+1:n))
                scalar = 2.0_dp*sum(z(i,j+1:n)*v(j+1:n))
                z(i,j+1:n) = z(i,j+1:n)-scalar*conjg(v(j+1:n))
            end do
            t(j+2:n,j) = cmplx(0.0_dp,0.0_dp,dp)
        end do
        m = n
        iter = 0
        do while (m > 1)
            threshold = 32.0_dp*epsilon(1.0_dp)*(abs(t(m-1,m-1))+abs(t(m,m))+1.0_dp)
            if (abs(t(m,m-1)) <= threshold) then
                t(m,m-1) = cmplx(0.0_dp,0.0_dp,dp)
                m = m-1
                cycle
            end if
            iter = iter+1
            if (iter > max(200,200*n)) return
            alpha = t(m-1,m-1)
            beta = t(m-1,m)
            gamma = t(m,m-1)
            delta = t(m,m)
            mean = 0.5_dp*(alpha+delta)
            disc = sqrt((0.5_dp*(alpha-delta))**2+beta*gamma)
            mu = mean+disc
            if (abs(mean-disc-delta) < abs(mu-delta)) mu = mean-disc
            block(:m,:m) = t(:m,:m)
            do i = 1,m
                block(i,i) = block(i,i)-mu
            end do
            call complex_qr(block(:m,:m),u(:m,:m),r(:m,:m))
            t(:m,:m) = matmul(r(:m,:m),u(:m,:m))
            do i = 1,m
                t(i,i) = t(i,i)+mu
            end do
            if (m < n) then
                updated(:m,m+1:n) = matmul(conjg(transpose(u(:m,:m))),t(:m,m+1:n))
                t(:m,m+1:n) = updated(:m,m+1:n)
            end if
            tmp(:,1:m) = matmul(z(:,1:m),u(:m,:m))
            z(:,1:m) = tmp(:,1:m)
            do j = 1,m-2
                t(j+2:m,j) = cmplx(0.0_dp,0.0_dp,dp)
            end do
            if (.not. all(finite_complex(t))) return
        end do
        do i = 1,n
            values(i) = t(i,i)
            y = cmplx(0.0_dp,0.0_dp,dp)
            y(i) = cmplx(1.0_dp,0.0_dp,dp)
            do j = i-1,1,-1
                scalar = t(j,j)-values(i)
                threshold = 100.0_dp*epsilon(1.0_dp)* &
                    max(1.0_dp,abs(values(i)),abs(t(j,j)))
                if (abs(scalar) < threshold) scalar = cmplx(threshold,0.0_dp,dp)
                y(j) = -sum(t(j,j+1:i)*y(j+1:i))/scalar
                if (abs(y(j)) > 1.0e100_dp) y(1:i) = y(1:i)*1.0e-100_dp
            end do
            col = matmul(z,y)
            magnitude = cnorm(col)
            if (magnitude <= tiny(magnitude) .or. .not. ieee_is_finite(magnitude)) return
            vectors(:,i) = col/magnitude
        end do
        info = complex_eigs_success
    end subroutine schur_eigensystem
end module scifort_sparse_linalg_complex_eigs
