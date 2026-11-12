! SPDX-License-Identifier: MIT
! Nonsymmetric Schur invariant-subspace locking and coupled Ritz extraction.
program test_sparse_arnoldi_schur
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo
    use scifort_sparse_linalg, only : eigs, eigs_arnoldi, sparse_linalg_invalid_input
    implicit none
    integer, parameter :: n = 64
    type(csr_matrix) :: a
    real(dp) :: dense(n,n), test_norm, baseline(2)
    complex(dp), allocatable :: values(:), vectors(:,:), values2(:), vectors2(:,:)
    real(dp), allocatable :: errors(:), errors2(:)
    integer :: i, j, info, cycles, soft_count, hard_count, products

    ! Nonnormal unsymmetric problem with a nonzero upper-right coupling:
    ! a projection-only right-eigenvector deflation is not a valid substitute
    ! for retaining the triangular Schur coupling in Rayleigh--Ritz.
    dense = 0.0_dp
    do i = 1,n
        dense(i,i) = 0.6_dp*cos(0.24_dp*real(i,dp))
        if (i < n) then
            dense(i,i+1) = 0.07_dp
            dense(i+1,i) = -0.025_dp
        end if
    end do
    dense(1,1) = 9.0_dp
    dense(2,2) = 6.0_dp
    dense(1,3) = 5.0_dp
    dense(2,4) = -3.0_dp
    dense(1,6) = 7.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LR',ncv=6, &
        maxiter=60,tol=1.0e-8_dp,schur_locking=.true., &
        schur_locked_count=hard_count,locked_count=soft_count, &
        iterations=cycles,matvecs=products,residual_norms=errors)
    call require(info == 0, 'nonnormal Schur convergence')
    call require(hard_count >= 1 .and. hard_count <= 2, 'hard-locking occurred')
    call require(soft_count == 2 .and. cycles > 1, 'Schur restart diagnostics')
    call require(products >= 6, 'operator products counted')
    call check_eigenpairs(dense,values,vectors,errors,1.0e-6_dp)
    call require(abs(real(values(1),dp)-9.0_dp) < 0.02_dp, 'dominant nonnormal eigenvalue')
    call require(abs(real(values(2),dp)-6.0_dp) < 0.02_dp, 'secondary nonnormal eigenvalue')
    baseline = real(values,dp)

    ! Method dispatch matches the direct opt-in interface.
    call eigs(a,2,values2,vectors2,info,which='LR',ncv=6,maxiter=60, &
        tol=1.0e-8_dp,method='arnoldi_schur',residual_norms=errors2)
    call require(info == 0, 'Schur dispatch convergence')
    call check_eigenpairs(dense,values2,vectors2,errors2,1.0e-6_dp)
    call require(maxval(abs(real(values2,dp)-baseline)) < 1.0e-6_dp, 'dispatch agreement')

    ! A complex spectrum that is not a symmetric eigenproblem.
    dense(1:4,1:4) = 0.0_dp
    dense(1,2) = -9.0_dp
    dense(2,1) = 9.0_dp
    dense(3,4) = -5.0_dp
    dense(4,3) = 5.0_dp
    dense(1,6) = 3.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LI',tol=1.0e-8_dp, &
        ncv=6,maxiter=70,schur_locking=.true.,schur_locked_count=hard_count, &
        residual_norms=errors)
    call require(info == 0, 'complex Schur convergence')
    call require(hard_count >= 1, 'complex Schur hard lock')
    call check_eigenpairs(dense,values,vectors,errors,1.0e-6_dp)
    call require(abs(aimag(values(1))-9.0_dp) < 0.02_dp, 'largest imaginary eigenvalue')
    call require(abs(aimag(values(2))-5.0_dp) < 0.02_dp, 'second imaginary eigenvalue')

    ! Three well-separated desired eigenpairs with a narrow search space.
    ! Two vectors become hard-locked while the third continues iterating:
    ! this exercises a nontrivial 2x2 complex Schur rotation and the upper-
    ! triangular coupling to the remaining Krylov search directions.
    dense = 0.0_dp
    do i = 1,n
        dense(i,i) = 0.6_dp*cos(0.24_dp*real(i,dp))
        if (i < n) then
            dense(i,i+1) = 0.07_dp
            dense(i+1,i) = -0.025_dp
        end if
    end do
    dense(1,1) = 9.0_dp
    dense(2,2) = 7.0_dp
    dense(3,3) = 5.0_dp
    dense(1,6) = 1.1_dp
    dense(2,7) = -0.8_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,3,values,vectors,info,which='LR',ncv=5, &
        maxiter=75,tol=1.0e-8_dp,schur_locking=.true., &
        schur_locked_count=hard_count,iterations=cycles,residual_norms=errors)
    call require(info == 0, 'multiple Schur-lock convergence')
    call require(hard_count >= 2 .and. cycles >= 3, 'Schur subspace dimension two')
    call check_eigenpairs(dense,values,vectors,errors,1.0e-6_dp)
    call require(abs(real(values(1),dp)-9.0_dp) < 0.02_dp .and. &
        abs(real(values(2),dp)-7.0_dp) < 0.02_dp .and. &
        abs(real(values(3),dp)-5.0_dp) < 0.02_dp, 'three locked-spectrum eigenvalues')

    ! Complex shift-invert with an interior eigenvalue.  The locked space is
    ! invariant under the inverse operator as well as under the original A.
    dense = 0.0_dp
    do i = 1,n
        dense(i,i) = 0.6_dp*cos(0.24_dp*real(i,dp))
        if (i < n) then
            dense(i,i+1) = 0.07_dp
            dense(i+1,i) = -0.025_dp
        end if
    end do
    dense(1,1) = 1.0_dp
    dense(2,2) = 3.0_dp
    dense(3,3) = 6.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LM',tol=1.0e-8_dp, &
        ncv=10,maxiter=80,sigma=cmplx(0.98_dp,0.0_dp,dp), &
        schur_locking=.true.,schur_locked_count=hard_count,residual_norms=errors)
    call require(info == 0, 'shifted Schur convergence')
    call require(hard_count >= 1, 'shifted Schur locking')
    call check_eigenpairs(dense,values,vectors,errors,1.0e-6_dp)
    call require(abs(real(values(1),dp)-1.0_dp) < 0.02_dp, 'shifted interior eigenvalue')

    call eigs(a,2,values2,vectors2,info,which='LM',ncv=10,maxiter=80, &
        tol=1.0e-8_dp,sigma=cmplx(0.98_dp,0.0_dp,dp), &
        method='arnoldi_schur',residual_norms=errors2)
    call require(info == 0, 'shifted Schur method dispatch')
    call check_eigenpairs(dense,values2,vectors2,errors2,1.0e-6_dp)
    call require(maxval(abs(values2-values)) < 1.0e-5_dp, 'shifted dispatch agreement')

    ! Sparse large-n regression: keep only O(n*ncv) eigensolver storage.
    block
        integer, parameter :: ns = 900
        integer :: rows(3*ns-2), cols(3*ns-2), t
        real(dp) :: entries(3*ns-2)
        do t = 1,ns
            rows(t) = t
            cols(t) = t
            entries(t) = 0.3_dp*cos(real(t,dp)*0.12_dp)
        end do
        entries(1) = 9.0_dp
        entries(2) = 7.0_dp
        do t = 1,ns-1
            rows(ns+t) = t
            cols(ns+t) = t+1
            entries(ns+t) = 0.03_dp
            rows(2*ns-1+t) = t+1
            cols(2*ns-1+t) = t
            entries(2*ns-1+t) = -0.018_dp
        end do
        a = csr_from_coo(ns,ns,rows,cols,entries)
        call eigs_arnoldi(a,2,values,vectors,info,which='LR',ncv=8, &
            maxiter=60,tol=1.0e-8_dp,schur_locking=.true., &
            schur_locked_count=hard_count,residual_norms=errors)
        call require(info == 0, '900-row CSR Schur convergence')
        call require(maxval(errors) < 1.0e-6_dp, 'large sparse residual')
    end block

    ! A one-cycle budget must not silently certify unconverged eigenpairs.
    dense(1,1) = 9.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LR',ncv=6, &
        maxiter=1,tol=1.0e-12_dp,schur_locking=.true., &
        schur_locked_count=hard_count,iterations=cycles)
    call require(info > 0 .and. cycles == 1, 'restart budget exhaustion')
    call require(hard_count >= 0, 'nonconvergence locking diagnostic')

    call eigs_arnoldi(a,2,values,vectors,info,ncv=2,schur_locking=.true., &
        schur_locked_count=hard_count)
    call require(info == sparse_linalg_invalid_input .and. hard_count == 0, 'invalid subspace')
    call eigs_arnoldi(a,2,values,vectors,info,schur_locking=.true., &
        v0=[(cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp),i=1,n)], &
        schur_locked_count=hard_count)
    call require(info == sparse_linalg_invalid_input .and. hard_count == 0, 'nonfinite guess')
    call eigs(a,2,values,vectors,info,method='invalid')
    call require(info == sparse_linalg_invalid_input, 'invalid dispatch')
    print '(a)', 'test_sparse_arnoldi_schur: PASS'
contains
    subroutine check_eigenpairs(d, eigvals, eigvecs, residuals, limit)
        real(dp), intent(in) :: d(:,:) !! Original dense matrix for independent validation.
        complex(dp), intent(in) :: eigvals(:) !! Returned complex eigenvalues.
        complex(dp), intent(in) :: eigvecs(:,:) !! Returned right eigenvectors.
        real(dp), intent(in) :: residuals(:) !! Reported true residuals.
        real(dp), intent(in) :: limit !! Maximum allowed true residual.
        integer :: col
        real(dp) :: actual
        call require(size(eigvals) == size(residuals), 'residual shape')
        do col = 1,size(eigvals)
            actual = sqrt(sum(abs(matmul(d,eigvecs(:,col))-eigvals(col)*eigvecs(:,col))**2))
            call require(actual < limit, 'independent true eigenpair residual')
            call require(abs(actual-residuals(col)) < 1.0e-10_dp, 'residual diagnostics')
            call require(abs(sqrt(sum(abs(eigvecs(:,col))**2))-1.0_dp) < 1.0e-10_dp, &
                'unit eigenvectors')
        end do
    end subroutine check_eigenpairs

    subroutine require(condition,message)
        logical, intent(in) :: condition !! Test assertion.
        character(len=*), intent(in) :: message !! Test description.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_arnoldi_schur
