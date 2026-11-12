! SPDX-License-Identifier: MIT
! Block restarted Arnoldi: exact and nearly repeated nonsymmetric spectra.
program test_sparse_arnoldi_block
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo
    use scifort_sparse_linalg, only : eigs, eigs_arnoldi, sparse_linalg_invalid_input
    implicit none
    integer, parameter :: n = 48
    type(csr_matrix) :: a
    real(dp) :: dense(n,n)
    complex(dp), allocatable :: values(:), vectors(:,:), dispatch_values(:), dispatch_vectors(:,:)
    real(dp), allocatable :: errors(:), dispatch_errors(:)
    integer :: i, info, cycles, nlocked

    call background_matrix(dense)
    dense(1,1) = 9.0_dp
    dense(2,2) = 9.0_dp
    dense(3,3) = 8.9_dp
    dense(1,2) = 0.0_dp  ! Semisimple double eigenvalue, not a Jordan block.
    dense(1,5) = 0.15_dp
    dense(2,7) = -0.25_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LR',ncv=10, &
        maxiter=90,tol=1.0e-8_dp,block_size=2,schur_locking=.true., &
        residual_norms=errors,iterations=cycles,schur_locked_count=nlocked)
    call require(info == 0 .and. cycles > 0, 'double eigenvalue convergence')
    call check_eigenpairs(dense,values,vectors,errors,2.0e-7_dp)
    call require(maxval(abs(values-cmplx(9.0_dp,0.0_dp,dp))) < 2.0e-4_dp, &
        'both repeated eigenvalues recovered')
    call independent_pair(vectors(:,1),vectors(:,2))

    ! The public dispatcher must forward block_size to either Arnoldi variant.
    call eigs(a,2,dispatch_values,dispatch_vectors,info,which='LR',ncv=10, &
        maxiter=90,tol=1.0e-8_dp,method='arnoldi_schur',block_size=2, &
        residual_norms=dispatch_errors)
    call require(info == 0, 'Schur block dispatch')
    call check_eigenpairs(dense,dispatch_values,dispatch_vectors,dispatch_errors,2.0e-7_dp)
    call require(maxval(abs(dispatch_values-values)) < 2.0e-5_dp, 'direct and dispatch agree')
    call eigs(a,2,dispatch_values,dispatch_vectors,info,which='LR',ncv=10, &
        maxiter=90,tol=1.0e-8_dp,method='arnoldi',block_size=2)
    call require(info == 0, 'soft block dispatch')
    call require(maxval(abs(dispatch_values-cmplx(9.0_dp,0.0_dp,dp))) < 2.0e-4_dp, &
        'soft block repeats')

    ! A cluster of three nearby eigenvalues with nonnormal right coupling.
    call background_matrix(dense)
    dense(1,1) = 9.0_dp
    dense(2,2) = 8.999_dp
    dense(3,3) = 8.998_dp
    dense(1,2) = 0.0_dp
    dense(2,3) = 0.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,3,values,vectors,info,which='LR',ncv=12, &
        maxiter=80,tol=1.0e-9_dp,block_size=3,schur_locking=.true.,residual_norms=errors)
    call require(info == 0, 'clustered three eigenvalues')
    call check_eigenpairs(dense,values,vectors,errors,2.0e-7_dp)
    call require(maxval(abs(real(values,dp)-[9.0_dp,8.999_dp,8.998_dp])) < 2.0e-4_dp, &
        'clustered real eigenvalues match')

    ! Two independent real 2x2 rotation blocks give a repeated complex
    ! eigenvalue +8i.  Treat the two linearly independent eigenvectors.
    call background_matrix(dense)
    dense(1:4,1:4) = 0.0_dp
    dense(1,2) = -8.0_dp
    dense(2,1) = 8.0_dp
    dense(3,4) = -8.0_dp
    dense(4,3) = 8.0_dp
    dense(1,7) = 0.5_dp
    dense(3,6) = 0.25_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LI',ncv=12, &
        maxiter=80,tol=1.0e-9_dp,block_size=2,schur_locking=.true.,residual_norms=errors)
    call require(info == 0, 'repeated complex eigenvalues')
    call check_eigenpairs(dense,values,vectors,errors,2.0e-7_dp)
    call require(maxval(abs(values-cmplx(0.0_dp,8.0_dp,dp))) < 2.0e-4_dp, &
        'both imaginary eigenvalues recovered')
    call independent_pair(vectors(:,1),vectors(:,2))

    ! Interior double eigenvalue: complex inverse solves must still verify
    ! residuals against the original, unshifted real matrix.
    call background_matrix(dense)
    dense(1,1) = 1.0_dp
    dense(2,2) = 1.0_dp
    dense(3,3) = 1.25_dp
    dense(1,2) = 0.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,which='LM',ncv=12, &
        maxiter=70,tol=1.0e-9_dp,sigma=cmplx(1.025_dp,0.0_dp,dp), &
        block_size=2,schur_locking=.true.,residual_norms=errors)
    call require(info == 0, 'shift-invert repeated eigenvalues')
    call check_eigenpairs(dense,values,vectors,errors,2.0e-7_dp)
    call require(maxval(abs(values-cmplx(1.0_dp,0.0_dp,dp))) < 2.0e-4_dp, &
        'both shifted eigenvalues recovered')
    call independent_pair(vectors(:,1),vectors(:,2))

    ! A partly converged nonnormal problem exercises frozen Schur coupling
    ! and multiple residual directions on a subsequent restart.
    call background_matrix(dense)
    dense(1,1) = 9.0_dp
    dense(2,2) = 7.0_dp
    dense(3,3) = 5.0_dp
    dense(1,4) = 1.2_dp
    dense(2,5) = -0.8_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,3,values,vectors,info,which='LR',ncv=9, &
        maxiter=80,tol=1.0e-9_dp,block_size=2,schur_locking=.true., &
        residual_norms=errors,iterations=cycles,schur_locked_count=nlocked)
    call require(info == 0 .and. cycles >= 2 .and. nlocked >= 1, 'block Schur lock')
    call check_eigenpairs(dense,values,vectors,errors,2.0e-7_dp)

    ! Explicit input validation.  Subspace iteration does not accept this
    ! distinct Arnoldi control, avoiding silent ignored arguments.
    call eigs_arnoldi(a,2,values,vectors,info,ncv=9,block_size=0)
    call require(info == sparse_linalg_invalid_input, 'zero block width')
    call eigs_arnoldi(a,2,values,vectors,info,ncv=9,block_size=3)
    call require(info == sparse_linalg_invalid_input, 'block wider than k')
    call eigs_arnoldi(a,2,values,vectors,info,ncv=4,block_size=2)
    call require(info == sparse_linalg_invalid_input, 'insufficient restart slots')
    call eigs(a,2,values,vectors,info,ncv=9,block_size=2)
    call require(info == sparse_linalg_invalid_input, 'subspace rejects block_size')
    call eigs_arnoldi(a,2,values,vectors,info,ncv=9,block_size=2,maxiter=0)
    call require(info > 0, 'zero budget is nonconvergence')

    ! Full search spaces are accepted without attempting an impossible
    ! additional restart when ncv == n.
    block
        real(dp) :: small(5,5)
        integer :: j
        small = 0.0_dp
        do j=1,5
            small(j,j)=real(j,dp)
        end do
        small(1,1)=5.0_dp
        small(2,2)=5.0_dp
        a = csr_from_dense(small)
        call eigs_arnoldi(a,2,values,vectors,info,ncv=5,block_size=2, &
            maxiter=3,tol=1.0e-8_dp)
        call require(info == 0, 'full block subspace')
        call check_eigenpairs(small,values,vectors,[0.0_dp,0.0_dp],1.0e-7_dp, .false.)
    end block
    ! Large sparse semisimple eigenspace: exercise O(n*ncv), not dense n*n,
    ! storage for the search itself.  The independent dense references above
    ! provide correctness checks for smaller nonnormal and shifted matrices.
    block
        integer, parameter :: ns = 1400
        integer :: rows(2*ns-1), cols(2*ns-1), j
        real(dp) :: entries(2*ns-1)
        do j = 1,ns
            rows(j) = j
            cols(j) = j
            entries(j) = 0.3_dp*cos(real(j,dp)*0.12_dp)
            if (j < ns) then
                rows(ns+j) = j
                cols(ns+j) = j+1
                entries(ns+j) = 0.04_dp
            end if
        end do
        entries(1) = 9.0_dp
        entries(2) = 9.0_dp
        entries(ns+1) = 0.0_dp
        a = csr_from_coo(ns,ns,rows,cols,entries)
        call eigs_arnoldi(a,2,values,vectors,info,which='LR',block_size=2, &
            ncv=10,tol=1.0e-8_dp,maxiter=60,schur_locking=.true., &
            residual_norms=errors)
        call require(info == 0, '1400-row sparse repeated eigenspace')
        call require(maxval(abs(values-cmplx(9.0_dp,0.0_dp,dp))) < 1.0e-4_dp, &
            'large repeated eigenvalues')
        call require(maxval(errors) < 2.0e-7_dp, 'large sparse true residuals')
        call independent_pair(vectors(:,1),vectors(:,2))
    end block
    print '(a)', 'test_sparse_arnoldi_block: PASS'
contains
    subroutine background_matrix(d)
        real(dp), intent(out) :: d(:,:) !! Nonsymmetric upper-bidiagonal baseline.
        integer :: j
        d = 0.0_dp
        do j=1,size(d,1)
            d(j,j)=0.2_dp*cos(real(j,dp)*0.15_dp)
            if (j < size(d,1)) d(j,j+1)=0.02_dp
        end do
    end subroutine background_matrix

    subroutine check_eigenpairs(d, eigvals, eigvecs, residuals, limit, check_report)
        real(dp), intent(in) :: d(:,:) !! Original dense matrix.
        complex(dp), intent(in) :: eigvals(:) !! Eigenvalues returned by SciFort.
        complex(dp), intent(in) :: eigvecs(:,:) !! Normalized eigenvectors.
        real(dp), intent(in) :: residuals(:) !! Reported residuals.
        real(dp), intent(in) :: limit !! Maximum independent residual.
        logical, intent(in), optional :: check_report !! Compare reported and actual norms.
        logical :: compare
        real(dp) :: actual
        integer :: j
        compare=.true.
        if (present(check_report)) compare=check_report
        call require(size(eigvals) == size(eigvecs,2), 'eigenvector shape')
        do j=1,size(eigvals)
            actual = sqrt(sum(abs(matmul(d,eigvecs(:,j))-eigvals(j)*eigvecs(:,j))**2))
            call require(actual < limit, 'independently recomputed residual')
            if (compare) call require(abs(actual-residuals(j)) < 1.0e-10_dp, 'reported residual')
            call require(abs(sqrt(sum(abs(eigvecs(:,j))**2))-1.0_dp) < 1.0e-10_dp, &
                'unit eigenvector')
        end do
    end subroutine check_eigenpairs

    subroutine independent_pair(v1,v2)
        complex(dp), intent(in) :: v1(:) !! First normalized right eigenvector.
        complex(dp), intent(in) :: v2(:) !! Second normalized right eigenvector.
        ! Do not require orthogonality: nonsymmetric matrices need not have
        ! mutually orthogonal right eigenvectors.
        call require(1.0_dp-abs(dot_product(v1,v2))**2 > 1.0e-2_dp, &
            'two independent eigenvectors')
    end subroutine independent_pair

    subroutine require(ok, message)
        logical, intent(in) :: ok !! Condition under test.
        character(len=*), intent(in) :: message !! Failure description.
        if (.not. ok) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_arnoldi_block
