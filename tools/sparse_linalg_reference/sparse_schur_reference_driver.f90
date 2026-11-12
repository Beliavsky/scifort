! SPDX-License-Identifier: MIT
! Column-major dense input; Schur-locked Arnoldi outputs for independent tests.
program sparse_schur_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigs_arnoldi
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: dense(:,:), residuals(:)
    complex(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
    integer :: n, k, p, cap, info, cycles, wanted, frozen, products, j, i
    character(len=2) :: mode
    real(dp) :: tolerance, sigma_real, sigma_imag

    read(*,*) n,k,p,cap,mode,tolerance,sigma_real,sigma_imag
    allocate(dense(n,n))
    read(*,*) dense
    a = csr_from_dense(dense)
    if (sigma_real /= 0.0_dp .or. sigma_imag /= 0.0_dp) then
        call eigs_arnoldi(a,k,eigenvalues,eigenvectors,info,which=mode, &
            tol=tolerance,ncv=p,maxiter=cap, &
            schur_locking=.true.,schur_locked_count=frozen, &
            locked_count=wanted,matvecs=products,iterations=cycles, &
            residual_norms=residuals, &
            sigma=cmplx(sigma_real,sigma_imag,dp))
    else
        call eigs_arnoldi(a,k,eigenvalues,eigenvectors,info,which=mode, &
            tol=tolerance,ncv=p,maxiter=cap, &
            schur_locking=.true.,schur_locked_count=frozen, &
            locked_count=wanted,matvecs=products,iterations=cycles, &
            residual_norms=residuals)
    end if
    write(*,'(5(i0,1x))') info,cycles,wanted,frozen,products
    do j=1,size(eigenvalues)
        write(*,'(3(es26.17,1x))') real(eigenvalues(j),dp),aimag(eigenvalues(j)),residuals(j)
        do i=1,size(eigenvectors,1)
            write(*,'(2(es26.17,1x))') real(eigenvectors(i,j),dp),aimag(eigenvectors(i,j))
        end do
    end do
end program sparse_schur_reference_driver
