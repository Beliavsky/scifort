! SPDX-License-Identifier: MIT
! Matrix input in Fortran column-major order; output complex eigenpairs.
program sparse_arnoldi_block_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigs_arnoldi
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: dense(:,:), errors(:)
    complex(dp), allocatable :: values(:), vectors(:,:)
    integer :: n, k, p, cap, width, info, cycles, locked, wanted, products, i, j, locking
    character(len=2) :: mode
    real(dp) :: tolerance, sr, si

    read(*,*) n, k, p, cap, width, mode, tolerance, sr, si, locking
    allocate(dense(n,n))
    read(*,*) dense
    a = csr_from_dense(dense)
    if (sr /= 0.0_dp .or. si /= 0.0_dp) then
        call eigs_arnoldi(a,k,values,vectors,info,which=mode,tol=tolerance, &
            maxiter=cap,ncv=p,block_size=width,schur_locking=locking/=0, &
            iterations=cycles,locked_count=wanted,schur_locked_count=locked, &
            matvecs=products,residual_norms=errors,sigma=cmplx(sr,si,dp))
    else
        call eigs_arnoldi(a,k,values,vectors,info,which=mode,tol=tolerance, &
            maxiter=cap,ncv=p,block_size=width,schur_locking=locking/=0, &
            iterations=cycles,locked_count=wanted,schur_locked_count=locked, &
            matvecs=products,residual_norms=errors)
    end if
    write(*,'(5(i0,1x))') info,cycles,wanted,locked,products
    do j=1,size(values)
        write(*,'(3(es26.17,1x))') real(values(j),dp),aimag(values(j)),errors(j)
        do i=1,size(vectors,1)
            write(*,'(2(es26.17,1x))') real(vectors(i,j),dp),aimag(vectors(i,j))
        end do
    end do
end program sparse_arnoldi_block_reference_driver
