! SPDX-License-Identifier: MIT
program sparse_lanczos_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigsh_lanczos
    implicit none
    integer :: cases, ci, n, k, p, cap, info, cycles, locked, matvecs, j
    real(dp) :: tolerance
    real(dp), allocatable :: dense(:,:), values(:), vectors(:,:), errors(:)
    type(csr_matrix) :: a
    character(len=2) :: mode

    read(*,*) cases
    do ci = 1,cases
        read(*,*) n,k,p,cap,mode,tolerance
        allocate(dense(n,n))
        read(*,*) dense
        a = csr_from_dense(dense)
        call eigsh_lanczos(a,k,values,vectors,info,which=mode, &
            tol=tolerance,maxiter=cap,ncv=p,iterations=cycles, &
            locked_count=locked,matvecs=matvecs,residual_norms=errors)
        write(*,'(a,1x,i0)') 'CASE',ci
        write(*,'(4(i0,1x))') info,cycles,locked,matvecs
        do j = 1,k
            write(*,'(*(es26.17,1x))') values(j),errors(j),vectors(:,j)
        end do
        deallocate(dense)
    end do
end program sparse_lanczos_reference_driver
