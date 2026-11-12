! SPDX-License-Identifier: MIT
program sparse_lanczos_locked_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense
    use scifort_sparse_linalg, only : eigsh_lanczos
    implicit none
    integer :: cases, ci, n, k, p, cap, info, cycles, locked, matvecs, j, block_size, shifted
    real(dp) :: tolerance, sigma
    real(dp), allocatable :: dense(:,:), values(:), vectors(:,:), errors(:)
    type(csr_matrix) :: a
    character(len=2) :: mode

    read(*,*) cases
    do ci = 1,cases
        read(*,*) n,k,p,cap,mode,tolerance,block_size,shifted,sigma
        allocate(dense(n,n))
        read(*,*) dense
        a = csr_from_dense(dense)
        if (shifted == 1) then
            call eigsh_lanczos(a,k,values,vectors,info,which=mode, &
                tol=tolerance,maxiter=cap,ncv=p,iterations=cycles, &
                locked_count=locked,matvecs=matvecs,residual_norms=errors, &
                hard_locking=.true.,block_size=block_size,sigma=sigma)
        else
            call eigsh_lanczos(a,k,values,vectors,info,which=mode, &
                tol=tolerance,maxiter=cap,ncv=p,iterations=cycles, &
                locked_count=locked,matvecs=matvecs,residual_norms=errors, &
                hard_locking=.true.,block_size=block_size)
        end if
        write(*,'(a,1x,i0)') 'CASE',ci
        write(*,'(4(i0,1x))') info,cycles,locked,matvecs
        do j = 1,k
            write(*,'(*(es26.17,1x))') values(j),errors(j),vectors(:,j)
        end do
        deallocate(dense)
    end do
end program sparse_lanczos_locked_reference_driver
