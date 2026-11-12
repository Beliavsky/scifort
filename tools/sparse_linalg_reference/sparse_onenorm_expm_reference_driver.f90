! SPDX-License-Identifier: MIT
! Read COO test cases and return norm estimates, witnesses and exp(A).
program sparse_onenorm_expm_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_to_dense
    use scifort_sparse_linalg, only : onenormest, expm
    implicit none
    type(csr_matrix) :: a, e
    integer, allocatable :: rows(:), columns(:)
    real(dp), allocatable :: values(:), v(:), w(:), output(:,:)
    real(dp) :: estimate
    integer :: ncases, n, nnz, width, limit, chunk, info_norm, info_exp
    integer :: i, k, case_id

    read(*,*) ncases
    do case_id = 1,ncases
        read(*,*) n,nnz,width,limit,chunk
        allocate(rows(nnz),columns(nnz),values(nnz))
        do k = 1,nnz
            read(*,*) rows(k),columns(k),values(k)
        end do
        a = csr_from_coo(n,n,rows,columns,values)
        estimate = onenormest(a,t=width,itmax=limit,v=v,w=w,info=info_norm)
        e = expm(a,info_exp,block_size=chunk)
        if (info_exp == 0) then
            output = csr_to_dense(e)
        else
            allocate(output(0,0))
        end if
        write(*,'(i0,1x,i0,1x,es25.16e3)') info_norm,info_exp,estimate
        write(*,'(*(es25.16e3,1x))') v
        write(*,'(*(es25.16e3,1x))') w
        write(*,'(*(es25.16e3,1x))') reshape(output,[size(output)])
        deallocate(rows,columns,values)
    end do
end program sparse_onenorm_expm_reference_driver
