! SPDX-License-Identifier: MIT
! Machine-readable reference driver for the action of the CSR exponential.
program sparse_expm_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo
    use scifort_sparse_linalg, only : expm_multiply, expm_multiply_interval
    implicit none
    integer :: count_cases, ci, n, nrhs, nnz, num, j, k, i, ip, info, info_grid
    integer :: include_last, have_trace, status, nprod, nprod_grid
    integer, allocatable :: rows(:), cols(:)
    real(dp), allocatable :: values(:), b(:,:), y(:,:), grid(:,:,:)
    real(dp) :: time, start, stop, trace_value
    type(csr_matrix) :: a

    read(*,*) count_cases
    do ci = 1,count_cases
        read(*,*) n,nrhs,nnz,time,start,stop,num,include_last,have_trace,trace_value
        allocate(rows(nnz),cols(nnz),values(nnz),b(n,nrhs))
        do k = 1,nnz
            read(*,*) rows(k),cols(k),values(k)
        end do
        do j = 1,nrhs
            do i = 1,n
                read(*,*) b(i,j)
            end do
        end do
        a = csr_from_coo(n,n,rows,cols,values,status)
        if (status /= 0) error stop 'invalid COO in exponential reference input'
        if (have_trace == 1) then
            call expm_multiply(a,b,y,info,t=time,tracea=trace_value,products=nprod)
            call expm_multiply_interval(a,b,start,stop,num,grid,info_grid, &
                endpoint=include_last == 1,tracea=trace_value,products=nprod_grid)
        else
            call expm_multiply(a,b,y,info,t=time,products=nprod)
            call expm_multiply_interval(a,b,start,stop,num,grid,info_grid, &
                endpoint=include_last == 1,products=nprod_grid)
        end if
        write(*,'(a,1x,i0)') 'CASE',ci
        write(*,'(4(i0,1x))') info,info_grid,nprod,nprod_grid
        if (info /= 0 .or. info_grid /= 0) then
            print '(a)', 'FAILED'
        else
            do j = 1,nrhs
                write(*,'(*(es30.17e3,1x))') y(:,j)
            end do
            do ip = 1,num
                do j = 1,nrhs
                    write(*,'(*(es30.17e3,1x))') grid(:,j,ip)
                end do
            end do
        end if
        deallocate(rows,cols,values,b)
    end do
end program sparse_expm_reference_driver
