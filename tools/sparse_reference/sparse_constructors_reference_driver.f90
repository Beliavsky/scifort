! SPDX-License-Identifier: MIT
! Batch reference driver for scipy.sparse diagonal and random construction.
program sparse_constructors_reference_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_seed
    use scifort_sparse, only : csr_matrix, csr_to_coo, sparse_success
    use scifort_sparse_constructors, only : sparse_diagonal, csr_diags_array, csr_random
    implicit none
    type(csr_matrix) :: a
    type(rng_state) :: rng
    type(sparse_diagonal), allocatable :: diagonals(:)
    integer, allocatable :: offsets(:), rows(:), columns(:)
    real(dp), allocatable :: values(:)
    character(len=1) :: command, dist
    integer :: ierr, stat, m, n, ndiags, t, count_t, status, i
    integer(int64) :: seed
    real(dp) :: density

    do
        read(*,*,iostat=ierr) command
        if (ierr /= 0) exit
        select case (command)
        case ('D')
            read(*,*) m,n,ndiags
            allocate(diagonals(ndiags),offsets(ndiags))
            do t = 1, ndiags
                read(*,*) offsets(t), count_t
                allocate(diagonals(t)%values(count_t))
                if (count_t > 0) read(*,*) diagonals(t)%values
            end do
            if (m == -1 .and. n == -1) then
                a = csr_diags_array(diagonals,offsets,status=status)
            else
                a = csr_diags_array(diagonals,offsets,m,n,status)
            end if
            deallocate(diagonals,offsets)
        case ('R')
            read(*,*) m,n,density,seed,dist
            call rng_seed(rng,seed)
            if (dist == 'N') then
                a = csr_random(m,n,density,rng,status,'normal')
            else
                a = csr_random(m,n,density,rng,status)
            end if
        case default
            error stop 'unrecognized constructors reference command'
        end select
        if (status /= sparse_success) then
            write(*,'(a,1x,i0)') 'ERR',status
            cycle
        end if
        call csr_to_coo(a,rows,columns,values,stat)
        if (stat /= sparse_success) error stop 'unexpected invalid CSR'
        write(*,'(a,3(1x,i0))') 'OK',a%n_rows,a%n_cols,size(values)
        do i = 1, size(values)
            write(*,'(2(i0,1x),es26.17e3)') rows(i),columns(i),values(i)
        end do
    end do
end program sparse_constructors_reference_driver
