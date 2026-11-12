! SPDX-License-Identifier: MIT
! Line-oriented cross-language driver for sparse assembly and indexing.
program sparse_assembly_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_hstack, csr_vstack, &
        csr_block_diag, csr_kron, csr_select, csr_slice, valid_csr, sparse_success
    implicit none
    type(csr_matrix), allocatable :: blocks(:)
    type(csr_matrix) :: a, b, c
    integer, allocatable :: row_ids(:), col_ids(:)
    integer :: cases, case_num, op, count, i, info, m, n
    integer :: r0, r1, c0, c1, rs, cs

    read(*,*) cases
    do case_num = 1, cases
        read(*,*) op
        info = 1
        select case (op)
        case (1,2,3)
            read(*,*) count
            allocate(blocks(count))
            do i = 1, count
                call input_csr(blocks(i))
            end do
            select case (op)
            case (1)
                c = csr_hstack(blocks,info)
            case (2)
                c = csr_vstack(blocks,info)
            case (3)
                c = csr_block_diag(blocks,info)
            end select
            deallocate(blocks)
        case (4)
            call input_csr(a)
            call input_csr(b)
            c = csr_kron(a,b,info)
        case (5)
            call input_csr(a)
            read(*,*) m,n
            allocate(row_ids(m),col_ids(n))
            do i = 1, m
                read(*,*) row_ids(i)
            end do
            do i = 1, n
                read(*,*) col_ids(i)
            end do
            c = csr_select(a,row_ids,col_ids,info)
            deallocate(row_ids,col_ids)
        case (6)
            call input_csr(a)
            read(*,*) r0,r1,c0,c1,rs,cs
            c = csr_slice(a,r0,r1,c0,c1,info,row_step=rs,col_step=cs)
        case default
            error stop 'invalid reference operation'
        end select
        if (info /= sparse_success) then
            write(*,'(i0,1x,i0,1x,i0,1x,i0)') info,0,0,0
        else
            if (.not. valid_csr(c)) error stop 'invalid output CSR'
            write(*,'(i0,1x,i0,1x,i0,1x,i0)') info,c%n_rows,c%n_cols,size(c%data)
            do i = 1, c%n_rows
                do n = c%indptr(i), c%indptr(i+1)-1
                    write(*,'(i0,1x,i0,1x,es26.17e3)') i,c%indices(n),c%data(n)
                end do
            end do
        end if
    end do
contains
    subroutine input_csr(mat)
        type(csr_matrix), intent(out) :: mat !! Parsed canonical CSR matrix.
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: nr, nc, nz, k, ierr

        read(*,*) nr,nc,nz
        allocate(rows(nz),cols(nz),vals(nz))
        do k = 1, nz
            read(*,*) rows(k),cols(k),vals(k)
        end do
        mat = csr_from_coo(nr,nc,rows,cols,vals,ierr)
        if (ierr /= sparse_success) error stop 'invalid reference source CSR'
    end subroutine input_csr
end program sparse_assembly_reference_driver
