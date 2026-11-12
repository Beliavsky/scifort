! SPDX-License-Identifier: MIT
! Line-oriented driver for Python/SciPy block/assignment reference comparisons.
program sparse_block_assignment_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_bmat, csr_assign, &
        csr_assign_coo, csr_set, valid_csr, sparse_success
    implicit none
    type(csr_matrix) :: a, c, source
    type(csr_matrix), allocatable :: blocks(:,:)
    integer, allocatable :: heights(:), widths(:), rr(:), cc(:)
    real(dp), allocatable :: vv(:), dense(:,:)
    real(dp) :: val
    integer :: cases, case_no, op, nr, nc, m, k, i, j, r, col, info, flag

    read(*,*) cases
    do case_no = 1, cases
        read(*,*) op
        info = 1
        select case(op)
        case(1)
            read(*,*) nr,nc
            allocate(blocks(nr,nc),heights(nr),widths(nc))
            do i = 1, nr
                read(*,*) heights(i)
            end do
            do j = 1, nc
                read(*,*) widths(j)
            end do
            do j = 1, nc
                do i = 1, nr
                    read(*,*) flag
                    if (flag == 1) call input_csr(blocks(i,j))
                end do
            end do
            c = csr_bmat(blocks,info,row_sizes=heights,col_sizes=widths)
            deallocate(blocks,heights,widths)
        case(2)
            call input_csr(a)
            read(*,*) m
            allocate(rr(m),cc(m),vv(m))
            do k = 1, m
                read(*,*) rr(k),cc(k),vv(k)
            end do
            call csr_assign_coo(a,rr,cc,vv,info)
            c = a
            deallocate(rr,cc,vv)
        case(3,4,5)
            call input_csr(a)
            read(*,*) nr,nc
            allocate(rr(nr),cc(nc))
            do i = 1, nr
                read(*,*) rr(i)
            end do
            do j = 1, nc
                read(*,*) cc(j)
            end do
            select case(op)
            case(3)
                allocate(dense(nr,nc))
                do j = 1, nc
                    do i = 1, nr
                        read(*,*) dense(i,j)
                    end do
                end do
                call csr_assign(a,rr,cc,dense,info)
                deallocate(dense)
            case(4)
                read(*,*) val
                call csr_assign(a,rr,cc,val,info)
            case(5)
                call input_csr(source)
                call csr_assign(a,rr,cc,source,info)
            end select
            c = a
            deallocate(rr,cc)
        case(6)
            call input_csr(a)
            read(*,*) r,col,val
            call csr_set(a,r,col,val,info)
            c = a
        case default
            error stop 'unrecognized operation'
        end select
        if (info /= sparse_success) then
            write(*,'(4(i0,1x))') info,0,0,0
        else
            if (.not. valid_csr(c)) error stop 'noncanonical output'
            write(*,'(4(i0,1x))') info,c%n_rows,c%n_cols,size(c%data)
            do i = 1, c%n_rows
                do k = c%indptr(i),c%indptr(i+1)-1
                    write(*,'(i0,1x,i0,1x,es26.17e3)') i,c%indices(k),c%data(k)
                end do
            end do
        end if
    end do
contains
    subroutine input_csr(mat)
        type(csr_matrix), intent(out) :: mat !! New canonical matrix.
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: height, width, nz, k, ierr
        read(*,*) height,width,nz
        allocate(rows(nz),cols(nz),vals(nz))
        do k = 1,nz
            read(*,*) rows(k),cols(k),vals(k)
        end do
        mat = csr_from_coo(height,width,rows,cols,vals,ierr)
        if (ierr /= sparse_success) error stop 'source input invalid'
    end subroutine input_csr
end program sparse_block_assignment_reference_driver
