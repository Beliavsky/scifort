! SPDX-License-Identifier: MIT
! Text-stream reference driver for CSR deletion, reshape and sparse scatter.
program sparse_deletion_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_to_coo, csr_assign, &
        csr_delete_entries, csr_delete_rows, csr_delete_columns, csr_reshape, &
        sparse_success
    implicit none
    integer :: op, nrows, ncols, nz, nr, nc, newrows, newcols, srcnz
    integer :: status, i, ios
    integer, allocatable :: r(:), c(:), selrows(:), selcols(:), sr(:), sc(:)
    real(dp), allocatable :: v(:), sv(:)
    type(csr_matrix) :: a, b, source

    do
        read(*,*,iostat=ios) op,nrows,ncols,nz,nr,nc,newrows,newcols,srcnz
        if (ios < 0) exit
        if (ios /= 0) error stop 'bad case header'
        allocate(r(nz),c(nz),v(nz),selrows(nr),selcols(nc),sr(srcnz),sc(srcnz),sv(srcnz))
        do i = 1, nz
            read(*,*) r(i),c(i),v(i)
        end do
        do i = 1, nr
            read(*,*) selrows(i)
        end do
        do i = 1, nc
            read(*,*) selcols(i)
        end do
        do i = 1, srcnz
            read(*,*) sr(i),sc(i),sv(i)
        end do
        a = csr_from_coo(nrows,ncols,r,c,v,status)
        if (status /= sparse_success) error stop 'bad initial COO'
        select case (op)
        case (1)
            call csr_delete_entries(a,selrows,selcols,status)
        case (2)
            call csr_delete_rows(a,selrows,status)
        case (3)
            call csr_delete_columns(a,selcols,status)
        case (4)
            b = csr_reshape(a,newrows,newcols,status)
            if (status == sparse_success) a = b
        case (5)
            source = csr_from_coo(nr,nc,sr,sc,sv,status)
            if (status /= sparse_success) error stop 'bad source COO'
            call csr_assign(a,selrows,selcols,source,status)
        case default
            error stop 'unknown opcode'
        end select
        call csr_to_coo(a,r,c,v,ios)
        print '(a,4(1x,i0))', 'RESULT',status,a%n_rows,a%n_cols,size(v)
        do i = 1, size(v)
            print '(i0,1x,i0,1x,es26.17)',r(i),c(i),v(i)
        end do
        deallocate(r,c,v,selrows,selcols,sr,sc,sv)
    end do
end program sparse_deletion_reference_driver
