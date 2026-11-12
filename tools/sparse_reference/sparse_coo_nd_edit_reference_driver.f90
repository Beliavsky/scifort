! SPDX-License-Identifier: MIT
! Streaming reference driver for COO N-D boolean indexing and scatter edits.
program sparse_coo_nd_edit_reference_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, coo_nd_to_dense_flat
    use scifort_sparse_coo_nd_edit, only : coo_nd_mask_select, coo_nd_mask_axis, &
        coo_nd_scatter, coo_nd_mask_assign
    implicit none
    type(coo_nd_array) :: a, b
    logical, allocatable :: mask(:)
    integer, allocatable :: points(:,:)
    real(dp), allocatable :: values(:), flat(:)
    integer :: cases, id, mode, nmask, ax, nupdates, add_flag, zeros_flag, limit, &
        status, k
    real(dp) :: scalar

    read(*,*) cases
    do id = 1, cases
        read(*,*) mode
        call read_array(a)
        select case (mode)
        case (1)
            read(*,*) nmask
            allocate(mask(nmask))
            if (nmask > 0) read(*,*) mask
            b = coo_nd_mask_select(a,mask,status)
            deallocate(mask)
        case (2)
            read(*,*) ax,nmask
            allocate(mask(nmask))
            if (nmask > 0) read(*,*) mask
            b = coo_nd_mask_axis(a,ax,mask,status)
            deallocate(mask)
        case (3,4)
            read(*,*) nupdates,zeros_flag,limit
            allocate(points(size(a%shape),nupdates),values(nupdates))
            do k = 1, nupdates
                read(*,*) points(:,k),values(k)
            end do
            call coo_nd_scatter(a,points,values,status,accumulate=(mode==4), &
                keep_zeros=(zeros_flag/=0),max_nnz=int(limit,int64))
            b = a
            deallocate(points,values)
        case (5)
            read(*,*) nmask,add_flag,zeros_flag,scalar,limit
            allocate(mask(nmask))
            if (nmask > 0) read(*,*) mask
            call coo_nd_mask_assign(a,mask,scalar,status,accumulate=(add_flag/=0), &
                keep_zeros=(zeros_flag/=0),max_nnz=int(limit,int64))
            b = a
            deallocate(mask)
        case default
            error stop 'unknown operation'
        end select
        if (status /= sparse_success) then
            print '(i0)', status
            cycle
        end if
        flat = coo_nd_to_dense_flat(b,status)
        if (status /= sparse_success) error stop 'dense conversion failed'
        write(*,'(i0,1x,i0,1x,i0)',advance='no') status,size(b%shape),size(b%data)
        do ax = 1,size(b%shape)
            write(*,'(1x,i0)',advance='no') b%shape(ax)
        end do
        do k = 1,size(flat)
            write(*,'(1x,es25.16e3)',advance='no') flat(k)
        end do
        write(*,*)
    end do
contains
    subroutine read_array(z)
        type(coo_nd_array), intent(out) :: z !! Input COO array.
        integer, allocatable :: dims(:), coords(:,:)
        real(dp), allocatable :: vals(:)
        integer :: nd, n, j, info

        read(*,*) nd,n
        allocate(dims(nd),coords(nd,n),vals(n))
        read(*,*) dims
        do j = 1,n
            read(*,*) coords(:,j),vals(j)
        end do
        z = coo_nd_from_triplets(dims,coords,vals,info)
        if (info /= sparse_success) error stop 'bad COO input'
    end subroutine read_array
end program sparse_coo_nd_edit_reference_driver
