! SPDX-License-Identifier: MIT
! Stream reference cases to the N-dimensional COO extension.
program sparse_coo_nd_ops_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, coo_nd_to_dense_flat
    use scifort_sparse_coo_nd_ops, only : coo_nd_sum_all, coo_nd_sum_axes, &
        coo_nd_broadcast_to, coo_nd_add, coo_nd_multiply, coo_nd_tensordot
    implicit none
    type(coo_nd_array) :: a, b, c
    real(dp), allocatable :: dense(:)
    real(dp) :: scalar
    integer, allocatable :: axes(:), axes_b(:), target(:)
    integer :: cases, case_id, mode, n, status, keep, j

    read(*,*) cases
    do case_id = 1, cases
        read(*,*) mode
        call read_array(a)
        select case (mode)
        case (1)
            read(*,*) n, keep
            allocate(axes(n))
            if (n > 0) read(*,*) axes
            c = coo_nd_sum_axes(a, axes,keep/=0,status)
        case (2)
            read(*,*) n
            allocate(target(n))
            read(*,*) target
            c = coo_nd_broadcast_to(a,target,status)
        case (3,4)
            call read_array(b)
            if (mode == 3) then
                c = coo_nd_add(a,b,status)
            else
                c = coo_nd_multiply(a,b,status)
            end if
        case (5)
            call read_array(b)
            read(*,*) n
            allocate(axes(n),axes_b(n))
            if (n > 0) then
                read(*,*) axes
                read(*,*) axes_b
            end if
            c = coo_nd_tensordot(a,b,axes,axes_b,status)
        case (6)
            scalar = coo_nd_sum_all(a,status)
        case default
            error stop 'unexpected operation'
        end select
        if (status /= sparse_success) then
            write(*,'(i0)') status
        else if (mode == 6) then
            write(*,'(i0,1x,i0,1x,es25.16e3)') status,0,scalar
        else
            dense = coo_nd_to_dense_flat(c,status)
            if (status /= sparse_success) error stop 'dense conversion failed'
            write(*,'(i0,1x,i0)',advance='no') status,size(c%shape)
            do j = 1,size(c%shape)
                write(*,'(1x,i0)',advance='no') c%shape(j)
            end do
            do j = 1,size(dense)
                write(*,'(1x,es25.16e3)',advance='no') dense(j)
            end do
            write(*,*)
        end if
        if (allocated(axes)) deallocate(axes)
        if (allocated(axes_b)) deallocate(axes_b)
        if (allocated(target)) deallocate(target)
    end do
contains
    subroutine read_array(z)
        type(coo_nd_array), intent(out) :: z !! Read a new COO array.
        integer, allocatable :: dims(:), coords(:,:)
        real(dp), allocatable :: values(:)
        integer :: nr, nnz, k, ierr

        read(*,*) nr, nnz
        allocate(dims(nr),coords(nr,nnz),values(nnz))
        read(*,*) dims
        do k = 1,nnz
            read(*,*) coords(:,k),values(k)
        end do
        z = coo_nd_from_triplets(dims,coords,values,ierr)
        if (ierr /= sparse_success) error stop 'invalid case input'
    end subroutine read_array
end program sparse_coo_nd_ops_reference_driver
