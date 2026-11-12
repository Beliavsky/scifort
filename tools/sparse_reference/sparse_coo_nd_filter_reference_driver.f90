! SPDX-License-Identifier: MIT
! Streaming validation driver for sparse COO N-D comparisons and masks.
program sparse_coo_nd_filter_reference_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, valid_coo_nd
    use scifort_sparse_coo_nd_filter, only : coo_nd_compare_scalar, &
        coo_nd_filter_sparse, coo_nd_filter_points, coo_nd_filter_value, &
        coo_nd_count_nonzero
    implicit none
    type(coo_nd_array) :: a, b, mask
    integer, allocatable :: points(:,:)
    integer :: cases, caseid, mode, info, npoints, j, ax, inverse_flag
    integer(int64) :: n
    character(16) :: op
    real(dp) :: threshold

    read(*,*) cases
    do caseid = 1, cases
        read(*,*) mode
        call read_coo(a)
        select case (mode)
        case (1)
            read(*,*) op, threshold
            b = coo_nd_compare_scalar(a,trim(op),threshold,info)
        case (2)
            call read_coo(mask)
            read(*,*) inverse_flag
            b = coo_nd_filter_sparse(a,mask,info,complement=(inverse_flag/=0))
        case (3)
            read(*,*) npoints, inverse_flag
            allocate(points(size(a%shape),npoints))
            do j = 1, npoints
                read(*,*) points(:,j)
            end do
            b = coo_nd_filter_points(a,points,info,complement=(inverse_flag/=0))
            deallocate(points)
        case (4)
            read(*,*) op, threshold
            b = coo_nd_filter_value(a,trim(op),threshold,info)
        case (5)
            n = coo_nd_count_nonzero(a,info)
            print '(i0,1x,i0,1x,i0)', info,-2,n
            cycle
        case default
            error stop 'invalid test case'
        end select
        if (info /= sparse_success .or. .not. valid_coo_nd(b)) then
            print '(i0,1x,i0)', info, -1
            cycle
        end if
        write(*,'(i0,1x,i0,1x,i0,1x)',advance='no') info,size(b%shape),size(b%data)
        do ax = 1, size(b%shape)
            write(*,'(i0,1x)',advance='no') b%shape(ax)
        end do
        do j = 1, size(b%data)
            do ax = 1, size(b%shape)
                write(*,'(i0,1x)',advance='no') b%coords(ax,j)
            end do
            write(*,'(es25.16e3,1x)',advance='no') b%data(j)
        end do
        write(*,*)
    end do
contains
    subroutine read_coo(z)
        type(coo_nd_array), intent(out) :: z !! COO input.
        integer, allocatable :: dims(:), coords(:,:)
        real(dp), allocatable :: values(:)
        integer :: ndim, nnz, i, j, stat

        read(*,*) ndim,nnz
        allocate(dims(ndim),coords(ndim,nnz),values(nnz))
        read(*,*) dims
        do i = 1, nnz
            read(*,*) coords(:,i),values(i)
        end do
        z = coo_nd_from_triplets(dims,coords,values,stat)
        if (stat /= sparse_success) error stop 'invalid input COO'
    end subroutine read_coo
end program sparse_coo_nd_filter_reference_driver
