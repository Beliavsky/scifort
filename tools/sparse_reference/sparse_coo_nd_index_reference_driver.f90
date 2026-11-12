! SPDX-License-Identifier: MIT
! Streaming SciPy/NumPy reference driver for N-D COO indexing and slicing.
program sparse_coo_nd_index_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, coo_nd_to_dense_flat
    use scifort_sparse_coo_nd_index, only : coo_nd_selector, coo_nd_all, &
        coo_nd_at, coo_nd_slice, coo_nd_indices, coo_nd_select, coo_nd_gather
    implicit none
    type(coo_nd_array) :: a, b
    type(coo_nd_selector), allocatable :: selectors(:)
    integer, allocatable :: points(:,:), idx(:)
    real(dp), allocatable :: values(:), flattened(:)
    integer :: cases, id, mode, nnz, ndim, ax, kind, first, last, stride, &
        has_first, has_last, count_idx, k, status, keepdims

    read(*,*) cases
    do id = 1,cases
        read(*,*) mode
        call read_array(a)
        ndim = size(a%shape)
        select case (mode)
        case (1)
            read(*,*) keepdims
            allocate(selectors(ndim))
            do ax = 1,ndim
                read(*,*) kind,first,last,stride,has_first,has_last,count_idx
                select case (kind)
                case (0)
                    selectors(ax) = coo_nd_all()
                case (1)
                    selectors(ax) = coo_nd_at(first)
                case (2)
                    if (has_first /= 0) then
                        if (has_last /= 0) then
                            selectors(ax) = coo_nd_slice(first,last,stride)
                        else
                            selectors(ax) = coo_nd_slice(first=first,stride=stride)
                        end if
                    else if (has_last /= 0) then
                        selectors(ax) = coo_nd_slice(last=last,stride=stride)
                    else
                        selectors(ax) = coo_nd_slice(stride=stride)
                    end if
                case (3)
                    allocate(idx(count_idx))
                    if (count_idx > 0) read(*,*) idx
                    selectors(ax) = coo_nd_indices(idx)
                    deallocate(idx)
                case default
                    error stop 'unexpected selector'
                end select
            end do
            b = coo_nd_select(a,selectors,status,keepdims/=0)
            deallocate(selectors)
        case (2)
            read(*,*) nnz
            allocate(points(ndim,nnz))
            do k = 1,nnz
                read(*,*) points(:,k)
            end do
            b = coo_nd_gather(a,points,status)
            deallocate(points)
        case default
            error stop 'unexpected operation'
        end select
        if (status /= sparse_success) then
            write(*,'(i0)') status
            cycle
        end if
        flattened = coo_nd_to_dense_flat(b,status)
        if (status /= sparse_success) error stop 'dense conversion failed'
        write(*,'(i0,1x,i0)',advance='no') status,size(b%shape)
        do ax = 1,size(b%shape)
            write(*,'(1x,i0)',advance='no') b%shape(ax)
        end do
        do k = 1,size(flattened)
            write(*,'(1x,es25.16e3)',advance='no') flattened(k)
        end do
        write(*,*)
    end do
contains
    subroutine read_array(z)
        type(coo_nd_array), intent(out) :: z !! Input COO array.
        integer, allocatable :: dims(:), coords(:,:)
        real(dp), allocatable :: vals(:)
        integer :: nd, n, j, info

        read(*,*) nd, n
        allocate(dims(nd),coords(nd,n),vals(n))
        read(*,*) dims
        do j = 1,n
            read(*,*) coords(:,j),vals(j)
        end do
        z = coo_nd_from_triplets(dims,coords,vals,info)
        if (info /= sparse_success) error stop 'bad reference input'
    end subroutine read_array
end program sparse_coo_nd_index_reference_driver
