! SPDX-License-Identifier: MIT
! Streaming randomized reference driver for complement-aware N-D Boolean masks.
program sparse_coo_nd_bool_reference_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, valid_coo_nd
    use scifort_sparse_coo_nd_bool
    implicit none

    type(coo_nd_bool_mask) :: a, b, result_mask
    type(coo_nd_array) :: reals, result_coo
    integer, allocatable :: dims(:), points(:,:), left(:,:), right(:,:)
    integer :: nt, t, mode, ndim, na, nb, da, db, info, i, j, k
    integer(int64) :: ntrue
    character(16) :: op
    real(dp) :: threshold

    read(*,*) nt
    do t = 1, nt
        read(*,*) mode, ndim
        allocate(dims(ndim))
        read(*,*) dims
        read(*,*) na, nb, da, db
        allocate(left(ndim,na),right(ndim,nb))
        do i = 1, na
            read(*,*) left(:,i)
        end do
        do i = 1, nb
            read(*,*) right(:,i)
        end do
        a = coo_nd_bool_from_points(dims,left,info,default_value=da/=0)
        if (info /= sparse_success) error stop 'left point-mask construction'
        b = coo_nd_bool_from_points(dims,right,info,default_value=db/=0)
        if (info /= sparse_success) error stop 'right point-mask construction'
        select case (mode)
        case (1)
            result_mask = coo_nd_bool_and(a,b,info)
        case (2)
            result_mask = coo_nd_bool_or(a,b,info)
        case (3)
            result_mask = coo_nd_bool_xor(a,b,info)
        case (4)
            result_mask = coo_nd_bool_difference(a,b,info)
        case (5)
            result_mask = coo_nd_bool_not(a,info)
        case (6)
            call read_coo(reals)
            read(*,*) op, threshold
            result_mask = coo_nd_bool_compare_scalar(reals,trim(op),threshold,info)
        case (7)
            call read_coo(reals)
            result_coo = coo_nd_bool_filter(reals,a,info)
            call write_coo(result_coo,info)
            call release_all()
            cycle
        case (8)
            ntrue = coo_nd_bool_count_true(a,info)
            print '(i0,1x,i0,1x,i0)', info,-2,ntrue
            call release_all()
            cycle
        case (9)
            result_coo = coo_nd_bool_to_coo(a,info)
            call write_coo(result_coo,info)
            call release_all()
            cycle
        case (10)
            call read_coo(reals)
            result_mask = coo_nd_bool_from_coo(reals,info,complement=db/=0)
        case default
            error stop 'unknown mode'
        end select
        if (info /= sparse_success .or. .not. valid_coo_nd_bool_mask(result_mask)) then
            print '(i0,1x,i0)', info, -1
        else
            write(*,'(i0,1x,i0,1x,i0,1x,i0,1x)',advance='no') &
                info,ndim,size(result_mask%coords,2),merge(1,0,result_mask%default_value)
            do i = 1, ndim
                write(*,'(i0,1x)',advance='no') dims(i)
            end do
            do k = 1, size(result_mask%coords,2)
                do i = 1, ndim
                    write(*,'(i0,1x)',advance='no') result_mask%coords(i,k)
                end do
            end do
            write(*,*)
        end if
        call release_all()
    end do
contains
    subroutine release_all()
        deallocate(dims,left,right)
    end subroutine release_all

    subroutine read_coo(z)
        type(coo_nd_array), intent(out) :: z !! Read N-D real COO array.
        integer, allocatable :: ds(:), coords(:,:)
        real(dp), allocatable :: values(:)
        integer :: nd, nnz, ii, st

        read(*,*) nd,nnz
        allocate(ds(nd),coords(nd,nnz),values(nnz))
        read(*,*) ds
        do ii = 1, nnz
            read(*,*) coords(:,ii),values(ii)
        end do
        z = coo_nd_from_triplets(ds,coords,values,st)
        if (st /= sparse_success) error stop 'invalid COO input'
    end subroutine read_coo

    subroutine write_coo(z, st)
        type(coo_nd_array), intent(in) :: z !! COO result.
        integer, intent(in) :: st !! Operation status.
        integer :: ii, jj

        if (st /= sparse_success .or. .not. valid_coo_nd(z)) then
            print '(i0,1x,i0)', st,-1
            return
        end if
        write(*,'(i0,1x,i0,1x,i0,1x)',advance='no') st,size(z%shape),size(z%data)
        do ii = 1, size(z%shape)
            write(*,'(i0,1x)',advance='no') z%shape(ii)
        end do
        do jj = 1, size(z%data)
            do ii = 1, size(z%shape)
                write(*,'(i0,1x)',advance='no') z%coords(ii,jj)
            end do
            write(*,'(es25.16e3,1x)',advance='no') z%data(jj)
        end do
        write(*,*)
    end subroutine write_coo
end program sparse_coo_nd_bool_reference_driver
