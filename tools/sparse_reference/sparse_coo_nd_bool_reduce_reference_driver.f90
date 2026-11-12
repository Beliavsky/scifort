! SPDX-License-Identifier: MIT
! Read randomized cases, evaluate sparse Boolean operations, and emit
! flattened values (C-order) for independent Python/NumPy comparison.
program sparse_coo_nd_bool_reduce_reference_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd_bool, only : coo_nd_bool_mask, coo_nd_bool_from_points, &
        coo_nd_bool_get
    use scifort_sparse_coo_nd_bool_reduce
    implicit none

    type(coo_nd_bool_mask) :: a, b, out
    type(coo_nd_bool_count_array) :: counts
    integer, allocatable :: dims(:), target(:), points(:,:), coord(:)
    integer :: trials, t, mode, rank, rank_b, npts, ax, keep, op, info, i
    integer(int64) :: n, index, work, v
    logical :: def, full_any, full_all, val

    read(*,*) trials
    do t = 1, trials
        read(*,*) mode, rank
        allocate(dims(rank))
        read(*,*) dims
        read(*,*) npts, keep
        def = keep /= 0
        allocate(points(rank,npts))
        do i = 1, npts
            read(*,*) points(:,i)
        end do
        a = coo_nd_bool_from_points(dims,points,info,default_value=def)
        if (info /= sparse_success) error stop 'invalid randomized input mask'
        deallocate(dims,points)
        select case(mode)
        case (1)
            read(*,*) rank_b
            allocate(target(rank_b))
            read(*,*) target
            out = coo_nd_bool_broadcast_to(a,target,info)
            if (info /= sparse_success) error stop 'broadcast operation failure'
            call emit_mask(out)
            deallocate(target)
        case (2,3,4)
            read(*,*) ax, keep
            select case(mode)
            case(2)
                out = coo_nd_bool_any_axis(a,ax,info,keepdims=keep/=0)
                if (info /= sparse_success) error stop 'axis any operation failure'
                call emit_mask(out)
            case(3)
                out = coo_nd_bool_all_axis(a,ax,info,keepdims=keep/=0)
                if (info /= sparse_success) error stop 'axis all operation failure'
                call emit_mask(out)
            case(4)
                counts = coo_nd_bool_count_axis(a,ax,info,keepdims=keep/=0)
                if (info /= sparse_success) error stop 'axis count operation failure'
                call emit_counts(counts)
            end select
        case(5)
            read(*,*) rank_b
            allocate(dims(rank_b))
            read(*,*) dims
            read(*,*) npts, keep
            def = keep /= 0
            allocate(points(rank_b,npts))
            do i = 1, npts
                read(*,*) points(:,i)
            end do
            b = coo_nd_bool_from_points(dims,points,info,default_value=def)
            if (info /= sparse_success) error stop 'invalid second mask'
            deallocate(dims,points)
            read(*,*) op
            select case(op)
            case(1)
                out = coo_nd_bool_and_broadcast(a,b,info)
            case(2)
                out = coo_nd_bool_or_broadcast(a,b,info)
            case(3)
                out = coo_nd_bool_xor_broadcast(a,b,info)
            case(4)
                out = coo_nd_bool_difference_broadcast(a,b,info)
            end select
            if (info /= sparse_success) error stop 'binary broadcast failure'
            call emit_mask(out)
        case(6)
            full_any = coo_nd_bool_any(a,info)
            if (info /= sparse_success) error stop 'full any failure'
            full_all = coo_nd_bool_all(a,info)
            if (info /= sparse_success) error stop 'full all failure'
            write(*,'(2(I0,1X))') merge(1,0,full_any),merge(1,0,full_all)
        case default
            error stop 'unknown operation code'
        end select
    end do

contains

    subroutine emit_mask(mask)
        type(coo_nd_bool_mask), intent(in) :: mask !! Evaluated sparse mask.
        integer, allocatable :: pt(:)
        integer(int64) :: extent, offset, remainder
        integer :: j, status
        logical :: flag

        extent = 1_int64
        do j = 1, size(mask%shape)
            extent = extent*int(mask%shape(j),int64)
        end do
        write(*,'(I0,1X)',advance='no') size(mask%shape)
        do j = 1, size(mask%shape)
            write(*,'(I0,1X)',advance='no') mask%shape(j)
        end do
        write(*,'(I0,1X)',advance='no') extent
        allocate(pt(size(mask%shape)))
        do offset = 0_int64, extent-1_int64
            remainder = offset
            do j = size(mask%shape),1,-1
                pt(j) = int(mod(remainder,int(mask%shape(j),int64)))+1
                remainder = remainder/int(mask%shape(j),int64)
            end do
            flag = coo_nd_bool_get(mask,pt,status)
            if (status /= sparse_success) error stop 'invalid mask lookup'
            write(*,'(I0,1X)',advance='no') merge(1,0,flag)
        end do
        write(*,*)
    end subroutine emit_mask

    subroutine emit_counts(counts)
        type(coo_nd_bool_count_array), intent(in) :: counts !! Evaluated count array.
        integer, allocatable :: pt(:)
        integer(int64) :: extent, offset, remainder, value
        integer :: j, status

        extent = 1_int64
        do j = 1, size(counts%shape)
            extent = extent*int(counts%shape(j),int64)
        end do
        write(*,'(I0,1X)',advance='no') size(counts%shape)
        do j = 1, size(counts%shape)
            write(*,'(I0,1X)',advance='no') counts%shape(j)
        end do
        write(*,'(I0,1X)',advance='no') extent
        allocate(pt(size(counts%shape)))
        do offset = 0_int64, extent-1_int64
            remainder = offset
            do j = size(counts%shape),1,-1
                pt(j) = int(mod(remainder,int(counts%shape(j),int64)))+1
                remainder = remainder/int(counts%shape(j),int64)
            end do
            value = coo_nd_bool_count_get(counts,pt,status)
            if (status /= sparse_success) error stop 'invalid count lookup'
            write(*,'(I0,1X)',advance='no') value
        end do
        write(*,*)
    end subroutine emit_counts

end program sparse_coo_nd_bool_reduce_reference_driver
