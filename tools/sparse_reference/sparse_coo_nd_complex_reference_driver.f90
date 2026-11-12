! SPDX-License-Identifier: MIT
! Streaming complex128 N-D COO reference driver for Python/SciPy comparisons.
program sparse_coo_nd_complex_reference_driver
    use, intrinsic :: iso_fortran_env, only : input_unit, int64
    use scifort_kinds, only : dp
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd_complex
    use scifort_sparse_coo_nd_complex_ops
    implicit none
    type(coo_nd_complex_array) :: a, b, c
    integer, allocatable :: dims(:), coords(:,:), axes(:), axes_b(:), target(:)
    complex(dp), allocatable :: values(:), flat(:)
    complex(dp) :: z
    integer :: trials, trial, mode, n, nnz, i, j, rank_a, k, keep, ierr
    integer(int64) :: cap
    character(len=1) :: order
    real(dp) :: x, y

    read(input_unit,*) trials
    do trial = 1,trials
        read(input_unit,*) mode
        call read_array(a)
        ierr = sparse_success
        select case(mode)
        case(1)
            c = a
        case(2)
            c = coo_nd_complex_conjugate(a,ierr)
        case(3)
            read(input_unit,*) rank_a
            allocate(axes(rank_a))
            read(input_unit,*) axes
            c = coo_nd_complex_transpose(a,axes,ierr)
        case(4)
            read(input_unit,*) rank_a,order
            allocate(target(rank_a))
            read(input_unit,*) target
            c = coo_nd_complex_reshape(a,target,order,ierr)
        case(5)
            read(input_unit,*) rank_a
            allocate(target(rank_a))
            read(input_unit,*) target
            c = coo_nd_complex_broadcast_to(a,target,ierr,max_nnz=200000_int64)
        case(6,7,9,10)
            call read_array(b)
            select case(mode)
            case(6)
                c = coo_nd_complex_add(a,b,ierr,max_nnz=200000_int64)
            case(7)
                c = coo_nd_complex_multiply(a,b,ierr,max_nnz=200000_int64)
            case(9)
                read(input_unit,*) k
                allocate(axes(k),axes_b(k))
                if (k > 0) then
                    read(input_unit,*) axes
                    read(input_unit,*) axes_b
                end if
                c = coo_nd_complex_tensordot(a,b,axes,axes_b,ierr,max_nnz=200000_int64)
            case(10)
                z = coo_nd_complex_vdot(a,b,ierr)
            end select
        case(8)
            read(input_unit,*) k,keep
            allocate(axes(k))
            if (k > 0) read(input_unit,*) axes
            c = coo_nd_complex_sum_axes(a,axes,keep/=0,ierr)
        case(11)
            c = coo_nd_complex_hermitian(a,ierr)
        case(12)
            z = coo_nd_complex_sum_all(a,ierr)
        case(13)
            read(input_unit,*) order
            flat = coo_nd_complex_to_dense_flat(a,ierr,order)
            if (ierr == sparse_success) c = coo_nd_complex_from_dense_flat(a%shape,flat,ierr,order)
        case default
            error stop 'invalid complex reference driver operation'
        end select
        if (mode == 10 .or. mode == 12) then
            write(*,'(i0,1x,i0,2(1x,es25.16e3))') ierr,0,real(z,dp),aimag(z)
        else
            if (ierr == sparse_success) then
                flat = coo_nd_complex_to_dense_flat(c,ierr,'C',max_elements=200000_int64)
            end if
            if (ierr /= sparse_success) then
                write(*,'(i0,1x,i0)') ierr,0
            else
                write(*,'(i0,1x,i0)',advance='no') ierr,size(c%shape)
                do i = 1,size(c%shape)
                    write(*,'(1x,i0)',advance='no') c%shape(i)
                end do
                do i = 1,size(flat)
                    write(*,'(2(1x,es25.16e3))',advance='no') real(flat(i),dp),aimag(flat(i))
                end do
                write(*,*)
            end if
        end if
        if (allocated(axes)) deallocate(axes)
        if (allocated(axes_b)) deallocate(axes_b)
        if (allocated(target)) deallocate(target)
    end do
contains
    subroutine read_array(t)
        type(coo_nd_complex_array), intent(out) :: t !! New complex sparse array.
        integer, allocatable :: shape_a(:), coord_a(:,:)
        complex(dp), allocatable :: value_a(:)
        integer :: rank_t, n_t, row, col, rc
        real(dp) :: re, im

        read(input_unit,*) rank_t,n_t
        allocate(shape_a(rank_t),coord_a(rank_t,n_t),value_a(n_t))
        read(input_unit,*) shape_a
        do col = 1,n_t
            read(input_unit,*) coord_a(:,col),re,im
            value_a(col) = cmplx(re,im,kind=dp)
        end do
        t = coo_nd_complex_from_triplets(shape_a,coord_a,value_a,rc)
        if (rc /= sparse_success) error stop 'invalid test tensor'
    end subroutine read_array
end program sparse_coo_nd_complex_reference_driver
