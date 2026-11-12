! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Heterogeneous block diagonal assembly matching scipy.linalg.block_diag.
! Use as_block(value) for scalars, vectors, matrices or one-batch tensors.
! Output is provided by the caller to control type and allocation.
module scifort_linalg_block_diag
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private
    integer, parameter :: kind_integer = 1, kind_real = 2, kind_complex = 3

    type, public :: linalg_block
        private
        integer :: kind_code = 0
        integer :: rows = 0, cols = 0, batches = 1
        logical :: explicitly_batched = .false.
        integer(int64), allocatable :: integers(:,:,:)
        real(dp), allocatable :: reals(:,:,:)
        complex(dp), allocatable :: complexes(:,:,:)
    end type linalg_block

    public :: as_block, block_diag

    interface as_block
        module procedure as_block_i32_0
        module procedure as_block_i32_1
        module procedure as_block_i32_2
        module procedure as_block_i32_3
        module procedure as_block_i64_0
        module procedure as_block_i64_1
        module procedure as_block_i64_2
        module procedure as_block_i64_3
        module procedure as_block_r_0
        module procedure as_block_r_1
        module procedure as_block_r_2
        module procedure as_block_r_3
        module procedure as_block_c_0
        module procedure as_block_c_1
        module procedure as_block_c_2
        module procedure as_block_c_3
    end interface as_block

    interface block_diag
        module procedure block_diag_i32_2
        module procedure block_diag_i32_3
        module procedure block_diag_i64_2
        module procedure block_diag_i64_3
        module procedure block_diag_r_2
        module procedure block_diag_r_3
        module procedure block_diag_c_2
        module procedure block_diag_c_3
    end interface block_diag

contains

    pure function as_block_i32_0(a) result(block)
        integer(int32), intent(in) :: a !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%rows = 1
        block%cols = 1
        allocate(block%integers(1,1,1))
        block%integers(1,1,1) = int(a,int64)
    end function as_block_i32_0

    pure function as_block_i32_1(a) result(block)
        integer(int32), intent(in) :: a(:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%rows = 1
        block%cols = size(a)
        allocate(block%integers(1,1,size(a)))
        block%integers(1,1,:) = int(a,int64)
    end function as_block_i32_1

    pure function as_block_i32_2(a) result(block)
        integer(int32), intent(in) :: a(:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%rows = size(a,1)
        block%cols = size(a,2)
        allocate(block%integers(1,size(a,1),size(a,2)))
        block%integers(1,:,:) = int(a,int64)
    end function as_block_i32_2

    pure function as_block_i32_3(a) result(block)
        integer(int32), intent(in) :: a(:,:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%explicitly_batched = .true.
        block%batches = size(a,1)
        block%rows = size(a,2)
        block%cols = size(a,3)
        allocate(block%integers(size(a,1),size(a,2),size(a,3)))
        block%integers = int(a,int64)
    end function as_block_i32_3

    pure function as_block_i64_0(a) result(block)
        integer(int64), intent(in) :: a !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%rows = 1
        block%cols = 1
        allocate(block%integers(1,1,1))
        block%integers(1,1,1) = a
    end function as_block_i64_0

    pure function as_block_i64_1(a) result(block)
        integer(int64), intent(in) :: a(:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%rows = 1
        block%cols = size(a)
        allocate(block%integers(1,1,size(a)))
        block%integers(1,1,:) = a
    end function as_block_i64_1

    pure function as_block_i64_2(a) result(block)
        integer(int64), intent(in) :: a(:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%rows = size(a,1)
        block%cols = size(a,2)
        allocate(block%integers(1,size(a,1),size(a,2)))
        block%integers(1,:,:) = a
    end function as_block_i64_2

    pure function as_block_i64_3(a) result(block)
        integer(int64), intent(in) :: a(:,:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 1
        block%explicitly_batched = .true.
        block%batches = size(a,1)
        block%rows = size(a,2)
        block%cols = size(a,3)
        allocate(block%integers(size(a,1),size(a,2),size(a,3)))
        block%integers = a
    end function as_block_i64_3

    pure function as_block_r_0(a) result(block)
        real(dp), intent(in) :: a !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 2
        block%rows = 1
        block%cols = 1
        allocate(block%reals(1,1,1))
        block%reals(1,1,1) = a
    end function as_block_r_0

    pure function as_block_r_1(a) result(block)
        real(dp), intent(in) :: a(:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 2
        block%rows = 1
        block%cols = size(a)
        allocate(block%reals(1,1,size(a)))
        block%reals(1,1,:) = a
    end function as_block_r_1

    pure function as_block_r_2(a) result(block)
        real(dp), intent(in) :: a(:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 2
        block%rows = size(a,1)
        block%cols = size(a,2)
        allocate(block%reals(1,size(a,1),size(a,2)))
        block%reals(1,:,:) = a
    end function as_block_r_2

    pure function as_block_r_3(a) result(block)
        real(dp), intent(in) :: a(:,:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 2
        block%explicitly_batched = .true.
        block%batches = size(a,1)
        block%rows = size(a,2)
        block%cols = size(a,3)
        allocate(block%reals(size(a,1),size(a,2),size(a,3)))
        block%reals = a
    end function as_block_r_3

    pure function as_block_c_0(a) result(block)
        complex(dp), intent(in) :: a !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 3
        block%rows = 1
        block%cols = 1
        allocate(block%complexes(1,1,1))
        block%complexes(1,1,1) = a
    end function as_block_c_0

    pure function as_block_c_1(a) result(block)
        complex(dp), intent(in) :: a(:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 3
        block%rows = 1
        block%cols = size(a)
        allocate(block%complexes(1,1,size(a)))
        block%complexes(1,1,:) = a
    end function as_block_c_1

    pure function as_block_c_2(a) result(block)
        complex(dp), intent(in) :: a(:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 3
        block%rows = size(a,1)
        block%cols = size(a,2)
        allocate(block%complexes(1,size(a,1),size(a,2)))
        block%complexes(1,:,:) = a
    end function as_block_c_2

    pure function as_block_c_3(a) result(block)
        complex(dp), intent(in) :: a(:,:,:) !! numeric scalar, vector, matrix, or batch
        type(linalg_block) :: block !! tagged block for assembly

        block%kind_code = 3
        block%explicitly_batched = .true.
        block%batches = size(a,1)
        block%rows = size(a,2)
        block%cols = size(a,3)
        allocate(block%complexes(size(a,1),size(a,2),size(a,3)))
        block%complexes = a
    end function as_block_c_3

    pure subroutine validate_blocks(blocks, output_kind, output_rank, rows, cols, batches, status)
        type(linalg_block), intent(in) :: blocks(:) !! tagged input matrices
        integer, intent(in) :: output_kind !! 1=int32, 2=int64, 3=real, 4=complex
        integer, intent(in) :: output_rank !! 2=unbatched, 3=batch
        integer, intent(out) :: rows !! total rows
        integer, intent(out) :: cols !! total columns
        integer, intent(out) :: batches !! broadcasted batch count
        integer, intent(out) :: status !! success or invalid_shape/type
        integer(int64) :: total_rows, total_cols
        integer :: i, bs

        rows = 0
        cols = 0
        batches = 1
        status = linalg_status_invalid_shape
        total_rows = 0_int64
        total_cols = 0_int64
        if (size(blocks) == 0) total_rows = 1_int64
        do i=1,size(blocks)
            if (blocks(i)%kind_code < kind_integer .or. blocks(i)%kind_code > kind_complex) return
            if (blocks(i)%rows < 0 .or. blocks(i)%cols < 0) return
            if (blocks(i)%kind_code == kind_real .and. output_kind <= 2) return
            if (blocks(i)%kind_code == kind_complex .and. output_kind <= 3) return
            if (blocks(i)%kind_code == kind_integer .and. output_kind == 1) then
                if (any(blocks(i)%integers > int(huge(0_int32),int64))) return
                if (any(blocks(i)%integers < -int(huge(0_int32),int64)-1_int64)) return
            end if
            if (output_rank == 2 .and. blocks(i)%explicitly_batched) return
            if (blocks(i)%explicitly_batched) then
                bs = blocks(i)%batches
                if (bs /= 1) then
                    if (batches == 1) then
                        batches = bs
                    else if (bs /= batches) then
                        return
                    end if
                end if
            end if
            total_rows = total_rows + int(blocks(i)%rows,int64)
            total_cols = total_cols + int(blocks(i)%cols,int64)
            if (total_rows > int(huge(rows),int64)) return
            if (total_cols > int(huge(cols),int64)) return
        end do
        rows = int(total_rows)
        cols = int(total_cols)
        status = linalg_status_success
    end subroutine validate_blocks

    pure subroutine block_diag_i32_2(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        integer(int32), intent(out) :: a(:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = 0_int32
        call validate_blocks(blocks,1,2,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= rows .or. size(a,2) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = int(blocks(i)%integers(block_batch,:,:),int32)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_i32_2

    pure subroutine block_diag_i32_3(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        integer(int32), intent(out) :: a(:,:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = 0_int32
        call validate_blocks(blocks,1,3,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= batches .or. size(a,2) /= rows .or. size(a,3) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = int(blocks(i)%integers(block_batch,:,:),int32)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_i32_3

    pure subroutine block_diag_i64_2(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        integer(int64), intent(out) :: a(:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = 0_int64
        call validate_blocks(blocks,2,2,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= rows .or. size(a,2) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = blocks(i)%integers(block_batch,:,:)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_i64_2

    pure subroutine block_diag_i64_3(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        integer(int64), intent(out) :: a(:,:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = 0_int64
        call validate_blocks(blocks,2,3,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= batches .or. size(a,2) /= rows .or. size(a,3) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = blocks(i)%integers(block_batch,:,:)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_i64_3

    pure subroutine block_diag_r_2(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        real(dp), intent(out) :: a(:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = 0.0_dp
        call validate_blocks(blocks,3,2,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= rows .or. size(a,2) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = real(blocks(i)%integers(block_batch,:,:),dp)
                case (kind_real)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = blocks(i)%reals(block_batch,:,:)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_r_2

    pure subroutine block_diag_r_3(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        real(dp), intent(out) :: a(:,:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = 0.0_dp
        call validate_blocks(blocks,3,3,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= batches .or. size(a,2) /= rows .or. size(a,3) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = real(blocks(i)%integers(block_batch,:,:),dp)
                case (kind_real)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = blocks(i)%reals(block_batch,:,:)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_r_3

    pure subroutine block_diag_c_2(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        complex(dp), intent(out) :: a(:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = cmplx(0.0_dp,0.0_dp,dp)
        call validate_blocks(blocks,4,2,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= rows .or. size(a,2) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = &
                        cmplx(real(blocks(i)%integers(block_batch,:,:),dp),0.0_dp,dp)
                case (kind_real)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = cmplx(blocks(i)%reals(block_batch,:,:),kind=dp)
                case (kind_complex)
                    a(row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = blocks(i)%complexes(block_batch,:,:)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_c_2

    pure subroutine block_diag_c_3(blocks, a, status)
        type(linalg_block), intent(in) :: blocks(:) !! list of tagged blocks
        complex(dp), intent(out) :: a(:,:,:) !! block-diagonal output
        integer, intent(out) :: status !! success or invalid_shape/type
        integer :: rows, cols, batches, i, q, row0, col0, block_batch

        a = cmplx(0.0_dp,0.0_dp,dp)
        call validate_blocks(blocks,4,3,rows,cols,batches,status)
        if (status /= linalg_status_success) return
        if (size(a,1) /= batches .or. size(a,2) /= rows .or. size(a,3) /= cols) then
            status = linalg_status_invalid_shape
            return
        end if

        row0 = 0
        col0 = 0
        do i=1,size(blocks)
            do q=1,batches
                block_batch = 1
                if (blocks(i)%explicitly_batched .and. blocks(i)%batches /= 1) block_batch = q
                select case (blocks(i)%kind_code)
                case (kind_integer)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = &
                        cmplx(real(blocks(i)%integers(block_batch,:,:),dp),0.0_dp,dp)
                case (kind_real)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = cmplx(blocks(i)%reals(block_batch,:,:),kind=dp)
                case (kind_complex)
                    a(q,row0+1:row0+blocks(i)%rows,col0+1:col0+blocks(i)%cols) = blocks(i)%complexes(block_batch,:,:)
                end select
            end do
            row0 = row0 + blocks(i)%rows
            col0 = col0 + blocks(i)%cols
        end do
    end subroutine block_diag_c_3

end module scifort_linalg_block_diag
