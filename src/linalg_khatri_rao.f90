! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Column-wise Kronecker/Khatri-Rao products, matching scipy.linalg.khatri_rao.
! First input row is outer, second input row is inner:
! c((i-1)*m+j,k) = a(i,k)*b(j,k).
module scifort_linalg_khatri_rao
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private
    public :: khatri_rao

    interface khatri_rao
        module procedure khatri_rao_i32_i32_22
        module procedure khatri_rao_i32_i32_33
        module procedure khatri_rao_i32_i32_23
        module procedure khatri_rao_i32_i32_32
        module procedure khatri_rao_i32_i64_22
        module procedure khatri_rao_i32_i64_33
        module procedure khatri_rao_i32_i64_23
        module procedure khatri_rao_i32_i64_32
        module procedure khatri_rao_i32_r_22
        module procedure khatri_rao_i32_r_33
        module procedure khatri_rao_i32_r_23
        module procedure khatri_rao_i32_r_32
        module procedure khatri_rao_i32_c_22
        module procedure khatri_rao_i32_c_33
        module procedure khatri_rao_i32_c_23
        module procedure khatri_rao_i32_c_32
        module procedure khatri_rao_i64_i32_22
        module procedure khatri_rao_i64_i32_33
        module procedure khatri_rao_i64_i32_23
        module procedure khatri_rao_i64_i32_32
        module procedure khatri_rao_i64_i64_22
        module procedure khatri_rao_i64_i64_33
        module procedure khatri_rao_i64_i64_23
        module procedure khatri_rao_i64_i64_32
        module procedure khatri_rao_i64_r_22
        module procedure khatri_rao_i64_r_33
        module procedure khatri_rao_i64_r_23
        module procedure khatri_rao_i64_r_32
        module procedure khatri_rao_i64_c_22
        module procedure khatri_rao_i64_c_33
        module procedure khatri_rao_i64_c_23
        module procedure khatri_rao_i64_c_32
        module procedure khatri_rao_r_i32_22
        module procedure khatri_rao_r_i32_33
        module procedure khatri_rao_r_i32_23
        module procedure khatri_rao_r_i32_32
        module procedure khatri_rao_r_i64_22
        module procedure khatri_rao_r_i64_33
        module procedure khatri_rao_r_i64_23
        module procedure khatri_rao_r_i64_32
        module procedure khatri_rao_r_r_22
        module procedure khatri_rao_r_r_33
        module procedure khatri_rao_r_r_23
        module procedure khatri_rao_r_r_32
        module procedure khatri_rao_r_c_22
        module procedure khatri_rao_r_c_33
        module procedure khatri_rao_r_c_23
        module procedure khatri_rao_r_c_32
        module procedure khatri_rao_c_i32_22
        module procedure khatri_rao_c_i32_33
        module procedure khatri_rao_c_i32_23
        module procedure khatri_rao_c_i32_32
        module procedure khatri_rao_c_i64_22
        module procedure khatri_rao_c_i64_33
        module procedure khatri_rao_c_i64_23
        module procedure khatri_rao_c_i64_32
        module procedure khatri_rao_c_r_22
        module procedure khatri_rao_c_r_33
        module procedure khatri_rao_c_r_23
        module procedure khatri_rao_c_r_32
        module procedure khatri_rao_c_c_22
        module procedure khatri_rao_c_c_33
        module procedure khatri_rao_c_c_23
        module procedure khatri_rao_c_c_32
    end interface khatri_rao

contains

    pure subroutine khatri_rao_i32_i32_22(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        integer(int32), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0_int32
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i32_i32_22

    pure subroutine khatri_rao_i32_i32_33(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        integer(int32), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int32
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_i32_33

    pure subroutine khatri_rao_i32_i32_23(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        integer(int32), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int32
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_i32_23

    pure subroutine khatri_rao_i32_i32_32(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        integer(int32), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int32
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_i32_32

    pure subroutine khatri_rao_i32_i64_22(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        integer(int64), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i32_i64_22

    pure subroutine khatri_rao_i32_i64_33(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_i64_33

    pure subroutine khatri_rao_i32_i64_23(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_i64_23

    pure subroutine khatri_rao_i32_i64_32(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_i64_32

    pure subroutine khatri_rao_i32_r_22(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i32_r_22

    pure subroutine khatri_rao_i32_r_33(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_r_33

    pure subroutine khatri_rao_i32_r_23(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_r_23

    pure subroutine khatri_rao_i32_r_32(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_r_32

    pure subroutine khatri_rao_i32_c_22(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i32_c_22

    pure subroutine khatri_rao_i32_c_33(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_c_33

    pure subroutine khatri_rao_i32_c_23(a, b, c, status)
        integer(int32), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_c_23

    pure subroutine khatri_rao_i32_c_32(a, b, c, status)
        integer(int32), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i32_c_32

    pure subroutine khatri_rao_i64_i32_22(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        integer(int64), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i64_i32_22

    pure subroutine khatri_rao_i64_i32_33(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_i32_33

    pure subroutine khatri_rao_i64_i32_23(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_i32_23

    pure subroutine khatri_rao_i64_i32_32(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_i32_32

    pure subroutine khatri_rao_i64_i64_22(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        integer(int64), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i64_i64_22

    pure subroutine khatri_rao_i64_i64_33(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_i64_33

    pure subroutine khatri_rao_i64_i64_23(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_i64_23

    pure subroutine khatri_rao_i64_i64_32(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        integer(int64), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0_int64
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_i64_32

    pure subroutine khatri_rao_i64_r_22(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i64_r_22

    pure subroutine khatri_rao_i64_r_33(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_r_33

    pure subroutine khatri_rao_i64_r_23(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_r_23

    pure subroutine khatri_rao_i64_r_32(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_r_32

    pure subroutine khatri_rao_i64_c_22(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_i64_c_22

    pure subroutine khatri_rao_i64_c_33(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_c_33

    pure subroutine khatri_rao_i64_c_23(a, b, c, status)
        integer(int64), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_c_23

    pure subroutine khatri_rao_i64_c_32(a, b, c, status)
        integer(int64), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_i64_c_32

    pure subroutine khatri_rao_r_i32_22(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_r_i32_22

    pure subroutine khatri_rao_r_i32_33(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_i32_33

    pure subroutine khatri_rao_r_i32_23(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_i32_23

    pure subroutine khatri_rao_r_i32_32(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_i32_32

    pure subroutine khatri_rao_r_i64_22(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_r_i64_22

    pure subroutine khatri_rao_r_i64_33(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_i64_33

    pure subroutine khatri_rao_r_i64_23(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_i64_23

    pure subroutine khatri_rao_r_i64_32(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_i64_32

    pure subroutine khatri_rao_r_r_22(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_r_r_22

    pure subroutine khatri_rao_r_r_33(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_r_33

    pure subroutine khatri_rao_r_r_23(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_r_23

    pure subroutine khatri_rao_r_r_32(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        real(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_r_32

    pure subroutine khatri_rao_r_c_22(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_r_c_22

    pure subroutine khatri_rao_r_c_33(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_c_33

    pure subroutine khatri_rao_r_c_23(a, b, c, status)
        real(dp), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_c_23

    pure subroutine khatri_rao_r_c_32(a, b, c, status)
        real(dp), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_r_c_32

    pure subroutine khatri_rao_c_i32_22(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_c_i32_22

    pure subroutine khatri_rao_c_i32_33(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_i32_33

    pure subroutine khatri_rao_c_i32_23(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        integer(int32), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_i32_23

    pure subroutine khatri_rao_c_i32_32(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int32), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_i32_32

    pure subroutine khatri_rao_c_i64_22(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_c_i64_22

    pure subroutine khatri_rao_c_i64_33(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_i64_33

    pure subroutine khatri_rao_c_i64_23(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        integer(int64), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_i64_23

    pure subroutine khatri_rao_c_i64_32(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        integer(int64), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_i64_32

    pure subroutine khatri_rao_c_r_22(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_c_r_22

    pure subroutine khatri_rao_c_r_33(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_r_33

    pure subroutine khatri_rao_c_r_23(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        real(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_r_23

    pure subroutine khatri_rao_c_r_32(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        real(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_r_32

    pure subroutine khatri_rao_c_c_22(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,1)
        k = size(a,2)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        if (size(c,1) /= n*m .or. size(c,2) /= k) return
        status = linalg_status_success
        do col=1,k
            do i=1,n
                do j=1,m
                    c((i-1)*m+j,col) = a(i,col)*b(j,col)
                end do
            end do
        end do
    end subroutine khatri_rao_c_c_22

    pure subroutine khatri_rao_c_c_33(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,2)
        k = size(a,3)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_c_33

    pure subroutine khatri_rao_c_c_23(a, b, c, status)
        complex(dp), intent(in) :: a(:,:) !! first factor
        complex(dp), intent(in) :: b(:,:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        m = size(b,2)
        k = size(a,2)
        if (size(b,3) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        bb = size(b,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (bb /= 1) ib = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(i,col)*b(ib,j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_c_23

    pure subroutine khatri_rao_c_c_32(a, b, c, status)
        complex(dp), intent(in) :: a(:,:,:) !! first factor
        complex(dp), intent(in) :: b(:,:) !! second factor
        complex(dp), intent(out) :: c(:,:,:) !! columnwise Kronecker result
        integer, intent(out) :: status !! success or invalid_shape
        integer :: n, m, k, i, j, col
        integer :: batch, ba, bb, ia, ib

        c = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,2)
        m = size(b,1)
        k = size(a,3)
        if (size(b,2) /= k) return
        if (int(n,int64)*int(m,int64) > int(huge(n),int64)) return
        ba = 1
        bb = 1
        ba = size(a,1)
        if (ba == 0 .or. bb == 0) return  ! scipy disallows empty batches
        if (ba /= bb .and. ba /= 1 .and. bb /= 1) return
        batch = max(ba,bb)
        if (size(c,1) /= batch .or. size(c,2) /= n*m .or. size(c,3) /= k) return
        status = linalg_status_success
        do batch=1,size(c,1)
            ia = 1
            ib = 1
            if (ba /= 1) ia = batch
        do col=1,k
            do i=1,n
                do j=1,m
                    c(batch,(i-1)*m+j,col) = a(ia,i,col)*b(j,col)
                end do
            end do
        end do
        end do
    end subroutine khatri_rao_c_c_32

end module scifort_linalg_khatri_rao
