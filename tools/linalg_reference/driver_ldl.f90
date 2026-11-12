! SPDX-License-Identifier: MIT
! Streamed factorization for independent SciPy numerical comparisons.
program driver_ldl
    use scifort_kinds, only : dp
    use scifort_linalg_ldl, only : ldl
    implicit none
    real(dp), allocatable :: ar(:,:), lr(:,:), dr(:,:)
    complex(dp), allocatable :: a(:,:), lu(:,:), d(:,:)
    integer, allocatable :: perm(:)
    integer :: n, family, orientation, status, i, j, ios
    real(dp) :: re, im

    do
        read(*,*,iostat=ios) n, family, orientation
        if (ios /= 0) exit
        allocate(a(n,n),lu(n,n),d(n,n),perm(n),ar(n,n),lr(n,n),dr(n,n))
        do j = 1, n
            do i = 1, n
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
                ar(i,j) = re
            end do
        end do
        select case(family)
        case(1)
            call ldl(ar,lr,dr,perm,status,lower=(orientation == 1))
            lu = cmplx(lr,kind=dp)
            d = cmplx(dr,kind=dp)
        case(2)
            call ldl(a,lu,d,perm,status,lower=(orientation == 1))
        case(3)
            call ldl(a,lu,d,perm,status,lower=(orientation == 1),hermitian=.false.)
        case default
            status = -1
        end select
        print '(i0)', status
        do j = 1, n
            do i = 1, n
                print '(4(es26.17e3,1x))', real(lu(i,j),dp),aimag(lu(i,j)), &
                    real(d(i,j),dp),aimag(d(i,j))
            end do
        end do
        do i = 1, n
            print '(i0)', perm(i)
        end do
        deallocate(a,lu,d,perm,ar,lr,dr)
    end do
end program driver_ldl
