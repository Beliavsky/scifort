! SPDX-License-Identifier: MIT
! Driver for independent scipy.linalg.matmul_toeplitz comparisons.
! Input: n m nrhs kind explicit_row; then complex c, r (if present), x.
! kind: 0 real/real, 1 real/complex, 2 complex/real, 3 complex/complex.
program driver_matmul_toeplitz
    use scifort_kinds, only : dp
    use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
    implicit none
    complex(dp), allocatable :: c(:)
    complex(dp), allocatable :: r(:)
    complex(dp), allocatable :: x(:,:)
    complex(dp), allocatable :: y(:,:)
    real(dp) :: real_part
    real(dp) :: imaginary_part
    integer :: n
    integer :: m
    integer :: nrhs
    integer :: kind
    integer :: explicit_row
    integer :: status
    integer :: ios
    integer :: i
    integer :: j

    do
        read(*,*,iostat=ios) n, m, nrhs, kind, explicit_row
        if (ios /= 0) exit
        if (n < 0) exit
        allocate(c(n),r(m),x(m,nrhs),y(n,nrhs))
        do i=1,n
            read(*,*) real_part,imaginary_part
            c(i) = cmplx(real_part,imaginary_part,dp)
        end do
        if (explicit_row /= 0) then
            do i=1,m
                read(*,*) real_part,imaginary_part
                r(i) = cmplx(real_part,imaginary_part,dp)
            end do
        end if
        do j=1,nrhs
            do i=1,m
                read(*,*) real_part,imaginary_part
                x(i,j) = cmplx(real_part,imaginary_part,dp)
            end do
        end do
        call dispatch_matrix()
        print '(i0)',status
        if (status == 0) then
            do j=1,nrhs
                do i=1,n
                    write(*,'(es27.17e3,1x,es27.17e3)') real(y(i,j),dp),aimag(y(i,j))
                end do
            end do
        end if
        deallocate(c,r,x,y)
    end do

contains

    ! Dispatch matrix generic: all four real/complex combinations.
    subroutine dispatch_matrix()
        real(dp), allocatable :: ry(:,:)
        integer :: s
        select case(kind)
        case(0)
            allocate(ry(n,nrhs))
            if (explicit_row /= 0) then
                call matmul_toeplitz(real(c,dp),real(x,dp),ry,s,real(r,dp))
            else
                call matmul_toeplitz(real(c,dp),real(x,dp),ry,s)
            end if
            y = cmplx(ry,kind=dp)
        case(1)
            if (explicit_row /= 0) then
                call matmul_toeplitz(real(c,dp),x,y,s,real(r,dp))
            else
                call matmul_toeplitz(real(c,dp),x,y,s)
            end if
        case(2)
            if (explicit_row /= 0) then
                call matmul_toeplitz(c,real(x,dp),y,s,r)
            else
                call matmul_toeplitz(c,real(x,dp),y,s)
            end if
        case(3)
            if (explicit_row /= 0) then
                call matmul_toeplitz(c,x,y,s,r)
            else
                call matmul_toeplitz(c,x,y,s)
            end if
        case default
            s = 1
            y = (0.0_dp,0.0_dp)
        end select
        status = s
    end subroutine dispatch_matrix

end program driver_matmul_toeplitz
