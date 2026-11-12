! SPDX-License-Identifier: MIT
! Input: count, then operation, datatype and order, then n*n (real,imag) values in row order.
! Output: status and n*n (real,imag) values in row order.
program driver_trigonometric
    use scifort_kinds, only : dp
    use scifort_linalg_trigonometric, only : sinm, cosm, tanm, sinhm, coshm, tanhm
    implicit none
    integer :: count, k, operation, datatype, n, i, j, status
    real(dp) :: ar, ai
    complex(dp), allocatable :: a(:,:), x(:,:)
    real(dp), allocatable :: ra(:,:), rx(:,:)

    read(*,*) count
    do k = 1, count
        read(*,*) operation, datatype, n
        allocate(a(n,n), x(n,n))
        do i = 1, n
            do j = 1, n
                read(*,*) ar, ai
                a(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        if (datatype == 1) then
            allocate(ra(n,n),rx(n,n))
            ra = real(a,dp)
            select case (operation)
            case (1)
                call sinm(ra,rx,status)
            case (2)
                call cosm(ra,rx,status)
            case (3)
                call tanm(ra,rx,status)
            case (4)
                call sinhm(ra,rx,status)
            case (5)
                call coshm(ra,rx,status)
            case (6)
                call tanhm(ra,rx,status)
            end select
            x = cmplx(rx,kind=dp)
            deallocate(ra,rx)
        else
            select case (operation)
            case (1)
                call sinm(a,x,status)
            case (2)
                call cosm(a,x,status)
            case (3)
                call tanm(a,x,status)
            case (4)
                call sinhm(a,x,status)
            case (5)
                call coshm(a,x,status)
            case (6)
                call tanhm(a,x,status)
            end select
        end if
        write(*,'(i0)') status
        do i = 1, n
            do j = 1, n
                write(*,'(es29.17e3,1x,es29.17e3)') real(x(i,j),dp), aimag(x(i,j))
            end do
        end do
        deallocate(a,x)
    end do
end program driver_trigonometric
