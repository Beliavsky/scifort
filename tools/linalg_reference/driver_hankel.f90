! SPDX-License-Identifier: MIT
! Independent comparison driver for SciFort hankel and matmul_hankel.
! Input: n m nrhs kind explicit_row emit_dense; complex c, r, x values.
! kind: 0=real/real; 1=real/complex; 2=complex/real; 3=complex/complex.
program driver_hankel
    use scifort_kinds, only : dp
    use scifort_linalg_hankel, only : hankel
    use scifort_linalg_matmul_hankel, only : matmul_hankel
    implicit none
    complex(dp), allocatable :: c(:)
    complex(dp), allocatable :: r(:)
    complex(dp), allocatable :: x(:,:)
    complex(dp), allocatable :: y(:,:)
    complex(dp), allocatable :: h(:,:)
    real(dp) :: real_part, imaginary_part
    integer :: n, m, nrhs, kind, explicit_row, emit_dense
    integer :: status, ios, i, j

    do
        read(*,*,iostat=ios) n,m,nrhs,kind,explicit_row,emit_dense
        if (ios /= 0 .or. n < 0) exit
        allocate(c(n),r(m),x(m,nrhs),y(n,nrhs),h(n,m))
        do i=1,n
            read(*,*) real_part, imaginary_part
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
        call compute_product()
        if (status == 0 .and. emit_dense /= 0) call compute_dense()
        print '(i0)', status
        if (status == 0) then
            do j=1,nrhs
                do i=1,n
                    call print_complex(y(i,j))
                end do
            end do
            if (emit_dense /= 0) then
                do j=1,m
                    do i=1,n
                        call print_complex(h(i,j))
                    end do
                end do
            end if
        end if
        deallocate(c,r,x,y,h)
    end do

contains

    subroutine print_complex(z)
        complex(dp), intent(in) :: z !! value to print
        write(*,'(es27.17e3,1x,es27.17e3)') real(z,dp),aimag(z)
    end subroutine print_complex

    subroutine compute_product()
        real(dp), allocatable :: yr(:,:)
        integer :: s
        select case(kind)
        case(0)
            allocate(yr(n,nrhs))
            if (explicit_row /= 0) then
                call matmul_hankel(real(c,dp),real(x,dp),yr,s,real(r,dp))
            else
                call matmul_hankel(real(c,dp),real(x,dp),yr,s)
            end if
            y = cmplx(yr,kind=dp)
        case(1)
            if (explicit_row /= 0) then
                call matmul_hankel(real(c,dp),x,y,s,real(r,dp))
            else
                call matmul_hankel(real(c,dp),x,y,s)
            end if
        case(2)
            if (explicit_row /= 0) then
                call matmul_hankel(c,real(x,dp),y,s,r)
            else
                call matmul_hankel(c,real(x,dp),y,s)
            end if
        case(3)
            if (explicit_row /= 0) then
                call matmul_hankel(c,x,y,s,r)
            else
                call matmul_hankel(c,x,y,s)
            end if
        case default
            s = 1
            y = (0.0_dp,0.0_dp)
        end select
        status = s
    end subroutine compute_product

    subroutine compute_dense()
        real(dp), allocatable :: hr(:,:)
        integer :: s
        if (kind < 2) then
            allocate(hr(n,m))
            if (explicit_row /= 0) then
                call hankel(real(c,dp),hr,s,real(r,dp))
            else
                call hankel(real(c,dp),hr,s)
            end if
            h = cmplx(hr,kind=dp)
        else
            if (explicit_row /= 0) then
                call hankel(c,h,s,r)
            else
                call hankel(c,h,s)
            end if
        end if
        status = s
    end subroutine compute_dense

end program driver_hankel
