! SPDX-License-Identifier: MIT
! Text-protocol driver for reproducible independent SciPy comparisons.
program driver_expm
    use scifort_kinds, only : dp
    use scifort_linalg_expm, only : expm
    implicit none
    real(dp), allocatable :: a(:,:), x(:,:)
    complex(dp), allocatable :: z(:,:), y(:,:)
    real(dp) :: real_part, imaginary_part
    integer :: cases, case_index, data_type, n, i, j, status, scalings

    read(*,*) cases
    do case_index=1,cases
        read(*,*) data_type, n
        allocate(a(n,n),x(n,n),z(n,n),y(n,n))
        do i=1,n
            do j=1,n
                read(*,*) real_part, imaginary_part
                a(i,j) = real_part
                z(i,j) = cmplx(real_part, imaginary_part, dp)
            end do
        end do
        scalings = 0
        select case (data_type)
        case (1)
            call expm(a,x,status,scalings)
        case (2)
            call expm(z,y,status,scalings)
        case default
            error stop 'unknown matrix data type'
        end select
        write(*,'(i0,1x,i0,1x,i0)') status, n, scalings
        do i=1,n
            do j=1,n
                if (data_type == 1) then
                    write(*,'(es24.16e3,1x,es24.16e3)') x(i,j), 0.0_dp
                else
                    write(*,'(es24.16e3,1x,es24.16e3)') real(y(i,j),dp), aimag(y(i,j))
                end if
            end do
        end do
        deallocate(a,x,z,y)
    end do
end program driver_expm
