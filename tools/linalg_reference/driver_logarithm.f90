! SPDX-License-Identifier: MIT
! Streamed real/complex matrix-logarithm reference driver.
program driver_logarithm
    use scifort_kinds, only : dp
    use scifort_linalg_logarithm, only : logm, fractional_matrix_power
    implicit none
    integer :: cases, k, mode, n, i, j, status
    real(dp) :: r, im, t, err
    complex(dp), allocatable :: a(:, :), x(:, :)
    real(dp), allocatable :: ar(:, :)

    read (*, *) cases
    do k = 1, cases
        read (*, *) mode, n, t
        allocate(a(n,n), x(n,n), ar(n,n))
        do j = 1, n
            do i = 1, n
                read (*, *) r, im
                a(i,j) = cmplx(r, im, dp)
            end do
        end do
        ar = real(a, dp)
        err = -1.0_dp
        select case(mode)
        case(1)
            call logm(ar, x, status, err)
        case(2)
            call logm(a, x, status, err)
        case(3)
            call fractional_matrix_power(ar, t, x, status)
        case(4)
            call fractional_matrix_power(a, t, x, status)
        case default
            error stop 'unknown logarithm mode'
        end select
        write(*,'(i0,1x,es24.15)') status, err
        do j = 1, n
            do i = 1, n
                write(*,'(es29.17e3,1x,es29.17e3)') real(x(i,j),dp), aimag(x(i,j))
            end do
        end do
        deallocate(a,x,ar)
    end do
end program driver_logarithm
