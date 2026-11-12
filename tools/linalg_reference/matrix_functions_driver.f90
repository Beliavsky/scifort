program matrix_functions_driver
    use scifort_kinds, only : dp
    use scifort_linalg_matrix_functions, only : sqrtm, signm
    implicit none
    integer :: cases, k, mode, n, i, j, status
    real(dp) :: r, im, res
    complex(dp), allocatable :: a(:, :), y(:, :)
    real(dp), allocatable :: ar(:, :), yr(:, :)
    read (*, *) cases
    do k = 1, cases
        read (*, *) mode, n
        allocate(a(n, n), y(n, n), ar(n, n), yr(n, n))
        do j = 1, n
            do i = 1, n
                read (*, *) r, im
                a(i, j) = cmplx(r, im, dp)
            end do
        end do
        res = -1.0_dp
        select case (mode)
        case (1)
            call sqrtm(a, y, status, res)
        case (2)
            call signm(a, y, status, residual=res)
        case (3)
            ar = real(a, dp)
            call sqrtm(ar, y, status, res)
        case (4)
            ar = real(a, dp)
            call signm(ar, yr, status, residual=res)
            y = cmplx(yr, kind=dp)
        case default
            error stop 'unknown test mode'
        end select
        write (*, '(i0,1x,es24.15)') status, res
        do j = 1, n
            do i = 1, n
                write (*, '(es26.17,1x,es26.17)') real(y(i, j), dp), aimag(y(i, j))
            end do
        end do
        deallocate(a, y, ar, yr)
    end do
end program matrix_functions_driver
