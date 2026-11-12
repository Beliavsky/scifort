! SPDX-License-Identifier: MIT
! Streamed numerical reference driver for the SciPy comparison harness.
program driver_matrix_equations
    use scifort_kinds, only : dp
    use scifort_linalg_matrix_equations, only : solve_sylvester, &
        solve_continuous_lyapunov, solve_discrete_lyapunov
    implicit none
    real(dp), allocatable :: ar(:, :), br(:, :), qr(:, :), xr(:, :)
    complex(dp), allocatable :: ac(:, :), bc(:, :), qc(:, :), xc(:, :)
    integer :: code, m, n, ios, status, i, j
    character(len=1) :: kind

    do
        read(*, *, iostat=ios) code, kind, m, n
        if (ios < 0) exit
        if (ios /= 0 .or. m <= 0 .or. n <= 0) error stop 'invalid reference header'
        if (code /= 1) then
            if (m /= n) error stop 'Lyapunov requires square input'
        end if
        if (kind == 'R') then
            allocate(ar(m, m), br(n, n), qr(m, n), xr(m, n))
            read(*, *, iostat=ios) ar
            if (ios /= 0) error stop 'read A real'
            if (code == 1) then
                read(*, *, iostat=ios) br
                if (ios /= 0) error stop 'read B real'
            end if
            read(*, *, iostat=ios) qr
            if (ios /= 0) error stop 'read Q real'
            select case (code)
            case (1)
                call solve_sylvester(ar, br, qr, xr, status)
            case (2)
                call solve_continuous_lyapunov(ar, qr, xr, status)
            case (3)
                call solve_discrete_lyapunov(ar, qr, xr, status)
            case default
                error stop 'invalid equation code'
            end select
            write(*, '(i0)') status
            do j = 1, n
                do i = 1, m
                    write(*, '(es26.17e3)') xr(i, j)
                end do
            end do
            deallocate(ar, br, qr, xr)
        else if (kind == 'C') then
            allocate(ac(m, m), bc(n, n), qc(m, n), xc(m, n))
            read(*, *, iostat=ios) ac
            if (ios /= 0) error stop 'read A complex'
            if (code == 1) then
                read(*, *, iostat=ios) bc
                if (ios /= 0) error stop 'read B complex'
            end if
            read(*, *, iostat=ios) qc
            if (ios /= 0) error stop 'read Q complex'
            select case (code)
            case (1)
                call solve_sylvester(ac, bc, qc, xc, status)
            case (2)
                call solve_continuous_lyapunov(ac, qc, xc, status)
            case (3)
                call solve_discrete_lyapunov(ac, qc, xc, status)
            case default
                error stop 'invalid equation code'
            end select
            write(*, '(i0)') status
            do j = 1, n
                do i = 1, m
                    write(*, '(es26.17e3,1x,es26.17e3)') real(xc(i, j), dp), aimag(xc(i, j))
                end do
            end do
            deallocate(ac, bc, qc, xc)
        else
            error stop 'invalid real or complex flag'
        end if
    end do
end program driver_matrix_equations
