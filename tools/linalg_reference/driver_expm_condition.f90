! SPDX-License-Identifier: MIT
! Driver for independent Python/SciPy expm_cond comparison.
program driver_expm_condition
    use scifort_kinds, only : dp
    use scifort_linalg_expm_condition, only : expm_cond
    implicit none
    real(dp), allocatable :: r(:, :)
    complex(dp), allocatable :: a(:, :)
    real(dp) :: re, im, condition
    integer :: cases, dtype, n, i, j, k, status, count

    read(*,*) cases
    do k = 1, cases
        read(*,*) dtype, n
        allocate(a(n,n),r(n,n))
        do i = 1, n
            do j = 1, n
                read(*,*) re, im
                a(i,j) = cmplx(re, im, dp)
            end do
        end do
        if (dtype == 1) then
            r = real(a, dp)
            call expm_cond(r, condition, status, iterations=count)
        else
            call expm_cond(a, condition, status, iterations=count)
        end if
        write(*,'(i0,1x,i0,1x,es29.17e3)') status, count, condition
        deallocate(a,r)
    end do
end program driver_expm_condition
