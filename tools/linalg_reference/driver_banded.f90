! SPDX-License-Identifier: MIT
! Independent input/output protocol for scipy.linalg banded comparisons.
program driver_banded
    use scifort_kinds, only : dp
    use scifort_linalg_banded, only : solve_banded, solveh_banded, cholesky_banded, cho_solve_banded
    implicit none
    integer :: cases, case_number, operation, data_kind, lower, upper, n, nrhs, lower_flag
    integer :: i, j, status
    logical :: use_lower
    real(dp) :: ar, ai
    real(dp), allocatable :: real_ab(:,:), real_b(:,:), real_x(:,:), real_factor(:,:)
    complex(dp), allocatable :: complex_ab(:,:), complex_b(:,:), complex_x(:,:), complex_factor(:,:)

    read(*,*) cases
    do case_number = 1, cases
        read(*,*) operation, data_kind, lower, upper, n, nrhs, lower_flag
        use_lower = lower_flag == 1
        allocate(complex_ab(lower+upper+1,n),complex_b(n,nrhs),complex_x(n,nrhs))
        allocate(real_ab(lower+upper+1,n),real_b(n,nrhs),real_x(n,nrhs))
        allocate(complex_factor(lower+upper+1,n),real_factor(lower+upper+1,n))
        do j=1,n
            do i=1,lower+upper+1
                read(*,*) ar,ai
                complex_ab(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        do j=1,nrhs
            do i=1,n
                read(*,*) ar,ai
                complex_b(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        real_ab = real(complex_ab,dp)
        real_b = real(complex_b,dp)
        if (data_kind == 1) then
            select case (operation)
            case (1)
                call solve_banded(real_ab,lower,upper,real_b,real_x,status)
            case (2)
                call solveh_banded(real_ab,real_b,real_x,status,use_lower)
            case (3)
                call cholesky_banded(real_ab,real_factor,status,use_lower)
                if (status == 0) call cho_solve_banded(real_factor,real_b,real_x,status,use_lower)
            end select
            complex_x = cmplx(real_x,kind=dp)
        else
            select case (operation)
            case (1)
                call solve_banded(complex_ab,lower,upper,complex_b,complex_x,status)
            case (2)
                call solveh_banded(complex_ab,complex_b,complex_x,status,use_lower)
            case (3)
                call cholesky_banded(complex_ab,complex_factor,status,use_lower)
                if (status == 0) call cho_solve_banded(complex_factor,complex_b,complex_x,status,use_lower)
            end select
        end if
        write(*,'(I0)') status
        if (status == 0) then
            do j=1,nrhs
                do i=1,n
                    write(*,'(ES26.17E3,1X,ES26.17E3)') real(complex_x(i,j),dp), aimag(complex_x(i,j))
                end do
            end do
        end if
        deallocate(complex_ab,complex_b,complex_x,complex_factor)
        deallocate(real_ab,real_b,real_x,real_factor)
    end do
end program driver_banded
