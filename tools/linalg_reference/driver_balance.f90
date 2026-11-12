! SPDX-License-Identifier: MIT
! Stream packed real/complex matrices and similarity metadata for reference checks.
program driver_balance
    use scifort_kinds, only : dp
    use scifort_linalg_balance, only : matrix_balance
    implicit none
    integer :: cases, k, n, i, j, kind_case, perm_flag, scale_flag, status
    real(dp) :: ar, ai
    real(dp), allocatable :: a(:, :), b(:, :), factors(:)
    complex(dp), allocatable :: ac(:, :), bc(:, :)
    integer, allocatable :: perm(:)

    read(*,*) cases
    do k = 1, cases
        read(*,*) kind_case, n, perm_flag, scale_flag
        allocate(a(n,n),b(n,n),ac(n,n),bc(n,n),factors(n),perm(n))
        do j = 1, n
            do i = 1, n
                read(*,*) ar, ai
                a(i,j) = ar
                ac(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        if (kind_case == 1) then
            call matrix_balance(a,b,status,permute=perm_flag /= 0,scale=scale_flag /= 0, &
                scale_factors=factors,permutation=perm)
            bc = cmplx(b,0.0_dp,dp)
        else if (kind_case == 2) then
            call matrix_balance(ac,bc,status,permute=perm_flag /= 0,scale=scale_flag /= 0, &
                scale_factors=factors,permutation=perm)
        else
            error stop 'unsupported matrix kind'
        end if
        write(*,'(I0)') status
        do j = 1, n
            do i = 1, n
                write(*,'(ES25.16E3,1X,ES25.16E3)') real(bc(i,j),dp),aimag(bc(i,j))
            end do
        end do
        do i = 1, n
            write(*,'(ES25.16E3,1X,I0)') factors(i),perm(i)
        end do
        deallocate(a,b,ac,bc,factors,perm)
    end do
end program driver_balance
