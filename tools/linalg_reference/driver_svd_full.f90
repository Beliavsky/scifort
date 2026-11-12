! SPDX-License-Identifier: MIT
! Streaming full-size SVD and singular-value-only SciPy comparison driver.
program driver_svd_full
    use scifort_kinds, only : dp
    use scifort_linalg_svd, only : svd, svdvals, diagsvd
    implicit none
    complex(dp), allocatable :: a(:, :), u(:, :), vh(:, :)
    real(dp), allocatable :: s(:), ss(:), d(:, :)
    real(dp) :: re, im
    integer :: trials, testcase, type_id, m, n, k, i, j, st, stvals, stdiag

    read(*,*) trials
    do testcase = 1, trials
        read(*,*) type_id, m, n
        k = min(m,n)
        allocate(a(m,n),u(m,m),vh(n,n),s(k),ss(k),d(m,n))
        do j = 1, n
            do i = 1, m
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (type_id == 1) then
            block
                real(dp), allocatable :: ar(:, :), ur(:, :), vhr(:, :)
                allocate(ar(m,n),ur(m,m),vhr(n,n))
                ar = real(a,dp)
                call svd(ar,ur,s,vhr,st,full_matrices=.true.)
                call svdvals(ar,ss,stvals)
                u = cmplx(ur,0.0_dp,dp)
                vh = cmplx(vhr,0.0_dp,dp)
            end block
        else
            call svd(a,u,s,vh,st,full_matrices=.true.)
            call svdvals(a,ss,stvals)
        end if
        call diagsvd(s,m,n,d,stdiag)
        print '(3(I0,1X))', st,stvals,stdiag
        do i = 1, k
            print '(ES26.17E3)', s(i)
        end do
        do i = 1, k
            print '(ES26.17E3)', ss(i)
        end do
        do j = 1, m
            do i = 1, m
                call print_complex(u(i,j))
            end do
        end do
        do j = 1, n
            do i = 1, n
                call print_complex(vh(i,j))
            end do
        end do
        do j = 1, n
            do i = 1, m
                print '(ES26.17E3)', d(i,j)
            end do
        end do
        deallocate(a,u,vh,s,ss,d)
    end do
contains
    subroutine print_complex(z)
        complex(dp), intent(in) :: z !! finite value to serialize
        print '(ES26.17E3,1X,ES26.17E3)', real(z,dp),aimag(z)
    end subroutine print_complex
end program driver_svd_full
