! SPDX-License-Identifier: MIT
! Streaming reference harness: svdvals, cond2 and matrix_rank.
program driver_svd_condition
    use scifort_kinds, only : dp
    use scifort_linalg_svd, only : svdvals, cond2, matrix_rank
    implicit none
    integer :: cases, kind_case, m, n, k, c, i, j
    integer :: st_s, st_c, st_r, st_rt, rank, rank_tol
    real(dp) :: re, im, tol, condition
    real(dp), allocatable :: ar(:, :), s(:)
    complex(dp), allocatable :: ac(:, :)

    read(*, *) cases
    do c = 1, cases
        read(*, *) kind_case, m, n, tol
        k = min(m,n)
        allocate(ar(m,n), ac(m,n), s(k))
        do j = 1, n
            do i = 1, m
                read(*, *) re, im
                ar(i,j) = re
                ac(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (kind_case == 1) then
            call svdvals(ar,s,st_s)
            call cond2(ar,condition,st_c)
            call matrix_rank(ar,rank,st_r)
            call matrix_rank(ar,rank_tol,st_rt,tol=tol)
        else
            call svdvals(ac,s,st_s)
            call cond2(ac,condition,st_c)
            call matrix_rank(ac,rank,st_r)
            call matrix_rank(ac,rank_tol,st_rt,tol=tol)
        end if
        write(*,'(6(I0,1X))') st_s,st_c,st_r,st_rt,rank,rank_tol
        write(*,'(ES25.16E3)') condition
        do i = 1, k
            write(*,'(ES25.16E3)') s(i)
        end do
        deallocate(ar,ac,s)
    end do
end program driver_svd_condition
