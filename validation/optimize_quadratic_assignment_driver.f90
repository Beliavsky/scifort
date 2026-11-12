! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Batch driver for validation/compare_optimize_quadratic_assignment_scipy.py.
! Records:
! n method(1=faq,2=2opt) maximize nseed nguess maxiter hasp0 tol
! n matrix rows for A, n rows for B, nseed pairs, nguess pairs,
! optionally (n-nseed) rows of P0. Terminate with n=-1.
program optimize_quadratic_assignment_driver
    use scifort_kinds, only : dp
    use scifort_optimize_quadratic_assignment, only : quadratic_assignment, qap_result
    implicit none

    integer :: n, method, want_max, nseed, nguess, maxiter, hasp0, ios, i
    real(dp) :: tol
    real(dp), allocatable :: a(:, :), b(:, :), p0(:, :)
    integer, allocatable :: seed(:, :), guess(:, :)
    type(qap_result) :: res

    do
        read(*, *, iostat=ios) n, method, want_max, nseed, nguess, maxiter, hasp0, tol
        if (ios /= 0) exit
        if (n == -1) exit
        allocate(a(n, n), b(n, n), seed(nseed, 2), guess(nguess, 2), &
            p0(n - nseed, n - nseed))
        do i = 1, n
            read(*, *) a(i, :)
        end do
        do i = 1, n
            read(*, *) b(i, :)
        end do
        do i = 1, nseed
            read(*, *) seed(i, :)
        end do
        do i = 1, nguess
            read(*, *) guess(i, :)
        end do
        if (hasp0 /= 0) then
            do i = 1, n - nseed
                read(*, *) p0(i, :)
            end do
        end if
        select case (method)
        case (1)
            if (hasp0 /= 0) then
                res = quadratic_assignment(a, b, method='faq', maximize=want_max /= 0, &
                    partial_match=seed, p0=p0, maxiter=maxiter, tol=tol)
            else
                res = quadratic_assignment(a, b, method='faq', maximize=want_max /= 0, &
                    partial_match=seed, maxiter=maxiter, tol=tol)
            end if
        case (2)
            res = quadratic_assignment(a, b, method='2opt', maximize=want_max /= 0, &
                partial_match=seed, partial_guess=guess, maxiter=maxiter)
        end select
        write(*, '(a,1x,i0,1x,i0,1x,es26.17e3)') 'RESULT', res%status, res%nit, res%fun
        write(*, '(a)', advance='no') 'PERM'
        do i = 1, size(res%col_ind)
            write(*, '(1x,i0)', advance='no') res%col_ind(i)
        end do
        print *
        deallocate(a, b, seed, guess, p0)
    end do
end program optimize_quadratic_assignment_driver
