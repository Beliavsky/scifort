! SPDX-License-Identifier: MIT
! Demonstrate SciFort's values-only SVD and spectral condition/rank estimates.
program linalg_svd_condition_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_svd, only : svdvals, cond2, matrix_rank
    implicit none
    real(dp) :: a(3,3), s(3), condition
    integer :: status, rank

    a = 0.0_dp
    a(1,1) = 4.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 0.25_dp
    call svdvals(a,s,status)
    if (status /= linalg_status_success) error stop 'svdvals failed'
    call cond2(a,condition,status)
    if (status /= linalg_status_success) error stop 'cond2 failed'
    call matrix_rank(a,rank,status)
    if (status /= linalg_status_success) error stop 'matrix_rank failed'

    print '(A,3(F8.3,1X))', 'singular values: ', s
    print '(A,F8.3)', '2-norm condition: ', condition
    print '(A,I0)', 'numerical rank: ', rank
end program linalg_svd_condition_example
