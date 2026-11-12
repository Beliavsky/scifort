! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reproduce the published SciPy 4x4 QAP example, then refine FAQ with 2-opt.
program optimize_quadratic_assignment_example
    use scifort_kinds, only : dp
    use scifort_optimize_quadratic_assignment, only : quadratic_assignment, qap_result
    implicit none
    real(dp) :: a(4, 4), b(4, 4)
    type(qap_result) :: faq, refined, seeded
    integer :: i, guess(4, 2)

    a(1, :) = [0.0_dp, 5.0_dp, 8.0_dp, 6.0_dp]
    a(2, :) = [5.0_dp, 0.0_dp, 5.0_dp, 1.0_dp]
    a(3, :) = [8.0_dp, 5.0_dp, 0.0_dp, 2.0_dp]
    a(4, :) = [6.0_dp, 1.0_dp, 2.0_dp, 0.0_dp]
    b(1, :) = [0.0_dp, 1.0_dp, 8.0_dp, 4.0_dp]
    b(2, :) = [1.0_dp, 0.0_dp, 5.0_dp, 2.0_dp]
    b(3, :) = [8.0_dp, 5.0_dp, 0.0_dp, 5.0_dp]
    b(4, :) = [4.0_dp, 2.0_dp, 5.0_dp, 0.0_dp]

    faq = quadratic_assignment(a, b, method='faq')
    if (.not. faq%success) error stop 'FAQ failed'
    do i = 1, 4
        guess(i, :) = [i, faq%col_ind(i)]
    end do
    refined = quadratic_assignment(a, b, method='2opt', partial_guess=guess)
    if (.not. refined%success) error stop '2opt failed'
    print '(a,f10.2,a,4(1x,i0))', 'FAQ objective: ', faq%fun, ' permutation:', faq%col_ind
    print '(a,f10.2,a,4(1x,i0))', '2opt objective:', refined%fun, ' permutation:', refined%col_ind
    if (abs(faq%fun - 178.0_dp) > 1.0e-10_dp) error stop 'unexpected FAQ score'
    if (abs(refined%fun - 176.0_dp) > 1.0e-10_dp) error stop 'unexpected 2opt score'

    ! Fixed A-node 1 -> B-node 2: the optimizer may change only other rows.
    seeded = quadratic_assignment(a, b, method='2opt', &
        partial_match=reshape([1, 2], [1, 2]))
    if (.not. seeded%success .or. seeded%col_ind(1) /= 2) error stop 'seeded match failed'
    print '(a,4(1x,i0))', 'seeded permutation:', seeded%col_ind
end program optimize_quadratic_assignment_example
