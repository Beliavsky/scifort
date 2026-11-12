! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_optimize_quadratic_assignment
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_quadratic_assignment, only : quadratic_assignment, &
        quadratic_assignment_score, qap_result, qap_success, qap_invalid_input
    implicit none

    real(dp) :: a4(4, 4), b4(4, 4), p_start(2, 2)
    real(dp), allocatable :: a(:, :), b(:, :)
    integer :: n, i, j, k, checks
    integer :: seeds(2, 2), known(4)
    logical :: maximize
    type(qap_result) :: r, fixed

    checks = 0
    a4(1, :) = [0.0_dp, 5.0_dp, 8.0_dp, 6.0_dp]
    a4(2, :) = [5.0_dp, 0.0_dp, 5.0_dp, 1.0_dp]
    a4(3, :) = [8.0_dp, 5.0_dp, 0.0_dp, 2.0_dp]
    a4(4, :) = [6.0_dp, 1.0_dp, 2.0_dp, 0.0_dp]
    b4(1, :) = [0.0_dp, 1.0_dp, 8.0_dp, 4.0_dp]
    b4(2, :) = [1.0_dp, 0.0_dp, 5.0_dp, 2.0_dp]
    b4(3, :) = [8.0_dp, 5.0_dp, 0.0_dp, 5.0_dp]
    b4(4, :) = [4.0_dp, 2.0_dp, 5.0_dp, 0.0_dp]

    ! SciPy documentation's 4x4 example: FAQ reaches 178; 2opt 176.
    r = quadratic_assignment(a4, b4, method='faq')
    call check(r%status == qap_success, 'example FAQ status')
    call check(abs(r%fun - 178.0_dp) < 1.0e-12_dp, 'example FAQ objective')
    call check(all(r%col_ind == [2, 1, 4, 3]), 'example FAQ permutation')
    call check(r%nit == 13, 'example FAQ iteration count')
    r = quadratic_assignment(a4, b4, method='2opt')
    call check(abs(r%fun - 176.0_dp) < 1.0e-12_dp, 'example 2opt objective')
    call check(all(r%col_ind == [2, 3, 4, 1]), 'example 2opt permutation')
    call check(r%nit > 0, 'example 2opt iterations')
    call check(abs(quadratic_assignment_score(a4, b4, r%col_ind) - r%fun) < 1.0e-12_dp, &
        'example 2opt independent score')

    seeds(1, :) = [4, 3]
    seeds(2, :) = [2, 4]
    do k = 1, 2
        maximize = k == 2
        r = quadratic_assignment(a4, b4, method='2opt', partial_match=seeds, maximize=maximize)
        call valid(r, a4, b4, 'seeded 2opt')
        call check(r%col_ind(4) == 3 .and. r%col_ind(2) == 4, 'fixed seeded pairs 2opt')
        r = quadratic_assignment(a4, b4, method='faq', partial_match=seeds, maximize=maximize)
        call valid(r, a4, b4, 'seeded FAQ')
        call check(r%col_ind(4) == 3 .and. r%col_ind(2) == 4, 'fixed seeded pairs FAQ')
    end do

    ! Complete fixed permutation must be returned without any optimization.
    known = [3, 1, 4, 2]
    fixed = quadratic_assignment(a4, b4, partial_match=reshape( &
        [1, 2, 3, 4, 3, 1, 4, 2], [4, 2]))
    call check(fixed%success .and. fixed%nit == 0, 'all pairs fixed')
    call check(all(fixed%col_ind == known), 'all fixed exact permutation')
    call check(abs(fixed%fun - quadratic_assignment_score(a4, b4, known)) < 1.0e-12_dp, &
        'all fixed exact score')

    ! Initial guessing works without fixing; maxiter=0 returns the completed guess.
    r = quadratic_assignment(a4, b4, method='2opt', maxiter=0, partial_guess=reshape( &
        [4, 1, 4, 3], [2, 2]))
    call check(r%success .and. r%nit == 0, 'no 2opt swaps')
    call check(all(r%col_ind == [3, 1, 2, 4]), 'partial guess completed')
    p_start(1, :) = [0.0_dp, 1.0_dp]
    p_start(2, :) = [1.0_dp, 0.0_dp]
    r = quadratic_assignment(a4, b4, partial_match=seeds, p0=p_start, maxiter=0)
    call check(r%success .and. r%nit == 0, 'no FAQ steps')
    call check(all(r%col_ind == [2, 4, 1, 3]), 'P0 free-block map')

    r = quadratic_assignment(a4, b4, method='FAQ')
    call check(r%success .and. abs(r%fun - 178.0_dp) < 1.0e-12_dp, &
        'case-insensitive method')
    r = quadratic_assignment(a4, b4, method='invalid')
    call check(r%status == qap_invalid_input, 'unknown method rejected')
    r = quadratic_assignment(a4, b4, method='faq', partial_guess=reshape([1, 1], [1, 2]))
    call check(r%status == qap_invalid_input, 'guess restricted to 2opt')
    r = quadratic_assignment(a4, b4, method='2opt', p0=p_start)
    call check(r%status == qap_invalid_input, 'P0 restricted to FAQ')
    r = quadratic_assignment(a4, b4, maxiter=-1)
    call check(r%status == qap_invalid_input, 'negative iterations rejected')
    r = quadratic_assignment(a4, b4, tol=-1.0_dp)
    call check(r%status == qap_invalid_input, 'negative tolerance rejected')
    r = quadratic_assignment(a4, b4, tol=ieee_value(0.0_dp, ieee_quiet_nan))
    call check(r%status == qap_invalid_input, 'NaN tolerance rejected')
    r = quadratic_assignment(a4(:, 1:3), b4)
    call check(r%status == qap_invalid_input, 'nonsquare A rejected')
    r = quadratic_assignment(a4, b4(1:3, 1:3))
    call check(r%status == qap_invalid_input, 'size mismatch rejected')
    r = quadratic_assignment(a4, b4, partial_match=reshape([1, 1, 2, 3], [2, 2]))
    call check(r%status == qap_invalid_input, 'duplicate seeded rows rejected')
    r = quadratic_assignment(a4, b4, partial_match=reshape([1, 2, 3, 3], [2, 2]))
    call check(r%status == qap_invalid_input, 'duplicate seeded columns rejected')
    r = quadratic_assignment(a4, b4, partial_match=reshape([1, 5], [1, 2]))
    call check(r%status == qap_invalid_input, 'seed out of bounds rejected')
    r = quadratic_assignment(a4, b4, partial_match=reshape([1, 2, 3], [1, 3]))
    call check(r%status == qap_invalid_input, 'seed shape rejected')
    r = quadratic_assignment(a4, b4, method='2opt', partial_guess=reshape( &
        [1, 2, 3, 3], [2, 2]))
    call check(r%status == qap_invalid_input, 'duplicate guessed columns rejected')
    r = quadratic_assignment(a4, b4, method='2opt', partial_guess=reshape( &
        [1, 1, 2, 3], [2, 2]))
    call check(r%status == qap_invalid_input, 'duplicate guessed rows rejected')
    r = quadratic_assignment(a4, b4, method='2opt', partial_match=reshape( &
        [1, 2], [1, 2]), partial_guess=reshape([1, 3], [1, 2]))
    call check(r%status == qap_invalid_input, 'conflicting guess and seed rejected')
    r = quadratic_assignment(a4, b4, p0=reshape([0.0_dp], [1, 1]))
    call check(r%status == qap_invalid_input, 'wrong-sized P0 rejected')
    r = quadratic_assignment(a4, b4, p0=spread( &
        [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], dim=2, ncopies=4))
    call check(r%status == qap_invalid_input, 'non-stochastic P0 rejected')
    r = quadratic_assignment(a4, b4, p0=reshape( &
        [0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 1.5_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [4, 4]))
    call check(r%status == qap_invalid_input, 'negative P0 entry rejected')
    call check(ieee_is_nan(quadratic_assignment_score(a4, b4, [1, 2, 2, 4])), &
        'duplicate score permutation rejected')
    call check(ieee_is_nan(quadratic_assignment_score(a4, b4, [1, 2, 3])), &
        'incorrect score permutation size rejected')

    block
        real(dp) :: zero(0, 0), small(1, 1)
        r = quadratic_assignment(zero, zero)
        call check(r%success .and. size(r%col_ind) == 0, 'empty problem')
        call check(abs(r%fun) < 1.0e-14_dp, 'empty score')
        small(1, 1) = -3.0_dp
        r = quadratic_assignment(small, small)
        call check(r%success .and. all(r%col_ind == [1]), 'single-node problem')
        call check(abs(r%fun - 9.0_dp) < 1.0e-14_dp, 'single-node score')
    end block

    block
        real(dp) :: bad(4, 4)
        bad = a4
        bad(2, 3) = ieee_value(0.0_dp, ieee_quiet_nan)
        r = quadratic_assignment(bad, b4)
        call check(r%status == qap_invalid_input, 'NaN matrix rejected')
    end block

    ! Exhaustive local-optimality checks use direct objective evaluations,
    ! independent of 2opt's O(n) analytical swap-delta formula.
    do n = 2, 10
        allocate(a(n, n), b(n, n))
        do j = 1, n
            do i = 1, n
                a(i, j) = real(mod(23 * i + 17 * j + 7 * i * j + 9 * i * i, 41) - 20, dp)
                b(i, j) = real(mod(11 * i + 13 * j + 5 * i * i * j, 43) - 21, dp)
            end do
        end do
        do k = 1, 2
            maximize = k == 2
            r = quadratic_assignment(a, b, method='2opt', maximize=maximize)
            call valid(r, a, b, '2opt generated')
            call check(r%nit >= 0 .and. r%nit <= 10000, '2opt accepted swap count')
            call check(local_optimal(a, b, r%col_ind, maximize), '2opt is pairwise locally optimal')
            fixed = quadratic_assignment(a, b, method='faq', maximize=maximize)
            call valid(fixed, a, b, 'FAQ generated')
            call check(fixed%nit >= 0 .and. fixed%nit <= 30, 'FAQ step count')
            call check(score_better_or_equal(r%fun, &
                quadratic_assignment_score(a, b, [(i, i=1,n)]), maximize), &
                '2opt no worse than identity')
        end do
        deallocate(a, b)
    end do
    print '(a,i0,a)', 'quadratic_assignment: ', checks, ' native assertions passed'

contains

    subroutine check(ok, description)
        logical, intent(in) :: ok !! whether assertion succeeds
        character(len=*), intent(in) :: description !! failed-case diagnostic
        if (.not. ok) then
            print '(a,a)', 'FAIL: ', description
            error stop 1
        end if
        checks = checks + 1
    end subroutine check

    subroutine valid(res, aa, bb, label)
        type(qap_result), intent(in) :: res !! solver output to inspect
        real(dp), intent(in) :: aa(:, :) !! original A matrix
        real(dp), intent(in) :: bb(:, :) !! original B matrix
        character(len=*), intent(in) :: label !! assertion label
        integer :: t
        call check(res%success .and. res%status == qap_success, label//' success')
        call check(size(res%col_ind) == size(aa, 1), label//' permutation size')
        do t = 1, size(res%col_ind)
            call check(count(res%col_ind == t) == 1, label//' permutation multiplicity')
        end do
        call check(abs(res%fun - quadratic_assignment_score(aa, bb, res%col_ind)) < &
            1.0e-10_dp * max(1.0_dp, abs(res%fun)), label//' objective')
    end subroutine valid

    pure function score_better_or_equal(x, y, want_max) result(ok)
        real(dp), intent(in) :: x !! candidate objective
        real(dp), intent(in) :: y !! comparison objective
        logical, intent(in) :: want_max !! maximize or minimize
        logical :: ok !! true if x is at least as good as y
        ok = x <= y + 1.0e-9_dp
        if (want_max) ok = x >= y - 1.0e-9_dp
    end function score_better_or_equal

    function local_optimal(aa, bb, perm, want_max) result(ok)
        real(dp), intent(in) :: aa(:, :) !! first matrix
        real(dp), intent(in) :: bb(:, :) !! second matrix
        integer, intent(in) :: perm(:) !! candidate permutation
        logical, intent(in) :: want_max !! maximize or minimize
        logical :: ok !! no pair swap improves the score
        integer :: i, j, t, trial(size(perm))
        real(dp) :: f0, f1
        ok = .true.
        f0 = quadratic_assignment_score(aa, bb, perm)
        do i = 1, size(perm) - 1
            do j = i + 1, size(perm)
                trial = perm
                t = trial(i)
                trial(i) = trial(j)
                trial(j) = t
                f1 = quadratic_assignment_score(aa, bb, trial)
                if (.not. score_better_or_equal(f0, f1, want_max)) ok = .false.
            end do
        end do
    end function local_optimal

end program test_optimize_quadratic_assignment
