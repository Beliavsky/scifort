! SPDX-License-Identifier: MIT
! Regression of the common elementwise batching engine at every rank through 15.
program test_optimize_elementwise_highrank
    use scifort_kinds, only : dp
    use scifort_optimize_root_scalar, only : bracket_root, bracket_root_result, &
        find_root, find_root_result, root_status_invalid_input
    use scifort_optimize_minimize_scalar, only : bracket_minimum, bracket_minimum_result, &
        find_minimum, find_minimum_result, minimize_status_invalid_input
    implicit none
    integer :: calls, failures, callback_calls

    failures = 0
    callback_calls = 0
    ! Rank 1: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2)
        type(bracket_root_result) :: root_b(2)
        type(bracket_minimum_result) :: min_b(2)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 1 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 1 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 1 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 1 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 1 invalid broadcast is rejected')
        call check(calls == 0, 'rank 1 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2), min_triples(3,2)
        type(find_root_result) :: roots(2)
        type(find_minimum_result) :: mins(2)
        root_pairs = 0.0_dp
        root_pairs(2,:) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:) = 1.0_dp
        min_triples(3,:) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 1 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 1 root values')
        call check(calls == sum(roots%nfev), 'rank 1 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 1 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 1 minimum values')
        call check(calls == sum(mins%nfev), 'rank 1 minimum refinement evaluations')
    end block
    ! Rank 2: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1)
        type(bracket_root_result) :: root_b(2,1)
        type(bracket_minimum_result) :: min_b(2,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 2 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 2 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 2 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 2 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 2 invalid broadcast is rejected')
        call check(calls == 0, 'rank 2 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1), min_triples(3,2,1)
        type(find_root_result) :: roots(2,1)
        type(find_minimum_result) :: mins(2,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :) = 1.0_dp
        min_triples(3,:, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 2 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 2 root values')
        call check(calls == sum(roots%nfev), 'rank 2 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 2 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 2 minimum values')
        call check(calls == sum(mins%nfev), 'rank 2 minimum refinement evaluations')
    end block
    ! Rank 3: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1)
        type(bracket_root_result) :: root_b(2,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 3 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 3 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 3 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 3 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 3 invalid broadcast is rejected')
        call check(calls == 0, 'rank 3 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1), min_triples(3,2,1,1)
        type(find_root_result) :: roots(2,1,1)
        type(find_minimum_result) :: mins(2,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :) = 1.0_dp
        min_triples(3,:, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 3 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 3 root values')
        call check(calls == sum(roots%nfev), 'rank 3 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 3 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 3 minimum values')
        call check(calls == sum(mins%nfev), 'rank 3 minimum refinement evaluations')
    end block
    ! Rank 4: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 4 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 4 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 4 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 4 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 4 invalid broadcast is rejected')
        call check(calls == 0, 'rank 4 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1), min_triples(3,2,1,1,1)
        type(find_root_result) :: roots(2,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 4 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 4 root values')
        call check(calls == sum(roots%nfev), 'rank 4 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 4 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 4 minimum values')
        call check(calls == sum(mins%nfev), 'rank 4 minimum refinement evaluations')
    end block
    ! Rank 5: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 5 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 5 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 5 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 5 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 5 invalid broadcast is rejected')
        call check(calls == 0, 'rank 5 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1), min_triples(3,2,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 5 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 5 root values')
        call check(calls == sum(roots%nfev), 'rank 5 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 5 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 5 minimum values')
        call check(calls == sum(mins%nfev), 'rank 5 minimum refinement evaluations')
    end block
    ! Rank 6: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 6 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 6 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 6 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 6 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 6 invalid broadcast is rejected')
        call check(calls == 0, 'rank 6 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1), min_triples(3,2,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 6 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 6 root values')
        call check(calls == sum(roots%nfev), 'rank 6 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 6 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 6 minimum values')
        call check(calls == sum(mins%nfev), 'rank 6 minimum refinement evaluations')
    end block
    ! Rank 7: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 7 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 7 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 7 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 7 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 7 invalid broadcast is rejected')
        call check(calls == 0, 'rank 7 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 7 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 7 root values')
        call check(calls == sum(roots%nfev), 'rank 7 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 7 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 7 minimum values')
        call check(calls == sum(mins%nfev), 'rank 7 minimum refinement evaluations')
    end block
    ! Rank 8: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 8 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 8 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 8 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 8 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 8 invalid broadcast is rejected')
        call check(calls == 0, 'rank 8 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 8 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 8 root values')
        call check(calls == sum(roots%nfev), 'rank 8 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 8 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 8 minimum values')
        call check(calls == sum(mins%nfev), 'rank 8 minimum refinement evaluations')
    end block
    ! Rank 9: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 9 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 9 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 9 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 9 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 9 invalid broadcast is rejected')
        call check(calls == 0, 'rank 9 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 9 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 9 root values')
        call check(calls == sum(roots%nfev), 'rank 9 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 9 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 9 minimum values')
        call check(calls == sum(mins%nfev), 'rank 9 minimum refinement evaluations')
    end block
    ! Rank 10: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 10 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 10 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 10 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 10 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 10 invalid broadcast is rejected')
        call check(calls == 0, 'rank 10 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 10 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 10 root values')
        call check(calls == sum(roots%nfev), 'rank 10 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 10 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 10 minimum values')
        call check(calls == sum(mins%nfev), 'rank 10 minimum refinement evaluations')
    end block
    ! Rank 11: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 11 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 11 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 11 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 11 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 11 invalid broadcast is rejected')
        call check(calls == 0, 'rank 11 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 11 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 11 root values')
        call check(calls == sum(roots%nfev), 'rank 11 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 11 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 11 minimum values')
        call check(calls == sum(mins%nfev), 'rank 11 minimum refinement evaluations')
    end block
    ! Rank 12: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 12 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 12 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 12 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 12 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 12 invalid broadcast is rejected')
        call check(calls == 0, 'rank 12 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 12 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 12 root values')
        call check(calls == sum(roots%nfev), 'rank 12 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 12 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 12 minimum values')
        call check(calls == sum(mins%nfev), 'rank 12 minimum refinement evaluations')
    end block
    ! Rank 13: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 13 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 13 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 13 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 13 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 13 invalid broadcast is rejected')
        call check(calls == 0, 'rank 13 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 13 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 13 root values')
        call check(calls == sum(roots%nfev), 'rank 13 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 13 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 13 minimum values')
        call check(calls == sum(mins%nfev), 'rank 13 minimum refinement evaluations')
    end block
    ! Rank 14: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 14 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 14 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 14 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 14 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 14 invalid broadcast is rejected')
        call check(calls == 0, 'rank 14 invalid broadcast does not call objective')
    end block
    block
        real(dp) :: root_pairs(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1), min_triples(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(find_minimum_result) :: mins(2,1,1,1,1,1,1,1,1,1,1,1,1,1)
        root_pairs = 0.0_dp
        root_pairs(2,:, :, :, :, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        min_triples = 0.0_dp
        min_triples(2,:, :, :, :, :, :, :, :, :, :, :, :, :, :) = 1.0_dp
        min_triples(3,:, :, :, :, :, :, :, :, :, :, :, :, :, :) = 3.0_dp
        calls = 0
        roots = find_root(linear_function, root_pairs, xatol=1.0e-12_dp)
        call check(all(roots%success), 'rank 14 root refinement')
        call check(all(abs(roots%x - 2.0_dp) < 1.0e-9_dp), &
            'rank 14 root values')
        call check(calls == sum(roots%nfev), 'rank 14 root refinement evaluations')
        calls = 0
        mins = find_minimum(quadratic_function, min_triples, xatol=1.0e-9_dp)
        call check(all(mins%success), 'rank 14 minimum refinement')
        call check(all(abs(mins%x - 1.0_dp) < 1.0e-8_dp), &
            'rank 14 minimum values')
        call check(calls == sum(mins%nfev), 'rank 14 minimum refinement evaluations')
    end block
    ! Rank 15: both bracket searches, and refinements through rank 14.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: root_b(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_minimum_result) :: min_b(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp], &
            xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(root_b%success), 'rank 15 root brackets')
        call check(calls == sum(root_b%nfev), 'rank 15 root bracket evaluations')
        calls = 0
        min_b = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
        call check(all(min_b%success), 'rank 15 minimum brackets')
        call check(calls == sum(min_b%nfev), 'rank 15 minimum bracket evaluations')
        calls = 0
        root_b = bracket_root(linear_function, starts, xr0=[1.0_dp, 2.0_dp, 3.0_dp])
        call check(all(root_b%status == root_status_invalid_input), &
            'rank 15 invalid broadcast is rejected')
        call check(calls == 0, 'rank 15 invalid broadcast does not call objective')
    end block
    ! Rank-15 limits supplied as a rank-14 trailing-axis array. A failure
    ! in one slab must not change the statuses or evaluation counts elsewhere.
    block
        real(dp) :: starts(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
        real(dp) :: limits(1,1,1,1,1,1,1,1,1,1,1,1,1,2)
        real(dp) :: full_limits(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
        type(bracket_root_result) :: roots(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2), r0
        type(bracket_minimum_result) :: minima(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2), m0
        integer :: i, j, batch_calls
        starts = 0.0_dp
        limits = 10.0_dp
        limits(1,1,1,1,1,1,1,1,1,1,1,1,1,2) = 0.5_dp
        calls = 0
        roots = bracket_root(linear_function, starts, xr0=1.0_dp, &
            xmin=-10.0_dp, xmax=limits)
        batch_calls = calls
        call check(batch_calls == sum(roots%nfev), &
            'rank 15 trailing-axis root accounting')
        calls = 0
        minima = bracket_minimum(quadratic_function, starts, &
            xl0=-1.0_dp, xr0=1.0_dp, xmin=-10.0_dp, xmax=limits)
        batch_calls = calls
        call check(batch_calls == sum(minima%nfev), &
            'rank 15 trailing-axis minimum accounting')
        do j = 1, 2
            do i = 1, 2
                r0 = bracket_root(linear_function, 0.0_dp, xr0=1.0_dp, &
                    xmin=-10.0_dp, xmax=limits(1,1,1,1,1,1,1,1,1,1,1,1,1,j))
                m0 = bracket_minimum(quadratic_function, 0.0_dp, &
                    xl0=-1.0_dp, xr0=1.0_dp, &
                    xmin=-10.0_dp, xmax=limits(1,1,1,1,1,1,1,1,1,1,1,1,1,j))
                call check(roots(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j)%status == r0%status, &
                    'rank 15 trailing broadcast root status')
                call check(roots(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j)%nfev == r0%nfev, &
                    'rank 15 trailing broadcast root evaluations')
                call check(minima(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j)%status == m0%status, &
                    'rank 15 trailing broadcast minimum status')
                call check(minima(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j)%nfev == m0%nfev, &
                    'rank 15 trailing broadcast minimum evaluations')
            end do
        end do
        ! Rank-15 optional argument, not just a rank-14 one: this exercises
        ! the highest select-rank branch and a nontrailing varying axis.
        full_limits = 10.0_dp
        full_limits(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1) = 0.5_dp
        calls = 0
        roots = bracket_root(linear_function, starts, xr0=1.0_dp, &
            xmin=-10.0_dp, xmax=full_limits)
        call check(calls == sum(roots%nfev), &
            'rank 15 full-shape argument evaluation counts')
        do j = 1, 2
            do i = 1, 2
                r0 = bracket_root(linear_function, 0.0_dp, xr0=1.0_dp, &
                    xmin=-10.0_dp, xmax=full_limits(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j))
                call check(roots(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j)%status == r0%status, &
                    'rank 15 full-shape argument status')
                call check(roots(i,1,1,1,1,1,1,1,1,1,1,1,1,1,j)%nfev == r0%nfev, &
                    'rank 15 full-shape argument work count')
            end do
        end do
    end block

    ! Rank-14 refinements have rank-15 leading-endpoint arrays. One of the
    ! four brackets is deliberately invalid; the other three must converge.
    block
        real(dp) :: rb(2,2,1,1,1,1,1,1,1,1,1,1,1,1,2), mb(3,2,1,1,1,1,1,1,1,1,1,1,1,1,2)
        real(dp) :: tol(1,1,1,1,1,1,1,1,1,1,1,1,1,2)
        type(find_root_result) :: roots(2,1,1,1,1,1,1,1,1,1,1,1,1,2), r0
        type(find_minimum_result) :: minima(2,1,1,1,1,1,1,1,1,1,1,1,1,2), m0
        integer :: i, j, batch_calls
        rb = 0.0_dp
        rb(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 3.0_dp
        mb = 0.0_dp
        mb(2,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 1.0_dp
        mb(3,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 3.0_dp
        tol = 1.0e-12_dp
        tol(1,1,1,1,1,1,1,1,1,1,1,1,1,2) = 1.0e-8_dp
        rb(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1) = 3.5_dp
        calls = 0
        roots = find_root(linear_function, rb, xatol=tol)
        batch_calls = calls
        call check(batch_calls == sum(roots%nfev), &
            'rank 14 trailing-axis root refinement accounting')
        calls = 0
        minima = find_minimum(quadratic_function, mb, xatol=tol)
        batch_calls = calls
        call check(batch_calls == sum(minima%nfev), &
            'rank 14 trailing-axis minimum refinement accounting')
        do j = 1, 2
            do i = 1, 2
                r0 = find_root(linear_function, &
                    rb(:,i,1,1,1,1,1,1,1,1,1,1,1,1,j), xatol=tol(1,1,1,1,1,1,1,1,1,1,1,1,1,j))
                m0 = find_minimum(quadratic_function, &
                    mb(:,i,1,1,1,1,1,1,1,1,1,1,1,1,j), xatol=tol(1,1,1,1,1,1,1,1,1,1,1,1,1,j))
                call check(roots(i,1,1,1,1,1,1,1,1,1,1,1,1,j)%status == r0%status, &
                    'rank 14 trailing broadcast root status')
                call check(roots(i,1,1,1,1,1,1,1,1,1,1,1,1,j)%nfev == r0%nfev, &
                    'rank 14 trailing broadcast root evaluations')
                call check(minima(i,1,1,1,1,1,1,1,1,1,1,1,1,j)%status == m0%status, &
                    'rank 14 trailing broadcast minimum status')
                call check(minima(i,1,1,1,1,1,1,1,1,1,1,1,1,j)%nfev == m0%nfev, &
                    'rank 14 trailing broadcast minimum evaluations')
            end do
        end do
        rb(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1) = 0.0_dp
        callback_calls = 0
        roots = find_root(linear_function, rb, callback=stop_root_callback)
        call check(all(roots%status == -4), 'rank 14 root callback termination')
        call check(callback_calls == size(roots), 'rank 14 callback independence')
        callback_calls = 0
        minima = find_minimum(quadratic_function, mb, callback=stop_minimum_callback)
        call check(all(minima%status == -4), 'rank 14 minimum callback termination')
        call check(callback_calls == size(minima), 'rank 14 minimum callback independence')
    end block

    ! Zero extent in a rank-15 batch must not enter index arithmetic or
    ! evaluate the objective.
    block
        real(dp) :: starts(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
        type(bracket_root_result) :: roots(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
        starts = 0.0_dp
        calls = 0
        roots = bracket_root(linear_function, starts, xmax=10.0_dp)
        call check(size(roots) == 0 .and. calls == 0, &
            'empty rank 15 batch is safe')
    end block

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_elementwise_highrank: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_elementwise_highrank: PASS'
contains
    function linear_function(x) result(y)
        real(dp), intent(in) :: x !! scalar position
        real(dp) :: y !! position relative to root
        calls = calls + 1
        y = x - 2.0_dp
    end function linear_function

    function quadratic_function(x) result(y)
        real(dp), intent(in) :: x !! scalar position
        real(dp) :: y !! objective value
        calls = calls + 1
        y = (x - 1.0_dp)**2
    end function quadratic_function

    subroutine stop_root_callback(state, stop)
        use scifort_optimize_root_scalar, only : find_root_result
        type(find_root_result), intent(in) :: state !! one root state
        logical, intent(out) :: stop !! request termination
        callback_calls = callback_calls + 1
        stop = .true.
    end subroutine stop_root_callback

    subroutine stop_minimum_callback(state, stop)
        use scifort_optimize_minimize_scalar, only : find_minimum_result
        type(find_minimum_result), intent(in) :: state !! one minimum state
        logical, intent(out) :: stop !! request termination
        callback_calls = callback_calls + 1
        stop = .true.
    end subroutine stop_minimum_callback

    subroutine check(yes, name)
        logical, intent(in) :: yes !! required condition
        character(len=*), intent(in) :: name !! diagnostic name
        if (yes) return
        failures = failures + 1
        print '(a)', 'FAIL: '//name
    end subroutine check
end program test_optimize_elementwise_highrank
