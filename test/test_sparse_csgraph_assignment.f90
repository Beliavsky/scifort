! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_sparse_csgraph_assignment
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_csgraph, only : csr_weighted_bipartite_graph, &
        weighted_bipartite_from_edges, min_weight_full_bipartite_matching, &
        minimum_weight_full_bipartite_matching, csgraph_success, &
        csgraph_invalid_input, csgraph_no_full_matching
    use scifort_optimize_assignment, only : linear_sum_assignment, assignment_status_success, &
        assignment_status_infeasible
    implicit none
    type(csr_weighted_bipartite_graph) :: graph
    integer, allocatable :: ri(:), ci(:), old_ri(:), old_ci(:), dr(:), dc(:)
    integer :: nr, nc, mask, mode, i, j, k, code, old_code, dense_code, checks, cases
    integer :: rows(9), cols(9)
    real(dp) :: weights(9), costs(3,3), total, old_total, best, infinity
    logical :: present_edge(3,3), used(3), maximum, found

    checks = 0
    cases = 0
    infinity = ieee_value(0.0_dp,ieee_positive_inf)
    ! Exhaust every sparsity pattern through 3x3, including empty partitions,
    ! Hall-deficient graphs, negative/zero costs, transposes, and ties.
    do nr = 0, 3
        do nc = 0, 3
            do mask = 0, 2**(nr*nc)-1
                k = 0
                costs = 0.0_dp
                present_edge = .false.
                do i = 1, nr
                    do j = 1, nc
                        if (.not. btest(mask,(i-1)*nc+j-1)) cycle
                        k = k+1
                        rows(k) = i
                        cols(k) = j
                        weights(k) = real(mod(7*i+11*j+mask,13)-6,dp)
                        costs(i,j) = weights(k)
                        present_edge(i,j) = .true.
                    end do
                end do
                graph = weighted_bipartite_from_edges(nr,nc,rows(:k),cols(:k),weights(:k),code)
                call require(code == csgraph_success,'constructor')
                do mode = 0, 1
                    maximum = mode == 1
                    found = .false.
                    used = .false.
                    best = 0.0_dp
                    call enumerate(1,0.0_dp)
                    call min_weight_full_bipartite_matching(graph,ri,ci,code,maximum,total)
                    call minimum_weight_full_bipartite_matching(graph,old_ri,old_ci, &
                        old_total,old_code,maximum)
                    call require(code == old_code,'original API status')
                    call require(size(ri) == size(old_ri),'original API size')
                    call require(all(ri == old_ri) .and. all(ci == old_ci),'original API pairs')
                    call require(total == old_total,'original API total')
                    if (found) then
                        call require(code == csgraph_success,'feasible status')
                        call require(size(ri) == min(nr,nc),'full cardinality')
                        call require(size(ci) == size(ri),'paired cardinality')
                        call require(total == best,'exhaustive optimum')
                        call require(all(ri(2:) > ri(:size(ri)-1)),'sorted unique rows')
                        do i = 1, size(ri)
                            call require(ri(i) >= 1 .and. ri(i) <= nr,'row range')
                            call require(ci(i) >= 1 .and. ci(i) <= nc,'column range')
                            call require(present_edge(ri(i),ci(i)),'selected edge exists')
                            call require(count(ci == ci(i)) == 1,'unique columns')
                        end do
                        old_total = 0.0_dp
                        do i = 1, size(ri)
                            old_total = old_total+costs(ri(i),ci(i))
                        end do
                        call require(total == old_total,'selected original costs')
                    else
                        call require(code == csgraph_no_full_matching,'infeasible status')
                        call require(size(ri) == 0 .and. size(ci) == 0,'empty failure outputs')
                        call require(total == 0.0_dp,'zero failure total')
                    end if
                    ! Independent dense assignment, with missing edges forbidden.
                    block
                        real(dp) :: dense(nr,nc)
                        dense = merge(costs(:nr,:nc), &
                            merge(-infinity,infinity,maximum),present_edge(:nr,:nc))
                        call linear_sum_assignment(dense,dr,dc,dense_code,maximum)
                        if (found) then
                            call require(dense_code == assignment_status_success,'dense status')
                            old_total = 0.0_dp
                            do i = 1, size(dr)
                                old_total = old_total+dense(dr(i),dc(i))
                            end do
                            call require(old_total == total,'dense assignment optimum')
                        else
                            call require(dense_code == assignment_status_infeasible, &
                                'dense infeasible')
                        end if
                    end block
                    cases = cases+1
                end do
            end do
        end do
    end do

    ! Default options and optional total omitted; parallel edges and stored
    ! zeros are deliberate SciFort extensions to SciPy's sparse semantics.
    graph = weighted_bipartite_from_edges(2,2,[1,1,1,2],[1,1,2,2], &
        [4.0_dp,0.0_dp,-3.0_dp,2.0_dp],code)
    call min_weight_full_bipartite_matching(graph,ri,ci,code)
    call require(code == 0 .and. all(ci == [1,2]),'default options and parallel zero')
    call min_weight_full_bipartite_matching(graph,ri,ci,code,total_weight=total,maximize=.true.)
    call require(code == 0 .and. total == 6.0_dp,'keyword maximum and parallel edges')
    graph%data(1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call min_weight_full_bipartite_matching(graph,ri,ci,code,total_weight=total)
    call require(code == csgraph_invalid_input,'nonfinite cost rejected')
    call require(size(ri) == 0 .and. size(ci) == 0 .and. total == 0.0_dp,'invalid outputs')
    print '(a,i0,a,i0,a)', 'Sparse assignment: ',checks,' checks, ', &
        cases,' exhaustive cases passed.'

contains

    recursive subroutine enumerate(level, value)
        implicit none
        integer, intent(in) :: level !! Next vertex on the smaller side, starting at one.
        real(dp), intent(in) :: value !! Accumulated exact integer-valued cost.
        integer :: partner, r, c
        if (level > min(nr,nc)) then
            if (.not. found) then
                best = value
            else if (maximum) then
                best = max(best,value)
            else
                best = min(best,value)
            end if
            found = .true.
            return
        end if
        do partner = 1, max(nr,nc)
            if (used(partner)) cycle
            r = level
            c = partner
            if (nr > nc) then
                r = partner
                c = level
            end if
            if (.not. present_edge(r,c)) cycle
            used(partner) = .true.
            call enumerate(level+1,value+costs(r,c))
            used(partner) = .false.
        end do
    end subroutine enumerate

    subroutine require(ok, label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
end program test_sparse_csgraph_assignment
