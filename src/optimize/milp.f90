! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dense, bounded-size mixed-integer linear programming by branch-and-bound.
! Each node solves the original SciFort two-phase-simplex LP relaxation.
! This is an independent implementation, not a binding to HiGHS.
module scifort_optimize_milp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_linprog, only : linprog, linprog_result, &
        linprog_optimal, linprog_iteration_limit, linprog_infeasible, &
        linprog_unbounded
    implicit none
    private

    integer, parameter, public :: milp_optimal = 0
    integer, parameter, public :: milp_node_limit = 1
    integer, parameter, public :: milp_infeasible = 2
    integer, parameter, public :: milp_unbounded = 3
    integer, parameter, public :: milp_numerical = 4
    integer, parameter, public :: milp_invalid = -1

    type, public :: milp_result
        real(dp), allocatable :: x(:) !! incumbent solution, if available
        real(dp) :: fun = huge(1.0_dp) !! incumbent c^T*x or HUGE if absent
        real(dp) :: mip_dual_bound = -huge(1.0_dp) !! lower bound on global optimum
        real(dp) :: mip_gap = huge(1.0_dp) !! relative primal-dual gap, if known
        integer :: mip_node_count = 0 !! number of solved LP relaxations
        integer :: nit = 0 !! total simplex pivot count
        integer :: status = milp_invalid !! 0 optimal/gap, 1 limit, 2 infeasible, 3 unbounded, 4 numerical, -1 invalid
        logical :: success = .false. !! true iff optimal or requested gap achieved
        character(len=180) :: message = 'invalid mixed-integer programming input' !! termination explanation
    end type milp_result

    type :: bb_node
        real(dp), allocatable :: lower(:) !! effective node lower bounds
        real(dp), allocatable :: upper(:) !! effective node upper bounds
        real(dp) :: bound = -huge(1.0_dp) !! valid inherited LP lower bound
    end type bb_node

    public :: milp

contains

    function milp(c, integrality, a, b_l, b_u, lower, upper, node_limit, &
            mip_rel_gap, tol, lp_maxiter, time_limit) result(res)
        real(dp), intent(in) :: c(:) !! objective coefficients, length n
        integer, intent(in), optional :: integrality(:) !! n entries: 0 continuous, 1 integer, 2 semi-continuous, 3 semi-integer
        real(dp), intent(in), optional :: a(:, :) !! m-by-n constraint matrix, dense
        real(dp), intent(in), optional :: b_l(:) !! m lower row bounds, -infinity permitted
        real(dp), intent(in), optional :: b_u(:) !! m upper row bounds, +infinity permitted
        real(dp), intent(in), optional :: lower(:) !! n variable lower bounds, default zero
        real(dp), intent(in), optional :: upper(:) !! n variable upper bounds, default +infinity
        integer, intent(in), optional :: node_limit !! maximum LP relaxations, default 10000
        real(dp), intent(in), optional :: mip_rel_gap !! stopping gap, default 1.0e-4
        real(dp), intent(in), optional :: tol !! feasibility/LP tolerance, default 1.0e-9
        integer, intent(in), optional :: lp_maxiter !! maximum simplex pivots per node, default solver-dependent
        real(dp), intent(in), optional :: time_limit !! approximate CPU-seconds limit, checked between nodes
        type(milp_result) :: res !! incumbent and branch-and-bound diagnostics

        type(bb_node), allocatable :: pending(:), expanded(:)
        type(bb_node) :: current
        type(linprog_result) :: lp
        real(dp), allocatable :: original_lo(:), original_hi(:), root_lo(:), root_hi(:)
        real(dp), allocatable :: aub(:, :), bub(:), aeq(:, :), beq(:), candidate(:)
        real(dp), allocatable :: child_lo(:), child_hi(:)
        integer, allocatable :: itype(:)
        real(dp) :: infinity, tolerance, integral_tol, gap_target, bound, &
            elapsed, start_cpu, split, value, score, best_score
        integer :: n, m, nu, ne, j, i, chosen, idx, active, limit, lp_limit, &
            selected_kind
        logical :: found, valid, has_integer

        n = size(c)
        if (.not. all(ieee_is_finite(c))) return
        if (present(a) .neqv. present(b_l)) return
        if (present(a) .neqv. present(b_u)) return
        m = 0
        if (present(a)) then
            m = size(a, 1)
            if (size(a, 2) /= n .or. size(b_l) /= m .or. size(b_u) /= m) return
            if (.not. all(ieee_is_finite(a))) return
            if (any(ieee_is_nan(b_l)) .or. any(ieee_is_nan(b_u))) return
        end if
        if (present(lower)) then
            if (size(lower) /= n) return
            if (any(ieee_is_nan(lower))) return
        end if
        if (present(upper)) then
            if (size(upper) /= n) return
            if (any(ieee_is_nan(upper))) return
        end if
        allocate(itype(n))
        itype = 0
        if (present(integrality)) then
            if (size(integrality) /= n) return
            itype = integrality
            if (any(itype < 0) .or. any(itype > 3)) return
        end if
        has_integer = any(itype /= 0)
        tolerance = 1.0e-9_dp
        if (present(tol)) tolerance = tol
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0e-2_dp) return
        integral_tol = max(1.0e-7_dp, 10.0_dp*tolerance)
        gap_target = 1.0e-4_dp
        if (present(mip_rel_gap)) gap_target = mip_rel_gap
        if (.not. ieee_is_finite(gap_target)) return
        if (gap_target < 0.0_dp) return
        limit = 10000
        if (present(node_limit)) limit = node_limit
        if (limit < 0) return
        lp_limit = max(500, 40*(n+m+1))
        if (present(lp_maxiter)) lp_limit = lp_maxiter
        if (lp_limit < 0) return
        if (present(time_limit)) then
            if (.not. ieee_is_finite(time_limit)) return
            if (time_limit < 0.0_dp) return
        end if

        infinity = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(original_lo(n), original_hi(n), root_lo(n), root_hi(n))
        original_lo = 0.0_dp
        original_hi = infinity
        if (present(lower)) original_lo = lower
        if (present(upper)) original_hi = upper
        if (any(original_lo > original_hi)) return
        if (any(original_lo >= infinity) .or. any(original_hi <= -infinity)) return
        if (present(a)) then
            if (any(b_l > b_u)) return
            if (any(b_l >= infinity) .or. any(b_u <= -infinity)) return
        end if

        root_lo = original_lo
        root_hi = original_hi
        do j = 1, n
            if (itype(j) == 2 .or. itype(j) == 3) then
                ! For semi-variables the disjunction is {0} union [lower,upper].
                root_lo(j) = min(0.0_dp, original_lo(j))
                root_hi(j) = max(0.0_dp, original_hi(j))
            end if
        end do

        ! Convert finite two-sided row bounds to equality and <= rows.
        nu = 0
        ne = 0
        if (present(a)) then
            do i = 1, m
                if (b_l(i) >= b_u(i) .and. ieee_is_finite(b_l(i))) then
                    ne = ne + 1
                else
                    if (ieee_is_finite(b_l(i))) nu = nu + 1
                    if (ieee_is_finite(b_u(i))) nu = nu + 1
                end if
            end do
        end if
        allocate(aub(nu,n), bub(nu), aeq(ne,n), beq(ne))
        nu = 0
        ne = 0
        if (present(a)) then
            do i = 1, m
                if (b_l(i) >= b_u(i) .and. ieee_is_finite(b_l(i))) then
                    ne = ne + 1
                    aeq(ne,:) = a(i,:)
                    beq(ne) = b_l(i)
                else
                    if (ieee_is_finite(b_u(i))) then
                        nu = nu + 1
                        aub(nu,:) = a(i,:)
                        bub(nu) = b_u(i)
                    end if
                    if (ieee_is_finite(b_l(i))) then
                        nu = nu + 1
                        aub(nu,:) = -a(i,:)
                        bub(nu) = -b_l(i)
                    end if
                end if
            end do
        end if

        ! A branch-and-bound queue ordered by inherited lower bounds.
        ! Every unexplored child has a valid bound from its parent relaxation.
        allocate(pending(32), candidate(n), child_lo(n), child_hi(n))
        active = 1
        pending(1)%lower = root_lo
        pending(1)%upper = root_hi
        pending(1)%bound = -huge(1.0_dp)
        found = .false.
        res%status = milp_node_limit
        res%message = 'mixed-integer search stopped by node limit'
        call cpu_time(start_cpu)
        do while (active > 0)
            call update_bound_gap()
            if (found .and. res%mip_gap <= gap_target) then
                res%status = milp_optimal
                res%success = .true.
                res%message = 'requested relative mixed-integer gap reached'
                return
            end if
            if (res%mip_node_count >= limit) exit
            if (present(time_limit)) then
                call cpu_time(elapsed)
                if (elapsed - start_cpu >= time_limit) then
                    res%message = 'mixed-integer search stopped by CPU time limit'
                    exit
                end if
            end if

            idx = 1
            do i = 2, active
                if (pending(i)%bound < pending(idx)%bound) idx = i
            end do
            current = pending(idx)
            pending(idx) = pending(active)
            active = active - 1
            ! Discard the popped node's superseded allocatable buffers.
            if (allocated(pending(active+1)%lower)) deallocate(pending(active+1)%lower)
            if (allocated(pending(active+1)%upper)) deallocate(pending(active+1)%upper)

            if (found .and. current%bound >= res%fun - &
                    8.0_dp*tolerance*max(1.0_dp,abs(res%fun))) cycle
            lp = linprog(c, a_ub=aub, b_ub=bub, a_eq=aeq, b_eq=beq, &
                lower=current%lower, upper=current%upper, tol=tolerance, maxiter=lp_limit)
            res%mip_node_count = res%mip_node_count + 1
            res%nit = res%nit + lp%nit
            select case (lp%status)
            case (linprog_infeasible)
                cycle
            case (linprog_iteration_limit)
                res%status = milp_node_limit
                res%message = 'LP relaxation reached the simplex iteration limit'
                ! Still include the popped node's inherited lower bound.
                call push_node(current%lower, current%upper, current%bound)
                exit
            case (linprog_unbounded)
                if (.not. has_integer) then
                    res%status = milp_unbounded
                    res%message = 'continuous LP objective is unbounded below'
                else
                    res%status = milp_numerical
                    res%message = 'unbounded LP relaxation; mixed-integer status unresolved'
                end if
                return
            case (linprog_optimal)
                continue
            case default
                res%status = milp_numerical
                res%message = 'LP relaxation failed numerically'
                return
            end select
            if (.not. allocated(lp%x)) then
                res%status = milp_numerical
                res%message = 'LP relaxation returned no solution'
                return
            end if
            if (size(lp%x) /= n .or. .not. all(ieee_is_finite(lp%x)) .or. &
                    .not. ieee_is_finite(lp%fun)) then
                res%status = milp_numerical
                res%message = 'LP relaxation returned a nonfinite solution'
                return
            end if
            if (.not. relaxation_feasible(lp%x, current%lower, current%upper)) then
                res%status = milp_numerical
                res%message = 'LP relaxation failed independent feasibility check'
                return
            end if
            bound = lp%fun
            if (found .and. bound >= res%fun - &
                    8.0_dp*tolerance*max(1.0_dp,abs(res%fun))) cycle

            ! First branch on forbidden semi-variable gaps, then integrality.
            chosen = 0
            selected_kind = 0
            do j = 1, n
                if (itype(j) /= 2 .and. itype(j) /= 3) cycle
                value = lp%x(j)
                if (abs(value) <= integral_tol) cycle
                if (value >= original_lo(j) - integral_tol .and. &
                    value <= original_hi(j) + integral_tol) cycle
                chosen = j
                selected_kind = 1
                exit
            end do
            if (chosen == 0) then
                best_score = integral_tol
                do j = 1, n
                    if (itype(j) /= 1 .and. itype(j) /= 3) cycle
                    score = abs(lp%x(j) - anint(lp%x(j)))
                    if (score > best_score) then
                        best_score = score
                        chosen = j
                        selected_kind = 2
                    end if
                end do
            end if

            if (chosen == 0) then
                ! Correct near-integral values and check original constraints.
                candidate = lp%x
                do j = 1, n
                    if (itype(j) == 1 .or. itype(j) == 3) candidate(j) = anint(candidate(j))
                    if ((itype(j) == 2 .or. itype(j) == 3) .and. &
                        abs(candidate(j)) <= integral_tol) candidate(j) = 0.0_dp
                end do
                valid = original_feasible(candidate)
                if (.not. valid) then
                    ! Rounding at tolerance may upset a tight equality: split
                    ! the largest residual instead of accepting a false incumbent.
                    best_score = -1.0_dp
                    do j = 1, n
                        if (itype(j) /= 1 .and. itype(j) /= 3) cycle
                        score = abs(lp%x(j) - anint(lp%x(j)))
                        if (score > best_score) then
                            best_score = score
                            chosen = j
                            selected_kind = 2
                        end if
                    end do
                    if (chosen == 0 .or. best_score <= epsilon(1.0_dp)) then
                        res%status = milp_numerical
                        res%message = 'integer rounding failed independent feasibility check'
                        return
                    end if
                else
                    value = dot_product(c,candidate)
                    if (.not. ieee_is_finite(value)) then
                        res%status = milp_numerical
                        res%message = 'nonfinite integer objective'
                        return
                    end if
                    if (.not. found .or. value < res%fun) then
                        res%x = candidate
                        res%fun = value
                        found = .true.
                    end if
                    cycle
                end if
            end if

            j = chosen
            if (selected_kind == 1) then
                ! Disjoint semi-variable alternatives: x=0 or x in its
                ! original [lower,upper] interval, even if lower > 0.
                child_lo = current%lower
                child_hi = current%upper
                child_lo(j) = max(child_lo(j), 0.0_dp)
                child_hi(j) = min(child_hi(j), 0.0_dp)
                call push_node(child_lo, child_hi, bound)
                child_lo = current%lower
                child_hi = current%upper
                child_lo(j) = max(child_lo(j), original_lo(j))
                child_hi(j) = min(child_hi(j), original_hi(j))
                call push_node(child_lo, child_hi, bound)
            else
                ! REAL-valued floor avoids integer overflow for large bounds.
                split = aint(lp%x(j))
                if (split > lp%x(j)) split = split - 1.0_dp
                child_lo = current%lower
                child_hi = current%upper
                child_hi(j) = min(child_hi(j), split)
                call push_node(child_lo, child_hi, bound)
                child_lo = current%lower
                child_hi = current%upper
                child_lo(j) = max(child_lo(j), split + 1.0_dp)
                call push_node(child_lo, child_hi, bound)
            end if
        end do

        if (active == 0) then
            if (found) then
                res%status = milp_optimal
                res%success = .true.
                res%message = 'mixed-integer optimum certified by branch-and-bound'
                res%mip_dual_bound = res%fun
                res%mip_gap = 0.0_dp
            else
                res%status = milp_infeasible
                res%message = 'mixed-integer constraints are infeasible'
            end if
        else
            call update_bound_gap()
        end if

    contains

        subroutine push_node(l, u, parent_bound)
            real(dp), intent(in) :: l(:) !! new node lower variable bounds
            real(dp), intent(in) :: u(:) !! new node upper variable bounds
            real(dp), intent(in) :: parent_bound !! valid LP lower bound inherited from parent
            if (any(l > u)) return
            if (found) then
                if (parent_bound >= res%fun - 8.0_dp*tolerance*max(1.0_dp,abs(res%fun))) return
            end if
            if (active == size(pending)) then
                allocate(expanded(2*size(pending)))
                expanded(1:active) = pending(1:active)
                call move_alloc(expanded, pending)
            end if
            active = active + 1
            pending(active)%lower = l
            pending(active)%upper = u
            pending(active)%bound = parent_bound
        end subroutine push_node

        subroutine update_bound_gap()
            real(dp) :: dual
            integer :: t
            if (active == 0) then
                if (found) then
                    res%mip_dual_bound = res%fun
                    res%mip_gap = 0.0_dp
                end if
                return
            end if
            dual = pending(1)%bound
            do t = 2, active
                dual = min(dual, pending(t)%bound)
            end do
            if (found) dual = min(dual, res%fun)
            res%mip_dual_bound = dual
            res%mip_gap = huge(1.0_dp)
            if (found .and. dual > -0.5_dp*huge(1.0_dp)) then
                res%mip_gap = max(0.0_dp, (res%fun-dual)/max(1.0_dp,abs(res%fun)))
            end if
        end subroutine update_bound_gap

        logical function relaxation_feasible(x, l, u) result(ok)
            real(dp), intent(in) :: x(:) !! candidate relaxation solution
            real(dp), intent(in) :: l(:) !! current node lower bounds
            real(dp), intent(in) :: u(:) !! current node upper bounds
            real(dp), allocatable :: v(:)
            real(dp) :: permitted, row_magnitude
            integer :: t
            ok = .false.
            do t = 1, n
                if (ieee_is_finite(l(t))) then
                    permitted = feasibility_error(l(t), abs(x(t)))
                    if (x(t) < l(t)-permitted) return
                end if
                if (ieee_is_finite(u(t))) then
                    permitted = feasibility_error(u(t), abs(x(t)))
                    if (x(t) > u(t)+permitted) return
                end if
            end do
            if (present(a)) then
                v = matmul(a,x)
                if (.not. all(ieee_is_finite(v))) return
                do t = 1, m
                    row_magnitude = sum(abs(a(t,:)*x))
                    if (.not. ieee_is_finite(row_magnitude)) return
                    if (ieee_is_finite(b_l(t))) then
                        permitted = feasibility_error(b_l(t), row_magnitude)
                        if (v(t) < b_l(t)-permitted) return
                    end if
                    if (ieee_is_finite(b_u(t))) then
                        permitted = feasibility_error(b_u(t), row_magnitude)
                        if (v(t) > b_u(t)+permitted) return
                    end if
                end do
            end if
            ok = .true.
        end function relaxation_feasible

        real(dp) function feasibility_error(reference, magnitude) result(permitted)
            real(dp), intent(in) :: reference !! finite constraint/bound right-hand side
            real(dp), intent(in) :: magnitude !! scale of the computed left-hand side
            ! An absolute user tolerance plus a rounding allowance prevents
            ! accepting order-one violations of constraints around 1.0e12.
            permitted = max(32.0_dp*tolerance, &
                128.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(reference), magnitude))
        end function feasibility_error

        logical function original_feasible(x) result(ok)
            real(dp), intent(in) :: x(:) !! candidate solution after integer rounding
            real(dp), allocatable :: v(:)
            real(dp) :: permitted, row_magnitude
            integer :: t
            ok = .false.
            if (.not. all(ieee_is_finite(x))) return
            do t = 1, n
                if (itype(t) == 2 .or. itype(t) == 3) then
                    if (abs(x(t)) <= 0.0_dp) cycle
                end if
                if (ieee_is_finite(original_lo(t))) then
                    permitted = feasibility_error(original_lo(t), abs(x(t)))
                    if (x(t) < original_lo(t)-permitted) return
                end if
                if (ieee_is_finite(original_hi(t))) then
                    permitted = feasibility_error(original_hi(t), abs(x(t)))
                    if (x(t) > original_hi(t)+permitted) return
                end if
                if (itype(t) == 1 .or. itype(t) == 3) then
                    if (abs(x(t) - anint(x(t))) > integral_tol) return
                end if
            end do
            if (present(a)) then
                v = matmul(a,x)
                if (.not. all(ieee_is_finite(v))) return
                do t = 1, m
                    row_magnitude = sum(abs(a(t,:)*x))
                    if (.not. ieee_is_finite(row_magnitude)) return
                    if (ieee_is_finite(b_l(t))) then
                        permitted = feasibility_error(b_l(t), row_magnitude)
                        if (v(t) < b_l(t)-permitted) return
                    end if
                    if (ieee_is_finite(b_u(t))) then
                        permitted = feasibility_error(b_u(t), row_magnitude)
                        if (v(t) > b_u(t)+permitted) return
                    end if
                end do
            end if
            ok = .true.
        end function original_feasible
    end function milp
end module scifort_optimize_milp
