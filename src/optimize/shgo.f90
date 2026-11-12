! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent, deterministic bounded simplicial global search.
! The vertices of successive dyadic Kuhn/Freudenthal triangulations form the
! zero-dimensional sublevel-set complex. A local minimum of the vertex graph
! is a birth of a connected component (H_0), and supplies a polishing start.
! This is a deliberately restricted SHGO variant, not a port of SciPy's
! algorithm or its convergence guarantees. Optional linear/nonlinear
! constraint bounds use feasibility-first graph ranking and a constrained
! local-polishing callback (SLSQP through the public adapter).
! Optional adaptive refinement subdivides selected dyadic cells, using
! cached vertices and the edges of each cell's Kuhn simplices. With
! conforming=.true. in two dimensions, every leaf is triangulated as a fan
! about a sampled interior point. Boundary segments include all hanging
! vertices. In three dimensions each cube is tetrahedralized by coning
! conforming shared face triangles to its cached interior center. Face-center
! vertices and dyadic subface splits provide matching coarse-fine interfaces.
! The default adaptive graph and dimensions greater than three are nonconforming.
! Optional homology-growth stopping uses the number of sampled graph minima
! (H_0 births) as a heuristic rank estimate. It is not a certificate of
! convergence or a complete persistent-homology calculation.
!
! Geometric basis: Freudenthal (1942), triangulation of the n-cube.
! Optimization basis: Endres, Sandrock and Focke (2018), A simplicial homology
! algorithm for Lipschitz optimisation, J. Global Optim. 72, 181-217.
!
module scifort_optimize_shgo
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: shgo_completed = 0
    integer, parameter, public :: shgo_maxfun = 1
    integer, parameter, public :: shgo_target = 2
    integer, parameter, public :: shgo_callback_stop = 3
    integer, parameter, public :: shgo_mesh_limit = 4
    integer, parameter, public :: shgo_hgr_stable = 5
    integer, parameter, public :: shgo_infeasible = 6
    integer, parameter, public :: shgo_invalid = -1
    integer, parameter, public :: shgo_numerical = -2
    integer, parameter :: maximum_dimension = 8
    integer, parameter :: maximum_level = 20
    ! One more dyadic bit supports cell centers at the deepest permitted level.
    integer, parameter :: lattice_level = maximum_level+1

    type, public :: shgo_result
        real(dp), allocatable :: x(:) !! best point, length n
        real(dp) :: fun = huge(1.0_dp) !! best finite objective found
        real(dp), allocatable :: xl(:,:) !! distinct local minima, shape (n,nminima)
        real(dp), allocatable :: funl(:) !! associated minima, ascending order
        integer :: nfev = 0 !! actual calls to the user's objective, including local search
        integer :: nlfev = 0 !! calls made during local searches
        integer :: ncev = 0 !! constraint-vector evaluations, including local endpoint checks
        integer :: ncjev = 0 !! user nonlinear Jacobian evaluations
        integer :: nrestore = 0 !! attempted feasibility-restoration starts
        integer :: nrestored = 0 !! feasible points accepted through restoration
        real(dp) :: maxcv = 0.0_dp !! maximum component constraint violation at returned x
        logical :: feasible = .false. !! returned x satisfies all constraint bounds
        integer :: nit = 0 !! completed mesh refinement levels
        integer :: nlocal = 0 !! local minimizer invocations
        integer :: nminima = 0 !! stored distinct local minima
        integer :: nmesh = 0 !! sampled mesh vertices (last complete uniform mesh)
        integer :: nleaf = 0 !! active leaf cells, adaptive mode only
        integer :: nrefined = 0 !! number of adaptively subdivided cells
        integer :: ntriangles = 0 !! triangles in last validated conforming 2D mesh
        integer :: ntetrahedra = 0 !! tetrahedra in last validated conforming 3D mesh
        logical :: conforming_valid = .false. !! last built conforming mesh passed topology checks
        integer :: hgr_estimate = -1 !! graph minima on most recent completed mesh, -1 unknown
        integer :: hgr_streak = 0 !! consecutive rounds with unchanged rank and no material gain
        integer :: status = shgo_invalid !! termination status constant
        logical :: success = .false. !! at least one finite feasible minimum was found
        character(len=160) :: message = 'invalid input' !! termination explanation
    end type shgo_result

    ! Dyadic lattice keys permit exact reuse of vertices on adjacent and
    ! nested cells, independent of their floating-point coordinate values.
    type :: adaptive_mesh
        integer :: dim = 0
        integer :: cap = 0
        integer :: npoint = 0
        integer :: ncell = 0
        integer(int64), allocatable :: keys(:,:) !! integer lattice coordinates
        integer(int64), allocatable :: origins(:,:) !! low corner of each leaf
        integer(int64), allocatable :: spans(:) !! cell side length in lattice units
        integer, allocatable :: slots(:) !! open-addressed vertex lookup
        real(dp), allocatable :: points(:,:) !! mesh vertex coordinates
        real(dp), allocatable :: values(:) !! cached graph ranking values
        real(dp), allocatable :: raw(:) !! raw feasible objective, huge for infeasible vertices
    end type adaptive_mesh

    abstract interface
        function shgo_objective(x) result(f)
            import dp
            real(dp), intent(in) :: x(:) !! vector of n bounded coordinates
            real(dp) :: f !! scalar objective, nonfinite values treated as infeasible
        end function shgo_objective

        subroutine shgo_local_minimizer(fun, lower, upper, start, maxfev, x, value, nfev, success)
            import dp, shgo_objective
            procedure(shgo_objective) :: fun !! budget-guarded objective
            real(dp), intent(in) :: lower(:) !! inclusive lower bounds, length n
            real(dp), intent(in) :: upper(:) !! inclusive upper bounds, length n
            real(dp), intent(in) :: start(:) !! feasible starting point, length n
            integer, intent(in) :: maxfev !! remaining objective-call budget, >= 1
            real(dp), intent(out) :: x(:) !! minimizer, length n
            real(dp), intent(out) :: value !! objective at the minimizer
            integer, intent(out) :: nfev !! local solver's reported evaluation count
            logical, intent(out) :: success !! whether the returned solution is usable
        end subroutine shgo_local_minimizer

        subroutine shgo_constraint_function(x, c)
            import dp
            real(dp), intent(in) :: x(:) !! candidate vector, length nvar
            real(dp), intent(out) :: c(:) !! vector-valued nonlinear constraints, length m
        end subroutine shgo_constraint_function

        subroutine shgo_constraint_jacobian(x, a)
            import dp
            real(dp), intent(in) :: x(:) !! bounded optimization variables (nvar)
            real(dp), intent(out) :: a(:,:) !! nonlinear Jacobian (m,nvar), one row per component
        end subroutine shgo_constraint_jacobian

        subroutine shgo_callback(x, value, stop)
            import dp
            real(dp), intent(in) :: x(:) !! best point after a refinement iteration
            real(dp), intent(in) :: value !! best finite objective so far
            logical, intent(out) :: stop !! request clean early termination
        end subroutine shgo_callback
    end interface

    ! Independent nonlinear constraint blocks, analogous to SciPy's sequence
    ! of NonlinearConstraint objects. Callbacks are borrowed: their procedure
    ! targets must remain valid until shgo_core/shgo returns.
    type, public :: shgo_nonlinear_constraint
        procedure(shgo_constraint_function), pointer, nopass :: fun => null() !! vector-valued constraints
        procedure(shgo_constraint_jacobian), pointer, nopass :: jac => null() !! optional analytic rows
        real(dp), allocatable :: lb(:) !! lower limits, length equal to the callback output
        real(dp), allocatable :: ub(:) !! upper limits, matching lb
    end type shgo_nonlinear_constraint

    public :: shgo_core, shgo_objective, shgo_local_minimizer, shgo_callback
    public :: shgo_constraint_function, shgo_constraint_jacobian
    public :: shgo_constraint_blocks_valid

contains

    function shgo_constraint_blocks_valid(constraints) result(valid)
        type(shgo_nonlinear_constraint), intent(in) :: constraints(:) !! independent constraint records
        logical :: valid !! all procedures, shapes and scalar bounds are valid
        integer :: k
        valid = .false.
        do k = 1, size(constraints)
            if (.not. associated(constraints(k)%fun)) return
            if (.not. allocated(constraints(k)%lb)) return
            if (.not. allocated(constraints(k)%ub)) return
            if (size(constraints(k)%lb) /= size(constraints(k)%ub)) return
            if (any(ieee_is_nan(constraints(k)%lb)) .or. &
                    any(ieee_is_nan(constraints(k)%ub))) return
            if (any(constraints(k)%lb > constraints(k)%ub)) return
            if (any((.not. ieee_is_finite(constraints(k)%lb)) .and. &
                    constraints(k)%lb > 0.0_dp)) return
            if (any((.not. ieee_is_finite(constraints(k)%ub)) .and. &
                    constraints(k)%ub < 0.0_dp)) return
        end do
        valid = .true.
    end function shgo_constraint_blocks_valid


    function shgo_core(fun, bounds, n, iters, maxfun, max_mesh_points, local_iter, &
            minimize_every_iter, f_min, f_tol, merge_tol, minimizer, callback, &
            adaptive, adaptive_fraction, hgr_stability_iters, hgr_min_iters, &
            hgr_f_tol, conforming, linear_a, linear_lb, linear_ub, nonlinear, &
            nonlinear_lb, nonlinear_ub, constraint_tol, nonlinear_jac, &
            restore_feasibility, restoration_starts, restoration_maxiter, &
            constraints) result(res)
        procedure(shgo_objective) :: fun !! scalar bounded objective function
        real(dp), intent(in) :: bounds(:,:) !! finite bounds, shape (nvar,2), lower < upper
        integer, intent(in), optional :: n !! minimum vertices in first full grid, default 3**nvar
        integer, intent(in), optional :: iters !! number of nested dyadic grid levels, default 1
        integer, intent(in), optional :: maxfun !! strict call cap, default 20000
        integer, intent(in), optional :: max_mesh_points !! per-grid vertex cap, default 100000
        integer, intent(in), optional :: local_iter !! maximum local starts per level, default all
        logical, intent(in), optional :: minimize_every_iter !! polish on every level, default true
        real(dp), intent(in), optional :: f_min !! optional known global minimum
        real(dp), intent(in), optional :: f_tol !! relative target accuracy, default 1e-8
        real(dp), intent(in), optional :: merge_tol !! normalized dedup threshold, default 1e-4
        procedure(shgo_local_minimizer), optional :: minimizer !! optional bounded local search
        procedure(shgo_callback), optional :: callback !! end-of-level callback
        logical, intent(in), optional :: adaptive !! use selective dyadic cell refinement
        real(dp), intent(in), optional :: adaptive_fraction !! fraction of leaves refined per round, default 0.25
        integer, intent(in), optional :: hgr_stability_iters !! stable rank comparisons needed, 0 disables
        integer, intent(in), optional :: hgr_min_iters !! earliest completed round for HGR stop, default 3
        real(dp), intent(in), optional :: hgr_f_tol !! max relative best-value improvement for stability, 1e-8
        logical, intent(in), optional :: conforming !! conforming adaptive triangulation in 2D/3D, default false
        real(dp), intent(in), optional :: linear_a(:,:) !! linear constraint matrix (m,nvar)
        real(dp), intent(in), optional :: linear_lb(:) !! lower limits on A*x (m); -infinity allowed
        real(dp), intent(in), optional :: linear_ub(:) !! upper limits on A*x (m); +infinity allowed
        procedure(shgo_constraint_function), optional :: nonlinear !! vector nonlinear callback
        real(dp), intent(in), optional :: nonlinear_lb(:) !! lower nonlinear limits (m)
        real(dp), intent(in), optional :: nonlinear_ub(:) !! upper nonlinear limits (m)
        real(dp), intent(in), optional :: constraint_tol !! absolute component tolerance, default 1e-8
        procedure(shgo_constraint_jacobian), optional :: nonlinear_jac !! derivatives of nonlinear(x,c)
        logical, intent(in), optional :: restore_feasibility !! opt-in projected feasibility recovery
        integer, intent(in), optional :: restoration_starts !! starts per mesh level, default 3
        integer, intent(in), optional :: restoration_maxiter !! maximum projection iterations, default 25
        type(shgo_nonlinear_constraint), intent(in), optional :: constraints(:) !! independent nonlinear blocks
        type(shgo_result) :: res !! point, local minima and diagnostics

        integer :: dim, requested, iterations, eval_limit, mesh_limit, local_limit
        integer :: level0, level, iteration, q, prev_q, grid_size, j, k, idx
        integer :: digit, digit_index, previous_index, old_stride, ncand, cand_limit
        integer :: locally_reported, nlocal_stored, local_capacity
        integer :: old_mesh_size
        integer :: required_streak, minimum_rounds, previous_rank
        integer :: nlin, nnonlin
        type(adaptive_mesh) :: mesh
        integer, allocatable :: triangles(:,:), tetrahedra(:,:)
        logical :: too_many_tetrahedra
        logical :: triangulation_ok
        integer, allocatable :: candidate(:)
        real(dp), allocatable :: previous(:), previous_raw(:), values(:), raw_values(:)
        real(dp), allocatable :: point(:), local_point(:), linear_values(:), nonlinear_values(:)
        real(dp), allocatable :: all_nonlinear_lb(:), all_nonlinear_ub(:)
        real(dp), allocatable :: local_points(:,:), local_values(:), target(:), widths(:)
        real(dp) :: fx, fcenter, fcenter_raw, sampled_raw, f_target, target_tol
        real(dp) :: duplicate_tol, trial_value, cv_tol, point_cv, best_violation
        logical :: do_every, has_target, done, stop_request, local_ok, accepted, local_budget_hit
        logical :: adaptive_mode, conforming_mode, has_constraints, local_feasible
        logical :: saw_feasible_constraint_point, restoration_enabled
        integer :: restore_starts, restore_limit, restoration_index, restored_seed_count
        real(dp) :: adapt_fraction
        real(dp) :: hgr_tol, best_before_round
        logical :: hgr_enabled, tentative_hgr_stop, have_candidates
        integer(int64) :: total, needed

        dim = size(bounds, 1)
        if (dim < 1 .or. dim > maximum_dimension .or. size(bounds,2) /= 2) then
            res%message = 'bounds must be (n,2), with 1 <= n <= 8'
            return
        end if
        if (.not. all(ieee_is_finite(bounds))) then
            res%message = 'bounds must be finite'
            return
        end if
        if (any(bounds(:,2) <= bounds(:,1))) then
            res%message = 'each upper bound must exceed its lower bound'
            return
        end if
        do j = 1, dim
            if (bounds(j,1) < 0.0_dp .and. bounds(j,2) > 0.0_dp) then
                if (bounds(j,2) > huge(1.0_dp) + bounds(j,1)) then
                    res%message = 'bound widths must not overflow'
                    return
                end if
            end if
        end do
        allocate(widths(dim))
        widths = bounds(:,2) - bounds(:,1)
        if (.not. all(ieee_is_finite(widths))) then
            res%message = 'bound widths must not overflow'
            return
        end if

        iterations = 1
        if (present(iters)) iterations = iters
        eval_limit = 20000
        if (present(maxfun)) eval_limit = maxfun
        mesh_limit = 100000
        if (present(max_mesh_points)) mesh_limit = max_mesh_points
        local_limit = huge(1)
        if (present(local_iter)) local_limit = local_iter
        requested = 0
        if (present(n)) requested = n
        duplicate_tol = 1.0e-4_dp
        if (present(merge_tol)) duplicate_tol = merge_tol
        target_tol = 1.0e-8_dp
        if (present(f_tol)) target_tol = f_tol
        has_target = present(f_min)
        f_target = 0.0_dp
        if (has_target) f_target = f_min
        do_every = .true.
        if (present(minimize_every_iter)) do_every = minimize_every_iter
        adaptive_mode = .false.
        if (present(adaptive)) adaptive_mode = adaptive
        conforming_mode = .false.
        if (present(conforming)) conforming_mode = conforming
        adapt_fraction = 0.25_dp
        if (present(adaptive_fraction)) adapt_fraction = adaptive_fraction
        required_streak = 0
        if (present(hgr_stability_iters)) required_streak = hgr_stability_iters
        minimum_rounds = 3
        if (present(hgr_min_iters)) minimum_rounds = hgr_min_iters
        hgr_tol = 1.0e-8_dp
        if (present(hgr_f_tol)) hgr_tol = hgr_f_tol
        hgr_enabled = required_streak > 0
        restoration_enabled = .false.
        if (present(restore_feasibility)) restoration_enabled = restore_feasibility
        restore_starts = 3
        if (present(restoration_starts)) restore_starts = restoration_starts
        restore_limit = 25
        if (present(restoration_maxiter)) restore_limit = restoration_maxiter
        previous_rank = -1
        cv_tol = 1.0e-8_dp
        if (present(constraint_tol)) cv_tol = constraint_tol
        if (.not. ieee_is_finite(cv_tol)) then
            res%message = 'constraint_tol must be finite'
            return
        end if
        if (cv_tol < 0.0_dp) then
            res%message = 'constraint_tol must be nonnegative'
            return
        end if
        if ((present(linear_a) .neqv. present(linear_lb)) .or. &
                (present(linear_a) .neqv. present(linear_ub))) then
            res%message = 'linear_a, linear_lb and linear_ub must be supplied together'
            return
        end if
        if (present(constraints)) then
            if (present(nonlinear) .or. present(nonlinear_jac) .or. &
                    present(nonlinear_lb) .or. present(nonlinear_ub)) then
                res%message = 'constraints cannot be combined with legacy nonlinear arguments'
                return
            end if
            if (.not. shgo_constraint_blocks_valid(constraints)) then
                res%message = 'invalid nonlinear constraint block or bounds'
                return
            end if
        end if
        if (present(nonlinear_jac) .and. .not. present(nonlinear)) then
            res%message = 'nonlinear_jac requires nonlinear constraints'
            return
        end if
        if (restore_starts < 0 .or. restore_limit < 1) then
            res%message = 'invalid feasibility-restoration limits'
            return
        end if
        if ((present(nonlinear) .neqv. present(nonlinear_lb)) .or. &
                (present(nonlinear) .neqv. present(nonlinear_ub))) then
            res%message = 'nonlinear and its bound arrays must be supplied together'
            return
        end if
        nlin = 0
        nnonlin = 0
        if (present(linear_a)) then
            nlin = size(linear_a,1)
            if (size(linear_a,2) /= dim .or. size(linear_lb) /= nlin .or. &
                    size(linear_ub) /= nlin) then
                res%message = 'linear constraint array shape mismatch'
                return
            end if
            if (.not. all(ieee_is_finite(linear_a)) .or. &
                    any(ieee_is_nan(linear_lb)) .or. any(ieee_is_nan(linear_ub))) then
                res%message = 'invalid linear constraint coefficients or bounds'
                return
            end if
            if (any(linear_lb > linear_ub) .or. &
                    any((.not. ieee_is_finite(linear_lb)) .and. linear_lb > 0.0_dp) .or. &
                    any((.not. ieee_is_finite(linear_ub)) .and. linear_ub < 0.0_dp)) then
                res%message = 'invalid linear constraint coefficients or bounds'
                return
            end if
            allocate(linear_values(nlin))
        end if
        if (present(constraints)) then
            do j = 1, size(constraints)
                if (size(constraints(j)%lb) > huge(1)-nnonlin) then
                    res%message = 'too many nonlinear constraint components'
                    return
                end if
                nnonlin = nnonlin + size(constraints(j)%lb)
            end do
            allocate(all_nonlinear_lb(nnonlin), all_nonlinear_ub(nnonlin))
            k = 0
            do j = 1, size(constraints)
                idx = size(constraints(j)%lb)
                all_nonlinear_lb(k+1:k+idx) = constraints(j)%lb
                all_nonlinear_ub(k+1:k+idx) = constraints(j)%ub
                k = k+idx
            end do
            allocate(nonlinear_values(nnonlin))
        else if (present(nonlinear)) then
            nnonlin = size(nonlinear_lb)
            if (nnonlin /= size(nonlinear_ub)) then
                res%message = 'nonlinear constraint bound lengths differ'
                return
            end if
            if (any(ieee_is_nan(nonlinear_lb)) .or. any(ieee_is_nan(nonlinear_ub))) then
                res%message = 'invalid nonlinear constraint bounds'
                return
            end if
            if (any(nonlinear_lb > nonlinear_ub) .or. &
                    any((.not. ieee_is_finite(nonlinear_lb)) .and. nonlinear_lb > 0.0_dp) .or. &
                    any((.not. ieee_is_finite(nonlinear_ub)) .and. nonlinear_ub < 0.0_dp)) then
                res%message = 'invalid nonlinear constraint bounds'
                return
            end if
            allocate(nonlinear_values(nnonlin))
            all_nonlinear_lb = nonlinear_lb
            all_nonlinear_ub = nonlinear_ub
        end if
        has_constraints = nlin > 0 .or. nnonlin > 0
        best_violation = huge(1.0_dp)
        saw_feasible_constraint_point = .false.

        if (conforming_mode .and. (.not. adaptive_mode .or. (dim /= 2 .and. dim /= 3))) then
            res%message = 'conforming mode requires adaptive=.true. and 2 or 3 variables'
            return
        end if
        if (iterations < 1 .or. eval_limit < 1 .or. mesh_limit < 1 .or. &
                local_limit < 0 .or. requested < 0 .or. required_streak < 0 .or. &
                minimum_rounds < 2) then
            res%message = 'invalid iteration, sampling, local or evaluation limit'
            return
        end if
        if (.not. ieee_is_finite(duplicate_tol) .or. .not. ieee_is_finite(target_tol) .or. &
                .not. ieee_is_finite(f_target) .or. .not. ieee_is_finite(adapt_fraction) .or. &
                .not. ieee_is_finite(hgr_tol)) then
            res%message = 'target and tolerance parameters must be finite'
            return
        end if
        if (adapt_fraction <= 0.0_dp .or. adapt_fraction > 1.0_dp) then
            res%message = 'adaptive_fraction must be in (0,1]'
            return
        end if
        if (duplicate_tol < 0.0_dp .or. target_tol < 0.0_dp .or. hgr_tol < 0.0_dp) then
            res%message = 'tolerances must be nonnegative'
            return
        end if

        ! For this structured simplicial variant n selects the first dyadic
        ! grid level; unlike SciPy it is not an exact sampling count.
        level0 = 1
        do while (level0 < maximum_level)
            q = 2**level0 + 1
            total = capped_grid_size(q,dim,int(requested,int64))
            if (total >= int(requested, int64)) exit
            level0 = level0 + 1
        end do
        if (level0 + int(iterations, int64) - 1_int64 > maximum_level) then
            res%message = 'requested dyadic refinement level is too large'
            return
        end if

        allocate(res%x(dim), point(dim), local_point(dim), target(dim))
        local_capacity = min(32, eval_limit)
        allocate(local_points(dim,local_capacity), local_values(local_capacity))
        nlocal_stored = 0
        target = bounds(:,1) + 0.5_dp * widths
        res%x = target
        fcenter = counted_objective(target)
        fcenter_raw = sampled_raw
        res%status = shgo_completed
        res%message = 'requested refinement levels completed'
        if (target_reached()) then
            res%status = shgo_target
            res%message = 'known target objective reached'
            goto 900
        end if

        prev_q = 0
        old_mesh_size = 0
        if (adaptive_mode) allocate(candidate(min(mesh_limit,eval_limit)))
        do iteration = 1, iterations
            best_before_round = res%fun
            tentative_hgr_stop = .false.
            have_candidates = .false.
            level = level0 + iteration - 1
            q = 2**level + 1
            if (adaptive_mode) then
                if (iteration == 1) then
                    call initialize_adaptive_mesh(q)
                else
                    call refine_adaptive_mesh()
                end if
                if (res%status /= shgo_completed) exit
                if (conforming_mode) then
                    if (dim == 2) then
                        call build_conforming_triangles(mesh, triangles, triangulation_ok)
                        if (triangulation_ok) then
                            res%ntriangles = size(triangles,2)
                            call find_triangle_candidates(mesh, triangles, candidate, ncand)
                        end if
                    else
                        call build_conforming_tetrahedra(mesh, mesh_limit, tetrahedra, &
                            triangulation_ok, too_many_tetrahedra)
                        if (too_many_tetrahedra) then
                            res%status = shgo_mesh_limit
                            res%message = 'conforming tetrahedral mesh exceeds simplex capacity'
                            exit
                        end if
                        if (triangulation_ok) then
                            res%ntetrahedra = size(tetrahedra,2)
                            call find_tetrahedral_candidates(mesh, tetrahedra, candidate, ncand)
                        end if
                    end if
                    if (.not. triangulation_ok) then
                        res%status = shgo_numerical
                        res%message = 'conforming mesh failed face, edge or volume validation'
                        exit
                    end if
                    res%conforming_valid = .true.
                else
                    call find_adaptive_candidates(mesh, candidate, ncand)
                end if
                call sort_candidates(candidate, ncand, mesh%values)
                have_candidates = .true.
            else
                total = capped_grid_size(q,dim,int(mesh_limit,int64))
                if (total > int(mesh_limit, int64)) then
                    res%status = shgo_mesh_limit
                    res%message = 'next simplicial grid exceeds max_mesh_points'
                    exit
                end if
                grid_size = int(total)
                ! All old vertices (and the original center) are reused exactly.
                needed = total - int(old_mesh_size, int64)
                if (old_mesh_size == 0) needed = total - 1_int64
                if (needed > int(eval_limit - res%nfev, int64)) then
                    res%status = shgo_maxfun
                    res%message = 'insufficient objective budget for next complete mesh'
                    exit
                end if
                allocate(values(grid_size), raw_values(grid_size), candidate(grid_size))
                do idx = 1, grid_size
                    digit_index = idx - 1
                    old_stride = 1
                    previous_index = 1
                    done = old_mesh_size > 0
                    accepted = .true.
                    do j = 1, dim
                        digit = mod(digit_index,q)
                        digit_index = digit_index / q
                        point(j) = bounds(j,1) + (real(digit, dp) / real(q-1, dp))*widths(j)
                        if (mod(digit,2) /= 0) done = .false.
                        if (digit /= (q-1)/2) accepted = .false.
                        if (old_mesh_size > 0) then
                            previous_index = previous_index + (digit/2)*old_stride
                            old_stride = old_stride * prev_q
                        end if
                    end do
                    if (done) then
                        fx = previous(previous_index)
                        sampled_raw = previous_raw(previous_index)
                    else if (old_mesh_size == 0 .and. accepted) then
                        fx = fcenter
                        sampled_raw = fcenter_raw
                    else
                        fx = counted_objective(point)
                    end if
                    values(idx) = fx
                    raw_values(idx) = sampled_raw
                    if (target_reached()) then
                        res%status = shgo_target
                        res%message = 'known target objective reached'
                        goto 900
                    end if
                end do
                res%nit = iteration
                res%nmesh = grid_size

            end if
            if (hgr_enabled .or. do_every .or. iteration == iterations) then
                if (.not. adaptive_mode) then
                    ncand = 0
                    do idx = 1, grid_size
                        if (.not. valid_value(values(idx))) cycle
                        if (.not. is_graph_minimum(idx, q, dim, values)) cycle
                        ncand = ncand + 1
                        candidate(ncand) = idx
                    end do
                    call sort_candidates(candidate, ncand, values)
                end if
                have_candidates = .true.
                ! A potential stop must still perform any deferred local
                ! polishing before the stopping criterion is accepted.
                if (hgr_enabled .and. ncand > 0 .and. ncand == previous_rank) then
                    tentative_hgr_stop = res%hgr_streak+1 >= required_streak .and. &
                        iteration >= minimum_rounds .and. &
                        no_material_improvement(best_before_round,res%fun)
                end if
            end if
            if (restoration_enabled .and. has_constraints .and. have_candidates) then
                ! Start with graph minima, then draw additional independent
                ! vertices from the cache. On a flat violation plateau all
                ! boundary vertices may be one connected graph minimum; using
                ! only that representative can miss nonzero derivatives.
                restored_seed_count = min(ncand,restore_starts)
                do restoration_index = 1, restored_seed_count
                    idx = candidate(restoration_index)
                    if (adaptive_mode) then
                        point = mesh%points(:,idx)
                    else
                        call point_from_index(idx,q,bounds,widths,point)
                    end if
                    call restore_candidate(point)
                    if (target_reached()) then
                        res%status = shgo_target
                        res%message = 'known target objective reached'
                        goto 900
                    end if
                end do
                if (adaptive_mode) then
                    grid_size = mesh%npoint
                end if
                do idx = 1, grid_size
                    if (restored_seed_count >= restore_starts) exit
                    if (any(candidate(1:min(ncand,restore_starts)) == idx)) cycle
                    if (adaptive_mode) then
                        point = mesh%points(:,idx)
                    else
                        call point_from_index(idx,q,bounds,widths,point)
                    end if
                    call restore_candidate(point)
                    restored_seed_count = restored_seed_count+1
                    if (target_reached()) then
                        res%status = shgo_target
                        res%message = 'known target objective reached'
                        goto 900
                    end if
                end do
            end if
            if (do_every .or. iteration == iterations .or. tentative_hgr_stop) then
                ! Unpolished coarse-grid minima are not necessarily stationary
                ! on a refined complex. Rebuild their pool at each level.
                if (.not. present(minimizer)) nlocal_stored = 0
                cand_limit = min(ncand, local_limit)
                do k = 1, cand_limit
                    idx = candidate(k)
                    if (adaptive_mode) then
                        point = mesh%points(:,idx)
                        fx = mesh%raw(idx)
                    else
                        call point_from_index(idx, q, bounds, widths, point)
                        fx = raw_values(idx)
                    end if
                    if (present(minimizer)) then
                        if (res%nfev >= eval_limit) then
                            res%status = shgo_maxfun
                            res%message = 'objective budget exhausted during local refinement'
                            exit
                        end if
                        local_point = point
                        trial_value = fx
                        locally_reported = 0
                        local_ok = .false.
                        res%nlocal = res%nlocal + 1
                        local_budget_hit = .false.
                        call minimizer(local_objective, bounds(:,1), bounds(:,2), &
                            point, eval_limit-res%nfev, local_point, trial_value, &
                            locally_reported, local_ok)
                        local_feasible = .true.
                        point_cv = 0.0_dp
                        if (has_constraints .and. local_ok .and. valid_value(trial_value)) then
                            ! Fortran logical expressions do not short circuit.
                            ! Do not compare a NaN return against the bounds.
                            if (all(ieee_is_finite(local_point))) then
                                if (all(local_point >= bounds(:,1)) .and. &
                                        all(local_point <= bounds(:,2))) then
                                    call evaluate_violation(local_point, point_cv)
                                    local_feasible = point_cv <= cv_tol
                                    call record_violation(local_point,point_cv)
                                    if (local_feasible) then
                                        ! Never trust a custom minimizer's reported
                                        ! objective without an actual evaluation.
                                        ! The extra call is budget guarded.
                                        if (res%nfev < eval_limit) then
                                            trial_value = local_objective(local_point)
                                        else
                                            trial_value = huge(1.0_dp)
                                        end if
                                    end if
                                end if
                            end if
                        end if
                        if (.not. local_feasible .or. .not. local_ok .or. &
                                .not. valid_value(trial_value) .or. &
                                .not. all(ieee_is_finite(local_point))) then
                            local_point = point
                            trial_value = fx
                        else if (any(local_point < bounds(:,1)) .or. &
                                any(local_point > bounds(:,2))) then
                            local_point = point
                            trial_value = fx
                        else if (trial_value > fx) then
                            local_point = point
                            trial_value = fx
                        end if
                        ! Count actual objective calls, never a solver-provided estimate.
                        if (local_budget_hit) then
                            res%status = shgo_maxfun
                            res%message = 'local solver attempted calls beyond maxfun'
                        end if
                    else
                        local_point = point
                        trial_value = fx
                    end if
                    call add_minimum(local_point, trial_value)
                    if (trial_value < res%fun .or. .not. valid_value(res%fun)) then
                        res%fun = trial_value
                        res%x = local_point
                        res%maxcv = point_cv
                        res%feasible = .true.
                    end if
                    if (target_reached()) then
                        res%status = shgo_target
                        res%message = 'known target objective reached'
                        exit
                    end if
                end do
                if (res%status /= shgo_completed) exit
            end if
            if (have_candidates) then
                res%hgr_estimate = ncand
                if (hgr_enabled) then
                    if (ncand > 0 .and. ncand == previous_rank .and. &
                            no_material_improvement(best_before_round, res%fun)) then
                        res%hgr_streak = res%hgr_streak + 1
                    else
                        res%hgr_streak = 0
                    end if
                    previous_rank = ncand
                end if
            end if
            if (present(callback) .and. valid_value(res%fun)) then
                stop_request = .false.
                call callback(res%x, res%fun, stop_request)
                if (stop_request) then
                    res%status = shgo_callback_stop
                    res%message = 'callback requested termination'
                    exit
                end if
            end if
            if (hgr_enabled) then
                if (res%hgr_streak >= required_streak .and. iteration >= minimum_rounds) then
                    res%status = shgo_hgr_stable
                    res%message = 'sampled graph-minimum count stabilized'
                    exit
                end if
            end if
            if (iteration == iterations) exit
            if (.not. adaptive_mode) then
                call move_alloc(values, previous)
                call move_alloc(raw_values, previous_raw)
                old_mesh_size = grid_size
                prev_q = q
                deallocate(candidate)
            end if
        end do
900     continue
        ! Keep the best actual evaluation in the sorted result pool, even
        ! when a custom local solver aborts after visiting an improved point.
        if (valid_value(res%fun)) then
            call add_minimum(res%x, res%fun)
        end if
        if (valid_value(res%fun)) then
            res%success = .true.
            res%feasible = .true.
        else if (has_constraints) then
            res%maxcv = best_violation
            if (res%status == shgo_completed .or. res%status == shgo_hgr_stable) then
                if (saw_feasible_constraint_point) then
                    res%status = shgo_numerical
                    res%message = 'all feasible objective evaluations were nonfinite'
                else
                    res%status = shgo_infeasible
                    res%message = 'no feasible point found'
                end if
            end if
        else
            res%status = shgo_numerical
            res%message = 'objective returned no finite feasible values'
        end if
        call sort_minima(local_points, local_values, nlocal_stored)
        res%nminima = nlocal_stored
        allocate(res%xl(dim,nlocal_stored), res%funl(nlocal_stored))
        if (nlocal_stored > 0) then
            res%xl = local_points(:,1:nlocal_stored)
            res%funl = local_values(1:nlocal_stored)
        end if

    contains

        subroutine restore_candidate(seed)
            ! Damped projected Gauss-Newton for the squared active constraint
            ! residual. This is a feasibility heuristic, not a global solver.
            real(dp), intent(in) :: seed(:) !! mesh vertex supplying initial coordinates
            real(dp) :: x(dim), trial(dim), residual(nlin+nnonlin), proposed(nlin+nnonlin)
            real(dp) :: jacobian(nlin+nnonlin,dim), scaled(nlin+nnonlin,dim)
            real(dp) :: normal(dim,dim), rhs(dim), delta(dim)
            real(dp) :: cv, new_cv, penalty, new_penalty, alpha, damping, scale
            real(dp) :: value
            integer :: t, j, l
            logical :: ok, accepted

            if (res%nfev >= eval_limit) return
            res%nrestore = res%nrestore + 1
            x = seed
            call constraint_residual(x,cv,residual,ok)
            if (.not. ok) return
            if (cv <= cv_tol) return ! a feasible mesh vertex is already sampled
            do t = 1, restore_limit
                call restoration_jacobian(x,jacobian,ok,nonlinear_jac)
                if (.not. ok) return
                do j = 1, dim
                    scaled(:,j) = jacobian(:,j)*widths(j)
                end do
                if (.not. all(ieee_is_finite(scaled))) return
                normal = matmul(transpose(scaled),scaled)
                rhs = matmul(transpose(scaled),residual)
                if (.not. all(ieee_is_finite(normal))) return
                if (.not. all(ieee_is_finite(rhs))) return
                damping = 1.0e-8_dp*max(1.0_dp,maxval(abs(normal)))
                do j = 1, dim
                    normal(j,j) = normal(j,j) + damping
                end do
                call shgo_solve_dense(normal,rhs,delta,ok)
                if (.not. ok) return
                scale = maxval(abs(delta))
                if (.not. ieee_is_finite(scale) .or. scale == 0.0_dp) return
                if (scale > 0.5_dp) delta = delta*(0.5_dp/scale)
                penalty = sum(residual*residual)
                if (.not. ieee_is_finite(penalty)) return
                accepted = .false.
                alpha = 1.0_dp
                do l = 1, 16
                    trial = max(bounds(:,1),min(bounds(:,2),x + alpha*widths*delta))
                    if (all(trial == x)) exit
                    call constraint_residual(trial,new_cv,proposed,ok)
                    if (ok) then
                        new_penalty = sum(proposed*proposed)
                        if (ieee_is_finite(new_penalty)) then
                            if (new_penalty < penalty*(1.0_dp-1.0e-4_dp*alpha) .or. &
                                    new_cv <= cv_tol) then
                                accepted = .true.
                                exit
                            end if
                        end if
                    end if
                    alpha = 0.5_dp*alpha
                end do
                if (.not. accepted) return
                x = trial
                cv = new_cv
                residual = proposed
                if (cv <= cv_tol) then
                    if (res%nfev >= eval_limit) return
                    value = counted_objective(x)
                    if (valid_value(sampled_raw)) then
                        res%nrestored = res%nrestored + 1
                        call add_minimum(x,sampled_raw)
                    end if
                    return
                end if
            end do
        end subroutine restore_candidate

        subroutine constraint_residual(x,cv,residual,ok)
            real(dp), intent(in) :: x(:) !! bounded point for feasibility evaluation
            real(dp), intent(out) :: cv !! max violation across components
            real(dp), intent(out) :: residual(:) !! signed corrections for active bounds
            logical, intent(out) :: ok !! all component evaluations are finite
            real(dp) :: v
            integer :: k
            call evaluate_violation(x,cv)
            call record_violation(x,cv)
            ok = ieee_is_finite(cv) .and. cv < huge(1.0_dp)/4.0_dp
            residual = 0.0_dp
            if (.not. ok) return
            do k = 1, nlin+nnonlin
                if (k <= nlin) then
                    v = linear_values(k)
                    if (ieee_is_finite(linear_lb(k))) then
                        if (linear_lb(k)-v > cv_tol) residual(k) = linear_lb(k)-v
                    end if
                    if (ieee_is_finite(linear_ub(k))) then
                        if (v-linear_ub(k) > cv_tol) residual(k) = linear_ub(k)-v
                    end if
                else
                    v = nonlinear_values(k-nlin)
                    if (ieee_is_finite(all_nonlinear_lb(k-nlin))) then
                        if (all_nonlinear_lb(k-nlin)-v > cv_tol) residual(k) = all_nonlinear_lb(k-nlin)-v
                    end if
                    if (ieee_is_finite(all_nonlinear_ub(k-nlin))) then
                        if (v-all_nonlinear_ub(k-nlin) > cv_tol) residual(k) = all_nonlinear_ub(k-nlin)-v
                    end if
                end if
            end do
            ok = all(ieee_is_finite(residual))
        end subroutine constraint_residual

        subroutine restoration_jacobian(x,jacobian,ok,nonlinear_jac)
            procedure(shgo_constraint_jacobian), optional :: nonlinear_jac !! legacy Jacobian callback
            real(dp), intent(in) :: x(:) !! bounded point for Jacobian evaluation
            real(dp), intent(out) :: jacobian(:,:) !! original constraint derivatives (m,nvar)
            logical, intent(out) :: ok !! finite, valid computed derivatives
            integer :: block, offset, m
            ok = .true.
            jacobian = 0.0_dp
            if (nlin > 0) jacobian(1:nlin,:) = linear_a
            if (nnonlin == 0) return
            if (present(constraints)) then
                offset = nlin
                do block = 1, size(constraints)
                    m = size(constraints(block)%lb)
                    if (m == 0) cycle
                    if (associated(constraints(block)%jac)) then
                        call constraints(block)%jac(x,jacobian(offset+1:offset+m,:))
                        res%ncjev = res%ncjev + 1
                    else
                        call finite_difference_block(x,jacobian(offset+1:offset+m,:), &
                            constraints(block)%fun,ok)
                        if (.not. ok) return
                    end if
                    offset = offset + m
                end do
            else if (present(nonlinear_jac)) then
                call nonlinear_jac(x,jacobian(nlin+1:,:))
                res%ncjev = res%ncjev + 1
            else
                call finite_difference_block(x,jacobian(nlin+1:,:),nonlinear,ok)
                if (.not. ok) return
            end if
            ok = all(ieee_is_finite(jacobian))
        end subroutine restoration_jacobian

        subroutine finite_difference_block(x, jacobian, evaluate, ok)
            real(dp), intent(in) :: x(:) !! bound-constrained point
            real(dp), intent(out) :: jacobian(:,:) !! derivatives for one nonlinear block
            procedure(shgo_constraint_function) :: evaluate !! nonlinear constraint callback
            logical, intent(out) :: ok !! all computed values and derivatives finite
            real(dp) :: original(size(jacobian,1)), shifted(size(jacobian,1))
            real(dp) :: perturbed(dim), h
            integer :: j
            ok = .false.
            call evaluate(x,original)
            res%ncev = res%ncev + 1
            if (.not. all(ieee_is_finite(original))) return
            do j = 1, dim
                perturbed = x
                h = max(1.0e-6_dp*widths(j), &
                    sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(j))))
                perturbed(j) = min(bounds(j,2),x(j)+h)
                if (perturbed(j) == x(j)) perturbed(j) = max(bounds(j,1),x(j)-h)
                if (perturbed(j) == x(j)) return
                call evaluate(perturbed,shifted)
                res%ncev = res%ncev + 1
                if (.not. all(ieee_is_finite(shifted))) return
                jacobian(:,j) = (shifted-original)/(perturbed(j)-x(j))
            end do
            ok = all(ieee_is_finite(jacobian))
        end subroutine finite_difference_block

        subroutine evaluate_nonlinear(x, values, nonlinear)
            procedure(shgo_constraint_function), optional :: nonlinear !! legacy nonlinear callback
            real(dp), intent(in) :: x(:) !! optimization variables
            real(dp), intent(out) :: values(:) !! flattened nonlinear constraint values
            integer :: block, offset, m
            if (present(constraints)) then
                offset = 0
                do block = 1, size(constraints)
                    m = size(constraints(block)%lb)
                    if (m > 0) call constraints(block)%fun(x,values(offset+1:offset+m))
                    offset = offset+m
                end do
            else
                call nonlinear(x,values)
            end if
        end subroutine evaluate_nonlinear

        subroutine evaluate_violation(x, cv)
            real(dp), intent(in) :: x(:) !! point to check against linear/nonlinear bounds
            real(dp), intent(out) :: cv !! largest component violation (0 when feasible)
            integer :: l
            cv = 0.0_dp
            if (.not. has_constraints) return
            res%ncev = res%ncev + 1
            if (nlin > 0) then
                linear_values = matmul(linear_a,x)
                do l = 1, nlin
                    call accumulate_cv(linear_values(l),linear_lb(l),linear_ub(l),cv)
                end do
            end if
            if (nnonlin > 0) then
                call evaluate_nonlinear(x,nonlinear_values,nonlinear)
                do l = 1, nnonlin
                    call accumulate_cv(nonlinear_values(l),all_nonlinear_lb(l),all_nonlinear_ub(l),cv)
                end do
            end if
        end subroutine evaluate_violation

        subroutine accumulate_cv(value, lower, upper, cv)
            real(dp), intent(in) :: value !! constraint component value
            real(dp), intent(in) :: lower !! allowed lower limit
            real(dp), intent(in) :: upper !! allowed upper limit
            real(dp), intent(inout) :: cv !! accumulated largest violation
            if (.not. ieee_is_finite(value)) then
                cv = huge(1.0_dp)
                return
            end if
            if (ieee_is_finite(lower)) cv = max(cv, max(0.0_dp,lower-value))
            if (ieee_is_finite(upper)) cv = max(cv, max(0.0_dp,value-upper))
        end subroutine accumulate_cv

        subroutine record_violation(x,cv)
            real(dp), intent(in) :: x(:) !! evaluated candidate
            real(dp), intent(in) :: cv !! maximum constraint violation
            if (.not. has_constraints) return
            if (size(x) /= dim) return
            if (.not. all(ieee_is_finite(x))) return
            if (any(x < bounds(:,1)) .or. any(x > bounds(:,2))) return
            if (cv < best_violation) then
                best_violation = cv
                if (.not. valid_value(res%fun)) then
                    res%x = x
                    res%maxcv = cv
                end if
            end if
        end subroutine record_violation

        subroutine initialize_adaptive_mesh(q_start)
            integer, intent(in) :: q_start !! points on each axis of initial mesh
            integer :: index, flat, j, digit, cell_count, side, axis
            integer(int64) :: nvertices, nleaves, spacing, nfaces, conforming_count
            integer(int64) :: key(dim)
            real(dp) :: sample(dim), sample_value

            nvertices = capped_grid_size(q_start, dim, int(mesh_limit,int64))
            if (nvertices > int(mesh_limit,int64)) then
                res%status = shgo_mesh_limit
                res%message = 'initial adaptive grid exceeds max_mesh_points'
                return
            end if
            if (nvertices > int(eval_limit,int64)) then
                res%status = shgo_maxfun
                res%message = 'insufficient budget for initial adaptive grid'
                return
            end if
            nleaves = capped_grid_size(q_start-1, dim, nvertices)
            if (nleaves > nvertices) then
                res%status = shgo_mesh_limit
                res%message = 'initial adaptive leaf count exceeds mesh capacity'
                return
            end if
            if (conforming_mode) then
                ! Initial cell centers and (in 3D) unique shared face centers.
                ! Preflight the entire first mesh before sampling any vertex.
                nfaces = 0_int64
                if (dim == 3) nfaces = 3_int64*int(q_start,int64)* &
                    int(q_start-1,int64)**2
                conforming_count = nvertices + nleaves + nfaces
                if (conforming_count > int(mesh_limit,int64)) then
                    res%status = shgo_mesh_limit
                    res%message = 'initial conforming mesh exceeds max_mesh_points'
                    return
                end if
                if (conforming_count > int(eval_limit,int64)) then
                    res%status = shgo_maxfun
                    res%message = 'insufficient budget for initial conforming mesh'
                    return
                end if
            end if
            call setup_adaptive_mesh(mesh, dim, min(mesh_limit,eval_limit))
            spacing = 2_int64**(lattice_level-level0)
            key = 2_int64**(lattice_level-1)
            call insert_adaptive_vertex(mesh, key, target, fcenter, fcenter_raw)
            grid_size = int(nvertices)
            do index = 1, grid_size
                flat = index-1
                do j = 1, dim
                    digit = mod(flat,q_start)
                    flat = flat/q_start
                    key(j) = int(digit,int64)*spacing
                end do
                if (lookup_adaptive_vertex(mesh,key) /= 0) cycle
                sample = bounds(:,1) + real(key,dp) / real(2_int64**lattice_level,dp)*widths
                sample_value = counted_objective(sample)
                call insert_adaptive_vertex(mesh,key,sample,sample_value,sampled_raw)
            end do
            cell_count = int(nleaves)
            do index = 1, cell_count
                flat = index-1
                do j = 1, dim
                    digit = mod(flat,q_start-1)
                    flat = flat/(q_start-1)
                    mesh%origins(j,index) = int(digit,int64)*spacing
                end do
                mesh%spans(index) = spacing
            end do
            mesh%ncell = cell_count
            if (conforming_mode) then
                do index = 1, mesh%ncell
                    key = mesh%origins(:,index) + mesh%spans(index)/2_int64
                    sample = bounds(:,1) + real(key,dp)/real(2_int64**lattice_level,dp)*widths
                    sample_value = counted_objective(sample)
                    call insert_adaptive_vertex(mesh,key,sample,sample_value,sampled_raw)
                end do
                if (dim == 3) then
                    do index = 1, mesh%ncell
                        do axis = 1, dim
                            do side = 0, 1
                                key = mesh%origins(:,index) + mesh%spans(index)/2_int64
                                key(axis) = mesh%origins(axis,index) + &
                                    int(side,int64)*mesh%spans(index)
                                if (lookup_adaptive_vertex(mesh,key) /= 0) cycle
                                sample = bounds(:,1) + real(key,dp)/real(2_int64**lattice_level,dp)*widths
                                sample_value = counted_objective(sample)
                                call insert_adaptive_vertex(mesh,key,sample,sample_value,sampled_raw)
                            end do
                        end do
                    end do
                end if
            end if
            res%nmesh = mesh%npoint
            res%nleaf = mesh%ncell
            res%nit = 1
            if (target_reached()) then
                res%status = shgo_target
                res%message = 'known target objective reached'
            end if
        end subroutine initialize_adaptive_mesh

        subroutine refine_adaptive_mesh()
            ! Split a deterministic fraction of the most promising active
            ! leaves. Each complete split adds a 3**dim dyadic subgrid;
            ! vertices shared with neighbors are fetched from the hash cache.
            integer, allocatable :: order(:)
            real(dp), allocatable :: scores(:)
            integer :: original_cells, selected, rank, leaf, vertex, ix
            integer :: digit, j, child, mask, index, nfinite, axis, side, nstaged, staged_index
            real(dp) :: avg, observed(2**dim+1)
            integer(int64) :: parent(dim), key(dim), child_step
            integer(int64) :: missing
            integer(int64) :: staged(dim,3**dim + 2**dim*(1+2*dim))
            real(dp) :: sample(dim), sample_value

            original_cells = mesh%ncell
            allocate(order(original_cells),scores(original_cells))
            do leaf = 1, original_cells
                order(leaf) = leaf
                scores(leaf) = huge(1.0_dp)
                nfinite = 0
                do mask = 0, 2**dim-1
                    key = mesh%origins(:,leaf)
                    do j = 1, dim
                        if (btest(mask,j-1)) key(j) = key(j)+mesh%spans(leaf)
                    end do
                    vertex = lookup_adaptive_vertex(mesh,key)
                    if (vertex <= 0) error stop 'SHGO: missing adaptive vertex'
                    if (.not. valid_value(mesh%values(vertex))) cycle
                    nfinite = nfinite+1
                    observed(nfinite) = mesh%values(vertex)
                    scores(leaf) = min(scores(leaf),mesh%values(vertex))
                end do
                if (conforming_mode) then
                    key = mesh%origins(:,leaf) + mesh%spans(leaf)/2_int64
                    vertex = lookup_adaptive_vertex(mesh,key)
                    if (vertex <= 0) error stop 'SHGO: missing conforming cell center'
                    if (valid_value(mesh%values(vertex))) then
                        nfinite = nfinite+1
                        observed(nfinite) = mesh%values(vertex)
                        scores(leaf) = min(scores(leaf),mesh%values(vertex))
                    end if
                end if
                if (nfinite > 0) then
                    ! A small average-corner penalty resolves tied minima
                    ! between cells sharing a best vertex. Summing divided
                    ! terms avoids overflow when objective values are large.
                    avg = sum(observed(1:nfinite)/real(nfinite,dp))
                    scores(leaf) = 0.875_dp*scores(leaf) + 0.125_dp*avg
                end if
            end do
            call sort_candidates(order,original_cells,scores)
            selected = max(1, ceiling(adapt_fraction*real(original_cells,dp)))
            do rank = 1, selected
                leaf = order(rank)
                child_step = mesh%spans(leaf)/2_int64
                if (child_step < 1_int64 .or. (conforming_mode .and. child_step < 2_int64)) then
                    res%status = shgo_mesh_limit
                    res%message = 'adaptive dyadic refinement depth exceeded'
                    return
                end if
                if (mesh%ncell > mesh%cap - (2**dim-1)) then
                    res%status = shgo_mesh_limit
                    res%message = 'adaptive leaf-cell capacity exceeded'
                    return
                end if
                parent = mesh%origins(:,leaf)
                nstaged = 0
                do index = 0, 3**dim-1
                    ix = index
                    key = parent
                    do j = 1, dim
                        digit = mod(ix,3)
                        ix = ix/3
                        key(j) = key(j) + int(digit,int64)*child_step
                    end do
                    call stage_adaptive_vertex(mesh,key,staged,nstaged)
                end do
                if (conforming_mode) then
                    do child = 0, 2**dim-1
                        key = parent + child_step/2_int64
                        do j = 1, dim
                            if (btest(child,j-1)) key(j) = key(j)+child_step
                        end do
                        call stage_adaptive_vertex(mesh,key,staged,nstaged)
                        if (dim == 3) then
                            do axis = 1, 3
                                do side = 0, 1
                                    key = parent + child_step/2_int64
                                    do j = 1, dim
                                        if (btest(child,j-1)) key(j) = key(j)+child_step
                                    end do
                                    key(axis) = parent(axis)+int(side,int64)*child_step
                                    if (btest(child,axis-1)) key(axis) = key(axis)+child_step
                                    call stage_adaptive_vertex(mesh,key,staged,nstaged)
                                end do
                            end do
                        end if
                    end do
                end if
                missing = int(nstaged,int64)
                if (missing > int(eval_limit-res%nfev,int64)) then
                    res%status = shgo_maxfun
                    res%message = 'insufficient budget for next adaptive cell split'
                    return
                end if
                if (missing > int(mesh_limit-mesh%npoint,int64)) then
                    res%status = shgo_mesh_limit
                    res%message = 'adaptive vertices exceed max_mesh_points'
                    return
                end if
                do staged_index = 1, nstaged
                    key = staged(:,staged_index)
                    sample = bounds(:,1) + real(key,dp)/real(2_int64**lattice_level,dp)*widths
                    sample_value = counted_objective(sample)
                    call insert_adaptive_vertex(mesh,key,sample,sample_value,sampled_raw)
                end do
                mesh%spans(leaf) = child_step
                do child = 1, 2**dim-1
                    mesh%ncell = mesh%ncell + 1
                    mesh%origins(:,mesh%ncell) = parent
                    do j = 1, dim
                        if (btest(child,j-1)) then
                            mesh%origins(j,mesh%ncell) = parent(j)+child_step
                        end if
                    end do
                    mesh%spans(mesh%ncell) = child_step
                end do
                res%nrefined = res%nrefined+1
                res%nleaf = mesh%ncell
                res%nmesh = mesh%npoint
                if (conforming_mode) then
                    ! A partially committed round is not yet triangulated.
                    res%conforming_valid = .false.
                    res%ntriangles = 0
                    res%ntetrahedra = 0
                end if
                if (target_reached()) then
                    res%status = shgo_target
                    res%message = 'known target objective reached'
                    return
                end if
            end do
            res%nit = iteration
        end subroutine refine_adaptive_mesh

        function counted_objective(x) result(value)
            real(dp), intent(in) :: x(:) !! point within rectangular bounds
            real(dp) :: value !! graph rank; equals raw f for unconstrained runs
            real(dp) :: cv
            value = huge(1.0_dp)
            sampled_raw = huge(1.0_dp)
            cv = 0.0_dp
            if (has_constraints) then
                call evaluate_violation(x,cv)
                call record_violation(x,cv)
                if (cv > cv_tol) then
                    ! Feasibility-first mesh ranking: all feasible ranks are
                    ! below -1; infeasible ranks are in [0,1]. This admits
                    ! local restoration starts even for equality surfaces
                    ! containing no sampled mesh vertex.
                    if (ieee_is_finite(cv)) then
                        if (cv > huge(1.0_dp)/2.0_dp) then
                            value = 1.0_dp
                        else
                            value = cv/(1.0_dp+cv)
                        end if
                    end if
                    return
                end if
                saw_feasible_constraint_point = .true.
            end if
            if (res%nfev >= eval_limit) return
            res%nfev = res%nfev + 1
            sampled_raw = fun(x)
            if (.not. ieee_is_finite(sampled_raw)) sampled_raw = huge(1.0_dp)
            if (.not. valid_value(sampled_raw)) return
            if (has_constraints) then
                value = -2.0_dp + (2.0_dp/acos(-1.0_dp))*atan(sampled_raw)
            else
                value = sampled_raw
            end if
            if (sampled_raw < res%fun .or. .not. valid_value(res%fun)) then
                res%fun = sampled_raw
                res%x = x
                res%maxcv = cv
                res%feasible = .true.
            end if
        end function counted_objective

        function local_objective(x) result(value)
            real(dp), intent(in) :: x(:) !! feasible box point (constraints may be violated)
            real(dp) :: value !! raw objective for constrained SLSQP or other polisher
            real(dp) :: cv
            value = huge(1.0_dp)
            if (size(x) /= dim) return
            if (.not. all(ieee_is_finite(x))) return
            if (any(x < bounds(:,1)) .or. any(x > bounds(:,2))) return
            if (res%nfev >= eval_limit) then
                local_budget_hit = .true.
                return
            end if
            res%nlfev = res%nlfev + 1
            if (.not. has_constraints) then
                value = counted_objective(x)
                return
            end if
            call evaluate_violation(x,cv)
            call record_violation(x,cv)
            res%nfev = res%nfev + 1
            value = fun(x)
            if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
            if (cv <= cv_tol .and. valid_value(value)) then
                if (value < res%fun .or. .not. valid_value(res%fun)) then
                    res%fun = value
                    res%x = x
                    res%maxcv = cv
                    res%feasible = .true.
                end if
            end if
        end function local_objective

        logical function target_reached() result(hit)
            real(dp) :: scale, tolerance
            hit = has_target .and. valid_value(res%fun)
            if (.not. hit) return
            scale = max(1.0_dp,abs(f_target))
            if (target_tol > huge(1.0_dp)/scale) return
            tolerance = target_tol*scale
            if (f_target > 0.0_dp) then
                if (tolerance > huge(1.0_dp)-f_target) return
            end if
            hit = res%fun <= f_target+tolerance
        end function target_reached

        logical function no_material_improvement(previous_best, current_best) result(stable)
            real(dp), intent(in) :: previous_best !! best finite value before current round
            real(dp), intent(in) :: current_best !! best finite value after local polishing
            real(dp) :: scale
            stable = .false.
            if (.not. valid_value(previous_best)) return
            if (.not. valid_value(current_best)) return
            ! Normalize before subtraction to avoid overflow on wide ranges.
            scale = max(1.0_dp,abs(previous_best),abs(current_best))
            stable = previous_best/scale-current_best/scale <= hgr_tol
        end function no_material_improvement

        subroutine add_minimum(x, value)
            real(dp), intent(in) :: x(:) !! candidate local minimum
            real(dp), intent(in) :: value !! finite candidate function value
            real(dp), allocatable :: grown_points(:,:), grown_values(:)
            integer :: l, new_capacity
            if (.not. valid_value(value)) return
            do l = 1, nlocal_stored
                if (maxval(abs((x-local_points(:,l))/widths)) <= duplicate_tol) then
                    if (value < local_values(l)) then
                        local_points(:,l) = x
                        local_values(l) = value
                    end if
                    return
                end if
            end do
            if (nlocal_stored == local_capacity) then
                new_capacity = min(eval_limit, max(2*local_capacity,1))
                if (new_capacity <= local_capacity) return
                allocate(grown_points(dim,new_capacity),grown_values(new_capacity))
                grown_points(:,1:nlocal_stored) = local_points(:,1:nlocal_stored)
                grown_values(1:nlocal_stored) = local_values(1:nlocal_stored)
                call move_alloc(grown_points,local_points)
                call move_alloc(grown_values,local_values)
                local_capacity = new_capacity
            end if
            nlocal_stored = nlocal_stored + 1
            local_points(:,nlocal_stored) = x
            local_values(nlocal_stored) = value
        end subroutine add_minimum

    end function shgo_core

    subroutine shgo_solve_dense(matrix, rhs, x, ok)
        ! Pivoted Gaussian elimination for the small (nvar <= 8) damped
        ! normal equations used only during optional feasibility restoration.
        real(dp), intent(in) :: matrix(:,:) !! symmetric regularized coefficient matrix
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! solution, or zero on failure
        logical, intent(out) :: ok !! successful finite linear solve
        real(dp) :: a(size(rhs),size(rhs)), b(size(rhs)), tmp(size(rhs)), pivot, factor
        integer :: n, i, k, r
        n = size(rhs)
        a = matrix
        b = rhs
        x = 0.0_dp
        ok = .false.
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        do k = 1, n
            r = k-1+maxloc(abs(a(k:n,k)),dim=1)
            pivot = abs(a(r,k))
            if (pivot <= tiny(1.0_dp)) return
            if (r /= k) then
                tmp = a(k,:)
                a(k,:) = a(r,:)
                a(r,:) = tmp
                pivot = b(k)
                b(k) = b(r)
                b(r) = pivot
            end if
            do i = k+1, n
                factor = a(i,k)/a(k,k)
                a(i,k:n) = a(i,k:n)-factor*a(k,k:n)
                b(i) = b(i)-factor*b(k)
            end do
        end do
        do i = n, 1, -1
            x(i) = (b(i)-sum(a(i,i+1:n)*x(i+1:n)))/a(i,i)
        end do
        ok = all(ieee_is_finite(x))
    end subroutine shgo_solve_dense

    subroutine setup_adaptive_mesh(mesh, dim, cap)
        type(adaptive_mesh), intent(out) :: mesh !! new vertex/leaf store
        integer, intent(in) :: dim !! optimization dimension
        integer, intent(in) :: cap !! upper bound on sampled vertices
        mesh%dim = dim
        mesh%cap = cap
        allocate(mesh%keys(dim,cap),mesh%origins(dim,cap),mesh%spans(cap))
        allocate(mesh%points(dim,cap),mesh%values(cap),mesh%raw(cap))
        allocate(mesh%slots(2*cap+1))
        mesh%slots = 0
    end subroutine setup_adaptive_mesh

    pure integer function adaptive_bucket(mesh, key) result(bucket)
        type(adaptive_mesh), intent(in) :: mesh !! initialized hash table
        integer(int64), intent(in) :: key(:) !! dyadic lattice coordinates
        integer(int64) :: hash, modulus
        integer :: j
        modulus = int(size(mesh%slots),int64)
        hash = 0_int64
        do j = 1, size(key)
            hash = modulo(hash*1009_int64 + modulo(key(j),modulus) + &
                int(j,int64), modulus)
        end do
        bucket = int(hash)+1
    end function adaptive_bucket

    pure integer function lookup_adaptive_vertex(mesh, key) result(index)
        type(adaptive_mesh), intent(in) :: mesh !! initialized vertex store
        integer(int64), intent(in) :: key(:) !! exact dyadic vertex key
        integer :: bucket, stored, i
        bucket = adaptive_bucket(mesh,key)
        index = 0
        do i = 1, size(mesh%slots)
            stored = mesh%slots(bucket)
            if (stored == 0) return
            if (all(mesh%keys(:,stored) == key)) then
                index = stored
                return
            end if
            bucket = mod(bucket,size(mesh%slots))+1
        end do
    end function lookup_adaptive_vertex

    subroutine insert_adaptive_vertex(mesh, key, point, value, raw)
        type(adaptive_mesh), intent(inout) :: mesh !! vertex cache
        integer(int64), intent(in) :: key(:) !! unique dyadic lattice coordinate
        real(dp), intent(in) :: point(:) !! physical coordinate
        real(dp), intent(in) :: value !! vertex ranking score
        real(dp), intent(in), optional :: raw !! actual feasible objective or huge sentinel
        integer :: bucket
        if (mesh%npoint >= mesh%cap) error stop 'SHGO: adaptive vertex capacity'
        bucket = adaptive_bucket(mesh,key)
        do while (mesh%slots(bucket) /= 0)
            if (all(mesh%keys(:,mesh%slots(bucket)) == key)) then
                error stop 'SHGO: duplicate adaptive vertex insertion'
            end if
            bucket = mod(bucket,size(mesh%slots))+1
        end do
        mesh%npoint = mesh%npoint+1
        mesh%keys(:,mesh%npoint) = key
        mesh%points(:,mesh%npoint) = point
        mesh%values(mesh%npoint) = value
        mesh%raw(mesh%npoint) = value
        if (present(raw)) mesh%raw(mesh%npoint) = raw
        mesh%slots(bucket) = mesh%npoint
    end subroutine insert_adaptive_vertex

    subroutine stage_adaptive_vertex(mesh,key,staged,count)
        ! Deduplicate prospective samples against both the committed mesh and
        ! other vertices in this as-yet-uncommitted split. No objective is
        ! called until the whole split has passed the budget/memory preflight.
        type(adaptive_mesh), intent(in) :: mesh !! committed vertex cache
        integer(int64), intent(in) :: key(:) !! prospective lattice coordinate
        integer(int64), intent(inout) :: staged(:,:) !! pending unique coordinates
        integer, intent(inout) :: count !! pending vertex count
        integer :: i
        if (lookup_adaptive_vertex(mesh,key) /= 0) return
        do i = 1, count
            if (all(staged(:,i) == key)) return
        end do
        if (count >= size(staged,2)) error stop 'SHGO: staging capacity'
        count = count + 1
        staged(:,count) = key
    end subroutine stage_adaptive_vertex

    subroutine build_conforming_tetrahedra(mesh, mesh_limit, tetrahedra, valid, over_limit)
        ! A common dyadic subdivision is computed for each geometric face.
        ! All cached edge/face vertices are reused from either side of a
        ! coarse-fine interface. Boundary polygons are triangulated using a
        ! sampled face center, or exact-integer ear clipping when a recursive
        ! subface has no center sample. Each surface triangle is coned to the
        ! sampled interior center of its owning leaf cube.
        type(adaptive_mesh), intent(in) :: mesh !! three-dimensional cached adaptive leaves
        integer, intent(in) :: mesh_limit !! caller's vertex memory limit
        integer, allocatable, intent(out) :: tetrahedra(:,:) !! vertex ids, shape (4,ntetra)
        logical, intent(out) :: valid !! true if all faces/volumes and incidences are valid
        logical, intent(out) :: over_limit !! exceeded the implementation simplex safety limit
        integer, allocatable :: polygon(:), working(:)
        integer, allocatable :: face_a(:), face_b(:), face_c(:), face_count(:)
        integer :: leaf, axis, side, uaxis, vaxis, center, phase, total
        integer :: tetra_limit, t, f, nslots, slot, count_exterior
        integer :: ids(3), a, b, c, tmp
        integer(int64) :: low(3), center_key(3), face_area2, span
        integer(int64) :: d1(3), d2(3), d3(3), det, extent
        logical :: ok

        valid = .false.
        over_limit = .false.
        if (mesh%dim /= 3 .or. mesh%ncell < 1 .or. mesh%npoint < 1) return
        tetra_limit = int(min(300000_int64,24_int64*int(mesh_limit,int64)))
        allocate(polygon(mesh%npoint),working(mesh%npoint))
        extent = 2_int64**lattice_level
        ok = .true.
        ! The first pass counts simplices before allocating their connectivity.
        ! The second pass validates their signed determinants and face areas.
        do phase = 1, 2
            total = 0
            if (phase == 2) allocate(tetrahedra(4,t))
            do leaf = 1, mesh%ncell
                center_key = mesh%origins(:,leaf)+mesh%spans(leaf)/2_int64
                center = lookup_adaptive_vertex(mesh,center_key)
                if (center == 0) return
                span = mesh%spans(leaf)
                do axis = 1, 3
                    uaxis = mod(axis,3)+1
                    vaxis = mod(axis+1,3)+1
                    do side = 0, 1
                        low = mesh%origins(:,leaf)
                        low(axis) = low(axis) + int(side,int64)*span
                        face_area2 = 0_int64
                        call face_patch(low,span)
                        if (.not. ok) return
                        if (face_area2 /= 2_int64*span*span) return
                    end do
                end do
            end do
            if (phase == 1) t = total
        end do
        if (total /= size(tetrahedra,2)) return

        ! Each triangle face must have either two incident tetrahedra or be
        ! an exterior boundary face of the original rectangular domain.
        nslots = 8*total+1
        allocate(face_a(nslots),face_b(nslots),face_c(nslots),face_count(nslots))
        face_count = 0
        do t = 1, total
            do f = 1, 4
                ! Face opposite vertex f.
                select case (f)
                case (1)
                    ids = tetrahedra(2:4,t)
                case (2)
                    ids = tetrahedra([1,3,4],t)
                case (3)
                    ids = tetrahedra([1,2,4],t)
                case (4)
                    ids = tetrahedra(1:3,t)
                end select
                if (ids(2) < ids(1)) then
                    tmp = ids(1); ids(1) = ids(2); ids(2) = tmp
                end if
                if (ids(3) < ids(2)) then
                    tmp = ids(2); ids(2) = ids(3); ids(3) = tmp
                end if
                if (ids(2) < ids(1)) then
                    tmp = ids(1); ids(1) = ids(2); ids(2) = tmp
                end if
                slot = int(modulo(int(ids(1),int64)*1000003_int64 + &
                    int(ids(2),int64)*9176_int64 + int(ids(3),int64)*13007_int64, &
                    int(nslots,int64)))+1
                do while (face_count(slot) /= 0)
                    if (face_a(slot) == ids(1) .and. face_b(slot) == ids(2) .and. &
                            face_c(slot) == ids(3)) exit
                    slot = mod(slot,nslots)+1
                end do
                if (face_count(slot) == 0) then
                    face_a(slot) = ids(1)
                    face_b(slot) = ids(2)
                    face_c(slot) = ids(3)
                end if
                face_count(slot) = face_count(slot)+1
                if (face_count(slot) > 2) return
            end do
        end do
        count_exterior = 0
        do slot = 1, nslots
            if (face_count(slot) /= 1) cycle
            a = face_a(slot); b = face_b(slot); c = face_c(slot)
            ok = .false.
            do axis = 1, 3
                if ((mesh%keys(axis,a) == 0_int64 .and. &
                        mesh%keys(axis,b) == 0_int64 .and. mesh%keys(axis,c) == 0_int64) .or. &
                        (mesh%keys(axis,a) == extent .and. &
                        mesh%keys(axis,b) == extent .and. mesh%keys(axis,c) == extent)) ok = .true.
            end do
            if (.not. ok) return
            count_exterior = count_exterior+1
        end do
        valid = count_exterior >= 12

    contains

        recursive subroutine face_patch(patch_low,patch_span)
            integer(int64), intent(in) :: patch_low(3) !! low u/v corner, face axis fixed
            integer(int64), intent(in) :: patch_span !! dyadic patch side length
            integer(int64) :: mid(3), child_low(3), half
            integer :: i, j, nbound, face_center, first, last
            logical :: can_split
            if (.not. ok) return
            half = patch_span/2_int64
            if (half >= 1_int64 .and. mod(patch_span,2_int64) == 0_int64) then
                can_split = .true.
                ! A full 2-by-2 split requires all four edge midpoints and
                ! the patch center. If only some edge midpoints exist, the
                ! polygon fan/ear-cut below handles them as hanging vertices.
                do i = 1, 5
                    mid = patch_low
                    select case (i)
                    case (1)
                        mid(uaxis) = mid(uaxis)+half
                    case (2)
                        mid(uaxis) = mid(uaxis)+patch_span
                        mid(vaxis) = mid(vaxis)+half
                    case (3)
                        mid(uaxis) = mid(uaxis)+half
                        mid(vaxis) = mid(vaxis)+patch_span
                    case (4)
                        mid(vaxis) = mid(vaxis)+half
                    case (5)
                        mid(uaxis) = mid(uaxis)+half
                        mid(vaxis) = mid(vaxis)+half
                    end select
                    if (lookup_adaptive_vertex(mesh,mid) == 0) can_split = .false.
                end do
                if (can_split) then
                    do j = 0, 1
                        do i = 0, 1
                            child_low = patch_low
                            child_low(uaxis) = child_low(uaxis)+int(i,int64)*half
                            child_low(vaxis) = child_low(vaxis)+int(j,int64)*half
                            call face_patch(child_low,half)
                            if (.not. ok) return
                        end do
                    end do
                    return
                end if
            end if
            call patch_boundary(patch_low,patch_span,nbound)
            if (.not. ok .or. nbound < 4) then
                ok = .false.
                return
            end if
            mid = patch_low
            mid(uaxis) = mid(uaxis)+half
            mid(vaxis) = mid(vaxis)+half
            face_center = 0
            if (2_int64*half == patch_span) face_center = lookup_adaptive_vertex(mesh,mid)
            if (face_center > 0) then
                do i = 1, nbound
                    call emit_face_triangle(face_center,polygon(i),polygon(mod(i,nbound)+1))
                    if (.not. ok) return
                end do
            else
                ! Ordinary polygon triangulation when the patch is only a
                ! portion of a larger leaf face and has no sampled center.
                working(1:nbound) = polygon(1:nbound)
                first = nbound
                do while (first > 3)
                    last = 0
                    do i = 1, first
                        if (.not. is_ear(working(1:first),i)) cycle
                        last = i
                        exit
                    end do
                    if (last == 0) then
                        ok = .false.
                        return
                    end if
                    call emit_face_triangle(working(mod(last+first-2,first)+1), &
                        working(last),working(mod(last,first)+1))
                    if (.not. ok) return
                    working(last:first-1) = working(last+1:first)
                    first = first-1
                end do
                call emit_face_triangle(working(1),working(2),working(3))
            end if
        end subroutine face_patch

        subroutine patch_boundary(patch_low,patch_span,nbound)
            integer(int64), intent(in) :: patch_low(3) !! low corner on face
            integer(int64), intent(in) :: patch_span !! side length
            integer, intent(out) :: nbound !! count of oriented boundary vertices
            integer(int64) :: uend(3), vend(3), diagonal(3)
            uend = patch_low
            vend = patch_low
            uend(uaxis) = uend(uaxis)+patch_span
            vend(vaxis) = vend(vaxis)+patch_span
            diagonal = uend
            diagonal(vaxis) = diagonal(vaxis)+patch_span
            nbound = 0
            call trace_segment(patch_low,uend,nbound)
            call trace_segment(uend,diagonal,nbound)
            call trace_segment(diagonal,vend,nbound)
            call trace_segment(vend,patch_low,nbound)
        end subroutine patch_boundary

        recursive subroutine trace_segment(first,last,nbound)
            integer(int64), intent(in) :: first(3) !! inclusive first lattice key
            integer(int64), intent(in) :: last(3) !! excluded final key
            integer, intent(inout) :: nbound !! boundary point count
            integer(int64) :: middle(3)
            integer :: found
            if (.not. ok) return
            middle = (first+last)/2_int64
            if (all(2_int64*middle == first+last)) then
                found = lookup_adaptive_vertex(mesh,middle)
                if (found > 0) then
                    call trace_segment(first,middle,nbound)
                    call trace_segment(middle,last,nbound)
                    return
                end if
            end if
            found = lookup_adaptive_vertex(mesh,first)
            if (found <= 0 .or. nbound >= size(polygon)) then
                ok = .false.
                return
            end if
            nbound = nbound+1
            polygon(nbound) = found
        end subroutine trace_segment

        pure integer(int64) function orient(a,b,c) result(cross)
            integer, intent(in) :: a,b,c !! cached face vertex indices
            integer(int64) :: ab_u, ab_v, ac_u, ac_v
            ab_u = mesh%keys(uaxis,b)-mesh%keys(uaxis,a)
            ab_v = mesh%keys(vaxis,b)-mesh%keys(vaxis,a)
            ac_u = mesh%keys(uaxis,c)-mesh%keys(uaxis,a)
            ac_v = mesh%keys(vaxis,c)-mesh%keys(vaxis,a)
            cross = ab_u*ac_v-ab_v*ac_u
        end function orient

        logical function is_ear(boundary,i) result(yes)
            integer, intent(in) :: boundary(:) !! oriented polygon ring
            integer, intent(in) :: i !! proposed ear tip
            integer :: m, prev, next, j, a, b, c, p
            m = size(boundary)
            prev = mod(i+m-2,m)+1
            next = mod(i,m)+1
            a = boundary(prev); b = boundary(i); c = boundary(next)
            yes = orient(a,b,c) > 0_int64
            if (.not. yes) return
            do j = 1, m
                if (j == prev .or. j == i .or. j == next) cycle
                p = boundary(j)
                if (orient(a,b,p) >= 0_int64 .and. orient(b,c,p) >= 0_int64 .and. &
                        orient(c,a,p) >= 0_int64) then
                    yes = .false.
                    return
                end if
            end do
        end function is_ear

        subroutine emit_face_triangle(a,b,c)
            integer, intent(in) :: a,b,c !! CCW face vertex ids in local (u,v)
            integer(int64) :: area2
            integer :: pa,pb,pc
            if (.not. ok) return
            area2 = orient(a,b,c)
            if (area2 <= 0_int64) then
                ok = .false.
                return
            end if
            face_area2 = face_area2+area2
            total = total+1
            if (total > tetra_limit) then
                over_limit = .true.
                ok = .false.
                return
            end if
            if (phase /= 2) return
            pa = a; pb = b; pc = c
            if (side == 0) then
                pb = c; pc = b
            end if
            tetrahedra(:,total) = [center,pa,pb,pc]
            d1 = mesh%keys(:,pa)-mesh%keys(:,center)
            d2 = mesh%keys(:,pb)-mesh%keys(:,center)
            d3 = mesh%keys(:,pc)-mesh%keys(:,center)
            det = d1(1)*(d2(2)*d3(3)-d2(3)*d3(2)) - &
                d1(2)*(d2(1)*d3(3)-d2(3)*d3(1)) + &
                d1(3)*(d2(1)*d3(2)-d2(2)*d3(1))
            if (det /= span/2_int64*area2 .or. det <= 0_int64) ok = .false.
        end subroutine emit_face_triangle

    end subroutine build_conforming_tetrahedra

    subroutine find_tetrahedral_candidates(mesh,tetrahedra,candidates,ncand)
        ! Reuse the plateau-aware triangular 1-skeleton implementation by
        ! representing each tetrahedron as its four triangular faces.
        type(adaptive_mesh), intent(in) :: mesh !! cached simplex vertices
        integer, intent(in) :: tetrahedra(:,:) !! four vertex indices per simplex
        integer, intent(out) :: candidates(:) !! plateau representatives
        integer, intent(out) :: ncand !! number of sampled graph minima
        integer, allocatable :: faces(:,:)
        integer :: t, offset
        allocate(faces(3,4*size(tetrahedra,2)))
        do t = 1, size(tetrahedra,2)
            offset = 4*(t-1)
            faces(:,offset+1) = tetrahedra(2:4,t)
            faces(:,offset+2) = tetrahedra([1,3,4],t)
            faces(:,offset+3) = tetrahedra([1,2,4],t)
            faces(:,offset+4) = tetrahedra(1:3,t)
        end do
        call find_triangle_candidates(mesh,faces,candidates,ncand)
    end subroutine find_tetrahedral_candidates

    subroutine build_conforming_triangles(mesh, triangles, valid)
        ! Fan triangulation of every adaptive square. Recursively finding
        ! dyadic edge midpoints in the vertex cache ensures that a coarse
        ! cell uses exactly the same face segments as its refined neighbors.
        type(adaptive_mesh), intent(in) :: mesh !! two-dimensional leaf mesh and cached samples
        integer, allocatable, intent(out) :: triangles(:,:) !! vertex ids, shape (3,ntriangles)
        logical, intent(out) :: valid !! true only for a connected, conforming planar mesh
        integer, allocatable :: polygon(:), edge_u(:), edge_v(:), edge_count(:)
        integer :: leaf, nbound, total, t, i, center, nslots, slot, u, v
        integer :: a, b, c, e, nboundary
        integer(int64) :: corner_a(2), corner_b(2), center_key(2)
        integer(int64) :: area2, triangle_area2, extent
        logical :: boundary_ok
        valid = .false.
        boundary_ok = .true.
        if (mesh%dim /= 2) return
        if (mesh%npoint < 1 .or. mesh%ncell < 1) return
        allocate(polygon(mesh%npoint))
        total = 0
        do leaf = 1, mesh%ncell
            call leaf_boundary(leaf, polygon, nbound)
            if (.not. boundary_ok) return
            if (nbound < 4 .or. nbound > huge(total)-total) return
            total = total + nbound
        end do
        if (total > (huge(total)-1)/6) return
        allocate(triangles(3,total))
        t = 0
        area2 = 0_int64
        do leaf = 1, mesh%ncell
            call leaf_boundary(leaf, polygon, nbound)
            if (.not. boundary_ok) return
            center_key = mesh%origins(:,leaf)+mesh%spans(leaf)/2_int64
            center = lookup_adaptive_vertex(mesh,center_key)
            if (center == 0) return
            do i = 1, nbound
                t = t+1
                triangles(1,t) = center
                triangles(2,t) = polygon(i)
                triangles(3,t) = polygon(mod(i,nbound)+1)
                a = triangles(1,t)
                b = triangles(2,t)
                c = triangles(3,t)
                corner_a = mesh%keys(:,b)-mesh%keys(:,a)
                corner_b = mesh%keys(:,c)-mesh%keys(:,a)
                triangle_area2 = corner_a(1)*corner_b(2)-corner_a(2)*corner_b(1)
                if (triangle_area2 <= 0_int64) return
                area2 = area2+triangle_area2
            end do
        end do
        extent = 2_int64**lattice_level
        if (area2 /= 2_int64*extent*extent) return

        ! Every simplex edge must occur twice in the interior, once on the
        ! exterior rectangle, and never three or more times. This detects
        ! missing hanging-node subdivisions as well as overlaps/gaps.
        nslots = 6*total+1
        allocate(edge_u(nslots),edge_v(nslots),edge_count(nslots))
        edge_count = 0
        do t = 1, total
            do e = 1, 3
                a = triangles(e,t)
                b = triangles(mod(e,3)+1,t)
                u = min(a,b)
                v = max(a,b)
                slot = int(modulo(int(u,int64)*1000003_int64+ &
                    int(v,int64)*9176_int64,int(nslots,int64)))+1
                do while (edge_count(slot) /= 0)
                    if (edge_u(slot) == u .and. edge_v(slot) == v) exit
                    slot = mod(slot,nslots)+1
                end do
                if (edge_count(slot) == 0) then
                    edge_u(slot) = u
                    edge_v(slot) = v
                end if
                edge_count(slot) = edge_count(slot)+1
                if (edge_count(slot) > 2) return
            end do
        end do
        nboundary = 0
        do slot = 1, nslots
            if (edge_count(slot) /= 1) cycle
            a = edge_u(slot)
            b = edge_v(slot)
            if (.not. ((mesh%keys(1,a) == 0_int64 .and. &
                    mesh%keys(1,b) == 0_int64) .or. &
                    (mesh%keys(1,a) == extent .and. &
                    mesh%keys(1,b) == extent) .or. &
                    (mesh%keys(2,a) == 0_int64 .and. &
                    mesh%keys(2,b) == 0_int64) .or. &
                    (mesh%keys(2,a) == extent .and. &
                    mesh%keys(2,b) == extent))) return
            nboundary = nboundary+1
        end do
        valid = nboundary >= 4
    contains
        subroutine leaf_boundary(cell, ids, nids)
            integer, intent(in) :: cell !! index of the current adaptive leaf
            integer, intent(out) :: ids(:) !! counterclockwise boundary vertex ids
            integer, intent(out) :: nids !! number of boundary vertices
            integer(int64) :: low(2), high(2), akey(2), bkey(2)
            integer :: side
            low = mesh%origins(:,cell)
            high = low + mesh%spans(cell)
            nids = 0
            do side = 1, 4
                select case (side)
                case (1)
                    akey = low
                    bkey = [high(1),low(2)]
                case (2)
                    akey = [high(1),low(2)]
                    bkey = high
                case (3)
                    akey = high
                    bkey = [low(1),high(2)]
                case (4)
                    akey = [low(1),high(2)]
                    bkey = low
                end select
                call trace_edge(akey,bkey,ids,nids)
            end do
        end subroutine leaf_boundary

        recursive subroutine trace_edge(first,last,ids,nids)
            integer(int64), intent(in) :: first(2) !! inclusive start lattice point
            integer(int64), intent(in) :: last(2) !! excluded end lattice point
            integer, intent(inout) :: ids(:) !! polygon vertex ids being appended
            integer, intent(inout) :: nids !! number of entries written
            integer(int64) :: middle(2)
            integer :: vertex
            if (.not. boundary_ok) return
            middle = (first+last)/2_int64
            if (all(2_int64*middle == first+last)) then
                vertex = lookup_adaptive_vertex(mesh,middle)
                if (vertex > 0) then
                    call trace_edge(first,middle,ids,nids)
                    call trace_edge(middle,last,ids,nids)
                    return
                end if
            end if
            vertex = lookup_adaptive_vertex(mesh,first)
            if (vertex <= 0 .or. nids >= size(ids)) then
                boundary_ok = .false.
                return
            end if
            nids = nids+1
            ids(nids) = vertex
        end subroutine trace_edge
    end subroutine build_conforming_triangles

    subroutine find_triangle_candidates(mesh, triangles, candidates, ncand)
        ! The conforming simplicial 1-skeleton uses every triangle edge.
        ! Equal-valued connected vertices form a single graph minimum.
        type(adaptive_mesh), intent(in) :: mesh !! sampled simplicial mesh
        integer, intent(in) :: triangles(:,:) !! triangular vertex connectivity
        integer, intent(out) :: candidates(:) !! minima representative vertex indices
        integer, intent(out) :: ncand !! number of distinct minimum plateaus
        integer, allocatable :: parent(:), height(:), representative(:)
        logical, allocatable :: has_lower_neighbor(:)
        integer :: i, t, a, b, root, edge
        allocate(parent(mesh%npoint),height(mesh%npoint))
        allocate(representative(mesh%npoint),has_lower_neighbor(mesh%npoint))
        do i = 1, mesh%npoint
            parent(i) = i
        end do
        height = 0
        representative = 0
        has_lower_neighbor = .false.
        do t = 1, size(triangles,2)
            do edge = 1, 3
                a = triangles(edge,t)
                b = triangles(mod(edge,3)+1,t)
                if (mesh%values(a) == mesh%values(b)) call unite(a,b)
            end do
        end do
        do t = 1, size(triangles,2)
            do edge = 1, 3
                a = triangles(edge,t)
                b = triangles(mod(edge,3)+1,t)
                if (mesh%values(a) > mesh%values(b)) then
                    has_lower_neighbor(find_root(a)) = .true.
                else if (mesh%values(b) > mesh%values(a)) then
                    has_lower_neighbor(find_root(b)) = .true.
                end if
            end do
        end do
        do i = 1, mesh%npoint
            if (.not. valid_value(mesh%values(i))) cycle
            root = find_root(i)
            if (has_lower_neighbor(root)) cycle
            if (representative(root) == 0) representative(root) = i
        end do
        ncand = 0
        do i = 1, mesh%npoint
            if (parent(i) /= i) cycle
            if (representative(i) == 0) cycle
            ncand = ncand+1
            candidates(ncand) = representative(i)
        end do
    contains
        integer function find_root(vertex) result(root)
            integer, intent(in) :: vertex !! graph vertex index
            integer :: current, next_vertex
            root = vertex
            do while (parent(root) /= root)
                root = parent(root)
            end do
            current = vertex
            do while (parent(current) /= root)
                next_vertex = parent(current)
                parent(current) = root
                current = next_vertex
            end do
        end function find_root
        subroutine unite(vertex_a, vertex_b)
            integer, intent(in) :: vertex_a !! first graph vertex
            integer, intent(in) :: vertex_b !! second graph vertex
            integer :: root_a, root_b
            root_a = find_root(vertex_a)
            root_b = find_root(vertex_b)
            if (root_a == root_b) return
            if (height(root_a) < height(root_b)) then
                parent(root_a) = root_b
            else if (height(root_a) > height(root_b)) then
                parent(root_b) = root_a
            else
                parent(root_b) = root_a
                height(root_a) = height(root_a)+1
            end if
        end subroutine unite
    end subroutine find_triangle_candidates

    subroutine find_adaptive_candidates(mesh, candidates, ncand)
        type(adaptive_mesh), intent(in) :: mesh !! leaf-cell vertex-edge graph
        integer, intent(out) :: candidates(:) !! qualifying vertex indices
        integer, intent(out) :: ncand !! number of graph local minima
        integer, allocatable :: parent(:), height(:), representative(:)
        logical, allocatable :: has_lower_neighbor(:)
        integer :: ids(2**mesh%dim), leaf, a, b, j, mask, va, vb, idx, root
        integer(int64) :: key(mesh%dim)

        ! Equal-value plateaus are single connected components in the
        ! sublevel graph. Vertex-insertion order is not a valid tie breaker:
        ! refining cells removes some old edges and could otherwise create
        ! spurious minima on a constant function.
        allocate(parent(mesh%npoint),height(mesh%npoint))
        allocate(representative(mesh%npoint),has_lower_neighbor(mesh%npoint))
        do idx = 1, mesh%npoint
            parent(idx) = idx
        end do
        height = 0
        representative = 0
        has_lower_neighbor = .false.
        ! First identify connected plateaus of equal-valued vertices.
        do leaf = 1, mesh%ncell
            do mask = 0, 2**mesh%dim-1
                key = mesh%origins(:,leaf)
                do j = 1, mesh%dim
                    if (btest(mask,j-1)) key(j) = key(j)+mesh%spans(leaf)
                end do
                ids(mask+1) = lookup_adaptive_vertex(mesh,key)
                if (ids(mask+1) == 0) error stop 'SHGO: incomplete adaptive mesh'
            end do
            ! Kuhn simplex edges: vertex bit masks are comparable under
            ! coordinatewise order, with deterministic tie breaking.
            do a = 0, 2**mesh%dim-2
                va = ids(a+1)
                do b = a+1, 2**mesh%dim-1
                    if (iand(a,b) /= a .and. iand(a,b) /= b) cycle
                    vb = ids(b+1)
                    if (mesh%values(va) == mesh%values(vb)) call unite(va,vb)
                end do
            end do
        end do

        ! A plateau is a local minimum only if none of its boundary edges
        ! goes downhill. A single representative is returned per plateau.
        do leaf = 1, mesh%ncell
            do mask = 0, 2**mesh%dim-1
                key = mesh%origins(:,leaf)
                do j = 1, mesh%dim
                    if (btest(mask,j-1)) key(j) = key(j)+mesh%spans(leaf)
                end do
                ids(mask+1) = lookup_adaptive_vertex(mesh,key)
            end do
            do a = 0, 2**mesh%dim-2
                va = ids(a+1)
                do b = a+1, 2**mesh%dim-1
                    if (iand(a,b) /= a .and. iand(a,b) /= b) cycle
                    vb = ids(b+1)
                    if (mesh%values(va) > mesh%values(vb)) then
                        has_lower_neighbor(find_root(va)) = .true.
                    else if (mesh%values(vb) > mesh%values(va)) then
                        has_lower_neighbor(find_root(vb)) = .true.
                    end if
                end do
            end do
        end do
        do idx = 1, mesh%npoint
            if (.not. valid_value(mesh%values(idx))) cycle
            root = find_root(idx)
            if (has_lower_neighbor(root)) cycle
            if (representative(root) == 0) representative(root) = idx
        end do
        ncand = 0
        do idx = 1, mesh%npoint
            if (parent(idx) /= idx) cycle
            if (representative(idx) == 0) cycle
            ncand = ncand+1
            candidates(ncand) = representative(idx)
        end do

    contains

        integer function find_root(vertex) result(root)
            integer, intent(in) :: vertex !! vertex index in the plateau forest
            integer :: current, next_vertex
            root = vertex
            do while (parent(root) /= root)
                root = parent(root)
            end do
            current = vertex
            do while (parent(current) /= root)
                next_vertex = parent(current)
                parent(current) = root
                current = next_vertex
            end do
        end function find_root

        subroutine unite(vertex_a, vertex_b)
            integer, intent(in) :: vertex_a !! first equal-value graph vertex
            integer, intent(in) :: vertex_b !! second equal-value graph vertex
            integer :: root_a, root_b
            root_a = find_root(vertex_a)
            root_b = find_root(vertex_b)
            if (root_a == root_b) return
            if (height(root_a) < height(root_b)) then
                parent(root_a) = root_b
            else if (height(root_a) > height(root_b)) then
                parent(root_b) = root_a
            else
                parent(root_b) = root_a
                height(root_a) = height(root_a)+1
            end if
        end subroutine unite
    end subroutine find_adaptive_candidates

    pure function capped_grid_size(q, dim, cap) result(total)
        integer, intent(in) :: q !! number of grid points per dimension
        integer, intent(in) :: dim !! number of dimensions
        integer(int64), intent(in) :: cap !! maximum permitted product
        integer(int64) :: total !! true vertex count if <= cap, otherwise cap+1
        integer :: j
        total = 1_int64
        do j = 1, dim
            if (total > cap/int(q,int64)) then
                total = cap + 1_int64
                return
            end if
            total = total * int(q,int64)
        end do
    end function capped_grid_size

    pure logical function valid_value(value) result(ok)
        real(dp), intent(in) :: value !! objective sample
        ok = .false.
        ! Fortran does not require short-circuit .and.; never compare a NaN.
        if (.not. ieee_is_finite(value)) return
        ok = value < huge(1.0_dp)
    end function valid_value

    pure subroutine point_from_index(index, q, bounds, widths, x)
        integer, intent(in) :: index !! one-based flattened grid position
        integer, intent(in) :: q !! points along each coordinate axis, >= 3
        real(dp), intent(in) :: bounds(:,:) !! search domain, shape (n,2)
        real(dp), intent(in) :: widths(:) !! interval lengths, shape (n)
        real(dp), intent(out) :: x(:) !! physical point, shape (n)
        integer :: j, flat, digit
        flat = index-1
        do j = 1, size(x)
            digit = mod(flat,q)
            flat = flat/q
            x(j) = bounds(j,1) + (real(digit,dp)/real(q-1,dp))*widths(j)
        end do
    end subroutine point_from_index

    pure logical function is_graph_minimum(index, q, dim, values) result(is_min)
        integer, intent(in) :: index !! one-based flattened grid vertex
        integer, intent(in) :: q !! odd grid width, >= 3
        integer, intent(in) :: dim !! number of axes in [1,8]
        real(dp), intent(in) :: values(:) !! values at all grid vertices
        integer :: digits(dim), steps(dim), flat, j, stride, mask, direction, other
        logical :: in_grid
        is_min = .true.
        flat = index-1
        stride = 1
        do j = 1, dim
            digits(j) = mod(flat,q)
            steps(j) = stride
            flat = flat/q
            stride = stride*q
        end do
        ! Kuhn simplex edges join comparable coordinatewise vertices whose
        ! coordinate differences are 0 or all +1 or all -1.
        do direction = -1, 1, 2
            do mask = 1, 2**dim-1
                other = index
                in_grid = .true.
                do j = 1, dim
                    if (.not. btest(mask,j-1)) cycle
                    if (digits(j)+direction < 0 .or. digits(j)+direction >= q) then
                        in_grid = .false.
                        exit
                    end if
                    other = other + direction * steps(j)
                end do
                if (.not. in_grid) cycle
                if (values(other) < values(index) .or. &
                        (values(other) == values(index) .and. other < index)) then
                    is_min = .false.
                    return
                end if
            end do
        end do
    end function is_graph_minimum

    recursive subroutine sort_candidates(indices, count, values)
        integer, intent(inout) :: indices(:) !! candidate indices to sort in place
        integer, intent(in) :: count !! number of valid indices
        real(dp), intent(in) :: values(:) !! vertex objective values
        if (count > 1) call qsort(1,count)
    contains
        recursive subroutine qsort(left,right)
            integer, intent(in) :: left, right !! inclusive sort range
            integer :: i, j, pivot, temp
            if (left >= right) return
            i = left
            j = right
            pivot = indices((left+right)/2)
            do
                do while (vertex_less(indices(i),pivot))
                    i = i+1
                end do
                do while (vertex_less(pivot,indices(j)))
                    j = j-1
                end do
                if (i > j) exit
                temp = indices(i)
                indices(i) = indices(j)
                indices(j) = temp
                i = i+1
                j = j-1
                if (i > j) exit
            end do
            if (left < j) call qsort(left,j)
            if (i < right) call qsort(i,right)
        end subroutine qsort
        pure logical function vertex_less(a,b)
            integer, intent(in) :: a,b !! flattened vertex indices
            vertex_less = values(a) < values(b)
            if (values(a) == values(b)) vertex_less = a < b
        end function vertex_less
    end subroutine sort_candidates

    subroutine sort_minima(points, values, count)
        real(dp), intent(inout) :: points(:,:) !! stored local-minimum columns
        real(dp), intent(inout) :: values(:) !! objective values to sort
        integer, intent(in) :: count !! number of populated columns
        real(dp) :: temp_point(size(points,1)), temp_value
        integer :: i,j
        do i = 2, count
            temp_value = values(i)
            temp_point = points(:,i)
            j = i-1
            do while (j >= 1)
                if (values(j) <= temp_value) exit
                values(j+1) = values(j)
                points(:,j+1) = points(:,j)
                j = j-1
            end do
            values(j+1) = temp_value
            points(:,j+1) = temp_point
        end do
    end subroutine sort_minima

end module scifort_optimize_shgo
