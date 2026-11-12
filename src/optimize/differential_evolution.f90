! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Differential evolution, with algorithmic attribution to the Storn-Price
! DE/rand/1/bin method and its MIT-licensed Fortran presentation by
! Ian Morrall (2017), https://github.com/ian-mmm/differential-evolution_f95 .
! This module is a new implementation with SciFort interfaces, local RNG,
! robust donor sampling, bounds repair, convergence and accounting.
! For original MIT license and attribution see
! LICENSES/ian-morrall-differential-evolution-MIT.txt .
!
! Public interface analogous to scipy.optimize.differential_evolution.
! In bounds(n,2), column 1 is lower and column 2 is upper.
!
module scifort_optimize_differential_evolution
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: de_status_converged = 0
    integer, parameter, public :: de_status_maxiter = 1
    integer, parameter, public :: de_status_callback = 2
    integer, parameter, public :: de_status_maxfun = 3
    integer, parameter, public :: de_status_infeasible = 4
    integer, parameter, public :: de_status_invalid = -1
    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64
    real(dp), parameter :: bad_energy = huge(1.0_dp)/16.0_dp

    type, public :: differential_evolution_result
        real(dp), allocatable :: x(:) !! best feasible point, shape (n)
        real(dp) :: fun = huge(1.0_dp) !! objective at best point
        real(dp), allocatable :: population(:, :) !! final points, shape (n,np)
        real(dp), allocatable :: population_energies(:) !! final objective values, shape (np)
        integer :: nit = 0 !! completed generations
        integer :: nfev = 0 !! number of objective evaluations
        integer :: ncev = 0 !! constraint evaluations, includes infeasible candidates
        real(dp) :: maxcv = 0.0_dp !! largest constraint violation at returned x
        logical :: feasible = .true. !! whether returned point satisfies constraints
        real(dp), allocatable :: population_constraint_violation(:) !! maximum violation per member
        integer :: status = de_status_invalid !! termination code
        logical :: success = .false. !! converged
        character(len=120) :: message = 'invalid input' !! termination description
    end type differential_evolution_result

    abstract interface
        function de_objective(x) result(value)
            import dp
            real(dp), intent(in) :: x(:) !! candidate vector, size n
            real(dp) :: value !! objective to minimize
        end function de_objective

        subroutine de_callback(x, fun, generation, stop)
            import dp
            real(dp), intent(in) :: x(:) !! current best point, shape (n)
            real(dp), intent(in) :: fun !! current best objective value
            integer, intent(in) :: generation !! completed generations
            logical, intent(out) :: stop !! request early termination
        end subroutine de_callback

        subroutine de_constraint_function(x, c)
            import dp
            real(dp), intent(in) :: x(:) !! candidate vector, shape (n)
            real(dp), intent(out) :: c(:) !! nonlinear constraint values, shape (m)
        end subroutine de_constraint_function

        subroutine de_polisher(func, x0, lower, upper, x, fun, nfev, success)
            import dp, de_objective
            procedure(de_objective) :: func !! objective to minimize
            real(dp), intent(in) :: x0(:) !! best starting point, shape (n)
            real(dp), intent(in) :: lower(:) !! bounds (integer coordinates are fixed)
            real(dp), intent(in) :: upper(:) !! bounds (integer coordinates are fixed)
            real(dp), intent(out) :: x(:) !! final point, shape (n)
            real(dp), intent(out) :: fun !! objective at returned point
            integer, intent(out) :: nfev !! objective evaluations made by polisher
            logical, intent(out) :: success !! true if polisher returned valid point
        end subroutine de_polisher
    end interface

    public :: differential_evolution, de_objective, de_callback
    public :: de_constraint_function, de_polisher

contains

    function differential_evolution(func, bounds, strategy, maxiter, popsize, tol, atol, &
                                    mutation, mutation_range, recombination, seed, init, &
                                    updating, x0, maxfun, callback, integrality, &
                                    linear_a, linear_lb, linear_ub, nonlinear, nonlinear_lb, &
                                    nonlinear_ub, constraint_tol, polisher) result(res)
        procedure(de_objective) :: func !! real scalar objective function
        real(dp), intent(in) :: bounds(:, :) !! finite lower/upper bounds, shape (n,2)
        character(len=*), optional, intent(in) :: strategy !! best1bin/rand1bin/etc.
        integer, optional, intent(in) :: maxiter !! maximum generations; default 1000
        integer, optional, intent(in) :: popsize !! population multiplier; default 15
        real(dp), optional, intent(in) :: tol !! relative spread tolerance; default .01
        real(dp), optional, intent(in) :: atol !! absolute spread tolerance; default 0
        real(dp), optional, intent(in) :: mutation !! constant mutation scale in [0,2)
        real(dp), optional, intent(in) :: mutation_range(2) !! dither interval [fmin,fmax), within [0,2)
        real(dp), optional, intent(in) :: recombination !! crossover probability in [0,1]
        integer, optional, intent(in) :: seed !! local RNG seed; deterministic when supplied
        character(len=*), optional, intent(in) :: init !! 'latinhypercube' or 'random'
        character(len=*), optional, intent(in) :: updating !! 'immediate' or 'deferred'
        real(dp), optional, intent(in) :: x0(:) !! initial point, must be within bounds
        integer, optional, intent(in) :: maxfun !! maximum objective evaluations
        procedure(de_callback), optional :: callback !! called after each completed generation
        logical, optional, intent(in) :: integrality(:) !! true for integer coordinates, size n
        real(dp), optional, intent(in) :: linear_a(:, :) !! linear constraint matrix, (m,n)
        real(dp), optional, intent(in) :: linear_lb(:) !! lower bounds on A*x, size m
        real(dp), optional, intent(in) :: linear_ub(:) !! upper bounds on A*x, size m
        procedure(de_constraint_function), optional :: nonlinear !! nonlinear constraint callback
        real(dp), optional, intent(in) :: nonlinear_lb(:) !! lower nonlinear bounds, size m
        real(dp), optional, intent(in) :: nonlinear_ub(:) !! upper nonlinear bounds, size m
        real(dp), optional, intent(in) :: constraint_tol !! absolute feasibility tolerance, default 1e-10
        procedure(de_polisher), optional :: polisher !! optional local bounded minimizer
        type(differential_evolution_result) :: res

        real(dp), allocatable :: pop(:, :), next_pop(:, :), energy(:), next_energy(:)
        real(dp), allocatable :: donor(:), trial(:), lo(:), hi(:), width(:)
        real(dp), allocatable :: int_lo(:), int_hi(:), violation(:), next_violation(:)
        real(dp), allocatable :: max_violation(:), next_max_violation(:)
        real(dp), allocatable :: linear_values(:), nonlinear_values(:), polish_x(:)
        logical, allocatable :: integer_mask(:), feasible(:), next_feasible(:)
        integer, allocatable :: perm(:)
        integer :: n, nfree, np, niter_limit, multiplier, eval_limit, iteration, i, j
        integer :: best_idx, mode, crossover_type, donor_count, first, gene, other
        integer :: picked(5), p, chosen
        integer(int64) :: state, clock_seed
        real(dp) :: convergence_tol, convergence_atol, cross, fmin, fmax, factor
        real(dp) :: u, candidate_fun, candidate_total_violation, candidate_maxcv, cv_tol
        real(dp) :: polished_fun
        integer :: polished_evals, nlinear, nnonlinear
        character(len=32) :: alg, init_method, update_method
        logical :: use_deferred, stop_requested, candidate_feasible, take_trial
        logical :: polished_success, has_constraints

        n = size(bounds, 1)
        if (n <= 0 .or. size(bounds, 2) /= 2) then
            res%message = 'bounds must have shape (n,2), n >= 1'
            return
        end if
        if (.not. all(ieee_is_finite(bounds))) then
            res%message = 'all bounds must be finite'
            return
        end if
        if (any(bounds(:, 1) > bounds(:, 2))) then
            res%message = 'lower bounds must not exceed upper bounds'
            return
        end if
        cv_tol = 1.0e-10_dp
        if (present(constraint_tol)) cv_tol = constraint_tol
        if (.not. ieee_is_finite(cv_tol)) then
            res%message = 'constraint_tol must be finite'
            return
        end if
        if (cv_tol < 0.0_dp) then
            res%message = 'constraint_tol must be nonnegative'
            return
        end if
        allocate(integer_mask(n), int_lo(n), int_hi(n))
        integer_mask = .false.
        if (present(integrality)) then
            if (size(integrality) /= n) then
                res%message = 'integrality must have length n'
                return
            end if
            integer_mask = integrality
        end if
        int_lo = bounds(:,1)
        int_hi = bounds(:,2)
        do j = 1,n
            if (.not. integer_mask(j)) cycle
            if (max(abs(bounds(j,1)),abs(bounds(j,2))) > 9007199254740991.0_dp) then
                res%message = 'integer bounds exceed exactly representable range'
                return
            end if
            int_lo(j) = real(ceiling(bounds(j,1),kind=int64),dp)
            int_hi(j) = real(floor(bounds(j,2),kind=int64),dp)
            if (int_lo(j) > int_hi(j)) then
                res%message = 'integral parameter has no allowed integer values'
                return
            end if
        end do
        nlinear = 0
        nnonlinear = 0
        if ((present(linear_a) .neqv. present(linear_lb)) .or. &
            (present(linear_a) .neqv. present(linear_ub))) then
            res%message = 'linear_a, linear_lb and linear_ub must be supplied together'
            return
        end if
        if ((present(nonlinear) .neqv. present(nonlinear_lb)) .or. &
            (present(nonlinear) .neqv. present(nonlinear_ub))) then
            res%message = 'nonlinear and both nonlinear bound arrays must be supplied together'
            return
        end if
        if (present(linear_a)) then
            nlinear = size(linear_a,1)
            if (size(linear_a,2) /= n .or. size(linear_lb) /= nlinear .or. &
                size(linear_ub) /= nlinear) then
                res%message = 'linear constraint matrix or bound shape mismatch'
                return
            end if
            if (.not. all(ieee_is_finite(linear_a))) then
                res%message = 'linear constraint coefficients must be finite'
                return
            end if
            if (any(ieee_is_nan(linear_lb)) .or. any(ieee_is_nan(linear_ub)) .or. &
                any(linear_lb > linear_ub)) then
                res%message = 'invalid linear constraint bounds'
                return
            end if
        end if
        if (present(nonlinear)) then
            nnonlinear = size(nonlinear_lb)
            if (size(nonlinear_ub) /= nnonlinear) then
                res%message = 'nonlinear constraint bound length mismatch'
                return
            end if
            if (any(ieee_is_nan(nonlinear_lb)) .or. any(ieee_is_nan(nonlinear_ub)) .or. &
                any(nonlinear_lb > nonlinear_ub)) then
                res%message = 'invalid nonlinear constraint bounds'
                return
            end if
        end if
        has_constraints = nlinear > 0 .or. nnonlinear > 0
        if (present(linear_a)) allocate(linear_values(nlinear))
        if (present(nonlinear)) allocate(nonlinear_values(nnonlinear))
        nfree = count(bounds(:, 1) < bounds(:, 2))
        multiplier = 15
        if (present(popsize)) multiplier = popsize
        if (multiplier < 1 .or. multiplier > huge(1)/max(1,nfree)) then
            res%message = 'invalid or overflowing population multiplier'
            return
        end if
        np = max(6, multiplier*max(1,nfree))
        niter_limit = 1000
        if (present(maxiter)) niter_limit = maxiter
        if (niter_limit < 0) then
            res%message = 'maxiter must be nonnegative'
            return
        end if
        eval_limit = huge(1)
        if (present(maxfun)) eval_limit = maxfun
        if (eval_limit < np .and. nfree > 0) then
            res%message = 'maxfun must permit the initial population'
            return
        end if
        if (eval_limit < 1) then
            res%message = 'maxfun must be positive'
            return
        end if
        convergence_tol = 0.01_dp
        if (present(tol)) convergence_tol = tol
        convergence_atol = 0.0_dp
        if (present(atol)) convergence_atol = atol
        if (.not. ieee_is_finite(convergence_tol) .or. .not. ieee_is_finite(convergence_atol)) then
            res%message = 'tol and atol must be finite'
            return
        end if
        if (min(convergence_tol, convergence_atol) < 0.0_dp) then
            res%message = 'tol and atol must be nonnegative'
            return
        end if
        cross = 0.7_dp
        if (present(recombination)) cross = recombination
        if (.not. ieee_is_finite(cross)) then
            res%message = 'recombination must be finite'
            return
        end if
        if (cross < 0.0_dp .or. cross > 1.0_dp) then
            res%message = 'recombination must be between zero and one'
            return
        end if
        fmin = 0.5_dp
        fmax = 1.0_dp
        if (present(mutation)) then
            fmin = mutation
            fmax = mutation
        end if
        if (present(mutation_range)) then
            if (present(mutation)) then
                res%message = 'choose mutation or mutation_range, not both'
                return
            end if
            fmin = mutation_range(1)
            fmax = mutation_range(2)
        end if
        if (.not. ieee_is_finite(fmin) .or. .not. ieee_is_finite(fmax)) then
            res%message = 'mutation parameters must be finite'
            return
        end if
        if (fmin < 0.0_dp .or. fmax >= 2.0_dp .or. fmin > fmax) then
            res%message = 'mutation scale/range must satisfy 0 <= fmin <= fmax < 2'
            return
        end if

        alg = 'best1bin'
        if (present(strategy)) alg = trim(adjustl(strategy))
        crossover_type = 0
        if (index(trim(alg), 'bin') == len_trim(alg) - 2) crossover_type = 1
        if (index(trim(alg), 'exp') == len_trim(alg) - 2) crossover_type = 2
        if (crossover_type == 0) then
            res%message = 'unknown crossover strategy'
            return
        end if
        mode = 0
        select case (trim(alg(:len_trim(alg) - 3)))
        case ('best1')
            mode = 1
            donor_count = 2
        case ('rand1')
            mode = 2
            donor_count = 3
        case ('best2')
            mode = 3
            donor_count = 4
        case ('rand2')
            mode = 4
            donor_count = 5
        case ('randtobest1')
            mode = 5
            donor_count = 3
        case ('currenttobest1')
            mode = 6
            donor_count = 2
        case default
            res%message = 'unknown differential evolution strategy'
            return
        end select
        init_method = 'latinhypercube'
        if (present(init)) init_method = trim(adjustl(init))
        if (init_method /= 'latinhypercube' .and. init_method /= 'random') then
            res%message = 'init must be latinhypercube or random'
            return
        end if
        update_method = 'immediate'
        if (present(updating)) update_method = trim(adjustl(updating))
        if (update_method /= 'immediate' .and. update_method /= 'deferred') then
            res%message = 'updating must be immediate or deferred'
            return
        end if
        use_deferred = update_method == 'deferred'
        if (present(x0)) then
            if (size(x0) /= n) then
                res%message = 'x0 must have length n'
                return
            end if
            if (.not. all(ieee_is_finite(x0))) then
                res%message = 'x0 must be finite'
                return
            end if
            if (any(x0 < bounds(:,1)) .or. any(x0 > bounds(:,2))) then
                res%message = 'x0 must lie within bounds'
                return
            end if
        end if

        if (present(seed)) then
            state = modulo(int(seed,int64), rng_modulus-1_int64) + 1_int64
        else
            call system_clock(count=clock_seed)
            state = modulo(clock_seed, rng_modulus-1_int64) + 1_int64
        end if
        allocate(lo(n), hi(n), width(n), donor(n), trial(n))
        lo = bounds(:,1)
        hi = bounds(:,2)
        width = hi - lo
        if (.not. all(ieee_is_finite(width))) then
            res%message = 'bound widths must be finite'
            return
        end if
        if (nfree == 0) then
            allocate(res%x(n), res%population(n,1), res%population_energies(1), &
                     res%population_constraint_violation(1))
            res%x = lo
            call evaluate_violation(lo,candidate_total_violation,candidate_maxcv,nonlinear)
            res%maxcv = candidate_maxcv
            res%feasible = candidate_maxcv <= cv_tol
            res%fun = bad_energy
            if (res%feasible) then
                res%fun = sanitize_value(func(lo))
                res%nfev = 1
            end if
            res%population(:,1) = lo
            res%population_energies(1) = res%fun
            res%population_constraint_violation(1) = res%maxcv
            res%success = res%feasible .and. res%fun < bad_energy
            if (res%success) then
                res%status = de_status_converged
                res%message = 'all variables fixed by bounds'
            else if (.not. res%feasible) then
                res%status = de_status_infeasible
                res%message = 'fixed parameters violate constraints'
            else
                res%message = 'objective returned a nonfinite value'
            end if
            return
        end if
        allocate(pop(n,np), energy(np), next_pop(n,np), next_energy(np), perm(np))
        do j = 1,n
            if (hi(j) <= lo(j)) then
                pop(j,:) = lo(j)
            else if (init_method == 'random') then
                do i = 1,np
                    call uniform(state, u)
                    pop(j,i) = lo(j) + u*width(j)
                end do
            else
                perm = [(i, i=1,np)]
                do i = np,2,-1
                    call randint(state,i,p)
                    other = perm(i)
                    perm(i) = perm(p)
                    perm(p) = other
                end do
                do i = 1,np
                    call uniform(state, u)
                    pop(j,i) = lo(j) + (real(perm(i)-1,dp)+u)/real(np,dp)*width(j)
                end do
            end if
        end do
        if (present(x0)) pop(:,1) = x0
        do i = 1,np
            call enforce_integrality(pop(:,i))
        end do
        allocate(violation(np), next_violation(np), max_violation(np), &
                 next_max_violation(np), feasible(np), next_feasible(np))
        do i = 1,np
            call evaluate_violation(pop(:,i), violation(i), max_violation(i), nonlinear)
            feasible(i) = max_violation(i) <= cv_tol
            energy(i) = bad_energy
            if (feasible(i)) then
                energy(i) = sanitize_value(func(pop(:,i)))
                res%nfev = res%nfev + 1
            end if
        end do
        best_idx = choose_best(energy,feasible,violation,max_violation)
        res%status = de_status_maxiter
        res%message = 'maximum generations reached'
        if (res%nfev >= eval_limit .and. niter_limit > 0) then
            res%status = de_status_maxfun
            res%message = 'maximum function evaluations reached'
        end if

        do iteration = 1,niter_limit
            if (res%nfev >= eval_limit) exit
            call uniform(state,u)
            factor = fmin + (fmax-fmin)*u
            if (use_deferred) then
                next_pop = pop
                next_energy = energy
                next_violation = violation
                next_max_violation = max_violation
                next_feasible = feasible
            end if
            do i = 1,np
                if (res%nfev >= eval_limit) exit
                ! Distinct donor members, excluding the candidate.
                do j = 1,donor_count
                    do
                        call randint(state,np,chosen)
                        if (chosen == i) cycle
                        if (any(picked(1:j-1) == chosen)) cycle
                        exit
                    end do
                    picked(j) = chosen
                end do
                select case (mode)
                case (1)
                    donor = pop(:,best_idx) + factor*(pop(:,picked(1)) - pop(:,picked(2)))
                case (2)
                    donor = pop(:,picked(1)) + factor*(pop(:,picked(2)) - pop(:,picked(3)))
                case (3)
                    donor = pop(:,best_idx) + factor*(pop(:,picked(1)) + pop(:,picked(2)) - &
                        pop(:,picked(3)) - pop(:,picked(4)))
                case (4)
                    donor = pop(:,picked(1)) + factor*(pop(:,picked(2)) + pop(:,picked(3)) - &
                        pop(:,picked(4)) - pop(:,picked(5)))
                case (5)
                    donor = pop(:,picked(1)) + factor*(pop(:,best_idx) - pop(:,picked(1)) + &
                        pop(:,picked(2)) - pop(:,picked(3)))
                case (6)
                    donor = pop(:,i) + factor*(pop(:,best_idx) - pop(:,i) + &
                        pop(:,picked(1)) - pop(:,picked(2)))
                end select
                trial = pop(:,i)
                call randint(state,n,first)
                if (crossover_type == 1) then
                    do j = 1,n
                        call uniform(state,u)
                        if (u < cross .or. j == first) trial(j) = donor(j)
                    end do
                else
                    gene = first
                    do j = 1,n
                        trial(gene) = donor(gene)
                        if (j == n) exit
                        call uniform(state,u)
                        if (u >= cross) exit
                        gene = mod(gene,n)+1
                    end do
                end if
                do j = 1,n
                    if (.not. ieee_is_finite(trial(j))) then
                        call uniform(state,u)
                        trial(j) = lo(j) + u*width(j)
                    else if (trial(j) < lo(j) .or. trial(j) > hi(j)) then
                        call uniform(state,u)
                        trial(j) = lo(j) + u*width(j)
                    end if
                    if (hi(j) <= lo(j)) trial(j) = lo(j)
                end do
                call enforce_integrality(trial)
                call evaluate_violation(trial, candidate_total_violation, candidate_maxcv, nonlinear)
                candidate_feasible = candidate_maxcv <= cv_tol
                candidate_fun = bad_energy
                if (candidate_feasible) then
                    if (res%nfev >= eval_limit) exit
                    candidate_fun = sanitize_value(func(trial))
                    res%nfev = res%nfev + 1
                end if
                take_trial = .false.
                if (candidate_feasible .and. .not. feasible(i)) then
                    take_trial = .true.
                else if (candidate_feasible .and. feasible(i)) then
                    take_trial = candidate_fun <= energy(i)
                else if (.not. candidate_feasible .and. .not. feasible(i)) then
                    take_trial = candidate_total_violation < violation(i)
                    if (candidate_total_violation <= violation(i)) then
                        take_trial = candidate_maxcv < max_violation(i)
                    end if
                end if
                if (take_trial) then
                    if (use_deferred) then
                        next_pop(:,i) = trial
                        next_energy(i) = candidate_fun
                        next_violation(i) = candidate_total_violation
                        next_max_violation(i) = candidate_maxcv
                        next_feasible(i) = candidate_feasible
                    else
                        pop(:,i) = trial
                        energy(i) = candidate_fun
                        violation(i) = candidate_total_violation
                        max_violation(i) = candidate_maxcv
                        feasible(i) = candidate_feasible
                        best_idx = choose_best(energy,feasible,violation,max_violation)
                    end if
                end if
            end do
            if (use_deferred) then
                pop = next_pop
                energy = next_energy
                violation = next_violation
                max_violation = next_max_violation
                feasible = next_feasible
            end if
            best_idx = choose_best(energy,feasible,violation,max_violation)
            res%nit = iteration
            if (present(callback)) then
                call callback(pop(:,best_idx),energy(best_idx),iteration,stop_requested)
                if (stop_requested) then
                    res%status = de_status_callback
                    res%message = 'stopped by callback'
                    exit
                end if
            end if
            if (all(feasible) .and. converged(energy,convergence_tol,convergence_atol)) then
                res%status = de_status_converged
                res%success = .true.
                res%message = 'population energies converged'
                exit
            end if
            if (res%nfev >= eval_limit) then
                res%status = de_status_maxfun
                res%message = 'maximum function evaluations reached'
                exit
            end if
        end do
        allocate(res%x(n),res%population(n,np),res%population_energies(np), &
                 res%population_constraint_violation(np))
        res%x = pop(:,best_idx)
        res%fun = energy(best_idx)
        res%maxcv = max_violation(best_idx)
        res%feasible = feasible(best_idx)
        res%population = pop
        res%population_energies = energy
        res%population_constraint_violation = max_violation
        if (.not. res%feasible) then
            res%success = .false.
            res%status = de_status_infeasible
            res%message = 'no feasible solution found'
        end if
        if (present(polisher)) then
            if (res%feasible .and. res%fun < bad_energy .and. .not. has_constraints) then
                if (any(.not. integer_mask)) then
                    allocate(polish_x(n))
                    do j = 1,n
                        if (integer_mask(j)) then
                            int_lo(j) = res%x(j)
                            int_hi(j) = res%x(j)
                        end if
                    end do
                    do j = 1,n
                        if (.not. integer_mask(j)) then
                            int_lo(j) = lo(j)
                            int_hi(j) = hi(j)
                        end if
                    end do
                    call polisher(func,res%x,int_lo,int_hi,polish_x,polished_fun, &
                                  polished_evals,polished_success)
                    res%nfev = res%nfev + max(0,polished_evals)
                    if (polished_success .and. all(ieee_is_finite(polish_x)) .and. &
                        ieee_is_finite(polished_fun)) then
                        if (all(polish_x >= int_lo) .and. all(polish_x <= int_hi)) then
                            if (polished_fun < res%fun) then
                                res%x = polish_x
                                res%fun = polished_fun
                            end if
                        end if
                    end if
                end if
            end if
        end if
    contains
        subroutine enforce_integrality(x)
            real(dp), intent(inout) :: x(:) !! candidate, shape (n)
            integer :: jdim
            do jdim = 1,n
                if (.not. integer_mask(jdim)) cycle
                x(jdim) = max(int_lo(jdim),min(int_hi(jdim),anint(x(jdim))))
            end do
        end subroutine enforce_integrality

        subroutine evaluate_violation(x, total, maximum, nonlinear)
            real(dp), intent(in) :: x(:) !! candidate, shape (n)
            real(dp), intent(out) :: total !! sum of positive component violations
            real(dp), intent(out) :: maximum !! maximum positive component violation
            procedure(de_constraint_function), optional :: nonlinear !! nonlinear constraint callback
            integer :: jj
            total = 0.0_dp
            maximum = 0.0_dp
            if (.not. has_constraints) return
            res%ncev = res%ncev + 1
            if (present(linear_a)) then
                linear_values = matmul(linear_a,x)
                do jj = 1,nlinear
                    call accumulate_violation(linear_values(jj),linear_lb(jj),linear_ub(jj),total,maximum)
                end do
            end if
            if (present(nonlinear)) then
                call nonlinear(x,nonlinear_values)
                do jj = 1,nnonlinear
                    call accumulate_violation(nonlinear_values(jj),nonlinear_lb(jj),nonlinear_ub(jj),total,maximum)
                end do
            end if
        end subroutine evaluate_violation
    end function differential_evolution

    subroutine accumulate_violation(value,lower,upper,total,maximum)
        real(dp), intent(in) :: value !! constraint value
        real(dp), intent(in) :: lower !! lower permitted bound
        real(dp), intent(in) :: upper !! upper permitted bound
        real(dp), intent(inout) :: total !! total violation
        real(dp), intent(inout) :: maximum !! maximum violation
        real(dp) :: deviation
        if (.not. ieee_is_finite(value)) then
            deviation = bad_energy
        else
            deviation = max(0.0_dp,lower-value,value-upper)
            if (.not. ieee_is_finite(deviation)) deviation = bad_energy
        end if
        maximum = max(maximum,deviation)
        total = min(bad_energy,total+min(bad_energy,deviation))
    end subroutine accumulate_violation

    integer function choose_best(energy,feasible,violation,max_violation) result(idx)
        real(dp), intent(in) :: energy(:) !! population objective values
        logical, intent(in) :: feasible(:) !! constraint feasibility flags
        real(dp), intent(in) :: violation(:) !! total constraint violation
        real(dp), intent(in) :: max_violation(:) !! maximum component violation
        integer :: i
        idx = 1
        do i = 2,size(energy)
            if (feasible(i) .and. .not. feasible(idx)) then
                idx = i
            else if (feasible(i) .and. feasible(idx)) then
                if (energy(i) < energy(idx)) idx = i
            else if (.not. feasible(i) .and. .not. feasible(idx)) then
                if (violation(i) < violation(idx)) then
                    idx = i
                else if (violation(i) <= violation(idx)) then
                    if (max_violation(i) < max_violation(idx)) idx = i
                end if
            end if
        end do
    end function choose_best

    subroutine uniform(state,u)
        integer(int64), intent(inout) :: state !! local Park-Miller RNG state
        real(dp), intent(out) :: u !! uniform variate strictly between zero and one
        state = modulo(rng_multiplier*state,rng_modulus)
        u = real(state,dp)/real(rng_modulus,dp)
    end subroutine uniform

    subroutine randint(state,upper,j)
        integer(int64), intent(inout) :: state !! local RNG state
        integer, intent(in) :: upper !! upper integer bound, >=1
        integer, intent(out) :: j !! uniformly selected index in [1,upper]
        real(dp) :: u
        call uniform(state,u)
        j = min(upper,1+int(u*real(upper,dp)))
    end subroutine randint

    real(dp) function sanitize_value(value)
        real(dp), intent(in) :: value !! raw objective value
        sanitize_value = bad_energy
        if (.not. ieee_is_finite(value)) return
        if (value >= bad_energy) return
        sanitize_value = value
    end function sanitize_value

    logical function converged(values,tol,atol)
        real(dp), intent(in) :: values(:) !! objective values for population
        real(dp), intent(in) :: tol !! relative population-energy spread
        real(dp), intent(in) :: atol !! absolute population-energy spread
        real(dp) :: scale, mean_scaled, std_scaled
        integer :: n
        converged = .false.
        if (any(values >= bad_energy)) return
        n = size(values)
        if (n == 0) return
        scale = maxval(abs(values))
        if (scale <= 0.0_dp) then
            converged = .true.
            return
        end if
        mean_scaled = sum(values/scale)/real(n,dp)
        std_scaled = sqrt(sum((values/scale - mean_scaled)**2)/real(n,dp))
        converged = std_scaled <= atol/scale + tol*abs(mean_scaled)
    end function converged

end module scifort_optimize_differential_evolution
