! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Local-polishing adapter for deterministic simplicial global search.
! Reuses SciFort's bounded Nelder-Mead and Powell optimizers.
module scifort_optimize_shgo_minimize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, shgo_objective, &
        shgo_local_minimizer, shgo_callback, shgo_constraint_function, &
        shgo_constraint_jacobian, shgo_nonlinear_constraint, &
        shgo_constraint_blocks_valid
    use scifort_optimize_minimize, only : nelder_mead, powell, slsqp, multivariate_minimize_result
    implicit none
    private
    public :: shgo

contains

    function shgo(fun, bounds, n, iters, maxfun, max_mesh_points, method, &
            local_maxiter, local_iter, minimize_every_iter, f_min, f_tol, &
            merge_tol, minimizer, callback, adaptive, adaptive_fraction, &
            hgr_stability_iters, hgr_min_iters, hgr_f_tol, conforming, &
            linear_a, linear_lb, linear_ub, nonlinear, nonlinear_lb, &
            nonlinear_ub, constraint_tol, nonlinear_jac, &
            restore_feasibility, restoration_starts, restoration_maxiter, &
            constraints) result(res)
        procedure(shgo_objective) :: fun !! scalar objective on the search rectangle
        real(dp), intent(in) :: bounds(:,:) !! finite bounds, shape (nvar,2)
        integer, intent(in), optional :: n !! minimum first-mesh vertices, default 3**nvar
        integer, intent(in), optional :: iters !! dyadic refinement levels, default 1
        integer, intent(in), optional :: maxfun !! strict objective-call budget, default 20000
        integer, intent(in), optional :: max_mesh_points !! grid memory limit, default 100000
        character(len=*), intent(in), optional :: method !! 'nelder-mead' (default) or 'powell'
        integer, intent(in), optional :: local_maxiter !! local solver limit, default 200
        integer, intent(in), optional :: local_iter !! number of local starts per level, default all
        logical, intent(in), optional :: minimize_every_iter !! polish each level, default true
        real(dp), intent(in), optional :: f_min !! optional known global minimum
        real(dp), intent(in), optional :: f_tol !! relative target tolerance, default 1e-8
        real(dp), intent(in), optional :: merge_tol !! normalized duplicate threshold, default 1e-4
        procedure(shgo_local_minimizer), optional :: minimizer !! override built-in local search
        procedure(shgo_callback), optional :: callback !! iteration-end/early-stop hook
        logical, intent(in), optional :: adaptive !! selective dyadic leaf refinement
        real(dp), intent(in), optional :: adaptive_fraction !! fraction of leaves, default 0.25
        integer, intent(in), optional :: hgr_stability_iters !! consecutive stable comparisons, default 0
        integer, intent(in), optional :: hgr_min_iters !! minimum round before stopping, default 3
        real(dp), intent(in), optional :: hgr_f_tol !! relative improvement tolerance, default 1e-8
        logical, intent(in), optional :: conforming !! conforming 2D/3D adaptive mesh, default false
        real(dp), intent(in), optional :: linear_a(:,:) !! linear constraint matrix (m,nvar)
        real(dp), intent(in), optional :: linear_lb(:) !! lower limits of A*x
        real(dp), intent(in), optional :: linear_ub(:) !! upper limits of A*x
        procedure(shgo_constraint_function), optional :: nonlinear !! vector constraint callback
        real(dp), intent(in), optional :: nonlinear_lb(:) !! lower nonlinear limits
        real(dp), intent(in), optional :: nonlinear_ub(:) !! upper nonlinear limits
        real(dp), intent(in), optional :: constraint_tol !! absolute component feasibility tolerance
        procedure(shgo_constraint_jacobian), optional :: nonlinear_jac !! derivatives of nonlinear(x,c)
        logical, intent(in), optional :: restore_feasibility !! project mesh candidates onto constraints
        integer, intent(in), optional :: restoration_starts !! selected starts per round, default 3
        integer, intent(in), optional :: restoration_maxiter !! projection iteration limit, default 25
        type(shgo_nonlinear_constraint), intent(in), optional :: constraints(:) !! separate nonlinear blocks
        type(shgo_result) :: res !! global minimum, all local minima, and counters

        character(len=24) :: local_method
        integer :: local_limit, j, ascii_code, meq_local, mineq_local, nc_local, ncjev_local
        integer :: nnonlinear, block, m, offset
        real(dp), allocatable :: linear_buffer(:), nonlinear_buffer(:), nonlinear_jac_buffer(:,:)
        real(dp), allocatable :: all_nonlinear_lb(:), all_nonlinear_ub(:)
        logical :: have_constraints, use_constraint_jac

        ! Do not inspect unallocated block bounds or dereference disassociated
        ! callbacks before validating the records (including mixed legacy API).
        if (present(constraints)) then
            if (present(nonlinear) .or. present(nonlinear_lb) .or. &
                    present(nonlinear_ub) .or. present(nonlinear_jac)) then
                res%message = 'constraints cannot be combined with legacy nonlinear arguments'
                return
            end if
            if (.not. shgo_constraint_blocks_valid(constraints)) then
                res%message = 'invalid nonlinear constraint block or bounds'
                return
            end if
        end if
        nnonlinear = 0
        if (present(constraints)) then
            do block = 1, size(constraints)
                if (size(constraints(block)%lb) > huge(1)-nnonlinear) then
                    res%message = 'too many nonlinear constraint components'
                    return
                end if
                nnonlinear = nnonlinear+size(constraints(block)%lb)
            end do
        else if (present(nonlinear_lb)) then
            nnonlinear = size(nonlinear_lb)
        end if
        have_constraints = .false.
        if (present(linear_a)) have_constraints = have_constraints .or. size(linear_a,1) > 0
        have_constraints = have_constraints .or. nnonlinear > 0
        local_method = 'nelder-mead'
        if (have_constraints) local_method = 'slsqp'
        meq_local = 0
        mineq_local = 0
        nc_local = 0
        ncjev_local = 0
        if (have_constraints) then
            if (present(linear_a) .and. present(linear_lb) .and. present(linear_ub)) then
                allocate(linear_buffer(size(linear_a,1)))
                call count_components(linear_lb,linear_ub,meq_local,mineq_local)
            end if
            if (present(constraints)) then
                allocate(nonlinear_buffer(nnonlinear), nonlinear_jac_buffer(nnonlinear,size(bounds,1)))
                allocate(all_nonlinear_lb(nnonlinear), all_nonlinear_ub(nnonlinear))
                offset = 0
                do block = 1, size(constraints)
                    m = size(constraints(block)%lb)
                    all_nonlinear_lb(offset+1:offset+m) = constraints(block)%lb
                    all_nonlinear_ub(offset+1:offset+m) = constraints(block)%ub
                    offset = offset+m
                end do
                call count_components(all_nonlinear_lb,all_nonlinear_ub,meq_local,mineq_local)
            else if (present(nonlinear) .and. present(nonlinear_lb) .and. present(nonlinear_ub)) then
                allocate(nonlinear_buffer(nnonlinear))
                allocate(all_nonlinear_lb(nnonlinear),all_nonlinear_ub(nnonlinear))
                all_nonlinear_lb = nonlinear_lb
                all_nonlinear_ub = nonlinear_ub
                if (present(nonlinear_jac)) then
                    allocate(nonlinear_jac_buffer(nnonlinear,size(bounds,1)))
                end if
                call count_components(all_nonlinear_lb,all_nonlinear_ub,meq_local,mineq_local)
            end if
        end if
        use_constraint_jac = .not. allocated(nonlinear_buffer) .or. present(nonlinear_jac)
        if (present(constraints)) use_constraint_jac = .true.
        if (present(method)) local_method = trim(adjustl(method))
        do j = 1, len_trim(local_method)
            ascii_code = iachar(local_method(j:j))
            if (ascii_code >= iachar('A') .and. ascii_code <= iachar('Z')) then
                local_method(j:j) = achar(ascii_code + 32)
            end if
        end do
        local_limit = 200
        if (present(local_maxiter)) local_limit = local_maxiter
        if (local_limit <= 0) then
            res%message = 'local_maxiter must be positive'
            return
        end if
        if (.not. present(minimizer)) then
            select case (trim(local_method))
            case ('nelder-mead', 'nelder_mead', 'powell')
                if (have_constraints) then
                    res%message = 'constrained local polishing requires slsqp or custom minimizer'
                    return
                end if
            case ('slsqp')
                continue
            case default
                res%message = 'unsupported local method; use nelder-mead, powell, slsqp'
                return
            end select
        end if

        if (present(minimizer)) then
            res = shgo_core(fun, bounds, n=n, iters=iters, maxfun=maxfun, &
                max_mesh_points=max_mesh_points, local_iter=local_iter, &
                minimize_every_iter=minimize_every_iter, f_min=f_min, f_tol=f_tol, &
                merge_tol=merge_tol, minimizer=minimizer, callback=callback, &
                adaptive=adaptive, adaptive_fraction=adaptive_fraction, &
                hgr_stability_iters=hgr_stability_iters, hgr_min_iters=hgr_min_iters, &
                hgr_f_tol=hgr_f_tol, conforming=conforming, &
                linear_a=linear_a, linear_lb=linear_lb, linear_ub=linear_ub, &
                nonlinear=nonlinear, nonlinear_lb=nonlinear_lb, &
                nonlinear_ub=nonlinear_ub, constraint_tol=constraint_tol, &
                nonlinear_jac=nonlinear_jac, restore_feasibility=restore_feasibility, &
                restoration_starts=restoration_starts, restoration_maxiter=restoration_maxiter, &
                constraints=constraints)
        else
            res = shgo_core(fun, bounds, n=n, iters=iters, maxfun=maxfun, &
                max_mesh_points=max_mesh_points, local_iter=local_iter, &
                minimize_every_iter=minimize_every_iter, f_min=f_min, f_tol=f_tol, &
                merge_tol=merge_tol, minimizer=local_search, callback=callback, &
                adaptive=adaptive, adaptive_fraction=adaptive_fraction, &
                hgr_stability_iters=hgr_stability_iters, hgr_min_iters=hgr_min_iters, &
                hgr_f_tol=hgr_f_tol, conforming=conforming, &
                linear_a=linear_a, linear_lb=linear_lb, linear_ub=linear_ub, &
                nonlinear=nonlinear, nonlinear_lb=nonlinear_lb, &
                nonlinear_ub=nonlinear_ub, constraint_tol=constraint_tol, &
                nonlinear_jac=nonlinear_jac, restore_feasibility=restore_feasibility, &
                restoration_starts=restoration_starts, restoration_maxiter=restoration_maxiter, &
                constraints=constraints)
        end if
        res%ncev = res%ncev + nc_local
        res%ncjev = res%ncjev + ncjev_local

    contains

        subroutine count_components(lb,ub,ne,ni)
            real(dp), intent(in) :: lb(:) !! lower bounds for one constraint block
            real(dp), intent(in) :: ub(:) !! upper bounds for one constraint block
            integer, intent(inout) :: ne !! equality count
            integer, intent(inout) :: ni !! inequality count
            integer :: k
            if (size(lb) /= size(ub)) return
            do k = 1, size(lb)
                if (ieee_is_finite(lb(k)) .and. ieee_is_finite(ub(k))) then
                    if (lb(k) == ub(k)) then
                        ne = ne+1
                        cycle
                    end if
                end if
                if (ieee_is_finite(lb(k))) ni = ni+1
                if (ieee_is_finite(ub(k))) ni = ni+1
            end do
        end subroutine count_components

        subroutine evaluate_nonlinear(x, values, nonlinear)
            procedure(shgo_constraint_function), optional :: nonlinear !! legacy nonlinear callback
            real(dp), intent(in) :: x(:) !! candidate local point
            real(dp), intent(out) :: values(:) !! concatenated nonlinear constraint values
            integer :: k, offset, m
            if (present(constraints)) then
                offset = 0
                do k = 1, size(constraints)
                    m = size(constraints(k)%lb)
                    if (m > 0) call constraints(k)%fun(x,values(offset+1:offset+m))
                    offset = offset+m
                end do
            else
                call nonlinear(x,values)
            end if
        end subroutine evaluate_nonlinear

        subroutine evaluate_nonlinear_jac(x, a, nonlinear_jac)
            procedure(shgo_constraint_jacobian), optional :: nonlinear_jac !! legacy Jacobian callback
            ! Fill only missing blocks with finite differences. Supplying a
            ! Jacobian for one block does not force users to provide the others.
            real(dp), intent(in) :: x(:) !! bounded local point
            real(dp), intent(out) :: a(:,:) !! concatenated nonlinear Jacobian
            integer :: k, offset, m
            if (present(constraints)) then
                offset = 0
                do k = 1, size(constraints)
                    m = size(constraints(k)%lb)
                    if (m == 0) cycle
                    if (associated(constraints(k)%jac)) then
                        call constraints(k)%jac(x,a(offset+1:offset+m,:))
                        ncjev_local = ncjev_local+1
                    else
                        call finite_difference_block(x,a(offset+1:offset+m,:),constraints(k)%fun)
                    end if
                    offset = offset+m
                end do
            else
                call nonlinear_jac(x,a)
                ncjev_local = ncjev_local+1
            end if
        end subroutine evaluate_nonlinear_jac

        subroutine finite_difference_block(x, a, evaluate)
            ! Bound-aware one-sided differences. Each extra callback call is
            ! reflected in the constraint-evaluation diagnostics.
            real(dp), intent(in) :: x(:) !! local bounded coordinates
            real(dp), intent(out) :: a(:,:) !! block Jacobian
            procedure(shgo_constraint_function) :: evaluate !! block callback
            real(dp) :: base(size(a,1)), shifted(size(a,1)), perturbed(size(x))
            real(dp) :: step
            integer :: j
            a = 0.0_dp
            call evaluate(x,base)
            nc_local = nc_local+1
            if (.not. all(ieee_is_finite(base))) then
                a = huge(1.0_dp)
                return
            end if
            do j = 1, size(x)
                perturbed = x
                step = max(1.0e-6_dp*(bounds(j,2)-bounds(j,1)), &
                    sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(j))))
                perturbed(j) = min(bounds(j,2),x(j)+step)
                if (perturbed(j) == x(j)) perturbed(j) = max(bounds(j,1),x(j)-step)
                if (perturbed(j) == x(j)) then
                    a = huge(1.0_dp)
                    return
                end if
                call evaluate(perturbed,shifted)
                nc_local = nc_local+1
                if (.not. all(ieee_is_finite(shifted))) then
                    a = huge(1.0_dp)
                    return
                end if
                a(:,j) = (shifted-base)/(perturbed(j)-x(j))
            end do
        end subroutine finite_difference_block

        subroutine local_constraints(x,c)
            real(dp), intent(in) :: x(:) !! SLSQP current point
            real(dp), intent(out) :: c(:) !! equalities, then inequalities >=0
            integer :: idx, k, block
            nc_local = nc_local+1
            if (allocated(linear_buffer)) linear_buffer = matmul(linear_a,x)
            if (allocated(nonlinear_buffer)) call evaluate_nonlinear(x,nonlinear_buffer,nonlinear)
            idx = 0
            do block = 1, 2
                if (block == 1) then
                    if (.not. allocated(linear_buffer)) cycle
                    do k = 1, size(linear_buffer)
                        if (ieee_is_finite(linear_lb(k)) .and. &
                                ieee_is_finite(linear_ub(k))) then
                            if (linear_lb(k) == linear_ub(k)) then
                                idx = idx+1
                                c(idx) = linear_buffer(k)-linear_lb(k)
                            end if
                        end if
                    end do
                else
                    if (.not. allocated(nonlinear_buffer)) cycle
                    do k = 1, size(nonlinear_buffer)
                        if (ieee_is_finite(all_nonlinear_lb(k)) .and. &
                                ieee_is_finite(all_nonlinear_ub(k))) then
                            if (all_nonlinear_lb(k) == all_nonlinear_ub(k)) then
                                idx = idx+1
                                c(idx) = nonlinear_buffer(k)-all_nonlinear_lb(k)
                            end if
                        end if
                    end do
                end if
            end do
            if (allocated(linear_buffer)) then
                do k = 1, size(linear_buffer)
                    if (ieee_is_finite(linear_lb(k)) .and. ieee_is_finite(linear_ub(k))) then
                        if (linear_lb(k) == linear_ub(k)) cycle
                    end if
                    if (ieee_is_finite(linear_lb(k))) then
                        idx = idx+1
                        c(idx) = linear_buffer(k)-linear_lb(k)
                    end if
                    if (ieee_is_finite(linear_ub(k))) then
                        idx = idx+1
                        c(idx) = linear_ub(k)-linear_buffer(k)
                    end if
                end do
            end if
            if (allocated(nonlinear_buffer)) then
                do k = 1, size(nonlinear_buffer)
                    if (ieee_is_finite(all_nonlinear_lb(k)) .and. ieee_is_finite(all_nonlinear_ub(k))) then
                        if (all_nonlinear_lb(k) == all_nonlinear_ub(k)) cycle
                    end if
                    if (ieee_is_finite(all_nonlinear_lb(k))) then
                        idx = idx+1
                        c(idx) = nonlinear_buffer(k)-all_nonlinear_lb(k)
                    end if
                    if (ieee_is_finite(all_nonlinear_ub(k))) then
                        idx = idx+1
                        c(idx) = all_nonlinear_ub(k)-nonlinear_buffer(k)
                    end if
                end do
            end if
        end subroutine local_constraints

        subroutine local_constraints_jac(x,a)
            ! The SLSQP Jacobian uses the exact same row ordering and bound
            ! signs as local_constraints: equalities, then inequalities >= 0.
            real(dp), intent(in) :: x(:) !! SLSQP point
            real(dp), intent(out) :: a(:,:) !! canonical Jacobian (meq+mineq,nvar)
            integer :: idx, k
            a = 0.0_dp
            if (allocated(nonlinear_jac_buffer)) then
                call evaluate_nonlinear_jac(x,nonlinear_jac_buffer,nonlinear_jac)
            end if
            idx = 0
            if (allocated(linear_buffer)) then
                do k = 1, size(linear_buffer)
                    if (ieee_is_finite(linear_lb(k)) .and. ieee_is_finite(linear_ub(k))) then
                        if (linear_lb(k) == linear_ub(k)) then
                            idx = idx+1
                            a(idx,:) = linear_a(k,:)
                        end if
                    end if
                end do
            end if
            if (allocated(nonlinear_buffer)) then
                do k = 1, size(nonlinear_buffer)
                    if (ieee_is_finite(all_nonlinear_lb(k)) .and. &
                            ieee_is_finite(all_nonlinear_ub(k))) then
                        if (all_nonlinear_lb(k) == all_nonlinear_ub(k)) then
                            idx = idx+1
                            a(idx,:) = nonlinear_jac_buffer(k,:)
                        end if
                    end if
                end do
            end if
            if (allocated(linear_buffer)) then
                do k = 1, size(linear_buffer)
                    if (ieee_is_finite(linear_lb(k)) .and. ieee_is_finite(linear_ub(k))) then
                        if (linear_lb(k) == linear_ub(k)) cycle
                    end if
                    if (ieee_is_finite(linear_lb(k))) then
                        idx = idx+1
                        a(idx,:) = linear_a(k,:)
                    end if
                    if (ieee_is_finite(linear_ub(k))) then
                        idx = idx+1
                        a(idx,:) = -linear_a(k,:)
                    end if
                end do
            end if
            if (allocated(nonlinear_buffer)) then
                do k = 1, size(nonlinear_buffer)
                    if (ieee_is_finite(all_nonlinear_lb(k)) .and. ieee_is_finite(all_nonlinear_ub(k))) then
                        if (all_nonlinear_lb(k) == all_nonlinear_ub(k)) cycle
                    end if
                    if (ieee_is_finite(all_nonlinear_lb(k))) then
                        idx = idx+1
                        a(idx,:) = nonlinear_jac_buffer(k,:)
                    end if
                    if (ieee_is_finite(all_nonlinear_ub(k))) then
                        idx = idx+1
                        a(idx,:) = -nonlinear_jac_buffer(k,:)
                    end if
                end do
            end if
        end subroutine local_constraints_jac

        subroutine local_search(objective, lower, upper, start, maxfev, x, value, nfev, success)
            procedure(shgo_objective) :: objective !! budget-guarded scalar objective
            real(dp), intent(in) :: lower(:) !! box lower bounds
            real(dp), intent(in) :: upper(:) !! box upper bounds
            real(dp), intent(in) :: start(:) !! feasible initial vertex
            integer, intent(in) :: maxfev !! strict local evaluation budget
            real(dp), intent(out) :: x(:) !! refined feasible point
            real(dp), intent(out) :: value !! finite refined objective
            integer, intent(out) :: nfev !! objective evaluations by local method
            logical, intent(out) :: success !! returned a finite feasible local solution
            type(multivariate_minimize_result) :: local

            x = start
            value = huge(1.0_dp)
            nfev = 0
            success = .false.
            if (maxfev < 1) return
            select case (trim(local_method))
            case ('nelder-mead', 'nelder_mead')
                local = nelder_mead(objective, start, lower=lower, upper=upper, &
                    maxiter=local_limit, maxfev=maxfev)
            case ('powell')
                local = powell(objective, start, lower=lower, upper=upper, &
                    maxiter=local_limit, maxfev=maxfev)
            case ('slsqp')
                if (meq_local+mineq_local > 0) then
                    if (use_constraint_jac) then
                        local = slsqp(objective, start, lower=lower, upper=upper, &
                            constraints=local_constraints, constraint_jac=local_constraints_jac, &
                            meq=meq_local, mineq=mineq_local, &
                            maxiter=local_limit, ftol=1.0e-10_dp)
                    else
                        local = slsqp(objective, start, lower=lower, upper=upper, &
                            constraints=local_constraints, meq=meq_local, mineq=mineq_local, &
                            maxiter=local_limit, ftol=1.0e-10_dp)
                    end if
                else
                    local = slsqp(objective, start, lower=lower, upper=upper, &
                        maxiter=local_limit, ftol=1.0e-10_dp)
                end if
            end select
            nfev = local%nfev
            if (.not. allocated(local%x)) return
            if (size(local%x) /= size(start)) return
            if (.not. ieee_is_finite(local%fun)) return
            if (.not. all(ieee_is_finite(local%x))) return
            if (local%fun >= huge(1.0_dp)) return
            if (any(local%x < lower) .or. any(local%x > upper)) return
            x = local%x
            value = local%fun
            ! A finite improvement is usable even if maxiter stopped the local solver.
            success = .true.
        end subroutine local_search

    end function shgo

end module scifort_optimize_shgo_minimize
