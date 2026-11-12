! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Multivariate minimization for the scipy.optimize translation.
!
! The Nelder-Mead, modified Powell, BFGS, and Polak-Ribiere+ conjugate-gradient
! algorithms, defaults, stopping rules, and accounting are adapted from SciPy
! 1.17.0 scipy/optimize/_optimize.py and scipy/optimize/_minimize.py; the
! dogleg/trust-ncg/trust-exact/trust-krylov controllers and their trust-region subproblems are
! adapted from SciPy 1.17.0 scipy/optimize/_trustregion.py,
! scipy/optimize/_trustregion_dogleg.py, scipy/optimize/_trustregion_ncg.py,
! scipy/optimize/_trustregion_exact.py, and scipy/optimize/_trustregion_krylov.py
! under SciPy's BSD-3-Clause license. The native trust-krylov reduced solver
! follows the GLTR/TRLIB algorithm and forcing rules used by SciPy. Powell's one-dimensional searches reuse
! SciFort's attributed scalar Brent implementation. COBYLA wrapper semantics,
! constraint conversion, bounds-only start projection, and stopping controls follow
! SciPy 1.17.0; its PRIMA-derived numerical kernel is in cobyla_prima.f90.
! The trust-constr path follows SciPy 1.17.0's Byrd-Omojokun equality SQP,
! projected-CG QP subproblem, and trust-region interior-point structure.
! Inequalities and finite bounds are converted
! to positive-slack barrier equalities solved by the dense equality SQP kernel.
! L-BFGS-B and SLSQP are original SciFort adapters around pinned external FPM
! numerical kernels; see CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.

module scifort_optimize_minimize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_linalg, only : cholesky_lower, linalg_status_success, solve_lower, &
        symmetric_eigen_jacobi
    use lbfgsb_module, only : setulb
    use slsqp_kinds, only : slsqp_wp => wp
    use slsqp_module, only : slsqp_solver
    use scifort_optimize_cobyla_prima, only : cobyla_damaging_rounding, &
        cobyla_ftarget_achieved, cobyla_maxfun_reached, cobyla_maxtr_reached, &
        cobyla_small_tr_radius, prima_cobyla
    use scifort_optimize_minimize_scalar, only : bounded, brent, minimize_result
    use scifort_optimize_wolfe, only : scalar_search_wolfe1, scalar_search_wolfe2
    implicit none
    private

    real(dp), parameter :: default_nm_xatol = 1.0e-4_dp
    real(dp), parameter :: default_nm_fatol = 1.0e-4_dp
    real(dp), parameter :: default_powell_xtol = 1.0e-4_dp
    real(dp), parameter :: default_powell_ftol = 1.0e-4_dp
    real(dp), parameter :: default_gtol = 1.0e-5_dp
    real(dp), parameter :: default_fd_eps = sqrt(epsilon(1.0_dp))

    type, public :: multivariate_minimize_result
        real(dp), allocatable :: x(:) !! final minimizer estimate
        real(dp) :: fun = 0.0_dp !! objective value at x
        real(dp), allocatable :: jac(:) !! final gradient when available
        real(dp), allocatable :: constr(:) !! final trust-constr values, user block then bounds block
        real(dp), allocatable :: constr_jac(:, :) !! matching row-stacked trust-constr constraint Jacobians
        real(dp), allocatable :: lagrangian_grad(:) !! final trust-constr Lagrangian gradient when available
        real(dp), allocatable :: v(:) !! final trust-constr working Lagrange multipliers when available
        real(dp), allocatable :: hess_inv(:, :) !! final BFGS inverse-Hessian estimate
        real(dp), allocatable :: final_simplex(:, :) !! Nelder-Mead vertices, one per row
        real(dp), allocatable :: final_simplex_fun(:) !! objective values at simplex vertices
        real(dp), allocatable :: direc(:, :) !! final Powell direction set, one per row
        integer :: nit = 0 !! completed outer iterations
        integer :: nfev = 0 !! objective evaluations
        integer :: njev = 0 !! gradient evaluations
        integer :: nhev = 0 !! Hessian or Hessian-vector evaluations
        integer :: cg_niter = 0 !! cumulative trust-constr projected-CG iterations
        integer :: cg_stop_cond = 0 !! final trust-constr projected-CG stop condition
        integer, allocatable :: constr_nfev(:) !! calls by trust-constr block: user constraints, then bounds
        integer, allocatable :: constr_njev(:) !! Jacobian calls in the same trust-constr block order
        integer, allocatable :: constr_nhev(:) !! Hessian calls in the same trust-constr block order
        integer :: status = -1 !! SciPy-style method-specific termination status
        logical :: success = .false. !! true when the method converged
        character(len=96) :: message = 'invalid input' !! termination description
        real(dp) :: maxcv = 0.0_dp !! maximum constraint violation when available
        real(dp) :: constr_violation = 0.0_dp !! SciPy-style trust-constr constraint violation
        real(dp) :: optimality = 0.0_dp !! final trust-constr Lagrangian-gradient infinity norm
        real(dp) :: constr_penalty = 0.0_dp !! final trust-constr merit penalty when available
        real(dp) :: tr_radius = 0.0_dp !! final trust-constr trust-region radius when available
        real(dp) :: barrier_parameter = 0.0_dp !! final trust-constr barrier parameter when applicable
        real(dp) :: barrier_tolerance = 0.0_dp !! final trust-constr barrier tolerance when applicable
        character(len=32) :: method = '' !! method actually used
    end type multivariate_minimize_result

    abstract interface
        function multivariate_objective_function(x) result(fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp) :: fx !! scalar objective value
        end function multivariate_objective_function

        subroutine multivariate_gradient_function(x, grad)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(out) :: grad(:) !! objective gradient at x
        end subroutine multivariate_gradient_function

        subroutine multivariate_hessian_function(x, hess)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(out) :: hess(:, :) !! objective Hessian at x
        end subroutine multivariate_hessian_function

        subroutine multivariate_hessian_product_function(x, p, hp)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(in) :: p(:) !! vector multiplied by the Hessian
            real(dp), intent(out) :: hp(:) !! Hessian-vector product
        end subroutine multivariate_hessian_product_function

        subroutine multivariate_constraint_function(x, c)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(out) :: c(:) !! equality constraints first, then inequalities >= 0
        end subroutine multivariate_constraint_function

        subroutine multivariate_constraint_jacobian(x, a)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(out) :: a(:, :) !! constraint Jacobian, one constraint per row
        end subroutine multivariate_constraint_jacobian

        subroutine multivariate_constraint_hessian(x, v, hess)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(in) :: v(:) !! Lagrange multipliers for the constraints
            real(dp), intent(out) :: hess(:, :) !! Hessian of dot_product(v, constraints(x))
        end subroutine multivariate_constraint_hessian

        subroutine multivariate_constraint_hessian_product(x, v, p, hp)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(in) :: v(:) !! multipliers for this constraint block
            real(dp), intent(in) :: p(:) !! vector multiplied by the weighted constraint Hessian
            real(dp), intent(out) :: hp(:) !! Hessian-vector product
        end subroutine multivariate_constraint_hessian_product

        subroutine trust_constr_iteration_callback(x, nit, stop)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! physical parameter vector at the current iteration
            integer, intent(in) :: nit !! iteration number, including the initial state as 1
            logical, intent(out) :: stop !! set true to terminate the optimization
        end subroutine trust_constr_iteration_callback
    end interface

    ! Independent SciPy-style constraint blocks. Allocating %matrix specifies
    ! a linear constraint matrix*x, otherwise %fun supplies nonlinear values.
    ! A bound array defines each block's output count; a missing endpoint is unbounded.
    type, public :: nonlinear_constraint
        procedure(multivariate_constraint_function), pointer, nopass :: fun => null() !! vector-valued constraint callback
        procedure(multivariate_constraint_jacobian), pointer, nopass :: jac => null() !! optional analytic Jacobian
        procedure(multivariate_constraint_hessian), pointer, nopass :: hess => null() !! optional weighted Hessian
        procedure(multivariate_constraint_hessian_product), pointer, nopass :: hessp => null() !! weighted Hessian product
        real(dp), allocatable :: matrix(:, :) !! linear constraint matrix, one output per row
        real(dp), allocatable :: lower(:) !! component lower bounds; absent means -infinity
        real(dp), allocatable :: upper(:) !! component upper bounds; absent means +infinity
        logical, allocatable :: keep_feasible(:) !! retain feasibility of inequality components
        character(len=8) :: fd_scheme = '2-point' !! numerical Jacobian: '2-point' or '3-point'
        character(len=8) :: hess_fd_scheme = '' !! explicit numerical Hessian: '2-point', '3-point', or none
        character(len=8) :: hess_update = '' !! optional per-block Hessian update: 'bfgs' or 'sr1'
        real(dp), allocatable :: finite_diff_rel_step(:) !! signed relative steps; zero uses automatic step
        logical, allocatable :: finite_diff_jac_sparsity(:, :) !! output-by-variable FD Jacobian pattern
    end type nonlinear_constraint

    public :: multivariate_objective_function
    public :: multivariate_gradient_function
    public :: multivariate_hessian_function
    public :: multivariate_hessian_product_function
    public :: multivariate_constraint_function
    public :: multivariate_constraint_jacobian
    public :: multivariate_constraint_hessian
    public :: multivariate_constraint_hessian_product
    public :: trust_constr_iteration_callback
    public :: linear_constraint
    public :: minimize
    public :: nelder_mead
    public :: powell
    public :: bfgs
    public :: cg
    public :: newton_cg
    public :: dogleg
    public :: trust_ncg
    public :: trust_exact
    public :: trust_krylov
    public :: trust_constr
    public :: lbfgsb
    public :: slsqp
    public :: tnc
    public :: cobyla

contains

    function linear_constraint(matrix, lower, upper, keep_feasible) result(block)
        ! Build a constant-matrix block compatible with constraint_objects=.
        ! Shape and bounds validation occurs in trust_constr_objects.
        real(dp), intent(in) :: matrix(:, :) !! linear constraint matrix, one row per output
        real(dp), intent(in), optional :: lower(:) !! lower bounds for matrix*x
        real(dp), intent(in), optional :: upper(:) !! upper bounds for matrix*x
        logical, intent(in), optional :: keep_feasible(:) !! retained inequality components
        type(nonlinear_constraint) :: block !! resulting linear constraint block

        allocate(block%matrix(size(matrix, 1), size(matrix, 2)))
        block%matrix = matrix
        if (present(lower)) block%lower = lower
        if (present(upper)) block%upper = upper
        if (present(keep_feasible)) block%keep_feasible = keep_feasible
    end function linear_constraint

    function minimize(fun, x0, method, jac, hess, hessp, tol, maxiter, maxfev, &
            xatol, fatol, xtol, ftol, gtol, eps, adaptive, initial_simplex, direc, &
            xrtol, c1, c2, hess_inv0, lower, upper, maxfun, maxcor, constraints, &
            constraint_jac, meq, mineq, initial_trust_radius, max_trust_radius, eta, &
            subproblem_maxiter, inexact, max_cg_iter, rhobeg, catol, f_target, &
            barrier_tol, initial_barrier_parameter, initial_barrier_tolerance, &
            initial_constr_penalty, initial_tr_radius, constraint_hess, &
            constraint_lower, constraint_upper, constraint_keep_feasible, &
            bounds_keep_feasible, constraint_objects, objective_fd_scheme, &
            objective_finite_diff_rel_step, objective_hess_update, callback, &
            objective_hess_fd_scheme) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        character(len=*), intent(in), optional :: method !! selected SciPy-style minimization method
        procedure(multivariate_gradient_function), optional :: jac !! analytic gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic Hessian
        procedure(multivariate_hessian_product_function), optional :: hessp !! Hessian product
        real(dp), intent(in), optional :: tol !! generic method-specific tolerance
        integer, intent(in), optional :: maxiter !! maximum outer iterations
        integer, intent(in), optional :: maxfev !! maximum objective evaluations
        real(dp), intent(in), optional :: xatol !! Nelder-Mead absolute parameter tolerance
        real(dp), intent(in), optional :: fatol !! Nelder-Mead absolute objective tolerance
        real(dp), intent(in), optional :: xtol !! Powell/Newton-CG/TNC/trust-constr step tolerance
        real(dp), intent(in), optional :: ftol !! Powell/TNC relative objective tolerance
        real(dp), intent(in), optional :: gtol !! gradient tolerance for supported methods
        real(dp), intent(in), optional :: eps !! forward-difference absolute step
        logical, intent(in), optional :: adaptive !! use dimension-adaptive Nelder-Mead constants
        real(dp), intent(in), optional :: initial_simplex(:, :) !! SciPy-shape (n+1,n) simplex
        real(dp), intent(in), optional :: direc(:, :) !! Powell direction rows, shape (n,n)
        real(dp), intent(in), optional :: xrtol !! BFGS relative step tolerance
        real(dp), intent(in), optional :: c1 !! strong-Wolfe sufficient-decrease constant
        real(dp), intent(in), optional :: c2 !! strong-Wolfe curvature constant
        real(dp), intent(in), optional :: hess_inv0(:, :) !! initial BFGS inverse Hessian
        real(dp), intent(in), optional :: lower(:) !! lower variable bounds; -Inf/NaN means unbounded
        real(dp), intent(in), optional :: upper(:) !! upper variable bounds; +Inf/NaN means unbounded
        integer, intent(in), optional :: maxfun !! L-BFGS-B/TNC maximum objective evaluations
        integer, intent(in), optional :: maxcor !! L-BFGS-B correction-history length
        procedure(multivariate_constraint_function), optional :: constraints !! combined constraints
        procedure(multivariate_constraint_jacobian), optional :: constraint_jac !! analytic constraint Jacobian
        integer, intent(in), optional :: meq !! number of leading equality constraints
        integer, intent(in), optional :: mineq !! number of trailing inequality constraints
        real(dp), intent(in), optional :: constraint_lower(:) !! trust-constr lower bounds on user constraint values
        real(dp), intent(in), optional :: constraint_upper(:) !! trust-constr upper bounds on user constraint values
        logical, intent(in), optional :: constraint_keep_feasible(:) !! keep bounded user constraint components feasible
        logical, intent(in), optional :: bounds_keep_feasible(:) !! keep variable-bound components feasible
        type(nonlinear_constraint), intent(in), optional :: constraint_objects(:) !! independent constraint blocks
        character(len=*), intent(in), optional :: objective_fd_scheme !! trust-constr objective '2-point' or '3-point'
        real(dp), intent(in), optional :: objective_finite_diff_rel_step(:) !! per-variable objective relative steps
        character(len=*), intent(in), optional :: objective_hess_update !! objective Hessian BFGS or SR1 update
        character(len=*), intent(in), optional :: objective_hess_fd_scheme !! numerical Hessian from analytic jac
        procedure(trust_constr_iteration_callback), optional :: callback !! trust-constr iteration/termination hook
        real(dp), intent(in), optional :: initial_trust_radius !! initial trust-region radius
        real(dp), intent(in), optional :: max_trust_radius !! maximum trust-region radius
        real(dp), intent(in), optional :: eta !! trust-region step acceptance threshold
        integer, intent(in), optional :: subproblem_maxiter !! trust-exact inner iteration limit
        logical, intent(in), optional :: inexact !! trust-krylov inexact GLTR forcing sequence
        integer, intent(in), optional :: max_cg_iter !! TNC inner CG cap; negative uses default
        real(dp), intent(in), optional :: rhobeg !! COBYLA initial trust-region radius
        real(dp), intent(in), optional :: catol !! COBYLA absolute constraint-violation tolerance
        real(dp), intent(in), optional :: f_target !! COBYLA objective target for early termination
        real(dp), intent(in), optional :: initial_constr_penalty !! trust-constr initial merit penalty
        real(dp), intent(in), optional :: initial_tr_radius !! SciPy trust-constr initial trust radius alias
        real(dp), intent(in), optional :: barrier_tol !! trust-constr barrier termination threshold
        real(dp), intent(in), optional :: initial_barrier_parameter !! trust-constr initial barrier parameter
        real(dp), intent(in), optional :: initial_barrier_tolerance !! trust-constr initial barrier subproblem tolerance
        procedure(multivariate_constraint_hessian), optional :: constraint_hess !! analytic constraint Hessian
        type(multivariate_minimize_result) :: res

        character(len=16) :: meth
        real(dp) :: barrier_tol_use
        real(dp) :: fatol_use
        real(dp) :: ftol_use
        real(dp) :: gtol_use
        real(dp) :: xatol_use
        real(dp) :: xtol_use

        meth = 'bfgs'
        if (present(constraint_objects)) then
            meth = 'trust-constr'
        else if (present(constraints) .or. present(meq) .or. present(mineq)) then
            meth = 'slsqp'
        else if (present(lower) .or. present(upper)) then
            meth = 'l-bfgs-b'
        end if
        if (present(method)) meth = lowercase(trim(method))

        if (present(constraint_objects) .and. trim(meth) /= 'trust-constr' .and. &
                trim(meth) /= 'trust_constr') then
            call init_result(res, trim(meth), size(x0))
            call finish_result(res, -1, .false., 'constraint_objects require trust-constr')
            return
        end if
        if (present(objective_fd_scheme) .and. trim(meth) /= 'trust-constr' .and. &
                trim(meth) /= 'trust_constr') then
            call init_result(res, trim(meth), size(x0))
            call finish_result(res, -1, .false., &
                'objective_fd_scheme is available only for trust-constr')
            return
        end if
        if (present(objective_finite_diff_rel_step) .and. trim(meth) /= 'trust-constr' .and. &
                trim(meth) /= 'trust_constr') then
            call init_result(res, trim(meth), size(x0))
            call finish_result(res, -1, .false., &
                'objective_finite_diff_rel_step is available only for trust-constr')
            return
        end if
        if (present(callback) .and. trim(meth) /= 'trust-constr' .and. &
                trim(meth) /= 'trust_constr') then
            call finish_result(res, -1, .false., 'callback is available only for trust-constr')
            return
        end if
        if (present(objective_hess_fd_scheme) .and. trim(meth) /= 'trust-constr' .and. &
                trim(meth) /= 'trust_constr') then
            call init_result(res, trim(meth), size(x0))
            call finish_result(res, -1, .false., &
                'objective_hess_fd_scheme is available only for trust-constr')
            return
        end if
        if (present(objective_hess_update) .and. trim(meth) /= 'trust-constr' .and. &
                trim(meth) /= 'trust_constr') then
            call init_result(res, trim(meth), size(x0))
            call finish_result(res, -1, .false., &
                'objective_hess_update is available only for trust-constr')
            return
        end if
        select case (trim(meth))
        case ('nelder-mead', 'nelder_mead', 'simplex')
            xatol_use = default_nm_xatol
            fatol_use = default_nm_fatol
            if (present(tol)) then
                xatol_use = tol
                fatol_use = tol
            end if
            if (present(xatol)) xatol_use = xatol
            if (present(fatol)) fatol_use = fatol
            res = nelder_mead(fun, x0, xatol=xatol_use, fatol=fatol_use, &
                maxiter=maxiter, maxfev=maxfev, adaptive=adaptive, &
                initial_simplex=initial_simplex, lower=lower, upper=upper)
        case ('powell')
            xtol_use = default_powell_xtol
            ftol_use = default_powell_ftol
            if (present(tol)) then
                xtol_use = tol
                ftol_use = tol
            end if
            if (present(xtol)) xtol_use = xtol
            if (present(ftol)) ftol_use = ftol
            res = powell(fun, x0, xtol=xtol_use, ftol=ftol_use, maxiter=maxiter, &
                maxfev=maxfev, direc=direc, lower=lower, upper=upper)
        case ('bfgs')
            gtol_use = default_gtol
            if (present(tol)) gtol_use = tol
            if (present(gtol)) gtol_use = gtol
            if (present(jac)) then
                res = bfgs(fun, x0, jac=jac, gtol=gtol_use, eps=eps, maxiter=maxiter, &
                    xrtol=xrtol, c1=c1, c2=c2, hess_inv0=hess_inv0)
            else
                res = bfgs(fun, x0, gtol=gtol_use, eps=eps, maxiter=maxiter, &
                    xrtol=xrtol, c1=c1, c2=c2, hess_inv0=hess_inv0)
            end if
        case ('cg')
            gtol_use = default_gtol
            if (present(tol)) gtol_use = tol
            if (present(gtol)) gtol_use = gtol
            if (present(jac)) then
                res = cg(fun, x0, jac=jac, gtol=gtol_use, eps=eps, maxiter=maxiter, &
                    c1=c1, c2=c2)
            else
                res = cg(fun, x0, gtol=gtol_use, eps=eps, maxiter=maxiter, &
                    c1=c1, c2=c2)
            end if
        case ('newton-cg', 'newton_cg')
            if (.not. present(jac)) then
                call init_result(res, 'newton-cg', size(x0))
                return
            end if
            if (present(hess)) then
                res = newton_cg(fun, x0, jac, hess=hess, xtol=xtol, eps=eps, &
                    maxiter=maxiter, c1=c1, c2=c2)
            else if (present(hessp)) then
                res = newton_cg(fun, x0, jac, hessp=hessp, xtol=xtol, eps=eps, &
                    maxiter=maxiter, c1=c1, c2=c2)
            else
                res = newton_cg(fun, x0, jac, xtol=xtol, eps=eps, maxiter=maxiter, &
                    c1=c1, c2=c2)
            end if
        case ('dogleg')
            gtol_use = 1.0e-4_dp
            if (present(tol)) gtol_use = tol
            if (present(gtol)) gtol_use = gtol
            if (.not. present(jac) .or. .not. present(hess)) then
                call init_result(res, 'dogleg', size(x0))
                return
            end if
            res = dogleg(fun, x0, jac, hess, &
                initial_trust_radius=initial_trust_radius, &
                max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                maxiter=maxiter)
        case ('trust-ncg', 'trust_ncg')
            gtol_use = 1.0e-4_dp
            if (present(tol)) gtol_use = tol
            if (present(gtol)) gtol_use = gtol
            if (.not. present(jac) .or. (.not. present(hess) .and. .not. present(hessp))) then
                call init_result(res, 'trust-ncg', size(x0))
                return
            end if
            if (present(hessp)) then
                if (present(hess)) then
                    res = trust_ncg(fun, x0, jac, hess=hess, hessp=hessp, &
                        initial_trust_radius=initial_trust_radius, &
                        max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                        maxiter=maxiter)
                else
                    res = trust_ncg(fun, x0, jac, hessp=hessp, &
                        initial_trust_radius=initial_trust_radius, &
                        max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                        maxiter=maxiter)
                end if
            else
                res = trust_ncg(fun, x0, jac, hess=hess, &
                    initial_trust_radius=initial_trust_radius, &
                    max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                    maxiter=maxiter)
            end if
        case ('trust-krylov', 'trust_krylov')
            gtol_use = 1.0e-4_dp
            if (present(tol)) gtol_use = tol
            if (present(gtol)) gtol_use = gtol
            if (.not. present(jac) .or. (.not. present(hess) .and. .not. present(hessp))) then
                call init_result(res, 'trust-krylov', size(x0))
                return
            end if
            if (present(hessp)) then
                if (present(hess)) then
                    res = trust_krylov(fun, x0, jac, hess=hess, hessp=hessp, &
                        initial_trust_radius=initial_trust_radius, &
                        max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                        maxiter=maxiter, inexact=inexact)
                else
                    res = trust_krylov(fun, x0, jac, hessp=hessp, &
                        initial_trust_radius=initial_trust_radius, &
                        max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                        maxiter=maxiter, inexact=inexact)
                end if
            else
                res = trust_krylov(fun, x0, jac, hess=hess, &
                    initial_trust_radius=initial_trust_radius, &
                    max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                    maxiter=maxiter, inexact=inexact)
            end if
        case ('trust-exact', 'trust_exact')
            gtol_use = 1.0e-4_dp
            if (present(tol)) gtol_use = tol
            if (present(gtol)) gtol_use = gtol
            if (.not. present(jac) .or. .not. present(hess)) then
                call init_result(res, 'trust-exact', size(x0))
                return
            end if
            res = trust_exact(fun, x0, jac, hess, &
                initial_trust_radius=initial_trust_radius, &
                max_trust_radius=max_trust_radius, eta=eta, gtol=gtol_use, &
                maxiter=maxiter, subproblem_maxiter=subproblem_maxiter)
        case ('trust-constr', 'trust_constr')
            gtol_use = 1.0e-8_dp
            xtol_use = 1.0e-8_dp
            barrier_tol_use = 1.0e-8_dp
            if (present(tol)) then
                gtol_use = tol
                xtol_use = tol
                barrier_tol_use = tol
            end if
            if (present(gtol)) gtol_use = gtol
            if (present(xtol)) xtol_use = xtol
            if (present(barrier_tol)) barrier_tol_use = barrier_tol
            res = trust_constr(fun, x0, jac=jac, hess=hess, lower=lower, upper=upper, &
                constraints=constraints, constraint_jac=constraint_jac, meq=meq, &
                mineq=mineq, gtol=gtol_use, xtol=xtol_use, eps=eps, maxiter=maxiter, &
                initial_trust_radius=initial_trust_radius, &
                initial_constr_penalty=initial_constr_penalty, barrier_tol=barrier_tol_use, &
                initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                initial_tr_radius=initial_tr_radius, constraint_hess=constraint_hess, &
                constraint_lower=constraint_lower, constraint_upper=constraint_upper, &
                constraint_keep_feasible=constraint_keep_feasible, &
                bounds_keep_feasible=bounds_keep_feasible, constraint_objects=constraint_objects, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, hessp=hessp, &
                objective_hess_update=objective_hess_update, callback=callback, &
                objective_hess_fd_scheme=objective_hess_fd_scheme)
        case ('cobyla')
            res = cobyla(fun, x0, lower=lower, upper=upper, constraints=constraints, &
                meq=meq, mineq=mineq, rhobeg=rhobeg, tol=tol, maxiter=maxiter, &
                catol=catol, f_target=f_target)
        case ('tnc')
            ftol_use = 0.0_dp
            xtol_use = sqrt(epsilon(1.0_dp))
            gtol_use = 1.0e-2_dp * sqrt(sqrt(epsilon(1.0_dp)))
            if (present(tol)) then
                ftol_use = tol
                xtol_use = tol
                gtol_use = tol
            end if
            if (present(ftol)) ftol_use = ftol
            if (present(xtol)) xtol_use = xtol
            if (present(gtol)) gtol_use = gtol
            if (present(jac)) then
                res = tnc(fun, x0, jac=jac, lower=lower, upper=upper, ftol=ftol_use, &
                    xtol=xtol_use, gtol=gtol_use, eps=eps, maxfun=maxfun, &
                    max_cg_iter=max_cg_iter)
            else
                res = tnc(fun, x0, lower=lower, upper=upper, ftol=ftol_use, &
                    xtol=xtol_use, gtol=gtol_use, eps=eps, maxfun=maxfun, &
                    max_cg_iter=max_cg_iter)
            end if
        case ('l-bfgs-b', 'lbfgsb', 'l_bfgs_b')
            gtol_use = default_gtol
            ftol_use = 2.2204460492503131e-9_dp
            if (present(tol)) then
                gtol_use = tol
                ftol_use = tol
            end if
            if (present(gtol)) gtol_use = gtol
            if (present(ftol)) ftol_use = ftol
            if (present(jac)) then
                res = lbfgsb(fun, x0, jac=jac, lower=lower, upper=upper, &
                    ftol=ftol_use, gtol=gtol_use, eps=eps, maxiter=maxiter, &
                    maxfun=maxfun, maxcor=maxcor)
            else
                res = lbfgsb(fun, x0, lower=lower, upper=upper, ftol=ftol_use, &
                    gtol=gtol_use, eps=eps, maxiter=maxiter, maxfun=maxfun, &
                    maxcor=maxcor)
            end if
        case ('slsqp')
            ftol_use = 1.0e-6_dp
            if (present(tol)) ftol_use = tol
            if (present(ftol)) ftol_use = ftol
            if (present(jac)) then
                if (present(constraint_jac)) then
                    res = slsqp(fun, x0, jac=jac, lower=lower, upper=upper, &
                        constraints=constraints, constraint_jac=constraint_jac, &
                        meq=meq, mineq=mineq, ftol=ftol_use, eps=eps, maxiter=maxiter)
                else
                    res = slsqp(fun, x0, jac=jac, lower=lower, upper=upper, &
                        constraints=constraints, meq=meq, mineq=mineq, ftol=ftol_use, &
                        eps=eps, maxiter=maxiter)
                end if
            else
                if (present(constraint_jac)) then
                    res = slsqp(fun, x0, lower=lower, upper=upper, &
                        constraints=constraints, constraint_jac=constraint_jac, &
                        meq=meq, mineq=mineq, ftol=ftol_use, eps=eps, maxiter=maxiter)
                else
                    res = slsqp(fun, x0, lower=lower, upper=upper, &
                        constraints=constraints, meq=meq, mineq=mineq, ftol=ftol_use, &
                        eps=eps, maxiter=maxiter)
                end if
            end if
        case default
            call init_result(res, meth, size(x0))
        end select
    end function minimize

    function nelder_mead(fun, x0, xatol, fatol, maxiter, maxfev, adaptive, &
            initial_simplex, lower, upper) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        real(dp), intent(in), optional :: xatol !! absolute simplex-coordinate tolerance
        real(dp), intent(in), optional :: fatol !! absolute simplex-objective tolerance
        integer, intent(in), optional :: maxiter !! maximum simplex iterations
        integer, intent(in), optional :: maxfev !! maximum objective evaluations
        logical, intent(in), optional :: adaptive !! dimension-adaptive coefficients
        real(dp), intent(in), optional :: initial_simplex(:, :) !! shape (n+1,n)
        real(dp), intent(in), optional :: lower(:) !! lower bounds; nonfinite means unbounded
        real(dp), intent(in), optional :: upper(:) !! upper bounds; nonfinite means unbounded
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: sim(:, :)
        real(dp), allocatable :: fsim(:)
        real(dp), allocatable :: lb(:)
        real(dp), allocatable :: ub(:)
        real(dp), allocatable :: xstart(:)
        real(dp), allocatable :: xbar(:)
        real(dp), allocatable :: xr(:)
        real(dp), allocatable :: xe(:)
        real(dp), allocatable :: xc(:)
        real(dp) :: chi
        real(dp) :: fatol_use
        real(dp) :: fxc
        real(dp) :: fxe
        real(dp) :: fxr
        real(dp) :: psi
        real(dp) :: rho
        real(dp) :: sigma
        real(dp) :: xatol_use
        integer :: i
        integer :: j
        integer :: maxfev_use
        integer :: maxiter_use
        integer :: n
        integer, allocatable :: nbd(:)
        logical :: adaptive_use
        logical :: doshrink
        logical :: invalid_bounds

        n = size(x0)
        call init_result(res, 'nelder-mead', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        allocate(lb(n), ub(n), nbd(n), xstart(n))
        call setup_box_bounds(n, lower, upper, lb, ub, nbd, invalid_bounds)
        if (invalid_bounds) return
        xstart = x0
        call clip_to_box(xstart, lb, ub, nbd)

        xatol_use = default_nm_xatol
        if (present(xatol)) xatol_use = xatol
        fatol_use = default_nm_fatol
        if (present(fatol)) fatol_use = fatol
        if (.not. ieee_is_finite(xatol_use) .or. xatol_use < 0.0_dp) return
        if (.not. ieee_is_finite(fatol_use) .or. fatol_use < 0.0_dp) return

        call resolve_limits(n, 200, maxiter, maxfev, maxiter_use, maxfev_use)
        if (maxiter_use < 0 .or. maxfev_use < 0) return

        adaptive_use = .false.
        if (present(adaptive)) adaptive_use = adaptive
        rho = 1.0_dp
        if (adaptive_use) then
            chi = 1.0_dp + 2.0_dp / real(n, dp)
            psi = 0.75_dp - 1.0_dp / (2.0_dp * real(n, dp))
            sigma = 1.0_dp - 1.0_dp / real(n, dp)
        else
            chi = 2.0_dp
            psi = 0.5_dp
            sigma = 0.5_dp
        end if

        allocate(sim(n + 1, n), fsim(n + 1), xbar(n), xr(n), xe(n), xc(n))
        if (present(initial_simplex)) then
            if (size(initial_simplex, 1) /= n + 1 .or. size(initial_simplex, 2) /= n) return
            sim = initial_simplex
            if (.not. all(ieee_is_finite(sim))) return
        else
            sim(1, :) = xstart
            do i = 1, n
                sim(i + 1, :) = xstart
                if (abs(xstart(i)) > 0.0_dp) then
                    sim(i + 1, i) = 1.05_dp * xstart(i)
                else
                    sim(i + 1, i) = 2.5e-4_dp
                end if
            end do
        end if

        if (any(nbd /= 0)) then
            do j = 1, n + 1
                do i = 1, n
                    if ((nbd(i) == 2 .or. nbd(i) == 3) .and. sim(j, i) > ub(i)) then
                        sim(j, i) = 2.0_dp * ub(i) - sim(j, i)
                    end if
                end do
                call clip_to_box(sim(j, :), lb, ub, nbd)
            end do
        end if

        fsim = huge(1.0_dp)
        do i = 1, n + 1
            if (res%nfev >= maxfev_use) exit
            fsim(i) = fun(sim(i, :))
            res%nfev = res%nfev + 1
        end do
        call sort_simplex(sim, fsim)
        res%nit = 1

        do while (res%nfev < maxfev_use .and. res%nit < maxiter_use)
            if (max_simplex_dx(sim) <= xatol_use .and. &
                    maxval(abs(fsim(1) - fsim(2:))) <= fatol_use) exit

            xbar = sum(sim(1:n, :), dim=1) / real(n, dp)
            xr = (1.0_dp + rho) * xbar - rho * sim(n + 1, :)
            call clip_to_box(xr, lb, ub, nbd)
            if (.not. eval_nm(fun, xr, maxfev_use, res%nfev, fxr)) exit
            doshrink = .false.

            if (fxr < fsim(1)) then
                xe = (1.0_dp + rho * chi) * xbar - rho * chi * sim(n + 1, :)
                call clip_to_box(xe, lb, ub, nbd)
                if (.not. eval_nm(fun, xe, maxfev_use, res%nfev, fxe)) exit
                if (fxe < fxr) then
                    sim(n + 1, :) = xe
                    fsim(n + 1) = fxe
                else
                    sim(n + 1, :) = xr
                    fsim(n + 1) = fxr
                end if
            else if (fxr < fsim(n)) then
                sim(n + 1, :) = xr
                fsim(n + 1) = fxr
            else
                if (fxr < fsim(n + 1)) then
                    xc = (1.0_dp + psi * rho) * xbar - psi * rho * sim(n + 1, :)
                    call clip_to_box(xc, lb, ub, nbd)
                    if (.not. eval_nm(fun, xc, maxfev_use, res%nfev, fxc)) exit
                    if (fxc <= fxr) then
                        sim(n + 1, :) = xc
                        fsim(n + 1) = fxc
                    else
                        doshrink = .true.
                    end if
                else
                    xc = (1.0_dp - psi) * xbar + psi * sim(n + 1, :)
                    call clip_to_box(xc, lb, ub, nbd)
                    if (.not. eval_nm(fun, xc, maxfev_use, res%nfev, fxc)) exit
                    if (fxc < fsim(n + 1)) then
                        sim(n + 1, :) = xc
                        fsim(n + 1) = fxc
                    else
                        doshrink = .true.
                    end if
                end if

                if (doshrink) then
                    do j = 2, n + 1
                        sim(j, :) = sim(1, :) + sigma * (sim(j, :) - sim(1, :))
                        call clip_to_box(sim(j, :), lb, ub, nbd)
                        if (.not. eval_nm(fun, sim(j, :), maxfev_use, res%nfev, fsim(j))) exit
                    end do
                    if (j <= n + 1) exit
                end if
            end if

            res%nit = res%nit + 1
            call sort_simplex(sim, fsim)
        end do

        call sort_simplex(sim, fsim)
        res%x = sim(1, :)
        res%fun = fsim(1)
        allocate(res%final_simplex(n + 1, n), res%final_simplex_fun(n + 1))
        res%final_simplex = sim
        res%final_simplex_fun = fsim

        if (res%nfev >= maxfev_use) then
            call finish_result(res, 1, .false., 'Maximum number of function evaluations has been exceeded.')
        else if (res%nit >= maxiter_use) then
            call finish_result(res, 2, .false., 'Maximum number of iterations has been exceeded.')
        else if (ieee_is_nan(res%fun) .or. any(ieee_is_nan(res%x))) then
            call finish_result(res, 3, .false., 'NaN result encountered.')
        else
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if
    end function nelder_mead

    function powell(fun, x0, xtol, ftol, maxiter, maxfev, direc, lower, upper) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        real(dp), intent(in), optional :: xtol !! relative line-search tolerance
        real(dp), intent(in), optional :: ftol !! relative objective tolerance
        integer, intent(in), optional :: maxiter !! maximum outer iterations
        integer, intent(in), optional :: maxfev !! maximum objective evaluations
        real(dp), intent(in), optional :: direc(:, :) !! initial direction rows, shape (n,n)
        real(dp), intent(in), optional :: lower(:) !! lower bounds; nonfinite means unbounded
        real(dp), intent(in), optional :: upper(:) !! upper bounds; nonfinite means unbounded
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: dirs(:, :)
        real(dp), allocatable :: direction(:)
        real(dp), allocatable :: lb(:)
        real(dp), allocatable :: ub(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x1(:)
        real(dp), allocatable :: x2(:)
        real(dp) :: bnd
        real(dp) :: delta
        real(dp) :: ftol_use
        real(dp) :: fval
        real(dp) :: fx
        real(dp) :: fx2
        real(dp) :: temp
        real(dp) :: t
        real(dp) :: xtol_use
        integer :: bigind
        integer :: i
        integer :: maxfev_use
        integer :: maxiter_use
        integer :: n
        integer, allocatable :: nbd(:)
        logical :: invalid_bounds
        logical :: line_ok

        n = size(x0)
        call init_result(res, 'powell', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        allocate(lb(n), ub(n), nbd(n))
        call setup_box_bounds(n, lower, upper, lb, ub, nbd, invalid_bounds)
        if (invalid_bounds) return

        xtol_use = default_powell_xtol
        if (present(xtol)) xtol_use = xtol
        ftol_use = default_powell_ftol
        if (present(ftol)) ftol_use = ftol
        if (.not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(ftol_use) .or. ftol_use < 0.0_dp) return
        call resolve_limits(n, 1000, maxiter, maxfev, maxiter_use, maxfev_use)
        if (maxiter_use < 0 .or. maxfev_use < 0) return

        allocate(x(n), x1(n), x2(n), direction(n), dirs(n, n))
        x = x0
        dirs = 0.0_dp
        do i = 1, n
            dirs(i, i) = 1.0_dp
        end do
        if (present(direc)) then
            if (size(direc, 1) /= n .or. size(direc, 2) /= n) return
            if (.not. all(ieee_is_finite(direc))) return
            dirs = direc
        end if

        if (.not. eval_powell(fun, x, maxfev_use, res%nfev, fval)) then
            res%x = x
            res%fun = ieee_value(0.0_dp, ieee_quiet_nan)
            call finish_result(res, 1, .false., 'Maximum number of function evaluations has been exceeded.')
            return
        end if
        x1 = x
        res%nit = 0

        do
            fx = fval
            bigind = 1
            delta = 0.0_dp
            do i = 1, n
                direction = dirs(i, :)
                fx2 = fval
                call powell_line_search(fun, x, direction, xtol_use * 100.0_dp, &
                    maxfev_use, res%nfev, fval, line_ok, lb, ub, nbd)
                if (.not. line_ok) exit
                if (fx2 - fval > delta) then
                    delta = fx2 - fval
                    bigind = i
                end if
            end do
            if (i <= n) exit

            res%nit = res%nit + 1
            bnd = ftol_use * (abs(fx) + abs(fval)) + 1.0e-20_dp
            if (2.0_dp * (fx - fval) <= bnd) exit
            if (res%nfev >= maxfev_use) exit
            if (res%nit >= maxiter_use) exit
            if (ieee_is_nan(fx) .and. ieee_is_nan(fval)) exit

            direction = x - x1
            x1 = x
            call powell_extrapolated_point(x, direction, lb, ub, nbd, x2)
            if (.not. eval_powell(fun, x2, maxfev_use, res%nfev, fx2)) exit

            if (fx > fx2) then
                t = 2.0_dp * (fx + fx2 - 2.0_dp * fval)
                temp = fx - fval - delta
                t = t * temp * temp
                temp = fx - fx2
                t = t - delta * temp * temp
                if (t < 0.0_dp) then
                    call powell_line_search(fun, x, direction, xtol_use * 100.0_dp, &
                        maxfev_use, res%nfev, fval, line_ok, lb, ub, nbd)
                    if (.not. line_ok) exit
                    if (any(abs(direction) > 0.0_dp)) then
                        dirs(bigind, :) = dirs(n, :)
                        dirs(n, :) = direction
                    end if
                end if
            end if
        end do

        res%x = x
        res%fun = fval
        allocate(res%direc(n, n))
        res%direc = dirs
        if (.not. point_in_box(res%x, lb, ub, nbd)) then
            call finish_result(res, 4, .false., 'The result is outside of the provided bounds.')
        else if (res%nfev >= maxfev_use) then
            call finish_result(res, 1, .false., 'Maximum number of function evaluations has been exceeded.')
        else if (res%nit >= maxiter_use) then
            call finish_result(res, 2, .false., 'Maximum number of iterations has been exceeded.')
        else if (ieee_is_nan(res%fun) .or. any(ieee_is_nan(res%x))) then
            call finish_result(res, 3, .false., 'NaN result encountered.')
        else
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if
    end function powell

    function bfgs(fun, x0, jac, gtol, eps, maxiter, xrtol, c1, c2, hess_inv0) &
            result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic gradient
        real(dp), intent(in), optional :: gtol !! max-norm gradient tolerance
        real(dp), intent(in), optional :: eps !! forward-difference absolute step
        integer, intent(in), optional :: maxiter !! maximum BFGS iterations
        real(dp), intent(in), optional :: xrtol !! relative step tolerance
        real(dp), intent(in), optional :: c1 !! strong-Wolfe sufficient-decrease constant
        real(dp), intent(in), optional :: c2 !! strong-Wolfe curvature constant
        real(dp), intent(in), optional :: hess_inv0(:, :) !! initial inverse Hessian
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: a1(:, :)
        real(dp), allocatable :: a2(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: gnew(:)
        real(dp), allocatable :: h(:, :)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xnew(:)
        real(dp), allocatable :: y(:)
        real(dp) :: alpha
        real(dp) :: c1_use
        real(dp) :: c2_use
        real(dp) :: eps_use
        real(dp) :: fval
        real(dp) :: gnorm
        real(dp) :: gtol_use
        real(dp) :: old_old_fval
        real(dp) :: rhok
        real(dp) :: rhok_inv
        real(dp) :: xrtol_use
        integer :: i
        integer :: maxiter_use
        integer :: n
        logical :: line_ok

        real(dp), allocatable :: cache_x(:)
        real(dp), allocatable :: cache_g(:)
        real(dp) :: cache_f
        logical :: cache_f_valid
        logical :: cache_g_valid

        n = size(x0)
        call init_result(res, 'bfgs', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        call gradient_options(n, gtol, eps, maxiter, c1, c2, gtol_use, eps_use, &
            maxiter_use, c1_use, c2_use)
        if (maxiter_use < 0) return
        xrtol_use = 0.0_dp
        if (present(xrtol)) xrtol_use = xrtol
        if (.not. ieee_is_finite(xrtol_use) .or. xrtol_use < 0.0_dp) return

        allocate(x(n), xnew(n), g(n), gnew(n), p(n), s(n), y(n), h(n, n), &
            a1(n, n), a2(n, n), cache_x(n), cache_g(n))
        x = x0
        h = 0.0_dp
        do i = 1, n
            h(i, i) = 1.0_dp
        end do
        if (present(hess_inv0)) then
            if (size(hess_inv0, 1) /= n .or. size(hess_inv0, 2) /= n) return
            if (.not. all(ieee_is_finite(hess_inv0))) return
            h = hess_inv0
        end if
        cache_f_valid = .false.
        cache_g_valid = .false.

        call get_fun(x, fval)
        call get_grad(x, g)
        if (.not. ieee_is_finite(fval) .or. .not. all(ieee_is_finite(g))) then
            res%x = x
            res%fun = fval
            res%jac = g
            call finish_result(res, 3, .false., 'NaN result encountered.')
            return
        end if
        old_old_fval = fval + norm2(g) / 2.0_dp
        gnorm = maxval(abs(g))

        do while (gnorm > gtol_use .and. res%nit < maxiter_use)
            p = -matmul(h, g)
            call vector_wolfe_search(p, g, fval, old_old_fval, c1_use, c2_use, &
                alpha, line_ok)
            if (.not. line_ok) then
                call finish_result(res, 2, .false., &
                    'Desired error not necessarily achieved due to precision loss.')
                exit
            end if

            s = alpha * p
            xnew = x + s
            old_old_fval = fval
            call get_fun(xnew, fval)
            call get_grad(xnew, gnew)
            y = gnew - g
            x = xnew
            g = gnew
            res%nit = res%nit + 1
            gnorm = maxval(abs(g))
            if (gnorm <= gtol_use) exit
            if (alpha * norm2(p) <= xrtol_use * (xrtol_use + norm2(x))) exit
            if (.not. ieee_is_finite(fval)) then
                call finish_result(res, 2, .false., &
                    'Desired error not necessarily achieved due to precision loss.')
                exit
            end if

            rhok_inv = dot_product(y, s)
            if (same_real(rhok_inv, 0.0_dp)) then
                rhok = 1000.0_dp
            else
                rhok = 1.0_dp / rhok_inv
            end if
            call identity_minus_outer(s, y, rhok, a1)
            call identity_minus_outer(y, s, rhok, a2)
            h = matmul(a1, matmul(h, a2)) + rhok * outer_product(s, s)
        end do

        res%x = x
        res%fun = fval
        res%jac = g
        allocate(res%hess_inv(n, n))
        res%hess_inv = h
        if (res%status == 2) return
        if (res%nit >= maxiter_use .and. gnorm > gtol_use) then
            call finish_result(res, 1, .false., 'Maximum number of iterations has been exceeded.')
        else if (.not. ieee_is_finite(gnorm) .or. .not. ieee_is_finite(fval) .or. &
                any(.not. ieee_is_finite(x))) then
            call finish_result(res, 3, .false., 'NaN result encountered.')
        else
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine get_fun(xv, fv)
            real(dp), intent(in) :: xv(:) !! point to evaluate
            real(dp), intent(out) :: fv !! objective value

            if (cache_f_valid .and. same_vector(xv, cache_x)) then
                fv = cache_f
                return
            end if
            fv = fun(xv)
            res%nfev = res%nfev + 1
            cache_x = xv
            cache_f = fv
            cache_f_valid = .true.
            cache_g_valid = .false.
        end subroutine get_fun

        subroutine get_grad(xv, gv)
            real(dp), intent(in) :: xv(:) !! point to differentiate
            real(dp), intent(out) :: gv(:) !! objective gradient

            real(dp), allocatable :: xp(:)
            real(dp) :: f0
            real(dp) :: fp
            integer :: k

            if (cache_g_valid .and. same_vector(xv, cache_x)) then
                gv = cache_g
                return
            end if
            if (present(jac)) then
                call invoke_gradient(jac, xv, gv)
                res%njev = res%njev + 1
                if (.not. cache_f_valid .or. .not. same_vector(xv, cache_x)) then
                    cache_x = xv
                    cache_f_valid = .false.
                end if
            else
                call get_fun(xv, f0)
                allocate(xp(n))
                do k = 1, n
                    xp = xv
                    xp(k) = xp(k) + eps_use
                    fp = fun(xp)
                    res%nfev = res%nfev + 1
                    gv(k) = (fp - f0) / eps_use
                end do
                res%njev = res%njev + 1
                cache_x = xv
                cache_f = f0
                cache_f_valid = .true.
            end if
            cache_g = gv
            cache_g_valid = .true.
        end subroutine get_grad

        subroutine vector_wolfe_search(direction, gbase, phi0, old_phi0, &
                c1v, c2v, step, ok)
            real(dp), intent(in) :: direction(:) !! search direction
            real(dp), intent(in) :: gbase(:) !! gradient at current point
            real(dp), intent(in) :: phi0 !! objective at current point
            real(dp), intent(in) :: old_phi0 !! objective at preceding outer point
            real(dp), intent(in) :: c1v !! Armijo constant
            real(dp), intent(in) :: c2v !! curvature constant
            real(dp), intent(out) :: step !! accepted step
            logical, intent(out) :: ok !! true when a Wolfe step is found

            real(dp) :: phi1
            real(dp) :: dphi0

            dphi0 = dot_product(gbase, direction)
            call scalar_search_wolfe1(line_phi, line_derphi, phi0, step, phi1, ok, &
                derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, amin=1.0e-100_dp, &
                old_phi0=old_phi0)
            if (.not. ok) then
                call scalar_search_wolfe2(line_phi, line_derphi, phi0, step, phi1, ok, &
                    derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, &
                    old_phi0=old_phi0)
            end if
        end subroutine vector_wolfe_search

        function line_phi(step) result(value)
            real(dp), intent(in) :: step !! trial line step
            real(dp) :: value !! objective along current line
            real(dp) :: xt(n)

            xt = x + step * p
            call get_fun(xt, value)
        end function line_phi

        function line_derphi(step) result(value)
            real(dp), intent(in) :: step !! trial line step
            real(dp) :: value !! directional derivative
            real(dp) :: gt(n)
            real(dp) :: xt(n)

            xt = x + step * p
            call get_grad(xt, gt)
            value = dot_product(gt, p)
        end function line_derphi

    end function bfgs

    function cg(fun, x0, jac, gtol, eps, maxiter, c1, c2) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic gradient
        real(dp), intent(in), optional :: gtol !! max-norm gradient tolerance
        real(dp), intent(in), optional :: eps !! forward-difference absolute step
        integer, intent(in), optional :: maxiter !! maximum CG iterations
        real(dp), intent(in), optional :: c1 !! strong-Wolfe sufficient-decrease constant
        real(dp), intent(in), optional :: c2 !! strong-Wolfe curvature constant
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: g(:)
        real(dp), allocatable :: gnew(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: pnew(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xnew(:)
        real(dp), allocatable :: y(:)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: c1_use
        real(dp) :: c2_use
        real(dp) :: delta
        real(dp) :: eps_use
        real(dp) :: fval
        real(dp) :: gnorm
        real(dp) :: gtol_use
        real(dp) :: old_old_fval
        integer :: maxiter_use
        integer :: n
        logical :: descent_ok
        logical :: line_ok

        real(dp), allocatable :: cache_x(:)
        real(dp), allocatable :: cache_g(:)
        real(dp) :: cache_f
        logical :: cache_f_valid
        logical :: cache_g_valid

        n = size(x0)
        call init_result(res, 'cg', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        call gradient_options(n, gtol, eps, maxiter, c1, c2, gtol_use, eps_use, &
            maxiter_use, c1_use, c2_use, default_c2=0.4_dp)
        if (maxiter_use < 0) return

        allocate(x(n), xnew(n), g(n), gnew(n), p(n), pnew(n), y(n), &
            cache_x(n), cache_g(n))
        x = x0
        cache_f_valid = .false.
        cache_g_valid = .false.
        call get_fun(x, fval)
        call get_grad(x, g)
        if (.not. ieee_is_finite(fval) .or. .not. all(ieee_is_finite(g))) then
            res%x = x
            res%fun = fval
            res%jac = g
            call finish_result(res, 3, .false., 'NaN result encountered.')
            return
        end if
        old_old_fval = fval + norm2(g) / 2.0_dp
        p = -g
        gnorm = maxval(abs(g))

        do while (gnorm > gtol_use .and. res%nit < maxiter_use)
            delta = dot_product(g, g)
            call vector_wolfe_search(p, g, fval, old_old_fval, c1_use, c2_use, &
                alpha, line_ok)
            if (.not. line_ok) then
                call finish_result(res, 2, .false., &
                    'Desired error not necessarily achieved due to precision loss.')
                exit
            end if

            xnew = x + alpha * p
            old_old_fval = fval
            call get_fun(xnew, fval)
            call get_grad(xnew, gnew)
            y = gnew - g
            beta = max(0.0_dp, dot_product(y, gnew) / delta)
            pnew = -gnew + beta * p
            gnorm = maxval(abs(gnew))
            descent_ok = gnorm <= gtol_use .or. &
                dot_product(pnew, gnew) <= -0.01_dp * dot_product(gnew, gnew)
            if (.not. descent_ok) then
                pnew = -gnew
            end if
            x = xnew
            g = gnew
            p = pnew
            res%nit = res%nit + 1
        end do

        res%x = x
        res%fun = fval
        res%jac = g
        if (res%status == 2) return
        if (res%nit >= maxiter_use .and. gnorm > gtol_use) then
            call finish_result(res, 1, .false., 'Maximum number of iterations has been exceeded.')
        else if (.not. ieee_is_finite(gnorm) .or. .not. ieee_is_finite(fval) .or. &
                any(.not. ieee_is_finite(x))) then
            call finish_result(res, 3, .false., 'NaN result encountered.')
        else
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine get_fun(xv, fv)
            real(dp), intent(in) :: xv(:) !! point to evaluate
            real(dp), intent(out) :: fv !! objective value

            if (cache_f_valid .and. same_vector(xv, cache_x)) then
                fv = cache_f
                return
            end if
            fv = fun(xv)
            res%nfev = res%nfev + 1
            cache_x = xv
            cache_f = fv
            cache_f_valid = .true.
            cache_g_valid = .false.
        end subroutine get_fun

        subroutine get_grad(xv, gv)
            real(dp), intent(in) :: xv(:) !! point to differentiate
            real(dp), intent(out) :: gv(:) !! objective gradient

            real(dp), allocatable :: xp(:)
            real(dp) :: f0
            real(dp) :: fp
            integer :: k

            if (cache_g_valid .and. same_vector(xv, cache_x)) then
                gv = cache_g
                return
            end if
            if (present(jac)) then
                call invoke_gradient(jac, xv, gv)
                res%njev = res%njev + 1
                if (.not. cache_f_valid .or. .not. same_vector(xv, cache_x)) then
                    cache_x = xv
                    cache_f_valid = .false.
                end if
            else
                call get_fun(xv, f0)
                allocate(xp(n))
                do k = 1, n
                    xp = xv
                    xp(k) = xp(k) + eps_use
                    fp = fun(xp)
                    res%nfev = res%nfev + 1
                    gv(k) = (fp - f0) / eps_use
                end do
                res%njev = res%njev + 1
                cache_x = xv
                cache_f = f0
                cache_f_valid = .true.
            end if
            cache_g = gv
            cache_g_valid = .true.
        end subroutine get_grad

        subroutine vector_wolfe_search(direction, gbase, phi0, old_phi0, &
                c1v, c2v, step, ok)
            real(dp), intent(in) :: direction(:) !! search direction
            real(dp), intent(in) :: gbase(:) !! gradient at current point
            real(dp), intent(in) :: phi0 !! objective at current point
            real(dp), intent(in) :: old_phi0 !! objective at preceding outer point
            real(dp), intent(in) :: c1v !! Armijo constant
            real(dp), intent(in) :: c2v !! curvature constant
            real(dp), intent(out) :: step !! accepted step
            logical, intent(out) :: ok !! true when a Wolfe step is found

            real(dp) :: phi1
            real(dp) :: dphi0

            dphi0 = dot_product(gbase, direction)
            call scalar_search_wolfe1(line_phi, line_derphi, phi0, step, phi1, ok, &
                derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, amin=1.0e-100_dp, &
                old_phi0=old_phi0)
            if (ok) ok = cg_extra_condition(step, phi1)
            if (.not. ok) then
                call scalar_search_wolfe2(line_phi, line_derphi, phi0, step, phi1, ok, &
                    derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, &
                    old_phi0=old_phi0, extra_condition=cg_extra_condition)
            end if
        end subroutine vector_wolfe_search

        logical function cg_extra_condition(step, phi_value) result(accept)
            real(dp), intent(in) :: step !! candidate strong-Wolfe step
            real(dp), intent(in) :: phi_value !! objective at candidate step

            real(dp) :: beta_trial
            real(dp) :: gt(n)
            real(dp) :: pt(n)
            real(dp) :: xt(n)
            real(dp) :: yt(n)

            if (.not. ieee_is_finite(phi_value)) then
                accept = .false.
                return
            end if
            xt = x + step * p
            call get_grad(xt, gt)
            yt = gt - g
            beta_trial = max(0.0_dp, dot_product(yt, gt) / dot_product(g, g))
            pt = -gt + beta_trial * p
            accept = maxval(abs(gt)) <= gtol_use .or. &
                dot_product(pt, gt) <= -0.01_dp * dot_product(gt, gt)
        end function cg_extra_condition

        function line_phi(step) result(value)
            real(dp), intent(in) :: step !! trial line step
            real(dp) :: value !! objective along current line
            real(dp) :: xt(n)

            xt = x + step * p
            call get_fun(xt, value)
        end function line_phi

        function line_derphi(step) result(value)
            real(dp), intent(in) :: step !! trial line step
            real(dp) :: value !! directional derivative
            real(dp) :: gt(n)
            real(dp) :: xt(n)

            xt = x + step * p
            call get_grad(xt, gt)
            value = dot_product(gt, p)
        end function line_derphi

    end function cg

    function newton_cg(fun, x0, jac, hess, hessp, xtol, eps, maxiter, c1, c2) &
            result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function) :: jac !! analytic gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic Hessian
        procedure(multivariate_hessian_product_function), optional :: hessp !! Hessian product
        real(dp), intent(in), optional :: xtol !! average relative step tolerance
        real(dp), intent(in), optional :: eps !! finite-difference Hessian-product step
        integer, intent(in), optional :: maxiter !! maximum Newton iterations
        real(dp), intent(in), optional :: c1 !! strong-Wolfe sufficient-decrease constant
        real(dp), intent(in), optional :: c2 !! strong-Wolfe curvature constant
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: ap(:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: gbase(:)
        real(dp), allocatable :: gshift(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: update(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xsup(:)
        real(dp) :: alpha
        real(dp) :: alphai
        real(dp) :: beta
        real(dp) :: c1_use
        real(dp) :: c2_use
        real(dp) :: curv
        real(dp) :: dri0
        real(dp) :: dri1
        real(dp) :: eps_use
        real(dp) :: eta
        real(dp) :: fval
        real(dp) :: maggrad
        real(dp) :: line_fval
        real(dp) :: old_fval
        real(dp) :: old_old_fval
        real(dp) :: termcond
        real(dp) :: xtol_use
        integer :: i
        integer :: k2
        integer :: maxiter_use
        integer :: n
        logical :: have_old_old
        logical :: inner_ok
        logical :: line_ok

        real(dp), allocatable :: cache_x(:)
        real(dp), allocatable :: cache_g(:)
        real(dp) :: cache_f
        logical :: cache_f_valid
        logical :: cache_g_valid

        n = size(x0)
        call init_result(res, 'newton-cg', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        xtol_use = 1.0e-5_dp
        if (present(xtol)) xtol_use = xtol
        eps_use = default_fd_eps
        if (present(eps)) eps_use = eps
        maxiter_use = 200 * n
        if (present(maxiter)) maxiter_use = maxiter
        c1_use = 1.0e-4_dp
        if (present(c1)) c1_use = c1
        c2_use = 0.9_dp
        if (present(c2)) c2_use = c2
        if (.not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp) return
        if (maxiter_use < 0) return
        if (.not. ieee_is_finite(c1_use) .or. .not. ieee_is_finite(c2_use)) return
        if (c1_use <= 0.0_dp .or. c1_use >= c2_use .or. c2_use >= 1.0_dp) return

        allocate(x(n), g(n), b(n), r(n), p(n), xsup(n), ap(n), update(n), &
            gbase(n), gshift(n), cache_x(n), cache_g(n))
        if (present(hess)) allocate(a(n, n))
        x = x0
        cache_f_valid = .false.
        cache_g_valid = .false.
        have_old_old = .false.
        update = huge(1.0_dp)

        call get_fun(x, fval)
        if (.not. ieee_is_finite(fval)) then
            res%x = x
            res%fun = fval
            call finish_result(res, 3, .false., 'NaN result encountered.')
            return
        end if

        do while (sum(abs(update)) > real(n, dp) * xtol_use)
            if (res%nit >= maxiter_use) then
                call finish_result(res, 1, .false., &
                    'Maximum number of iterations has been exceeded.')
                exit
            end if

            call get_grad(x, g)
            if (.not. all(ieee_is_finite(g))) then
                call finish_result(res, 3, .false., 'NaN result encountered.')
                exit
            end if
            b = -g
            maggrad = sum(abs(b))
            eta = min(0.5_dp, sqrt(maggrad))
            termcond = eta * maggrad
            xsup = 0.0_dp
            r = -b
            p = -r
            i = 0
            dri0 = dot_product(r, r)
            inner_ok = .false.

            if (present(hess)) then
                call hess(x, a)
                res%nhev = res%nhev + 1
                if (.not. all(ieee_is_finite(a))) then
                    call finish_result(res, 3, .false., 'NaN result encountered.')
                    exit
                end if
            end if

            do k2 = 1, 20 * n
                if (sum(abs(r)) <= termcond) then
                    inner_ok = .true.
                    exit
                end if
                if (present(hess)) then
                    ap = matmul(a, p)
                else if (present(hessp)) then
                    call hessp(x, p, ap)
                    res%nhev = res%nhev + 1
                else
                    call get_grad(x, gbase)
                    call get_grad(x + eps_use * p, gshift)
                    ap = (gshift - gbase) / eps_use
                end if
                if (.not. all(ieee_is_finite(ap))) then
                    call finish_result(res, 3, .false., 'NaN result encountered.')
                    inner_ok = .true.
                    exit
                end if
                curv = dot_product(p, ap)
                if (curv >= 0.0_dp .and. curv <= 3.0_dp * epsilon(1.0_dp)) then
                    inner_ok = .true.
                    exit
                else if (curv < 0.0_dp) then
                    if (i > 0) then
                        inner_ok = .true.
                        exit
                    end if
                    xsup = dri0 / (-curv) * b
                    inner_ok = .true.
                    exit
                end if
                alphai = dri0 / curv
                xsup = xsup + alphai * p
                r = r + alphai * ap
                dri1 = dot_product(r, r)
                beta = dri1 / dri0
                p = -r + beta * p
                i = i + 1
                dri0 = dri1
                if (k2 == 20 * n) inner_ok = .false.
            end do
            if (res%status == 3) exit
            if (.not. inner_ok) then
                call finish_result(res, 3, .false., &
                    'CG iterations did not converge; Hessian may not be positive definite.')
                exit
            end if

            p = xsup
            old_fval = fval
            call vector_wolfe_search(p, g, fval, old_old_fval, have_old_old, &
                c1_use, c2_use, alpha, line_fval, line_ok)
            if (.not. line_ok) then
                call finish_result(res, 2, .false., &
                    'Desired error not necessarily achieved due to precision loss.')
                exit
            end if
            old_old_fval = old_fval
            have_old_old = .true.
            fval = line_fval
            update = alpha * p
            x = x + update
            res%nit = res%nit + 1
        end do

        res%x = x
        res%fun = fval
        if (allocated(g)) then
            allocate(res%jac(n))
            res%jac = g
        end if
        if (res%status >= 0) return
        if (.not. ieee_is_finite(fval) .or. any(.not. ieee_is_finite(x)) .or. &
                any(ieee_is_nan(update))) then
            call finish_result(res, 3, .false., 'NaN result encountered.')
        else
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine get_fun(xv, fv)
            real(dp), intent(in) :: xv(:) !! point to evaluate
            real(dp), intent(out) :: fv !! objective value

            if (cache_f_valid .and. same_vector(xv, cache_x)) then
                fv = cache_f
                return
            end if
            fv = fun(xv)
            res%nfev = res%nfev + 1
            cache_x = xv
            cache_f = fv
            cache_f_valid = .true.
            cache_g_valid = .false.
        end subroutine get_fun

        subroutine get_grad(xv, gv)
            real(dp), intent(in) :: xv(:) !! point to differentiate
            real(dp), intent(out) :: gv(:) !! objective gradient

            if (cache_g_valid .and. same_vector(xv, cache_x)) then
                gv = cache_g
                return
            end if
            call jac(xv, gv)
            res%njev = res%njev + 1
            if (.not. cache_f_valid .or. .not. same_vector(xv, cache_x)) then
                cache_x = xv
                cache_f_valid = .false.
            end if
            cache_g = gv
            cache_g_valid = .true.
        end subroutine get_grad

        subroutine vector_wolfe_search(direction, gbase0, phi0, previous_phi0, &
                have_previous, c1v, c2v, step, phi1, ok)
            real(dp), intent(in) :: direction(:) !! search direction
            real(dp), intent(in) :: gbase0(:) !! gradient at current point
            real(dp), intent(in) :: phi0 !! objective at current point
            real(dp), intent(in) :: previous_phi0 !! preceding objective, when available
            logical, intent(in) :: have_previous !! true when previous_phi0 is valid
            real(dp), intent(in) :: c1v !! Armijo constant
            real(dp), intent(in) :: c2v !! curvature constant
            real(dp), intent(out) :: step !! accepted step
            real(dp), intent(out) :: phi1 !! objective at accepted step
            logical, intent(out) :: ok !! true when a Wolfe step is found

            real(dp) :: dphi0

            dphi0 = dot_product(gbase0, direction)
            if (have_previous) then
                call scalar_search_wolfe1(line_phi, line_derphi, phi0, step, phi1, ok, &
                    derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, &
                    amin=1.0e-100_dp, old_phi0=previous_phi0)
            else
                call scalar_search_wolfe1(line_phi, line_derphi, phi0, step, phi1, ok, &
                    derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, &
                    amin=1.0e-100_dp)
            end if
            if (.not. ok) then
                if (have_previous) then
                    call scalar_search_wolfe2(line_phi, line_derphi, phi0, step, phi1, ok, &
                        derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp, &
                        old_phi0=previous_phi0)
                else
                    call scalar_search_wolfe2(line_phi, line_derphi, phi0, step, phi1, ok, &
                        derphi0=dphi0, c1=c1v, c2=c2v, amax=1.0e100_dp)
                end if
            end if
        end subroutine vector_wolfe_search

        function line_phi(step) result(value)
            real(dp), intent(in) :: step !! trial line step
            real(dp) :: value !! objective along current line
            real(dp) :: xt(n)

            xt = x + step * p
            call get_fun(xt, value)
        end function line_phi

        function line_derphi(step) result(value)
            real(dp), intent(in) :: step !! trial line step
            real(dp) :: value !! directional derivative
            real(dp) :: gt(n)
            real(dp) :: xt(n)

            xt = x + step * p
            call get_grad(xt, gt)
            value = dot_product(gt, p)
        end function line_derphi

    end function newton_cg

    function dogleg(fun, x0, jac, hess, initial_trust_radius, max_trust_radius, &
            eta, gtol, maxiter) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function) :: jac !! analytic gradient
        procedure(multivariate_hessian_function) :: hess !! analytic positive-definite Hessian
        real(dp), intent(in), optional :: initial_trust_radius !! initial trust-region radius
        real(dp), intent(in), optional :: max_trust_radius !! maximum trust-region radius
        real(dp), intent(in), optional :: eta !! step acceptance threshold in [0, 0.25)
        real(dp), intent(in), optional :: gtol !! Euclidean gradient-norm tolerance
        integer, intent(in), optional :: maxiter !! maximum trust-region iterations
        type(multivariate_minimize_result) :: res

        real(dp) :: actual_reduction
        real(dp) :: eta_use
        real(dp) :: f_proposed
        real(dp) :: fval
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: g_proposed(:)
        real(dp) :: gtol_use
        real(dp), allocatable :: hmat(:, :)
        real(dp) :: initial_radius_use
        real(dp) :: max_radius_use
        integer :: maxiter_use
        integer :: n
        real(dp), allocatable :: p(:)
        real(dp) :: predicted_reduction
        real(dp) :: predicted_value
        real(dp) :: rho
        real(dp) :: trust_radius
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_proposed(:)
        logical :: f_valid
        logical :: hess_valid
        logical :: hits_boundary
        logical :: ok

        n = size(x0)
        call init_result(res, 'dogleg', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        initial_radius_use = 1.0_dp
        if (present(initial_trust_radius)) initial_radius_use = initial_trust_radius
        max_radius_use = 1000.0_dp
        if (present(max_trust_radius)) max_radius_use = max_trust_radius
        eta_use = 0.15_dp
        if (present(eta)) eta_use = eta
        gtol_use = 1.0e-4_dp
        if (present(gtol)) gtol_use = gtol
        maxiter_use = 200 * n
        if (present(maxiter)) maxiter_use = maxiter

        if (.not. ieee_is_finite(initial_radius_use) .or. initial_radius_use <= 0.0_dp) return
        if (.not. ieee_is_finite(max_radius_use) .or. max_radius_use <= 0.0_dp) return
        if (initial_radius_use >= max_radius_use) return
        if (.not. ieee_is_finite(eta_use) .or. eta_use < 0.0_dp .or. eta_use >= 0.25_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return
        if (maxiter_use < 0) return

        allocate(x(n), x_proposed(n), g(n), g_proposed(n), p(n), hmat(n, n))
        x = x0
        trust_radius = initial_radius_use
        f_valid = .false.
        hess_valid = .false.

        call jac(x, g)
        res%njev = res%njev + 1
        if (.not. all(ieee_is_finite(g))) then
            res%x = x
            res%jac = g
            call finish_result(res, 3, .false., 'A linalg error occurred, such as a non-psd Hessian.')
            return
        end if

        do while (norm2(g) >= gtol_use)
            call solve_subproblem(p, hits_boundary, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            call model_value(p, predicted_value, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            x_proposed = x + p
            f_proposed = fun(x_proposed)
            res%nfev = res%nfev + 1
            if (.not. ieee_is_finite(f_proposed)) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            actual_reduction = fval - f_proposed
            predicted_reduction = fval - predicted_value
            if (predicted_reduction <= 0.0_dp) then
                call finish_result(res, 2, .false., &
                    'A bad approximation caused failure to predict improvement.')
                exit
            end if
            rho = actual_reduction / predicted_reduction

            if (rho < 0.25_dp) then
                trust_radius = 0.25_dp * trust_radius
            else if (rho > 0.75_dp .and. hits_boundary) then
                trust_radius = min(2.0_dp * trust_radius, max_radius_use)
            end if

            if (rho > eta_use) then
                x = x_proposed
                fval = f_proposed
                f_valid = .true.
                hess_valid = .false.
                call jac(x, g_proposed)
                res%njev = res%njev + 1
                if (.not. all(ieee_is_finite(g_proposed))) then
                    call finish_result(res, 3, .false., &
                        'A linalg error occurred, such as a non-psd Hessian.')
                    exit
                end if
                g = g_proposed
            end if

            res%nit = res%nit + 1
            if (norm2(g) < gtol_use) then
                call finish_result(res, 0, .true., 'Optimization terminated successfully.')
                exit
            end if
            if (res%nit >= maxiter_use) then
                call finish_result(res, 1, .false., &
                    'Maximum number of iterations has been exceeded.')
                exit
            end if
        end do

        if (.not. f_valid) then
            fval = fun(x)
            res%nfev = res%nfev + 1
            f_valid = .true.
        end if
        res%x = x
        res%fun = fval
        res%jac = g
        if (res%status < 0) then
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine ensure_fun(ok_eval)
            logical, intent(out) :: ok_eval !! true when the cached objective is finite

            if (.not. f_valid) then
                fval = fun(x)
                res%nfev = res%nfev + 1
                f_valid = .true.
            end if
            ok_eval = ieee_is_finite(fval)
        end subroutine ensure_fun

        subroutine ensure_hessian(ok_eval)
            logical, intent(out) :: ok_eval !! true when the cached Hessian is finite

            if (.not. hess_valid) then
                call hess(x, hmat)
                res%nhev = res%nhev + 1
                hess_valid = .true.
            end if
            ok_eval = all(ieee_is_finite(hmat))
        end subroutine ensure_hessian

        subroutine model_value(step, value, ok_eval)
            real(dp), intent(in) :: step(:) !! trial trust-region step
            real(dp), intent(out) :: value !! quadratic model at the trial step
            logical, intent(out) :: ok_eval !! true when the model is finite

            logical :: h_ok

            call ensure_fun(ok_eval)
            if (.not. ok_eval) return
            call ensure_hessian(h_ok)
            if (.not. h_ok) then
                ok_eval = .false.
                return
            end if
            value = fval + dot_product(g, step) &
                + 0.5_dp * dot_product(step, matmul(hmat, step))
            ok_eval = ieee_is_finite(value)
        end subroutine model_value

        subroutine solve_subproblem(step, boundary, ok_solve)
            real(dp), intent(out) :: step(:) !! proposed dog-leg trust-region step
            logical, intent(out) :: boundary !! true when step lies on trust boundary
            logical, intent(out) :: ok_solve !! true when the subproblem solve succeeds

            real(dp) :: bg(n)
            real(dp) :: denominator
            real(dp) :: d(n)
            integer :: info
            real(dp) :: l(n, n)
            real(dp) :: p_best(n)
            real(dp) :: p_u(n)
            real(dp) :: p_u_norm
            real(dp) :: scale
            real(dp) :: ta
            real(dp) :: tb
            real(dp) :: y(n)
            integer :: i
            logical :: h_ok

            call ensure_hessian(h_ok)
            if (.not. h_ok) then
                ok_solve = .false.
                step = 0.0_dp
                boundary = .false.
                return
            end if

            call cholesky_lower(hmat, l, info)
            if (info /= linalg_status_success) then
                ok_solve = .false.
                step = 0.0_dp
                boundary = .false.
                return
            end if
            call solve_lower(l, -g, y, info)
            if (info /= linalg_status_success) then
                ok_solve = .false.
                step = 0.0_dp
                boundary = .false.
                return
            end if
            p_best = 0.0_dp
            do i = n, 1, -1
                denominator = l(i, i)
                if (.not. ieee_is_finite(denominator) .or. abs(denominator) <= tiny(1.0_dp)) then
                    ok_solve = .false.
                    step = 0.0_dp
                    boundary = .false.
                    return
                end if
                p_best(i) = y(i)
                if (i < n) p_best(i) = p_best(i) &
                    - dot_product(l(i + 1:n, i), p_best(i + 1:n))
                p_best(i) = p_best(i) / denominator
            end do
            if (.not. all(ieee_is_finite(p_best))) then
                ok_solve = .false.
                step = 0.0_dp
                boundary = .false.
                return
            end if

            if (norm2(p_best) < trust_radius) then
                step = p_best
                boundary = .false.
                ok_solve = .true.
                return
            end if

            bg = matmul(hmat, g)
            denominator = dot_product(g, bg)
            if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_dp) then
                ok_solve = .false.
                step = 0.0_dp
                boundary = .false.
                return
            end if
            p_u = -(dot_product(g, g) / denominator) * g
            p_u_norm = norm2(p_u)
            if (.not. ieee_is_finite(p_u_norm) .or. p_u_norm <= 0.0_dp) then
                ok_solve = .false.
                step = 0.0_dp
                boundary = .false.
                return
            end if
            if (p_u_norm >= trust_radius) then
                scale = trust_radius / p_u_norm
                step = scale * p_u
                boundary = .true.
                ok_solve = .true.
                return
            end if

            d = p_best - p_u
            call boundary_intersections(p_u, d, trust_radius, ta, tb, ok_solve)
            if (.not. ok_solve) then
                step = 0.0_dp
                boundary = .false.
                return
            end if
            step = p_u + tb * d
            boundary = .true.
        end subroutine solve_subproblem

        subroutine boundary_intersections(z, d, radius, ta, tb, ok_intersection)
            real(dp), intent(in) :: z(:) !! line origin
            real(dp), intent(in) :: d(:) !! line direction
            real(dp), intent(in) :: radius !! trust-region radius
            real(dp), intent(out) :: ta !! smaller line-sphere parameter
            real(dp), intent(out) :: tb !! larger line-sphere parameter
            logical, intent(out) :: ok_intersection !! true when real intersections exist

            real(dp) :: a
            real(dp) :: aux
            real(dp) :: b
            real(dp) :: c
            real(dp) :: discriminant
            real(dp) :: root_disc
            real(dp) :: tmp

            a = dot_product(d, d)
            b = 2.0_dp * dot_product(z, d)
            c = dot_product(z, z) - radius**2
            discriminant = b * b - 4.0_dp * a * c
            if (a <= 0.0_dp .or. discriminant < 0.0_dp .or. &
                    .not. ieee_is_finite(discriminant)) then
                ok_intersection = .false.
                ta = 0.0_dp
                tb = 0.0_dp
                return
            end if

            root_disc = sqrt(max(0.0_dp, discriminant))
            aux = b + sign(root_disc, b)
            if (abs(aux) > tiny(1.0_dp)) then
                ta = -aux / (2.0_dp * a)
                tb = -2.0_dp * c / aux
            else
                ta = (-b - root_disc) / (2.0_dp * a)
                tb = (-b + root_disc) / (2.0_dp * a)
            end if
            if (ta > tb) then
                tmp = ta
                ta = tb
                tb = tmp
            end if
            ok_intersection = ieee_is_finite(ta) .and. ieee_is_finite(tb)
        end subroutine boundary_intersections

    end function dogleg

    function trust_ncg(fun, x0, jac, hess, hessp, initial_trust_radius, &
            max_trust_radius, eta, gtol, maxiter) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function) :: jac !! analytic gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic Hessian
        procedure(multivariate_hessian_product_function), optional :: hessp !! Hessian product
        real(dp), intent(in), optional :: initial_trust_radius !! initial trust-region radius
        real(dp), intent(in), optional :: max_trust_radius !! maximum trust-region radius
        real(dp), intent(in), optional :: eta !! step acceptance threshold in [0, 0.25)
        real(dp), intent(in), optional :: gtol !! Euclidean gradient-norm tolerance
        integer, intent(in), optional :: maxiter !! maximum trust-region iterations
        type(multivariate_minimize_result) :: res

        real(dp) :: actual_reduction
        real(dp) :: eta_use
        real(dp) :: f_proposed
        real(dp) :: fval
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: g_proposed(:)
        real(dp) :: gtol_use
        real(dp), allocatable :: hmat(:, :)
        real(dp) :: initial_radius_use
        real(dp) :: max_radius_use
        integer :: maxiter_use
        integer :: n
        real(dp), allocatable :: p(:)
        real(dp) :: predicted_reduction
        real(dp) :: predicted_value
        real(dp) :: rho
        real(dp) :: trust_radius
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_proposed(:)
        logical :: f_valid
        logical :: hess_valid
        logical :: hits_boundary
        logical :: ok

        n = size(x0)
        call init_result(res, 'trust-ncg', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        if (.not. present(hess) .and. .not. present(hessp)) return
        ! SciPy 1.17 ScalarFunction contributes one Hessian evaluation to the
        ! reported count when trust-ncg is driven by a Hessian product callback.
        if (present(hessp)) res%nhev = 1

        initial_radius_use = 1.0_dp
        if (present(initial_trust_radius)) initial_radius_use = initial_trust_radius
        max_radius_use = 1000.0_dp
        if (present(max_trust_radius)) max_radius_use = max_trust_radius
        eta_use = 0.15_dp
        if (present(eta)) eta_use = eta
        gtol_use = 1.0e-4_dp
        if (present(gtol)) gtol_use = gtol
        maxiter_use = 200 * n
        if (present(maxiter)) maxiter_use = maxiter

        if (.not. ieee_is_finite(initial_radius_use) .or. initial_radius_use <= 0.0_dp) return
        if (.not. ieee_is_finite(max_radius_use) .or. max_radius_use <= 0.0_dp) return
        if (initial_radius_use >= max_radius_use) return
        if (.not. ieee_is_finite(eta_use) .or. eta_use < 0.0_dp .or. eta_use >= 0.25_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return
        if (maxiter_use < 0) return

        allocate(x(n), x_proposed(n), g(n), g_proposed(n), p(n), hmat(n, n))
        x = x0
        trust_radius = initial_radius_use
        f_valid = .false.
        hess_valid = .false.

        call jac(x, g)
        res%njev = res%njev + 1
        if (.not. all(ieee_is_finite(g))) then
            res%x = x
            res%jac = g
            call finish_result(res, 3, .false., 'A linalg error occurred, such as a non-psd Hessian.')
            return
        end if

        do while (norm2(g) >= gtol_use)
            call solve_subproblem(p, hits_boundary, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            call model_value(p, predicted_value, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            x_proposed = x + p
            f_proposed = fun(x_proposed)
            res%nfev = res%nfev + 1
            if (.not. ieee_is_finite(f_proposed)) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            actual_reduction = fval - f_proposed
            predicted_reduction = fval - predicted_value
            if (predicted_reduction <= 0.0_dp) then
                call finish_result(res, 2, .false., &
                    'A bad approximation caused failure to predict improvement.')
                exit
            end if
            rho = actual_reduction / predicted_reduction

            if (rho < 0.25_dp) then
                trust_radius = 0.25_dp * trust_radius
            else if (rho > 0.75_dp .and. hits_boundary) then
                trust_radius = min(2.0_dp * trust_radius, max_radius_use)
            end if

            if (rho > eta_use) then
                x = x_proposed
                fval = f_proposed
                f_valid = .true.
                hess_valid = .false.
                call jac(x, g_proposed)
                res%njev = res%njev + 1
                if (.not. all(ieee_is_finite(g_proposed))) then
                    call finish_result(res, 3, .false., &
                        'A linalg error occurred, such as a non-psd Hessian.')
                    exit
                end if
                g = g_proposed
            end if

            res%nit = res%nit + 1
            if (norm2(g) < gtol_use) then
                call finish_result(res, 0, .true., 'Optimization terminated successfully.')
                exit
            end if
            if (res%nit >= maxiter_use) then
                call finish_result(res, 1, .false., &
                    'Maximum number of iterations has been exceeded.')
                exit
            end if
        end do

        if (.not. f_valid) then
            fval = fun(x)
            res%nfev = res%nfev + 1
            f_valid = .true.
        end if
        res%x = x
        res%fun = fval
        res%jac = g
        if (res%status < 0) then
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine ensure_fun(ok_eval)
            logical, intent(out) :: ok_eval !! true when the cached objective is finite

            if (.not. f_valid) then
                fval = fun(x)
                res%nfev = res%nfev + 1
                f_valid = .true.
            end if
            ok_eval = ieee_is_finite(fval)
        end subroutine ensure_fun

        subroutine apply_hessian(v, hv, ok_eval)
            real(dp), intent(in) :: v(:) !! Hessian-product direction
            real(dp), intent(out) :: hv(:) !! Hessian-vector product
            logical, intent(out) :: ok_eval !! true when the product is finite

            if (present(hessp)) then
                call invoke_hessian_product(hessp, x, v, hv)
                res%nhev = res%nhev + 1
            else
                if (.not. hess_valid) then
                    call invoke_hessian(hess, x, hmat)
                    res%nhev = res%nhev + 1
                    hess_valid = .true.
                end if
                hv = matmul(hmat, v)
            end if
            ok_eval = all(ieee_is_finite(hv))
        end subroutine apply_hessian

        subroutine model_value(step, value, ok_eval)
            real(dp), intent(in) :: step(:) !! trial trust-region step
            real(dp), intent(out) :: value !! quadratic model at the trial step
            logical, intent(out) :: ok_eval !! true when the model is finite

            real(dp) :: hp(n)
            logical :: h_ok

            call ensure_fun(ok_eval)
            if (.not. ok_eval) return
            call apply_hessian(step, hp, h_ok)
            if (.not. h_ok) then
                ok_eval = .false.
                return
            end if
            value = fval + dot_product(g, step) + 0.5_dp * dot_product(step, hp)
            ok_eval = ieee_is_finite(value)
        end subroutine model_value

        subroutine solve_subproblem(step, boundary, ok_solve)
            real(dp), intent(out) :: step(:) !! proposed trust-region step
            logical, intent(out) :: boundary !! true when step lies on trust boundary
            logical, intent(out) :: ok_solve !! true when the subproblem solve succeeds

            real(dp) :: alpha
            real(dp) :: beta
            real(dp) :: bd(n)
            real(dp) :: d(n)
            real(dp) :: dbd
            real(dp) :: gmag
            logical :: h_ok
            real(dp) :: pa(n)
            real(dp) :: pb(n)
            real(dp) :: qa
            real(dp) :: qb
            real(dp) :: r(n)
            real(dp) :: r_next(n)
            real(dp) :: r_next_squared
            real(dp) :: r_squared
            real(dp) :: ta
            real(dp) :: tb
            real(dp) :: tolerance
            real(dp) :: z(n)
            real(dp) :: z_next(n)

            gmag = norm2(g)
            tolerance = min(0.5_dp, sqrt(gmag)) * gmag
            z = 0.0_dp
            step = 0.0_dp
            boundary = .false.
            ok_solve = .true.
            if (gmag < tolerance) return

            r = g
            d = -r
            do
                call apply_hessian(d, bd, h_ok)
                if (.not. h_ok) then
                    ok_solve = .false.
                    return
                end if
                dbd = dot_product(d, bd)
                if (dbd <= 0.0_dp) then
                    call boundary_intersections(z, d, trust_radius, ta, tb, ok_solve)
                    if (.not. ok_solve) return
                    pa = z + ta * d
                    pb = z + tb * d
                    call model_value(pa, qa, ok_solve)
                    if (.not. ok_solve) return
                    call model_value(pb, qb, ok_solve)
                    if (.not. ok_solve) return
                    if (qa < qb) then
                        step = pa
                    else
                        step = pb
                    end if
                    boundary = .true.
                    return
                end if

                r_squared = dot_product(r, r)
                alpha = r_squared / dbd
                z_next = z + alpha * d
                if (norm2(z_next) >= trust_radius) then
                    call boundary_intersections(z, d, trust_radius, ta, tb, ok_solve)
                    if (.not. ok_solve) return
                    step = z + tb * d
                    boundary = .true.
                    return
                end if

                r_next = r + alpha * bd
                r_next_squared = dot_product(r_next, r_next)
                if (sqrt(r_next_squared) < tolerance) then
                    step = z_next
                    return
                end if

                beta = r_next_squared / r_squared
                d = -r_next + beta * d
                z = z_next
                r = r_next
            end do
        end subroutine solve_subproblem

        subroutine boundary_intersections(z, d, radius, ta, tb, ok_intersection)
            real(dp), intent(in) :: z(:) !! line origin
            real(dp), intent(in) :: d(:) !! line direction
            real(dp), intent(in) :: radius !! trust-region radius
            real(dp), intent(out) :: ta !! smaller line-sphere parameter
            real(dp), intent(out) :: tb !! larger line-sphere parameter
            logical, intent(out) :: ok_intersection !! true when real intersections exist

            real(dp) :: a
            real(dp) :: aux
            real(dp) :: b
            real(dp) :: c
            real(dp) :: discriminant
            real(dp) :: root_disc
            real(dp) :: tmp

            a = dot_product(d, d)
            b = 2.0_dp * dot_product(z, d)
            c = dot_product(z, z) - radius**2
            discriminant = b * b - 4.0_dp * a * c
            if (a <= 0.0_dp .or. discriminant < 0.0_dp .or. &
                    .not. ieee_is_finite(discriminant)) then
                ok_intersection = .false.
                ta = 0.0_dp
                tb = 0.0_dp
                return
            end if

            root_disc = sqrt(max(0.0_dp, discriminant))
            aux = b + sign(root_disc, b)
            if (abs(aux) > tiny(1.0_dp)) then
                ta = -aux / (2.0_dp * a)
                tb = -2.0_dp * c / aux
            else
                ta = (-b - root_disc) / (2.0_dp * a)
                tb = (-b + root_disc) / (2.0_dp * a)
            end if
            if (ta > tb) then
                tmp = ta
                ta = tb
                tb = tmp
            end if
            ok_intersection = ieee_is_finite(ta) .and. ieee_is_finite(tb)
        end subroutine boundary_intersections

    end function trust_ncg

    function trust_krylov(fun, x0, jac, hess, hessp, initial_trust_radius, &
            max_trust_radius, eta, gtol, maxiter, inexact) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function) :: jac !! analytic gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic Hessian
        procedure(multivariate_hessian_product_function), optional :: hessp !! Hessian product
        real(dp), intent(in), optional :: initial_trust_radius !! initial trust-region radius
        real(dp), intent(in), optional :: max_trust_radius !! maximum trust-region radius
        real(dp), intent(in), optional :: eta !! step acceptance threshold in [0, 0.25)
        real(dp), intent(in), optional :: gtol !! Euclidean gradient-norm tolerance
        integer, intent(in), optional :: maxiter !! maximum trust-region iterations
        logical, intent(in), optional :: inexact !! use SciPy's default inexact GLTR forcing sequence
        type(multivariate_minimize_result) :: res

        real(dp) :: actual_reduction
        logical :: basis_exhausted
        real(dp) :: eta_use
        real(dp) :: f_proposed
        real(dp) :: fval
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: g_proposed(:)
        real(dp) :: gtol_use
        real(dp), allocatable :: hmat(:, :)
        real(dp), allocatable :: hq(:, :)
        real(dp) :: initial_radius_use
        logical :: inexact_use
        integer :: k_built
        real(dp) :: max_radius_use
        integer :: maxiter_use
        integer :: n
        real(dp), allocatable :: p(:)
        real(dp) :: predicted_reduction
        real(dp) :: predicted_value
        real(dp), allocatable :: q(:, :)
        real(dp) :: rho
        real(dp), allocatable :: tmat(:, :)
        real(dp) :: trust_radius
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_proposed(:)
        logical :: hess_valid
        logical :: hits_boundary
        logical :: ok

        n = size(x0)
        call init_result(res, 'trust-krylov', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        if (.not. present(hess) .and. .not. present(hessp)) return
        ! SciPy's ScalarFunction/TRLib wrapper contributes one bookkeeping
        ! Hessian evaluation when trust-krylov is driven by hessp.
        if (present(hessp)) res%nhev = 1

        initial_radius_use = 1.0_dp
        if (present(initial_trust_radius)) initial_radius_use = initial_trust_radius
        max_radius_use = 1000.0_dp
        if (present(max_trust_radius)) max_radius_use = max_trust_radius
        eta_use = 0.15_dp
        if (present(eta)) eta_use = eta
        gtol_use = 1.0e-4_dp
        if (present(gtol)) gtol_use = gtol
        maxiter_use = 200 * n
        if (present(maxiter)) maxiter_use = maxiter
        inexact_use = .true.
        if (present(inexact)) inexact_use = inexact

        if (.not. ieee_is_finite(initial_radius_use) .or. initial_radius_use <= 0.0_dp) return
        if (.not. ieee_is_finite(max_radius_use) .or. max_radius_use <= 0.0_dp) return
        if (initial_radius_use >= max_radius_use) return
        if (.not. ieee_is_finite(eta_use) .or. eta_use < 0.0_dp .or. eta_use >= 0.25_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return
        if (maxiter_use < 0) return

        allocate(x(n), x_proposed(n), g(n), g_proposed(n), p(n), hmat(n, n))
        allocate(q(n, n), hq(n, n), tmat(n, n))
        x = x0
        trust_radius = initial_radius_use
        hess_valid = .false.

        fval = fun(x)
        res%nfev = res%nfev + 1
        call jac(x, g)
        res%njev = res%njev + 1
        if (.not. ieee_is_finite(fval) .or. .not. all(ieee_is_finite(g))) then
            res%x = x
            res%fun = fval
            res%jac = g
            call finish_result(res, 3, .false., &
                'A linalg error occurred, such as a non-psd Hessian.')
            return
        end if
        call reset_krylov_basis()

        do while (norm2(g) >= gtol_use)
            call solve_subproblem(p, hits_boundary, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            call model_value(p, predicted_value, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            x_proposed = x + p
            f_proposed = fun(x_proposed)
            res%nfev = res%nfev + 1
            call jac(x_proposed, g_proposed)
            res%njev = res%njev + 1
            if (.not. ieee_is_finite(f_proposed) .or. .not. all(ieee_is_finite(g_proposed))) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            actual_reduction = fval - f_proposed
            predicted_reduction = fval - predicted_value
            if (predicted_reduction <= 0.0_dp) then
                call finish_result(res, 2, .false., &
                    'A bad approximation caused failure to predict improvement.')
                exit
            end if
            rho = actual_reduction / predicted_reduction

            if (rho < 0.25_dp) then
                trust_radius = 0.25_dp * trust_radius
            else if (rho > 0.75_dp .and. hits_boundary) then
                trust_radius = min(2.0_dp * trust_radius, max_radius_use)
            end if

            if (rho > eta_use) then
                x = x_proposed
                fval = f_proposed
                g = g_proposed
                hess_valid = .false.
                call reset_krylov_basis()
            end if

            res%nit = res%nit + 1
            if (norm2(g) < gtol_use) then
                call finish_result(res, 0, .true., 'Optimization terminated successfully.')
                exit
            end if
            if (res%nit >= maxiter_use) then
                call finish_result(res, 1, .false., &
                    'Maximum number of iterations has been exceeded.')
                exit
            end if
        end do

        res%x = x
        res%fun = fval
        res%jac = g
        if (res%status < 0) then
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine reset_krylov_basis()
            real(dp) :: gmag

            q = 0.0_dp
            hq = 0.0_dp
            tmat = 0.0_dp
            k_built = 0
            basis_exhausted = .false.
            gmag = norm2(g)
            if (gmag > 0.0_dp) q(:, 1) = g / gmag
        end subroutine reset_krylov_basis

        subroutine apply_hessian(v, hv, ok_eval)
            real(dp), intent(in) :: v(:) !! Hessian-product direction
            real(dp), intent(out) :: hv(:) !! Hessian-vector product
            logical, intent(out) :: ok_eval !! true when the product is finite

            if (present(hessp)) then
                call invoke_hessian_product(hessp, x, v, hv)
                res%nhev = res%nhev + 1
            else
                if (.not. hess_valid) then
                    call invoke_hessian(hess, x, hmat)
                    res%nhev = res%nhev + 1
                    hess_valid = .true.
                end if
                hv = matmul(hmat, v)
            end if
            ok_eval = all(ieee_is_finite(hv))
        end subroutine apply_hessian

        subroutine model_value(step, value, ok_eval)
            real(dp), intent(in) :: step(:) !! trial trust-region step
            real(dp), intent(out) :: value !! quadratic model at the trial step
            logical, intent(out) :: ok_eval !! true when the model is finite

            real(dp) :: hp(n)
            logical :: h_ok

            call apply_hessian(step, hp, h_ok)
            if (.not. h_ok) then
                ok_eval = .false.
                return
            end if
            value = fval + dot_product(g, step) + 0.5_dp * dot_product(step, hp)
            ok_eval = ieee_is_finite(value)
        end subroutine model_value

        subroutine solve_subproblem(step, boundary, ok_solve)
            real(dp), intent(out) :: step(:) !! proposed trust-region step
            logical, intent(out) :: boundary !! true when step lies on trust boundary
            logical, intent(out) :: ok_solve !! true when the subproblem solve succeeds

            integer :: i
            integer :: j
            real(dp) :: gmag
            logical :: h_ok
            real(dp) :: lambda_shift
            real(dp) :: residual(n)
            real(dp) :: residual_norm
            real(dp) :: stop_b
            real(dp) :: stop_i
            real(dp) :: threshold
            real(dp) :: w(n)
            real(dp), allocatable :: y(:)

            gmag = norm2(g)
            step = 0.0_dp
            boundary = .false.
            ok_solve = .true.
            if (gmag <= tiny(1.0_dp)) return

            if (inexact_use) then
                stop_i = min(0.5_dp, gmag) * gmag
                stop_b = max(1.0e-6_dp, min(0.5_dp, sqrt(gmag))) * gmag
            else
                stop_i = 1.0e-8_dp * gmag
                stop_b = 1.0e-6_dp * gmag
            end if

            do
                if (k_built > 0) then
                    allocate(y(k_built))
                    call solve_reduced_trust_region(k_built, gmag, trust_radius, &
                        y, lambda_shift, boundary, ok_solve)
                    if (.not. ok_solve) return
                    step = matmul(q(:, 1:k_built), y)
                    residual = g + matmul(hq(:, 1:k_built), y) + lambda_shift * step
                    residual_norm = norm2(residual)
                    deallocate(y)
                    if (residual_norm <= merge(stop_b, stop_i, boundary)) return
                    if (k_built >= n .or. basis_exhausted) return
                end if

                j = k_built + 1
                call apply_hessian(q(:, j), hq(:, j), h_ok)
                if (.not. h_ok) then
                    ok_solve = .false.
                    return
                end if
                do i = 1, j
                    tmat(i, j) = dot_product(q(:, i), hq(:, j))
                    tmat(j, i) = tmat(i, j)
                end do
                k_built = j

                allocate(y(k_built))
                call solve_reduced_trust_region(k_built, gmag, trust_radius, &
                    y, lambda_shift, boundary, ok_solve)
                if (.not. ok_solve) return
                step = matmul(q(:, 1:k_built), y)
                residual = g + matmul(hq(:, 1:k_built), y) + lambda_shift * step
                residual_norm = norm2(residual)
                deallocate(y)

                if (k_built < n) then
                    w = hq(:, k_built)
                    w = w - matmul(q(:, 1:k_built), &
                        matmul(transpose(q(:, 1:k_built)), w))
                    w = w - matmul(q(:, 1:k_built), &
                        matmul(transpose(q(:, 1:k_built)), w))
                    threshold = 64.0_dp * epsilon(1.0_dp) * &
                        max(1.0_dp, norm2(hq(:, k_built)))
                    if (norm2(w) > threshold) then
                        q(:, k_built + 1) = w / norm2(w)
                    else
                        basis_exhausted = .true.
                    end if
                end if

                if (residual_norm <= merge(stop_b, stop_i, boundary)) return
                if (k_built >= n .or. basis_exhausted) return
            end do
        end subroutine solve_subproblem

        subroutine solve_reduced_trust_region(k, gmag, radius, y, lambda_shift, &
                boundary, ok_solve)
            integer, intent(in) :: k !! current Krylov dimension
            real(dp), intent(in) :: gmag !! norm of the full-space gradient
            real(dp), intent(in) :: radius !! trust-region radius
            real(dp), intent(out) :: y(:) !! reduced-space trust-region step
            real(dp), intent(out) :: lambda_shift !! trust-region multiplier
            logical, intent(out) :: boundary !! true when the trust boundary is active
            logical, intent(out) :: ok_solve !! true on a finite reduced solve

            real(dp) :: cbar(k)
            real(dp) :: denom(k)
            real(dp) :: denom_tol
            real(dp) :: eig_scale
            real(dp) :: eig_tol
            real(dp) :: evals(k)
            real(dp) :: evecs(k, k)
            real(dp) :: hi
            integer :: i
            integer :: idx_zero
            integer :: iter
            real(dp) :: lin_tol
            real(dp) :: lo
            real(dp) :: lower
            real(dp) :: mid
            real(dp) :: nmid
            real(dp) :: nhi
            real(dp) :: norm_particular
            integer :: status_linalg
            real(dp) :: tred(k, k)
            logical :: zero_compatible
            real(dp) :: ybar(k)

            y = 0.0_dp
            lambda_shift = 0.0_dp
            boundary = .false.
            ok_solve = .false.
            tred = tmat(1:k, 1:k)
            call symmetric_eigen_jacobi(tred, evals, evecs, status_linalg)
            if (status_linalg /= linalg_status_success) return
            if (.not. all(ieee_is_finite(evals)) .or. &
                    .not. all(ieee_is_finite(evecs))) return

            cbar = gmag * evecs(1, :)
            eig_scale = max(1.0_dp, maxval(abs(evals)))
            eig_tol = 64.0_dp * epsilon(1.0_dp) * eig_scale

            if (minval(evals) > eig_tol) then
                ybar = -cbar / evals
                if (norm2(ybar) <= radius) then
                    y = matmul(evecs, ybar)
                    ok_solve = all(ieee_is_finite(y))
                    return
                end if
            end if

            lower = max(0.0_dp, -evals(1))
            denom = evals + lower
            denom_tol = 256.0_dp * epsilon(1.0_dp) * &
                max(1.0_dp, eig_scale, abs(lower))
            lin_tol = 256.0_dp * epsilon(1.0_dp) * max(1.0_dp, gmag)
            zero_compatible = .true.
            idx_zero = 0
            do i = 1, k
                if (abs(denom(i)) <= denom_tol) then
                    if (idx_zero == 0) idx_zero = i
                    if (abs(cbar(i)) > lin_tol) zero_compatible = .false.
                end if
            end do

            if (zero_compatible .and. idx_zero > 0) then
                ybar = 0.0_dp
                do i = 1, k
                    if (abs(denom(i)) > denom_tol) ybar(i) = -cbar(i) / denom(i)
                end do
                norm_particular = norm2(ybar)
                if (norm_particular <= radius * (1.0_dp + 64.0_dp * epsilon(1.0_dp))) then
                    if (lower <= eig_tol) then
                        y = matmul(evecs, ybar)
                        ok_solve = all(ieee_is_finite(y))
                        return
                    end if
                    ybar(idx_zero) = sqrt(max(0.0_dp, radius * radius - &
                        norm_particular * norm_particular))
                    y = matmul(evecs, ybar)
                    lambda_shift = lower
                    boundary = .true.
                    ok_solve = all(ieee_is_finite(y))
                    return
                end if
            end if

            lo = lower + max(1.0e-14_dp, epsilon(1.0_dp) * &
                max(1.0_dp, eig_scale, abs(lower)))
            hi = max(1.0_dp, 2.0_dp * lo, gmag / radius - evals(1))
            do iter = 1, 100
                nhi = shifted_solution_norm(hi, evals, cbar, k)
                if (ieee_is_finite(nhi) .and. nhi <= radius) exit
                hi = max(2.0_dp * hi, hi + 1.0_dp)
            end do
            if (iter > 100) return

            do iter = 1, 120
                mid = 0.5_dp * (lo + hi)
                nmid = shifted_solution_norm(mid, evals, cbar, k)
                if (.not. ieee_is_finite(nmid)) then
                    lo = mid
                else if (nmid > radius) then
                    lo = mid
                else
                    hi = mid
                end if
                if (hi - lo <= 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, hi)) exit
            end do

            lambda_shift = hi
            do i = 1, k
                denom(i) = evals(i) + lambda_shift
                if (abs(denom(i)) <= tiny(1.0_dp)) return
                ybar(i) = -cbar(i) / denom(i)
            end do
            y = matmul(evecs, ybar)
            boundary = .true.
            ok_solve = all(ieee_is_finite(y)) .and. ieee_is_finite(lambda_shift)
        end subroutine solve_reduced_trust_region

        function shifted_solution_norm(shift, evals, cbar, k) result(value)
            real(dp), intent(in) :: shift !! spectral shift
            real(dp), intent(in) :: evals(:) !! reduced Hessian eigenvalues
            real(dp), intent(in) :: cbar(:) !! reduced linear term in eigenbasis
            integer, intent(in) :: k !! active reduced dimension
            real(dp) :: value

            real(dp) :: denom
            integer :: i
            real(dp) :: sumsq

            sumsq = 0.0_dp
            do i = 1, k
                denom = evals(i) + shift
                if (abs(denom) <= tiny(1.0_dp)) then
                    value = huge(1.0_dp)
                    return
                end if
                sumsq = sumsq + (cbar(i) / denom)**2
            end do
            value = sqrt(sumsq)
        end function shifted_solution_norm

    end function trust_krylov

    function trust_exact(fun, x0, jac, hess, initial_trust_radius, max_trust_radius, &
            eta, gtol, maxiter, subproblem_maxiter) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function) :: jac !! analytic gradient
        procedure(multivariate_hessian_function) :: hess !! analytic Hessian
        real(dp), intent(in), optional :: initial_trust_radius !! initial trust-region radius
        real(dp), intent(in), optional :: max_trust_radius !! maximum trust-region radius
        real(dp), intent(in), optional :: eta !! step acceptance threshold in [0, 0.25)
        real(dp), intent(in), optional :: gtol !! Euclidean gradient-norm tolerance
        integer, intent(in), optional :: maxiter !! maximum trust-region iterations
        integer, intent(in), optional :: subproblem_maxiter !! maximum damping iterations per subproblem
        type(multivariate_minimize_result) :: res

        real(dp) :: actual_reduction
        real(dp) :: eta_use
        real(dp) :: f_proposed
        real(dp) :: fval
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: g_proposed(:)
        real(dp) :: gtol_use
        real(dp), allocatable :: hmat(:, :)
        real(dp), allocatable :: h_proposed(:, :)
        real(dp) :: initial_radius_use
        real(dp) :: max_radius_use
        integer :: maxiter_use
        integer :: n
        real(dp), allocatable :: p(:)
        real(dp) :: predicted_reduction
        real(dp) :: predicted_value
        real(dp) :: previous_lambda_lb
        real(dp) :: previous_trust_radius
        real(dp) :: rho
        integer :: subproblem_maxiter_use
        real(dp) :: trust_radius
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_proposed(:)
        logical :: f_valid
        logical :: have_previous_lambda_lb
        logical :: hess_valid
        logical :: hits_boundary
        logical :: ok

        n = size(x0)
        call init_result(res, 'trust-exact', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        initial_radius_use = 1.0_dp
        if (present(initial_trust_radius)) initial_radius_use = initial_trust_radius
        max_radius_use = 1000.0_dp
        if (present(max_trust_radius)) max_radius_use = max_trust_radius
        eta_use = 0.15_dp
        if (present(eta)) eta_use = eta
        gtol_use = 1.0e-4_dp
        if (present(gtol)) gtol_use = gtol
        maxiter_use = 200 * n
        if (present(maxiter)) maxiter_use = maxiter
        subproblem_maxiter_use = 25
        if (present(subproblem_maxiter)) subproblem_maxiter_use = subproblem_maxiter

        if (.not. ieee_is_finite(initial_radius_use) .or. initial_radius_use <= 0.0_dp) return
        if (.not. ieee_is_finite(max_radius_use) .or. max_radius_use <= 0.0_dp) return
        if (initial_radius_use >= max_radius_use) return
        if (.not. ieee_is_finite(eta_use) .or. eta_use < 0.0_dp .or. eta_use >= 0.25_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return
        if (maxiter_use < 0 .or. subproblem_maxiter_use < 0) return

        allocate(x(n), x_proposed(n), g(n), g_proposed(n), p(n))
        allocate(hmat(n, n), h_proposed(n, n))
        x = x0
        trust_radius = initial_radius_use
        f_valid = .false.
        hess_valid = .false.
        previous_trust_radius = -1.0_dp
        previous_lambda_lb = 0.0_dp
        have_previous_lambda_lb = .false.

        call jac(x, g)
        res%njev = res%njev + 1
        if (.not. all(ieee_is_finite(g))) then
            res%x = x
            res%jac = g
            call finish_result(res, 3, .false., &
                'A linalg error occurred, such as a non-psd Hessian.')
            return
        end if

        do while (norm2(g) >= gtol_use)
            call solve_subproblem(p, hits_boundary, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            call model_value(p, predicted_value, ok)
            if (.not. ok) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            x_proposed = x + p
            f_proposed = fun(x_proposed)
            res%nfev = res%nfev + 1
            call hess(x_proposed, h_proposed)
            res%nhev = res%nhev + 1
            if (.not. ieee_is_finite(f_proposed) .or. &
                    .not. all(ieee_is_finite(h_proposed))) then
                call finish_result(res, 3, .false., &
                    'A linalg error occurred, such as a non-psd Hessian.')
                exit
            end if

            actual_reduction = fval - f_proposed
            predicted_reduction = fval - predicted_value
            if (predicted_reduction <= 0.0_dp) then
                call finish_result(res, 2, .false., &
                    'A bad approximation caused failure to predict improvement.')
                exit
            end if
            rho = actual_reduction / predicted_reduction

            if (rho < 0.25_dp) then
                trust_radius = 0.25_dp * trust_radius
            else if (rho > 0.75_dp .and. hits_boundary) then
                trust_radius = min(2.0_dp * trust_radius, max_radius_use)
            end if

            if (rho > eta_use) then
                x = x_proposed
                fval = f_proposed
                f_valid = .true.
                hmat = h_proposed
                hess_valid = .true.
                previous_trust_radius = -1.0_dp
                previous_lambda_lb = 0.0_dp
                have_previous_lambda_lb = .false.
                call jac(x, g_proposed)
                res%njev = res%njev + 1
                if (.not. all(ieee_is_finite(g_proposed))) then
                    call finish_result(res, 3, .false., &
                        'A linalg error occurred, such as a non-psd Hessian.')
                    exit
                end if
                g = g_proposed
            end if

            res%nit = res%nit + 1
            if (norm2(g) < gtol_use) then
                call finish_result(res, 0, .true., 'Optimization terminated successfully.')
                exit
            end if
            if (res%nit >= maxiter_use) then
                call finish_result(res, 1, .false., &
                    'Maximum number of iterations has been exceeded.')
                exit
            end if
        end do

        if (.not. f_valid) then
            fval = fun(x)
            res%nfev = res%nfev + 1
            f_valid = .true.
        end if
        if (.not. hess_valid) then
            call hess(x, hmat)
            res%nhev = res%nhev + 1
            hess_valid = .true.
        end if
        res%x = x
        res%fun = fval
        res%jac = g
        if (res%status < 0) then
            call finish_result(res, 0, .true., 'Optimization terminated successfully.')
        end if

    contains

        subroutine ensure_fun(ok_eval)
            logical, intent(out) :: ok_eval !! true when the cached objective is finite

            if (.not. f_valid) then
                fval = fun(x)
                res%nfev = res%nfev + 1
                f_valid = .true.
            end if
            ok_eval = ieee_is_finite(fval)
        end subroutine ensure_fun

        subroutine ensure_hessian(ok_eval)
            logical, intent(out) :: ok_eval !! true when the cached Hessian is finite

            if (.not. hess_valid) then
                call hess(x, hmat)
                res%nhev = res%nhev + 1
                hess_valid = .true.
            end if
            ok_eval = all(ieee_is_finite(hmat))
        end subroutine ensure_hessian

        subroutine model_value(step, value, ok_eval)
            real(dp), intent(in) :: step(:) !! trial trust-region step
            real(dp), intent(out) :: value !! quadratic model at the trial step
            logical, intent(out) :: ok_eval !! true when the model is finite

            logical :: h_ok

            call ensure_fun(ok_eval)
            if (.not. ok_eval) return
            call ensure_hessian(h_ok)
            if (.not. h_ok) then
                ok_eval = .false.
                return
            end if
            value = fval + dot_product(g, step) &
                + 0.5_dp * dot_product(step, matmul(hmat, step))
            ok_eval = ieee_is_finite(value)
        end subroutine model_value

        subroutine solve_subproblem(step, boundary, ok_solve)
            real(dp), intent(out) :: step(:) !! nearly exact trust-region step
            logical, intent(out) :: boundary !! true when step lies on trust boundary
            logical, intent(out) :: ok_solve !! true when the subproblem solve succeeds

            real(dp) :: a_shift(n, n)
            logical :: already_factorized
            real(dp) :: close_to_zero
            real(dp) :: delta
            real(dp) :: delta_lambda
            real(dp) :: diagonal_min
            real(dp) :: gersh_lb
            real(dp) :: gersh_ub
            real(dp) :: gmag
            real(dp) :: hess_fro
            real(dp) :: hess_inf
            logical :: h_ok
            integer :: i
            integer :: info
            integer :: inner_iter
            real(dp) :: lambda_current
            real(dp) :: lambda_lb
            real(dp) :: lambda_new
            real(dp) :: lambda_ub
            real(dp) :: p_norm
            logical :: p_valid
            real(dp) :: quadratic_term
            real(dp) :: relative_error
            real(dp) :: row_sum
            real(dp) :: s_min
            real(dp) :: step_len
            real(dp) :: ta
            real(dp) :: tb
            real(dp) :: u(n, n)
            real(dp) :: u_new(n, n)
            real(dp) :: v(n)
            real(dp) :: v_norm
            real(dp) :: w(n)
            real(dp) :: w_norm
            real(dp) :: z_min(n)

            call ensure_hessian(h_ok)
            if (.not. h_ok .or. subproblem_maxiter_use == 0) then
                step = 0.0_dp
                boundary = .false.
                ok_solve = .false.
                return
            end if

            gmag = norm2(g)
            diagonal_min = hmat(1, 1)
            gersh_lb = huge(1.0_dp)
            gersh_ub = -huge(1.0_dp)
            hess_inf = 0.0_dp
            do i = 1, n
                diagonal_min = min(diagonal_min, hmat(i, i))
                row_sum = sum(abs(hmat(i, :)))
                hess_inf = max(hess_inf, row_sum)
                gersh_lb = min(gersh_lb, hmat(i, i) + abs(hmat(i, i)) - row_sum)
                gersh_ub = max(gersh_ub, hmat(i, i) - abs(hmat(i, i)) + row_sum)
            end do
            hess_fro = sqrt(sum(hmat * hmat))
            close_to_zero = real(n, dp) * epsilon(1.0_dp) * hess_inf

            lambda_ub = max(0.0_dp, gmag / trust_radius &
                + min(-gersh_lb, hess_fro, hess_inf))
            lambda_lb = max(0.0_dp, -diagonal_min, gmag / trust_radius &
                - min(gersh_ub, hess_fro, hess_inf))
            if (trust_radius < previous_trust_radius .and. have_previous_lambda_lb) then
                lambda_lb = max(previous_lambda_lb, lambda_lb)
            end if
            if (lambda_lb <= 0.0_dp) then
                lambda_current = 0.0_dp
            else
                lambda_current = max(sqrt(lambda_lb * lambda_ub), &
                    lambda_lb + 0.01_dp * (lambda_ub - lambda_lb))
            end if

            boundary = .true.
            already_factorized = .false.
            p_valid = .false.
            do inner_iter = 1, subproblem_maxiter_use
                if (already_factorized) then
                    already_factorized = .false.
                else
                    a_shift = hmat
                    do i = 1, n
                        a_shift(i, i) = a_shift(i, i) + lambda_current
                    end do
                    call cholesky_upper_partial(a_shift, u, info)
                end if

                if (info == 0 .and. gmag > close_to_zero) then
                    call solve_cholesky_upper(u, -g, step, info)
                    if (info /= 0 .or. .not. all(ieee_is_finite(step))) then
                        ok_solve = .false.
                        return
                    end if
                    p_valid = .true.
                    p_norm = norm2(step)
                    if (p_norm <= trust_radius .and. lambda_current <= 0.0_dp) then
                        boundary = .false.
                        exit
                    end if

                    call solve_upper_transpose(u, step, w, info)
                    if (info /= 0) then
                        ok_solve = .false.
                        return
                    end if
                    w_norm = norm2(w)
                    if (.not. ieee_is_finite(w_norm) .or. w_norm <= tiny(1.0_dp)) then
                        ok_solve = .false.
                        return
                    end if
                    delta_lambda = (p_norm / w_norm)**2 &
                        * (p_norm - trust_radius) / trust_radius
                    lambda_new = lambda_current + delta_lambda

                    if (p_norm < trust_radius) then
                        call estimate_smallest_singular_value(u, s_min, z_min, info)
                        if (info /= 0) then
                            ok_solve = .false.
                            return
                        end if
                        call boundary_intersections(step, z_min, trust_radius, ta, tb, ok_solve)
                        if (.not. ok_solve) return
                        if (abs(ta) <= abs(tb)) then
                            step_len = ta
                        else
                            step_len = tb
                        end if
                        quadratic_term = dot_product(step, matmul(a_shift, step))
                        if (quadratic_term + lambda_current * trust_radius**2 <= 0.0_dp) then
                            ok_solve = .false.
                            return
                        end if
                        relative_error = step_len**2 * s_min**2 &
                            / (quadratic_term + lambda_current * trust_radius**2)
                        if (relative_error <= 0.2_dp) then
                            step = step + step_len * z_min
                            exit
                        end if

                        lambda_ub = lambda_current
                        lambda_lb = max(lambda_lb, lambda_current - s_min**2)
                        a_shift = hmat
                        do i = 1, n
                            a_shift(i, i) = a_shift(i, i) + lambda_new
                        end do
                        call cholesky_upper_partial(a_shift, u_new, info)
                        if (info == 0) then
                            lambda_current = lambda_new
                            u = u_new
                            already_factorized = .true.
                        else
                            lambda_lb = max(lambda_lb, lambda_new)
                            lambda_current = damping_update(lambda_lb, lambda_ub)
                        end if
                    else
                        relative_error = abs(p_norm - trust_radius) / trust_radius
                        if (relative_error <= 0.1_dp) exit
                        lambda_lb = lambda_current
                        lambda_current = lambda_new
                    end if
                else if (info == 0 .and. gmag <= close_to_zero) then
                    if (lambda_current <= 0.0_dp) then
                        step = 0.0_dp
                        p_valid = .true.
                        boundary = .false.
                        exit
                    end if
                    call estimate_smallest_singular_value(u, s_min, z_min, info)
                    if (info /= 0) then
                        ok_solve = .false.
                        return
                    end if
                    step_len = trust_radius
                    step = step_len * z_min
                    p_valid = .true.
                    if (step_len**2 * s_min**2 &
                            <= 0.2_dp * lambda_current * trust_radius**2) exit
                    lambda_ub = lambda_current
                    lambda_lb = max(lambda_lb, lambda_current - s_min**2)
                    lambda_current = damping_update(lambda_lb, lambda_ub)
                else
                    call singular_leading_submatrix(a_shift, u, info, delta, v, ok_solve)
                    if (.not. ok_solve) return
                    v_norm = norm2(v)
                    if (.not. ieee_is_finite(v_norm) .or. v_norm <= tiny(1.0_dp)) then
                        ok_solve = .false.
                        return
                    end if
                    lambda_lb = max(lambda_lb, lambda_current + delta / v_norm**2)
                    lambda_current = damping_update(lambda_lb, lambda_ub)
                end if
            end do

            previous_lambda_lb = lambda_lb
            previous_trust_radius = trust_radius
            have_previous_lambda_lb = .true.
            ok_solve = p_valid .and. all(ieee_is_finite(step))
        end subroutine solve_subproblem

        subroutine cholesky_upper_partial(a, u, info)
            real(dp), intent(in) :: a(:, :) !! symmetric matrix; upper triangle is read
            real(dp), intent(out) :: u(:, :) !! partial upper factor
            integer, intent(out) :: info !! zero on success, first failing pivot otherwise

            real(dp) :: diagonal
            integer :: i
            integer :: j
            integer :: m
            real(dp) :: work

            m = size(a, 1)
            u = 0.0_dp
            do j = 1, m
                diagonal = a(j, j)
                if (j > 1) diagonal = diagonal - dot_product(u(1:j - 1, j), u(1:j - 1, j))
                if (.not. ieee_is_finite(diagonal) .or. diagonal <= 0.0_dp) then
                    info = j
                    return
                end if
                u(j, j) = sqrt(diagonal)
                do i = j + 1, m
                    work = a(j, i)
                    if (j > 1) work = work - dot_product(u(1:j - 1, j), u(1:j - 1, i))
                    u(j, i) = work / u(j, j)
                end do
            end do
            info = 0
        end subroutine cholesky_upper_partial

        subroutine solve_cholesky_upper(u, b, solution, info)
            real(dp), intent(in) :: u(:, :) !! upper factor with A = U^T U
            real(dp), intent(in) :: b(:) !! right-hand side
            real(dp), intent(out) :: solution(:) !! solution of A x = b
            integer, intent(out) :: info !! zero on success

            real(dp) :: y(size(b))

            call solve_upper_transpose(u, b, y, info)
            if (info /= 0) return
            call solve_upper(u, y, solution, info)
        end subroutine solve_cholesky_upper

        subroutine solve_upper(u, b, solution, info)
            real(dp), intent(in) :: u(:, :) !! nonsingular upper-triangular matrix
            real(dp), intent(in) :: b(:) !! right-hand side
            real(dp), intent(out) :: solution(:) !! solution of U x = b
            integer, intent(out) :: info !! zero on success

            integer :: i
            integer :: m
            real(dp) :: rhs

            m = size(b)
            solution = 0.0_dp
            do i = m, 1, -1
                if (.not. ieee_is_finite(u(i, i)) .or. abs(u(i, i)) <= tiny(1.0_dp)) then
                    info = i
                    return
                end if
                rhs = b(i)
                if (i < m) rhs = rhs - dot_product(u(i, i + 1:m), solution(i + 1:m))
                solution(i) = rhs / u(i, i)
            end do
            info = 0
        end subroutine solve_upper

        subroutine solve_upper_transpose(u, b, solution, info)
            real(dp), intent(in) :: u(:, :) !! nonsingular upper-triangular matrix
            real(dp), intent(in) :: b(:) !! right-hand side
            real(dp), intent(out) :: solution(:) !! solution of U^T x = b
            integer, intent(out) :: info !! zero on success

            integer :: i
            integer :: m
            real(dp) :: rhs

            m = size(b)
            solution = 0.0_dp
            do i = 1, m
                if (.not. ieee_is_finite(u(i, i)) .or. abs(u(i, i)) <= tiny(1.0_dp)) then
                    info = i
                    return
                end if
                rhs = b(i)
                if (i > 1) rhs = rhs - dot_product(u(1:i - 1, i), solution(1:i - 1))
                solution(i) = rhs / u(i, i)
            end do
            info = 0
        end subroutine solve_upper_transpose

        subroutine estimate_smallest_singular_value(u, s_min, z_min, info)
            real(dp), intent(in) :: u(:, :) !! nonsingular upper-triangular matrix
            real(dp), intent(out) :: s_min !! estimated smallest singular value
            real(dp), intent(out) :: z_min(:) !! estimated right singular vector
            integer, intent(out) :: info !! zero on success

            integer :: i
            integer :: k
            integer :: m
            real(dp) :: pwork(size(z_min))
            real(dp) :: score_minus
            real(dp) :: score_plus
            real(dp) :: v(size(z_min))
            real(dp) :: v_norm
            real(dp) :: w(size(z_min))
            real(dp) :: w_minus
            real(dp) :: w_norm
            real(dp) :: w_plus

            m = size(z_min)
            pwork = 0.0_dp
            w = 0.0_dp
            do k = 1, m
                if (.not. ieee_is_finite(u(k, k)) .or. abs(u(k, k)) <= tiny(1.0_dp)) then
                    info = k
                    return
                end if
                w_plus = (1.0_dp - pwork(k)) / u(k, k)
                w_minus = (-1.0_dp - pwork(k)) / u(k, k)
                score_plus = abs(w_plus)
                score_minus = abs(w_minus)
                do i = k + 1, m
                    score_plus = score_plus + abs(pwork(i) + u(k, i) * w_plus)
                    score_minus = score_minus + abs(pwork(i) + u(k, i) * w_minus)
                end do
                if (score_plus >= score_minus) then
                    w(k) = w_plus
                else
                    w(k) = w_minus
                end if
                do i = k + 1, m
                    pwork(i) = pwork(i) + u(k, i) * w(k)
                end do
            end do

            call solve_upper(u, w, v, info)
            if (info /= 0) return
            v_norm = norm2(v)
            w_norm = norm2(w)
            if (.not. ieee_is_finite(v_norm) .or. v_norm <= tiny(1.0_dp) .or. &
                    .not. ieee_is_finite(w_norm)) then
                info = 1
                return
            end if
            s_min = w_norm / v_norm
            z_min = v / v_norm
            info = 0
        end subroutine estimate_smallest_singular_value

        subroutine singular_leading_submatrix(a, u, k, delta, v, ok_singular)
            real(dp), intent(in) :: a(:, :) !! shifted symmetric matrix
            real(dp), intent(in) :: u(:, :) !! incomplete upper Cholesky factor
            integer, intent(in) :: k !! first non-positive-definite leading minor
            real(dp), intent(out) :: delta !! diagonal shift making the leading block singular
            real(dp), intent(out) :: v(:) !! associated null vector
            logical, intent(out) :: ok_singular !! true on success

            integer :: info_local

            if (k < 1 .or. k > size(v)) then
                delta = 0.0_dp
                v = 0.0_dp
                ok_singular = .false.
                return
            end if
            delta = -a(k, k)
            if (k > 1) delta = delta + sum(u(1:k - 1, k)**2)
            v = 0.0_dp
            v(k) = 1.0_dp
            if (k > 1) then
                call solve_upper(u(1:k - 1, 1:k - 1), -u(1:k - 1, k), &
                    v(1:k - 1), info_local)
                if (info_local /= 0) then
                    ok_singular = .false.
                    return
                end if
            end if
            ok_singular = ieee_is_finite(delta) .and. all(ieee_is_finite(v))
        end subroutine singular_leading_submatrix

        subroutine boundary_intersections(z, d, radius, ta, tb, ok_intersection)
            real(dp), intent(in) :: z(:) !! line origin
            real(dp), intent(in) :: d(:) !! line direction
            real(dp), intent(in) :: radius !! trust-region radius
            real(dp), intent(out) :: ta !! smaller line-sphere parameter
            real(dp), intent(out) :: tb !! larger line-sphere parameter
            logical, intent(out) :: ok_intersection !! true when real intersections exist

            real(dp) :: a
            real(dp) :: aux
            real(dp) :: b
            real(dp) :: c
            real(dp) :: discriminant
            real(dp) :: root_disc
            real(dp) :: tmp

            a = dot_product(d, d)
            b = 2.0_dp * dot_product(z, d)
            c = dot_product(z, z) - radius**2
            discriminant = b * b - 4.0_dp * a * c
            if (a <= 0.0_dp .or. discriminant < 0.0_dp .or. &
                    .not. ieee_is_finite(discriminant)) then
                ok_intersection = .false.
                ta = 0.0_dp
                tb = 0.0_dp
                return
            end if

            root_disc = sqrt(max(0.0_dp, discriminant))
            aux = b + sign(root_disc, b)
            if (abs(aux) > tiny(1.0_dp)) then
                ta = -aux / (2.0_dp * a)
                tb = -2.0_dp * c / aux
            else
                ta = (-b - root_disc) / (2.0_dp * a)
                tb = (-b + root_disc) / (2.0_dp * a)
            end if
            if (ta > tb) then
                tmp = ta
                ta = tb
                tb = tmp
            end if
            ok_intersection = ieee_is_finite(ta) .and. ieee_is_finite(tb)
        end subroutine boundary_intersections

        function damping_update(lambda_low, lambda_high) result(value)
            real(dp), intent(in) :: lambda_low !! lower damping bound
            real(dp), intent(in) :: lambda_high !! upper damping bound
            real(dp) :: value !! safeguarded next damping value

            value = max(sqrt(abs(lambda_low * lambda_high)), &
                lambda_low + 0.01_dp * (lambda_high - lambda_low))
        end function damping_update

    end function trust_exact

    function cobyla(fun, x0, lower, upper, constraints, meq, mineq, rhobeg, tol, &
            maxiter, catol, f_target) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        real(dp), intent(in), optional :: lower(:) !! lower bounds; nonfinite values mean unbounded
        real(dp), intent(in), optional :: upper(:) !! upper bounds; nonfinite values mean unbounded
        procedure(multivariate_constraint_function), optional :: constraints !! equality then >=0 inequalities
        integer, intent(in), optional :: meq !! number of leading equality constraints
        integer, intent(in), optional :: mineq !! number of trailing inequality constraints
        real(dp), intent(in), optional :: rhobeg !! initial trust-region radius, default 1
        real(dp), intent(in), optional :: tol !! lower bound on trust-region radius, default 1e-4
        integer, intent(in), optional :: maxiter !! maximum objective evaluations, default 1000
        real(dp), intent(in), optional :: catol !! absolute constraint-violation tolerance
        real(dp), intent(in), optional :: f_target !! target objective value for early termination
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: cfinal(:)
        real(dp), allocatable :: cuser(:)
        real(dp), allocatable :: reduced_x(:)
        real(dp), allocatable :: full_x(:)
        integer, allocatable :: free_indices(:)
        logical, allocatable :: fixed(:)
        real(dp), allocatable :: lb(:)
        real(dp), allocatable :: ub(:)
        real(dp), allocatable :: xfinal(:)
        real(dp), allocatable :: xstart(:)
        real(dp) :: catol_use
        real(dp) :: f_target_use
        real(dp) :: ffinal
        real(dp) :: maxcv
        real(dp) :: rhobeg_use
        real(dp) :: tol_use
        integer :: i
        integer :: info
        integer :: maxiter_use
        integer :: meq_use
        integer :: mineq_use
        integer :: mtotal
        integer :: muser
        integer :: n
        integer :: nbound
        integer :: nf
        integer :: nit
        integer, allocatable :: nbd(:)
        logical :: invalid_bounds
        logical :: success

        n = size(x0)
        call init_result(res, 'cobyla', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        meq_use = 0
        if (present(meq)) meq_use = meq
        mineq_use = 0
        if (present(mineq)) mineq_use = mineq
        if (meq_use < 0 .or. mineq_use < 0) return
        muser = meq_use + mineq_use
        if (muser > 0 .and. .not. present(constraints)) return
        if (muser == 0 .and. present(constraints)) return

        rhobeg_use = 1.0_dp
        if (present(rhobeg)) rhobeg_use = rhobeg
        tol_use = 1.0e-4_dp
        if (present(tol)) tol_use = tol
        catol_use = sqrt(epsilon(1.0_dp))
        if (present(catol)) catol_use = catol
        f_target_use = -huge(1.0_dp)
        if (present(f_target)) f_target_use = f_target
        maxiter_use = 1000
        if (present(maxiter)) maxiter_use = maxiter
        if (.not. ieee_is_finite(rhobeg_use) .or. rhobeg_use <= 0.0_dp) return
        if (.not. ieee_is_finite(tol_use) .or. tol_use <= 0.0_dp) return
        if (.not. ieee_is_finite(catol_use) .or. catol_use < 0.0_dp) return
        if (maxiter_use < 1) return
        rhobeg_use = max(rhobeg_use, tol_use)

        allocate(lb(n), ub(n), nbd(n))
        call setup_box_bounds(n, lower, upper, lb, ub, nbd, invalid_bounds)
        if (invalid_bounds) return
        fixed = nbd == 2 .and. lb == ub
        free_indices = pack([(i, i = 1, n)], .not. fixed)
        nbound = count((nbd == 1 .or. nbd == 2) .and. .not. fixed) + &
            count((nbd == 2 .or. nbd == 3) .and. .not. fixed)
        mtotal = nbound + 2 * meq_use + mineq_use
        allocate(cfinal(mtotal), cuser(muser), xfinal(n), xstart(n))
        xstart = x0
        ! SciPy/PRIMA projects bound-only starting points onto the feasible box
        ! before entering COBYLA.  Keep nonlinear-constraint starts unchanged.
        if (muser == 0) then
            do i = 1, n
                select case (nbd(i))
                case (1)
                    xstart(i) = max(xstart(i), lb(i))
                case (2)
                    xstart(i) = min(max(xstart(i), lb(i)), ub(i))
                case (3)
                    xstart(i) = min(xstart(i), ub(i))
                end select
            end do
        end if

        where (fixed) xstart = lb
        allocate(full_x(n), reduced_x(size(free_indices)))
        if (size(free_indices) == 0) then
            xfinal = xstart
            call eval_cobyla(xfinal, ffinal, cfinal)
            maxcv = max(0.0_dp, maxval(cfinal))
            nf = 1
            nit = 0
            info = cobyla_small_tr_radius
        else
            call prima_cobyla(eval_reduced_cobyla, xstart(free_indices), mtotal, &
                rhobeg_use, tol_use, maxiter_use, catol_use, f_target_use, reduced_x, &
                ffinal, cfinal, maxcv, nf, info, nit)
            xfinal = xstart
            xfinal(free_indices) = reduced_x
        end if

        res%x = xfinal
        res%fun = ffinal
        res%maxcv = maxcv
        res%nfev = nf
        res%nit = nit
        success = (info == cobyla_small_tr_radius .or. info == cobyla_ftarget_achieved) .and. &
            maxcv <= catol_use
        if (maxcv > catol_use) then
            call finish_result(res, info, .false., &
                'Did not converge to a solution satisfying the constraints.')
        else
            select case (info)
            case (cobyla_small_tr_radius)
                call finish_result(res, info, success, &
                    'Trust region radius reaches its lower bound.')
            case (cobyla_ftarget_achieved)
                call finish_result(res, info, success, 'Target function value achieved.')
            case (cobyla_maxfun_reached)
                call finish_result(res, info, .false., &
                    'Maximum number of function evaluations has been exceeded.')
            case (cobyla_damaging_rounding)
                call finish_result(res, info, .false., 'Rounding errors are becoming damaging.')
            case (cobyla_maxtr_reached)
                call finish_result(res, info, .false., &
                    'Maximum number of trust-region iterations has been exceeded.')
            case default
                call finish_result(res, info, .false., 'COBYLA terminated without convergence.')
            end select
        end if

    contains

        subroutine eval_reduced_cobyla(xv, fv, cv)
            implicit none
            real(dp), intent(in) :: xv(:) !! free coordinates requested by the kernel
            real(dp), intent(out) :: fv !! objective value in the full parameter space
            real(dp), intent(out) :: cv(:) !! constraints evaluated with fixed coordinates restored

            full_x = xstart
            full_x(free_indices) = xv
            call eval_cobyla(full_x, fv, cv)
        end subroutine eval_reduced_cobyla

        subroutine eval_cobyla(xv, fv, cv)
            real(dp), intent(in) :: xv(:) !! point requested by the COBYLA kernel
            real(dp), intent(out) :: fv !! objective value at xv
            real(dp), intent(out) :: cv(:) !! PRIMA-style constraints, each <= 0

            integer :: k
            integer :: pos

            fv = fun(xv)
            pos = 0
            do k = 1, n
                if (fixed(k)) cycle
                if (nbd(k) == 1 .or. nbd(k) == 2) then
                    pos = pos + 1
                    cv(pos) = lb(k) - xv(k)
                end if
            end do
            do k = 1, n
                if (fixed(k)) cycle
                if (nbd(k) == 2 .or. nbd(k) == 3) then
                    pos = pos + 1
                    cv(pos) = xv(k) - ub(k)
                end if
            end do
            if (muser > 0) then
                call invoke_constraints(constraints, xv, cuser)
                do k = 1, meq_use
                    pos = pos + 1
                    cv(pos) = cuser(k)
                    pos = pos + 1
                    cv(pos) = -cuser(k)
                end do
                do k = 1, mineq_use
                    pos = pos + 1
                    cv(pos) = -cuser(meq_use + k)
                end do
            end if
        end subroutine eval_cobyla

    end function cobyla

    function tnc(fun, x0, jac, lower, upper, ftol, xtol, gtol, eps, maxfun, &
            max_cg_iter) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic gradient
        real(dp), intent(in), optional :: lower(:) !! lower bounds; nonfinite values mean unbounded
        real(dp), intent(in), optional :: upper(:) !! upper bounds; nonfinite values mean unbounded
        real(dp), intent(in), optional :: ftol !! relative objective-reduction tolerance
        real(dp), intent(in), optional :: xtol !! scaled-step Euclidean-norm tolerance
        real(dp), intent(in), optional :: gtol !! scaled projected-gradient Euclidean-norm tolerance
        real(dp), intent(in), optional :: eps !! absolute finite-difference step, default 1e-8
        integer, intent(in), optional :: maxfun !! maximum objective evaluations
        integer, intent(in), optional :: max_cg_iter !! inner CG cap; negative uses default
        type(multivariate_minimize_result) :: res

        real(dp), parameter :: armijo = 1.0e-4_dp
        real(dp), parameter :: min_alpha = 16.0_dp * epsilon(1.0_dp)
        integer, parameter :: max_line_iter = 60
        integer :: cg_iter
        integer :: i
        integer :: j
        integer :: line_iter
        integer :: max_cg_use
        integer :: maxfun_use
        integer :: n
        integer, allocatable :: nbd(:)
        logical :: accepted
        logical :: all_constant
        logical :: invalid
        real(dp) :: accuracy
        real(dp) :: alpha
        real(dp) :: alpha_max
        real(dp) :: beta
        real(dp) :: curvature
        real(dp) :: eps_use
        real(dp) :: f_old
        real(dp) :: f_trial
        real(dp) :: ftol_use
        real(dp) :: gnorm
        real(dp) :: gtol_use
        real(dp) :: pnorm
        real(dp) :: pg_norm
        real(dp) :: residual_norm
        real(dp) :: rz_new
        real(dp) :: rz_old
        real(dp) :: sigma
        real(dp) :: slope
        real(dp) :: step_norm
        real(dp) :: xtol_use
        real(dp) :: yhy
        real(dp) :: ynorm
        real(dp) :: ys
        real(dp) :: snorm
        real(dp), allocatable :: direction(:)
        real(dp), allocatable :: g_scaled(:)
        real(dp), allocatable :: g_old(:)
        real(dp), allocatable :: g_trial(:)
        real(dp), allocatable :: hess_direction(:)
        real(dp), allocatable :: hinv(:, :)
        real(dp), allocatable :: hy(:)
        real(dp), allocatable :: lb(:)
        real(dp), allocatable :: offset(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: pg(:)
        real(dp), allocatable :: precond(:)
        real(dp), allocatable :: residual(:)
        real(dp), allocatable :: scale(:)
        real(dp), allocatable :: svec(:)
        real(dp), allocatable :: ub(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_old(:)
        real(dp), allocatable :: x_trial(:)
        real(dp), allocatable :: yvec(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: z_old(:)
        real(dp), allocatable :: z_trial(:)

        n = size(x0)
        call init_result(res, 'tnc', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        allocate(lb(n), ub(n), nbd(n), x(n), x_old(n), x_trial(n), z(n), z_trial(n))
        allocate(scale(n), offset(n), g_scaled(n), g_old(n), g_trial(n), pg(n), p(n), &
            residual(n), direction(n), hess_direction(n), hinv(n, n), hy(n), svec(n), &
            precond(n), yvec(n), z_old(n))

        call setup_box_bounds(n, lower, upper, lb, ub, nbd, invalid)
        if (invalid) then
            call finish_result(res, -1, .false., 'Infeasible or invalid bounds')
            return
        end if

        x = x0
        call clip_to_box(x, lb, ub, nbd)

        eps_use = 1.0e-8_dp
        if (present(eps)) eps_use = eps
        if (.not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp) return

        accuracy = sqrt(epsilon(1.0_dp))
        ftol_use = 0.0_dp
        if (present(ftol)) ftol_use = max(0.0_dp, ftol)
        xtol_use = sqrt(epsilon(1.0_dp))
        if (present(xtol)) then
            if (xtol >= 0.0_dp) xtol_use = xtol
        end if
        gtol_use = 1.0e-2_dp * sqrt(accuracy)
        if (present(gtol)) then
            if (gtol >= 0.0_dp) gtol_use = gtol
        end if
        if (.not. ieee_is_finite(ftol_use) .or. .not. ieee_is_finite(xtol_use) .or. &
                .not. ieee_is_finite(gtol_use)) return

        maxfun_use = max(100, 10 * n)
        if (present(maxfun)) maxfun_use = maxfun
        if (maxfun_use < 1) then
            res%x = x
            call finish_result(res, 3, .false., &
                'Maximum number of function evaluations reached')
            return
        end if

        max_cg_use = max(1, min(50, n / 2))
        if (present(max_cg_iter)) then
            if (max_cg_iter >= 0) max_cg_use = min(max_cg_iter, n)
        end if

        call setup_tnc_scaling(x, lb, ub, nbd, scale, offset)
        z = (x - offset) / scale
        hinv = 0.0_dp
        do i = 1, n
            hinv(i, i) = 1.0_dp
        end do

        call evaluate_fg(x, res%fun, g_trial)
        if (.not. ieee_is_finite(res%fun) .or. .not. all(ieee_is_finite(g_trial))) then
            res%x = x
            res%jac = g_trial
            call finish_result(res, 4, .false., 'Nonfinite objective or gradient')
            return
        end if
        all_constant = .true.
        do i = 1, n
            if (nbd(i) /= 2 .or. .not. same_real(lb(i), ub(i))) then
                all_constant = .false.
                exit
            end if
        end do
        if (all_constant) then
            res%x = x
            res%jac = g_trial
            call finish_result(res, 5, .false., &
                'All lower bounds are equal to the upper bounds')
            return
        end if
        g_scaled = g_trial * scale

        do
            call projected_scaled_gradient(x, g_scaled, scale, lb, ub, nbd, pg)
            pg_norm = sqrt(max(0.0_dp, dot_product(pg, pg)))
            gnorm = pg_norm
            if (gnorm <= gtol_use) then
                call finish_result(res, 0, .true., 'Local minimum reached (|pg| ~= 0)')
                exit
            end if
            if (res%nfev >= maxfun_use) then
                call finish_result(res, 3, .false., &
                    'Maximum number of function evaluations reached')
                exit
            end if

            if (present(jac)) then
                residual = -pg
                p = 0.0_dp
                if (max_cg_use == 0) then
                    p = residual
                else
                    precond = matmul(hinv, residual)
                    where (same_active_component(x, g_scaled, scale, lb, ub, nbd))
                        precond = 0.0_dp
                    end where
                    rz_old = dot_product(residual, precond)
                    if (rz_old <= 0.0_dp) then
                        precond = residual
                        rz_old = dot_product(residual, residual)
                    end if
                    direction = precond
                    do cg_iter = 1, max_cg_use
                        pnorm = sqrt(max(0.0_dp, dot_product(direction, direction)))
                        if (same_real(pnorm, 0.0_dp)) exit

                        sigma = sqrt(accuracy) * (1.0_dp + sqrt(dot_product(z, z))) / pnorm
                        alpha_max = scaled_step_to_bounds(x, direction, scale, lb, ub, nbd)
                        if (ieee_is_finite(alpha_max)) then
                            sigma = min(sigma, 0.5_dp * max(0.0_dp, alpha_max))
                        end if
                        if (sigma <= min_alpha) exit

                        z_trial = z + sigma * direction
                        x_trial = offset + scale * z_trial
                        call clip_to_box(x_trial, lb, ub, nbd)
                        call evaluate_gradient_only(x_trial, f_trial, g_trial)
                        if (.not. all(ieee_is_finite(g_trial))) exit

                        hess_direction = (g_trial * scale - g_scaled) / sigma
                        where (same_active_component(x, g_scaled, scale, lb, ub, nbd))
                            hess_direction = 0.0_dp
                        end where
                        curvature = dot_product(direction, hess_direction)
                        if (curvature <= 1.0e-14_dp * max(1.0_dp, pnorm * pnorm)) then
                            if (cg_iter == 1) p = precond
                            exit
                        end if

                        alpha = rz_old / curvature
                        p = p + alpha * direction
                        residual = residual - alpha * hess_direction
                        residual_norm = dot_product(residual, residual)
                        if (sqrt(max(0.0_dp, residual_norm)) <= &
                                min(0.5_dp, sqrt(max(pg_norm, epsilon(1.0_dp)))) * pg_norm) exit
                        precond = matmul(hinv, residual)
                        where (same_active_component(x, g_scaled, scale, lb, ub, nbd))
                            precond = 0.0_dp
                        end where
                        rz_new = dot_product(residual, precond)
                        if (rz_new <= 0.0_dp) exit
                        beta = rz_new / rz_old
                        direction = precond + beta * direction
                        rz_old = rz_new
                    end do
                end if
            else
                p = -matmul(hinv, pg)
                where (same_active_component(x, g_scaled, scale, lb, ub, nbd))
                    p = 0.0_dp
                end where
            end if

            if (dot_product(g_scaled, p) >= 0.0_dp .or. &
                    same_real(dot_product(p, p), 0.0_dp)) then
                p = -pg
                hinv = 0.0_dp
                do i = 1, n
                    hinv(i, i) = 1.0_dp
                end do
            end if

            alpha_max = scaled_step_to_bounds(x, p, scale, lb, ub, nbd)
            alpha = 1.0_dp
            if (ieee_is_finite(alpha_max)) alpha = min(alpha, alpha_max)
            if (alpha <= 0.0_dp) then
                p = -pg
                alpha_max = scaled_step_to_bounds(x, p, scale, lb, ub, nbd)
                alpha = 1.0_dp
                if (ieee_is_finite(alpha_max)) alpha = min(alpha, alpha_max)
            end if
            if (alpha <= 0.0_dp) then
                call finish_result(res, 6, .false., 'Unable to progress')
                exit
            end if

            x_old = x
            z_old = z
            g_old = g_scaled
            f_old = res%fun
            slope = dot_product(g_scaled, p)
            accepted = .false.
            do line_iter = 1, max_line_iter
                z_trial = z + alpha * p
                x_trial = offset + scale * z_trial
                call clip_to_box(x_trial, lb, ub, nbd)
                if (same_vector(x_trial, x_old)) then
                    alpha = 0.5_dp * alpha
                    if (alpha <= min_alpha) exit
                    cycle
                end if

                f_trial = fun(x_trial)
                res%nfev = res%nfev + 1
                if (ieee_is_finite(f_trial)) then
                    if (f_trial <= f_old + armijo * alpha * slope) then
                        accepted = .true.
                        exit
                    end if
                end if
                if (res%nfev >= maxfun_use) exit
                alpha = 0.5_dp * alpha
                if (alpha <= min_alpha) exit
            end do

            if (.not. accepted) then
                if (res%nfev >= maxfun_use) then
                    call finish_result(res, 3, .false., &
                        'Maximum number of function evaluations reached')
                else
                    call finish_result(res, 4, .false., 'Linear search failed')
                end if
                exit
            end if

            x = x_trial
            z = (x - offset) / scale
            res%fun = f_trial
            call evaluate_gradient_from_base(x, res%fun, g_trial)
            g_scaled = g_trial * scale
            svec = z - z_old
            yvec = g_scaled - g_old
            ys = dot_product(yvec, svec)
            snorm = sqrt(max(0.0_dp, dot_product(svec, svec)))
            ynorm = sqrt(max(0.0_dp, dot_product(yvec, yvec)))
            if (ys > 1.0e-12_dp * snorm * ynorm .and. ys > 0.0_dp) then
                hy = matmul(hinv, yvec)
                yhy = dot_product(yvec, hy)
                do j = 1, n
                    do i = 1, n
                        hinv(i, j) = hinv(i, j) + &
                            (1.0_dp + yhy / ys) * svec(i) * svec(j) / ys - &
                            (svec(i) * hy(j) + hy(i) * svec(j)) / ys
                    end do
                end do
            else
                hinv = 0.0_dp
                do i = 1, n
                    hinv(i, i) = 1.0_dp
                end do
            end if
            res%nit = res%nit + 1

            if (ftol_use > 0.0_dp) then
                if (abs(f_old - res%fun) <= ftol_use * &
                        max(1.0_dp, abs(f_old), abs(res%fun))) then
                    call finish_result(res, 1, .true., &
                        'Converged (|f_n-f_(n-1)| ~= 0)')
                    exit
                end if
            end if

            step_norm = sqrt(max(0.0_dp, &
                dot_product((x - x_old) / scale, (x - x_old) / scale)))
            if (step_norm <= xtol_use) then
                call finish_result(res, 2, .true., &
                    'Converged (|x_n-x_(n-1)| ~= 0)')
                exit
            end if

            if (res%nfev >= maxfun_use) then
                call finish_result(res, 3, .false., &
                    'Maximum number of function evaluations reached')
                exit
            end if
        end do

        res%x = x
        res%jac = g_scaled / scale

    contains

        subroutine evaluate_fg(xv, fv, gv)
            real(dp), intent(in) :: xv(:) !! point to evaluate
            real(dp), intent(out) :: fv !! objective value
            real(dp), intent(out) :: gv(:) !! objective gradient

            fv = fun(xv)
            res%nfev = res%nfev + 1
            call evaluate_gradient_from_base(xv, fv, gv)
        end subroutine evaluate_fg

        subroutine evaluate_gradient_only(xv, fv, gv)
            real(dp), intent(in) :: xv(:) !! point at which a gradient is needed
            real(dp), intent(out) :: fv !! objective value when numerical differences require it
            real(dp), intent(out) :: gv(:) !! objective gradient

            if (present(jac)) then
                fv = res%fun
                call invoke_gradient(jac, xv, gv)
                res%njev = res%njev + 1
            else
                fv = fun(xv)
                res%nfev = res%nfev + 1
                call evaluate_gradient_from_base(xv, fv, gv)
            end if
        end subroutine evaluate_gradient_only

        subroutine evaluate_gradient_from_base(xv, fv, gv)
            real(dp), intent(in) :: xv(:) !! point at which the base objective is known
            real(dp), intent(in) :: fv !! objective value at xv
            real(dp), intent(out) :: gv(:) !! objective gradient

            real(dp) :: fp
            real(dp) :: h
            real(dp) :: lower_i
            real(dp) :: upper_i
            real(dp) :: xp(n)
            integer :: j

            if (present(jac)) then
                call invoke_gradient(jac, xv, gv)
                res%njev = res%njev + 1
                return
            end if

            xp = xv
            do j = 1, n
                lower_i = 0.0_dp
                upper_i = 0.0_dp
                if (nbd(j) == 1 .or. nbd(j) == 2) lower_i = lb(j)
                if (nbd(j) == 2 .or. nbd(j) == 3) upper_i = ub(j)
                call bounded_fd_step(xv(j), lower_i, upper_i, nbd(j), eps_use, h)
                if (same_real(h, 0.0_dp)) then
                    gv(j) = 0.0_dp
                else
                    xp = xv
                    xp(j) = xv(j) + h
                    fp = fun(xp)
                    res%nfev = res%nfev + 1
                    gv(j) = (fp - fv) / h
                end if
            end do
            res%njev = res%njev + 1
        end subroutine evaluate_gradient_from_base

    end function tnc

    recursive function trust_constr(fun, x0, jac, hess, lower, upper, constraints, &
            constraint_jac, meq, mineq, gtol, xtol, eps, maxiter, &
            initial_trust_radius, barrier_tol, initial_barrier_parameter, &
            initial_barrier_tolerance, initial_constr_penalty, initial_tr_radius, &
            constraint_hess, constraint_lower, constraint_upper, &
            constraint_keep_feasible, bounds_keep_feasible, constraint_objects, objective_fd_scheme, &
            objective_finite_diff_rel_step, hessp, objective_hess_update, callback, &
            objective_hess_fd_scheme) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic objective gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic objective Hessian
        procedure(multivariate_hessian_product_function), optional :: hessp !! objective Hessian-vector product
        real(dp), intent(in), optional :: lower(:) !! lower variable bounds; nonfinite means unbounded
        real(dp), intent(in), optional :: upper(:) !! upper variable bounds; nonfinite means unbounded
        procedure(multivariate_constraint_function), optional :: constraints !! equalities then inequalities >= 0
        procedure(multivariate_constraint_jacobian), optional :: constraint_jac !! combined constraint Jacobian
        integer, intent(in), optional :: meq !! number of leading equality constraints
        integer, intent(in), optional :: mineq !! number of trailing inequality constraints
        real(dp), intent(in), optional :: constraint_lower(:) !! lower bounds on user constraint values
        real(dp), intent(in), optional :: constraint_upper(:) !! upper bounds on user constraint values
        logical, intent(in), optional :: constraint_keep_feasible(:) !! keep bounded user constraint components feasible
        logical, intent(in), optional :: bounds_keep_feasible(:) !! keep variable-bound components feasible
        type(nonlinear_constraint), intent(in), optional :: constraint_objects(:) !! independent constraint blocks
        character(len=*), intent(in), optional :: objective_fd_scheme !! trust-constr objective '2-point' or '3-point'
        real(dp), intent(in), optional :: objective_finite_diff_rel_step(:) !! per-variable objective relative steps
        character(len=*), intent(in), optional :: objective_hess_update !! objective Hessian BFGS or SR1 update
        character(len=*), intent(in), optional :: objective_hess_fd_scheme !! 2-point or 3-point Hessian of analytic jac
        procedure(trust_constr_iteration_callback), optional :: callback !! iteration callback, may stop the solve
        real(dp), intent(in), optional :: gtol !! Lagrangian-gradient/feasibility tolerance
        real(dp), intent(in), optional :: xtol !! trust-region-radius tolerance
        real(dp), intent(in), optional :: eps !! finite-difference step
        integer, intent(in), optional :: maxiter !! maximum total SQP iterations
        real(dp), intent(in), optional :: initial_trust_radius !! legacy SciFort initial trust-region radius spelling
        real(dp), intent(in), optional :: initial_constr_penalty !! initial merit-function constraint penalty
        real(dp), intent(in), optional :: initial_tr_radius !! SciPy trust-constr initial trust radius spelling
        real(dp), intent(in), optional :: barrier_tol !! threshold for terminating the barrier sequence
        real(dp), intent(in), optional :: initial_barrier_parameter !! initial logarithmic-barrier parameter
        real(dp), intent(in), optional :: initial_barrier_tolerance !! initial barrier-subproblem tolerance
        procedure(multivariate_constraint_hessian), optional :: constraint_hess !! analytic constraint Hessian
        type(multivariate_minimize_result) :: res

        real(dp), parameter :: barrier_decay = 0.2_dp
        integer, parameter :: max_barrier_stages = 32
        integer :: i
        integer :: j
        integer :: mbar
        integer :: mfix
        integer :: meq_use
        integer :: mineq_use
        integer :: morig
        integer :: mraw
        integer :: constraint_eq_count
        integer :: constraint_upper_count
        integer :: constraint_lower_count
        integer :: constraint_interval_count
        integer :: mwork
        integer :: n
        integer :: nwork
        integer :: remaining
        integer :: slack
        integer :: stage
        integer :: total_cg_niter
        integer :: total_nfev
        integer :: total_nhev
        integer :: total_nit
        integer :: total_njev
        integer :: kcache
        integer :: maxiter_use
        integer :: hess_fd_njev
        integer, allocatable :: hess_fd_nbd(:)
        integer :: user_constr_nfev
        integer :: user_constr_njev
        integer :: user_constr_nhev
        integer, allocatable :: nbd(:)
        integer, allocatable :: strict_fd_nbd(:)
        integer, allocatable :: constraint_eq_idx(:)
        integer, allocatable :: constraint_upper_idx(:)
        integer, allocatable :: constraint_lower_idx(:)
        integer, allocatable :: constraint_interval_idx(:)
        logical :: hess_fd_three_point
        logical :: hess_fd_invalid
        logical :: has_bounds_block
        logical :: has_constraint_bounds
        logical, allocatable :: canonical_keep_feasible(:)
        logical, allocatable :: barrier_keep_feasible(:)
        logical, allocatable :: bound_keep_feasible_use(:)
        logical, allocatable :: raw_constraint_keep_feasible(:)
        logical :: raw_constraint_cache_valid
        logical :: raw_constraint_jac_cache_valid
        logical :: invalid_bounds
        logical :: final_success
        logical :: optimality_ok
        logical :: objective_three_point
        real(dp) :: barrier_parameter
        real(dp) :: barrier_tolerance
        real(dp) :: barrier_tol_use
        real(dp) :: eps_use
        real(dp) :: objective_fd_step_use
        real(dp) :: gtol_use
        real(dp) :: initial_barrier_parameter_use
        real(dp) :: initial_barrier_tolerance_use
        real(dp) :: initial_constr_penalty_use
        real(dp) :: initial_radius_use
        real(dp) :: maxcv
        real(dp) :: barrier_optimality
        real(dp) :: cached_barrier_fun
        real(dp) :: original_fun
        real(dp) :: hess_fd_step
        real(dp) :: fd_step
        real(dp) :: fd_request
        real(dp) :: fd_value
        real(dp) :: stage_radius
        real(dp) :: xtol_use
        real(dp), allocatable :: abase(:, :)
        real(dp), allocatable :: cached_barrier_grad(:)
        real(dp), allocatable :: cached_barrier_hess(:, :)
        real(dp), allocatable :: cached_barrier_constr(:)
        real(dp), allocatable :: cached_barrier_jac(:, :)
        real(dp), allocatable :: barrier_lagrangian_grad(:)
        real(dp), allocatable :: cbase(:)
        real(dp), allocatable :: c_fd(:)
        real(dp), allocatable :: x_fd(:)
        real(dp), allocatable :: hess_fd_lb(:)
        real(dp), allocatable :: hess_fd_ub(:)
        real(dp), allocatable :: lb(:)
        real(dp), allocatable :: ub(:)
        real(dp), allocatable :: strict_fd_lb(:)
        real(dp), allocatable :: strict_fd_ub(:)
        real(dp), allocatable :: constraint_lb(:)
        real(dp), allocatable :: constraint_ub(:)
        real(dp), allocatable :: raw_constraint_cache(:)
        real(dp), allocatable :: raw_constraint_cache_x(:)
        real(dp), allocatable :: raw_constraint_jac_cache(:, :)
        real(dp), allocatable :: raw_constraint_jac_cache_x(:)
        real(dp), allocatable :: stage_barrier_hess(:, :)
        real(dp), allocatable :: initial_constraint_hess(:, :)
        real(dp), allocatable :: zero_constraint_multipliers(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: z_fd(:)
        type(multivariate_minimize_result) :: inner

        n = size(x0)
        call init_result(res, 'trust-constr', n)
        if (present(objective_hess_update)) then
            select case (lowercase(trim(objective_hess_update)))
            case ('bfgs', 'sr1')
                if (present(hess) .or. present(hessp)) then
                    call finish_result(res, -1, .false., &
                        'objective_hess_update conflicts with exact objective Hessian')
                    return
                end if
            case default
                call finish_result(res, -1, .false., 'Invalid objective_hess_update')
                return
            end select
        end if
        if (present(objective_hess_fd_scheme)) then
            ! Differentiate the analytic objective gradient columnwise, using
            ! feasible stencils for finite variable bounds. Reenter the common
            ! dense SQP controller with an explicit Hessian callback so all
            ! existing barrier, constraint-object, and diagnostics paths apply.
            if (objective_hess_fd_scheme /= '2-point' .and. &
                    objective_hess_fd_scheme /= '3-point') then
                call finish_result(res, -1, .false., 'Invalid objective_hess_fd_scheme')
                return
            end if
            if (present(hess) .or. present(hessp) .or. present(objective_hess_update)) then
                call finish_result(res, -1, .false., &
                    'objective_hess_fd_scheme conflicts with objective Hessian')
                return
            end if
            if (.not. present(jac)) then
                call finish_result(res, -1, .false., &
                    'objective_hess_fd_scheme requires analytic jac')
                return
            end if
            if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
            if (present(objective_finite_diff_rel_step)) then
                if (size(objective_finite_diff_rel_step) /= n) then
                    call finish_result(res, -1, .false., &
                        'Invalid objective_finite_diff_rel_step size')
                    return
                end if
                if (.not. all(ieee_is_finite(objective_finite_diff_rel_step))) then
                    call finish_result(res, -1, .false., &
                        'Nonfinite objective_finite_diff_rel_step')
                    return
                end if
            end if
            allocate(hess_fd_lb(n), hess_fd_ub(n), hess_fd_nbd(n))
            call setup_box_bounds(n, lower, upper, hess_fd_lb, hess_fd_ub, &
                hess_fd_nbd, hess_fd_invalid)
            if (hess_fd_invalid) then
                call finish_result(res, -1, .false., 'Invalid trust-constr bounds')
                return
            end if
            hess_fd_three_point = objective_hess_fd_scheme == '3-point'
            hess_fd_step = sqrt(epsilon(1.0_dp))
            if (hess_fd_three_point) hess_fd_step = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
            hess_fd_njev = 0
            res = trust_constr(fun, x0, jac=jac, hess=finite_difference_objective_hessian, &
                lower=lower, upper=upper, constraints=constraints, &
                constraint_jac=constraint_jac, meq=meq, mineq=mineq, &
                gtol=gtol, xtol=xtol, eps=eps, maxiter=maxiter, &
                initial_trust_radius=initial_trust_radius, barrier_tol=barrier_tol, &
                initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                initial_constr_penalty=initial_constr_penalty, &
                initial_tr_radius=initial_tr_radius, constraint_hess=constraint_hess, &
                constraint_lower=constraint_lower, constraint_upper=constraint_upper, &
                constraint_keep_feasible=constraint_keep_feasible, &
                bounds_keep_feasible=bounds_keep_feasible, constraint_objects=constraint_objects, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, &
                callback=callback)
            res%njev = res%njev + hess_fd_njev
            ! SciPy counts numerical Hessians as gradient work (njev), not
            ! exact Hessian callback work (nhev).
            res%nhev = 0
            return
        end if
        if (present(hessp) .and. .not. present(hess)) then
            ! The dense SQP kernel needs an explicit Hessian. Assemble it from
            ! n Hessian-vector products at each Hessian evaluation. Forwarding
            ! through the regular analytic-Hessian path also preserves its
            ! barrier, nonlinear-block, and evaluation-counting behavior.
            ! As in SciPy, an explicit hess takes precedence over hessp.
            res = trust_constr(fun, x0, jac=jac, hess=assembled_objective_hessian, &
                lower=lower, upper=upper, constraints=constraints, &
                constraint_jac=constraint_jac, meq=meq, mineq=mineq, &
                gtol=gtol, xtol=xtol, eps=eps, maxiter=maxiter, &
                initial_trust_radius=initial_trust_radius, barrier_tol=barrier_tol, &
                initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                initial_constr_penalty=initial_constr_penalty, &
                initial_tr_radius=initial_tr_radius, constraint_hess=constraint_hess, &
                constraint_lower=constraint_lower, constraint_upper=constraint_upper, &
                constraint_keep_feasible=constraint_keep_feasible, &
                bounds_keep_feasible=bounds_keep_feasible, constraint_objects=constraint_objects, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, &
                objective_hess_update=objective_hess_update, &
                callback=callback)
            return
        end if
        user_constr_nfev = 0
        user_constr_njev = 0
        user_constr_nhev = 0
        has_bounds_block = present(lower) .or. present(upper)
        has_constraint_bounds = present(constraint_lower) .or. present(constraint_upper)
        raw_constraint_cache_valid = .false.
        raw_constraint_jac_cache_valid = .false.
        objective_three_point = .false.
        if (present(objective_fd_scheme)) then
            if (objective_fd_scheme /= '2-point' .and. objective_fd_scheme /= '3-point') then
                call finish_result(res, -1, .false., 'Invalid trust-constr objective_fd_scheme')
                return
            end if
            objective_three_point = objective_fd_scheme == '3-point'
        end if
        if (present(objective_finite_diff_rel_step)) then
            if (size(objective_finite_diff_rel_step) /= n) then
                call finish_result(res, -1, .false., &
                    'Invalid objective_finite_diff_rel_step size')
                return
            end if
            if (.not. all(ieee_is_finite(objective_finite_diff_rel_step))) then
                call finish_result(res, -1, .false., &
                    'Nonfinite objective_finite_diff_rel_step')
                return
            end if
        end if
        if (present(constraint_objects)) then
            if (present(constraints) .or. present(constraint_jac) .or. &
                    present(constraint_hess) .or. present(meq) .or. present(mineq) .or. &
                    present(constraint_lower) .or. present(constraint_upper) .or. &
                    present(constraint_keep_feasible)) then
                call finish_result(res, -1, .false., &
                    'constraint_objects cannot be combined with legacy constraints')
                return
            end if
            res = trust_constr_objects(fun, x0, constraint_objects, jac=jac, hess=hess, &
                lower=lower, upper=upper, gtol=gtol, xtol=xtol, eps=eps, &
                maxiter=maxiter, initial_trust_radius=initial_trust_radius, &
                initial_tr_radius=initial_tr_radius, initial_constr_penalty=initial_constr_penalty, &
                barrier_tol=barrier_tol, initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                bounds_keep_feasible=bounds_keep_feasible, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, &
                objective_hess_update=objective_hess_update, &
                callback=callback)
            return
        end if
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        meq_use = 0
        if (present(meq)) meq_use = meq
        mineq_use = 0
        if (present(mineq)) mineq_use = mineq
        if (meq_use < 0 .or. mineq_use < 0) return

        if (has_constraint_bounds) then
            if (present(meq) .or. present(mineq)) then
                call finish_result(res, -1, .false., &
                    'constraint_lower/constraint_upper cannot be combined with meq/mineq')
                return
            end if
            if (.not. present(constraints)) then
                call finish_result(res, -1, .false., &
                    'constraint_lower/constraint_upper require a constraint callback')
                return
            end if
            if (present(constraint_lower)) then
                mraw = size(constraint_lower)
            else
                mraw = size(constraint_upper)
            end if
            if (mraw <= 0) then
                call finish_result(res, -1, .false., 'Empty trust-constr constraint bounds')
                return
            end if
            if (present(constraint_lower)) then
                if (size(constraint_lower) /= mraw .or. any(ieee_is_nan(constraint_lower))) then
                    call finish_result(res, -1, .false., 'Invalid trust-constr constraint lower bounds')
                    return
                end if
            end if
            if (present(constraint_upper)) then
                if (size(constraint_upper) /= mraw .or. any(ieee_is_nan(constraint_upper))) then
                    call finish_result(res, -1, .false., 'Invalid trust-constr constraint upper bounds')
                    return
                end if
            end if
            allocate(constraint_lb(mraw), constraint_ub(mraw))
            constraint_lb = -ieee_value(0.0_dp, ieee_positive_inf)
            constraint_ub = ieee_value(0.0_dp, ieee_positive_inf)
            if (present(constraint_lower)) constraint_lb = constraint_lower
            if (present(constraint_upper)) constraint_ub = constraint_upper
            do i = 1, mraw
                if (constraint_lb(i) > constraint_ub(i)) then
                    call finish_result(res, -1, .false., 'Invalid trust-constr constraint interval')
                    return
                end if
            end do

            constraint_eq_count = 0
            constraint_upper_count = 0
            constraint_lower_count = 0
            constraint_interval_count = 0
            do i = 1, mraw
                if (ieee_is_finite(constraint_lb(i)) .and. ieee_is_finite(constraint_ub(i))) then
                    if (constraint_lb(i) <= constraint_ub(i) .and. &
                            constraint_ub(i) <= constraint_lb(i)) then
                        constraint_eq_count = constraint_eq_count + 1
                    else
                        constraint_interval_count = constraint_interval_count + 1
                    end if
                else if (ieee_is_finite(constraint_ub(i))) then
                    constraint_upper_count = constraint_upper_count + 1
                else if (ieee_is_finite(constraint_lb(i))) then
                    constraint_lower_count = constraint_lower_count + 1
                end if
            end do
            allocate(constraint_eq_idx(constraint_eq_count), &
                constraint_upper_idx(constraint_upper_count), &
                constraint_lower_idx(constraint_lower_count), &
                constraint_interval_idx(constraint_interval_count))
            constraint_eq_count = 0
            constraint_upper_count = 0
            constraint_lower_count = 0
            constraint_interval_count = 0
            do i = 1, mraw
                if (ieee_is_finite(constraint_lb(i)) .and. ieee_is_finite(constraint_ub(i))) then
                    if (constraint_lb(i) <= constraint_ub(i) .and. &
                            constraint_ub(i) <= constraint_lb(i)) then
                        constraint_eq_count = constraint_eq_count + 1
                        constraint_eq_idx(constraint_eq_count) = i
                    else
                        constraint_interval_count = constraint_interval_count + 1
                        constraint_interval_idx(constraint_interval_count) = i
                    end if
                else if (ieee_is_finite(constraint_ub(i))) then
                    constraint_upper_count = constraint_upper_count + 1
                    constraint_upper_idx(constraint_upper_count) = i
                else if (ieee_is_finite(constraint_lb(i))) then
                    constraint_lower_count = constraint_lower_count + 1
                    constraint_lower_idx(constraint_lower_count) = i
                end if
            end do
            meq_use = constraint_eq_count
            mineq_use = constraint_upper_count + constraint_lower_count + &
                2 * constraint_interval_count
            morig = meq_use + mineq_use
        else
            morig = meq_use + mineq_use
            mraw = morig
            if (morig > 0 .and. .not. present(constraints)) then
                call finish_result(res, -1, .false., &
                    'trust-constr constraint counts require a constraint callback')
                return
            end if
            if (morig == 0 .and. (present(constraints) .or. present(constraint_jac))) then
                call finish_result(res, -1, .false., &
                    'Zero constraint count cannot have constraint callbacks')
                return
            end if
        end if

        if (present(constraint_keep_feasible) .and. .not. has_constraint_bounds) then
            call finish_result(res, -1, .false., &
                'constraint_keep_feasible requires constraint_lower/constraint_upper')
            return
        end if
        if (present(constraint_keep_feasible)) then
            if (size(constraint_keep_feasible) /= mraw) then
                call finish_result(res, -1, .false., 'Invalid constraint_keep_feasible size')
                return
            end if
        end if
        allocate(canonical_keep_feasible(morig))
        canonical_keep_feasible = .false.
        if (has_constraint_bounds .and. mraw > 0) then
            allocate(raw_constraint_keep_feasible(mraw))
            raw_constraint_keep_feasible = .false.
            if (present(constraint_keep_feasible)) then
                raw_constraint_keep_feasible = constraint_keep_feasible
            end if
            i = constraint_eq_count
            do j = 1, constraint_upper_count
                i = i + 1
                canonical_keep_feasible(i) = &
                    raw_constraint_keep_feasible(constraint_upper_idx(j))
            end do
            do j = 1, constraint_lower_count
                i = i + 1
                canonical_keep_feasible(i) = &
                    raw_constraint_keep_feasible(constraint_lower_idx(j))
            end do
            do j = 1, constraint_interval_count
                i = i + 1
                canonical_keep_feasible(i) = &
                    raw_constraint_keep_feasible(constraint_interval_idx(j))
            end do
            do j = 1, constraint_interval_count
                i = i + 1
                canonical_keep_feasible(i) = &
                    raw_constraint_keep_feasible(constraint_interval_idx(j))
            end do
        end if

        if (present(bounds_keep_feasible) .and. .not. has_bounds_block) then
            call finish_result(res, -1, .false., &
                'bounds_keep_feasible requires lower/upper bounds')
            return
        end if
        allocate(bound_keep_feasible_use(n))
        bound_keep_feasible_use = .false.
        if (present(bounds_keep_feasible)) then
            if (size(bounds_keep_feasible) /= n) then
                call finish_result(res, -1, .false., 'Invalid bounds_keep_feasible size')
                return
            end if
            bound_keep_feasible_use = bounds_keep_feasible
        end if

        if (mraw > 0) then
            allocate(raw_constraint_cache(mraw), raw_constraint_cache_x(n), &
                raw_constraint_jac_cache(mraw, n), raw_constraint_jac_cache_x(n))
            raw_constraint_cache = 0.0_dp
            raw_constraint_cache_x = 0.0_dp
            raw_constraint_jac_cache = 0.0_dp
            raw_constraint_jac_cache_x = 0.0_dp
        end if

        gtol_use = 1.0e-8_dp
        if (present(gtol)) gtol_use = gtol
        xtol_use = 1.0e-8_dp
        if (present(xtol)) xtol_use = xtol
        eps_use = default_fd_eps
        if (present(eps)) eps_use = eps
        objective_fd_step_use = eps_use
        if (objective_three_point .and. .not. present(eps)) then
            objective_fd_step_use = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
        end if
        maxiter_use = 1000
        if (present(maxiter)) maxiter_use = maxiter
        if (present(initial_trust_radius) .and. present(initial_tr_radius)) then
            call finish_result(res, -1, .false., &
                'Specify only one of initial_trust_radius and initial_tr_radius')
            return
        end if
        initial_radius_use = 1.0_dp
        if (present(initial_trust_radius)) initial_radius_use = initial_trust_radius
        if (present(initial_tr_radius)) initial_radius_use = initial_tr_radius
        initial_constr_penalty_use = 1.0_dp
        if (present(initial_constr_penalty)) initial_constr_penalty_use = initial_constr_penalty
        barrier_tol_use = 1.0e-8_dp
        if (present(barrier_tol)) barrier_tol_use = barrier_tol
        initial_barrier_parameter_use = 0.1_dp
        if (present(initial_barrier_parameter)) then
            initial_barrier_parameter_use = initial_barrier_parameter
        end if
        initial_barrier_tolerance_use = 0.1_dp
        if (present(initial_barrier_tolerance)) then
            initial_barrier_tolerance_use = initial_barrier_tolerance
        end if
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp .or. &
                .not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp .or. &
                .not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp .or. &
                .not. ieee_is_finite(initial_radius_use) .or. initial_radius_use <= 0.0_dp .or. &
                .not. ieee_is_finite(initial_constr_penalty_use) .or. &
                .not. ieee_is_finite(barrier_tol_use) .or. barrier_tol_use < 0.0_dp .or. &
                .not. ieee_is_finite(initial_barrier_parameter_use) .or. &
                initial_barrier_parameter_use <= 0.0_dp .or. &
                .not. ieee_is_finite(initial_barrier_tolerance_use) .or. &
                initial_barrier_tolerance_use <= 0.0_dp .or. &
                maxiter_use < 1) return

        allocate(lb(n), ub(n), nbd(n))
        call setup_box_bounds(n, lower, upper, lb, ub, nbd, invalid_bounds)
        if (invalid_bounds) then
            call finish_result(res, -1, .false., 'Invalid trust-constr bounds')
            return
        end if

        mfix = 0
        mbar = mineq_use
        do i = 1, n
            select case (nbd(i))
            case (1, 3)
                mbar = mbar + 1
            case (2)
                if (ub(i) <= lb(i)) then
                    mfix = mfix + 1
                else
                    mbar = mbar + 2
                end if
            end select
        end do

        allocate(barrier_keep_feasible(mbar))
        barrier_keep_feasible = .false.
        slack = 0
        do j = 1, mineq_use
            slack = slack + 1
            if (morig > 0) barrier_keep_feasible(slack) = &
                canonical_keep_feasible(meq_use + j)
        end do
        do i = 1, n
            select case (nbd(i))
            case (1)
                slack = slack + 1
                barrier_keep_feasible(slack) = bound_keep_feasible_use(i)
                if (bound_keep_feasible_use(i) .and. x0(i) < lb(i)) then
                    call finish_result(res, -1, .false., &
                        'Initial point violates a lower bound with keep_feasible')
                    return
                end if
            case (3)
                slack = slack + 1
                barrier_keep_feasible(slack) = bound_keep_feasible_use(i)
                if (bound_keep_feasible_use(i) .and. x0(i) > ub(i)) then
                    call finish_result(res, -1, .false., &
                        'Initial point violates an upper bound with keep_feasible')
                    return
                end if
            case (2)
                if (ub(i) > lb(i)) then
                    slack = slack + 1
                    barrier_keep_feasible(slack) = bound_keep_feasible_use(i)
                    slack = slack + 1
                    barrier_keep_feasible(slack) = bound_keep_feasible_use(i)
                    if (bound_keep_feasible_use(i) .and. &
                            (x0(i) < lb(i) .or. x0(i) > ub(i))) then
                        call finish_result(res, -1, .false., &
                            'Initial point violates an interval bound with keep_feasible')
                        return
                    end if
                end if
            end select
        end do

        if (present(constraint_hess) .and. morig <= 0) then
            call finish_result(res, -1, .false., &
                'constraint_hess requires nonlinear constraints')
            return
        end if

        if (mbar == 0 .and. mfix == 0) then
            if (meq_use == 0) then
                ! Without finite bounds or user constraints, trust-constr is
                ! still a valid unconstrained trust-region method.  The SQP
                ! kernel handles the zero-row Jacobian as an identity
                ! null-space projection, with no constraint evaluations.
                res = trust_constr_equality_sqp(fun, x0, jac=jac, hess=hess, &
                    meq=0, gtol=gtol_use, xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                    maxiter=maxiter_use, initial_trust_radius=initial_radius_use, &
                    initial_constr_penalty=initial_constr_penalty_use, &
                    callback=callback)
                return
            end if
            if (has_constraint_bounds) then
                if (present(constraint_jac)) then
                    if (present(constraint_hess)) then
                        res = trust_constr_equality_sqp(fun, x0, jac=jac, hess=hess, &
                            constraints=canonical_constraint_callback, &
                            constraint_jac=canonical_constraint_jacobian_callback, &
                            constraint_hess=canonical_constraint_hessian_callback, &
                            meq=meq_use, mineq=0, gtol=gtol_use, xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=maxiter_use, initial_trust_radius=initial_radius_use, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            callback=callback)
                    else
                        res = trust_constr_equality_sqp(fun, x0, jac=jac, hess=hess, &
                            constraints=canonical_constraint_callback, &
                            constraint_jac=canonical_constraint_jacobian_callback, &
                            meq=meq_use, mineq=0, gtol=gtol_use, xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=maxiter_use, initial_trust_radius=initial_radius_use, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            callback=callback)
                    end if
                else
                    if (present(constraint_hess)) then
                        res = trust_constr_equality_sqp(fun, x0, jac=jac, hess=hess, &
                            constraints=canonical_constraint_callback, &
                            constraint_hess=canonical_constraint_hessian_callback, &
                            meq=meq_use, mineq=0, gtol=gtol_use, xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=maxiter_use, initial_trust_radius=initial_radius_use, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            callback=callback)
                    else
                        res = trust_constr_equality_sqp(fun, x0, jac=jac, hess=hess, &
                            constraints=canonical_constraint_callback, meq=meq_use, mineq=0, &
                            gtol=gtol_use, xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=maxiter_use, &
                            initial_trust_radius=initial_radius_use, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            callback=callback)
                    end if
                end if
                inner = res
                call populate_constraint_diagnostics(res%x, inner)
            else
                res = trust_constr_equality_sqp(fun, x0, jac=jac, hess=hess, &
                    constraints=constraints, constraint_jac=constraint_jac, meq=meq_use, &
                    mineq=0, gtol=gtol_use, xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                    maxiter=maxiter_use, initial_trust_radius=initial_radius_use, &
                    initial_constr_penalty=initial_constr_penalty_use, &
                    constraint_hess=constraint_hess, &
                    callback=callback)
            end if
            if ((res%status == 1 .or. res%status == 2) .and. &
                    res%constr_violation > gtol_use) then
                call finish_result(res, 4, .false., "Constraint violation exceeds 'gtol'")
            end if
            return
        end if

        mwork = meq_use + mfix + mbar
        nwork = n + mbar
        ! The Gram-matrix pseudoinverse also handles overdetermined systems
        ! with redundant equality rows.  A full-rank constraint system is
        ! not required; infeasible residuals must still prevent convergence.

        allocate(z(nwork), z_fd(nwork), cbase(morig), c_fd(morig), x_fd(n), abase(morig, n), &
            cached_barrier_grad(nwork), cached_barrier_constr(mwork), cached_barrier_jac(mwork, nwork), &
            barrier_lagrangian_grad(n), strict_fd_nbd(nwork), strict_fd_lb(nwork), strict_fd_ub(nwork))
        strict_fd_nbd = 0
        strict_fd_lb = -ieee_value(0.0_dp, ieee_positive_inf)
        strict_fd_ub = ieee_value(0.0_dp, ieee_positive_inf)
        if (objective_three_point .and. mbar > 0) then
            ! The barrier logarithm is defined only for positive slacks.
            strict_fd_nbd(n + 1:nwork) = 1
            strict_fd_lb(n + 1:nwork) = 0.0_dp
        end if
        do i = 1, n
            if (bound_keep_feasible_use(i) .and. .not. &
                    (nbd(i) == 2 .and. ub(i) <= lb(i))) then
                strict_fd_nbd(i) = nbd(i)
                strict_fd_lb(i) = lb(i)
                strict_fd_ub(i) = ub(i)
            end if
        end do
        if (present(hess)) then
            allocate(cached_barrier_hess(nwork, nwork), stage_barrier_hess(nwork, nwork))
        end if
        z(1:n) = x0
        cbase = 0.0_dp
        if (morig > 0) call evaluate_original_constraints(x0, cbase)
        if (morig > 0) then
            if (.not. all(ieee_is_finite(cbase))) then
                call finish_result(res, -1, .false., 'Nonfinite initial constraints')
                return
            end if
        end if

        do j = 1, mineq_use
            if (barrier_keep_feasible(j) .and. cbase(meq_use + j) < 0.0_dp) then
                call finish_result(res, -1, .false., &
                    'Initial point violates a nonlinear constraint with keep_feasible')
                return
            end if
        end do

        slack = 0
        do j = 1, mineq_use
            slack = slack + 1
            if (barrier_keep_feasible(slack)) then
                z(n + slack) = max(cbase(meq_use + j), sqrt(epsilon(1.0_dp)))
            else
                z(n + slack) = max(1.0_dp, 1.5_dp * cbase(meq_use + j))
            end if
        end do
        do i = 1, n
            select case (nbd(i))
            case (1)
                slack = slack + 1
                if (barrier_keep_feasible(slack)) then
                    z(n + slack) = max(x0(i) - lb(i), sqrt(epsilon(1.0_dp)))
                else
                    z(n + slack) = max(1.0_dp, 1.5_dp * (x0(i) - lb(i)))
                end if
            case (3)
                slack = slack + 1
                if (barrier_keep_feasible(slack)) then
                    z(n + slack) = max(ub(i) - x0(i), sqrt(epsilon(1.0_dp)))
                else
                    z(n + slack) = max(1.0_dp, 1.5_dp * (ub(i) - x0(i)))
                end if
            case (2)
                if (ub(i) > lb(i)) then
                    slack = slack + 1
                    if (barrier_keep_feasible(slack)) then
                        z(n + slack) = max(x0(i) - lb(i), sqrt(epsilon(1.0_dp)))
                    else
                        z(n + slack) = max(1.0_dp, 1.5_dp * (x0(i) - lb(i)))
                    end if
                    slack = slack + 1
                    if (barrier_keep_feasible(slack)) then
                        z(n + slack) = max(ub(i) - x0(i), sqrt(epsilon(1.0_dp)))
                    else
                        z(n + slack) = max(1.0_dp, 1.5_dp * (ub(i) - x0(i)))
                    end if
                end if
            end select
        end do

        call assemble_barrier_constraints(z, cbase, cached_barrier_constr)
        abase = 0.0_dp
        if (morig > 0) then
            if (present(constraint_jac)) then
                call evaluate_original_constraint_jacobian(x0, abase)
            else
                do kcache = 1, n
                    fd_step = eps_use * max(1.0_dp, abs(x0(kcache)))
                    if (strict_fd_nbd(kcache) /= 0) then
                        call bounded_fd_step(x0(kcache), strict_fd_lb(kcache), &
                            strict_fd_ub(kcache), strict_fd_nbd(kcache), abs(fd_step), fd_step)
                    end if
                    x_fd = x0
                    x_fd(kcache) = x_fd(kcache) + fd_step
                    call evaluate_original_constraints(x_fd, c_fd)
                    abase(:, kcache) = (c_fd - cbase) / fd_step
                end do
            end if
        end if
        call assemble_barrier_jacobian(abase, cached_barrier_jac)
        if (present(constraint_hess)) then
            allocate(initial_constraint_hess(n, n), zero_constraint_multipliers(morig))
            zero_constraint_multipliers = 0.0_dp
            call evaluate_original_constraint_hessian(x0, zero_constraint_multipliers, &
                initial_constraint_hess)
            if (.not. all(ieee_is_finite(initial_constraint_hess))) then
                call finish_result(res, -1, .false., 'Nonfinite initial constraint Hessian')
                return
            end if
        end if

        barrier_parameter = initial_barrier_parameter_use
        barrier_tolerance = initial_barrier_tolerance_use
        stage_radius = initial_radius_use
        total_nit = 0
        total_cg_niter = 0
        cached_barrier_fun = barrier_objective(z)
        if (mbar > 0) then
            original_fun = cached_barrier_fun + barrier_parameter * sum(log(z(n + 1:nwork)))
        else
            original_fun = cached_barrier_fun
        end if
        barrier_optimality = huge(1.0_dp)
        barrier_lagrangian_grad = 0.0_dp
        total_nfev = 1
        if (present(jac)) then
            call barrier_gradient(z, cached_barrier_grad)
        else if (objective_three_point) then
            call three_point_gradient(barrier_objective, z, cached_barrier_fun, &
                objective_fd_step_use, cached_barrier_grad, kcache, &
                strict_fd_lb, strict_fd_ub, strict_fd_nbd, &
                objective_finite_diff_rel_step)
            total_nfev = total_nfev + kcache
        else
            do kcache = 1, nwork
                fd_step = eps_use * max(1.0_dp, abs(z(kcache)))
                if (present(objective_finite_diff_rel_step)) then
                    if (kcache <= n) fd_step = objective_relative_step( &
                        z(kcache), eps_use, objective_finite_diff_rel_step(kcache))
                end if
                if (strict_fd_nbd(kcache) /= 0) then
                    ! A user-supplied negative relative step remains negative
                    ! when it fits inside the retained variable bounds.
                    fd_request = fd_step
                    call bounded_signed_fd_step(z(kcache), strict_fd_lb(kcache), &
                        strict_fd_ub(kcache), strict_fd_nbd(kcache), fd_request, fd_step)
                end if
                z_fd = z
                z_fd(kcache) = z_fd(kcache) + fd_step
                fd_step = z_fd(kcache) - z(kcache)
                if (fd_step == 0.0_dp) then
                    cached_barrier_grad(kcache) = 0.0_dp
                    cycle
                end if
                fd_value = barrier_objective(z_fd)
                cached_barrier_grad(kcache) = (fd_value - cached_barrier_fun) / fd_step
                total_nfev = total_nfev + 1
            end do
        end if
        total_njev = 1
        total_nhev = 0
        if (present(hess)) then
            call barrier_hessian(z, cached_barrier_hess)
            total_nhev = 1
        end if
        final_success = .false.

        do stage = 1, max_barrier_stages
            remaining = maxiter_use - total_nit
            if (remaining <= 0) exit

            if (present(jac)) then
                if (present(hess)) then
                    if (present(constraint_hess)) then
                        inner = trust_constr_equality_sqp(barrier_objective, z, &
                            jac=barrier_gradient, hess=barrier_hessian, &
                            initial_objective_hess=cached_barrier_hess, &
                            final_objective_hess=stage_barrier_hess, &
                            constraints=barrier_constraints, &
                            constraint_jac=barrier_jacobian, &
                            meq=mwork, mineq=0, gtol=max(gtol_use, barrier_tolerance), &
                            xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=remaining, &
                            initial_trust_radius=stage_radius, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                            initial_constraints=cached_barrier_constr, &
                            initial_constraint_jac=cached_barrier_jac, &
                            fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                            fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                            enforce_tail=barrier_keep_feasible, &
                            constraint_hess=barrier_constraint_hessian, &
                            callback=barrier_iteration_callback)
                    else
                        inner = trust_constr_equality_sqp(barrier_objective, z, &
                            jac=barrier_gradient, hess=barrier_hessian, &
                            initial_objective_hess=cached_barrier_hess, &
                            final_objective_hess=stage_barrier_hess, &
                            constraints=barrier_constraints, &
                            constraint_jac=barrier_jacobian, &
                            meq=mwork, mineq=0, gtol=max(gtol_use, barrier_tolerance), &
                            xtol=xtol_use, eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=remaining, &
                            initial_trust_radius=stage_radius, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                            initial_constraints=cached_barrier_constr, &
                            initial_constraint_jac=cached_barrier_jac, &
                            fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                            fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                            enforce_tail=barrier_keep_feasible, &
                            approx_constraint_hess=morig > 0, &
                            callback=barrier_iteration_callback)
                    end if
                else if (present(constraint_hess)) then
                    inner = trust_constr_equality_sqp(barrier_objective, z, &
                        jac=barrier_gradient, constraints=barrier_constraints, &
                        constraint_jac=barrier_jacobian, meq=mwork, mineq=0, &
                        gtol=max(gtol_use, barrier_tolerance), xtol=xtol_use, &
                        eps=eps_use, objective_three_point=objective_three_point, &
                        objective_fd_step=objective_fd_step_use, &
                        objective_fd_rel_step=objective_finite_diff_rel_step, &
                        objective_hess_update=objective_hess_update, &
                        maxiter=remaining, &
                        initial_trust_radius=stage_radius, &
                        initial_constr_penalty=initial_constr_penalty_use, &
                        initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                        initial_constraints=cached_barrier_constr, &
                        initial_constraint_jac=cached_barrier_jac, &
                        fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                        fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                        enforce_tail=barrier_keep_feasible, &
                        constraint_hess=barrier_constraint_hessian, &
                        callback=barrier_iteration_callback)
                else
                    inner = trust_constr_equality_sqp(barrier_objective, z, &
                        jac=barrier_gradient, constraints=barrier_constraints, &
                        constraint_jac=barrier_jacobian, meq=mwork, mineq=0, &
                        gtol=max(gtol_use, barrier_tolerance), xtol=xtol_use, &
                        eps=eps_use, objective_three_point=objective_three_point, &
                        objective_fd_step=objective_fd_step_use, &
                        objective_fd_rel_step=objective_finite_diff_rel_step, &
                        objective_hess_update=objective_hess_update, &
                        maxiter=remaining, &
                        initial_trust_radius=stage_radius, &
                        initial_constr_penalty=initial_constr_penalty_use, &
                        initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                        initial_constraints=cached_barrier_constr, &
                        initial_constraint_jac=cached_barrier_jac, &
                        fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                        fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                        enforce_tail=barrier_keep_feasible, &
                        approx_constraint_hess=morig > 0, &
                        callback=barrier_iteration_callback)
                end if
            else
                if (present(hess)) then
                    if (present(constraint_hess)) then
                        inner = trust_constr_equality_sqp(barrier_objective, z, &
                            hess=barrier_hessian, &
                            initial_objective_hess=cached_barrier_hess, &
                            final_objective_hess=stage_barrier_hess, &
                            constraints=barrier_constraints, &
                            constraint_jac=barrier_jacobian, meq=mwork, mineq=0, &
                            gtol=max(gtol_use, barrier_tolerance), xtol=xtol_use, &
                            eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=remaining, &
                            initial_trust_radius=stage_radius, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                            initial_constraints=cached_barrier_constr, &
                            initial_constraint_jac=cached_barrier_jac, &
                            fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                            fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                            enforce_tail=barrier_keep_feasible, &
                            constraint_hess=barrier_constraint_hessian, &
                            callback=barrier_iteration_callback)
                    else
                        inner = trust_constr_equality_sqp(barrier_objective, z, &
                            hess=barrier_hessian, &
                            initial_objective_hess=cached_barrier_hess, &
                            final_objective_hess=stage_barrier_hess, &
                            constraints=barrier_constraints, &
                            constraint_jac=barrier_jacobian, meq=mwork, mineq=0, &
                            gtol=max(gtol_use, barrier_tolerance), xtol=xtol_use, &
                            eps=eps_use, &
                            objective_three_point=objective_three_point, &
                            objective_fd_step=objective_fd_step_use, &
                            objective_fd_rel_step=objective_finite_diff_rel_step, &
                            objective_hess_update=objective_hess_update, &
                            maxiter=remaining, &
                            initial_trust_radius=stage_radius, &
                            initial_constr_penalty=initial_constr_penalty_use, &
                            initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                            initial_constraints=cached_barrier_constr, &
                            initial_constraint_jac=cached_barrier_jac, &
                            fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                            fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                            enforce_tail=barrier_keep_feasible, &
                            approx_constraint_hess=morig > 0, &
                            callback=barrier_iteration_callback)
                    end if
                else if (present(constraint_hess)) then
                    inner = trust_constr_equality_sqp(barrier_objective, z, &
                        constraints=barrier_constraints, constraint_jac=barrier_jacobian, &
                        meq=mwork, mineq=0, gtol=max(gtol_use, barrier_tolerance), &
                        xtol=xtol_use, eps=eps_use, objective_three_point=objective_three_point, &
                        objective_fd_step=objective_fd_step_use, &
                        objective_fd_rel_step=objective_finite_diff_rel_step, &
                        objective_hess_update=objective_hess_update, &
                        maxiter=remaining, &
                        initial_trust_radius=stage_radius, &
                        initial_constr_penalty=initial_constr_penalty_use, &
                        initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                        initial_constraints=cached_barrier_constr, &
                        initial_constraint_jac=cached_barrier_jac, &
                        fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                        fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                        enforce_tail=barrier_keep_feasible, &
                        constraint_hess=barrier_constraint_hessian, &
                        callback=barrier_iteration_callback)
                else
                    inner = trust_constr_equality_sqp(barrier_objective, z, &
                        constraints=barrier_constraints, constraint_jac=barrier_jacobian, &
                        meq=mwork, mineq=0, gtol=max(gtol_use, barrier_tolerance), &
                        xtol=xtol_use, eps=eps_use, objective_three_point=objective_three_point, &
                        objective_fd_step=objective_fd_step_use, &
                        objective_fd_rel_step=objective_finite_diff_rel_step, &
                        objective_hess_update=objective_hess_update, &
                        maxiter=remaining, &
                        initial_trust_radius=stage_radius, &
                        initial_constr_penalty=initial_constr_penalty_use, &
                        initial_fun=cached_barrier_fun, initial_grad=cached_barrier_grad, &
                        initial_constraints=cached_barrier_constr, &
                        initial_constraint_jac=cached_barrier_jac, &
                        fd_lower=strict_fd_lb, fd_upper=strict_fd_ub, &
                        fd_nbd=strict_fd_nbd, use_soc=mbar > 0, positive_tail=mbar, &
                        enforce_tail=barrier_keep_feasible, &
                        approx_constraint_hess=morig > 0, &
                        callback=barrier_iteration_callback)
                end if
            end if

            total_nit = total_nit + inner%nit
            total_cg_niter = total_cg_niter + inner%cg_niter
            total_nfev = total_nfev + inner%nfev
            total_njev = total_njev + inner%njev
            total_nhev = total_nhev + inner%nhev
            if (allocated(inner%x)) z = inner%x
            if (allocated(inner%jac)) cached_barrier_grad = inner%jac
            if (allocated(inner%constr)) cached_barrier_constr = inner%constr
            if (allocated(inner%constr_jac)) cached_barrier_jac = inner%constr_jac
            if (present(hess)) cached_barrier_hess = stage_barrier_hess
            if (mbar > 0 .and. all(z(n + 1:nwork) > 0.0_dp) .and. ieee_is_finite(inner%fun)) then
                original_fun = inner%fun + barrier_parameter * sum(log(z(n + 1:nwork)))
            else if (mbar == 0) then
                original_fun = inner%fun
            else
                original_fun = fun(z(1:n))
                total_nfev = total_nfev + 1
            end if

            if (allocated(inner%constr)) then
                call barrier_original_constraint_violation(z, inner%constr, maxcv)
            else
                call original_constraint_violation(z(1:n), maxcv)
            end if
            optimality_ok = .false.
            barrier_optimality = huge(1.0_dp)
            if (allocated(inner%jac) .and. allocated(inner%v) .and. &
                    allocated(inner%constr_jac)) then
                call barrier_original_optimality(inner%jac, inner%v, inner%constr_jac, &
                    barrier_lagrangian_grad, barrier_optimality, optimality_ok)
            end if
            if (inner%status == 3) exit
            if (inner%status < 0 .or. .not. all(ieee_is_finite(z))) exit
            if (optimality_ok .and. barrier_optimality <= gtol_use .and. &
                    maxcv <= gtol_use .and. (morig > 0 .or. &
                    barrier_parameter <= max(gtol_use, barrier_tol_use))) then
                final_success = .true.
                inner%status = 1
                inner%message = '`gtol` termination condition is satisfied'
                exit
            end if
            if (barrier_parameter <= barrier_tol_use .and. inner%success .and. &
                    maxcv <= gtol_use) then
                final_success = .true.
                exit
            end if
            if (.not. inner%success .and. inner%status == 0) exit

            stage_radius = max(initial_radius_use, 5.0_dp * inner%tr_radius)
            barrier_parameter = barrier_parameter * barrier_decay
            barrier_tolerance = barrier_tolerance * barrier_decay
            if (mbar > 0) then
                cached_barrier_fun = original_fun - &
                    barrier_parameter * sum(log(z(n + 1:nwork)))
                cached_barrier_grad(n + 1:nwork) = -barrier_parameter / z(n + 1:nwork)
                if (present(hess)) then
                    do kcache = 1, mbar
                        cached_barrier_hess(n + kcache, n + kcache) = &
                            barrier_parameter / z(n + kcache)**2
                    end do
                end if
            else
                cached_barrier_fun = original_fun
            end if
        end do

        res%x = z(1:n)
        res%fun = original_fun
        if (allocated(inner%jac)) then
            res%jac = inner%jac(1:n)
        else
            allocate(res%jac(n))
            res%jac = 0.0_dp
        end if
        if (allocated(inner%constr)) then
            call barrier_original_constraint_violation(z, inner%constr, maxcv)
        else
            call original_constraint_violation(res%x, maxcv)
        end if
        res%maxcv = maxcv
        res%constr_violation = maxcv
        res%optimality = barrier_optimality
        res%lagrangian_grad = barrier_lagrangian_grad
        res%constr_penalty = inner%constr_penalty
        res%tr_radius = inner%tr_radius
        res%barrier_parameter = barrier_parameter
        res%barrier_tolerance = barrier_tolerance
        if (allocated(inner%v)) res%v = inner%v
        res%nit = total_nit
        res%nfev = total_nfev
        res%njev = total_njev
        res%nhev = total_nhev
        res%cg_niter = total_cg_niter
        res%cg_stop_cond = inner%cg_stop_cond
        res%method = 'tr_interior_point'
        call populate_constraint_diagnostics(z, inner)

        if (inner%status == 3) then
            call finish_result(res, 3, .false., trim(inner%message))
        else if (final_success) then
            call finish_result(res, inner%status, .true., trim(inner%message))
        else if (total_nit >= maxiter_use .or. inner%status == 0) then
            call finish_result(res, 0, .false., 'Maximum number of iterations reached')
        else if (inner%status < 0) then
            call finish_result(res, inner%status, .false., trim(inner%message))
        else
            call finish_result(res, 0, .false., 'Barrier sequence did not converge')
        end if

    contains

        subroutine barrier_iteration_callback(zv, nit, stop)
            real(dp), intent(in) :: zv(:) !! inner variables with appended barrier slacks
            integer, intent(in) :: nit !! iteration within the current barrier stage
            logical, intent(out) :: stop !! early termination requested by user

            stop = .false.
            if (present(callback)) then
                call invoke_iteration_callback(callback, zv(1:n), total_nit + nit, stop)
            end if
        end subroutine barrier_iteration_callback

        subroutine finite_difference_objective_hessian(xv, hmat)
            ! Two-point or nonuniform three-point differentiation of jac(x).
            ! Symmetrize the computed columns to control differencing noise.
            real(dp), intent(in) :: xv(:) !! physical parameter vector
            real(dp), intent(out) :: hmat(:, :) !! finite-difference Hessian

            integer :: col
            real(dp) :: step
            real(dp) :: h1
            real(dp) :: h2
            real(dp) :: xprobe(size(xv))
            real(dp) :: gbase(size(xv))
            real(dp) :: gfirst(size(xv))
            real(dp) :: gsecond(size(xv))

            hmat = 0.0_dp
            call invoke_gradient(jac, xv, gbase)
            hess_fd_njev = hess_fd_njev + 1
            do col = 1, size(xv)
                step = hess_fd_step * sign(1.0_dp, xv(col)) * max(1.0_dp, abs(xv(col)))
                if (present(objective_finite_diff_rel_step)) then
                    step = objective_relative_step(xv(col), hess_fd_step, &
                        objective_finite_diff_rel_step(col))
                end if
                if (hess_fd_three_point) then
                    call bounded_three_point_fd_steps(xv(col), hess_fd_lb(col), &
                        hess_fd_ub(col), hess_fd_nbd(col), abs(step), h1, h2)
                else
                    call bounded_signed_fd_step(xv(col), hess_fd_lb(col), hess_fd_ub(col), &
                        hess_fd_nbd(col), step, h1)
                    h2 = 0.0_dp
                end if
                if (h1 == 0.0_dp) cycle
                xprobe = xv
                xprobe(col) = xv(col) + h1
                if (hess_fd_nbd(col) == 1 .or. hess_fd_nbd(col) == 2) &
                    xprobe(col) = max(xprobe(col), hess_fd_lb(col))
                if (hess_fd_nbd(col) == 2 .or. hess_fd_nbd(col) == 3) &
                    xprobe(col) = min(xprobe(col), hess_fd_ub(col))
                h1 = xprobe(col) - xv(col)
                if (h1 == 0.0_dp) cycle
                call invoke_gradient(jac, xprobe, gfirst)
                hess_fd_njev = hess_fd_njev + 1
                if (hess_fd_three_point .and. h2 /= 0.0_dp) then
                    xprobe = xv
                    xprobe(col) = xv(col) + h2
                    if (hess_fd_nbd(col) == 1 .or. hess_fd_nbd(col) == 2) &
                        xprobe(col) = max(xprobe(col), hess_fd_lb(col))
                    if (hess_fd_nbd(col) == 2 .or. hess_fd_nbd(col) == 3) &
                        xprobe(col) = min(xprobe(col), hess_fd_ub(col))
                    h2 = xprobe(col) - xv(col)
                    if (h2 /= 0.0_dp .and. h2 /= h1) then
                        call invoke_gradient(jac, xprobe, gsecond)
                        hess_fd_njev = hess_fd_njev + 1
                        hmat(:, col) = (h2 * ((gfirst - gbase) / h1) - &
                            h1 * ((gsecond - gbase) / h2)) / (h2 - h1)
                    else
                        hmat(:, col) = (gfirst - gbase) / h1
                    end if
                else
                    hmat(:, col) = (gfirst - gbase) / h1
                end if
            end do
            hmat = 0.5_dp * (hmat + transpose(hmat))
        end subroutine finite_difference_objective_hessian

        subroutine assembled_objective_hessian(xv, hmat)
            real(dp), intent(in) :: xv(:) !! coordinates for the objective Hessian
            real(dp), intent(out) :: hmat(:, :) !! dense Hessian assembled from hessp

            integer :: col
            real(dp) :: direction(size(xv))
            real(dp) :: product(size(xv))

            hmat = 0.0_dp
            do col = 1, size(xv)
                direction = 0.0_dp
                direction(col) = 1.0_dp
                call invoke_hessian_product(hessp, xv, direction, product)
                hmat(:, col) = product
            end do
            ! The optimization problem defines a symmetric Hessian. Average
            ! tiny inconsistencies in separately evaluated products before
            ! passing it to the symmetric dense SQP solver.
            hmat = 0.5_dp * (hmat + transpose(hmat))
        end subroutine assembled_objective_hessian

        subroutine evaluate_original_constraints(xv, cv)
            real(dp), intent(in) :: xv(:) !! original variables at which to evaluate constraints
            real(dp), intent(out) :: cv(:) !! canonical constraint residuals

            real(dp) :: raw_cv(mraw)

            if (has_constraint_bounds) then
                call invoke_constraints(constraints, xv, raw_cv)
                user_constr_nfev = user_constr_nfev + 1
                raw_constraint_cache = raw_cv
                raw_constraint_cache_x = xv
                raw_constraint_cache_valid = .true.
                call canonicalize_constraint_values(raw_cv, cv)
            else
                call invoke_constraints(constraints, xv, cv)
                user_constr_nfev = user_constr_nfev + 1
                if (mraw > 0) then
                    raw_constraint_cache = cv
                    raw_constraint_cache_x = xv
                    raw_constraint_cache_valid = .true.
                end if
            end if
        end subroutine evaluate_original_constraints

        subroutine evaluate_original_constraint_jacobian(xv, av)
            real(dp), intent(in) :: xv(:) !! original variables at which to evaluate the Jacobian
            real(dp), intent(out) :: av(:, :) !! canonical constraint Jacobian

            real(dp) :: raw_av(mraw, n)

            if (has_constraint_bounds) then
                call invoke_constraint_jacobian(constraint_jac, xv, raw_av)
                user_constr_njev = user_constr_njev + 1
                raw_constraint_jac_cache = raw_av
                raw_constraint_jac_cache_x = xv
                raw_constraint_jac_cache_valid = .true.
                call canonicalize_constraint_jacobian(raw_av, av)
            else
                call invoke_constraint_jacobian(constraint_jac, xv, av)
                user_constr_njev = user_constr_njev + 1
                if (mraw > 0) then
                    raw_constraint_jac_cache = av
                    raw_constraint_jac_cache_x = xv
                    raw_constraint_jac_cache_valid = .true.
                end if
            end if
        end subroutine evaluate_original_constraint_jacobian

        subroutine evaluate_original_constraint_hessian(xv, vv, hv)
            real(dp), intent(in) :: xv(:) !! original variables at which to evaluate the Hessian
            real(dp), intent(in) :: vv(:) !! canonical constraint multipliers
            real(dp), intent(out) :: hv(:, :) !! weighted original constraint Hessian

            real(dp) :: raw_v(mraw)

            if (has_constraint_bounds) then
                call canonicalize_constraint_multipliers(vv, raw_v)
                call invoke_constraint_hessian(constraint_hess, xv, raw_v, hv)
            else
                call invoke_constraint_hessian(constraint_hess, xv, vv, hv)
            end if
            user_constr_nhev = user_constr_nhev + 1
        end subroutine evaluate_original_constraint_hessian

        subroutine canonicalize_constraint_values(raw_cv, cv)
            real(dp), intent(in) :: raw_cv(:) !! raw user constraint values
            real(dp), intent(out) :: cv(:) !! equalities followed by >= 0 inequalities

            integer :: k
            integer :: row

            cv = 0.0_dp
            row = 0
            do k = 1, constraint_eq_count
                row = row + 1
                cv(row) = raw_cv(constraint_eq_idx(k)) - constraint_lb(constraint_eq_idx(k))
            end do
            do k = 1, constraint_upper_count
                row = row + 1
                cv(row) = constraint_ub(constraint_upper_idx(k)) - raw_cv(constraint_upper_idx(k))
            end do
            do k = 1, constraint_lower_count
                row = row + 1
                cv(row) = raw_cv(constraint_lower_idx(k)) - constraint_lb(constraint_lower_idx(k))
            end do
            do k = 1, constraint_interval_count
                row = row + 1
                cv(row) = constraint_ub(constraint_interval_idx(k)) - raw_cv(constraint_interval_idx(k))
            end do
            do k = 1, constraint_interval_count
                row = row + 1
                cv(row) = raw_cv(constraint_interval_idx(k)) - constraint_lb(constraint_interval_idx(k))
            end do
        end subroutine canonicalize_constraint_values

        subroutine canonicalize_constraint_jacobian(raw_av, av)
            real(dp), intent(in) :: raw_av(:, :) !! raw user constraint Jacobian
            real(dp), intent(out) :: av(:, :) !! canonical constraint Jacobian

            integer :: k
            integer :: row

            av = 0.0_dp
            row = 0
            do k = 1, constraint_eq_count
                row = row + 1
                av(row, :) = raw_av(constraint_eq_idx(k), :)
            end do
            do k = 1, constraint_upper_count
                row = row + 1
                av(row, :) = -raw_av(constraint_upper_idx(k), :)
            end do
            do k = 1, constraint_lower_count
                row = row + 1
                av(row, :) = raw_av(constraint_lower_idx(k), :)
            end do
            do k = 1, constraint_interval_count
                row = row + 1
                av(row, :) = -raw_av(constraint_interval_idx(k), :)
            end do
            do k = 1, constraint_interval_count
                row = row + 1
                av(row, :) = raw_av(constraint_interval_idx(k), :)
            end do
        end subroutine canonicalize_constraint_jacobian

        subroutine canonicalize_constraint_multipliers(vv, raw_v)
            real(dp), intent(in) :: vv(:) !! canonical equality/inequality multipliers
            real(dp), intent(out) :: raw_v(:) !! multipliers for the raw user callback

            integer :: k
            integer :: row

            raw_v = 0.0_dp
            row = 0
            do k = 1, constraint_eq_count
                row = row + 1
                raw_v(constraint_eq_idx(k)) = vv(row)
            end do
            do k = 1, constraint_upper_count
                row = row + 1
                raw_v(constraint_upper_idx(k)) = -vv(row)
            end do
            do k = 1, constraint_lower_count
                row = row + 1
                raw_v(constraint_lower_idx(k)) = vv(row)
            end do
            do k = 1, constraint_interval_count
                row = row + 1
                raw_v(constraint_interval_idx(k)) = -vv(row)
            end do
            do k = 1, constraint_interval_count
                row = row + 1
                raw_v(constraint_interval_idx(k)) = raw_v(constraint_interval_idx(k)) + vv(row)
            end do
        end subroutine canonicalize_constraint_multipliers

        subroutine canonical_constraint_callback(xv, cv)
            real(dp), intent(in) :: xv(:) !! original variables
            real(dp), intent(out) :: cv(:) !! canonical equality residuals

            call evaluate_original_constraints(xv, cv)
        end subroutine canonical_constraint_callback

        subroutine canonical_constraint_jacobian_callback(xv, av)
            real(dp), intent(in) :: xv(:) !! original variables
            real(dp), intent(out) :: av(:, :) !! canonical constraint Jacobian

            call evaluate_original_constraint_jacobian(xv, av)
        end subroutine canonical_constraint_jacobian_callback

        subroutine canonical_constraint_hessian_callback(xv, vv, hv)
            real(dp), intent(in) :: xv(:) !! original variables
            real(dp), intent(in) :: vv(:) !! canonical multipliers
            real(dp), intent(out) :: hv(:, :) !! weighted raw constraint Hessian

            call evaluate_original_constraint_hessian(xv, vv, hv)
        end subroutine canonical_constraint_hessian_callback

        subroutine populate_constraint_diagnostics(zv, inner_result)
            real(dp), intent(in) :: zv(:) !! final variables followed by barrier slacks
            type(multivariate_minimize_result), intent(in) :: inner_result !! final barrier SQP state

            integer :: block
            integer :: diag_row
            integer :: fixed_rows
            integer :: k
            integer :: nblocks
            integer :: out_row
            integer :: row
            integer :: slack_index
            real(dp) :: canonical_values(morig)
            real(dp) :: canonical_jac(morig, n)

            if (allocated(res%constr_nfev)) deallocate(res%constr_nfev)
            if (allocated(res%constr_njev)) deallocate(res%constr_njev)
            if (allocated(res%constr_nhev)) deallocate(res%constr_nhev)
            if (allocated(res%constr)) deallocate(res%constr)
            if (allocated(res%constr_jac)) deallocate(res%constr_jac)

            nblocks = 0
            if (mraw > 0) nblocks = nblocks + 1
            if (has_bounds_block) nblocks = nblocks + 1
            allocate(res%constr_nfev(nblocks), res%constr_njev(nblocks), res%constr_nhev(nblocks))
            if (nblocks > 0) then
                res%constr_nfev = 0
                res%constr_njev = 0
                res%constr_nhev = 0
            end if

            block = 0
            if (mraw > 0) then
                block = block + 1
                res%constr_nfev(block) = user_constr_nfev
                res%constr_njev(block) = user_constr_njev
                res%constr_nhev(block) = user_constr_nhev
            end if
            if (has_bounds_block) block = block + 1

            out_row = merge(mraw, morig, has_constraint_bounds)
            allocate(res%constr(out_row + merge(n, 0, has_bounds_block)))
            allocate(res%constr_jac(out_row + merge(n, 0, has_bounds_block), n))
            res%constr = 0.0_dp
            res%constr_jac = 0.0_dp

            canonical_values = 0.0_dp
            canonical_jac = 0.0_dp
            if (morig > 0 .and. allocated(inner_result%constr) .and. &
                    allocated(inner_result%constr_jac)) then
                fixed_rows = 0
                do k = 1, n
                    if (nbd(k) == 2 .and. ub(k) <= lb(k)) fixed_rows = fixed_rows + 1
                end do
                if (meq_use > 0) then
                    canonical_values(1:meq_use) = inner_result%constr(1:meq_use)
                    canonical_jac(1:meq_use, :) = inner_result%constr_jac(1:meq_use, 1:n)
                end if
                diag_row = meq_use + fixed_rows
                slack_index = 0
                do k = 1, mineq_use
                    diag_row = diag_row + 1
                    slack_index = slack_index + 1
                    canonical_values(meq_use + k) = &
                        inner_result%constr(diag_row) + zv(n + slack_index)
                    canonical_jac(meq_use + k, :) = inner_result%constr_jac(diag_row, 1:n)
                end do
            end if

            if (has_constraint_bounds) then
                if (raw_constraint_cache_valid) res%constr(1:mraw) = raw_constraint_cache
                if (raw_constraint_jac_cache_valid) then
                    res%constr_jac(1:mraw, :) = raw_constraint_jac_cache
                end if
                row = 0
                do k = 1, constraint_eq_count
                    row = row + 1
                    res%constr(constraint_eq_idx(k)) = &
                        canonical_values(row) + constraint_lb(constraint_eq_idx(k))
                    res%constr_jac(constraint_eq_idx(k), :) = canonical_jac(row, :)
                end do
                do k = 1, constraint_upper_count
                    row = row + 1
                    res%constr(constraint_upper_idx(k)) = &
                        constraint_ub(constraint_upper_idx(k)) - canonical_values(row)
                    res%constr_jac(constraint_upper_idx(k), :) = -canonical_jac(row, :)
                end do
                do k = 1, constraint_lower_count
                    row = row + 1
                    res%constr(constraint_lower_idx(k)) = &
                        constraint_lb(constraint_lower_idx(k)) + canonical_values(row)
                    res%constr_jac(constraint_lower_idx(k), :) = canonical_jac(row, :)
                end do
                do k = 1, constraint_interval_count
                    row = row + 1
                    res%constr(constraint_interval_idx(k)) = &
                        constraint_ub(constraint_interval_idx(k)) - canonical_values(row)
                    res%constr_jac(constraint_interval_idx(k), :) = -canonical_jac(row, :)
                end do
                row = row + constraint_interval_count
            else if (morig > 0) then
                res%constr(1:morig) = canonical_values
                res%constr_jac(1:morig, :) = canonical_jac
            end if

            if (has_bounds_block) then
                res%constr(out_row + 1:out_row + n) = zv(1:n)
                do k = 1, n
                    res%constr_jac(out_row + k, k) = 1.0_dp
                end do
            end if
        end subroutine populate_constraint_diagnostics

        function barrier_objective(zv) result(fv)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp) :: fv

            ! Match SciPy's barrier trial guard: reject trial points outside
            ! strict variable bounds without evaluating the user objective.
            if (violates_strict_bounds(zv(1:n))) then
                fv = ieee_value(0.0_dp, ieee_positive_inf)
                return
            end if
            if (mbar > 0) then
                if (any(zv(n + 1:nwork) <= 0.0_dp)) then
                    fv = ieee_value(0.0_dp, ieee_positive_inf)
                    return
                end if
                fv = fun(zv(1:n)) - barrier_parameter * sum(log(zv(n + 1:nwork)))
            else
                fv = fun(zv(1:n))
            end if
        end function barrier_objective

        subroutine barrier_gradient(zv, gv)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp), intent(out) :: gv(:) !! barrier objective gradient

            call invoke_gradient(jac, zv(1:n), gv(1:n))
            if (mbar > 0) gv(n + 1:nwork) = -barrier_parameter / zv(n + 1:nwork)
        end subroutine barrier_gradient

        subroutine barrier_hessian(zv, hmat)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp), intent(out) :: hmat(:, :) !! barrier objective Hessian

            hmat = 0.0_dp
            call invoke_hessian(hess, zv(1:n), hmat(1:n, 1:n))
            do i = 1, mbar
                hmat(n + i, n + i) = barrier_parameter / zv(n + i)**2
            end do
        end subroutine barrier_hessian

        subroutine barrier_constraint_hessian(zv, vv, hmat)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp), intent(in) :: vv(:) !! equality-form barrier multipliers
            real(dp), intent(out) :: hmat(:, :) !! Hessian of weighted barrier constraints

            integer :: fixed_row
            integer :: k
            real(dp) :: hbase(n, n)
            real(dp) :: vorig(morig)

            vorig = 0.0_dp
            fixed_row = meq_use
            if (meq_use > 0) vorig(1:meq_use) = vv(1:meq_use)
            do k = 1, n
                if (nbd(k) == 2 .and. ub(k) <= lb(k)) fixed_row = fixed_row + 1
            end do
            if (mineq_use > 0) then
                vorig(meq_use + 1:morig) = vv(fixed_row + 1:fixed_row + mineq_use)
            end if

            call evaluate_original_constraint_hessian(zv(1:n), vorig, hbase)
            hmat = 0.0_dp
            hmat(1:n, 1:n) = hbase
        end subroutine barrier_constraint_hessian

        subroutine barrier_constraints(zv, cv)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp), intent(out) :: cv(:) !! equality-form barrier constraints

            real(dp) :: corig(morig)

            ! Trial points outside keep_feasible bounds must never reach the
            ! user's nonlinear constraint callbacks (SciPy's strict bounds).
            ! The paired barrier objective is +Inf, so this trial is rejected.
            if (violates_strict_bounds(zv(1:n))) then
                cv = 0.0_dp
                return
            end if
            corig = 0.0_dp
            if (morig > 0) call evaluate_original_constraints(zv(1:n), corig)
            call assemble_barrier_constraints(zv, corig, cv)
        end subroutine barrier_constraints

        logical function violates_strict_bounds(xv) result(violated)
            real(dp), intent(in) :: xv(:) !! original optimization variables

            integer :: k

            violated = .false.
            do k = 1, n
                if (.not. bound_keep_feasible_use(k)) cycle
                ! Fixed bounds are equalities, not strict inequalities.
                if (nbd(k) == 2 .and. lb(k) >= ub(k)) cycle
                if (nbd(k) == 1 .or. nbd(k) == 2) then
                    if (xv(k) < lb(k)) violated = .true.
                end if
                if (nbd(k) == 2 .or. nbd(k) == 3) then
                    if (xv(k) > ub(k)) violated = .true.
                end if
                if (violated) return
            end do
        end function violates_strict_bounds

        subroutine assemble_barrier_constraints(zv, corig, cv)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp), intent(in) :: corig(:) !! already evaluated original constraints
            real(dp), intent(out) :: cv(:) !! equality-form barrier constraints

            integer :: i_local
            integer :: j_local
            integer :: row_local
            integer :: slack_local

            cv = 0.0_dp
            row_local = 0
            if (meq_use > 0) then
                cv(1:meq_use) = corig(1:meq_use)
                row_local = meq_use
            end if
            do i_local = 1, n
                if (nbd(i_local) == 2 .and. ub(i_local) <= lb(i_local)) then
                    row_local = row_local + 1
                    cv(row_local) = zv(i_local) - lb(i_local)
                end if
            end do
            slack_local = 0
            do j_local = 1, mineq_use
                row_local = row_local + 1
                slack_local = slack_local + 1
                cv(row_local) = corig(meq_use + j_local) - zv(n + slack_local)
            end do
            do i_local = 1, n
                select case (nbd(i_local))
                case (1)
                    row_local = row_local + 1
                    slack_local = slack_local + 1
                    cv(row_local) = zv(i_local) - lb(i_local) - zv(n + slack_local)
                case (3)
                    row_local = row_local + 1
                    slack_local = slack_local + 1
                    cv(row_local) = ub(i_local) - zv(i_local) - zv(n + slack_local)
                case (2)
                    if (ub(i_local) > lb(i_local)) then
                        row_local = row_local + 1
                        slack_local = slack_local + 1
                        cv(row_local) = zv(i_local) - lb(i_local) - zv(n + slack_local)
                        row_local = row_local + 1
                        slack_local = slack_local + 1
                        cv(row_local) = ub(i_local) - zv(i_local) - zv(n + slack_local)
                    end if
                end select
            end do
        end subroutine assemble_barrier_constraints

        subroutine barrier_jacobian(zv, amat)
            real(dp), intent(in) :: zv(:) !! variables followed by positive slacks
            real(dp), intent(out) :: amat(:, :) !! equality-form barrier Jacobian

            integer :: k
            real(dp) :: h
            real(dp) :: corig(morig)
            real(dp) :: cp(morig)
            real(dp) :: xp(n)

            abase = 0.0_dp
            corig = 0.0_dp
            cp = 0.0_dp
            if (morig > 0) then
                if (present(constraint_jac)) then
                    call evaluate_original_constraint_jacobian(zv(1:n), abase)
                else
                    call evaluate_original_constraints(zv(1:n), corig)
                    do k = 1, n
                        h = eps_use * max(1.0_dp, abs(zv(k)))
                        if (strict_fd_nbd(k) /= 0) then
                            call bounded_fd_step(zv(k), strict_fd_lb(k), strict_fd_ub(k), &
                                strict_fd_nbd(k), abs(h), h)
                        end if
                        xp = zv(1:n)
                        xp(k) = xp(k) + h
                        call evaluate_original_constraints(xp, cp)
                        abase(:, k) = (cp - corig) / h
                    end do
                end if
            end if
            call assemble_barrier_jacobian(abase, amat)
        end subroutine barrier_jacobian

        subroutine assemble_barrier_jacobian(abase_in, amat)
            real(dp), intent(in) :: abase_in(:, :) !! already evaluated original constraint Jacobian
            real(dp), intent(out) :: amat(:, :) !! equality-form barrier Jacobian

            integer :: i_local
            integer :: j_local
            integer :: row_local
            integer :: slack_local

            amat = 0.0_dp
            row_local = 0
            if (meq_use > 0) then
                amat(1:meq_use, 1:n) = abase_in(1:meq_use, :)
                row_local = meq_use
            end if
            do i_local = 1, n
                if (nbd(i_local) == 2 .and. ub(i_local) <= lb(i_local)) then
                    row_local = row_local + 1
                    amat(row_local, i_local) = 1.0_dp
                end if
            end do
            slack_local = 0
            do j_local = 1, mineq_use
                row_local = row_local + 1
                slack_local = slack_local + 1
                amat(row_local, 1:n) = abase_in(meq_use + j_local, :)
                amat(row_local, n + slack_local) = -1.0_dp
            end do
            do i_local = 1, n
                select case (nbd(i_local))
                case (1)
                    row_local = row_local + 1
                    slack_local = slack_local + 1
                    amat(row_local, i_local) = 1.0_dp
                    amat(row_local, n + slack_local) = -1.0_dp
                case (3)
                    row_local = row_local + 1
                    slack_local = slack_local + 1
                    amat(row_local, i_local) = -1.0_dp
                    amat(row_local, n + slack_local) = -1.0_dp
                case (2)
                    if (ub(i_local) > lb(i_local)) then
                        row_local = row_local + 1
                        slack_local = slack_local + 1
                        amat(row_local, i_local) = 1.0_dp
                        amat(row_local, n + slack_local) = -1.0_dp
                        row_local = row_local + 1
                        slack_local = slack_local + 1
                        amat(row_local, i_local) = -1.0_dp
                        amat(row_local, n + slack_local) = -1.0_dp
                    end if
                end select
            end do
        end subroutine assemble_barrier_jacobian


        subroutine barrier_original_optimality(gbar, vv, abar, lagrangian_grad, optimality, ok)
            real(dp), intent(in) :: gbar(:) !! physical-coordinate barrier objective gradient
            real(dp), intent(in) :: vv(:) !! current equality multipliers from the inner SQP
            real(dp), intent(in) :: abar(:, :) !! cached equality-form barrier Jacobian
            real(dp), intent(out) :: lagrangian_grad(:) !! original-variable Lagrangian gradient
            real(dp), intent(out) :: optimality !! original-variable Lagrangian infinity norm
            logical, intent(out) :: ok !! true when the multiplier state is finite

            real(dp), allocatable :: lagw(:)

            ok = size(gbar) >= n .and. size(vv) == mwork .and. &
                size(abar, 1) == mwork .and. size(abar, 2) >= n
            if (.not. ok) then
                lagrangian_grad = 0.0_dp
                optimality = huge(1.0_dp)
                return
            end if
            allocate(lagw(n))
            lagw = gbar(1:n) + matmul(transpose(abar(:, 1:n)), vv)
            lagrangian_grad = lagw
            optimality = maxval(abs(lagw))
            ok = ieee_is_finite(optimality) .and. all(ieee_is_finite(vv)) .and. &
                all(ieee_is_finite(abar(:, 1:n)))
        end subroutine barrier_original_optimality

        subroutine barrier_original_constraint_violation(zv, barrier_c, violation)
            real(dp), intent(in) :: zv(:) !! current barrier variables
            real(dp), intent(in) :: barrier_c(:) !! cached equality-form barrier constraints
            real(dp), intent(out) :: violation !! maximum original constraint violation

            integer :: diag_row
            integer :: fixed_rows
            integer :: i_local
            integer :: k
            integer :: slack_index
            real(dp) :: corig(morig)

            violation = 0.0_dp
            corig = 0.0_dp
            if (morig > 0) then
                if (size(barrier_c) /= mwork) then
                    violation = huge(1.0_dp)
                    return
                end if
                if (meq_use > 0) corig(1:meq_use) = barrier_c(1:meq_use)
                fixed_rows = 0
                do i_local = 1, n
                    if (nbd(i_local) == 2 .and. ub(i_local) <= lb(i_local)) then
                        fixed_rows = fixed_rows + 1
                    end if
                end do
                diag_row = meq_use + fixed_rows
                slack_index = 0
                do k = 1, mineq_use
                    diag_row = diag_row + 1
                    slack_index = slack_index + 1
                    corig(meq_use + k) = barrier_c(diag_row) + zv(n + slack_index)
                end do
                if (.not. all(ieee_is_finite(corig))) then
                    violation = huge(1.0_dp)
                    return
                end if
                if (meq_use > 0) violation = max(violation, maxval(abs(corig(1:meq_use))))
                if (mineq_use > 0) violation = max(violation, &
                    maxval(max(0.0_dp, -corig(meq_use + 1:morig))))
            end if
            do i_local = 1, n
                select case (nbd(i_local))
                case (1)
                    violation = max(violation, max(0.0_dp, lb(i_local) - zv(i_local)))
                case (3)
                    violation = max(violation, max(0.0_dp, zv(i_local) - ub(i_local)))
                case (2)
                    violation = max(violation, max(0.0_dp, lb(i_local) - zv(i_local)))
                    violation = max(violation, max(0.0_dp, zv(i_local) - ub(i_local)))
                end select
            end do
        end subroutine barrier_original_constraint_violation

        subroutine original_constraint_violation(xv, violation)
            real(dp), intent(in) :: xv(:) !! original variables
            real(dp), intent(out) :: violation !! maximum original constraint violation

            real(dp) :: corig(morig)

            violation = 0.0_dp
            corig = 0.0_dp
            if (morig > 0) then
                call evaluate_original_constraints(xv, corig)
                if (.not. all(ieee_is_finite(corig))) then
                    violation = huge(1.0_dp)
                    return
                end if
                if (meq_use > 0) violation = max(violation, maxval(abs(corig(1:meq_use))))
                if (mineq_use > 0) violation = max(violation, &
                    maxval(max(0.0_dp, -corig(meq_use + 1:morig))))
            end if
            do i = 1, n
                select case (nbd(i))
                case (1)
                    violation = max(violation, max(0.0_dp, lb(i) - xv(i)))
                case (2)
                    violation = max(violation, max(0.0_dp, lb(i) - xv(i)))
                    violation = max(violation, max(0.0_dp, xv(i) - ub(i)))
                case (3)
                    violation = max(violation, max(0.0_dp, xv(i) - ub(i)))
                end select
            end do
        end subroutine original_constraint_violation

    end function trust_constr

    ! Assemble independent nonlinear constraint blocks into the existing
    ! vector-bound trust-constr machinery. Per-block raw values/Jacobians retain
    ! their user-specified order, while the solver performs the canonical sorting.
    recursive function trust_constr_objects(fun, x0, objects, jac, hess, lower, upper, &
            gtol, xtol, eps, maxiter, initial_trust_radius, initial_tr_radius, &
            initial_constr_penalty, barrier_tol, initial_barrier_parameter, &
            initial_barrier_tolerance, bounds_keep_feasible, objective_fd_scheme, &
            objective_finite_diff_rel_step, objective_hess_update, callback) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        type(nonlinear_constraint), intent(in) :: objects(:) !! independently bounded linear/nonlinear constraints
        procedure(multivariate_gradient_function), optional :: jac !! analytic objective gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic objective Hessian
        real(dp), intent(in), optional :: lower(:) !! variable lower bounds
        real(dp), intent(in), optional :: upper(:) !! variable upper bounds
        real(dp), intent(in), optional :: gtol !! optimality tolerance
        real(dp), intent(in), optional :: xtol !! trust radius tolerance
        real(dp), intent(in), optional :: eps !! numerical derivative step
        integer, intent(in), optional :: maxiter !! maximum total iterations
        real(dp), intent(in), optional :: initial_trust_radius !! legacy initial radius
        real(dp), intent(in), optional :: initial_tr_radius !! SciPy initial radius
        real(dp), intent(in), optional :: initial_constr_penalty !! initial merit penalty
        real(dp), intent(in), optional :: barrier_tol !! barrier termination threshold
        real(dp), intent(in), optional :: initial_barrier_parameter !! initial barrier value
        real(dp), intent(in), optional :: initial_barrier_tolerance !! initial barrier tolerance
        logical, intent(in), optional :: bounds_keep_feasible(:) !! retained variable bounds
        character(len=*), intent(in), optional :: objective_fd_scheme !! objective numerical gradient scheme
        real(dp), intent(in), optional :: objective_finite_diff_rel_step(:) !! per-variable objective relative steps
        character(len=*), intent(in), optional :: objective_hess_update !! objective Hessian BFGS or SR1 update
        procedure(trust_constr_iteration_callback), optional :: callback !! iteration/termination callback
        type(multivariate_minimize_result) :: res

        integer :: i
        integer :: k
        integer :: m
        integer :: n
        integer :: nb
        integer, allocatable :: first(:)
        integer, allocatable :: last(:)
        integer, allocatable :: block_nfev(:)
        integer, allocatable :: block_njev(:)
        integer, allocatable :: block_nhev(:)
        integer, allocatable :: fd_nbd(:)
        logical :: any_hess
        logical :: all_hess
        logical :: any_update
        logical :: invalid
        logical, allocatable :: cache_valid(:)
        logical, allocatable :: update_valid(:)
        logical, allocatable :: update_jac_ready(:)
        logical, allocatable :: update_first(:)
        real(dp) :: eps_use
        real(dp) :: infinity
        real(dp), allocatable :: lo(:)
        real(dp), allocatable :: hi(:)
        real(dp), allocatable :: bound_lo(:)
        real(dp), allocatable :: bound_hi(:)
        real(dp), allocatable :: cache_x(:, :)
        real(dp), allocatable :: cache_c(:)
        real(dp), allocatable :: update_x(:, :)
        real(dp), allocatable :: update_jac_x(:, :)
        real(dp), allocatable :: update_jac_work(:, :)
        real(dp), allocatable :: update_jac(:, :)
        real(dp), allocatable :: update_hess(:, :, :)
        logical, allocatable :: keep(:)

        n = size(x0)
        nb = size(objects)
        call init_result(res, 'trust-constr', n)
        if (n < 1 .or. .not. all(ieee_is_finite(x0))) return
        if (nb == 0) then
            ! SciPy treats constraints=[] as an unconstrained problem (or a
            ! bounds-only problem when variable bounds are supplied).
            res = trust_constr(fun, x0, jac=jac, hess=hess, lower=lower, upper=upper, &
                gtol=gtol, xtol=xtol, eps=eps, maxiter=maxiter, &
                initial_trust_radius=initial_trust_radius, initial_tr_radius=initial_tr_radius, &
                initial_constr_penalty=initial_constr_penalty, barrier_tol=barrier_tol, &
                initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                bounds_keep_feasible=bounds_keep_feasible, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, &
                objective_hess_update=objective_hess_update, &
                callback=callback)
            return
        end if
        eps_use = default_fd_eps
        if (present(eps)) eps_use = eps
        if (.not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp) then
            call finish_result(res, -1, .false., 'Invalid constraint finite-difference step')
            return
        end if
        infinity = ieee_value(0.0_dp, ieee_positive_inf)
        allocate(first(nb), last(nb), block_nfev(nb), block_njev(nb), block_nhev(nb))
        block_nfev = 0
        block_njev = 0
        block_nhev = 0
        m = 0
        any_hess = .false.
        all_hess = .true.
        any_update = .false.
        do i = 1, nb
            if (allocated(objects(i)%matrix)) then
                if (associated(objects(i)%fun) .or. associated(objects(i)%jac) .or. &
                        associated(objects(i)%hess) .or. associated(objects(i)%hessp)) then
                    call finish_result(res, -1, .false., 'Linear constraint block cannot have callbacks')
                    return
                end if
                if (size(objects(i)%matrix, 2) /= n .or. &
                        .not. all(ieee_is_finite(objects(i)%matrix))) then
                    call finish_result(res, -1, .false., 'Invalid linear constraint matrix')
                    return
                end if
                if (allocated(objects(i)%finite_diff_rel_step) .or. &
                        allocated(objects(i)%finite_diff_jac_sparsity) .or. &
                        objects(i)%fd_scheme /= '2-point' .or. &
                        objects(i)%hess_fd_scheme /= '' .or. objects(i)%hess_update /= '') then
                    call finish_result(res, -1, .false., 'Linear constraint block cannot use derivative strategies')
                    return
                end if
            else if (.not. associated(objects(i)%fun)) then
                call finish_result(res, -1, .false., 'Unassociated constraint block callback')
                return
            end if
            if (allocated(objects(i)%lower)) then
                k = size(objects(i)%lower)
            else if (allocated(objects(i)%upper)) then
                k = size(objects(i)%upper)
            else
                call finish_result(res, -1, .false., 'Constraint block requires lower or upper array')
                return
            end if
            if (k < 1) then
                call finish_result(res, -1, .false., 'Empty constraint block')
                return
            end if
            if (allocated(objects(i)%matrix)) then
                if (size(objects(i)%matrix, 1) /= k) then
                    call finish_result(res, -1, .false., 'Invalid linear constraint row count')
                    return
                end if
            end if
            if (allocated(objects(i)%lower)) then
                if (size(objects(i)%lower) /= k .or. any(ieee_is_nan(objects(i)%lower))) then
                    call finish_result(res, -1, .false., 'Invalid constraint block lower bounds')
                    return
                end if
            end if
            if (allocated(objects(i)%upper)) then
                if (size(objects(i)%upper) /= k .or. any(ieee_is_nan(objects(i)%upper))) then
                    call finish_result(res, -1, .false., 'Invalid constraint block upper bounds')
                    return
                end if
            end if
            if (allocated(objects(i)%keep_feasible)) then
                if (size(objects(i)%keep_feasible) /= k) then
                    call finish_result(res, -1, .false., 'Invalid constraint block keep_feasible size')
                    return
                end if
            end if
            if (objects(i)%fd_scheme /= '2-point' .and. objects(i)%fd_scheme /= '3-point') then
                call finish_result(res, -1, .false., 'Unknown constraint block finite-difference scheme')
                return
            end if
            if (objects(i)%hess_fd_scheme /= '' .and. &
                    objects(i)%hess_fd_scheme /= '2-point' .and. &
                    objects(i)%hess_fd_scheme /= '3-point') then
                call finish_result(res, -1, .false., 'Unknown constraint block Hessian finite-difference scheme')
                return
            end if
            if (objects(i)%hess_update /= '') then
                select case (lowercase(trim(objects(i)%hess_update)))
                case ('bfgs', 'sr1')
                    if (associated(objects(i)%hess) .or. associated(objects(i)%hessp) .or. &
                            objects(i)%hess_fd_scheme /= '') then
                        call finish_result(res, -1, .false., &
                            'Constraint block Hessian update conflicts with another Hessian')
                        return
                    end if
                case default
                    call finish_result(res, -1, .false., 'Unknown constraint block Hessian update')
                    return
                end select
            end if
            if (objects(i)%hess_fd_scheme /= '') then
                if (associated(objects(i)%hess) .or. associated(objects(i)%hessp)) then
                    call finish_result(res, -1, .false., &
                        'Constraint block cannot combine exact and finite-difference Hessians')
                    return
                end if
                if (.not. associated(objects(i)%jac)) then
                    call finish_result(res, -1, .false., &
                        'Finite-difference constraint Hessian requires analytic Jacobian')
                    return
                end if
            end if
            if (allocated(objects(i)%finite_diff_rel_step)) then
                if (size(objects(i)%finite_diff_rel_step) /= n) then
                    call finish_result(res, -1, .false., 'Invalid constraint block finite_diff_rel_step size')
                    return
                end if
                if (.not. all(ieee_is_finite(objects(i)%finite_diff_rel_step))) then
                    call finish_result(res, -1, .false., 'Invalid constraint block finite_diff_rel_step')
                    return
                end if
            end if
            if (allocated(objects(i)%finite_diff_jac_sparsity)) then
                if (size(objects(i)%finite_diff_jac_sparsity, 1) /= k .or. &
                        size(objects(i)%finite_diff_jac_sparsity, 2) /= n) then
                    call finish_result(res, -1, .false., &
                        'Invalid constraint block finite_diff_jac_sparsity shape')
                    return
                end if
            end if
            first(i) = m + 1
            m = m + k
            last(i) = m
            any_hess = any_hess .or. associated(objects(i)%hess) .or. &
                associated(objects(i)%hessp) .or. objects(i)%hess_fd_scheme /= '' .or. &
                objects(i)%hess_update /= ''
            all_hess = all_hess .and. &
                (allocated(objects(i)%matrix) .or. associated(objects(i)%hess) .or. &
                associated(objects(i)%hessp) .or. objects(i)%hess_fd_scheme /= '' .or. &
                objects(i)%hess_update /= '')
            any_update = any_update .or. objects(i)%hess_update /= ''
        end do
        ! An all-linear collection has a known zero constraint Hessian.
        ! Pass it to the exact-Hessian SQP path when the objective Hessian is
        ! available; do not start a spurious constraint BFGS approximation.
        ! When only some nonlinear blocks have analytic Hessians, missing
        ! Hessians are computed by differencing that block's Jacobian.
        if (all_hess .and. present(hess)) any_hess = .true.
        allocate(lo(m), hi(m), keep(m), cache_x(n, nb), cache_c(m), cache_valid(nb))
        lo = -infinity
        hi = infinity
        keep = .false.
        cache_valid = .false.
        if (any_update) then
            ! Separate updates for each scalar output allow multipliers to change
            ! without contaminating a secant pair with a multiplier change.
            allocate(update_x(n, nb), update_jac(n, m), update_hess(n, n, m), &
                update_jac_x(n, nb), update_jac_work(n, m), &
                update_valid(nb), update_jac_ready(nb), update_first(m))
            update_x = 0.0_dp
            update_jac_x = 0.0_dp
            update_jac_work = 0.0_dp
            update_jac_ready = .false.
            update_jac = 0.0_dp
            update_hess = 0.0_dp
            update_valid = .false.
            update_first = .true.
            do i = 1, m
                do k = 1, n
                    update_hess(k, k, i) = 1.0_dp
                end do
            end do
        end if
        do i = 1, nb
            if (allocated(objects(i)%lower)) lo(first(i):last(i)) = objects(i)%lower
            if (allocated(objects(i)%upper)) hi(first(i):last(i)) = objects(i)%upper
            if (allocated(objects(i)%keep_feasible)) then
                keep(first(i):last(i)) = objects(i)%keep_feasible
            end if
        end do
        if (any(lo > hi)) then
            call finish_result(res, -1, .false., 'Invalid constraint block interval')
            return
        end if

        ! Only bounds marked keep_feasible restrict finite-difference probes.
        allocate(bound_lo(n), bound_hi(n), fd_nbd(n))
        call setup_box_bounds(n, lower, upper, bound_lo, bound_hi, fd_nbd, invalid)
        if (invalid) then
            call finish_result(res, -1, .false., 'Invalid trust-constr variable bounds')
            return
        end if
        if (present(bounds_keep_feasible)) then
            if (size(bounds_keep_feasible) /= n) then
                call finish_result(res, -1, .false., 'Invalid bounds_keep_feasible size')
                return
            end if
            where (.not. bounds_keep_feasible) fd_nbd = 0
        else
            fd_nbd = 0
        end if

        if (any_hess) then
            res = trust_constr(fun, x0, jac=jac, hess=hess, lower=lower, upper=upper, &
                constraints=joined_constraints, constraint_jac=joined_jacobian, &
                constraint_hess=joined_hessian, constraint_lower=lo, constraint_upper=hi, &
                constraint_keep_feasible=keep, bounds_keep_feasible=bounds_keep_feasible, &
                gtol=gtol, xtol=xtol, eps=eps, maxiter=maxiter, &
                initial_trust_radius=initial_trust_radius, initial_tr_radius=initial_tr_radius, &
                initial_constr_penalty=initial_constr_penalty, barrier_tol=barrier_tol, &
                initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, &
                objective_hess_update=objective_hess_update, &
                callback=callback)
        else
            res = trust_constr(fun, x0, jac=jac, hess=hess, lower=lower, upper=upper, &
                constraints=joined_constraints, constraint_jac=joined_jacobian, &
                constraint_lower=lo, constraint_upper=hi, &
                constraint_keep_feasible=keep, bounds_keep_feasible=bounds_keep_feasible, &
                gtol=gtol, xtol=xtol, eps=eps, maxiter=maxiter, &
                initial_trust_radius=initial_trust_radius, initial_tr_radius=initial_tr_radius, &
                initial_constr_penalty=initial_constr_penalty, barrier_tol=barrier_tol, &
                initial_barrier_parameter=initial_barrier_parameter, &
                initial_barrier_tolerance=initial_barrier_tolerance, &
                objective_fd_scheme=objective_fd_scheme, &
                objective_finite_diff_rel_step=objective_finite_diff_rel_step, &
                objective_hess_update=objective_hess_update, &
                callback=callback)
        end if

        ! SciPy reports function/derivative counts independently for each object,
        ! with an additional zero-cost block for variable bounds, if present.
        if (allocated(res%constr_nfev)) then
            if (present(lower) .or. present(upper)) then
                res%constr_nfev = [block_nfev, 0]
                res%constr_njev = [block_njev, 0]
                res%constr_nhev = [block_nhev, 0]
            else
                res%constr_nfev = block_nfev
                res%constr_njev = block_njev
                res%constr_nhev = block_nhev
            end if
        end if

    contains

        subroutine joined_constraints(xv, cv)
            real(dp), intent(in) :: xv(:) !! point to evaluate
            real(dp), intent(out) :: cv(:) !! original values in user block order

            integer :: block

            do block = 1, nb
                if (allocated(objects(block)%matrix)) then
                    cv(first(block):last(block)) = matmul(objects(block)%matrix, xv)
                else
                    call objects(block)%fun(xv, cv(first(block):last(block)))
                    block_nfev(block) = block_nfev(block) + 1
                end if
                cache_x(:, block) = xv
                cache_c(first(block):last(block)) = cv(first(block):last(block))
                cache_valid(block) = .true.
            end do
        end subroutine joined_constraints

        subroutine joined_jacobian(xv, amat)
            real(dp), intent(in) :: xv(:) !! point to differentiate
            real(dp), intent(out) :: amat(:, :) !! block-row-stacked Jacobian
            call joined_jacobian_block(xv, amat, 0)
        end subroutine joined_jacobian

        subroutine joined_jacobian_block(xv, amat, only_block)
            real(dp), intent(in) :: xv(:) !! point to differentiate
            real(dp), intent(out) :: amat(:, :) !! block-row-stacked Jacobian
            integer, intent(in) :: only_block !! zero for all blocks; otherwise selected block

            integer :: block
            integer :: dim
            integer :: col
            real(dp) :: hstep
            real(dp) :: h1
            real(dp) :: h2
            real(dp) :: xprobe(n)
            real(dp) :: base(m)
            real(dp) :: probe(m)
            real(dp) :: probe2(m)

            amat = 0.0_dp
            do block = 1, nb
                if (only_block > 0 .and. block /= only_block) cycle
                if (allocated(objects(block)%matrix)) then
                    amat(first(block):last(block), :) = objects(block)%matrix
                    cycle
                end if
                if (associated(objects(block)%jac)) then
                    block_njev(block) = block_njev(block) + 1
                    call objects(block)%jac(xv, amat(first(block):last(block), :))
                    call cache_updated_jacobian(xv, amat, block)
                    cycle
                end if
                dim = last(block) - first(block) + 1
                if (cache_valid(block)) then
                    if (same_vector(xv, cache_x(:, block))) then
                        base(1:dim) = cache_c(first(block):last(block))
                    else
                        call objects(block)%fun(xv, base(1:dim))
                        block_nfev(block) = block_nfev(block) + 1
                        cache_x(:, block) = xv
                        cache_c(first(block):last(block)) = base(1:dim)
                        cache_valid(block) = .true.
                    end if
                else
                    call objects(block)%fun(xv, base(1:dim))
                    block_nfev(block) = block_nfev(block) + 1
                    cache_x(:, block) = xv
                    cache_c(first(block):last(block)) = base(1:dim)
                    cache_valid(block) = .true.
                end if
                if (allocated(objects(block)%finite_diff_jac_sparsity)) then
                    call colored_jacobian_block(xv, base(1:dim), &
                        amat(first(block):last(block), :), block)
                    call cache_updated_jacobian(xv, amat, block)
                    cycle
                end if
                do col = 1, n
                    hstep = eps_use * max(1.0_dp, abs(xv(col)))
                    if (objects(block)%fd_scheme == '3-point' .and. .not. present(eps)) then
                        hstep = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(xv(col)))
                    end if
                    if (allocated(objects(block)%finite_diff_rel_step)) then
                        hstep = objective_relative_step(xv(col), &
                            hstep / max(1.0_dp, abs(xv(col))), &
                            objects(block)%finite_diff_rel_step(col))
                    end if
                    if (objects(block)%fd_scheme == '3-point') then
                        call bounded_three_point_fd_steps(xv(col), bound_lo(col), bound_hi(col), &
                            fd_nbd(col), abs(hstep), h1, h2)
                        if (h1 == 0.0_dp .or. h2 == 0.0_dp .or. h1 == h2) cycle
                        xprobe = xv
                        xprobe(col) = xv(col) + h1
                        if (fd_nbd(col) == 1 .or. fd_nbd(col) == 2) then
                            xprobe(col) = max(xprobe(col), bound_lo(col))
                        end if
                        if (fd_nbd(col) == 2 .or. fd_nbd(col) == 3) then
                            xprobe(col) = min(xprobe(col), bound_hi(col))
                        end if
                        h1 = xprobe(col) - xv(col)
                        if (h1 == 0.0_dp) cycle
                        call objects(block)%fun(xprobe, probe(1:dim))
                        block_nfev(block) = block_nfev(block) + 1
                        xprobe(col) = xv(col) + h2
                        if (fd_nbd(col) == 1 .or. fd_nbd(col) == 2) then
                            xprobe(col) = max(xprobe(col), bound_lo(col))
                        end if
                        if (fd_nbd(col) == 2 .or. fd_nbd(col) == 3) then
                            xprobe(col) = min(xprobe(col), bound_hi(col))
                        end if
                        h2 = xprobe(col) - xv(col)
                        if (h2 == 0.0_dp .or. h2 == h1) then
                            ! A very narrow retained interval may admit only one
                            ! representable probe. Fall back to a first-order slope.
                            amat(first(block):last(block), col) = &
                                (probe(1:dim) - base(1:dim)) / h1
                            cycle
                        end if
                        call objects(block)%fun(xprobe, probe2(1:dim))
                        block_nfev(block) = block_nfev(block) + 1
                        ! Unequal actual floating-point steps also work for
                        ! shortened/one-sided probes near strict bounds.
                        amat(first(block):last(block), col) = &
                            (h2 * ((probe(1:dim) - base(1:dim)) / h1) - &
                            h1 * ((probe2(1:dim) - base(1:dim)) / h2)) / (h2 - h1)
                    else
                        if (fd_nbd(col) /= 0) then
                            h1 = hstep
                            call bounded_signed_fd_step(xv(col), bound_lo(col), bound_hi(col), &
                                fd_nbd(col), h1, hstep)
                        end if
                        if (hstep == 0.0_dp) cycle
                        xprobe = xv
                        xprobe(col) = xv(col) + hstep
                        hstep = xprobe(col) - xv(col)
                        if (hstep == 0.0_dp) cycle
                        call objects(block)%fun(xprobe, probe(1:dim))
                        block_nfev(block) = block_nfev(block) + 1
                        amat(first(block):last(block), col) = &
                            (probe(1:dim) - base(1:dim)) / hstep
                    end if
                end do
                call cache_updated_jacobian(xv, amat, block)
            end do

        end subroutine joined_jacobian_block

        subroutine cache_updated_jacobian(xv, amat, block)
            ! Reuse the Jacobian already evaluated at this iterate. As in
            ! SciPy, a quasi-Newton constraint Hessian does not require a
            ! second user-Jacobian evaluation at the same point.
            real(dp), intent(in) :: xv(:) !! evaluation point
            real(dp), intent(in) :: amat(:, :) !! row-stacked full Jacobian
            integer, intent(in) :: block !! constraint block to cache

            if (objects(block)%hess_update == '') return
            update_jac_x(:, block) = xv
            update_jac_work(:, first(block):last(block)) = &
                transpose(amat(first(block):last(block), :))
            update_jac_ready(block) = .true.
        end subroutine cache_updated_jacobian

        subroutine colored_jacobian_block(xv, base, deriv, block)
            real(dp), intent(in) :: xv(:) !! current point
            real(dp), intent(in) :: base(:) !! cached vector constraint value at current point
            real(dp), intent(out) :: deriv(:, :) !! sparse finite-difference Jacobian
            integer, intent(in) :: block !! constraint block to differentiate

            integer :: col
            integer :: group
            integer :: ngroup
            integer :: previous
            integer :: row
            integer :: group_id(n)
            real(dp) :: h1(n)
            real(dp) :: h2(n)
            real(dp) :: hstep
            real(dp) :: probe1(size(base))
            real(dp) :: probe2(size(base))
            real(dp) :: xprobe1(n)
            real(dp) :: xprobe2(n)
            logical :: first_needed
            logical :: second_needed
            logical :: conflicts

            ! Greedy coloring: columns in the same group have no common
            ! structurally nonzero rows, so their probes can be combined.
            group_id = 0
            ngroup = 0
            do col = 1, n
                if (.not. any(objects(block)%finite_diff_jac_sparsity(:, col))) cycle
                do group = 1, ngroup
                    conflicts = .false.
                    do previous = 1, col - 1
                        if (group_id(previous) /= group) cycle
                        if (any(objects(block)%finite_diff_jac_sparsity(:, previous) .and. &
                                objects(block)%finite_diff_jac_sparsity(:, col))) then
                            conflicts = .true.
                            exit
                        end if
                    end do
                    if (.not. conflicts) exit
                end do
                if (group > ngroup) ngroup = group
                group_id(col) = group
            end do

            deriv = 0.0_dp
            do group = 1, ngroup
                h1 = 0.0_dp
                h2 = 0.0_dp
                xprobe1 = xv
                xprobe2 = xv
                first_needed = .false.
                second_needed = .false.
                do col = 1, n
                    if (group_id(col) /= group) cycle
                    hstep = eps_use * max(1.0_dp, abs(xv(col)))
                    if (objects(block)%fd_scheme == '3-point' .and. .not. present(eps)) then
                        hstep = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * &
                            max(1.0_dp, abs(xv(col)))
                    end if
                    if (allocated(objects(block)%finite_diff_rel_step)) then
                        hstep = objective_relative_step(xv(col), &
                            hstep / max(1.0_dp, abs(xv(col))), &
                            objects(block)%finite_diff_rel_step(col))
                    end if

                    if (objects(block)%fd_scheme == '3-point') then
                        call bounded_three_point_fd_steps(xv(col), bound_lo(col), &
                            bound_hi(col), fd_nbd(col), abs(hstep), h1(col), h2(col))
                    else
                        if (fd_nbd(col) /= 0) then
                            h1(col) = hstep
                            call bounded_signed_fd_step(xv(col), bound_lo(col), bound_hi(col), &
                                fd_nbd(col), h1(col), hstep)
                        end if
                        h1(col) = hstep
                    end if
                    if (h1(col) == 0.0_dp) cycle
                    xprobe1(col) = xv(col) + h1(col)
                    if (fd_nbd(col) == 1 .or. fd_nbd(col) == 2) then
                        xprobe1(col) = max(xprobe1(col), bound_lo(col))
                    end if
                    if (fd_nbd(col) == 2 .or. fd_nbd(col) == 3) then
                        xprobe1(col) = min(xprobe1(col), bound_hi(col))
                    end if
                    h1(col) = xprobe1(col) - xv(col)
                    if (h1(col) == 0.0_dp) cycle
                    first_needed = .true.
                    if (objects(block)%fd_scheme /= '3-point') cycle
                    if (h2(col) == 0.0_dp) cycle
                    xprobe2(col) = xv(col) + h2(col)
                    if (fd_nbd(col) == 1 .or. fd_nbd(col) == 2) then
                        xprobe2(col) = max(xprobe2(col), bound_lo(col))
                    end if
                    if (fd_nbd(col) == 2 .or. fd_nbd(col) == 3) then
                        xprobe2(col) = min(xprobe2(col), bound_hi(col))
                    end if
                    h2(col) = xprobe2(col) - xv(col)
                    if (h2(col) /= 0.0_dp .and. h2(col) /= h1(col)) second_needed = .true.
                end do
                if (.not. first_needed) cycle
                call objects(block)%fun(xprobe1, probe1)
                block_nfev(block) = block_nfev(block) + 1
                if (second_needed) then
                    call objects(block)%fun(xprobe2, probe2)
                    block_nfev(block) = block_nfev(block) + 1
                end if
                do col = 1, n
                    if (group_id(col) /= group .or. h1(col) == 0.0_dp) cycle
                    do row = 1, size(base)
                        if (.not. objects(block)%finite_diff_jac_sparsity(row, col)) cycle
                        if (second_needed .and. h2(col) /= 0.0_dp .and. &
                                h2(col) /= h1(col)) then
                            deriv(row, col) = (h2(col) * ((probe1(row) - base(row)) / h1(col)) - &
                                h1(col) * ((probe2(row) - base(row)) / h2(col))) / &
                                (h2(col) - h1(col))
                        else
                            deriv(row, col) = (probe1(row) - base(row)) / h1(col)
                        end if
                    end do
                end do
            end do
        end subroutine colored_jacobian_block

        subroutine joined_hessian(xv, vv, hmat)
            real(dp), intent(in) :: xv(:) !! point to evaluate
            real(dp), intent(in) :: vv(:) !! multipliers in raw user order
            real(dp), intent(out) :: hmat(:, :) !! weighted sum of block Hessians

            real(dp) :: part(n, n)
            real(dp) :: direction(n)
            real(dp) :: jplus(m, n)
            real(dp) :: jminus(m, n)
            real(dp) :: jbase(m, n)
            real(dp) :: xprobe(n)
            real(dp) :: hstep
            real(dp) :: h1
            real(dp) :: h2
            logical :: base_valid
            logical :: two_point_hessian
            integer :: block
            integer :: col

            hmat = 0.0_dp
            do block = 1, nb
                if (allocated(objects(block)%matrix)) cycle
                if (associated(objects(block)%hess)) then
                    ! An explicitly supplied full Hessian takes precedence.
                    call objects(block)%hess(xv, vv(first(block):last(block)), part)
                else if (associated(objects(block)%hessp)) then
                    ! A Fortran counterpart of a SciPy LinearOperator Hessian:
                    ! assemble its columns by applying the product to e_j.
                    ! The dense SQP kernel needs a symmetric Hessian matrix.
                    do col = 1, n
                        direction = 0.0_dp
                        direction(col) = 1.0_dp
                        call objects(block)%hessp(xv, vv(first(block):last(block)), &
                            direction, part(:, col))
                    end do
                    part = 0.5_dp * (part + transpose(part))
                else if (objects(block)%hess_update /= '') then
                    call updated_constraint_hessian(xv, vv, block, part)
                else
                    ! Differentiate only this block's Jacobian, not the full
                    ! stacked Jacobian. This leaves the counts and caches for
                    ! unrelated blocks untouched. The larger outer step when
                    ! the Jacobian is numerical controls roundoff amplification.
                    part = 0.0_dp
                    two_point_hessian = objects(block)%hess_fd_scheme == '2-point'
                    base_valid = .false.
                    hstep = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
                    if (two_point_hessian) hstep = sqrt(epsilon(1.0_dp))
                    if (.not. associated(objects(block)%jac)) then
                        hstep = epsilon(1.0_dp)**0.25_dp
                    end if
                    do col = 1, n
                        h1 = hstep * max(1.0_dp, abs(xv(col)))
                        if (allocated(objects(block)%finite_diff_rel_step)) then
                            h1 = objective_relative_step(xv(col), hstep, &
                                objects(block)%finite_diff_rel_step(col))
                        end if
                        if (two_point_hessian) then
                            if (fd_nbd(col) /= 0) then
                                call bounded_signed_fd_step(xv(col), bound_lo(col), &
                                    bound_hi(col), fd_nbd(col), h1, h2)
                                h1 = h2
                            end if
                            xprobe = xv
                            xprobe(col) = xv(col) + h1
                            h1 = xprobe(col) - xv(col)
                            if (h1 == 0.0_dp) cycle
                            if (.not. base_valid) then
                                call joined_jacobian_block(xv, jbase, block)
                                base_valid = .true.
                            end if
                            call joined_jacobian_block(xprobe, jplus, block)
                            part(:, col) = matmul(transpose( &
                                jplus(first(block):last(block), :) - &
                                jbase(first(block):last(block), :)), &
                                vv(first(block):last(block))) / h1
                            cycle
                        end if
                        call bounded_three_point_fd_steps(xv(col), bound_lo(col), &
                            bound_hi(col), fd_nbd(col), &
                            abs(h1), h1, h2)
                        if (h1 == 0.0_dp .or. h2 == 0.0_dp) cycle
                        xprobe = xv
                        xprobe(col) = xv(col) + h1
                        h1 = xprobe(col) - xv(col)
                        if (h1 == 0.0_dp) cycle
                        call joined_jacobian_block(xprobe, jminus, block)
                        xprobe(col) = xv(col) + h2
                        h2 = xprobe(col) - xv(col)
                        if (h2 == h1 .or. h2 == 0.0_dp) cycle
                        call joined_jacobian_block(xprobe, jplus, block)
                        if (sign(1.0_dp, h1) /= sign(1.0_dp, h2)) then
                            ! Centered difference, allowing rounded steps.
                            part(:, col) = matmul(transpose( &
                                jplus(first(block):last(block), :) - &
                                jminus(first(block):last(block), :)), &
                                vv(first(block):last(block))) / (h2 - h1)
                        else
                            ! Unequal one-sided steps: differentiate the
                            ! quadratic interpolant through x, x+h1, x+h2.
                            if (.not. base_valid) then
                                call joined_jacobian_block(xv, jbase, block)
                                base_valid = .true.
                            end if
                            part(:, col) = matmul(transpose( &
                                (h2 / h1) * (jminus(first(block):last(block), :) - &
                                jbase(first(block):last(block), :)) - &
                                (h1 / h2) * (jplus(first(block):last(block), :) - &
                                jbase(first(block):last(block), :))), &
                                vv(first(block):last(block))) / (h2 - h1)
                        end if
                    end do
                    part = 0.5_dp * (part + transpose(part))
                end if
                hmat = hmat + part
                if (associated(objects(block)%hess) .or. &
                        associated(objects(block)%hessp)) then
                    block_nhev(block) = block_nhev(block) + 1
                end if
            end do
        end subroutine joined_hessian

        subroutine updated_constraint_hessian(xv, vv, block, part)
            ! Per-output quasi-Newton Hessian updates.  Each Hessian models one
            ! scalar constraint; current multipliers are applied after updating.
            ! This avoids a spurious secant term when multipliers vary between
            ! SQP iterations, and leaves unrelated constraint blocks untouched.
            real(dp), intent(in) :: xv(:) !! current physical optimization variables
            real(dp), intent(in) :: vv(:) !! multipliers in block order
            integer, intent(in) :: block !! nonlinear constraint block number
            real(dp), intent(out) :: part(:, :) !! weighted block Hessian

            real(dp), parameter :: min_denominator = 1.0e-8_dp
            real(dp) :: current_jac(last(block) - first(block) + 1, n)
            real(dp) :: stacked_jac(m, n)
            real(dp) :: step(n)
            real(dp) :: y(n)
            real(dp) :: bs(n)
            real(dp) :: correction(n)
            real(dp) :: denominator
            real(dp) :: sy
            real(dp) :: y2
            real(dp) :: scale
            real(dp) :: theta
            integer :: component
            integer :: j
            integer :: row

            if (update_jac_ready(block)) then
                if (all(xv == update_jac_x(:, block))) then
                    current_jac = transpose(update_jac_work(:, first(block):last(block)))
                else
                    call joined_jacobian_block(xv, stacked_jac, block)
                    current_jac = stacked_jac(first(block):last(block), :)
                end if
            else
                call joined_jacobian_block(xv, stacked_jac, block)
                current_jac = stacked_jac(first(block):last(block), :)
            end if
            part = 0.0_dp
            step = 0.0_dp
            if (update_valid(block)) step = xv - update_x(:, block)
            if (update_valid(block) .and. any(step /= 0.0_dp)) then
                do component = 1, size(current_jac, 1)
                    row = first(block) + component - 1
                    y = current_jac(component, :) - update_jac(:, row)
                    sy = dot_product(step, y)
                    y2 = dot_product(y, y)
                    if (update_first(row)) then
                        if (sy /= 0.0_dp .and. y2 > 0.0_dp) then
                            scale = y2 / abs(sy)
                            if (ieee_is_finite(scale) .and. scale > 0.0_dp) then
                                update_hess(:, :, row) = 0.0_dp
                                do j = 1, n
                                    update_hess(j, j, row) = scale
                                end do
                            end if
                        end if
                        update_first(row) = .false.
                    end if
                    bs = matmul(update_hess(:, :, row), step)
                    if (lowercase(trim(objects(block)%hess_update)) == 'sr1') then
                        correction = y - bs
                        denominator = dot_product(step, correction)
                        if (.not. ieee_is_finite(denominator)) cycle
                        if (abs(denominator) <= min_denominator * norm2(step) * &
                                norm2(correction)) cycle
                        do j = 1, n
                            update_hess(:, j, row) = update_hess(:, j, row) + &
                                correction * (correction(j) / denominator)
                        end do
                    else
                        ! Powell damping retains positive-definite BFGS updates
                        ! even when a nonlinear constraint has negative curvature.
                        denominator = dot_product(step, bs)
                        if (.not. ieee_is_finite(denominator) .or. &
                                denominator <= epsilon(1.0_dp) * norm2(step) * norm2(bs)) cycle
                        correction = y
                        if (sy < 0.2_dp * denominator) then
                            theta = 0.8_dp * denominator / (denominator - sy)
                            correction = theta * y + (1.0_dp - theta) * bs
                        end if
                        sy = dot_product(step, correction)
                        if (.not. ieee_is_finite(sy) .or. &
                                sy <= epsilon(1.0_dp) * norm2(step) * norm2(correction)) cycle
                        do j = 1, n
                            update_hess(:, j, row) = update_hess(:, j, row) - &
                                bs * (bs(j) / denominator) + correction * (correction(j) / sy)
                        end do
                    end if
                end do
            end if
            update_x(:, block) = xv
            do component = 1, size(current_jac, 1)
                row = first(block) + component - 1
                update_jac(:, row) = current_jac(component, :)
                part = part + vv(row) * update_hess(:, :, row)
            end do
            update_valid(block) = .true.
        end subroutine updated_constraint_hessian

    end function trust_constr_objects

    function trust_constr_equality_sqp(fun, x0, jac, hess, lower, upper, constraints, &
            constraint_jac, meq, mineq, gtol, xtol, eps, maxiter, &
            initial_trust_radius, initial_constr_penalty, constraint_hess, use_soc, &
            positive_tail, approx_constraint_hess, initial_fun, initial_grad, &
            initial_objective_hess, final_objective_hess, initial_constraints, &
            initial_constraint_jac, enforce_tail, fd_lower, fd_upper, fd_nbd, &
            objective_three_point, objective_fd_step, objective_fd_rel_step, &
            objective_hess_update, callback) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic objective gradient
        procedure(multivariate_hessian_function), optional :: hess !! analytic objective Hessian
        real(dp), intent(in), optional :: lower(:) !! reserved for the inequality/interior-point path
        real(dp), intent(in), optional :: upper(:) !! reserved for the inequality/interior-point path
        procedure(multivariate_constraint_function), optional :: constraints !! equality constraints
        procedure(multivariate_constraint_jacobian), optional :: constraint_jac !! equality Jacobian
        integer, intent(in), optional :: meq !! number of equality constraints
        integer, intent(in), optional :: mineq !! must be zero in the current equality-only path
        real(dp), intent(in), optional :: gtol !! Lagrangian-gradient/feasibility tolerance
        real(dp), intent(in), optional :: xtol !! trust-region-radius tolerance
        real(dp), intent(in), optional :: eps !! finite-difference step
        integer, intent(in), optional :: maxiter !! maximum SQP iterations
        real(dp), intent(in), optional :: initial_trust_radius !! initial trust-region radius
        real(dp), intent(in), optional :: initial_constr_penalty !! initial merit-function constraint penalty
        procedure(multivariate_constraint_hessian), optional :: constraint_hess !! Hessian of weighted constraints
        logical, intent(in), optional :: use_soc !! enable second-order correction for rejected SQP steps
        integer, intent(in), optional :: positive_tail !! trailing positive variables protected during SOC
        logical, intent(in), optional :: approx_constraint_hess !! use SciPy-style BFGS for missing constraint Hessian
        real(dp), intent(in), optional :: initial_fun !! cached objective value at the starting point
        real(dp), intent(in), optional :: initial_grad(:) !! cached objective gradient at the starting point
        real(dp), intent(in), optional :: initial_objective_hess(:, :) !! cached objective Hessian at the start
        real(dp), intent(out), optional :: final_objective_hess(:, :) !! final objective Hessian cache
        real(dp), intent(in), optional :: initial_constraints(:) !! cached equality values at the starting point
        real(dp), intent(in), optional :: initial_constraint_jac(:, :) !! cached equality Jacobian at the start
        logical, intent(in), optional :: enforce_tail(:) !! trailing barrier slacks tied to feasible inequalities
        real(dp), intent(in), optional :: fd_lower(:) !! strict lower bounds for numerical derivative probes
        real(dp), intent(in), optional :: fd_upper(:) !! strict upper bounds for numerical derivative probes
        integer, intent(in), optional :: fd_nbd(:) !! bound-type codes for strict numerical derivative probes
        logical, intent(in), optional :: objective_three_point !! central/one-sided second-order objective derivative
        real(dp), intent(in), optional :: objective_fd_step !! objective derivative absolute step scale
        real(dp), intent(in), optional :: objective_fd_rel_step(:) !! relative steps for physical variables
        character(len=*), intent(in), optional :: objective_hess_update !! BFGS or SR1 objective Hessian update
        procedure(trust_constr_iteration_callback), optional :: callback !! called at each SQP state
        type(multivariate_minimize_result) :: res

        real(dp), parameter :: penalty_factor = 0.3_dp
        real(dp), parameter :: large_ratio = 0.9_dp
        real(dp), parameter :: intermediary_ratio = 0.3_dp
        real(dp), parameter :: sufficient_ratio = 1.0e-8_dp
        real(dp), parameter :: soc_threshold = 0.1_dp
        real(dp), parameter :: trust_enlarge_large = 7.0_dp
        real(dp), parameter :: trust_enlarge_small = 2.0_dp
        real(dp), parameter :: max_trust_reduction = 0.5_dp
        real(dp), parameter :: min_trust_reduction = 0.1_dp
        real(dp), parameter :: normal_trust_factor = 0.8_dp
        real(dp), parameter :: boundary_parameter = 0.995_dp
        integer :: cg_step_niter
        integer :: cg_step_stop_cond
        integer :: i
        integer :: m
        integer :: maxiter_use
        integer :: n
        integer :: constraint_nfev
        integer :: constraint_njev
        integer :: constraint_nhev
        integer :: positive_tail_use
        integer :: linalg_status
        logical :: accepted
        logical :: approx_constraint_hess_use
        logical, allocatable :: enforce_tail_use(:)
        logical :: first_constraint_bfgs
        logical :: first_objective_sr1
        logical :: use_sr1
        logical :: projection_ok
        logical :: terminated
        logical :: callback_stop
        logical :: three_point_use
        logical :: use_cached_objective_hess
        logical :: use_soc_now
        logical :: soc_intersect
        real(dp) :: actual_reduction
        real(dp) :: constr_violation
        real(dp) :: delta
        real(dp) :: dnorm
        real(dp) :: eps_use
        real(dp) :: objective_fd_step_use
        real(dp) :: f
        real(dp) :: f_trial
        real(dp) :: gtol_use
        real(dp) :: initial_penalty_use
        real(dp) :: merit
        real(dp) :: merit_trial
        real(dp) :: optimality
        real(dp) :: penalty
        real(dp) :: predicted_reduction
        real(dp) :: previous_penalty
        real(dp) :: quadratic_model
        real(dp) :: reduction_ratio
        real(dp) :: reduction_ratio_soc
        real(dp) :: soc_scale
        real(dp) :: trust_reduction
        real(dp) :: actual_reduction_soc
        real(dp) :: f_soc
        real(dp) :: merit_soc
        real(dp) :: vpred
        real(dp) :: xtol_use
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: a_phys(:, :)
        real(dp), allocatable :: a_prev(:, :)
        real(dp), allocatable :: bmat(:, :)
        real(dp), allocatable :: bcon(:, :)
        real(dp), allocatable :: bs(:)
        real(dp), allocatable :: bscon(:)
        real(dp), allocatable :: cvec(:)
        real(dp), allocatable :: c_trial(:)
        real(dp), allocatable :: c_soc(:)
        real(dp), allocatable :: constr_hmat(:, :)
        real(dp), allocatable :: ct(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: dn(:)
        real(dp), allocatable :: dt(:)
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: evec(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: g_scaled(:)
        real(dp), allocatable :: g_trial(:)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: hmat(:, :)
        real(dp), allocatable :: h_scaled(:, :)
        real(dp), allocatable :: lag_grad(:)
        real(dp), allocatable :: lag_grad_trial(:)
        real(dp), allocatable :: last_objective_hess(:, :)
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: lambda_trial(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: scale_vec(:)
        real(dp), allocatable :: soc_step(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_trial(:)
        real(dp), allocatable :: x_soc(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: ycon(:)
        real(dp), allocatable :: ymod(:)
        real(dp) :: sbs
        real(dp) :: sy
        real(dp) :: theta

        n = size(x0)
        call init_result(res, 'equality_constrained_sqp', n)
        constraint_nfev = 0
        constraint_njev = 0
        constraint_nhev = 0
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        m = 0
        if (present(meq)) m = meq
        if (m < 0 .or. (m > 0 .and. .not. present(constraints))) then
            call finish_result(res, -1, .false., &
                'Equality-constrained trust-constr requires constraints when meq > 0')
            return
        end if
        ! An overdetermined equality system is accepted when its rows are
        ! compatible.  The Gram pseudoinverse handles redundant rows.
        if (present(mineq)) then
            if (mineq /= 0) then
                call finish_result(res, -1, .false., &
                    'Inequality trust-constr is not implemented yet')
                return
            end if
        end if
        if (present(lower) .or. present(upper)) then
            call finish_result(res, -1, .false., &
                'Bounded trust-constr is not implemented yet')
            return
        end if
        if ((present(fd_lower) .neqv. present(fd_upper)) .or. &
                (present(fd_lower) .neqv. present(fd_nbd))) then
            call finish_result(res, -1, .false., &
                'Incomplete finite-difference bound specification')
            return
        end if
        if (present(fd_nbd)) then
            if (size(fd_lower) /= n .or. size(fd_upper) /= n .or. size(fd_nbd) /= n) then
                call finish_result(res, -1, .false., &
                    'Invalid finite-difference bound size')
                return
            end if
        end if
        gtol_use = 1.0e-8_dp
        if (present(gtol)) gtol_use = gtol
        xtol_use = 1.0e-8_dp
        if (present(xtol)) xtol_use = xtol
        eps_use = default_fd_eps
        if (present(eps)) eps_use = eps
        three_point_use = .false.
        if (present(objective_three_point)) three_point_use = objective_three_point
        objective_fd_step_use = eps_use
        if (present(objective_fd_step)) objective_fd_step_use = objective_fd_step
        maxiter_use = 1000
        if (present(maxiter)) maxiter_use = maxiter
        delta = 1.0_dp
        if (present(initial_trust_radius)) delta = initial_trust_radius
        initial_penalty_use = 1.0_dp
        if (present(initial_constr_penalty)) initial_penalty_use = initial_constr_penalty
        use_soc_now = .true.
        if (present(use_soc)) use_soc_now = use_soc
        if (m == 0) use_soc_now = .false.
        positive_tail_use = 0
        if (present(positive_tail)) positive_tail_use = positive_tail
        allocate(enforce_tail_use(positive_tail_use))
        enforce_tail_use = .false.
        if (present(enforce_tail)) then
            if (size(enforce_tail) /= positive_tail_use) return
            enforce_tail_use = enforce_tail
        end if
        approx_constraint_hess_use = .false.
        if (present(approx_constraint_hess)) approx_constraint_hess_use = approx_constraint_hess
        first_constraint_bfgs = .true.
        first_objective_sr1 = .true.
        use_sr1 = .false.
        if (present(objective_hess_update)) &
            use_sr1 = lowercase(trim(objective_hess_update)) == 'sr1'
        if (use_sr1 .and. m > 0 .and. positive_tail_use == 0 .and. &
                .not. present(constraint_hess) .and. .not. present(hess)) then
            ! For SR1, the objective Hessian update must use objective
            ! gradient differences alone. Approximate missing constraint
            ! curvature separately instead of folding it into SR1.
            approx_constraint_hess_use = .true.
        end if
        if (present(objective_fd_rel_step)) then
            if (size(objective_fd_rel_step) > n) return
        end if
        if (present(initial_grad)) then
            if (size(initial_grad) /= n) return
        end if
        if (present(initial_objective_hess)) then
            if (size(initial_objective_hess, 1) /= n .or. &
                    size(initial_objective_hess, 2) /= n) return
        end if
        if (present(final_objective_hess)) then
            if (size(final_objective_hess, 1) /= n .or. &
                    size(final_objective_hess, 2) /= n) return
        end if
        if (present(initial_constraints)) then
            if (size(initial_constraints) /= m) return
        end if
        if (present(initial_constraint_jac)) then
            if (size(initial_constraint_jac, 1) /= m .or. &
                    size(initial_constraint_jac, 2) /= n) return
        end if
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp .or. &
                .not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp .or. &
                .not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp .or. &
                .not. ieee_is_finite(objective_fd_step_use) .or. objective_fd_step_use <= 0.0_dp .or. &
                .not. ieee_is_finite(delta) .or. delta <= 0.0_dp .or. &
                .not. ieee_is_finite(initial_penalty_use) .or. &
                maxiter_use < 1 .or. positive_tail_use < 0 .or. positive_tail_use > n .or. &
                (approx_constraint_hess_use .and. present(constraint_hess))) return

        allocate(x(n), x_trial(n), g(n), g_scaled(n), g_trial(n), hmat(n, n), &
            h_scaled(n, n), bmat(n, n), bcon(n, n), constr_hmat(n, n), cvec(m), c_trial(m), &
            c_soc(m), a(m, n), a_phys(m, n), a_prev(m, n), gram(m, m), eig(m), evec(m, m), &
            lambda(m), lambda_trial(m), lag_grad(n), lag_grad_trial(n), &
            last_objective_hess(n, n), &
            dn(n), dt(n), d(n), ct(n), s(n), scale_vec(n), soc_step(n), &
            x_soc(n), y(n), ycon(n), ymod(n), bs(n), bscon(n))

        x = x0
        if (positive_tail_use > 0) then
            if (any(x(n - positive_tail_use + 1:n) <= 0.0_dp)) then
                call finish_result(res, -1, .false., 'Positive SOC variables must start above zero')
                return
            end if
        end if
        if (present(initial_fun)) then
            f = initial_fun
            res%nfev = 0
        else
            f = fun(x)
            res%nfev = 1
        end if
        if (present(initial_grad)) then
            g = initial_grad
        else
            call evaluate_gradient(x, f, g)
        end if
        if (present(initial_constraints)) then
            cvec = initial_constraints
        else
            call evaluate_constraints(x, cvec)
        end if
        if (present(initial_constraint_jac)) then
            a_phys = initial_constraint_jac
        else
            call evaluate_constraint_jacobian(x, cvec, a_phys)
        end if
        call update_step_scaling(x, scale_vec)
        call scale_gradient(g, scale_vec, g_scaled)
        call scale_jacobian(a_phys, scale_vec, a)
        if (.not. ieee_is_finite(f) .or. .not. all(ieee_is_finite(g)) .or. &
                .not. all(ieee_is_finite(cvec)) .or. .not. all(ieee_is_finite(a))) then
            res%x = x
            res%fun = f
            res%jac = g
            call finish_result(res, -1, .false., 'Nonfinite initial objective or constraints')
            return
        end if

        use_cached_objective_hess = present(initial_objective_hess)
        last_objective_hess = 0.0_dp

        if (.not. present(hess)) then
            bmat = 0.0_dp
            do i = 1, n
                bmat(i, i) = 1.0_dp
            end do
            hmat = bmat
        end if

        bcon = 0.0_dp
        if (approx_constraint_hess_use) then
            do i = 1, n - positive_tail_use
                bcon(i, i) = 1.0_dp
            end do
            if (.not. present(hess)) hmat = hmat + bcon
        end if

        penalty = initial_penalty_use
        terminated = .false.
        call update_projection(a, projection_ok)
        if (.not. projection_ok) then
            res%x = x
            res%fun = f
            res%jac = g
            call finish_result(res, -1, .false., 'Constraint Jacobian projection failed')
            return
        end if
        call compute_multiplier(a, g_scaled, lambda)
        lag_grad = g + matmul(transpose(a_phys), lambda)
        if (present(constraint_hess) .and. .not. present(initial_constraints)) then
            lambda_trial = 0.0_dp
            call constraint_hess(x, lambda_trial, constr_hmat)
            constraint_nhev = constraint_nhev + 1
            if (.not. all(ieee_is_finite(constr_hmat))) then
                res%x = x
                res%fun = f
                res%jac = g
                call finish_result(res, -1, .false., 'Nonfinite initial constraint Hessian')
                return
            end if
        end if
        if (present(hess)) then
            call assemble_lagrangian_hessian(x, lambda, hmat, projection_ok)
            if (.not. projection_ok) then
                res%x = x
                res%fun = f
                res%jac = g
                call finish_result(res, -1, .false., 'Nonfinite initial Hessian')
                return
            end if
        else if (present(constraint_hess)) then
            ! With a BFGS objective Hessian and an exact constraint Hessian,
            ! add the latter at the current multipliers instead of updating
            ! BFGS with the whole (already curvature-corrected) Lagrangian.
            call constraint_hess(x, lambda, constr_hmat)
            constraint_nhev = constraint_nhev + 1
            hmat = bmat + constr_hmat
            if (.not. all(ieee_is_finite(hmat))) then
                call finish_result(res, -1, .false., 'Nonfinite initial constraint Hessian')
                return
            end if
        end if
        call scale_hessian(hmat, scale_vec, lambda, h_scaled)
        optimality = maxval(abs(scale_vec * lag_grad))
        constr_violation = 0.0_dp
        if (m > 0) constr_violation = maxval(abs(cvec))

        ! SciPy counts the initial stop-condition state as iteration one.  This
        ! matters for barrier restarts because each stage enters the equality
        ! SQP solver once before taking another trial step.
        res%nit = 1
        do
            if (present(callback)) then
                callback_stop = .false.
                call callback(x, res%nit, callback_stop)
                if (callback_stop) then
                    call finish_result(res, 3, .false., '`callback` function requested termination')
                    terminated = .true.
                    exit
                end if
            end if
            if (optimality <= gtol_use .and. constr_violation <= gtol_use) then
                call finish_result(res, 1, .true., '`gtol` termination condition is satisfied')
                terminated = .true.
                exit
            end if
            if (delta < xtol_use) then
                call finish_result(res, 2, .true., '`xtol` termination condition is satisfied')
                terminated = .true.
                exit
            end if
            if (res%nit >= maxiter_use) then
                call finish_result(res, 0, .false., 'Maximum number of iterations reached')
                terminated = .true.
                exit
            end if

            call normal_step(cvec, delta, dn)
            ct = matmul(h_scaled, dn) + g_scaled
            call tangential_step(ct, dn, delta, dt, cg_step_niter, cg_step_stop_cond)
            res%cg_niter = res%cg_niter + cg_step_niter
            res%cg_stop_cond = cg_step_stop_cond
            d = dn + dt
            dnorm = norm2(d)
            if (.not. all(ieee_is_finite(d)) .or. dnorm <= epsilon(1.0_dp)) then
                delta = min_trust_reduction * delta
                res%nit = res%nit + 1
                cycle
            end if

            quadratic_model = 0.5_dp * dot_product(d, matmul(h_scaled, d)) + &
                dot_product(g_scaled, d)
            vpred = norm2(cvec) - norm2(matmul(a, d) + cvec)
            vpred = max(1.0e-16_dp, vpred)
            previous_penalty = penalty
            if (quadratic_model > 0.0_dp) then
                penalty = max(penalty, quadratic_model / ((1.0_dp - penalty_factor) * vpred))
            end if
            predicted_reduction = -quadratic_model + penalty * vpred
            if (.not. ieee_is_finite(predicted_reduction) .or. predicted_reduction <= 0.0_dp) then
                predicted_reduction = max(1.0e-16_dp, abs(predicted_reduction))
            end if

            merit = f + penalty * norm2(cvec)
            x_trial = x + scale_vec * d
            if (any(enforce_tail_use)) then
                call evaluate_constraints(x_trial, c_trial)
                call enforce_feasible_tail(x_trial, c_trial)
                f_trial = fun(x_trial)
                res%nfev = res%nfev + 1
            else
                f_trial = fun(x_trial)
                res%nfev = res%nfev + 1
                call evaluate_constraints(x_trial, c_trial)
            end if
            if (.not. ieee_is_finite(f_trial) .or. .not. all(ieee_is_finite(c_trial))) then
                reduction_ratio = -huge(1.0_dp)
            else
                merit_trial = f_trial + penalty * norm2(c_trial)
                actual_reduction = merit - merit_trial
                reduction_ratio = actual_reduction / predicted_reduction
            end if

            if (use_soc_now .and. reduction_ratio < sufficient_ratio .and. &
                    norm2(dn) <= soc_threshold * norm2(dt) .and. &
                    all(ieee_is_finite(c_trial))) then
                call second_order_correction(c_trial, soc_step)
                call soc_slack_intersection(d, soc_step, soc_scale, soc_intersect)
                if (soc_intersect) then
                    x_soc = x + scale_vec * (d + soc_scale * soc_step)
                    if (any(enforce_tail_use)) then
                        call evaluate_constraints(x_soc, c_soc)
                        call enforce_feasible_tail(x_soc, c_soc)
                        f_soc = fun(x_soc)
                        res%nfev = res%nfev + 1
                    else
                        f_soc = fun(x_soc)
                        res%nfev = res%nfev + 1
                        call evaluate_constraints(x_soc, c_soc)
                    end if
                    if (ieee_is_finite(f_soc) .and. all(ieee_is_finite(c_soc))) then
                        merit_soc = f_soc + penalty * norm2(c_soc)
                        actual_reduction_soc = merit - merit_soc
                        reduction_ratio_soc = actual_reduction_soc / predicted_reduction
                        if (reduction_ratio_soc >= sufficient_ratio) then
                            x_trial = x_soc
                            f_trial = f_soc
                            c_trial = c_soc
                            reduction_ratio = reduction_ratio_soc
                        end if
                    end if
                end if
            end if

            if (reduction_ratio >= large_ratio) then
                delta = max(trust_enlarge_large * dnorm, delta)
            else if (reduction_ratio >= intermediary_ratio) then
                delta = max(trust_enlarge_small * dnorm, delta)
            else if (reduction_ratio < sufficient_ratio) then
                trust_reduction = (1.0_dp - sufficient_ratio) / &
                    max(1.0_dp - reduction_ratio, 1.0_dp)
                if (trust_reduction * dnorm >= max_trust_reduction * delta) then
                    delta = max_trust_reduction * delta
                else if (trust_reduction * dnorm >= min_trust_reduction * delta) then
                    delta = trust_reduction * dnorm
                else
                    delta = min_trust_reduction * delta
                end if
            end if

            accepted = reduction_ratio >= sufficient_ratio
            if (accepted) then
                if (approx_constraint_hess_use) a_prev = a_phys
                s = x_trial - x
                x = x_trial
                f = f_trial
                cvec = c_trial
                call evaluate_gradient(x, f, g_trial)
                call evaluate_constraint_jacobian(x, cvec, a_phys)
                call update_step_scaling(x, scale_vec)
                call scale_gradient(g_trial, scale_vec, g_scaled)
                call scale_jacobian(a_phys, scale_vec, a)
                if (.not. all(ieee_is_finite(g_trial)) .or. .not. all(ieee_is_finite(a))) then
                    call finish_result(res, -1, .false., 'Nonfinite derivative encountered')
                    terminated = .true.
                    exit
                end if
                call update_projection(a, projection_ok)
                if (.not. projection_ok) then
                    call finish_result(res, -1, .false., 'Constraint Jacobian projection failed')
                    terminated = .true.
                    exit
                end if
                call compute_multiplier(a, g_scaled, lambda_trial)
                lag_grad_trial = g_trial + matmul(transpose(a_phys), lambda_trial)

                if (approx_constraint_hess_use) then
                    ycon = matmul(transpose(a_phys - a_prev), lambda_trial)
                    call update_constraint_bfgs(s, ycon)
                end if

                if (present(hess)) then
                    call assemble_lagrangian_hessian(x, lambda_trial, hmat, projection_ok)
                    if (.not. projection_ok) then
                        call finish_result(res, -1, .false., 'Nonfinite Hessian encountered')
                        terminated = .true.
                        exit
                    end if
                else
                    if (use_sr1 .or. approx_constraint_hess_use .or. &
                            present(constraint_hess)) then
                        y = g_trial - g
                    else
                        y = lag_grad_trial - lag_grad
                    end if
                    if (use_sr1) then
                        call update_objective_sr1(s, y)
                    else
                        bs = matmul(bmat, s)
                        sbs = dot_product(s, bs)
                        sy = dot_product(s, y)
                        ymod = y
                        if (sbs > epsilon(1.0_dp) * max(1.0_dp, norm2(s) * norm2(bs))) then
                            if (sy < 0.2_dp * sbs) then
                                theta = 0.8_dp * sbs / max(sbs - sy, epsilon(1.0_dp))
                                ymod = theta * y + (1.0_dp - theta) * bs
                            end if
                            sy = dot_product(s, ymod)
                            if (sy > epsilon(1.0_dp) * max(1.0_dp, norm2(s) * norm2(ymod))) then
                                do i = 1, n
                                    bmat(:, i) = bmat(:, i) - bs * bs(i) / sbs + &
                                        ymod * ymod(i) / sy
                                end do
                            end if
                        end if
                    end if
                    hmat = bmat
                    if (approx_constraint_hess_use) then
                        hmat = hmat + bcon
                    else if (present(constraint_hess)) then
                        call constraint_hess(x, lambda_trial, constr_hmat)
                        constraint_nhev = constraint_nhev + 1
                        hmat = hmat + constr_hmat
                        if (.not. all(ieee_is_finite(hmat))) then
                            call finish_result(res, -1, .false., 'Nonfinite constraint Hessian encountered')
                            terminated = .true.
                            exit
                        end if
                    end if
                end if

                g = g_trial
                lambda = lambda_trial
                lag_grad = lag_grad_trial
                call scale_hessian(hmat, scale_vec, lambda, h_scaled)
                optimality = maxval(abs(scale_vec * lag_grad))
                constr_violation = 0.0_dp
                if (m > 0) constr_violation = maxval(abs(cvec))
            else
                penalty = previous_penalty
            end if
            res%nit = res%nit + 1
        end do

        if (.not. terminated) then
            call finish_result(res, 0, .false., 'Maximum number of iterations reached')
        end if
        res%x = x
        res%fun = f
        res%jac = g
        res%maxcv = constr_violation
        res%constr_violation = constr_violation
        res%optimality = optimality
        res%lagrangian_grad = lag_grad
        res%constr_penalty = penalty
        res%tr_radius = delta
        res%v = lambda
        res%constr = cvec
        res%constr_jac = a_phys
        if (m == 0) then
            allocate(res%constr_nfev(0), res%constr_njev(0), res%constr_nhev(0))
        else
            allocate(res%constr_nfev(1), res%constr_njev(1), res%constr_nhev(1))
            res%constr_nfev(1) = constraint_nfev
            res%constr_njev(1) = constraint_njev
            res%constr_nhev(1) = constraint_nhev
        end if
        if (present(final_objective_hess) .and. present(hess)) then
            final_objective_hess = last_objective_hess
        end if

    contains

        subroutine assemble_lagrangian_hessian(xv, vv, hv, ok)
            real(dp), intent(in) :: xv(:) !! evaluation point
            real(dp), intent(in) :: vv(:) !! current equality multipliers
            real(dp), intent(out) :: hv(:, :) !! Lagrangian Hessian
            logical, intent(out) :: ok !! true when all Hessian entries are finite

            if (use_cached_objective_hess) then
                hv = initial_objective_hess
                use_cached_objective_hess = .false.
            else
                call invoke_hessian(hess, xv, hv)
                res%nhev = res%nhev + 1
            end if
            last_objective_hess = hv
            if (present(constraint_hess)) then
                call invoke_constraint_hessian(constraint_hess, xv, vv, constr_hmat)
                constraint_nhev = constraint_nhev + 1
                hv = hv + constr_hmat
            else if (approx_constraint_hess_use) then
                hv = hv + bcon
            end if
            ok = all(ieee_is_finite(hv))
        end subroutine assemble_lagrangian_hessian

        subroutine update_objective_sr1(step_phys, delta_grad)
            real(dp), intent(in) :: step_phys(:) !! accepted objective step
            real(dp), intent(in) :: delta_grad(:) !! objective gradient difference

            real(dp), parameter :: min_denominator = 1.0e-8_dp
            real(dp) :: denominator
            real(dp) :: ys
            real(dp) :: y2
            real(dp) :: scale
            real(dp) :: correction(n)
            integer :: j

            if (first_objective_sr1) then
                ! SciPy's auto-scaled initial SR1 matrix uses ||y||^2/|s.y|.
                ys = abs(dot_product(step_phys, delta_grad))
                y2 = dot_product(delta_grad, delta_grad)
                if (ys > 0.0_dp .and. y2 > 0.0_dp) then
                    scale = y2 / ys
                    if (ieee_is_finite(scale)) then
                        bmat = 0.0_dp
                        do j = 1, n
                            bmat(j, j) = scale
                        end do
                    end if
                end if
                first_objective_sr1 = .false.
            end if

            correction = delta_grad - matmul(bmat, step_phys)
            denominator = dot_product(step_phys, correction)
            if (abs(denominator) <= min_denominator * norm2(step_phys) * norm2(correction)) return
            if (.not. ieee_is_finite(denominator)) return
            do j = 1, n
                bmat(:, j) = bmat(:, j) + correction * (correction(j) / denominator)
            end do
        end subroutine update_objective_sr1

        subroutine update_constraint_bfgs(step_phys, delta_grad)
            real(dp), intent(in) :: step_phys(:) !! accepted physical-coordinate displacement
            real(dp), intent(in) :: delta_grad(:) !! change in weighted constraint gradient

            integer :: j
            integer :: nbase
            real(dp) :: scale
            real(dp) :: sbs_local
            real(dp) :: sy_local
            real(dp) :: y2

            nbase = n - positive_tail_use
            if (nbase <= 0) return
            sy_local = dot_product(step_phys(1:nbase), delta_grad(1:nbase))
            y2 = dot_product(delta_grad(1:nbase), delta_grad(1:nbase))
            if (dot_product(step_phys(1:nbase), step_phys(1:nbase)) == 0.0_dp .or. &
                    y2 == 0.0_dp .or. sy_local == 0.0_dp) return

            if (first_constraint_bfgs) then
                scale = y2 / abs(sy_local)
                bcon = 0.0_dp
                do j = 1, nbase
                    bcon(j, j) = scale
                end do
                first_constraint_bfgs = .false.
            end if

            bscon = matmul(bcon, step_phys)
            sbs_local = dot_product(step_phys, bscon)
            if (sbs_local <= 0.0_dp) then
                scale = y2 / abs(sy_local)
                bcon = 0.0_dp
                do j = 1, nbase
                    bcon(j, j) = scale
                end do
                bscon = matmul(bcon, step_phys)
                sbs_local = dot_product(step_phys, bscon)
            end if
            if (sy_local <= 1.0e-8_dp * sbs_local) return

            do j = 1, nbase
                bcon(1:nbase, j) = bcon(1:nbase, j) - &
                    bscon(1:nbase) * bscon(j) / sbs_local + &
                    delta_grad(1:nbase) * delta_grad(j) / sy_local
            end do
        end subroutine update_constraint_bfgs

        subroutine enforce_feasible_tail(xv, cv)
            real(dp), intent(inout) :: xv(:) !! barrier variables, including trailing slacks
            real(dp), intent(inout) :: cv(:) !! equality residuals at xv

            integer :: k
            integer :: row
            integer :: tail_index

            if (positive_tail_use <= 0) return
            do k = 1, positive_tail_use
                if (.not. enforce_tail_use(k)) cycle
                row = m - positive_tail_use + k
                tail_index = n - positive_tail_use + k
                xv(tail_index) = xv(tail_index) + cv(row)
                cv(row) = 0.0_dp
            end do
        end subroutine enforce_feasible_tail

        subroutine evaluate_constraints(xv, cv)
            real(dp), intent(in) :: xv(:) !! point at which to evaluate the user constraints
            real(dp), intent(out) :: cv(:) !! user constraint values at xv

            if (m == 0) return
            call invoke_constraints(constraints, xv, cv)
            constraint_nfev = constraint_nfev + 1
        end subroutine evaluate_constraints

        subroutine evaluate_gradient(xv, fv, gv)
            real(dp), intent(in) :: xv(:) !! evaluation point
            real(dp), intent(in) :: fv !! already-known objective value at xv
            real(dp), intent(out) :: gv(:) !! objective gradient

            integer :: k
            real(dp) :: fp
            real(dp) :: h
            real(dp) :: requested_h
            real(dp) :: xp(n)

            if (present(jac)) then
                call invoke_gradient(jac, xv, gv)
            else if (three_point_use) then
                if (present(fd_nbd)) then
                    call three_point_gradient(fun, xv, fv, objective_fd_step_use, gv, k, &
                        fd_lower, fd_upper, fd_nbd, objective_fd_rel_step)
                else
                    call three_point_gradient(fun, xv, fv, objective_fd_step_use, gv, k, &
                        rel_step=objective_fd_rel_step)
                end if
                res%nfev = res%nfev + k
            else
                do k = 1, n
                    h = eps_use * max(1.0_dp, abs(xv(k)))
                    if (present(objective_fd_rel_step)) then
                        if (k <= size(objective_fd_rel_step)) then
                            h = objective_relative_step(xv(k), eps_use, objective_fd_rel_step(k))
                        end if
                    end if
                    if (present(fd_nbd)) then
                        if (fd_nbd(k) /= 0) then
                            requested_h = h
                            call bounded_signed_fd_step(xv(k), fd_lower(k), fd_upper(k), &
                                fd_nbd(k), requested_h, h)
                        end if
                    end if
                    xp = xv
                    xp(k) = xp(k) + h
                    h = xp(k) - xv(k)
                    if (h == 0.0_dp) then
                        gv(k) = 0.0_dp
                        cycle
                    end if
                    fp = fun(xp)
                    res%nfev = res%nfev + 1
                    gv(k) = (fp - fv) / h
                end do
            end if
            res%njev = res%njev + 1
        end subroutine evaluate_gradient

        subroutine evaluate_constraint_jacobian(xv, cv, av)
            real(dp), intent(in) :: xv(:) !! evaluation point
            real(dp), intent(in) :: cv(:) !! already-known constraints at xv
            real(dp), intent(out) :: av(:, :) !! equality constraint Jacobian

            integer :: k
            real(dp) :: h
            real(dp) :: cp(m)
            real(dp) :: xp(n)

            if (m == 0) return
            if (present(constraint_jac)) then
                call invoke_constraint_jacobian(constraint_jac, xv, av)
                constraint_njev = constraint_njev + 1
                return
            end if
            do k = 1, n
                h = eps_use * max(1.0_dp, abs(xv(k)))
                if (present(fd_nbd)) then
                    if (fd_nbd(k) /= 0) then
                        call bounded_fd_step(xv(k), fd_lower(k), fd_upper(k), &
                            fd_nbd(k), abs(h), h)
                    end if
                end if
                xp = xv
                xp(k) = xp(k) + h
                call evaluate_constraints(xp, cp)
                av(:, k) = (cp - cv) / h
            end do
        end subroutine evaluate_constraint_jacobian

        subroutine update_projection(av, ok)
            real(dp), intent(in) :: av(:, :) !! current equality Jacobian
            logical, intent(out) :: ok !! true when eigendecomposition succeeds

            if (m == 0) then
                ok = .true.
                return
            end if
            gram = matmul(av, transpose(av))
            call symmetric_eigen_jacobi(gram, eig, evec, linalg_status)
            ok = linalg_status == linalg_status_success .and. &
                all(ieee_is_finite(eig)) .and. all(ieee_is_finite(evec))
        end subroutine update_projection

        subroutine solve_gram(rhs, sol)
            real(dp), intent(in) :: rhs(:) !! right-hand side in constraint space
            real(dp), intent(out) :: sol(:) !! pseudoinverse solution

            integer :: k
            real(dp) :: cutoff
            real(dp) :: coeff(m)

            if (m == 0) then
                sol = 0.0_dp
                return
            end if
            coeff = matmul(transpose(evec), rhs)
            cutoff = epsilon(1.0_dp) * real(max(n, m), dp) * &
                max(1.0_dp, maxval(abs(eig)))
            do k = 1, m
                if (eig(k) > cutoff) then
                    coeff(k) = coeff(k) / eig(k)
                else
                    coeff(k) = 0.0_dp
                end if
            end do
            sol = matmul(evec, coeff)
        end subroutine solve_gram

        subroutine compute_multiplier(av, gv, vv)
            real(dp), intent(in) :: av(:, :) !! equality Jacobian
            real(dp), intent(in) :: gv(:) !! objective gradient
            real(dp), intent(out) :: vv(:) !! least-squares Lagrange multipliers

            real(dp) :: rhs(m)

            rhs = matmul(av, gv)
            call solve_gram(rhs, vv)
            vv = -vv
        end subroutine compute_multiplier

        subroutine project_null(v, pv)
            real(dp), intent(in) :: v(:) !! vector to project
            real(dp), intent(out) :: pv(:) !! null-space projection

            real(dp) :: rhs(m)
            real(dp) :: z(m)

            rhs = matmul(a, v)
            call solve_gram(rhs, z)
            pv = v - matmul(transpose(a), z)
        end subroutine project_null

        subroutine second_order_correction(cv, step)
            real(dp), intent(in) :: cv(:) !! nonlinear residual at the rejected trial point
            real(dp), intent(out) :: step(:) !! minimum-norm correction from the current Jacobian

            real(dp) :: z(m)

            call solve_gram(cv, z)
            step = -matmul(transpose(a), z)
        end subroutine second_order_correction

        subroutine soc_slack_intersection(base, direction, scale, intersect)
            real(dp), intent(in) :: base(:) !! rejected SQP step before SOC, in scaled coordinates
            real(dp), intent(in) :: direction(:) !! second-order correction direction
            real(dp), intent(out) :: scale !! largest feasible correction fraction in [0, 1]
            logical, intent(out) :: intersect !! true when the SOC segment remains slack-feasible

            integer :: k
            integer :: first_positive
            real(dp) :: lower_step
            real(dp) :: t
            real(dp) :: t_lower
            real(dp) :: t_upper

            scale = 0.0_dp
            intersect = .true.
            if (positive_tail_use <= 0) then
                scale = 1.0_dp
                return
            end if
            t_lower = 0.0_dp
            t_upper = 1.0_dp
            lower_step = -boundary_parameter
            first_positive = n - positive_tail_use + 1
            do k = first_positive, n
                if (direction(k) > 0.0_dp) then
                    t = (lower_step - base(k)) / direction(k)
                    t_lower = max(t_lower, t)
                else if (direction(k) < 0.0_dp) then
                    t = (lower_step - base(k)) / direction(k)
                    t_upper = min(t_upper, t)
                else if (base(k) < lower_step) then
                    intersect = .false.
                    return
                end if
            end do
            t_lower = max(0.0_dp, t_lower)
            t_upper = min(1.0_dp, t_upper)
            intersect = t_lower <= t_upper
            if (intersect) scale = max(0.0_dp, t_upper)
        end subroutine soc_slack_intersection

        subroutine update_step_scaling(xv, sv)
            real(dp), intent(in) :: xv(:) !! current physical variables
            real(dp), intent(out) :: sv(:) !! diagonal step scaling

            integer :: first_positive

            sv = 1.0_dp
            if (positive_tail_use <= 0) return
            first_positive = n - positive_tail_use + 1
            sv(first_positive:n) = xv(first_positive:n)
        end subroutine update_step_scaling

        subroutine scale_gradient(gv, sv, gwork)
            real(dp), intent(in) :: gv(:) !! physical-coordinate gradient
            real(dp), intent(in) :: sv(:) !! diagonal step scaling
            real(dp), intent(out) :: gwork(:) !! gradient in scaled step coordinates

            gwork = sv * gv
        end subroutine scale_gradient

        subroutine scale_jacobian(av, sv, awork)
            real(dp), intent(in) :: av(:, :) !! physical-coordinate Jacobian
            real(dp), intent(in) :: sv(:) !! diagonal step scaling
            real(dp), intent(out) :: awork(:, :) !! Jacobian in scaled step coordinates

            integer :: k

            do k = 1, n
                awork(:, k) = av(:, k) * sv(k)
            end do
        end subroutine scale_jacobian

        subroutine scale_hessian(hv, sv, vv, hwork)
            real(dp), intent(in) :: hv(:, :) !! physical-coordinate Lagrangian Hessian
            real(dp), intent(in) :: sv(:) !! diagonal step scaling
            real(dp), intent(in) :: vv(:) !! equality multipliers
            real(dp), intent(out) :: hwork(:, :) !! Hessian in scaled step coordinates

            integer :: first_multiplier
            integer :: first_positive
            integer :: j
            integer :: k
            real(dp) :: primal_dual

            do j = 1, n
                do k = 1, n
                    hwork(j, k) = sv(j) * hv(j, k) * sv(k)
                end do
            end do

            ! SciPy's barrier Hessian uses the primal-dual diagonal when the
            ! slack multiplier has the correct sign.  Our barrier equality is
            ! c(x) - s = 0, so its multiplier is negative at a KKT point.
            if (positive_tail_use > 0 .and. present(hess)) then
                first_positive = n - positive_tail_use + 1
                first_multiplier = m - positive_tail_use + 1
                do k = 0, positive_tail_use - 1
                    if (vv(first_multiplier + k) < 0.0_dp) then
                        primal_dual = -vv(first_multiplier + k) * sv(first_positive + k)
                        if (primal_dual > 0.0_dp .and. ieee_is_finite(primal_dual)) then
                            hwork(first_positive + k, first_positive + k) = primal_dual
                        end if
                    end if
                end do
            end if
        end subroutine scale_hessian

        logical function inside_step_box(step, lb_step, ub_step) result(ok)
            real(dp), intent(in) :: step(:) !! trial scaled step
            real(dp), intent(in) :: lb_step(:) !! lower step bounds
            real(dp), intent(in) :: ub_step(:) !! upper step bounds

            ok = all(step >= lb_step .and. step <= ub_step)
        end function inside_step_box

        subroutine reinforce_step_box(step, lb_step, ub_step)
            real(dp), intent(inout) :: step(:) !! scaled step clipped in place
            real(dp), intent(in) :: lb_step(:) !! lower step bounds
            real(dp), intent(in) :: ub_step(:) !! upper step bounds

            step = min(max(step, lb_step), ub_step)
        end subroutine reinforce_step_box

        subroutine box_sphere_upper(zv, direction, lb_step, ub_step, radius, &
                entire_line, t_upper, intersect)
            real(dp), intent(in) :: zv(:) !! line origin
            real(dp), intent(in) :: direction(:) !! line direction
            real(dp), intent(in) :: lb_step(:) !! lower box bounds
            real(dp), intent(in) :: ub_step(:) !! upper box bounds
            real(dp), intent(in) :: radius !! spherical trust radius
            logical, intent(in) :: entire_line !! allow all real line parameters
            real(dp), intent(out) :: t_upper !! largest feasible line parameter
            logical, intent(out) :: intersect !! true when box and sphere intersect

            integer :: k
            real(dp) :: aa
            real(dp) :: bb
            real(dp) :: cc
            real(dp) :: disc
            real(dp) :: root1
            real(dp) :: root2
            real(dp) :: t1
            real(dp) :: t2
            real(dp) :: t_lower
            real(dp) :: bound_guard

            bound_guard = 0.5_dp * huge(1.0_dp)
            if (entire_line) then
                t_lower = -huge(1.0_dp)
                t_upper = huge(1.0_dp)
            else
                t_lower = 0.0_dp
                t_upper = 1.0_dp
            end if

            do k = 1, size(zv)
                if (abs(direction(k)) <= tiny(1.0_dp)) then
                    if (zv(k) < lb_step(k) .or. zv(k) > ub_step(k)) then
                        intersect = .false.
                        t_upper = 0.0_dp
                        return
                    end if
                    cycle
                end if
                if (lb_step(k) > -bound_guard) then
                    t1 = (lb_step(k) - zv(k)) / direction(k)
                    if (direction(k) > 0.0_dp) then
                        t_lower = max(t_lower, t1)
                    else
                        t_upper = min(t_upper, t1)
                    end if
                end if
                if (ub_step(k) < bound_guard) then
                    t2 = (ub_step(k) - zv(k)) / direction(k)
                    if (direction(k) > 0.0_dp) then
                        t_upper = min(t_upper, t2)
                    else
                        t_lower = max(t_lower, t2)
                    end if
                end if
            end do
            if (t_lower > t_upper) then
                intersect = .false.
                t_upper = 0.0_dp
                return
            end if

            aa = dot_product(direction, direction)
            if (aa <= tiny(1.0_dp)) then
                intersect = norm2(zv) <= radius .and. t_lower <= t_upper
                if (.not. intersect) t_upper = 0.0_dp
                return
            end if
            bb = 2.0_dp * dot_product(zv, direction)
            cc = dot_product(zv, zv) - radius**2
            disc = bb**2 - 4.0_dp * aa * cc
            if (disc < 0.0_dp) then
                intersect = .false.
                t_upper = 0.0_dp
                return
            end if
            disc = sqrt(max(0.0_dp, disc))
            root1 = (-bb - disc) / (2.0_dp * aa)
            root2 = (-bb + disc) / (2.0_dp * aa)
            if (root1 > root2) then
                t1 = root1
                root1 = root2
                root2 = t1
            end if
            t_lower = max(t_lower, root1)
            t_upper = min(t_upper, root2)
            intersect = t_lower <= t_upper
            if (.not. intersect) t_upper = 0.0_dp
        end subroutine box_sphere_upper

        subroutine normal_step(cv, trust_radius, step)
            real(dp), intent(in) :: cv(:) !! equality residual
            real(dp), intent(in) :: trust_radius !! current trust radius
            real(dp), intent(out) :: step(:) !! normal feasibility step

            integer :: first_positive
            real(dp) :: alpha
            real(dp) :: denom
            real(dp) :: limit
            real(dp) :: resid1
            real(dp) :: resid2
            real(dp) :: zsolve(m)
            real(dp) :: ag(m)
            real(dp) :: cauchy(n)
            real(dp) :: gradient(n)
            real(dp) :: lb_step(n)
            real(dp) :: newton(n)
            real(dp) :: origin(n)
            real(dp) :: pdir(n)
            real(dp) :: x1(n)
            real(dp) :: x2(n)
            real(dp) :: ub_step(n)
            logical :: intersect

            limit = normal_trust_factor * trust_radius
            call solve_gram(cv, zsolve)
            newton = -matmul(transpose(a), zsolve)

            if (positive_tail_use <= 0) then
                step = newton
                if (norm2(step) > limit .and. norm2(step) > 0.0_dp) then
                    step = step * limit / norm2(step)
                end if
                return
            end if

            lb_step = -huge(1.0_dp)
            ub_step = huge(1.0_dp)
            first_positive = n - positive_tail_use + 1
            lb_step(first_positive:n) = -0.5_dp * boundary_parameter
            if (inside_step_box(newton, lb_step, ub_step) .and. norm2(newton) <= limit) then
                step = newton
                return
            end if

            gradient = matmul(transpose(a), cv)
            ag = matmul(a, gradient)
            denom = dot_product(ag, ag)
            if (denom <= tiny(1.0_dp)) then
                origin = 0.0_dp
                pdir = newton
                call box_sphere_upper(origin, pdir, lb_step, ub_step, limit, &
                    .false., alpha, intersect)
                if (intersect) then
                    step = origin + alpha * pdir
                else
                    step = 0.0_dp
                end if
                call reinforce_step_box(step, lb_step, ub_step)
                return
            end if
            cauchy = -dot_product(gradient, gradient) / denom * gradient
            origin = 0.0_dp

            pdir = newton - cauchy
            call box_sphere_upper(cauchy, pdir, lb_step, ub_step, limit, &
                .false., alpha, intersect)
            if (intersect) then
                x1 = cauchy + alpha * pdir
            else
                pdir = cauchy
                call box_sphere_upper(origin, pdir, lb_step, ub_step, limit, &
                    .false., alpha, intersect)
                if (intersect) then
                    x1 = alpha * pdir
                else
                    x1 = 0.0_dp
                end if
            end if

            pdir = newton
            call box_sphere_upper(origin, pdir, lb_step, ub_step, limit, &
                .false., alpha, intersect)
            if (intersect) then
                x2 = alpha * pdir
            else
                x2 = 0.0_dp
            end if
            call reinforce_step_box(x1, lb_step, ub_step)
            call reinforce_step_box(x2, lb_step, ub_step)
            resid1 = norm2(matmul(a, x1) + cv)
            resid2 = norm2(matmul(a, x2) + cv)
            if (resid1 < resid2) then
                step = x1
            else
                step = x2
            end if
        end subroutine normal_step

        subroutine tangential_step(c_t, normal, trust_radius, step, cg_niter, cg_stop_cond)
            real(dp), intent(in) :: c_t(:) !! tangential QP linear term
            real(dp), intent(in) :: normal(:) !! normal step
            real(dp), intent(in) :: trust_radius !! full trust radius
            real(dp), intent(out) :: step(:) !! tangential null-space step
            integer, intent(out) :: cg_niter !! projected-CG iterations used by this subproblem
            integer, intent(out) :: cg_stop_cond !! SciPy projected-CG termination code

            integer :: counter
            integer :: first_positive
            integer :: k
            integer :: max_cg
            real(dp) :: alpha
            real(dp) :: beta
            real(dp) :: curvature
            real(dp) :: rr
            real(dp) :: rr_new
            real(dp) :: tau
            real(dp) :: tol_cg
            real(dp) :: trust_t
            real(dp) :: candidate(n)
            real(dp) :: gproj(n)
            real(dp) :: gnew(n)
            real(dp) :: hp(n)
            real(dp) :: last_feasible(n)
            real(dp) :: lb_step(n)
            real(dp) :: p(n)
            real(dp) :: r(n)
            real(dp) :: rnew(n)
            real(dp) :: ub_step(n)
            logical :: intersect

            step = 0.0_dp
            cg_niter = 0
            cg_stop_cond = 0
            trust_t = sqrt(max(0.0_dp, trust_radius**2 - norm2(normal)**2))
            if (trust_t <= 1.0e-25_dp) then
                cg_stop_cond = 2
                return
            end if

            lb_step = -huge(1.0_dp)
            ub_step = huge(1.0_dp)
            if (positive_tail_use > 0) then
                first_positive = n - positive_tail_use + 1
                lb_step(first_positive:n) = -boundary_parameter - normal(first_positive:n)
            end if
            last_feasible = 0.0_dp
            counter = 0

            call project_null(c_t, r)
            call project_null(r, gproj)
            rr = dot_product(gproj, gproj)
            tol_cg = max(min(0.01_dp * sqrt(max(0.0_dp, rr)), 0.1_dp * rr), 1.0e-25_dp)
            p = -gproj
            max_cg = max(0, n - m)
            cg_stop_cond = 1
            do k = 1, max_cg
                if (rr < tol_cg) then
                    cg_stop_cond = 4
                    exit
                end if
                cg_niter = cg_niter + 1
                hp = matmul(h_scaled, p)
                curvature = dot_product(p, hp)
                if (curvature <= 0.0_dp) then
                    call box_sphere_upper(step, p, lb_step, ub_step, trust_t, &
                        .true., tau, intersect)
                    if (intersect) step = step + tau * p
                    call reinforce_step_box(step, lb_step, ub_step)
                    cg_stop_cond = 3
                    return
                end if

                alpha = rr / curvature
                candidate = step + alpha * p
                if (norm2(candidate) >= trust_t) then
                    call box_sphere_upper(step, alpha * p, lb_step, ub_step, trust_t, &
                        .false., tau, intersect)
                    if (intersect) step = step + tau * alpha * p
                    call reinforce_step_box(step, lb_step, ub_step)
                    cg_stop_cond = 2
                    return
                end if

                if (inside_step_box(candidate, lb_step, ub_step)) then
                    last_feasible = candidate
                    counter = 0
                else
                    counter = counter + 1
                end if
                if (counter > 0) then
                    call box_sphere_upper(step, alpha * p, lb_step, ub_step, trust_t, &
                        .false., tau, intersect)
                    if (intersect) then
                        last_feasible = step + tau * alpha * p
                        call reinforce_step_box(last_feasible, lb_step, ub_step)
                    end if
                end if
                if (counter > max_cg) exit

                rnew = r + alpha * hp
                call project_null(rnew, gnew)
                rr_new = dot_product(gnew, gnew)
                beta = rr_new / rr
                p = -gnew + beta * p
                step = candidate
                r = gnew
                gproj = gnew
                rr = rr_new
            end do
            if (.not. inside_step_box(step, lb_step, ub_step)) step = last_feasible
        end subroutine tangential_step


    end function trust_constr_equality_sqp

    function lbfgsb(fun, x0, jac, lower, upper, ftol, gtol, eps, maxiter, &
            maxfun, maxcor) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic gradient
        real(dp), intent(in), optional :: lower(:) !! lower bounds; nonfinite values mean unbounded
        real(dp), intent(in), optional :: upper(:) !! upper bounds; nonfinite values mean unbounded
        real(dp), intent(in), optional :: ftol !! relative objective-reduction tolerance
        real(dp), intent(in), optional :: gtol !! projected-gradient infinity-norm tolerance
        real(dp), intent(in), optional :: eps !! absolute forward-difference step, default 1e-8
        integer, intent(in), optional :: maxiter !! maximum completed L-BFGS-B iterations
        integer, intent(in), optional :: maxfun !! maximum objective evaluations, default 15000
        integer, intent(in), optional :: maxcor !! correction-history length, default 10
        type(multivariate_minimize_result) :: res

        character(len=60) :: csave
        character(len=60) :: task
        integer, allocatable :: iwa(:)
        integer, allocatable :: nbd(:)
        integer :: isave(44)
        integer :: maxcor_use
        integer :: maxfun_use
        integer :: maxiter_use
        integer :: n
        logical :: hit_fun_limit
        logical :: hit_iter_limit
        logical :: lsave(4)
        real(dp) :: dsave(29)
        real(dp) :: eps_use
        real(dp) :: factr
        real(dp) :: fval
        real(dp) :: ftol_use
        real(dp), allocatable :: g(:)
        real(dp) :: gtol_use
        real(dp), allocatable :: l(:)
        real(dp), allocatable :: u(:)
        real(dp), allocatable :: wa(:)
        real(dp), allocatable :: x(:)

        n = size(x0)
        call init_result(res, 'l-bfgs-b', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        maxcor_use = 10
        if (present(maxcor)) maxcor_use = maxcor
        maxiter_use = 15000
        if (present(maxiter)) maxiter_use = maxiter
        maxfun_use = 15000
        if (present(maxfun)) maxfun_use = maxfun
        ftol_use = 2.2204460492503131e-9_dp
        if (present(ftol)) ftol_use = ftol
        gtol_use = 1.0e-5_dp
        if (present(gtol)) gtol_use = gtol
        eps_use = 1.0e-8_dp
        if (present(eps)) eps_use = eps
        if (maxcor_use <= 0 .or. maxiter_use < 0 .or. maxfun_use < 0) return
        if (.not. ieee_is_finite(ftol_use) .or. ftol_use < 0.0_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp) return

        allocate(x(n), g(n), l(n), u(n), nbd(n), iwa(3 * n))
        allocate(wa(2 * maxcor_use * n + 5 * n + 11 * maxcor_use**2 + 8 * maxcor_use))
        x = x0
        call setup_lbfgsb_bounds(x, lower, upper, l, u, nbd, hit_iter_limit)
        if (hit_iter_limit) return

        fval = 0.0_dp
        g = 0.0_dp
        wa = 0.0_dp
        iwa = 0
        lsave = .false.
        isave = 0
        dsave = 0.0_dp
        csave = ''
        task = 'START'
        factr = ftol_use / epsilon(1.0_dp)
        hit_fun_limit = .false.
        hit_iter_limit = .false.

        do
            call setulb(n, maxcor_use, x, l, u, nbd, fval, g, factr, gtol_use, &
                wa, iwa, task, -1, csave, lsave, isave, dsave)
            if (task(1:2) == 'FG') then
                call evaluate_lbfgsb_fg(x, fval, g)
            else if (task(1:5) == 'NEW_X') then
                res%nit = res%nit + 1
                if (res%nit >= maxiter_use) then
                    hit_iter_limit = .true.
                    task = 'STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT'
                else if (res%nfev > maxfun_use) then
                    hit_fun_limit = .true.
                    task = 'STOP: TOTAL NO. OF F,G EVALUATIONS EXCEEDS LIMIT'
                end if
            else
                exit
            end if
        end do

        res%x = x
        res%fun = fval
        res%jac = g
        if (task(1:4) == 'CONV') then
            call finish_result(res, 0, .true., trim(task))
        else if (hit_fun_limit .or. hit_iter_limit .or. res%nfev > maxfun_use .or. &
                res%nit >= maxiter_use) then
            call finish_result(res, 1, .false., trim(task))
        else
            call finish_result(res, 2, .false., trim(task))
        end if

    contains

        subroutine evaluate_lbfgsb_fg(xv, fv, gv)
            real(dp), intent(in) :: xv(:) !! point requested by reverse communication
            real(dp), intent(out) :: fv !! objective value at xv
            real(dp), intent(out) :: gv(:) !! gradient at xv

            real(dp) :: fp
            real(dp) :: h
            real(dp), allocatable :: xp(:)
            integer :: k

            fv = fun(xv)
            res%nfev = res%nfev + 1
            if (present(jac)) then
                call invoke_gradient(jac, xv, gv)
                res%njev = res%njev + 1
                return
            end if

            allocate(xp(n))
            do k = 1, n
                call bounded_fd_step(xv(k), l(k), u(k), nbd(k), eps_use, h)
                if (same_real(h, 0.0_dp)) then
                    gv(k) = 0.0_dp
                else
                    xp = xv
                    xp(k) = xv(k) + h
                    fp = fun(xp)
                    res%nfev = res%nfev + 1
                    gv(k) = (fp - fv) / h
                end if
            end do
            res%njev = res%njev + 1
        end subroutine evaluate_lbfgsb_fg

    end function lbfgsb

    function slsqp(fun, x0, jac, lower, upper, constraints, constraint_jac, meq, &
            mineq, ftol, eps, maxiter) result(res)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        procedure(multivariate_gradient_function), optional :: jac !! analytic objective gradient
        real(dp), intent(in), optional :: lower(:) !! lower bounds; nonfinite values mean unbounded
        real(dp), intent(in), optional :: upper(:) !! upper bounds; nonfinite values mean unbounded
        procedure(multivariate_constraint_function), optional :: constraints !! combined constraint values
        procedure(multivariate_constraint_jacobian), optional :: constraint_jac !! combined constraint Jacobian
        integer, intent(in), optional :: meq !! number of equality constraints at the start of c
        integer, intent(in), optional :: mineq !! number of inequality constraints following equality constraints
        real(dp), intent(in), optional :: ftol !! SLSQP accuracy tolerance, default 1e-6
        real(dp), intent(in), optional :: eps !! absolute finite-difference step, default sqrt(epsilon)
        integer, intent(in), optional :: maxiter !! maximum SLSQP major iterations, default 100
        type(multivariate_minimize_result) :: res

        real(dp), allocatable :: cache_c(:)
        real(dp) :: cache_f
        real(dp), allocatable :: cache_g(:)
        real(dp), allocatable :: cache_x(:)
        real(dp) :: eps_use
        real(dp) :: ftol_use
        integer :: istat
        integer :: iterations
        integer :: maxiter_use
        integer :: meq_use
        integer :: mineq_use
        integer :: m
        integer :: n
        logical :: cache_c_valid
        logical :: cache_f_valid
        logical :: cache_g_valid
        logical :: ok
        type(slsqp_solver) :: solver
        character(len=96) :: final_message
        character(len=:), allocatable :: solver_message
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xl(:)
        real(dp), allocatable :: xu(:)

        n = size(x0)
        call init_result(res, 'slsqp', n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        meq_use = 0
        if (present(meq)) meq_use = meq
        mineq_use = 0
        if (present(mineq)) mineq_use = mineq
        if (meq_use < 0 .or. mineq_use < 0) return
        m = meq_use + mineq_use
        if (m > 0 .and. .not. present(constraints)) return
        if (m == 0 .and. present(constraints)) return
        if (present(constraint_jac) .and. .not. present(constraints)) return

        ftol_use = 1.0e-6_dp
        if (present(ftol)) ftol_use = ftol
        eps_use = default_fd_eps
        if (present(eps)) eps_use = eps
        maxiter_use = 100
        if (present(maxiter)) maxiter_use = maxiter
        if (.not. ieee_is_finite(ftol_use) .or. ftol_use < 0.0_dp) return
        if (.not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp) return
        if (maxiter_use < 0) return

        allocate(x(n), xl(n), xu(n), cache_x(n), cache_g(n), cache_c(m))
        x = x0
        call setup_slsqp_bounds(x, lower, upper, xl, xu, ok)
        if (.not. ok) return
        cache_f_valid = .false.
        cache_g_valid = .false.
        cache_c_valid = .false.

        call solver%initialize(n=n, m=m, meq=meq_use, max_iter=maxiter_use, &
            acc=ftol_use, f=eval_fc, g=eval_grad, xl=xl, xu=xu, status_ok=ok, &
            iprint=0)
        if (.not. ok) then
            call solver%destroy()
            return
        end if
        call solver%optimize(x, istat, iterations=iterations, status_message=solver_message)
        call solver%destroy()

        res%x = x
        res%nit = iterations
        ! The kernel increments its counter before checking the iteration limit.
        if (istat == 9) res%nit = min(iterations, maxiter_use)
        call final_values(x)
        if (allocated(solver_message)) then
            final_message = slsqp_status_message(istat, solver_message)
        else
            final_message = slsqp_status_message(istat)
        end if
        call finish_result(res, istat, istat == 0, trim(final_message))

    contains

        subroutine eval_fc(me, xv, fv, cv)
            class(slsqp_solver), intent(inout) :: me !! active SLSQP solver instance
            real(slsqp_wp), intent(in) :: xv(:) !! point requested by SLSQP
            real(slsqp_wp), intent(out) :: fv !! objective value at xv
            real(slsqp_wp), intent(out) :: cv(:) !! combined constraints at xv

            if (size(xv) < 0) call me%abort()
            fv = fun(xv)
            res%nfev = res%nfev + 1
            cache_x = xv
            cache_f = fv
            cache_f_valid = .true.
            cache_g_valid = .false.
            if (m > 0) then
                call invoke_constraints(constraints, xv, cv)
                cache_c = cv
                cache_c_valid = .true.
            else
                cache_c_valid = .true.
            end if
        end subroutine eval_fc

        subroutine eval_grad(me, xv, gv, av)
            class(slsqp_solver), intent(inout) :: me !! active SLSQP solver instance
            real(slsqp_wp), intent(in) :: xv(:) !! point requested by SLSQP
            real(slsqp_wp), intent(out) :: gv(:) !! objective gradient at xv
            real(slsqp_wp), intent(out) :: av(:, :) !! constraint Jacobian at xv

            if (size(xv) < 0) call me%abort()
            call objective_gradient(xv, gv)
            if (m > 0) call constraints_gradient(xv, av)
        end subroutine eval_grad

        subroutine objective_gradient(xv, gv)
            real(dp), intent(in) :: xv(:) !! point to differentiate
            real(dp), intent(out) :: gv(:) !! objective gradient at xv

            real(dp) :: f0
            real(dp) :: fp
            real(dp) :: h
            real(dp), allocatable :: xp(:)
            integer :: k

            if (present(jac)) then
                if (cache_f_valid .or. cache_c_valid) then
                    if (.not. same_vector(xv, cache_x)) then
                        cache_f_valid = .false.
                        cache_c_valid = .false.
                    end if
                end if
                call invoke_gradient(jac, xv, gv)
                res%njev = res%njev + 1
            else
                call get_cached_objective(xv, f0)
                allocate(xp(n))
                do k = 1, n
                    call slsqp_fd_step(xv(k), xl(k), xu(k), eps_use, h)
                    if (same_real(h, 0.0_dp)) then
                        gv(k) = 0.0_dp
                    else
                        xp = xv
                        xp(k) = xv(k) + h
                        fp = fun(xp)
                        res%nfev = res%nfev + 1
                        gv(k) = (fp - f0) / h
                    end if
                end do
                res%njev = res%njev + 1
            end if
            cache_x = xv
            cache_g = gv
            cache_g_valid = .true.
        end subroutine objective_gradient

        subroutine constraints_gradient(xv, av)
            real(dp), intent(in) :: xv(:) !! point to differentiate constraints
            real(dp), intent(out) :: av(:, :) !! constraint Jacobian at xv

            real(dp), allocatable :: c0(:)
            real(dp), allocatable :: cp(:)
            real(dp) :: h
            real(dp), allocatable :: xp(:)
            integer :: k

            if (present(constraint_jac)) then
                call invoke_constraint_jacobian(constraint_jac, xv, av)
                return
            end if
            allocate(c0(m), cp(m), xp(n))
            call get_cached_constraints(xv, c0)
            do k = 1, n
                call slsqp_fd_step(xv(k), xl(k), xu(k), eps_use, h)
                if (same_real(h, 0.0_dp)) then
                    av(:, k) = 0.0_dp
                else
                    xp = xv
                    xp(k) = xv(k) + h
                    call invoke_constraints(constraints, xp, cp)
                    av(:, k) = (cp - c0) / h
                end if
            end do
        end subroutine constraints_gradient

        subroutine get_cached_objective(xv, fv)
            real(dp), intent(in) :: xv(:) !! requested point
            real(dp), intent(out) :: fv !! objective value at requested point

            if (cache_f_valid .and. same_vector(xv, cache_x)) then
                fv = cache_f
            else
                fv = fun(xv)
                res%nfev = res%nfev + 1
                cache_x = xv
                cache_f = fv
                cache_f_valid = .true.
                cache_g_valid = .false.
                cache_c_valid = .false.
            end if
        end subroutine get_cached_objective

        subroutine get_cached_constraints(xv, cv)
            real(dp), intent(in) :: xv(:) !! requested point
            real(dp), intent(out) :: cv(:) !! combined constraints at requested point

            if (cache_c_valid .and. same_vector(xv, cache_x)) then
                cv = cache_c
            else
                call invoke_constraints(constraints, xv, cv)
                cache_x = xv
                cache_c = cv
                cache_c_valid = .true.
                cache_f_valid = .false.
                cache_g_valid = .false.
            end if
        end subroutine get_cached_constraints

        subroutine final_values(xv)
            real(dp), intent(in) :: xv(:) !! final SLSQP point

            real(dp), allocatable :: gv(:)

            allocate(gv(n))
            if (cache_f_valid .and. same_vector(xv, cache_x)) then
                res%fun = cache_f
            else
                res%fun = fun(xv)
                res%nfev = res%nfev + 1
            end if
            if (cache_g_valid .and. same_vector(xv, cache_x)) then
                gv = cache_g
            else
                call objective_gradient(xv, gv)
            end if
            res%jac = gv
        end subroutine final_values

    end function slsqp

    subroutine setup_lbfgsb_bounds(x, lower, upper, l, u, nbd, invalid)
        real(dp), intent(inout) :: x(:) !! starting point, clipped to finite bounds
        real(dp), intent(in), optional :: lower(:) !! optional lower bounds
        real(dp), intent(in), optional :: upper(:) !! optional upper bounds
        real(dp), intent(out) :: l(:) !! L-BFGS-B lower-bound array
        real(dp), intent(out) :: u(:) !! L-BFGS-B upper-bound array
        integer, intent(out) :: nbd(:) !! L-BFGS-B bound-type codes
        logical, intent(out) :: invalid !! true when bound shapes/order are invalid

        integer :: i
        logical :: has_lower
        logical :: has_upper

        invalid = .false.
        if (present(lower)) then
            if (size(lower) /= size(x)) then
                invalid = .true.
                return
            end if
        end if
        if (present(upper)) then
            if (size(upper) /= size(x)) then
                invalid = .true.
                return
            end if
        end if

        l = 0.0_dp
        u = 0.0_dp
        nbd = 0
        do i = 1, size(x)
            has_lower = .false.
            has_upper = .false.
            if (present(lower)) then
                if (ieee_is_finite(lower(i))) then
                    l(i) = lower(i)
                    has_lower = .true.
                end if
            end if
            if (present(upper)) then
                if (ieee_is_finite(upper(i))) then
                    u(i) = upper(i)
                    has_upper = .true.
                end if
            end if
            if (has_lower .and. has_upper) then
                if (l(i) > u(i)) then
                    invalid = .true.
                    return
                end if
                nbd(i) = 2
                x(i) = min(max(x(i), l(i)), u(i))
            else if (has_lower) then
                nbd(i) = 1
                x(i) = max(x(i), l(i))
            else if (has_upper) then
                nbd(i) = 3
                x(i) = min(x(i), u(i))
            end if
        end do
    end subroutine setup_lbfgsb_bounds

    subroutine three_point_gradient(fun, x, fx, step_scale, grad, nfev, lower, upper, nbd, rel_step)
        ! Centered or second-order one-sided objective differences. Reuse the
        ! known base value and track actual representable displacements.
        procedure(multivariate_objective_function) :: fun !! scalar objective callback
        real(dp), intent(in) :: x(:) !! differentiation point
        real(dp), intent(in) :: fx !! already evaluated objective at x
        real(dp), intent(in) :: step_scale !! positive finite-difference scale
        real(dp), intent(out) :: grad(:) !! numerical gradient
        integer, intent(out) :: nfev !! extra objective function evaluations
        real(dp), intent(in), optional :: lower(:) !! strict lower bounds
        real(dp), intent(in), optional :: upper(:) !! strict upper bounds
        integer, intent(in), optional :: nbd(:) !! bound codes
        real(dp), intent(in), optional :: rel_step(:) !! explicit relative steps for initial variables

        integer :: j
        real(dp) :: h1
        real(dp) :: h2
        real(dp) :: f1
        real(dp) :: f2
        real(dp) :: step
        real(dp) :: xp(size(x))

        grad = 0.0_dp
        nfev = 0
        do j = 1, size(x)
            step = step_scale * max(1.0_dp, abs(x(j)))
            if (present(rel_step)) then
                if (j <= size(rel_step)) then
                    step = abs(objective_relative_step(x(j), step_scale, rel_step(j)))
                end if
            end if
            if (present(nbd)) then
                call bounded_three_point_fd_steps(x(j), lower(j), upper(j), &
                    nbd(j), step, h1, h2)
            else
                h1 = -step
                h2 = step
            end if
            if (h1 == 0.0_dp .or. h2 == 0.0_dp .or. h1 == h2) cycle
            xp = x
            xp(j) = x(j) + h1
            if (present(nbd)) then
                if (nbd(j) == 1 .or. nbd(j) == 2) xp(j) = max(xp(j), lower(j))
                if (nbd(j) == 2 .or. nbd(j) == 3) xp(j) = min(xp(j), upper(j))
            end if
            h1 = xp(j) - x(j)
            if (h1 == 0.0_dp) cycle
            f1 = fun(xp)
            nfev = nfev + 1
            xp(j) = x(j) + h2
            if (present(nbd)) then
                if (nbd(j) == 1 .or. nbd(j) == 2) xp(j) = max(xp(j), lower(j))
                if (nbd(j) == 2 .or. nbd(j) == 3) xp(j) = min(xp(j), upper(j))
            end if
            h2 = xp(j) - x(j)
            if (h2 == 0.0_dp .or. h2 == h1) then
                grad(j) = (f1 - fx) / h1
                cycle
            end if
            f2 = fun(xp)
            nfev = nfev + 1
            grad(j) = (h2 * ((f1 - fx) / h1) - &
                h1 * ((f2 - fx) / h2)) / (h2 - h1)
        end do
    end subroutine three_point_gradient

    pure function objective_relative_step(x, default_scale, rel_step) result(h)
        ! Shared objective/constraint implementation of SciPy's
        ! _compute_absolute_step convention: explicit relative steps multiply
        ! abs(x), not max(1,abs(x)); their sign determines 2-point direction.
        ! For exact zero or sub-ULP steps, retain the baseline default scale.
        real(dp), intent(in) :: x !! physical coordinate
        real(dp), intent(in) :: default_scale !! fallback relative difference scale
        real(dp), intent(in) :: rel_step !! user relative step, zero permitted
        real(dp) :: h !! signed representable-step request

        h = rel_step * sign(1.0_dp, x) * abs(x)
        if (x + h == x) h = default_scale * sign(1.0_dp, x) * max(1.0_dp, abs(x))
    end function objective_relative_step

    subroutine bounded_three_point_fd_steps(x, lower, upper, nbd, eps, h1, h2)
        real(dp), intent(in) :: x !! coordinate to differentiate
        real(dp), intent(in) :: lower !! lower bound when nbd is 1 or 2
        real(dp), intent(in) :: upper !! upper bound when nbd is 2 or 3
        integer, intent(in) :: nbd !! box-bound type
        real(dp), intent(in) :: eps !! requested positive absolute step
        real(dp), intent(out) :: h1 !! first signed displacement
        real(dp), intent(out) :: h2 !! second signed displacement

        real(dp) :: backward_room
        real(dp) :: forward_room
        real(dp) :: h

        backward_room = huge(1.0_dp)
        forward_room = huge(1.0_dp)
        if (nbd == 1 .or. nbd == 2) backward_room = max(0.0_dp, x - lower)
        if (nbd == 2 .or. nbd == 3) forward_room = max(0.0_dp, upper - x)

        if (backward_room >= eps .and. forward_room >= eps) then
            ! Centered second-order stencil.
            h1 = -eps
            h2 = eps
        else if (forward_room >= backward_room .and. forward_room > 0.0_dp) then
            ! Second-order one-sided stencil (or shortened at narrow bounds).
            h = min(eps, forward_room / 2.0_dp)
            h1 = h
            h2 = 2.0_dp * h
        else if (backward_room > 0.0_dp) then
            h = min(eps, backward_room / 2.0_dp)
            h1 = -h
            h2 = -2.0_dp * h
        else
            ! Both bounds coincide, leaving no feasible displacement.
            h1 = 0.0_dp
            h2 = 0.0_dp
        end if
    end subroutine bounded_three_point_fd_steps

    subroutine bounded_signed_fd_step(x, lower, upper, nbd, requested, h)
        ! Two-point objective/constraint derivative stencils use SciPy's
        ! signed perturbation convention. Reverse only when bounds require it.
        real(dp), intent(in) :: x !! base coordinate
        real(dp), intent(in) :: lower !! lower endpoint when bounded
        real(dp), intent(in) :: upper !! upper endpoint when bounded
        integer, intent(in) :: nbd !! 0 none, 1 lower, 2 both, 3 upper
        real(dp), intent(in) :: requested !! signed desired displacement
        real(dp), intent(out) :: h !! signed feasible displacement

        real(dp) :: backward_room
        real(dp) :: forward_room
        real(dp) :: magnitude

        backward_room = huge(1.0_dp)
        forward_room = huge(1.0_dp)
        if (nbd == 1 .or. nbd == 2) backward_room = max(0.0_dp, x - lower)
        if (nbd == 2 .or. nbd == 3) forward_room = max(0.0_dp, upper - x)
        magnitude = abs(requested)
        if (requested < 0.0_dp) then
            if (backward_room >= magnitude) then
                h = -magnitude
            else if (forward_room >= magnitude) then
                h = magnitude
            else if (backward_room >= forward_room) then
                h = -backward_room
            else
                h = forward_room
            end if
        else
            if (forward_room >= magnitude) then
                h = magnitude
            else if (backward_room >= magnitude) then
                h = -magnitude
            else if (forward_room >= backward_room) then
                h = forward_room
            else
                h = -backward_room
            end if
        end if
    end subroutine bounded_signed_fd_step

    subroutine bounded_fd_step(x, lower, upper, nbd, eps, h)
        real(dp), intent(in) :: x !! scalar coordinate at which to differentiate
        real(dp), intent(in) :: lower !! finite lower bound when nbd is 1 or 2
        real(dp), intent(in) :: upper !! finite upper bound when nbd is 2 or 3
        integer, intent(in) :: nbd !! L-BFGS-B bound-type code
        real(dp), intent(in) :: eps !! requested absolute finite-difference step
        real(dp), intent(out) :: h !! feasible signed finite-difference step

        real(dp) :: backward_room
        real(dp) :: forward_room

        forward_room = huge(1.0_dp)
        backward_room = huge(1.0_dp)
        if (nbd == 2 .or. nbd == 3) forward_room = max(0.0_dp, upper - x)
        if (nbd == 1 .or. nbd == 2) backward_room = max(0.0_dp, x - lower)

        if (forward_room >= eps) then
            h = eps
        else if (backward_room >= eps) then
            h = -eps
        else if (forward_room >= backward_room .and. forward_room > 0.0_dp) then
            h = forward_room
        else if (backward_room > 0.0_dp) then
            h = -backward_room
        else
            h = 0.0_dp
        end if
    end subroutine bounded_fd_step

    subroutine setup_slsqp_bounds(x, lower, upper, xl, xu, ok)
        real(dp), intent(inout) :: x(:) !! starting point, clipped to finite bounds
        real(dp), intent(in), optional :: lower(:) !! optional lower bounds
        real(dp), intent(in), optional :: upper(:) !! optional upper bounds
        real(dp), intent(out) :: xl(:) !! SLSQP lower bounds; NaN means unbounded
        real(dp), intent(out) :: xu(:) !! SLSQP upper bounds; NaN means unbounded
        logical, intent(out) :: ok !! true when bound shapes/order are valid

        real(dp) :: nan_value
        integer :: i
        logical :: has_lower
        logical :: has_upper

        ok = .false.
        if (present(lower)) then
            if (size(lower) /= size(x)) return
        end if
        if (present(upper)) then
            if (size(upper) /= size(x)) return
        end if

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        xl = nan_value
        xu = nan_value
        do i = 1, size(x)
            has_lower = .false.
            has_upper = .false.
            if (present(lower)) then
                if (ieee_is_finite(lower(i))) then
                    xl(i) = lower(i)
                    has_lower = .true.
                end if
            end if
            if (present(upper)) then
                if (ieee_is_finite(upper(i))) then
                    xu(i) = upper(i)
                    has_upper = .true.
                end if
            end if
            if (has_lower .and. has_upper) then
                if (xl(i) > xu(i)) return
                x(i) = min(max(x(i), xl(i)), xu(i))
            else if (has_lower) then
                x(i) = max(x(i), xl(i))
            else if (has_upper) then
                x(i) = min(x(i), xu(i))
            end if
        end do
        ok = .true.
    end subroutine setup_slsqp_bounds

    subroutine slsqp_fd_step(x, lower, upper, eps, h)
        real(dp), intent(in) :: x !! scalar coordinate at which to differentiate
        real(dp), intent(in) :: lower !! lower bound or NaN when unbounded
        real(dp), intent(in) :: upper !! upper bound or NaN when unbounded
        real(dp), intent(in) :: eps !! requested absolute finite-difference step
        real(dp), intent(out) :: h !! feasible signed finite-difference step

        real(dp) :: backward_room
        real(dp) :: forward_room

        forward_room = huge(1.0_dp)
        backward_room = huge(1.0_dp)
        if (ieee_is_finite(upper)) forward_room = max(0.0_dp, upper - x)
        if (ieee_is_finite(lower)) backward_room = max(0.0_dp, x - lower)

        if (forward_room >= eps) then
            h = eps
        else if (backward_room >= eps) then
            h = -eps
        else if (forward_room >= backward_room .and. forward_room > 0.0_dp) then
            h = forward_room
        else if (backward_room > 0.0_dp) then
            h = -backward_room
        else
            h = 0.0_dp
        end if
    end subroutine slsqp_fd_step

    function slsqp_status_message(status, upstream_message) result(message)
        integer, intent(in) :: status !! SLSQP termination status
        character(len=*), intent(in), optional :: upstream_message !! backend status text
        character(len=96) :: message

        if (present(upstream_message)) then
            if (len_trim(upstream_message) > 0) then
                message = trim(upstream_message)
                return
            end if
        end if

        select case (status)
        case (-1)
            message = 'gradient evaluation required'
        case (0)
            message = 'optimization terminated successfully'
        case (2)
            message = 'more equality constraints than independent variables'
        case (3)
            message = 'more than 3*n iterations in LSQ subproblem'
        case (4)
            message = 'inequality constraints incompatible'
        case (5)
            message = 'singular matrix E in LSQ subproblem'
        case (6)
            message = 'singular matrix C in LSQ subproblem'
        case (7)
            message = 'rank-deficient equality constraint subproblem HFTI'
        case (8)
            message = 'positive directional derivative for linesearch'
        case (9)
            message = 'iteration limit reached'
        case default
            write(message, '(a,i0)') 'SLSQP terminated with status ', status
        end select
    end function slsqp_status_message

    subroutine setup_box_bounds(n, lower, upper, lb, ub, nbd, invalid)
        integer, intent(in) :: n !! problem dimension
        real(dp), intent(in), optional :: lower(:) !! optional lower bounds
        real(dp), intent(in), optional :: upper(:) !! optional upper bounds
        real(dp), intent(out) :: lb(:) !! finite lower-bound values when present
        real(dp), intent(out) :: ub(:) !! finite upper-bound values when present
        integer, intent(out) :: nbd(:) !! 0 none, 1 lower, 2 both, 3 upper
        logical, intent(out) :: invalid !! true when shapes/order are invalid

        integer :: i
        logical :: has_lower
        logical :: has_upper

        invalid = .false.
        if (size(lb) /= n .or. size(ub) /= n .or. size(nbd) /= n) then
            invalid = .true.
            return
        end if
        if (present(lower)) then
            if (size(lower) /= n) then
                invalid = .true.
                return
            end if
        end if
        if (present(upper)) then
            if (size(upper) /= n) then
                invalid = .true.
                return
            end if
        end if

        lb = 0.0_dp
        ub = 0.0_dp
        nbd = 0
        do i = 1, n
            has_lower = .false.
            has_upper = .false.
            if (present(lower)) then
                if (ieee_is_finite(lower(i))) then
                    lb(i) = lower(i)
                    has_lower = .true.
                end if
            end if
            if (present(upper)) then
                if (ieee_is_finite(upper(i))) then
                    ub(i) = upper(i)
                    has_upper = .true.
                end if
            end if
            if (has_lower .and. has_upper) then
                if (lb(i) > ub(i)) then
                    invalid = .true.
                    return
                end if
                nbd(i) = 2
            else if (has_lower) then
                nbd(i) = 1
            else if (has_upper) then
                nbd(i) = 3
            end if
        end do
    end subroutine setup_box_bounds

    subroutine clip_to_box(x, lb, ub, nbd)
        real(dp), intent(inout) :: x(:) !! vector clipped in place
        real(dp), intent(in) :: lb(:) !! lower-bound values
        real(dp), intent(in) :: ub(:) !! upper-bound values
        integer, intent(in) :: nbd(:) !! bound-type codes

        integer :: i

        do i = 1, size(x)
            if (nbd(i) == 1 .or. nbd(i) == 2) x(i) = max(x(i), lb(i))
            if (nbd(i) == 2 .or. nbd(i) == 3) x(i) = min(x(i), ub(i))
        end do
    end subroutine clip_to_box

    logical function point_in_box(x, lb, ub, nbd) result(ok)
        real(dp), intent(in) :: x(:) !! point to test
        real(dp), intent(in) :: lb(:) !! lower-bound values
        real(dp), intent(in) :: ub(:) !! upper-bound values
        integer, intent(in) :: nbd(:) !! bound-type codes

        integer :: i

        ok = .true.
        do i = 1, size(x)
            if ((nbd(i) == 1 .or. nbd(i) == 2) .and. x(i) < lb(i)) then
                ok = .false.
                return
            end if
            if ((nbd(i) == 2 .or. nbd(i) == 3) .and. x(i) > ub(i)) then
                ok = .false.
                return
            end if
        end do
    end function point_in_box

    subroutine setup_tnc_scaling(x, lb, ub, nbd, scale, offset)
        real(dp), intent(in) :: x(:) !! feasible starting point
        real(dp), intent(in) :: lb(:) !! finite lower bounds when present
        real(dp), intent(in) :: ub(:) !! finite upper bounds when present
        integer, intent(in) :: nbd(:) !! bound-type codes
        real(dp), intent(out) :: scale(:) !! TNC variable scales
        real(dp), intent(out) :: offset(:) !! TNC variable offsets

        integer :: i

        do i = 1, size(x)
            if (nbd(i) == 2) then
                scale(i) = ub(i) - lb(i)
                offset(i) = 0.5_dp * (ub(i) + lb(i))
            else
                scale(i) = 1.0_dp + abs(x(i))
                offset(i) = x(i)
            end if
            if (scale(i) <= 0.0_dp .or. .not. ieee_is_finite(scale(i))) scale(i) = 1.0_dp
        end do
    end subroutine setup_tnc_scaling

    subroutine projected_scaled_gradient(x, g_scaled, scale, lb, ub, nbd, pg)
        real(dp), intent(in) :: x(:) !! current point in original coordinates
        real(dp), intent(in) :: g_scaled(:) !! gradient with respect to scaled coordinates
        real(dp), intent(in) :: scale(:) !! positive variable scales
        real(dp), intent(in) :: lb(:) !! finite lower bounds when present
        real(dp), intent(in) :: ub(:) !! finite upper bounds when present
        integer, intent(in) :: nbd(:) !! bound-type codes
        real(dp), intent(out) :: pg(:) !! scaled projected gradient

        real(dp) :: bound_tol
        integer :: i

        pg = g_scaled
        do i = 1, size(x)
            bound_tol = 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)), scale(i))
            if ((nbd(i) == 1 .or. nbd(i) == 2) .and. &
                    x(i) <= lb(i) + bound_tol .and. g_scaled(i) > 0.0_dp) pg(i) = 0.0_dp
            if ((nbd(i) == 2 .or. nbd(i) == 3) .and. &
                    x(i) >= ub(i) - bound_tol .and. g_scaled(i) < 0.0_dp) pg(i) = 0.0_dp
        end do
    end subroutine projected_scaled_gradient

    function same_active_component(x, g_scaled, scale, lb, ub, nbd) result(active)
        real(dp), intent(in) :: x(:) !! current point in original coordinates
        real(dp), intent(in) :: g_scaled(:) !! gradient with respect to scaled coordinates
        real(dp), intent(in) :: scale(:) !! positive variable scales
        real(dp), intent(in) :: lb(:) !! finite lower bounds when present
        real(dp), intent(in) :: ub(:) !! finite upper bounds when present
        integer, intent(in) :: nbd(:) !! bound-type codes
        logical :: active(size(x))

        real(dp) :: bound_tol
        integer :: i

        active = .false.
        do i = 1, size(x)
            if (nbd(i) == 2 .and. same_real(lb(i), ub(i))) then
                active(i) = .true.
                cycle
            end if
            bound_tol = 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)), scale(i))
            if ((nbd(i) == 1 .or. nbd(i) == 2) .and. &
                    x(i) <= lb(i) + bound_tol .and. g_scaled(i) > 0.0_dp) active(i) = .true.
            if ((nbd(i) == 2 .or. nbd(i) == 3) .and. &
                    x(i) >= ub(i) - bound_tol .and. g_scaled(i) < 0.0_dp) active(i) = .true.
        end do
    end function same_active_component

    real(dp) function scaled_step_to_bounds(x, direction, scale, lb, ub, nbd) result(alpha_max)
        real(dp), intent(in) :: x(:) !! current point in original coordinates
        real(dp), intent(in) :: direction(:) !! direction in scaled coordinates
        real(dp), intent(in) :: scale(:) !! positive variable scales
        real(dp), intent(in) :: lb(:) !! finite lower bounds when present
        real(dp), intent(in) :: ub(:) !! finite upper bounds when present
        integer, intent(in) :: nbd(:) !! bound-type codes

        real(dp) :: dx
        real(dp) :: candidate
        integer :: i

        alpha_max = huge(1.0_dp)
        do i = 1, size(x)
            dx = scale(i) * direction(i)
            if (dx > 0.0_dp .and. (nbd(i) == 2 .or. nbd(i) == 3)) then
                candidate = max(0.0_dp, (ub(i) - x(i)) / dx)
                alpha_max = min(alpha_max, candidate)
            else if (dx < 0.0_dp .and. (nbd(i) == 1 .or. nbd(i) == 2)) then
                candidate = max(0.0_dp, (lb(i) - x(i)) / dx)
                alpha_max = min(alpha_max, candidate)
            end if
        end do
    end function scaled_step_to_bounds

    subroutine powell_line_bounds(x, direction, lb, ub, nbd, amin, amax, &
            lower_finite, upper_finite)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: direction(:) !! search direction
        real(dp), intent(in) :: lb(:) !! lower-bound values
        real(dp), intent(in) :: ub(:) !! upper-bound values
        integer, intent(in) :: nbd(:) !! bound-type codes
        real(dp), intent(out) :: amin !! lower scalar line bound when finite
        real(dp), intent(out) :: amax !! upper scalar line bound when finite
        logical, intent(out) :: lower_finite !! true when amin is finite
        logical, intent(out) :: upper_finite !! true when amax is finite

        real(dp) :: candidate
        integer :: i

        amin = 0.0_dp
        amax = 0.0_dp
        lower_finite = .false.
        upper_finite = .false.

        do i = 1, size(x)
            if (direction(i) == 0.0_dp) cycle
            if (nbd(i) == 1 .or. nbd(i) == 2) then
                candidate = (lb(i) - x(i)) / direction(i)
                if (direction(i) > 0.0_dp) then
                    if (.not. lower_finite .or. candidate > amin) amin = candidate
                    lower_finite = .true.
                else
                    if (.not. upper_finite .or. candidate < amax) amax = candidate
                    upper_finite = .true.
                end if
            end if
            if (nbd(i) == 2 .or. nbd(i) == 3) then
                candidate = (ub(i) - x(i)) / direction(i)
                if (direction(i) > 0.0_dp) then
                    if (.not. upper_finite .or. candidate < amax) amax = candidate
                    upper_finite = .true.
                else
                    if (.not. lower_finite .or. candidate > amin) amin = candidate
                    lower_finite = .true.
                end if
            end if
        end do

        if (lower_finite .and. upper_finite) then
            if (amax < amin) then
                amin = 0.0_dp
                amax = 0.0_dp
            end if
        end if
    end subroutine powell_line_bounds

    subroutine powell_extrapolated_point(x, direction, lb, ub, nbd, x2)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: direction(:) !! extrapolation direction
        real(dp), intent(in) :: lb(:) !! lower-bound values
        real(dp), intent(in) :: ub(:) !! upper-bound values
        integer, intent(in) :: nbd(:) !! bound-type codes
        real(dp), intent(out) :: x2(:) !! feasible extrapolated point where possible

        real(dp) :: amin
        real(dp) :: amax
        real(dp) :: alpha
        logical :: lower_finite
        logical :: upper_finite

        if (.not. any(nbd /= 0)) then
            x2 = x + direction
            return
        end if
        call powell_line_bounds(x, direction, lb, ub, nbd, amin, amax, &
            lower_finite, upper_finite)
        alpha = 1.0_dp
        if (upper_finite) alpha = min(alpha, amax)
        x2 = x + alpha * direction
    end subroutine powell_extrapolated_point

    subroutine powell_line_search(fun, x, direction, tol, maxfev, nfev, fval, ok, &
            lb, ub, nbd)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(inout) :: x(:) !! line-search base and returned point
        real(dp), intent(inout) :: direction(:) !! input direction and returned displacement
        real(dp), intent(in) :: tol !! Brent relative tolerance
        integer, intent(in) :: maxfev !! outer objective-evaluation cap
        integer, intent(inout) :: nfev !! outer objective-evaluation count
        real(dp), intent(inout) :: fval !! objective at x and returned point
        logical, intent(out) :: ok !! true when the line search produced a point
        real(dp), intent(in) :: lb(:) !! lower-bound values
        real(dp), intent(in) :: ub(:) !! upper-bound values
        integer, intent(in) :: nbd(:) !! bound-type codes

        type(minimize_result) :: line_res
        real(dp), allocatable :: base(:)
        real(dp), allocatable :: idir(:)
        real(dp) :: alpha
        real(dp) :: amin
        real(dp) :: amax
        real(dp) :: angle_min
        real(dp) :: angle_max
        real(dp), parameter :: half_pi = 0.5_dp * acos(-1.0_dp)
        logical :: exhausted
        logical :: lower_finite
        logical :: upper_finite

        ok = .true.
        if (.not. any(abs(direction) > 0.0_dp)) return
        if (nfev >= maxfev) then
            ok = .false.
            return
        end if
        allocate(base(size(x)), idir(size(x)))
        base = x
        idir = direction
        exhausted = .false.

        if (.not. any(nbd /= 0)) then
            line_res = brent(line_objective, xtol=tol)
            alpha = line_res%x
        else
            call powell_line_bounds(base, idir, lb, ub, nbd, amin, amax, &
                lower_finite, upper_finite)
            if (.not. lower_finite .and. .not. upper_finite) then
                line_res = brent(line_objective, xtol=tol)
                alpha = line_res%x
            else if (lower_finite .and. upper_finite) then
                line_res = bounded(line_objective, amin, amax, xatol=tol / 100.0_dp)
                alpha = line_res%x
            else
                angle_min = -half_pi
                angle_max = half_pi
                if (lower_finite) angle_min = atan(amin)
                if (upper_finite) angle_max = atan(amax)
                line_res = bounded(angle_objective, angle_min, angle_max, &
                    xatol=tol / 100.0_dp)
                alpha = tan(line_res%x)
            end if
        end if

        if (exhausted) then
            ok = .false.
            return
        end if
        direction = alpha * idir
        x = base + direction
        fval = line_res%fun

    contains

        function line_objective(a) result(value)
            real(dp), intent(in) :: a !! scalar line parameter
            real(dp) :: value !! objective along the Powell direction
            real(dp) :: xt(size(base))

            if (nfev >= maxfev) then
                exhausted = .true.
                value = huge(1.0_dp)
                return
            end if
            xt = base + a * idir
            value = fun(xt)
            nfev = nfev + 1
        end function line_objective

        function angle_objective(theta) result(value)
            real(dp), intent(in) :: theta !! transformed finite-angle parameter
            real(dp) :: value !! objective along the transformed Powell line

            value = line_objective(tan(theta))
        end function angle_objective

    end subroutine powell_line_search

    subroutine resolve_limits(n, multiplier, maxiter, maxfev, maxiter_use, maxfev_use)
        integer, intent(in) :: n !! problem dimension
        integer, intent(in) :: multiplier !! default per-dimension limit multiplier
        integer, intent(in), optional :: maxiter !! user iteration limit
        integer, intent(in), optional :: maxfev !! user objective-evaluation limit
        integer, intent(out) :: maxiter_use !! resolved iteration limit
        integer, intent(out) :: maxfev_use !! resolved evaluation limit

        integer, parameter :: unlimited = 1000000000

        if (.not. present(maxiter) .and. .not. present(maxfev)) then
            maxiter_use = n * multiplier
            maxfev_use = n * multiplier
        else if (.not. present(maxiter)) then
            maxiter_use = unlimited
            maxfev_use = maxfev
        else if (.not. present(maxfev)) then
            maxiter_use = maxiter
            maxfev_use = unlimited
        else
            maxiter_use = maxiter
            maxfev_use = maxfev
        end if
    end subroutine resolve_limits

    subroutine gradient_options(n, gtol, eps, maxiter, c1, c2, gtol_use, eps_use, &
            maxiter_use, c1_use, c2_use, default_c2)
        integer, intent(in) :: n !! problem dimension
        real(dp), intent(in), optional :: gtol !! requested gradient tolerance
        real(dp), intent(in), optional :: eps !! requested finite-difference step
        integer, intent(in), optional :: maxiter !! requested iteration limit
        real(dp), intent(in), optional :: c1 !! requested Armijo constant
        real(dp), intent(in), optional :: c2 !! requested curvature constant
        real(dp), intent(out) :: gtol_use !! resolved gradient tolerance
        real(dp), intent(out) :: eps_use !! resolved finite-difference step
        integer, intent(out) :: maxiter_use !! resolved iteration limit
        real(dp), intent(out) :: c1_use !! resolved Armijo constant
        real(dp), intent(out) :: c2_use !! resolved curvature constant
        real(dp), intent(in), optional :: default_c2 !! method-specific c2 default

        gtol_use = default_gtol
        if (present(gtol)) gtol_use = gtol
        eps_use = default_fd_eps
        if (present(eps)) eps_use = eps
        maxiter_use = n * 200
        if (present(maxiter)) maxiter_use = maxiter
        c1_use = 1.0e-4_dp
        if (present(c1)) c1_use = c1
        c2_use = 0.9_dp
        if (present(default_c2)) c2_use = default_c2
        if (present(c2)) c2_use = c2
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) maxiter_use = -1
        if (.not. ieee_is_finite(eps_use) .or. eps_use <= 0.0_dp) maxiter_use = -1
        if (.not. ieee_is_finite(c1_use) .or. .not. ieee_is_finite(c2_use)) maxiter_use = -1
        if (c1_use <= 0.0_dp .or. c1_use >= c2_use .or. c2_use >= 1.0_dp) maxiter_use = -1
    end subroutine gradient_options

    logical function eval_nm(fun, x, maxfev, nfev, fx) result(ok)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x(:) !! point to evaluate
        integer, intent(in) :: maxfev !! objective-evaluation cap
        integer, intent(inout) :: nfev !! evaluation counter
        real(dp), intent(out) :: fx !! objective value

        if (nfev >= maxfev) then
            ok = .false.
            fx = huge(1.0_dp)
            return
        end if
        fx = fun(x)
        nfev = nfev + 1
        ok = .true.
    end function eval_nm

    logical function eval_powell(fun, x, maxfev, nfev, fx) result(ok)
        procedure(multivariate_objective_function) :: fun !! scalar objective function
        real(dp), intent(in) :: x(:) !! point to evaluate
        integer, intent(in) :: maxfev !! objective-evaluation cap
        integer, intent(inout) :: nfev !! evaluation counter
        real(dp), intent(out) :: fx !! objective value

        ok = eval_nm(fun, x, maxfev, nfev, fx)
    end function eval_powell

    subroutine sort_simplex(sim, fsim)
        real(dp), intent(inout) :: sim(:, :) !! simplex rows to sort by objective
        real(dp), intent(inout) :: fsim(:) !! objective values corresponding to rows

        real(dp), allocatable :: row(:)
        real(dp) :: ftmp
        integer :: i
        integer :: j
        integer :: k

        allocate(row(size(sim, 2)))
        do i = 1, size(fsim) - 1
            k = i
            do j = i + 1, size(fsim)
                if (fsim(j) < fsim(k)) k = j
            end do
            if (k /= i) then
                ftmp = fsim(i)
                fsim(i) = fsim(k)
                fsim(k) = ftmp
                row = sim(i, :)
                sim(i, :) = sim(k, :)
                sim(k, :) = row
            end if
        end do
    end subroutine sort_simplex

    real(dp) function max_simplex_dx(sim) result(value)
        real(dp), intent(in) :: sim(:, :) !! simplex rows sorted best first

        integer :: j

        value = 0.0_dp
        do j = 2, size(sim, 1)
            value = max(value, maxval(abs(sim(j, :) - sim(1, :))))
        end do
    end function max_simplex_dx

    subroutine identity_minus_outer(a, b, scale, matrix)
        real(dp), intent(in) :: a(:) !! left outer-product vector
        real(dp), intent(in) :: b(:) !! right outer-product vector
        real(dp), intent(in) :: scale !! outer-product multiplier
        real(dp), intent(out) :: matrix(:, :) !! identity minus scaled outer product

        integer :: i
        integer :: j

        do j = 1, size(a)
            do i = 1, size(a)
                matrix(i, j) = -scale * a(i) * b(j)
            end do
            matrix(j, j) = matrix(j, j) + 1.0_dp
        end do
    end subroutine identity_minus_outer

    function outer_product(a, b) result(matrix)
        real(dp), intent(in) :: a(:) !! left vector
        real(dp), intent(in) :: b(:) !! right vector
        real(dp) :: matrix(size(a), size(b))

        integer :: i
        integer :: j

        do j = 1, size(b)
            do i = 1, size(a)
                matrix(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product

    logical function same_vector(a, b) result(equal)
        real(dp), intent(in) :: a(:) !! first vector
        real(dp), intent(in) :: b(:) !! second vector

        equal = size(a) == size(b)
        if (equal) equal = all(a <= b .and. a >= b)
    end function same_vector

    pure logical function same_real(a, b) result(equal)
        real(dp), intent(in) :: a !! first value to compare exactly
        real(dp), intent(in) :: b !! second value to compare exactly

        equal = a <= b .and. a >= b
    end function same_real

    subroutine init_result(res, method, n)
        type(multivariate_minimize_result), intent(out) :: res !! result to initialize
        character(len=*), intent(in) :: method !! selected solver method
        integer, intent(in) :: n !! parameter dimension

        allocate(res%x(max(0, n)), res%lagrangian_grad(max(0, n)))
        res%x = 0.0_dp
        res%lagrangian_grad = 0.0_dp
        res%fun = ieee_value(0.0_dp, ieee_quiet_nan)
        res%nit = 0
        res%nfev = 0
        res%njev = 0
        res%nhev = 0
        res%cg_niter = 0
        res%cg_stop_cond = 0
        res%maxcv = 0.0_dp
        res%constr_violation = 0.0_dp
        res%optimality = 0.0_dp
        res%constr_penalty = 0.0_dp
        res%tr_radius = 0.0_dp
        res%barrier_parameter = 0.0_dp
        res%barrier_tolerance = 0.0_dp
        res%status = -1
        res%success = .false.
        res%message = 'invalid input'
        res%method = method
    end subroutine init_result

    subroutine finish_result(res, status, success, message)
        type(multivariate_minimize_result), intent(inout) :: res !! result to finalize
        integer, intent(in) :: status !! method-specific termination code
        logical, intent(in) :: success !! true when optimization converged
        character(len=*), intent(in) :: message !! termination description

        res%status = status
        res%success = success
        res%message = message
    end subroutine finish_result

    pure function lowercase(text) result(out)
        character(len=*), intent(in) :: text !! text to convert to lowercase
        character(len=len(text)) :: out

        integer :: i
        integer :: code

        out = text
        do i = 1, len(text)
            code = iachar(out(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
        end do
    end function lowercase

    ! Dispatch optional callbacks outside internal procedures to retain explicit
    ! interfaces on compilers that diagnose host-associated optional calls.

    subroutine invoke_gradient(proc, x, grad)
        implicit none
        procedure(multivariate_gradient_function) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: grad(:) !! objective gradient at x

        call proc(x, grad)
    end subroutine invoke_gradient

    subroutine invoke_hessian(proc, x, hess)
        implicit none
        procedure(multivariate_hessian_function) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: hess(:, :) !! objective Hessian at x

        call proc(x, hess)
    end subroutine invoke_hessian

    subroutine invoke_hessian_product(proc, x, p, hp)
        implicit none
        procedure(multivariate_hessian_product_function) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: p(:) !! vector multiplied by the Hessian
        real(dp), intent(out) :: hp(:) !! Hessian-vector product

        call proc(x, p, hp)
    end subroutine invoke_hessian_product

    subroutine invoke_constraints(proc, x, c)
        implicit none
        procedure(multivariate_constraint_function) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! combined constraint values at x

        call proc(x, c)
    end subroutine invoke_constraints

    subroutine invoke_constraint_jacobian(proc, x, a)
        implicit none
        procedure(multivariate_constraint_jacobian) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: a(:, :) !! constraint Jacobian at x

        call proc(x, a)
    end subroutine invoke_constraint_jacobian

    subroutine invoke_constraint_hessian(proc, x, v, hess)
        implicit none
        procedure(multivariate_constraint_hessian) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(in) :: v(:) !! constraint multipliers
        real(dp), intent(out) :: hess(:, :) !! weighted constraint Hessian at x

        call proc(x, v, hess)
    end subroutine invoke_constraint_hessian

    subroutine invoke_iteration_callback(proc, x, nit, stop)
        implicit none
        procedure(trust_constr_iteration_callback) :: proc !! supplied callback
        real(dp), intent(in) :: x(:) !! current parameter vector
        integer, intent(in) :: nit !! current iteration number
        logical, intent(out) :: stop !! true to request termination

        call proc(x, nit, stop)
    end subroutine invoke_iteration_callback

end module scifort_optimize_minimize
