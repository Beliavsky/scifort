! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Multidimensional nonlinear-equation solving for the scipy.optimize translation.
!
! The MINPACK hybrid and Levenberg-Marquardt numerical kernels are provided by
! the community-maintained fortran-lang/minpack FPM dependency. This module is
! an original SciFort adapter that maps SciPy 1.17 root/fsolve defaults and
! result bookkeeping onto those public MINPACK interfaces.

module scifort_optimize_root
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use minpack_module, only : hybrd, hybrj, lmder, lmdif
    use scifort_optimize_df_sane, only : solve_df_sane
    use scifort_optimize_nonlin, only : solve_krylov_root, solve_nonlin_root, &
        vector_root_norm_function => nonlin_norm_function, &
        vector_root_linear_operator => nonlin_linear_operator, &
        vector_root_preconditioner_update => nonlin_preconditioner_update, &
        vector_root_inner_solver => nonlin_inner_solver, &
        vector_root_jacobian_product => nonlin_jacobian_product, &
        vector_root_base_diagonal => nonlin_base_diagonal
    use scifort_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: default_xtol = 1.49012e-8_dp
    real(dp), parameter :: default_ftol = 1.49012e-8_dp
    real(dp), parameter :: default_gtol = 0.0_dp
    real(dp), parameter :: default_factor = 100.0_dp

    type, public :: vector_root_result
        real(dp), allocatable :: x(:) !! final estimate of the solution vector
        real(dp), allocatable :: fun(:) !! residual vector evaluated at x
        real(dp), allocatable :: fjac(:, :) !! final MINPACK QR/Jacobian workspace
        real(dp), allocatable :: r(:) !! packed upper-triangular R for hybr methods
        real(dp), allocatable :: qtf(:) !! transformed residual vector Q^T f
        integer, allocatable :: ipvt(:) !! one-based MINPACK column pivots for lm
        integer :: nfev = 0 !! objective evaluations, including SciPy-style preflight calls
        integer :: njev = 0 !! Jacobian evaluations reported by MINPACK
        integer :: nit = 0 !! nonlinear iteration/check count when available
        integer :: status = 0 !! solver status; raw MINPACK status for hybr/lm
        logical :: success = .false. !! true for a SciPy-success termination status
        character(len=16) :: method = '' !! solver method actually used
        character(len=160) :: message = 'improper input parameters' !! termination message
    end type vector_root_result

    abstract interface
        subroutine vector_root_function(x, fvec)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent variable vector
            real(dp), intent(out) :: fvec(:) !! residual vector at x
        end subroutine vector_root_function

        subroutine vector_root_jacobian(x, fjac)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent variable vector
            real(dp), intent(out) :: fjac(:, :) !! residual-by-variable Jacobian
        end subroutine vector_root_jacobian
    end interface

    public :: anderson
    public :: broyden1
    public :: broyden2
    public :: diagbroyden
    public :: excitingmixing
    public :: fsolve
    public :: linearmixing
    public :: newton_krylov
    public :: root
    public :: vector_root_function
    public :: vector_root_jacobian
    public :: vector_root_norm_function
    public :: vector_root_linear_operator
    public :: vector_root_preconditioner_update
    public :: vector_root_inner_solver
    public :: vector_root_jacobian_product
    public :: vector_root_base_diagonal

contains

    function root(fun, x0, method, jac, m, tol, xtol, ftol, gtol, maxfev, &
            maxiter, band, eps, factor, diag, nit, fatol, xatol, line_search, &
            alpha, anderson_m, w0, alphamax, sane_m, sigma_eps, sigma_0, rdiff, &
            inner_method, inner_maxiter, outer_k, reduction_method, max_rank, &
            tol_norm, inner_m, inner_m_update, inner_solver, &
            inner_verify_residual, inner_tfqmr_right_preconditioned, &
            jacobian_diff_scheme, jacobian_step_mode, jacobian_product, &
            inner_check_residual, inner_orthogonalization, outer_prune_tol, &
            outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, &
            outer_image_basis_reuse, outer_candidate_cache_max_mib, inner_restart, &
            inner_restart_adaptive, inner_restart_max, inner_m_diag, &
            inner_m_diag_auto, inner_m_diag_update_interval, inner_m_diag_floor, &
            inner_m_diag_probes, inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed, inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol, inner_m_diag_probe_max, &
            inner_m_diag_rank_one, inner_m_diag_rank_one_right, &
            inner_m_diag_rank_two, inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch, inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two_discover, &
            inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base, jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive, inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth, inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost, inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost, inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        character(len=*), intent(in), optional :: method !! hybr or lm; default hybr
        procedure(vector_root_jacobian), optional :: jac !! analytic Jacobian callback
        integer, intent(in), optional :: m !! residual count for lm; default size(x0)
        real(dp), intent(in), optional :: tol !! generic SciPy root tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: ftol !! lm relative objective tolerance
        real(dp), intent(in), optional :: gtol !! lm orthogonality tolerance
        integer, intent(in), optional :: maxfev !! hybr function-evaluation cap
        integer, intent(in), optional :: maxiter !! lm function-evaluation cap
        integer, intent(in), optional :: band(:) !! hybr [subdiagonals, superdiagonals]
        real(dp), intent(in), optional :: eps !! finite-difference step control
        real(dp), intent(in), optional :: factor !! initial trust-region scale, (0.1,100)
        real(dp), intent(in), optional :: diag(:) !! positive variable scaling factors
        integer, intent(in), optional :: nit !! exact inexact-Newton iteration count
        real(dp), intent(in), optional :: fatol !! absolute inexact-Newton residual tolerance
        real(dp), intent(in), optional :: xatol !! absolute inexact-Newton step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial inexact-Newton Jacobian scaling
        integer, intent(in), optional :: anderson_m !! Anderson retained history, default 5
        real(dp), intent(in), optional :: w0 !! Anderson regularization, default 0.01
        real(dp), intent(in), optional :: alphamax !! exciting-mixing maximum scale
        integer, intent(in), optional :: sane_m !! DF-SANE nonmonotone window, default 10
        real(dp), intent(in), optional :: sigma_eps !! DF-SANE spectral safeguard
        real(dp), intent(in), optional :: sigma_0 !! DF-SANE initial spectral coefficient
        real(dp), intent(in), optional :: rdiff !! Krylov numerical-Jacobian step scale
        character(len=*), intent(in), optional :: inner_method !! lgmres, gmres, bicgstab, cgs, minres, or tfqmr
        integer, intent(in), optional :: inner_maxiter !! Krylov inner dimension, default 20
        integer, intent(in), optional :: outer_k !! retained LGMRES vectors, default 10
        character(len=*), intent(in), optional :: reduction_method !! Broyden rank reduction
        integer, intent(in), optional :: max_rank !! Broyden maximum inverse rank
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        procedure(vector_root_linear_operator), optional :: inner_m !! Krylov inverse-Jacobian action
        procedure(vector_root_preconditioner_update), optional :: inner_m_update !! accepted-step update
        procedure(vector_root_inner_solver), optional :: inner_solver !! custom Krylov linear solve
        logical, intent(in), optional :: inner_verify_residual !! verify original TFQMR residual
        logical, intent(in), optional :: inner_tfqmr_right_preconditioned !! use consistent right TFQMR
        character(len=*), intent(in), optional :: jacobian_diff_scheme !! Krylov 2-point or 3-point
        character(len=*), intent(in), optional :: jacobian_step_mode !! global (default), adaptive, or componentwise
        procedure(vector_root_jacobian_product), optional :: jacobian_product !! direct J(x)*vec callback
        logical, intent(in), optional :: inner_check_residual !! require true linear residual before accepting a step
        character(len=*), intent(in), optional :: inner_orthogonalization !! once (default), selective, or twice
        real(dp), intent(in), optional :: outer_prune_tol !! opt-in retained LGMRES dependence threshold
        real(dp), intent(in), optional :: outer_image_prune_tol !! current-Jacobian LGMRES dependence
        logical, intent(in), optional :: outer_image_prune_auto !! automatically infer image-space dependence cutoff
        logical, intent(in), optional :: outer_image_reuse !! reuse current-Jacobian pruning images in Arnoldi
        logical, intent(in), optional :: outer_candidate_image_reuse !! reconstruct candidate Jv from Arnoldi images
        logical, intent(in), optional :: outer_fifo_prune_fast !! use QR downdates for full solution-pruning history
        logical, intent(in), optional :: outer_image_basis_reuse !! reuse current-Jacobian image basis
        integer, intent(in), optional :: outer_candidate_cache_max_mib !! max raw Jv cache MiB; zero unlimited
        integer, intent(in), optional :: inner_restart !! GMRES/LGMRES restart subspace dimension
        logical, intent(in), optional :: inner_restart_adaptive !! grow LGMRES restart width on stagnation
        integer, intent(in), optional :: inner_restart_max !! maximum adaptive LGMRES restart width
        real(dp), intent(in), optional :: inner_m_diag(:) !! diagonal inverse-Jacobian preconditioner
        logical, intent(in), optional :: inner_m_diag_auto !! automatically estimate a diagonal inverse Jacobian
        integer, intent(in), optional :: inner_m_diag_update_interval !! accepted steps between refreshes; zero initial only
        real(dp), intent(in), optional :: inner_m_diag_floor !! relative floor for tiny estimated diagonal entries
        integer, intent(in), optional :: inner_m_diag_probes !! zero exact; positive Walsh probe count
        character(len=*), intent(in), optional :: inner_m_diag_probe_scheme !! walsh (default) or rademacher
        integer, intent(in), optional :: inner_m_diag_probe_seed !! reproducible Rademacher seed
        integer, intent(in), optional :: inner_m_diag_probe_max !! progressive Rademacher probe budget
        logical, intent(in), optional :: inner_m_diag_rank_one !! opt-in checked D + u*1^T inverse
        real(dp), intent(in), optional :: inner_m_diag_rank_one_right(:) !! fixed right coupling weights
        logical, intent(in), optional :: inner_m_diag_rank_one_discover !! learn unknown right weights
        logical, intent(in), optional :: inner_m_diag_rank_two !! checked rank-two Woodbury inverse
        real(dp), intent(in), optional :: inner_m_diag_rank_two_right(:,:) !! known n-by-2 right factors
        logical, intent(in), optional :: inner_m_diag_rank_two_sketch !! opt-in 3-probe diagonal estimate
        logical, intent(in), optional :: inner_m_diag_rank_two_discover !! learn unknown rank-two factors
        logical, intent(in), optional :: inner_m_diag_rank_two_stream !! O(n)-memory discovery
        real(dp), intent(in), optional :: inner_m_diag_rank_two_symmetric_base(:) !! known diagonal D in symmetric J=D+L
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_base(:) !! known D for nonsymmetric rank-two L
        procedure(vector_root_jacobian_product), optional :: jacobian_transpose_product !! analytic J(x)^T*v
        logical, intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_verify_transpose !! check two adjoint directions
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost !! J^T*v cost / Jv cost (> 0)
        procedure(vector_root_base_diagonal), optional :: inner_m_diag_rank_two_nonsymmetric_base_fun !! D(x) callback
        logical, intent(in), optional :: inner_m_diag_update_adaptive !! refresh early on weak residual reduction
        real(dp), intent(in), optional :: inner_m_diag_update_stall_ratio !! poor-step ratio threshold in (0,1]
        real(dp), intent(in), optional :: inner_m_diag_update_work_growth !! inner preconditioner work growth factor
        integer, intent(in), optional :: inner_m_diag_update_rebuild_cost !! rebuild cost in inner actions
        logical, intent(in), optional :: inner_m_diag_update_measured_cost !! count actual setup products
        real(dp), intent(in), optional :: inner_m_diag_update_transpose_cost !! J^T*v weight relative to Jv
        logical, intent(in), optional :: inner_m_diag_update_timed_cost !! time setup and inner solve
        real(dp), intent(in), optional :: inner_m_diag_update_timed_min_seconds !! minimum reliable timing interval
        logical, intent(in), optional :: inner_m_diag_update_wall_clock !! use wall seconds for timed refresh
        integer, intent(in), optional :: inner_m_diag_verify_samples !! diagonal spot checks; fall back if biased
        real(dp), intent(in), optional :: inner_m_diag_verify_rtol !! relative spot-check discrepancy limit
        type(vector_root_result) :: res

        character(len=16) :: meth
        integer :: m_use, n

        n = size(x0)
        meth = 'hybr'
        if (present(method)) meth = lowercase(trim(method))
        if (len_trim(meth) == 0) meth = 'hybr'
        m_use = n
        if (trim(meth) == 'lm' .and. present(m)) m_use = m

        call init_result(res, meth, n, m_use)
        if (n <= 0) return
        if (.not. all(ieee_is_finite(x0))) return

        if (present(inner_m) .or. present(inner_m_update) .or. present(inner_solver) .or. &
                present(inner_m_diag) .or. present(inner_m_diag_auto) .or. &
                present(inner_m_diag_update_interval) .or. present(inner_m_diag_floor) .or. &
                present(inner_m_diag_update_adaptive) .or. &
                present(inner_m_diag_update_stall_ratio) .or. &
                present(inner_m_diag_update_work_growth) .or. &
                present(inner_m_diag_update_rebuild_cost) .or. &
                present(inner_m_diag_update_measured_cost) .or. &
                present(inner_m_diag_update_transpose_cost) .or. &
                present(inner_m_diag_update_timed_cost) .or. &
                present(inner_m_diag_update_timed_min_seconds) .or. &
                present(inner_m_diag_update_wall_clock) .or. &
                present(inner_m_diag_probes) .or. present(inner_m_diag_probe_scheme) .or. &
                present(inner_m_diag_probe_seed) .or. present(inner_m_diag_verify_samples) .or. &
                present(inner_m_diag_verify_rtol) .or. &
                present(inner_m_diag_probe_max) .or. present(inner_m_diag_rank_one) .or. &
                present(inner_m_diag_rank_one_right) .or. &
                present(inner_m_diag_rank_one_discover) .or. &
                present(inner_m_diag_rank_two) .or. present(inner_m_diag_rank_two_right) .or. &
                present(inner_m_diag_rank_two_sketch) .or. &
                present(inner_m_diag_rank_two_discover) .or. &
                present(inner_m_diag_rank_two_stream) .or. &
                present(inner_m_diag_rank_two_symmetric_base) .or. &
                present(inner_m_diag_rank_two_nonsymmetric_base) .or. &
                present(inner_m_diag_rank_two_nonsymmetric_base_fun) .or. &
                present(jacobian_transpose_product) .or. &
                present(inner_m_diag_rank_two_nonsymmetric_verify_transpose) .or. &
                present(inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost)) then
            if (trim(meth) /= 'krylov') then
                res%message = 'inner Krylov hooks require method=krylov'
                return
            end if
        end if
        if (present(outer_image_basis_reuse)) then
            if (outer_image_basis_reuse .and. trim(meth) /= 'krylov') then
                res%message = 'outer_image_basis_reuse requires method=krylov'
                return
            end if
        end if
        if (present(outer_fifo_prune_fast)) then
            if (outer_fifo_prune_fast .and. trim(meth) /= 'krylov') then
                res%message = 'outer_fifo_prune_fast requires method=krylov'
                return
            end if
        end if
        if (present(outer_prune_tol) .and. trim(meth) /= 'krylov') then
            res%message = 'outer_prune_tol requires method=krylov'
            return
        end if
        if (present(outer_image_prune_tol) .and. trim(meth) /= 'krylov') then
            res%message = 'outer_image_prune_tol requires method=krylov'
            return
        end if
        if (present(outer_image_prune_auto)) then
            if (outer_image_prune_auto .and. trim(meth) /= 'krylov') then
                res%message = 'outer_image_prune_auto requires method=krylov'
                return
            end if
        end if
        if (present(outer_image_reuse)) then
            if (outer_image_reuse .and. trim(meth) /= 'krylov') then
                res%message = 'outer_image_reuse requires method=krylov'
                return
            end if
        end if
        if (present(outer_candidate_cache_max_mib) .and. trim(meth) /= 'krylov') then
            res%message = 'outer_candidate_cache_max_mib requires method=krylov'
            return
        end if
        if (present(outer_candidate_image_reuse)) then
            if (outer_candidate_image_reuse .and. trim(meth) /= 'krylov') then
                res%message = 'outer_candidate_image_reuse requires method=krylov'
                return
            end if
        end if
        if (present(inner_restart) .and. trim(meth) /= 'krylov') then
            res%message = 'inner_restart requires method=krylov'
            return
        end if
        if (present(inner_restart_adaptive) .and. trim(meth) /= 'krylov') then
            res%message = 'inner_restart_adaptive requires method=krylov'
            return
        end if
        if (present(inner_restart_max) .and. trim(meth) /= 'krylov') then
            res%message = 'inner_restart_max requires method=krylov'
            return
        end if
        if (present(inner_orthogonalization) .and. trim(meth) /= 'krylov') then
            res%message = 'inner_orthogonalization requires method=krylov'
            return
        end if
        if (present(inner_check_residual)) then
            if (inner_check_residual .and. trim(meth) /= 'krylov') then
                res%message = 'inner_check_residual requires method=krylov'
                return
            end if
        end if
        if (present(inner_verify_residual)) then
            if (inner_verify_residual .and. trim(meth) /= 'krylov') then
                res%message = 'inner_verify_residual requires method=krylov'
                return
            end if
        end if
        if (present(inner_tfqmr_right_preconditioned)) then
            if (inner_tfqmr_right_preconditioned .and. trim(meth) /= 'krylov') then
                res%message = 'inner_tfqmr_right_preconditioned requires method=krylov'
                return
            end if
        end if
        if (present(jacobian_product) .and. trim(meth) /= 'krylov') then
            res%message = 'jacobian_product requires method=krylov'
            return
        end if
        if (present(jacobian_step_mode) .and. trim(meth) /= 'krylov') then
            res%message = 'jacobian_step_mode requires method=krylov'
            return
        end if
        if (present(jacobian_diff_scheme) .and. trim(meth) /= 'krylov') then
            res%message = 'jacobian_diff_scheme requires method=krylov'
            return
        end if
        if (present(tol_norm)) then
            select case (trim(meth))
            case ('broyden1', 'broyden2', 'anderson', 'linearmixing', &
                    'diagbroyden', 'excitingmixing', 'krylov')
            case default
                res%message = 'tol_norm requires an inexact-Newton or Krylov method'
                return
            end select
        end if

        select case (trim(meth))
        case ('hybr')
            if (present(m)) then
                if (m /= n) return
            end if
            res = solve_hybr(fun, x0, jac, tol, xtol, maxfev, band, eps, factor, diag)
        case ('lm')
            res = solve_lm(fun, x0, jac, m, tol, xtol, ftol, gtol, maxiter, eps, &
                factor, diag)
        case ('broyden1', 'broyden2', 'anderson', 'linearmixing', 'diagbroyden', &
                'excitingmixing')
            res = solve_inexact(fun, x0, meth, tol, ftol, fatol, xtol, xatol, &
                maxiter, nit, line_search, alpha, anderson_m, w0, alphamax, &
                reduction_method, max_rank, tol_norm)
        case ('df-sane')
            res = solve_sane(fun, x0, tol, ftol, fatol, maxfev, sane_m, sigma_eps, &
                sigma_0, line_search)
        case ('krylov')
            res = solve_krylov(fun, x0, tol, ftol, fatol, xtol, xatol, maxiter, &
                nit, line_search, rdiff, inner_method, inner_maxiter, outer_k, &
                tol_norm, inner_m, inner_m_update, inner_solver, &
                inner_verify_residual, inner_tfqmr_right_preconditioned, &
                jacobian_diff_scheme, jacobian_step_mode, jacobian_product, &
                inner_check_residual, inner_orthogonalization, outer_prune_tol, &
                outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, &
            outer_image_basis_reuse, outer_candidate_cache_max_mib, inner_restart, &
                inner_restart_adaptive, inner_restart_max, inner_m_diag, &
            inner_m_diag_auto, inner_m_diag_update_interval, inner_m_diag_floor, &
            inner_m_diag_probes, inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed, inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol, inner_m_diag_probe_max, &
            inner_m_diag_rank_one, inner_m_diag_rank_one_right, &
            inner_m_diag_rank_two, inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch, inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two_discover, inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base, jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive, inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth, inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost, inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost, inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock)
        case default
            return
        end select
    end function root

    function fsolve(fun, x0, jac, xtol, maxfev, band, epsfcn, factor, diag) &
            result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        procedure(vector_root_jacobian), optional :: jac !! analytic Jacobian callback
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        integer, intent(in), optional :: maxfev !! function-evaluation cap
        integer, intent(in), optional :: band(:) !! [subdiagonals, superdiagonals]
        real(dp), intent(in), optional :: epsfcn !! finite-difference step control
        real(dp), intent(in), optional :: factor !! initial trust-region scale
        real(dp), intent(in), optional :: diag(:) !! positive variable scaling factors
        type(vector_root_result) :: res

        res = solve_hybr(fun, x0, jac, xtol=xtol, maxfev=maxfev, band=band, &
            eps=epsfcn, factor=factor, diag=diag)
        res%method = 'hybr'
    end function fsolve

    function broyden1(fun, x0, iter, maxiter, f_tol, f_rtol, x_tol, x_rtol, &
            line_search, alpha, reduction_method, max_rank, &
            tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of iterations to make
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial Jacobian scaling
        character(len=*), intent(in), optional :: reduction_method !! restart, simple, or svd
        integer, intent(in), optional :: max_rank !! maximum inverse-Jacobian rank
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        res = root(fun, x0, method='broyden1', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, alpha=alpha, reduction_method=reduction_method, &
            max_rank=max_rank, tol_norm=tol_norm)
    end function broyden1

    function broyden2(fun, x0, iter, maxiter, f_tol, f_rtol, x_tol, x_rtol, &
            line_search, alpha, reduction_method, max_rank, &
            tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of iterations to make
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial Jacobian scaling
        character(len=*), intent(in), optional :: reduction_method !! restart, simple, or svd
        integer, intent(in), optional :: max_rank !! maximum inverse-Jacobian rank
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        res = root(fun, x0, method='broyden2', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, alpha=alpha, reduction_method=reduction_method, &
            max_rank=max_rank, tol_norm=tol_norm)
    end function broyden2

    function anderson(fun, x0, iter, maxiter, f_tol, f_rtol, x_tol, x_rtol, &
            line_search, alpha, m, w0, tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of iterations to make
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial Jacobian scaling
        integer, intent(in), optional :: m !! retained Anderson history, default 5
        real(dp), intent(in), optional :: w0 !! Anderson regularization, default 0.01
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        res = root(fun, x0, method='anderson', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, alpha=alpha, anderson_m=m, w0=w0, &
            tol_norm=tol_norm)
    end function anderson

    function linearmixing(fun, x0, iter, maxiter, f_tol, f_rtol, x_tol, x_rtol, &
            line_search, alpha, tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of iterations to make
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! scalar mixing parameter
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        res = root(fun, x0, method='linearmixing', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, alpha=alpha, tol_norm=tol_norm)
    end function linearmixing

    function diagbroyden(fun, x0, iter, maxiter, f_tol, f_rtol, x_tol, x_rtol, &
            line_search, alpha, tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of iterations to make
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial diagonal scaling
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        res = root(fun, x0, method='diagbroyden', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, alpha=alpha, tol_norm=tol_norm)
    end function diagbroyden

    function excitingmixing(fun, x0, iter, maxiter, f_tol, f_rtol, x_tol, x_rtol, &
            line_search, alpha, alphamax, tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of iterations to make
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial diagonal mixing scale
        real(dp), intent(in), optional :: alphamax !! maximum diagonal mixing scale
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        res = root(fun, x0, method='excitingmixing', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, alpha=alpha, alphamax=alphamax, &
            tol_norm=tol_norm)
    end function excitingmixing

    function solve_hybr(fun, x0, jac, tol, xtol, maxfev, band, eps, factor, &
            diag) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        procedure(vector_root_jacobian), optional :: jac !! analytic Jacobian callback
        real(dp), intent(in), optional :: tol !! generic SciPy root tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        integer, intent(in), optional :: maxfev !! function-evaluation cap
        integer, intent(in), optional :: band(:) !! [subdiagonals, superdiagonals]
        real(dp), intent(in), optional :: eps !! finite-difference step control
        real(dp), intent(in), optional :: factor !! initial trust-region scale
        real(dp), intent(in), optional :: diag(:) !! positive variable scaling factors
        type(vector_root_result) :: res

        procedure(vector_root_function), pointer :: fun_ptr
        procedure(vector_root_jacobian), pointer :: jac_ptr
        real(dp), allocatable :: diag_use(:)
        real(dp), allocatable :: fvec(:)
        real(dp), allocatable :: fjac(:,:)
        real(dp), allocatable :: qtf(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: wa1(:), wa2(:), wa3(:), wa4(:)
        real(dp) :: eps_use, factor_use, xtol_use
        real(dp), allocatable :: preflight(:)
        integer :: info, lr, maxfev_use, ml, mode, mu, n, nfev, njev

        n = size(x0)
        call init_result(res, 'hybr', n, n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return

        xtol_use = default_xtol
        if (present(tol)) xtol_use = tol
        if (present(xtol)) xtol_use = xtol
        eps_use = epsilon(1.0_dp)
        if (present(eps)) eps_use = eps
        factor_use = default_factor
        if (present(factor)) factor_use = factor
        if (.not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(eps_use) .or. eps_use < 0.0_dp) return
        if (.not. ieee_is_finite(factor_use) .or. factor_use <= 0.1_dp .or. &
                factor_use > 100.0_dp) return

        maxfev_use = merge(100*(n + 1), 200*(n + 1), present(jac))
        if (present(maxfev)) then
            if (maxfev < 0) return
            if (maxfev > 0) maxfev_use = maxfev
        end if

        ml = n - 1
        mu = n - 1
        if (present(band) .and. .not. present(jac)) then
            if (size(band) /= 2) return
            if (any(band < 0)) return
            ml = band(1)
            mu = band(2)
        end if

        allocate(diag_use(n))
        diag_use = 1.0_dp
        mode = 1
        if (present(diag)) then
            if (size(diag) /= n) return
            if (.not. all(ieee_is_finite(diag)) .or. any(diag <= 0.0_dp)) return
            diag_use = diag
            mode = 2
        end if

        allocate(fvec(n), fjac(n,n), qtf(n), r(n*(n + 1)/2))
        allocate(wa1(n), wa2(n), wa3(n), wa4(n), preflight(n))
        res%x = x0
        fun_ptr => fun
        if (present(jac)) jac_ptr => jac

        ! SciPy's Python wrapper validates the callback before entering MINPACK.
        ! Two objective calls reproduce root/fsolve 1.17 work-count semantics.
        call fun_ptr(res%x, preflight)
        call fun_ptr(res%x, preflight)
        if (present(jac)) call jac_ptr(res%x, fjac)
        if (.not. all(ieee_is_finite(preflight))) then
            res%fun = preflight
            res%nfev = 2
            res%message = 'nonfinite residual encountered'
            return
        end if

        lr = size(r)
        nfev = 0
        njev = 0
        if (present(jac)) then
            call hybrj(hybrj_adapter, n, res%x, fvec, fjac, n, xtol_use, &
                maxfev_use, diag_use, mode, factor_use, 0, info, nfev, njev, &
                r, lr, qtf, wa1, wa2, wa3, wa4)
        else
            call hybrd(hybrd_adapter, n, res%x, fvec, xtol_use, maxfev_use, &
                ml, mu, eps_use, diag_use, mode, factor_use, 0, info, nfev, &
                fjac, n, r, lr, qtf, wa1, wa2, wa3, wa4)
        end if

        res%fun = fvec
        res%fjac = fjac
        res%r = r
        res%qtf = qtf
        res%nfev = nfev + 2
        res%njev = njev
        res%status = info
        res%success = info == 1
        res%message = hybr_message(info)

    contains

        subroutine hybrd_adapter(n_cb, x, f, iflag)
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current variable vector
            real(dp), intent(out) :: f(n_cb) !! residual vector
            integer, intent(inout) :: iflag !! MINPACK callback control flag

            if (iflag < 0) return
            call fun_ptr(x, f)
        end subroutine hybrd_adapter

        subroutine hybrj_adapter(n_cb, x, f, j, ldj, iflag)
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current variable vector
            real(dp), intent(inout) :: f(n_cb) !! residual vector
            integer, intent(in) :: ldj !! Jacobian leading dimension
            real(dp), intent(inout) :: j(ldj,n_cb) !! Jacobian workspace
            integer, intent(inout) :: iflag !! MINPACK callback control flag

            if (iflag == 1) then
                call fun_ptr(x, f)
            else if (iflag == 2) then
                call jac_ptr(x, j)
            end if
        end subroutine hybrj_adapter

    end function solve_hybr

    function solve_lm(fun, x0, jac, m, tol, xtol, ftol, gtol, maxiter, eps, &
            factor, diag) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        procedure(vector_root_jacobian), optional :: jac !! analytic Jacobian callback
        integer, intent(in), optional :: m !! residual count; default size(x0)
        real(dp), intent(in), optional :: tol !! generic SciPy root tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: ftol !! relative objective tolerance
        real(dp), intent(in), optional :: gtol !! orthogonality tolerance
        integer, intent(in), optional :: maxiter !! function-evaluation cap
        real(dp), intent(in), optional :: eps !! finite-difference step control
        real(dp), intent(in), optional :: factor !! initial trust-region scale
        real(dp), intent(in), optional :: diag(:) !! positive variable scaling factors
        type(vector_root_result) :: res

        procedure(vector_root_function), pointer :: fun_ptr
        procedure(vector_root_jacobian), pointer :: jac_ptr
        real(dp), allocatable :: diag_use(:), fvec(:), fjac_native(:,:)
        real(dp), allocatable :: preflight(:), qtf(:), wa1(:), wa2(:), wa3(:), wa4(:)
        real(dp) :: eps_use, factor_use, ftol_use, gtol_use, xtol_use
        integer, allocatable :: ipvt(:)
        integer :: info, maxiter_use, mode, m_use, n, nfev, njev

        n = size(x0)
        m_use = n
        if (present(m)) m_use = m
        call init_result(res, 'lm', n, m_use)
        if (n <= 0 .or. m_use < n) return
        if (.not. all(ieee_is_finite(x0))) return

        xtol_use = default_xtol
        if (present(tol)) xtol_use = tol
        if (present(xtol)) xtol_use = xtol
        ftol_use = default_ftol
        if (present(ftol)) ftol_use = ftol
        gtol_use = default_gtol
        if (present(gtol)) gtol_use = gtol
        eps_use = 0.0_dp
        if (present(eps)) eps_use = eps
        factor_use = default_factor
        if (present(factor)) factor_use = factor
        if (.not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(ftol_use) .or. ftol_use < 0.0_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(eps_use) .or. eps_use < 0.0_dp) return
        if (.not. ieee_is_finite(factor_use) .or. factor_use <= 0.1_dp .or. &
                factor_use > 100.0_dp) return

        maxiter_use = merge(100*(n + 1), 200*(n + 1), present(jac))
        if (present(maxiter)) then
            if (maxiter < 0) return
            if (maxiter > 0) maxiter_use = maxiter
        end if

        allocate(diag_use(n))
        diag_use = 1.0_dp
        mode = 1
        if (present(diag)) then
            if (size(diag) /= n) return
            if (.not. all(ieee_is_finite(diag)) .or. any(diag <= 0.0_dp)) return
            diag_use = diag
            mode = 2
        end if

        allocate(fvec(m_use), fjac_native(m_use,n), ipvt(n), qtf(n))
        allocate(wa1(n), wa2(n), wa3(n), wa4(m_use), preflight(m_use))
        res%x = x0
        fun_ptr => fun
        if (present(jac)) jac_ptr => jac

        call fun_ptr(res%x, preflight)
        call fun_ptr(res%x, preflight)
        if (present(jac)) call jac_ptr(res%x, fjac_native)
        if (.not. all(ieee_is_finite(preflight))) then
            res%fun = preflight
            res%nfev = 2
            res%message = 'nonfinite residual encountered'
            return
        end if

        nfev = 0
        njev = 0
        if (present(jac)) then
            call lmder(lmder_adapter, m_use, n, res%x, fvec, fjac_native, m_use, &
                ftol_use, xtol_use, gtol_use, maxiter_use, diag_use, mode, &
                factor_use, 0, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4)
        else
            call lmdif(lmdif_adapter, m_use, n, res%x, fvec, ftol_use, xtol_use, &
                gtol_use, maxiter_use, eps_use, diag_use, mode, factor_use, 0, &
                info, nfev, fjac_native, m_use, ipvt, qtf, wa1, wa2, wa3, wa4)
        end if

        res%fun = fvec
        res%fjac = transpose(fjac_native)
        res%ipvt = ipvt
        res%qtf = qtf
        res%nfev = nfev + 2
        res%njev = njev
        res%status = info
        res%success = info >= 1 .and. info <= 4
        res%message = lm_message(info)

    contains

        subroutine lmdif_adapter(m_cb, n_cb, x, f, iflag)
            integer, intent(in) :: m_cb !! MINPACK residual count
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current variable vector
            real(dp), intent(out) :: f(m_cb) !! residual vector
            integer, intent(inout) :: iflag !! MINPACK callback control flag

            if (iflag < 0) return
            call fun_ptr(x, f)
        end subroutine lmdif_adapter

        subroutine lmder_adapter(m_cb, n_cb, x, f, j, ldj, iflag)
            integer, intent(in) :: m_cb !! MINPACK residual count
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current variable vector
            real(dp), intent(inout) :: f(m_cb) !! residual vector
            integer, intent(in) :: ldj !! Jacobian leading dimension
            real(dp), intent(inout) :: j(ldj,n_cb) !! residual-by-variable Jacobian
            integer, intent(inout) :: iflag !! MINPACK callback control flag

            if (iflag == 1) then
                call fun_ptr(x, f)
            else if (iflag == 2) then
                call jac_ptr(x, j)
            end if
        end subroutine lmder_adapter

    end function solve_lm

    function newton_krylov(fun, x0, iter, rdiff, inner_method, inner_maxiter, &
            outer_k, maxiter, f_tol, f_rtol, x_tol, x_rtol, line_search, &
            tol_norm, inner_m, inner_m_update, inner_solver, &
            inner_verify_residual, inner_tfqmr_right_preconditioned, &
            jacobian_diff_scheme, jacobian_step_mode, jacobian_product, &
            inner_check_residual, inner_orthogonalization, outer_prune_tol, &
            outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, &
            outer_image_basis_reuse, outer_candidate_cache_max_mib, inner_restart, &
            inner_restart_adaptive, inner_restart_max, inner_m_diag, &
            inner_m_diag_auto, inner_m_diag_update_interval, inner_m_diag_floor, &
            inner_m_diag_probes, inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed, inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol, inner_m_diag_probe_max, &
            inner_m_diag_rank_one, inner_m_diag_rank_one_right, &
            inner_m_diag_rank_two, inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch, inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two_discover, &
            inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base, jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive, inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth, inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost, inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost, inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        integer, intent(in), optional :: iter !! exact number of nonlinear iterations
        real(dp), intent(in), optional :: rdiff !! numerical-Jacobian relative step
        character(len=*), intent(in), optional :: inner_method !! lgmres, gmres, bicgstab, cgs, minres, or tfqmr
        integer, intent(in), optional :: inner_maxiter !! inner Krylov dimension
        integer, intent(in), optional :: outer_k !! retained LGMRES augmentation vectors
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        real(dp), intent(in), optional :: f_tol !! absolute residual max-norm tolerance
        real(dp), intent(in), optional :: f_rtol !! relative residual tolerance
        real(dp), intent(in), optional :: x_tol !! absolute step max-norm tolerance
        real(dp), intent(in), optional :: x_rtol !! relative step tolerance
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        procedure(vector_root_linear_operator), optional :: inner_m !! Krylov preconditioner
        procedure(vector_root_preconditioner_update), optional :: inner_m_update !! accepted-step update
        procedure(vector_root_inner_solver), optional :: inner_solver !! custom inner linear solver
        logical, intent(in), optional :: inner_verify_residual !! verify original TFQMR residual
        logical, intent(in), optional :: inner_tfqmr_right_preconditioned !! use consistent right TFQMR
        character(len=*), intent(in), optional :: jacobian_diff_scheme !! Krylov 2-point or 3-point
        character(len=*), intent(in), optional :: jacobian_step_mode !! global (default), adaptive, or componentwise
        procedure(vector_root_jacobian_product), optional :: jacobian_product !! direct J(x)*vec callback
        logical, intent(in), optional :: inner_check_residual !! require true linear residual before accepting a step
        character(len=*), intent(in), optional :: inner_orthogonalization !! once (default), selective, or twice
        real(dp), intent(in), optional :: outer_prune_tol !! opt-in retained LGMRES dependence threshold
        real(dp), intent(in), optional :: outer_image_prune_tol !! current-Jacobian LGMRES dependence
        logical, intent(in), optional :: outer_image_prune_auto !! automatically infer image-space dependence cutoff
        logical, intent(in), optional :: outer_image_reuse !! reuse current-Jacobian pruning images in Arnoldi
        logical, intent(in), optional :: outer_candidate_image_reuse !! reconstruct candidate Jv from Arnoldi images
        logical, intent(in), optional :: outer_fifo_prune_fast !! use QR downdates for full solution-pruning history
        logical, intent(in), optional :: outer_image_basis_reuse !! reuse current-Jacobian image basis
        integer, intent(in), optional :: outer_candidate_cache_max_mib !! max raw Jv cache MiB; zero unlimited
        integer, intent(in), optional :: inner_restart !! GMRES/LGMRES restart subspace dimension
        logical, intent(in), optional :: inner_restart_adaptive !! grow LGMRES restart width on stagnation
        integer, intent(in), optional :: inner_restart_max !! maximum adaptive LGMRES restart width
        real(dp), intent(in), optional :: inner_m_diag(:) !! diagonal inverse-Jacobian preconditioner
        logical, intent(in), optional :: inner_m_diag_auto !! automatically estimate a diagonal inverse Jacobian
        integer, intent(in), optional :: inner_m_diag_update_interval !! accepted steps between refreshes; zero initial only
        real(dp), intent(in), optional :: inner_m_diag_floor !! relative floor for tiny estimated diagonal entries
        integer, intent(in), optional :: inner_m_diag_probes !! zero exact; positive Walsh probe count
        character(len=*), intent(in), optional :: inner_m_diag_probe_scheme !! walsh (default) or rademacher
        integer, intent(in), optional :: inner_m_diag_probe_seed !! reproducible Rademacher seed
        integer, intent(in), optional :: inner_m_diag_probe_max !! progressive Rademacher probe budget
        logical, intent(in), optional :: inner_m_diag_rank_one !! opt-in checked D + u*1^T inverse
        real(dp), intent(in), optional :: inner_m_diag_rank_one_right(:) !! fixed right coupling weights
        logical, intent(in), optional :: inner_m_diag_rank_one_discover !! learn unknown right weights
        logical, intent(in), optional :: inner_m_diag_rank_two !! checked rank-two Woodbury inverse
        real(dp), intent(in), optional :: inner_m_diag_rank_two_right(:,:) !! known n-by-2 right factors
        logical, intent(in), optional :: inner_m_diag_rank_two_sketch !! opt-in 3-probe diagonal estimate
        logical, intent(in), optional :: inner_m_diag_rank_two_discover !! learn unknown rank-two factors
        logical, intent(in), optional :: inner_m_diag_rank_two_stream !! O(n)-memory discovery
        real(dp), intent(in), optional :: inner_m_diag_rank_two_symmetric_base(:) !! known diagonal D in symmetric J=D+L
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_base(:) !! known D for nonsymmetric rank-two L
        procedure(vector_root_jacobian_product), optional :: jacobian_transpose_product !! analytic J(x)^T*v
        logical, intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_verify_transpose !! check two adjoint directions
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost !! J^T*v cost / Jv cost (> 0)
        procedure(vector_root_base_diagonal), optional :: inner_m_diag_rank_two_nonsymmetric_base_fun !! D(x) callback
        logical, intent(in), optional :: inner_m_diag_update_adaptive !! refresh early on weak residual reduction
        real(dp), intent(in), optional :: inner_m_diag_update_stall_ratio !! poor-step ratio threshold in (0,1]
        real(dp), intent(in), optional :: inner_m_diag_update_work_growth !! inner preconditioner work growth factor
        integer, intent(in), optional :: inner_m_diag_update_rebuild_cost !! rebuild cost in inner actions
        logical, intent(in), optional :: inner_m_diag_update_measured_cost !! count actual setup products
        real(dp), intent(in), optional :: inner_m_diag_update_transpose_cost !! J^T*v weight relative to Jv
        logical, intent(in), optional :: inner_m_diag_update_timed_cost !! time setup and inner solve
        real(dp), intent(in), optional :: inner_m_diag_update_timed_min_seconds !! minimum reliable timing interval
        logical, intent(in), optional :: inner_m_diag_update_wall_clock !! use wall seconds for timed refresh
        integer, intent(in), optional :: inner_m_diag_verify_samples !! diagonal spot checks; fall back if biased
        real(dp), intent(in), optional :: inner_m_diag_verify_rtol !! relative spot-check discrepancy limit
        type(vector_root_result) :: res

        res = root(fun, x0, method='krylov', maxiter=maxiter, nit=iter, &
            fatol=f_tol, ftol=f_rtol, xatol=x_tol, xtol=x_rtol, &
            line_search=line_search, rdiff=rdiff, inner_method=inner_method, &
            inner_maxiter=inner_maxiter, outer_k=outer_k, tol_norm=tol_norm, &
            inner_m=inner_m, inner_m_update=inner_m_update, &
            inner_solver=inner_solver, inner_verify_residual=inner_verify_residual, &
            inner_tfqmr_right_preconditioned=inner_tfqmr_right_preconditioned, &
            jacobian_diff_scheme=jacobian_diff_scheme, &
            jacobian_step_mode=jacobian_step_mode, &
            jacobian_product=jacobian_product, &
            inner_check_residual=inner_check_residual, &
            inner_orthogonalization=inner_orthogonalization, &
            outer_prune_tol=outer_prune_tol, &
            outer_image_prune_tol=outer_image_prune_tol, &
            outer_image_prune_auto=outer_image_prune_auto, &
            outer_image_reuse=outer_image_reuse, &
            outer_candidate_image_reuse=outer_candidate_image_reuse, &
            outer_fifo_prune_fast=outer_fifo_prune_fast, &
            outer_image_basis_reuse=outer_image_basis_reuse, &
            outer_candidate_cache_max_mib=outer_candidate_cache_max_mib, &
            inner_restart=inner_restart, inner_restart_adaptive=inner_restart_adaptive, &
            inner_restart_max=inner_restart_max, inner_m_diag=inner_m_diag, &
            inner_m_diag_auto=inner_m_diag_auto, &
            inner_m_diag_update_interval=inner_m_diag_update_interval, &
            inner_m_diag_floor=inner_m_diag_floor, &
            inner_m_diag_probes=inner_m_diag_probes, &
            inner_m_diag_probe_scheme=inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed=inner_m_diag_probe_seed, &
            inner_m_diag_verify_samples=inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol=inner_m_diag_verify_rtol, &
            inner_m_diag_probe_max=inner_m_diag_probe_max, &
            inner_m_diag_rank_one=inner_m_diag_rank_one, &
            inner_m_diag_rank_one_right=inner_m_diag_rank_one_right, &
            inner_m_diag_rank_one_discover=inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two=inner_m_diag_rank_two, &
            inner_m_diag_rank_two_right=inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch=inner_m_diag_rank_two_sketch, &
            inner_m_diag_rank_two_discover=inner_m_diag_rank_two_discover, &
            inner_m_diag_rank_two_stream=inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base=inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base=inner_m_diag_rank_two_nonsymmetric_base, &
            jacobian_transpose_product=jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose=inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun=inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive=inner_m_diag_update_adaptive, &
            inner_m_diag_update_stall_ratio=inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth=inner_m_diag_update_work_growth, &
            inner_m_diag_update_rebuild_cost=inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost=inner_m_diag_update_measured_cost, &
            inner_m_diag_update_transpose_cost=inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost=inner_m_diag_update_timed_cost, &
            inner_m_diag_update_timed_min_seconds=inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock=inner_m_diag_update_wall_clock)
    end function newton_krylov

    function solve_krylov(fun, x0, tol, ftol, fatol, xtol, xatol, maxiter, &
            iterations, line_search, rdiff, inner_method, inner_maxiter, outer_k, &
            tol_norm, inner_m, inner_m_update, inner_solver, &
            inner_verify_residual, inner_tfqmr_right_preconditioned, &
            jacobian_diff_scheme, jacobian_step_mode, jacobian_product, &
            inner_check_residual, inner_orthogonalization, outer_prune_tol, &
            outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, &
            outer_image_basis_reuse, outer_candidate_cache_max_mib, inner_restart, &
            inner_restart_adaptive, inner_restart_max, inner_m_diag, &
            inner_m_diag_auto, inner_m_diag_update_interval, inner_m_diag_floor, &
            inner_m_diag_probes, inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed, inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol, inner_m_diag_probe_max, &
            inner_m_diag_rank_one, inner_m_diag_rank_one_right, &
            inner_m_diag_rank_two, inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch, inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two_discover, &
            inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base, jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive, inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth, inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost, inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost, inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        real(dp), intent(in), optional :: tol !! generic root tolerance
        real(dp), intent(in), optional :: ftol !! relative residual tolerance
        real(dp), intent(in), optional :: fatol !! absolute residual tolerance
        real(dp), intent(in), optional :: xtol !! relative step tolerance
        real(dp), intent(in), optional :: xatol !! absolute step tolerance
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        integer, intent(in), optional :: iterations !! exact nonlinear iterations
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: rdiff !! numerical-Jacobian step scale
        character(len=*), intent(in), optional :: inner_method !! lgmres, gmres, bicgstab, cgs, minres, or tfqmr
        integer, intent(in), optional :: inner_maxiter !! inner Krylov dimension
        integer, intent(in), optional :: outer_k !! retained LGMRES vectors
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        procedure(vector_root_linear_operator), optional :: inner_m !! Krylov preconditioner
        procedure(vector_root_preconditioner_update), optional :: inner_m_update !! accepted-step update
        procedure(vector_root_inner_solver), optional :: inner_solver !! custom inner linear solver
        logical, intent(in), optional :: inner_verify_residual !! verify original TFQMR residual
        logical, intent(in), optional :: inner_tfqmr_right_preconditioned !! use consistent right TFQMR
        character(len=*), intent(in), optional :: jacobian_diff_scheme !! Krylov 2-point or 3-point
        character(len=*), intent(in), optional :: jacobian_step_mode !! global (default), adaptive, or componentwise
        procedure(vector_root_jacobian_product), optional :: jacobian_product !! direct J(x)*vec callback
        logical, intent(in), optional :: inner_check_residual !! require true linear residual before accepting a step
        character(len=*), intent(in), optional :: inner_orthogonalization !! once (default), selective, or twice
        real(dp), intent(in), optional :: outer_prune_tol !! opt-in retained LGMRES dependence threshold
        real(dp), intent(in), optional :: outer_image_prune_tol !! current-Jacobian LGMRES dependence
        logical, intent(in), optional :: outer_image_prune_auto !! automatically infer image-space dependence cutoff
        logical, intent(in), optional :: outer_image_reuse !! reuse current-Jacobian pruning images in Arnoldi
        logical, intent(in), optional :: outer_candidate_image_reuse !! reconstruct candidate Jv from Arnoldi images
        logical, intent(in), optional :: outer_fifo_prune_fast !! use QR downdates for full solution-pruning history
        logical, intent(in), optional :: outer_image_basis_reuse !! reuse current-Jacobian image basis
        integer, intent(in), optional :: outer_candidate_cache_max_mib !! max raw Jv cache MiB; zero unlimited
        integer, intent(in), optional :: inner_restart !! GMRES/LGMRES restart subspace dimension
        logical, intent(in), optional :: inner_restart_adaptive !! grow LGMRES restart width on stagnation
        integer, intent(in), optional :: inner_restart_max !! maximum adaptive LGMRES restart width
        real(dp), intent(in), optional :: inner_m_diag(:) !! diagonal inverse-Jacobian preconditioner
        logical, intent(in), optional :: inner_m_diag_auto !! automatically estimate a diagonal inverse Jacobian
        integer, intent(in), optional :: inner_m_diag_update_interval !! accepted steps between refreshes; zero initial only
        real(dp), intent(in), optional :: inner_m_diag_floor !! relative floor for tiny estimated diagonal entries
        integer, intent(in), optional :: inner_m_diag_probes !! zero exact; positive Walsh probe count
        character(len=*), intent(in), optional :: inner_m_diag_probe_scheme !! walsh (default) or rademacher
        integer, intent(in), optional :: inner_m_diag_probe_seed !! reproducible Rademacher seed
        integer, intent(in), optional :: inner_m_diag_probe_max !! progressive Rademacher probe budget
        logical, intent(in), optional :: inner_m_diag_rank_one !! opt-in checked D + u*1^T inverse
        real(dp), intent(in), optional :: inner_m_diag_rank_one_right(:) !! fixed right coupling weights
        logical, intent(in), optional :: inner_m_diag_rank_one_discover !! learn unknown right weights
        logical, intent(in), optional :: inner_m_diag_rank_two !! checked rank-two Woodbury inverse
        real(dp), intent(in), optional :: inner_m_diag_rank_two_right(:,:) !! known n-by-2 right factors
        logical, intent(in), optional :: inner_m_diag_rank_two_sketch !! opt-in 3-probe diagonal estimate
        logical, intent(in), optional :: inner_m_diag_rank_two_discover !! learn unknown rank-two factors
        logical, intent(in), optional :: inner_m_diag_rank_two_stream !! O(n)-memory discovery
        real(dp), intent(in), optional :: inner_m_diag_rank_two_symmetric_base(:) !! known diagonal D in symmetric J=D+L
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_base(:) !! known D for nonsymmetric rank-two L
        procedure(vector_root_jacobian_product), optional :: jacobian_transpose_product !! analytic J(x)^T*v
        logical, intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_verify_transpose !! check two adjoint directions
        real(dp), intent(in), optional :: inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost !! J^T*v cost / Jv cost (> 0)
        procedure(vector_root_base_diagonal), optional :: inner_m_diag_rank_two_nonsymmetric_base_fun !! D(x) callback
        logical, intent(in), optional :: inner_m_diag_update_adaptive !! refresh early on weak residual reduction
        real(dp), intent(in), optional :: inner_m_diag_update_stall_ratio !! poor-step ratio threshold in (0,1]
        real(dp), intent(in), optional :: inner_m_diag_update_work_growth !! inner preconditioner work growth factor
        integer, intent(in), optional :: inner_m_diag_update_rebuild_cost !! rebuild cost in inner actions
        logical, intent(in), optional :: inner_m_diag_update_measured_cost !! count actual setup products
        real(dp), intent(in), optional :: inner_m_diag_update_transpose_cost !! J^T*v weight relative to Jv
        logical, intent(in), optional :: inner_m_diag_update_timed_cost !! time setup and inner solve
        real(dp), intent(in), optional :: inner_m_diag_update_timed_min_seconds !! minimum reliable timing interval
        logical, intent(in), optional :: inner_m_diag_update_wall_clock !! use wall seconds for timed refresh
        integer, intent(in), optional :: inner_m_diag_verify_samples !! diagonal spot checks; fall back if biased
        real(dp), intent(in), optional :: inner_m_diag_verify_rtol !! relative spot-check discrepancy limit
        type(vector_root_result) :: res

        integer :: n

        n = size(x0)
        call init_result(res, 'krylov', n, n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        res%x = x0
        call solve_krylov_root(fun, x0, res%x, res%fun, res%nfev, res%nit, &
            res%status, res%success, res%message, tol, ftol, fatol, xtol, xatol, &
            maxiter, iterations, line_search, rdiff, inner_method, inner_maxiter, &
            outer_k, tol_norm, inner_m, inner_m_update, inner_solver, &
            inner_verify_residual, inner_tfqmr_right_preconditioned, &
            jacobian_diff_scheme, jacobian_step_mode, jacobian_product, &
            inner_check_residual, inner_orthogonalization, outer_prune_tol, &
                outer_image_prune_tol, outer_image_prune_auto, outer_image_reuse, &
            outer_candidate_image_reuse, outer_fifo_prune_fast, &
            outer_image_basis_reuse, outer_candidate_cache_max_mib, inner_restart, &
                inner_restart_adaptive, inner_restart_max, inner_m_diag, &
            inner_m_diag_auto, inner_m_diag_update_interval, inner_m_diag_floor, &
            inner_m_diag_probes, inner_m_diag_probe_scheme, &
            inner_m_diag_probe_seed, inner_m_diag_verify_samples, &
            inner_m_diag_verify_rtol, inner_m_diag_probe_max, &
            inner_m_diag_rank_one, inner_m_diag_rank_one_right, &
            inner_m_diag_rank_two, inner_m_diag_rank_two_right, &
            inner_m_diag_rank_two_sketch, inner_m_diag_rank_one_discover, &
            inner_m_diag_rank_two_discover, inner_m_diag_rank_two_stream, &
            inner_m_diag_rank_two_symmetric_base, &
            inner_m_diag_rank_two_nonsymmetric_base, jacobian_transpose_product, &
            inner_m_diag_rank_two_nonsymmetric_verify_transpose, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost, &
            inner_m_diag_rank_two_nonsymmetric_base_fun, &
            inner_m_diag_update_adaptive, inner_m_diag_update_stall_ratio, &
            inner_m_diag_update_work_growth, inner_m_diag_update_rebuild_cost, &
            inner_m_diag_update_measured_cost, inner_m_diag_update_transpose_cost, &
            inner_m_diag_update_timed_cost, inner_m_diag_update_timed_min_seconds, &
            inner_m_diag_update_wall_clock)
    end function solve_krylov

    function solve_sane(fun, x0, tol, ftol, fatol, maxfev, history, sigma_eps, &
            sigma_0, line_search) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        real(dp), intent(in), optional :: tol !! generic root tolerance, mapped to ftol
        real(dp), intent(in), optional :: ftol !! relative residual norm tolerance
        real(dp), intent(in), optional :: fatol !! absolute residual norm tolerance
        integer, intent(in), optional :: maxfev !! residual-evaluation cap
        integer, intent(in), optional :: history !! nonmonotone history length
        real(dp), intent(in), optional :: sigma_eps !! spectral safeguard
        real(dp), intent(in), optional :: sigma_0 !! initial spectral coefficient
        character(len=*), intent(in), optional :: line_search !! cruz or cheng
        type(vector_root_result) :: res

        integer :: n

        n = size(x0)
        call init_result(res, 'df-sane', n, n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        res%x = x0
        call solve_df_sane(fun, x0, res%x, res%fun, res%nfev, res%nit, res%status, &
            res%success, res%message, tol, ftol, fatol, maxfev, history, sigma_eps, &
            sigma_0, line_search)
    end function solve_sane

    function solve_inexact(fun, x0, method, tol, ftol, fatol, xtol, xatol, &
            maxiter, iterations, line_search, alpha, anderson_m, w0, alphamax, &
            reduction_method, max_rank, tol_norm) result(res)
        procedure(vector_root_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        character(len=*), intent(in) :: method !! inexact-Newton method name
        real(dp), intent(in), optional :: tol !! generic root tolerance
        real(dp), intent(in), optional :: ftol !! relative residual tolerance
        real(dp), intent(in), optional :: fatol !! absolute residual tolerance
        real(dp), intent(in), optional :: xtol !! relative step tolerance
        real(dp), intent(in), optional :: xatol !! absolute step tolerance
        integer, intent(in), optional :: maxiter !! maximum nonlinear iterations
        integer, intent(in), optional :: iterations !! exact iterations to make
        character(len=*), intent(in), optional :: line_search !! armijo, wolfe, or none
        real(dp), intent(in), optional :: alpha !! initial Jacobian scaling
        integer, intent(in), optional :: anderson_m !! Anderson history length
        real(dp), intent(in), optional :: w0 !! Anderson regularization
        real(dp), intent(in), optional :: alphamax !! exciting-mixing scale cap
        character(len=*), intent(in), optional :: reduction_method !! Broyden rank reduction
        integer, intent(in), optional :: max_rank !! Broyden maximum inverse rank
        procedure(vector_root_norm_function), optional :: tol_norm !! custom stopping norm
        type(vector_root_result) :: res

        integer :: n

        n = size(x0)
        call init_result(res, method, n, n)
        if (n <= 0 .or. .not. all(ieee_is_finite(x0))) return
        res%x = x0
        call solve_nonlin_root(fun, x0, method, res%x, res%fun, res%nfev, &
            res%nit, res%status, res%success, res%message, tol, ftol, fatol, &
            xtol, xatol, maxiter, iterations, line_search, alpha, anderson_m, &
            w0, alphamax, reduction_method, max_rank, tol_norm)
    end function solve_inexact

    subroutine init_result(res, method, n, m)
        type(vector_root_result), intent(out) :: res !! result object to initialize
        character(len=*), intent(in) :: method !! method name
        integer, intent(in) :: n !! variable count
        integer, intent(in) :: m !! residual count

        allocate(res%x(max(0,n)), res%fun(max(0,m)))
        res%x = 0.0_dp
        res%fun = 0.0_dp
        res%nfev = 0
        res%njev = 0
        res%nit = 0
        res%status = 0
        res%success = .false.
        res%method = method
        res%message = 'improper input parameters'
    end subroutine init_result

    function hybr_message(status) result(message)
        integer, intent(in) :: status !! MINPACK termination status
        character(len=160) :: message

        select case (status)
        case (1)
            message = 'the solution converged'
        case (2)
            message = 'maximum number of function evaluations reached'
        case (3)
            message = 'xtol is too small; no further improvement is possible'
        case (4)
            message = 'iteration is not making good progress over five Jacobian evaluations'
        case (5)
            message = 'iteration is not making good progress over ten iterations'
        case default
            message = 'improper input parameters'
        end select
    end function hybr_message

    function lm_message(status) result(message)
        integer, intent(in) :: status !! MINPACK termination status
        character(len=160) :: message

        select case (status)
        case (1)
            message = 'actual and predicted relative reductions are at most ftol'
        case (2)
            message = 'relative error between consecutive iterates is at most xtol'
        case (3)
            message = 'both ftol and xtol termination conditions are satisfied'
        case (4)
            message = 'residual is orthogonal to the Jacobian columns to machine precision'
        case (5)
            message = 'maximum number of function evaluations reached'
        case (6)
            message = 'ftol is too small; no further reduction is possible'
        case (7)
            message = 'xtol is too small; no further improvement is possible'
        case (8)
            message = 'gtol is too small; residual is orthogonal to Jacobian columns'
        case default
            message = 'improper input parameters'
        end select
    end function lm_message

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text !! source text
        character(len=len(text)) :: lower
        integer :: i, code

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module scifort_optimize_root
