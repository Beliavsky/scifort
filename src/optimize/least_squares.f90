! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Nonlinear least squares and curve fitting for the scipy.optimize translation.
!
! The Levenberg-Marquardt numerical kernels are provided by the community-
! maintained fortran-lang/minpack FPM dependency. This module supplies SciPy
! 1.17-style adapters, finite differences, result bookkeeping, covariance
! reconstruction, and a Fortran-native curve_fit interface.

module scifort_optimize_least_squares
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_value
    use, intrinsic :: iso_fortran_env, only : int64
    use minpack_module, only : lmder, lmdif
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, symmetric_eigen_jacobi
    use scifort_optimize_lsmr, only : lsmr_solve
    implicit none
    private

    real(dp), parameter :: leastsq_default_tol = 1.49012e-8_dp
    real(dp), parameter :: least_squares_default_tol = 1.0e-8_dp
    real(dp), parameter :: default_factor = 100.0_dp

    type, public :: least_squares_result
        real(dp), allocatable :: x(:) !! final parameter vector
        real(dp) :: cost = 0.0_dp !! one half of the residual sum of squares
        real(dp), allocatable :: fun(:) !! residual vector at x
        real(dp), allocatable :: jac(:, :) !! dense residual-by-parameter Jacobian, when available
        logical :: jac_is_operator = .false. !! true when the final Jacobian is matrix-free
        real(dp), allocatable :: grad(:) !! J^T f at x
        real(dp) :: optimality = 0.0_dp !! infinity norm of the gradient
        integer, allocatable :: active_mask(:) !! bound activity; zero for LM
        integer :: nfev = 0 !! residual evaluations reported by MINPACK
        integer :: njev = -1 !! Jacobian evaluations, or -1 for numerical Jacobians
        integer :: status = -1 !! SciPy least_squares-style termination status
        logical :: success = .false. !! true for a successful termination
        character(len=16) :: method = 'trf' !! least-squares method
        character(len=160) :: message = 'improper input parameters' !! termination message
    end type least_squares_result

    type, public :: leastsq_result
        real(dp), allocatable :: x(:) !! final parameter vector
        real(dp), allocatable :: cov_x(:, :) !! unscaled covariance approximation
        real(dp), allocatable :: fvec(:) !! residual vector at x
        real(dp), allocatable :: fjac(:, :) !! MINPACK QR workspace in SciPy orientation
        integer, allocatable :: ipvt(:) !! one-based MINPACK pivot indices
        real(dp), allocatable :: qtf(:) !! Q^T f
        integer :: nfev = 0 !! residual evaluation count
        integer :: njev = -1 !! Jacobian evaluations, or -1 when numerical
        integer :: status = 0 !! raw MINPACK ier value
        logical :: success = .false. !! true for ier 1 through 4
        character(len=160) :: message = 'Improper input parameters.' !! termination message
    end type leastsq_result

    type, public :: curve_fit_result
        real(dp), allocatable :: popt(:) !! fitted parameters
        real(dp), allocatable :: pcov(:, :) !! parameter covariance matrix
        real(dp), allocatable :: fvec(:) !! weighted residual vector at the solution
        real(dp) :: cost = 0.0_dp !! one half of the weighted residual sum of squares
        integer :: nfev = 0 !! residual evaluation count
        integer :: njev = -1 !! Jacobian evaluation count, or -1 when numerical
        integer :: status = 0 !! raw MINPACK ier value
        logical :: success = .false. !! true for ier 1 through 4
        character(len=160) :: message = 'Improper input parameters.' !! termination message
    end type curve_fit_result

    abstract interface
        subroutine least_squares_residual_function(x, fvec)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(out) :: fvec(:) !! residual vector
        end subroutine least_squares_residual_function

        subroutine least_squares_complex_residual_function(x, fvec)
            import dp
            implicit none
            complex(dp), intent(in) :: x(:) !! complex-step parameter vector
            complex(dp), intent(out) :: fvec(:) !! analytic residual continuation
        end subroutine least_squares_complex_residual_function

        subroutine least_squares_jacobian_function(x, fjac)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(out) :: fjac(:, :) !! residual-by-parameter Jacobian
        end subroutine least_squares_jacobian_function

        subroutine least_squares_jacobian_matvec_function(x, v, jv)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(in) :: v(:) !! parameter-space vector
            real(dp), intent(out) :: jv(:) !! Jacobian-vector product
        end subroutine least_squares_jacobian_matvec_function

        subroutine least_squares_jacobian_rmatvec_function(x, w, jtw)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! parameter vector
            real(dp), intent(in) :: w(:) !! residual-space vector
            real(dp), intent(out) :: jtw(:) !! transposed-Jacobian-vector product
        end subroutine least_squares_jacobian_rmatvec_function

        subroutine least_squares_linear_operator_function(v, y)
            import dp
            implicit none
            real(dp), intent(in) :: v(:) !! operator input vector
            real(dp), intent(out) :: y(:) !! operator output vector
        end subroutine least_squares_linear_operator_function

        subroutine curve_model_function(xdata, params, ymodel)
            import dp
            implicit none
            real(dp), intent(in) :: xdata(:) !! independent data
            real(dp), intent(in) :: params(:) !! model parameters
            real(dp), intent(out) :: ymodel(:) !! model values
        end subroutine curve_model_function

        subroutine curve_model_complex_function(xdata, params, ymodel)
            import dp
            implicit none
            real(dp), intent(in) :: xdata(:) !! independent data
            complex(dp), intent(in) :: params(:) !! complex-step model parameters
            complex(dp), intent(out) :: ymodel(:) !! analytic model continuation
        end subroutine curve_model_complex_function

        subroutine curve_model_jacobian(xdata, params, fjac)
            import dp
            implicit none
            real(dp), intent(in) :: xdata(:) !! independent data
            real(dp), intent(in) :: params(:) !! model parameters
            real(dp), intent(out) :: fjac(:, :) !! model-value Jacobian
        end subroutine curve_model_jacobian

        subroutine robust_loss_function(z, rho0, rho1, rho2)
            import dp
            implicit none
            real(dp), intent(in) :: z(:) !! squared residuals scaled by f_scale
            real(dp), intent(out) :: rho0(:) !! loss values rho(z)
            real(dp), intent(out) :: rho1(:) !! first derivatives with respect to z
            real(dp), intent(out) :: rho2(:) !! second derivatives with respect to z
        end subroutine robust_loss_function

        logical function least_squares_callback(x, fvec, cost, nit, nfev) result(stop)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! current parameter vector
            real(dp), intent(in) :: fvec(:) !! true residual vector at x
            real(dp), intent(in) :: cost !! current robust or linear cost
            integer, intent(in) :: nit !! completed outer iterations
            integer, intent(in) :: nfev !! residual evaluations reported by the solver
        end function least_squares_callback
    end interface

    public :: curve_fit
    public :: curve_model_complex_function
    public :: curve_model_function
    public :: curve_model_jacobian
    public :: least_squares
    public :: least_squares_callback
    public :: least_squares_complex_residual_function
    public :: least_squares_jacobian_function
    public :: least_squares_jacobian_matvec_function
    public :: least_squares_jacobian_rmatvec_function
    public :: least_squares_residual_function
    public :: robust_loss_function
    public :: leastsq

    interface least_squares
        module procedure least_squares_numeric_scale
        module procedure least_squares_jac_scale
    end interface least_squares

contains

    function least_squares_numeric_scale(fun, x0, m, jac, method, ftol, xtol, gtol, max_nfev, &
            x_scale, diff_step, lower, upper, loss, f_scale, loss_fun, diff_scheme, callback, &
            complex_fun, tr_solver, tr_regularize, jac_matvec, jac_rmatvec, jac_sparsity) result(res)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        integer, intent(in) :: m !! residual count
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        character(len=*), intent(in), optional :: method !! trf, dogbox, or lm; default trf
        real(dp), intent(in), optional :: ftol !! relative cost tolerance; default 1e-8
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance; default 1e-8
        real(dp), intent(in), optional :: gtol !! first-order tolerance; default 1e-8
        integer, intent(in), optional :: max_nfev !! maximum reported residual evaluations
        real(dp), intent(in), optional :: x_scale(:) !! positive variable scales
        real(dp), intent(in), optional :: diff_step !! relative finite-difference step
        real(dp), intent(in), optional :: lower(:) !! lower parameter bounds; default -infinity
        real(dp), intent(in), optional :: upper(:) !! upper parameter bounds; default +infinity
        character(len=*), intent(in), optional :: loss !! built-in robust loss name
        real(dp), intent(in), optional :: f_scale !! positive robust-loss scale; default 1
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss
        character(len=*), intent(in), optional :: diff_scheme !! finite-difference scheme
        procedure(least_squares_callback), optional :: callback !! return true to stop after an iteration
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! analytic continuation for cs
        character(len=*), intent(in), optional :: tr_solver !! TRF/dogbox trust solver: exact or lsmr
        logical, intent(in), optional :: tr_regularize !! regularize LSMR subproblem; default true
        procedure(least_squares_jacobian_matvec_function), optional :: jac_matvec !! compute J(x)*v without forming J
        procedure(least_squares_jacobian_rmatvec_function), optional :: jac_rmatvec !! compute transpose(J(x))*w
        logical, intent(in), optional :: jac_sparsity(:, :) !! structural nonzeros for grouped numerical differences
        type(least_squares_result) :: res

        res = least_squares_impl(fun, x0, m, jac, method, ftol, xtol, gtol, max_nfev, &
            x_scale, diff_step, lower, upper, loss, f_scale, loss_fun, diff_scheme, .false., &
            callback, complex_fun, tr_solver, tr_regularize, jac_matvec, jac_rmatvec, jac_sparsity)
    end function least_squares_numeric_scale

    function least_squares_jac_scale(fun, x0, m, jac, method, ftol, xtol, gtol, max_nfev, &
            x_scale, diff_step, lower, upper, loss, f_scale, loss_fun, diff_scheme, callback, &
            complex_fun, tr_solver, tr_regularize, jac_matvec, jac_rmatvec, jac_sparsity) result(res)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        integer, intent(in) :: m !! residual count
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        character(len=*), intent(in), optional :: method !! trf, dogbox, or lm; default trf
        real(dp), intent(in), optional :: ftol !! relative cost tolerance; default 1e-8
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance; default 1e-8
        real(dp), intent(in), optional :: gtol !! first-order tolerance; default 1e-8
        integer, intent(in), optional :: max_nfev !! maximum reported residual evaluations
        character(len=*), intent(in) :: x_scale !! Jacobian scaling selector; must be 'jac'
        real(dp), intent(in), optional :: diff_step !! relative finite-difference step
        real(dp), intent(in), optional :: lower(:) !! lower parameter bounds; default -infinity
        real(dp), intent(in), optional :: upper(:) !! upper parameter bounds; default +infinity
        character(len=*), intent(in), optional :: loss !! built-in robust loss name
        real(dp), intent(in), optional :: f_scale !! positive robust-loss scale; default 1
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss
        character(len=*), intent(in), optional :: diff_scheme !! finite-difference scheme
        procedure(least_squares_callback), optional :: callback !! return true to stop after an iteration
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! analytic continuation for cs
        character(len=*), intent(in), optional :: tr_solver !! TRF/dogbox trust solver: exact or lsmr
        logical, intent(in), optional :: tr_regularize !! regularize LSMR subproblem; default true
        procedure(least_squares_jacobian_matvec_function), optional :: jac_matvec !! matrix-free J(x)*v callback
        procedure(least_squares_jacobian_rmatvec_function), optional :: jac_rmatvec !! matrix-free transpose(J(x))*w callback
        logical, intent(in), optional :: jac_sparsity(:, :) !! structural nonzeros for grouped numerical differences
        type(least_squares_result) :: res

        if (trim(lowercase(x_scale)) /= 'jac') then
            call init_least_squares_result(res, size(x0), m)
            res%message = "character x_scale must be 'jac'"
            return
        end if
        res = least_squares_impl(fun, x0, m, jac, method, ftol, xtol, gtol, max_nfev, &
            diff_step=diff_step, lower=lower, upper=upper, loss=loss, f_scale=f_scale, &
            loss_fun=loss_fun, diff_scheme=diff_scheme, jac_scale=.true., callback=callback, &
            complex_fun=complex_fun, tr_solver=tr_solver, tr_regularize=tr_regularize, &
            jac_matvec=jac_matvec, jac_rmatvec=jac_rmatvec, jac_sparsity=jac_sparsity)
    end function least_squares_jac_scale

    function least_squares_impl(fun, x0, m, jac, method, ftol, xtol, gtol, max_nfev, &
            x_scale, diff_step, lower, upper, loss, f_scale, loss_fun, diff_scheme, jac_scale, &
            callback, complex_fun, tr_solver, tr_regularize, jac_matvec, jac_rmatvec, jac_sparsity) result(res)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        integer, intent(in) :: m !! residual count; positive for TRF/dogbox, >= n for LM
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        character(len=*), intent(in), optional :: method !! trf, dogbox, or lm; default trf
        real(dp), intent(in), optional :: ftol !! relative cost tolerance; default 1e-8
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance; default 1e-8
        real(dp), intent(in), optional :: gtol !! first-order tolerance; default 1e-8
        integer, intent(in), optional :: max_nfev !! maximum reported residual evaluations
        real(dp), intent(in), optional :: x_scale(:) !! positive variable scales
        real(dp), intent(in), optional :: diff_step !! relative finite-difference step
        real(dp), intent(in), optional :: lower(:) !! lower parameter bounds; default -infinity
        real(dp), intent(in), optional :: upper(:) !! upper parameter bounds; default +infinity
        character(len=*), intent(in), optional :: loss !! linear, soft_l1, huber, cauchy, or arctan
        real(dp), intent(in), optional :: f_scale !! positive soft-loss scale; default 1
        procedure(robust_loss_function), optional :: loss_fun !! custom rho(z), rho'(z), rho''(z)
        character(len=*), intent(in), optional :: diff_scheme !! 2-point, 3-point, or cs; default 2-point
        logical, intent(in) :: jac_scale !! use inverse Jacobian-column norms for variable scaling
        procedure(least_squares_callback), optional :: callback !! return true to stop TRF/dogbox iterations
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! analytic continuation for cs
        character(len=*), intent(in), optional :: tr_solver !! TRF/dogbox trust solver: exact or lsmr
        logical, intent(in), optional :: tr_regularize !! regularize LSMR subproblem; default true
        procedure(least_squares_jacobian_matvec_function), optional :: jac_matvec !! matrix-free J(x)*v callback
        procedure(least_squares_jacobian_rmatvec_function), optional :: jac_rmatvec !! matrix-free transpose(J(x))*w callback
        logical, intent(in), optional :: jac_sparsity(:, :) !! structural nonzeros for grouped numerical differences
        type(least_squares_result) :: res

        real(dp), allocatable :: lb(:), ub(:)
        real(dp) :: fscale_use, inf
        integer :: n
        character(len=16) :: diff_scheme_use, loss_use, meth, tr_solver_use
        logical :: operator_mode, tr_regularize_use

        n = size(x0)
        call init_least_squares_result(res, n, m)
        if (n <= 0 .or. m <= 0) return
        if (.not. all(ieee_is_finite(x0))) return

        operator_mode = present(jac_matvec) .or. present(jac_rmatvec)
        if (operator_mode .and. .not. (present(jac_matvec) .and. present(jac_rmatvec))) then
            res%message = 'matrix-free Jacobian requires both jac_matvec and jac_rmatvec'
            return
        end if
        if (operator_mode .and. present(jac)) then
            res%message = 'specify either jac or the matrix-free Jacobian callbacks, not both'
            return
        end if
        if (operator_mode .and. jac_scale) then
            res%message = "x_scale='jac' is not supported with a matrix-free Jacobian"
            return
        end if
        if (present(jac_sparsity)) then
            if (size(jac_sparsity, 1) /= m .or. size(jac_sparsity, 2) /= n) then
                res%message = 'jac_sparsity must have shape (m, size(x0))'
                return
            end if
            if (operator_mode) then
                res%message = 'jac_sparsity is incompatible with matrix-free Jacobian callbacks'
                return
            end if
        end if

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(lb(n), ub(n))
        lb = -inf
        ub = inf
        if (present(lower)) then
            if (size(lower) /= n .or. any(ieee_is_nan(lower))) then
                res%message = 'lower bounds must match x0 and contain no NaNs'
                return
            end if
            lb = lower
        end if
        if (present(upper)) then
            if (size(upper) /= n .or. any(ieee_is_nan(upper))) then
                res%message = 'upper bounds must match x0 and contain no NaNs'
                return
            end if
            ub = upper
        end if
        if (any(lb >= ub)) then
            res%message = 'each lower bound must be strictly less than its upper bound'
            return
        end if
        if (any(x0 < lb) .or. any(x0 > ub)) then
            res%message = 'initial parameter vector is outside the supplied bounds'
            return
        end if

        loss_use = 'linear'
        if (present(loss_fun)) then
            if (present(loss)) then
                res%message = 'specify either loss or loss_fun, not both'
                return
            end if
            loss_use = 'custom'
        else if (present(loss)) then
            loss_use = lowercase(trim(loss))
            select case (trim(loss_use))
            case ('linear', 'soft_l1', 'huber', 'cauchy', 'arctan')
            case default
                res%message = 'loss must be linear, soft_l1, huber, cauchy, or arctan'
                return
            end select
        end if
        fscale_use = 1.0_dp
        if (present(f_scale)) fscale_use = f_scale
        if (.not. ieee_is_finite(fscale_use) .or. fscale_use <= 0.0_dp) then
            res%message = 'f_scale must be finite and positive'
            return
        end if

        diff_scheme_use = '2-point'
        if (present(diff_scheme)) diff_scheme_use = lowercase(trim(diff_scheme))
        if (trim(diff_scheme_use) /= '2-point' .and. trim(diff_scheme_use) /= '3-point' .and. &
                trim(diff_scheme_use) /= 'cs') then
            res%message = 'diff_scheme must be ''2-point'', ''3-point'', or ''cs'''
            return
        end if
        if (trim(diff_scheme_use) == 'cs' .and. .not. present(jac) .and. &
                .not. operator_mode .and. .not. present(complex_fun)) then
            res%message = 'diff_scheme=''cs'' requires complex_fun when jac is absent'
            return
        end if

        tr_solver_use = 'exact'
        if (operator_mode) tr_solver_use = 'lsmr'
        if (present(tr_solver)) tr_solver_use = lowercase(trim(tr_solver))
        if (trim(tr_solver_use) /= 'exact' .and. trim(tr_solver_use) /= 'lsmr') then
            res%message = 'tr_solver must be ''exact'' or ''lsmr'''
            return
        end if
        tr_regularize_use = .true.
        if (present(tr_regularize)) tr_regularize_use = tr_regularize

        meth = 'trf'
        if (present(method)) meth = lowercase(trim(method))
        select case (trim(meth))
        case ('lm')
            if (present(jac_sparsity) .and. .not. present(jac)) then
                res%message = "method='lm' does not support jac_sparsity with numerical Jacobians"
                return
            end if
            if (operator_mode) then
                res%message = "method='lm' does not support matrix-free Jacobians"
                return
            end if
            if (trim(tr_solver_use) /= 'exact') then
                res%message = 'tr_solver is supported only for method=''trf'' or ''dogbox'''
                return
            end if
            if (trim(loss_use) /= 'linear') then
                res%message = 'method=''lm'' supports only loss=''linear'''
                return
            end if
            if (any(ieee_is_finite(lb)) .or. any(ieee_is_finite(ub))) then
                res%message = 'method=''lm'' does not support bounds'
                return
            end if
            res = least_squares_lm_core(fun, x0, m, jac=jac, ftol=ftol, xtol=xtol, &
                gtol=gtol, max_nfev=max_nfev, x_scale=x_scale, diff_step=diff_step, &
                diff_scheme=diff_scheme_use, jac_scale=jac_scale, complex_fun=complex_fun)
        case ('trf')
            if (present(jac_sparsity) .and. .not. present(jac)) then
                if (present(tr_solver) .and. trim(tr_solver_use) == 'exact') then
                    res%message = "tr_solver='exact' is incompatible with jac_sparsity"
                    return
                end if
                tr_solver_use = 'lsmr'
            end if
            if (operator_mode) then
                if (trim(tr_solver_use) /= 'lsmr') then
                    res%message = "tr_solver='exact' is incompatible with a matrix-free Jacobian"
                    return
                end if
                call least_squares_trf_operator_core(fun, x0, m, jac_matvec, jac_rmatvec, &
                    ftol, xtol, gtol, max_nfev, x_scale, lb, ub, loss_use, fscale_use, &
                    loss_fun, callback, res, tr_regularize_use)
            else
                call least_squares_trf_core(fun, x0, m, jac, ftol, xtol, gtol, max_nfev, &
                    x_scale, diff_step, lb, ub, loss_use, fscale_use, loss_fun, &
                    diff_scheme_use, jac_scale, callback, res, complex_fun, tr_solver_use, &
                    tr_regularize_use, jac_sparsity)
            end if
        case ('dogbox')
            if (operator_mode) then
                res%message = "matrix-free Jacobians are currently supported only with method='trf'"
                return
            end if
            if (present(jac_sparsity) .and. .not. present(jac)) then
                if (present(tr_solver) .and. trim(tr_solver_use) == 'exact') then
                    res%message = "tr_solver='exact' is incompatible with jac_sparsity"
                    return
                end if
                tr_solver_use = 'lsmr'
            end if
            call least_squares_dogbox_core(fun, x0, m, jac, ftol, xtol, gtol, max_nfev, &
                x_scale, diff_step, lb, ub, loss_use, fscale_use, loss_fun, &
                diff_scheme_use, jac_scale, callback, res, complex_fun, tr_solver_use, &
                jac_sparsity)
        case default
            res%message = 'method must be ''trf'', ''dogbox'', or ''lm'''
        end select
    end function least_squares_impl

    function least_squares_lm_core(fun, x0, m, jac, ftol, xtol, gtol, &
            max_nfev, x_scale, diff_step, diff_scheme, jac_scale, complex_fun) result(res)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        integer, intent(in) :: m !! residual count; must be at least size(x0)
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        real(dp), intent(in), optional :: ftol !! relative cost tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: gtol !! orthogonality tolerance
        integer, intent(in), optional :: max_nfev !! maximum MINPACK residual evaluations
        real(dp), intent(in), optional :: x_scale(:) !! positive variable scales
        real(dp), intent(in), optional :: diff_step !! relative finite-difference step
        character(len=*), intent(in) :: diff_scheme !! validated finite-difference scheme
        logical, intent(in) :: jac_scale !! use MINPACK internal Jacobian scaling
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! complex residual for cs
        type(least_squares_result) :: res

        procedure(least_squares_residual_function), pointer :: fun_ptr
        procedure(least_squares_jacobian_function), pointer :: jac_ptr
        real(dp), allocatable :: cache_f(:), cache_x(:), diag(:), fjac_native(:, :)
        real(dp), allocatable :: qtf(:), wa1(:), wa2(:), wa3(:), wa4(:)
        real(dp), allocatable :: jac_final(:, :)
        real(dp) :: ftol_use, gtol_use, xtol_use
        integer, allocatable :: ipvt(:)
        integer :: info, maxfev_use, mode, n, nfev, njev
        logical :: cache_valid

        n = size(x0)
        call init_least_squares_result(res, n, m)
        if (n <= 0 .or. m < n) return
        if (.not. all(ieee_is_finite(x0))) return

        res%method = 'lm'

        ftol_use = least_squares_default_tol
        xtol_use = least_squares_default_tol
        gtol_use = least_squares_default_tol
        if (present(ftol)) ftol_use = ftol
        if (present(xtol)) xtol_use = xtol
        if (present(gtol)) gtol_use = gtol
        if (.not. valid_lm_tolerance(ftol_use) .or. .not. valid_lm_tolerance(xtol_use) .or. &
                .not. valid_lm_tolerance(gtol_use)) then
            res%message = 'all LM tolerances must exceed machine epsilon'
            return
        end if

        maxfev_use = 100 * n
        if (present(max_nfev)) then
            if (max_nfev <= 0) return
            maxfev_use = max_nfev
        end if

        allocate(diag(n))
        diag = 1.0_dp
        mode = 1
        if (jac_scale .and. present(x_scale)) return
        if (present(x_scale)) then
            if (size(x_scale) /= n) return
            if (.not. all(ieee_is_finite(x_scale)) .or. any(x_scale <= 0.0_dp)) return
            diag = 1.0_dp / x_scale
            mode = 2
        end if
        if (present(diff_step)) then
            if (.not. ieee_is_finite(diff_step) .or. diff_step < 0.0_dp) return
        end if

        allocate(res%x(n), res%fun(m), fjac_native(m, n), ipvt(n), qtf(n))
        allocate(wa1(n), wa2(n), wa3(n), wa4(m), cache_x(n), cache_f(m))
        res%x = x0
        fun_ptr => fun
        if (present(jac)) jac_ptr => jac
        cache_valid = .false.
        nfev = 0
        njev = 0

        call lmder(lmder_adapter, m, n, res%x, res%fun, fjac_native, m, &
            ftol_use, xtol_use, gtol_use, maxfev_use, diag, mode, default_factor, &
            0, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4)

        allocate(jac_final(m, n))
        if (present(jac)) then
            call jac_ptr(res%x, jac_final)
            res%njev = njev
        else
            call numerical_jacobian(fun_ptr, res%x, res%fun, jac_final, diff_step, &
                diff_scheme, complex_fun)
            res%njev = -1
        end if
        res%jac = jac_final
        allocate(res%grad(n))
        res%grad = matmul(transpose(jac_final), res%fun)
        res%cost = 0.5_dp * dot_product(res%fun, res%fun)
        res%optimality = maxval(abs(res%grad))
        res%nfev = nfev
        res%status = least_squares_status(info)
        res%success = res%status > 0
        res%message = least_squares_message(res%status)

    contains

        subroutine lmder_adapter(m_cb, n_cb, x, f, j, ldj, iflag)
            integer, intent(in) :: m_cb !! MINPACK residual count
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current parameter vector
            real(dp), intent(inout) :: f(m_cb) !! current residual vector
            integer, intent(in) :: ldj !! Jacobian leading dimension
            real(dp), intent(inout) :: j(ldj, n_cb) !! residual Jacobian
            integer, intent(inout) :: iflag !! MINPACK callback selector

            if (iflag == 1) then
                if (cache_valid .and. same_vector(x, cache_x)) then
                    f = cache_f
                else
                    call fun_ptr(x, f)
                    cache_x = x
                    cache_f = f
                    cache_valid = .true.
                end if
            else if (iflag == 2) then
                if (present(jac)) then
                    call jac_ptr(x, j)
                else
                    if (.not. cache_valid .or. .not. same_vector(x, cache_x)) then
                        call fun_ptr(x, cache_f)
                        cache_x = x
                        cache_valid = .true.
                    end if
                    call numerical_jacobian(fun_ptr, x, cache_f, j, diff_step, &
                        diff_scheme, complex_fun)
                end if
            end if
        end subroutine lmder_adapter

    end function least_squares_lm_core

    subroutine least_squares_trf_core(fun, x0, m, jac, ftol, xtol, gtol, max_nfev, &
            x_scale, diff_step, lb, ub, loss_name, f_scale_use, loss_fun, diff_scheme, jac_scale, &
            callback, res, complex_fun, tr_solver, tr_regularize, jac_sparsity)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! feasible starting parameter vector
        integer, intent(in) :: m !! residual count
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        real(dp), intent(in), optional :: ftol !! relative cost tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: gtol !! scaled-gradient tolerance
        integer, intent(in), optional :: max_nfev !! maximum reported residual evaluations
        real(dp), intent(in), optional :: x_scale(:) !! positive variable scales
        real(dp), intent(in), optional :: diff_step !! relative finite-difference step
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(in) :: f_scale_use !! positive robust-loss scale
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss
        character(len=*), intent(in) :: diff_scheme !! validated finite-difference scheme
        logical, intent(in) :: jac_scale !! update scaling from Jacobian column norms
        procedure(least_squares_callback), optional :: callback !! return true to stop after an iteration
        type(least_squares_result), intent(out) :: res !! TRF result
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! complex residual for cs
        character(len=*), intent(in) :: tr_solver !! exact or lsmr trust-region solver
        logical, intent(in) :: tr_regularize !! regularize LSMR subproblem
        logical, intent(in), optional :: jac_sparsity(:, :) !! structural nonzeros for grouped differences

        real(dp), allocatable :: bmat(:, :), bsub(:, :), d(:), diag_h(:), dv(:), f(:), f_new(:)
        real(dp), allocatable :: f_true(:), g(:), g_h(:), jmat(:, :), j_h(:, :)
        real(dp), allocatable :: j_true(:, :), p(:), p_h(:), psub(:), sbasis(:, :)
        real(dp), allocatable :: scale(:), scale_inv(:), step(:), step_h(:), v(:), gsub(:)
        real(dp), allocatable :: x(:), x_new(:)
        real(dp) :: actual_reduction, alpha, alpha_new, cost, cost_new, delta, delta_new
        real(dp) :: ftol_use, g_norm, gtol_use, predicted_reduction, ratio
        real(dp) :: step_h_norm, step_norm, theta, xtol_use
        integer :: i, iteration, ksub, maxfev_use, n, nfev, njev, status, term_status
        logical :: ok

        n = size(x0)
        call init_least_squares_result(res, n, m)
        res%method = 'trf'
        if (n <= 0 .or. m <= 0) return

        ftol_use = least_squares_default_tol
        xtol_use = least_squares_default_tol
        gtol_use = least_squares_default_tol
        if (present(ftol)) ftol_use = ftol
        if (present(xtol)) xtol_use = xtol
        if (present(gtol)) gtol_use = gtol
        if (.not. valid_trf_tolerances(ftol_use, xtol_use, gtol_use)) then
            res%message = 'TRF tolerances must be nonnegative and at least one must exceed epsilon'
            return
        end if

        maxfev_use = 100 * n
        if (present(max_nfev)) then
            if (max_nfev <= 0) return
            maxfev_use = max_nfev
        end if
        if (present(diff_step)) then
            if (.not. ieee_is_finite(diff_step) .or. diff_step < 0.0_dp) return
        end if

        allocate(scale(n), scale_inv(n))
        scale = 1.0_dp
        if (jac_scale .and. present(x_scale)) return
        if (present(x_scale)) then
            if (size(x_scale) /= n) return
            if (.not. all(ieee_is_finite(x_scale)) .or. any(x_scale <= 0.0_dp)) return
            scale = x_scale
        end if
        scale_inv = 1.0_dp / scale

        allocate(x(n), x_new(n), f(m), f_new(m), f_true(m), jmat(m, n), j_true(m, n))
        allocate(g(n), v(n), dv(n))
        allocate(d(n), diag_h(n), g_h(n), j_h(m, n), bmat(n, n), p(n), p_h(n))
        allocate(step(n), step_h(n), sbasis(n, 2), bsub(2, 2), gsub(2), psub(2))
        x = x0
        if (any(ieee_is_finite(lb)) .or. any(ieee_is_finite(ub))) then
            call make_strictly_feasible(x, lb, ub, 1.0e-10_dp, x_new)
            x = x_new
        end if
        call fun(x, f_true)
        if (.not. all(ieee_is_finite(f_true))) then
            res%message = 'residual callback returned a nonfinite initial residual'
            return
        end if
        nfev = 1
        if (present(jac)) then
            call jac(x, j_true)
        else
            call numerical_jacobian_bounded(fun, x, f_true, j_true, diff_step, lb, ub, &
                diff_scheme, complex_fun, jac_sparsity)
        end if
        if (.not. all(ieee_is_finite(j_true))) then
            res%message = 'Jacobian callback returned a nonfinite value'
            return
        end if
        njev = 1
        call scale_robust_model(f_true, j_true, loss_name, f_scale_use, f, jmat, cost, loss_fun)
        g = matmul(transpose(jmat), f)
        if (jac_scale) call update_jac_scale(jmat, scale, scale_inv, .false.)

        call cl_scaling_vector(x, g, lb, ub, v, dv)
        do i = 1, n
            if (.not. same_real(dv(i), 0.0_dp)) v(i) = v(i) * scale_inv(i)
        end do
        delta = vector_norm2(x * scale_inv / sqrt(v))
        if (same_real(delta, 0.0_dp)) delta = 1.0_dp

        alpha = 0.0_dp
        status = -99
        iteration = 0
        actual_reduction = 0.0_dp
        step_norm = 0.0_dp
        g_norm = maxval(abs(g * v))

        do
            call cl_scaling_vector(x, g, lb, ub, v, dv)
            g_norm = maxval(abs(g * v))
            if (g_norm < gtol_use) status = 1
            if (status /= -99 .or. nfev >= maxfev_use) exit

            do i = 1, n
                if (.not. same_real(dv(i), 0.0_dp)) v(i) = v(i) * scale_inv(i)
            end do
            d = sqrt(v) * scale
            diag_h = g * dv * scale
            g_h = d * g
            do i = 1, n
                j_h(:, i) = jmat(:, i) * d(i)
            end do
            if (trim(tr_solver) == 'exact') then
                bmat = matmul(transpose(j_h), j_h)
                do i = 1, n
                    bmat(i, i) = bmat(i, i) + diag_h(i)
                end do
            else
                call build_lsmr_subspace(j_h, diag_h, g_h, f, delta, tr_regularize, &
                    sbasis, ksub, bsub, gsub, ok)
                if (.not. ok) then
                    res%message = 'LSMR TRF subspace construction failed'
                    return
                end if
            end if
            theta = max(0.995_dp, 1.0_dp - g_norm)

            actual_reduction = -1.0_dp
            do while (actual_reduction <= 0.0_dp .and. nfev < maxfev_use)
                if (trim(tr_solver) == 'exact') then
                    call solve_trust_region_dense(bmat, g_h, delta, alpha, p_h, alpha_new, ok)
                    alpha = alpha_new
                else
                    psub = 0.0_dp
                    call solve_trust_region_dense(bsub(1:ksub, 1:ksub), gsub(1:ksub), &
                        delta, 0.0_dp, psub(1:ksub), alpha_new, ok)
                    p_h = matmul(sbasis(:, 1:ksub), psub(1:ksub))
                end if
                if (.not. ok) then
                    res%message = 'TRF trust-region subproblem did not converge'
                    return
                end if
                p = d * p_h
                call select_trf_step(x, j_h, diag_h, g_h, p, p_h, d, delta, &
                    lb, ub, theta, step, step_h, predicted_reduction)
                call make_strictly_feasible(x + step, lb, ub, 0.0_dp, x_new)
                call fun(x_new, f_new)
                nfev = nfev + 1
                step_h_norm = vector_norm2(step_h)
                if (.not. all(ieee_is_finite(f_new))) then
                    delta = 0.25_dp * step_h_norm
                    cycle
                end if

                cost_new = robust_cost(f_new, loss_name, f_scale_use, loss_fun)
                actual_reduction = cost - cost_new
                call update_tr_radius(delta, actual_reduction, predicted_reduction, &
                    step_h_norm, step_h_norm > 0.95_dp * delta, delta_new, ratio)
                step_norm = vector_norm2(step)
                call check_ls_termination(actual_reduction, cost, step_norm, &
                    vector_norm2(x), ratio, ftol_use, xtol_use, term_status)
                if (term_status /= -99) then
                    status = term_status
                    exit
                end if
                if (trim(tr_solver) == 'exact' .and. delta_new > 0.0_dp) then
                    alpha = alpha * delta / delta_new
                end if
                delta = delta_new
            end do

            if (actual_reduction > 0.0_dp) then
                x = x_new
                f_true = f_new
                cost = cost_new
                if (present(jac)) then
                    call jac(x, j_true)
                else
                    call numerical_jacobian_bounded(fun, x, f_true, j_true, diff_step, lb, ub, &
                        diff_scheme, complex_fun, jac_sparsity)
                end if
                if (.not. all(ieee_is_finite(j_true))) then
                    res%message = 'Jacobian callback returned a nonfinite value'
                    return
                end if
                njev = njev + 1
                call scale_robust_model(f_true, j_true, loss_name, f_scale_use, f, jmat, cost, loss_fun)
                g = matmul(transpose(jmat), f)
                if (jac_scale) call update_jac_scale(jmat, scale, scale_inv, .true.)
            else
                step_norm = 0.0_dp
                actual_reduction = 0.0_dp
            end if

            iteration = iteration + 1
            if (present(callback)) then
                if (callback(x, f_true, cost, iteration, nfev)) then
                    status = -2
                    exit
                end if
            end if
        end do

        if (status == -99) status = 0
        res%x = x
        res%fun = f_true
        res%jac = jmat
        res%grad = g
        res%cost = cost
        res%optimality = g_norm
        call find_active_constraints(x, lb, ub, xtol_use, res%active_mask)
        res%nfev = nfev
        res%njev = njev
        res%status = status
        res%success = status > 0
        res%message = least_squares_message(status)
    end subroutine least_squares_trf_core

    subroutine least_squares_trf_operator_core(fun, x0, m, jac_matvec, jac_rmatvec, &
            ftol, xtol, gtol, max_nfev, x_scale, lb, ub, loss_name, f_scale_use, &
            loss_fun, callback, res, tr_regularize)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! feasible starting parameter vector
        integer, intent(in) :: m !! residual count
        procedure(least_squares_jacobian_matvec_function) :: jac_matvec !! compute J(x)*v
        procedure(least_squares_jacobian_rmatvec_function) :: jac_rmatvec !! compute transpose(J(x))*w
        real(dp), intent(in), optional :: ftol !! relative cost tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: gtol !! scaled-gradient tolerance
        integer, intent(in), optional :: max_nfev !! maximum reported residual evaluations
        real(dp), intent(in), optional :: x_scale(:) !! positive variable scales
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(in) :: f_scale_use !! positive robust-loss scale
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss
        procedure(least_squares_callback), optional :: callback !! return true to stop after an iteration
        type(least_squares_result), intent(out) :: res !! TRF result
        logical, intent(in) :: tr_regularize !! regularize LSMR subproblem

        real(dp), allocatable :: bsub(:, :), d(:), diag_h(:), dv(:), f(:), f_new(:)
        real(dp), allocatable :: f_true(:), g(:), g_h(:), gsub(:), op_work_m(:), op_work_n(:)
        real(dp), allocatable :: p(:), p_h(:), psub(:), row_scale(:), sbasis(:, :)
        real(dp), allocatable :: scale(:), scale_inv(:), step(:), step_h(:), v(:)
        real(dp), allocatable :: x(:), x_new(:)
        real(dp) :: actual_reduction, alpha_new, cost, cost_new, delta, delta_new
        real(dp) :: ftol_use, g_norm, gtol_use, predicted_reduction, ratio
        real(dp) :: step_h_norm, step_norm, theta, xtol_use
        integer :: i, iteration, ksub, maxfev_use, n, nfev, njev, status, term_status
        logical :: ok

        n = size(x0)
        call init_least_squares_result(res, n, m)
        res%method = 'trf'
        res%jac_is_operator = .true.
        if (n <= 0 .or. m <= 0) return

        ftol_use = least_squares_default_tol
        xtol_use = least_squares_default_tol
        gtol_use = least_squares_default_tol
        if (present(ftol)) ftol_use = ftol
        if (present(xtol)) xtol_use = xtol
        if (present(gtol)) gtol_use = gtol
        if (.not. valid_trf_tolerances(ftol_use, xtol_use, gtol_use)) then
            res%message = 'TRF tolerances must be nonnegative and at least one must exceed epsilon'
            return
        end if

        maxfev_use = 100 * n
        if (present(max_nfev)) then
            if (max_nfev <= 0) return
            maxfev_use = max_nfev
        end if

        allocate(scale(n), scale_inv(n))
        scale = 1.0_dp
        if (present(x_scale)) then
            if (size(x_scale) /= n) return
            if (.not. all(ieee_is_finite(x_scale)) .or. any(x_scale <= 0.0_dp)) return
            scale = x_scale
        end if
        scale_inv = 1.0_dp / scale

        allocate(x(n), x_new(n), f(m), f_new(m), f_true(m), row_scale(m))
        allocate(g(n), v(n), dv(n), d(n), diag_h(n), g_h(n))
        allocate(p(n), p_h(n), step(n), step_h(n), op_work_m(m), op_work_n(n))
        allocate(sbasis(n, 2), bsub(2, 2), gsub(2), psub(2))
        x = x0
        if (any(ieee_is_finite(lb)) .or. any(ieee_is_finite(ub))) then
            call make_strictly_feasible(x, lb, ub, 1.0e-10_dp, x_new)
            x = x_new
        end if
        call fun(x, f_true)
        if (.not. all(ieee_is_finite(f_true))) then
            res%message = 'residual callback returned a nonfinite initial residual'
            return
        end if
        nfev = 1
        njev = 1
        call scale_robust_operator_model(f_true, loss_name, f_scale_use, f, row_scale, &
            cost, loss_fun)
        call model_rmatvec(f, g)
        if (.not. all(ieee_is_finite(g))) then
            res%message = 'matrix-free Jacobian transpose product returned a nonfinite value'
            return
        end if

        call cl_scaling_vector(x, g, lb, ub, v, dv)
        do i = 1, n
            if (.not. same_real(dv(i), 0.0_dp)) v(i) = v(i) * scale_inv(i)
        end do
        delta = vector_norm2(x * scale_inv / sqrt(v))
        if (same_real(delta, 0.0_dp)) delta = 1.0_dp

        status = -99
        iteration = 0
        actual_reduction = 0.0_dp
        step_norm = 0.0_dp
        g_norm = maxval(abs(g * v))

        do
            call cl_scaling_vector(x, g, lb, ub, v, dv)
            g_norm = maxval(abs(g * v))
            if (g_norm < gtol_use) status = 1
            if (status /= -99 .or. nfev >= maxfev_use) exit

            do i = 1, n
                if (.not. same_real(dv(i), 0.0_dp)) v(i) = v(i) * scale_inv(i)
            end do
            d = sqrt(v) * scale
            diag_h = g * dv * scale
            g_h = d * g
            call build_lsmr_subspace_operator(j_h_matvec, j_h_rmatvec, m, n, diag_h, &
                g_h, f, delta, tr_regularize, sbasis, ksub, bsub, gsub, ok)
            if (.not. ok) then
                res%message = 'matrix-free LSMR TRF subspace construction failed'
                return
            end if
            theta = max(0.995_dp, 1.0_dp - g_norm)

            actual_reduction = -1.0_dp
            do while (actual_reduction <= 0.0_dp .and. nfev < maxfev_use)
                psub = 0.0_dp
                call solve_trust_region_dense(bsub(1:ksub, 1:ksub), gsub(1:ksub), &
                    delta, 0.0_dp, psub(1:ksub), alpha_new, ok)
                p_h = matmul(sbasis(:, 1:ksub), psub(1:ksub))
                if (.not. ok) then
                    res%message = 'TRF trust-region subproblem did not converge'
                    return
                end if
                p = d * p_h
                call select_trf_step_operator(x, j_h_matvec, m, diag_h, g_h, p, p_h, &
                    d, delta, lb, ub, theta, step, step_h, predicted_reduction)
                call make_strictly_feasible(x + step, lb, ub, 0.0_dp, x_new)
                call fun(x_new, f_new)
                nfev = nfev + 1
                step_h_norm = vector_norm2(step_h)
                if (.not. all(ieee_is_finite(f_new))) then
                    delta = 0.25_dp * step_h_norm
                    cycle
                end if

                cost_new = robust_cost(f_new, loss_name, f_scale_use, loss_fun)
                actual_reduction = cost - cost_new
                call update_tr_radius(delta, actual_reduction, predicted_reduction, &
                    step_h_norm, step_h_norm > 0.95_dp * delta, delta_new, ratio)
                step_norm = vector_norm2(step)
                call check_ls_termination(actual_reduction, cost, step_norm, &
                    vector_norm2(x), ratio, ftol_use, xtol_use, term_status)
                if (term_status /= -99) then
                    status = term_status
                    exit
                end if
                delta = delta_new
            end do

            if (actual_reduction > 0.0_dp) then
                x = x_new
                f_true = f_new
                cost = cost_new
                njev = njev + 1
                call scale_robust_operator_model(f_true, loss_name, f_scale_use, f, &
                    row_scale, cost, loss_fun)
                call model_rmatvec(f, g)
                if (.not. all(ieee_is_finite(g))) then
                    res%message = 'matrix-free Jacobian transpose product returned a nonfinite value'
                    return
                end if
            else
                step_norm = 0.0_dp
                actual_reduction = 0.0_dp
            end if

            iteration = iteration + 1
            if (present(callback)) then
                if (callback(x, f_true, cost, iteration, nfev)) then
                    status = -2
                    exit
                end if
            end if
        end do

        if (status == -99) status = 0
        res%x = x
        res%fun = f_true
        res%grad = g
        res%cost = cost
        res%optimality = g_norm
        call find_active_constraints(x, lb, ub, xtol_use, res%active_mask)
        res%nfev = nfev
        res%njev = njev
        res%status = status
        res%success = status > 0
        res%message = least_squares_message(status)

    contains

        subroutine model_matvec(vin, yout)
            real(dp), intent(in) :: vin(:) !! parameter-space vector
            real(dp), intent(out) :: yout(:) !! robust-scaled residual-space product

            call jac_matvec(x, vin, yout)
            yout = row_scale * yout
        end subroutine model_matvec

        subroutine model_rmatvec(win, yout)
            real(dp), intent(in) :: win(:) !! robust-scaled residual-space vector
            real(dp), intent(out) :: yout(:) !! parameter-space transpose product

            op_work_m = row_scale * win
            call jac_rmatvec(x, op_work_m, yout)
        end subroutine model_rmatvec

        subroutine j_h_matvec(vin, yout)
            real(dp), intent(in) :: vin(:) !! trust-region scaled parameter vector
            real(dp), intent(out) :: yout(:) !! residual-space product

            op_work_n = d * vin
            call model_matvec(op_work_n, yout)
        end subroutine j_h_matvec

        subroutine j_h_rmatvec(win, yout)
            real(dp), intent(in) :: win(:) !! residual-space vector
            real(dp), intent(out) :: yout(:) !! trust-region scaled transpose product

            call model_rmatvec(win, op_work_n)
            yout = d * op_work_n
        end subroutine j_h_rmatvec

    end subroutine least_squares_trf_operator_core

    subroutine least_squares_dogbox_core(fun, x0, m, jac, ftol, xtol, gtol, max_nfev, &
            x_scale, diff_step, lb, ub, loss_name, f_scale_use, loss_fun, diff_scheme, jac_scale, &
            callback, res, complex_fun, tr_solver, jac_sparsity)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! feasible starting parameter vector
        integer, intent(in) :: m !! residual count
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        real(dp), intent(in), optional :: ftol !! relative cost tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: gtol !! infinity-norm gradient tolerance
        integer, intent(in), optional :: max_nfev !! maximum reported residual evaluations
        real(dp), intent(in), optional :: x_scale(:) !! positive variable scales
        real(dp), intent(in), optional :: diff_step !! relative finite-difference step
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(in) :: f_scale_use !! positive robust-loss scale
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss
        character(len=*), intent(in) :: diff_scheme !! validated finite-difference scheme
        logical, intent(in) :: jac_scale !! update scaling from Jacobian column norms
        procedure(least_squares_callback), optional :: callback !! return true to stop after an iteration
        type(least_squares_result), intent(out) :: res !! dogbox result
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! complex residual for cs
        character(len=*), intent(in) :: tr_solver !! exact or lsmr Gauss-Newton solver
        logical, intent(in), optional :: jac_sparsity(:, :) !! structural nonzeros for grouped differences

        logical, allocatable :: active(:)
        integer, allocatable :: bound_hits(:), on_bound(:)
        real(dp), allocatable :: f(:), f_new(:), f_true(:), g(:), g_full(:)
        real(dp), allocatable :: jmat(:, :), j_true(:, :), jwork(:, :)
        real(dp), allocatable :: newton_step(:), scale(:), scale_inv(:), step(:)
        real(dp), allocatable :: tr_bounds(:), x(:), x_new(:)
        real(dp) :: a, actual_reduction, b, cost, cost_new, delta, delta_new, ftol_use
        real(dp) :: g_norm, gtol_use, predicted_reduction, ratio, step_h_norm
        real(dp) :: step_norm, xtol_use
        integer :: i, iteration, maxfev_use, n, nfev, njev, status, term_status
        logical :: ok, tr_hit

        n = size(x0)
        call init_least_squares_result(res, n, m)
        res%method = 'dogbox'
        if (n <= 0 .or. m <= 0) return

        ftol_use = least_squares_default_tol
        xtol_use = least_squares_default_tol
        gtol_use = least_squares_default_tol
        if (present(ftol)) ftol_use = ftol
        if (present(xtol)) xtol_use = xtol
        if (present(gtol)) gtol_use = gtol
        if (.not. valid_trf_tolerances(ftol_use, xtol_use, gtol_use)) then
            res%message = 'dogbox tolerances must be nonnegative and one must exceed epsilon'
            return
        end if
        maxfev_use = 100 * n
        if (present(max_nfev)) then
            if (max_nfev <= 0) return
            maxfev_use = max_nfev
        end if
        if (present(diff_step)) then
            if (.not. ieee_is_finite(diff_step) .or. diff_step < 0.0_dp) return
        end if

        allocate(scale(n), scale_inv(n))
        scale = 1.0_dp
        if (jac_scale .and. present(x_scale)) return
        if (present(x_scale)) then
            if (size(x_scale) /= n) return
            if (.not. all(ieee_is_finite(x_scale)) .or. any(x_scale <= 0.0_dp)) return
            scale = x_scale
        end if
        scale_inv = 1.0_dp / scale

        allocate(x(n), x_new(n), f(m), f_new(m), f_true(m), jmat(m, n), j_true(m, n), jwork(m, n))
        allocate(g(n), g_full(n), newton_step(n), step(n), tr_bounds(n))
        allocate(active(n), bound_hits(n), on_bound(n))
        x = x0
        call fun(x, f_true)
        if (.not. all(ieee_is_finite(f_true))) then
            res%message = 'residual callback returned a nonfinite initial residual'
            return
        end if
        nfev = 1
        if (present(jac)) then
            call jac(x, j_true)
        else
            call numerical_jacobian_bounded(fun, x, f_true, j_true, diff_step, lb, ub, &
                diff_scheme, complex_fun, jac_sparsity)
        end if
        if (.not. all(ieee_is_finite(j_true))) then
            res%message = 'Jacobian callback returned a nonfinite value'
            return
        end if
        njev = 1
        call scale_robust_model(f_true, j_true, loss_name, f_scale_use, f, jmat, cost, loss_fun)
        g = matmul(transpose(jmat), f)
        if (jac_scale) call update_jac_scale(jmat, scale, scale_inv, .false.)
        delta = maxval(abs(x * scale_inv))
        if (same_real(delta, 0.0_dp)) delta = 1.0_dp
        on_bound = 0
        do i = 1, n
            if (same_real(x(i), lb(i))) on_bound(i) = -1
            if (same_real(x(i), ub(i))) on_bound(i) = 1
        end do

        status = -99
        iteration = 0
        actual_reduction = 0.0_dp
        step_norm = 0.0_dp
        g_norm = maxval(abs(g))
        g_full = g

        do
            active = on_bound * g < 0.0_dp
            g_full = g
            where (active) g = 0.0_dp
            g_norm = maxval(abs(g))
            if (g_norm < gtol_use) status = 1
            if (status /= -99 .or. nfev >= maxfev_use) exit

            jwork = jmat
            do i = 1, n
                if (active(i)) jwork(:, i) = 0.0_dp
            end do
            if (trim(tr_solver) == 'exact') then
                call solve_linear_least_squares(jwork, -f, newton_step, ok)
                if (.not. ok) then
                    res%message = 'dense dogbox least-squares eigensolve did not converge'
                    return
                end if
            else
                call dogbox_lsmr_step(jmat, f, scale, active, newton_step, ok)
                if (.not. ok) then
                    res%message = 'dogbox LSMR least-squares solve returned a nonfinite step'
                    return
                end if
            end if
            call build_quadratic_1d(jwork, g, -g, a, b)

            actual_reduction = -1.0_dp
            do while (actual_reduction <= 0.0_dp .and. nfev < maxfev_use)
                tr_bounds = delta * scale
                where (active) tr_bounds = 0.0_dp
                call dogleg_step_dense(x, newton_step, g, a, b, tr_bounds, lb, ub, &
                    step, bound_hits, tr_hit)
                predicted_reduction = -evaluate_quadratic(jwork, g, step)
                x_new = min(max(x + step, lb), ub)
                call fun(x_new, f_new)
                nfev = nfev + 1
                step_h_norm = maxval(abs(step * scale_inv))
                if (.not. all(ieee_is_finite(f_new))) then
                    delta = 0.25_dp * step_h_norm
                    cycle
                end if
                cost_new = robust_cost(f_new, loss_name, f_scale_use, loss_fun)
                actual_reduction = cost - cost_new
                call update_tr_radius(delta, actual_reduction, predicted_reduction, &
                    step_h_norm, tr_hit, delta_new, ratio)
                delta = delta_new
                step_norm = vector_norm2(step)
                call check_ls_termination(actual_reduction, cost, step_norm, &
                    vector_norm2(x), ratio, ftol_use, xtol_use, term_status)
                if (term_status /= -99) then
                    status = term_status
                    exit
                end if
            end do

            if (actual_reduction > 0.0_dp) then
                do i = 1, n
                    if (.not. active(i)) on_bound(i) = bound_hits(i)
                end do
                x = x_new
                do i = 1, n
                    if (on_bound(i) == -1) x(i) = lb(i)
                    if (on_bound(i) == 1) x(i) = ub(i)
                end do
                f_true = f_new
                cost = cost_new
                if (present(jac)) then
                    call jac(x, j_true)
                else
                    call numerical_jacobian_bounded(fun, x, f_true, j_true, diff_step, lb, ub, &
                        diff_scheme, complex_fun, jac_sparsity)
                end if
                if (.not. all(ieee_is_finite(j_true))) then
                    res%message = 'Jacobian callback returned a nonfinite value'
                    return
                end if
                njev = njev + 1
                call scale_robust_model(f_true, j_true, loss_name, f_scale_use, f, jmat, cost, loss_fun)
                g = matmul(transpose(jmat), f)
                if (jac_scale) call update_jac_scale(jmat, scale, scale_inv, .true.)
            else
                step_norm = 0.0_dp
                actual_reduction = 0.0_dp
            end if

            iteration = iteration + 1
            if (present(callback)) then
                if (callback(x, f_true, cost, iteration, nfev)) then
                    status = -2
                    exit
                end if
            end if
        end do

        if (status == -99) status = 0
        res%x = x
        res%fun = f_true
        res%jac = jmat
        res%grad = g_full
        res%cost = cost
        res%optimality = g_norm
        res%active_mask = on_bound
        res%nfev = nfev
        res%njev = njev
        res%status = status
        res%success = status > 0
        res%message = least_squares_message(status)
    end subroutine least_squares_dogbox_core

    subroutine dogbox_lsmr_step(jmat, fvec, scale, active, step, ok)
        real(dp), intent(in) :: jmat(:, :) !! robust-scaled Jacobian
        real(dp), intent(in) :: fvec(:) !! robust-scaled residual vector
        real(dp), intent(in) :: scale(:) !! dogbox variable scale
        logical, intent(in) :: active(:) !! variables fixed by bound optimality
        real(dp), intent(out) :: step(:) !! approximate Gauss-Newton step
        logical, intent(out) :: ok !! true when LSMR returned a finite step

        real(dp), allocatable :: scaled_step(:)
        real(dp), allocatable :: work_n(:)
        integer :: m
        integer :: n

        m = size(jmat, 1)
        n = size(jmat, 2)
        ok = .false.
        step = 0.0_dp
        if (size(fvec) /= m .or. size(scale) /= n .or. size(active) /= n .or. &
                size(step) /= n) return

        allocate(scaled_step(n), work_n(n))
        call lsmr_solve(lsmr_matvec, lsmr_rmatvec, m, n, fvec, scaled_step)
        if (.not. all(ieee_is_finite(scaled_step))) return

        step = -scale * scaled_step
        where (active) step = 0.0_dp
        ok = all(ieee_is_finite(step))

    contains

        subroutine lsmr_matvec(v, y)
            real(dp), intent(in) :: v(:) !! scaled parameter-space vector
            real(dp), intent(out) :: y(:) !! residual-space product

            work_n = scale * v
            where (active) work_n = 0.0_dp
            y = matmul(jmat, work_n)
        end subroutine lsmr_matvec

        subroutine lsmr_rmatvec(w, y)
            real(dp), intent(in) :: w(:) !! residual-space vector
            real(dp), intent(out) :: y(:) !! scaled transposed-Jacobian product

            y = scale * matmul(transpose(jmat), w)
            where (active) y = 0.0_dp
        end subroutine lsmr_rmatvec

    end subroutine dogbox_lsmr_step

    function leastsq(fun, x0, m, jac, ftol, xtol, gtol, maxfev, epsfcn, factor, &
            diag) result(res)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting parameter vector
        integer, intent(in) :: m !! residual count; must be at least size(x0)
        procedure(least_squares_jacobian_function), optional :: jac !! analytic Jacobian
        real(dp), intent(in), optional :: ftol !! relative cost tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: gtol !! orthogonality tolerance
        integer, intent(in), optional :: maxfev !! maximum residual evaluations
        real(dp), intent(in), optional :: epsfcn !! lmdif finite-difference error estimate
        real(dp), intent(in), optional :: factor !! initial trust-region scale
        real(dp), intent(in), optional :: diag(:) !! positive variable scaling factors
        type(leastsq_result) :: res

        procedure(least_squares_residual_function), pointer :: fun_ptr
        procedure(least_squares_jacobian_function), pointer :: jac_ptr
        real(dp), allocatable :: diag_use(:), fjac_native(:, :), qtf(:)
        real(dp), allocatable :: wa1(:), wa2(:), wa3(:), wa4(:)
        real(dp) :: eps_use, factor_use, ftol_use, gtol_use, xtol_use
        integer, allocatable :: ipvt(:)
        integer :: info, maxfev_use, mode, n, nfev, njev
        logical :: cov_ok

        n = size(x0)
        call init_leastsq_result(res)
        if (n <= 0 .or. m < n) return
        if (.not. all(ieee_is_finite(x0))) return

        ftol_use = leastsq_default_tol
        xtol_use = leastsq_default_tol
        gtol_use = 0.0_dp
        if (present(ftol)) ftol_use = ftol
        if (present(xtol)) xtol_use = xtol
        if (present(gtol)) gtol_use = gtol
        if (.not. ieee_is_finite(ftol_use) .or. ftol_use < 0.0_dp) return
        if (.not. ieee_is_finite(xtol_use) .or. xtol_use < 0.0_dp) return
        if (.not. ieee_is_finite(gtol_use) .or. gtol_use < 0.0_dp) return

        eps_use = epsilon(1.0_dp)
        if (present(epsfcn)) eps_use = epsfcn
        if (.not. ieee_is_finite(eps_use) .or. eps_use < 0.0_dp) return
        factor_use = default_factor
        if (present(factor)) factor_use = factor
        if (.not. ieee_is_finite(factor_use) .or. factor_use <= 0.1_dp .or. &
                factor_use > 100.0_dp) return

        maxfev_use = merge(100 * (n + 1), 200 * (n + 1), present(jac))
        if (present(maxfev)) then
            if (maxfev < 0) return
            if (maxfev > 0) maxfev_use = maxfev
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

        allocate(res%x(n), res%fvec(m), fjac_native(m, n), ipvt(n), qtf(n))
        allocate(wa1(n), wa2(n), wa3(n), wa4(m))
        res%x = x0
        fun_ptr => fun
        if (present(jac)) jac_ptr => jac
        nfev = 0
        njev = 0

        if (present(jac)) then
            call lmder(lmder_adapter, m, n, res%x, res%fvec, fjac_native, m, &
                ftol_use, xtol_use, gtol_use, maxfev_use, diag_use, mode, factor_use, &
                0, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4)
            res%njev = njev
        else
            call lmdif(lmdif_adapter, m, n, res%x, res%fvec, ftol_use, xtol_use, &
                gtol_use, maxfev_use, eps_use, diag_use, mode, factor_use, 0, &
                info, nfev, fjac_native, m, ipvt, qtf, wa1, wa2, wa3, wa4)
            res%njev = -1
        end if

        res%fjac = transpose(fjac_native)
        res%ipvt = ipvt
        res%qtf = qtf
        res%nfev = nfev
        res%status = info
        res%success = info >= 1 .and. info <= 4
        res%message = leastsq_message(info, ftol_use, xtol_use, gtol_use, maxfev_use)
        if (res%success) then
            call covariance_from_minpack(fjac_native, ipvt, res%cov_x, cov_ok)
            if (.not. cov_ok .and. allocated(res%cov_x)) deallocate(res%cov_x)
        end if

    contains

        subroutine lmdif_adapter(m_cb, n_cb, x, f, iflag)
            integer, intent(in) :: m_cb !! MINPACK residual count
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current parameter vector
            real(dp), intent(out) :: f(m_cb) !! residual vector
            integer, intent(inout) :: iflag !! MINPACK callback selector

            if (iflag < 0) return
            call fun_ptr(x, f)
        end subroutine lmdif_adapter

        subroutine lmder_adapter(m_cb, n_cb, x, f, j, ldj, iflag)
            integer, intent(in) :: m_cb !! MINPACK residual count
            integer, intent(in) :: n_cb !! MINPACK variable count
            real(dp), intent(in) :: x(n_cb) !! current parameter vector
            real(dp), intent(inout) :: f(m_cb) !! residual vector
            integer, intent(in) :: ldj !! Jacobian leading dimension
            real(dp), intent(inout) :: j(ldj, n_cb) !! residual Jacobian
            integer, intent(inout) :: iflag !! MINPACK callback selector

            if (iflag == 1) then
                call fun_ptr(x, f)
            else if (iflag == 2) then
                call jac_ptr(x, j)
            end if
        end subroutine lmder_adapter

    end function leastsq

    function curve_fit(model, xdata, ydata, p0, jac, sigma, absolute_sigma, ftol, &
            xtol, gtol, maxfev, factor, diag, method, lower, upper, loss, f_scale, &
            diff_step, diff_scheme, complex_model) result(res)
        procedure(curve_model_function) :: model !! model values for xdata and parameters
        real(dp), intent(in) :: xdata(:) !! independent data vector
        real(dp), intent(in) :: ydata(:) !! observed dependent data vector
        real(dp), intent(in) :: p0(:) !! required initial parameter vector
        procedure(curve_model_jacobian), optional :: jac !! analytic model Jacobian
        real(dp), intent(in), optional :: sigma(:) !! positive one-sigma errors for ydata
        logical, intent(in), optional :: absolute_sigma !! retain absolute sigma covariance scale
        real(dp), intent(in), optional :: ftol !! relative cost tolerance
        real(dp), intent(in), optional :: xtol !! relative parameter tolerance
        real(dp), intent(in), optional :: gtol !! first-order tolerance
        integer, intent(in), optional :: maxfev !! maximum reported residual evaluations
        real(dp), intent(in), optional :: factor !! LM initial trust-region scale
        real(dp), intent(in), optional :: diag(:) !! positive LM parameter scales
        character(len=*), intent(in), optional :: method !! lm, trf, or dogbox
        real(dp), intent(in), optional :: lower(:) !! lower parameter bounds
        real(dp), intent(in), optional :: upper(:) !! upper parameter bounds
        character(len=*), intent(in), optional :: loss !! robust loss for TRF/dogbox
        real(dp), intent(in), optional :: f_scale !! positive robust-loss scale
        real(dp), intent(in), optional :: diff_step !! relative difference step for TRF/dogbox
        character(len=*), intent(in), optional :: diff_scheme !! 2-point, 3-point, or cs
        procedure(curve_model_complex_function), optional :: complex_model !! analytic continuation for cs
        type(curve_fit_result) :: res

        type(leastsq_result) :: lsres
        type(least_squares_result) :: trres
        procedure(curve_model_function), pointer :: model_ptr
        procedure(curve_model_complex_function), pointer :: complex_model_ptr
        procedure(curve_model_jacobian), pointer :: jac_ptr
        real(dp), allocatable :: lb(:), pcov_base(:, :), ub(:), weights(:), ywork(:)
        real(dp) :: fscale_use, inf, scale
        integer :: m, n
        logical :: absolute_sigma_use, cov_ok, has_bounds
        character(len=16) :: diff_scheme_use, loss_use, meth

        m = size(ydata)
        n = size(p0)
        if (m <= 0 .or. n <= 0 .or. size(xdata) /= m) return
        if (.not. all(ieee_is_finite(xdata)) .or. .not. all(ieee_is_finite(ydata)) .or. &
                .not. all(ieee_is_finite(p0))) return

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(lb(n), ub(n))
        lb = -inf
        ub = inf
        if (present(lower)) then
            if (size(lower) /= n .or. any(ieee_is_nan(lower))) then
                res%message = 'lower bounds must match p0 and contain no NaNs'
                return
            end if
            lb = lower
        end if
        if (present(upper)) then
            if (size(upper) /= n .or. any(ieee_is_nan(upper))) then
                res%message = 'upper bounds must match p0 and contain no NaNs'
                return
            end if
            ub = upper
        end if
        if (any(lb >= ub)) then
            res%message = 'each lower bound must be strictly less than its upper bound'
            return
        end if
        if (any(p0 < lb) .or. any(p0 > ub)) then
            res%message = 'initial parameter vector is outside the supplied bounds'
            return
        end if
        has_bounds = any(ieee_is_finite(lb)) .or. any(ieee_is_finite(ub))
        meth = 'lm'
        if (has_bounds) meth = 'trf'
        if (present(method)) meth = lowercase(trim(method))
        if (trim(meth) /= 'lm' .and. trim(meth) /= 'trf' .and. &
                trim(meth) /= 'dogbox') then
            res%message = 'curve_fit method must be ''lm'', ''trf'', or ''dogbox'''
            return
        end if
        if (trim(meth) == 'lm' .and. has_bounds) then
            res%message = 'curve_fit method=''lm'' does not support bounds'
            return
        end if
        if (trim(meth) == 'lm' .and. m < n) then
            res%message = 'curve_fit method=''lm'' requires observations >= parameters'
            return
        end if
        if (trim(meth) /= 'lm' .and. (present(factor) .or. present(diag))) then
            res%message = 'curve_fit factor and diag options are only supported by method=''lm'''
            return
        end if
        diff_scheme_use = '2-point'
        if (present(diff_scheme)) diff_scheme_use = lowercase(trim(diff_scheme))
        if (trim(diff_scheme_use) /= '2-point' .and. trim(diff_scheme_use) /= '3-point' .and. &
                trim(diff_scheme_use) /= 'cs') then
            res%message = 'curve_fit diff_scheme must be ''2-point'', ''3-point'', or ''cs'''
            return
        end if
        if (trim(meth) == 'lm' .and. trim(diff_scheme_use) /= '2-point') then
            res%message = 'curve_fit method=''lm'' supports only default two-point differences'
            return
        end if
        if (trim(diff_scheme_use) == 'cs' .and. .not. present(jac) .and. &
                .not. present(complex_model)) then
            res%message = 'curve_fit diff_scheme=''cs'' requires complex_model when jac is absent'
            return
        end if
        loss_use = 'linear'
        if (present(loss)) loss_use = lowercase(trim(loss))
        select case (trim(loss_use))
        case ('linear', 'soft_l1', 'huber', 'cauchy', 'arctan')
        case default
            res%message = 'curve_fit loss must be linear, soft_l1, huber, cauchy, or arctan'
            return
        end select
        fscale_use = 1.0_dp
        if (present(f_scale)) fscale_use = f_scale
        if (.not. ieee_is_finite(fscale_use) .or. fscale_use <= 0.0_dp) then
            res%message = 'curve_fit f_scale must be finite and positive'
            return
        end if
        if (trim(meth) == 'lm' .and. trim(loss_use) /= 'linear') then
            res%message = 'curve_fit method=''lm'' supports only loss=''linear'''
            return
        end if

        allocate(weights(m), ywork(m))
        weights = 1.0_dp
        if (present(sigma)) then
            if (size(sigma) /= m) return
            if (.not. all(ieee_is_finite(sigma)) .or. any(sigma <= 0.0_dp)) return
            weights = 1.0_dp / sigma
        end if
        absolute_sigma_use = .false.
        if (present(absolute_sigma)) absolute_sigma_use = absolute_sigma
        model_ptr => model
        if (present(complex_model)) complex_model_ptr => complex_model
        if (present(jac)) jac_ptr => jac

        if (trim(meth) == 'lm') then
            if (present(jac)) then
                lsres = leastsq(residual_adapter, p0, m, jac=jacobian_adapter, ftol=ftol, &
                    xtol=xtol, gtol=gtol, maxfev=maxfev, factor=factor, diag=diag)
            else
                lsres = leastsq(residual_adapter, p0, m, ftol=ftol, xtol=xtol, &
                    gtol=gtol, maxfev=maxfev, factor=factor, diag=diag)
            end if
            res%popt = lsres%x
            res%fvec = lsres%fvec
            res%cost = 0.5_dp * dot_product(lsres%fvec, lsres%fvec)
            res%nfev = lsres%nfev
            res%njev = lsres%njev
            res%status = lsres%status
            res%success = lsres%success
            res%message = lsres%message
            if (allocated(lsres%cov_x)) then
                pcov_base = lsres%cov_x
                cov_ok = .true.
            else
                cov_ok = .false.
            end if
        else
            if (present(jac)) then
                trres = least_squares(residual_adapter, p0, m, jac=jacobian_adapter, &
                    method=meth, ftol=ftol, xtol=xtol, gtol=gtol, max_nfev=maxfev, &
                    lower=lb, upper=ub, loss=loss_use, f_scale=fscale_use, &
                    diff_step=diff_step, diff_scheme=diff_scheme_use)
            else if (trim(diff_scheme_use) == 'cs') then
                trres = least_squares(residual_adapter, p0, m, method=meth, ftol=ftol, &
                    xtol=xtol, gtol=gtol, max_nfev=maxfev, lower=lb, upper=ub, &
                    loss=loss_use, f_scale=fscale_use, diff_step=diff_step, &
                    diff_scheme=diff_scheme_use, complex_fun=complex_residual_adapter)
            else
                trres = least_squares(residual_adapter, p0, m, method=meth, ftol=ftol, &
                    xtol=xtol, gtol=gtol, max_nfev=maxfev, lower=lb, upper=ub, &
                    loss=loss_use, f_scale=fscale_use, diff_step=diff_step, &
                    diff_scheme=diff_scheme_use)
            end if
            if (.not. allocated(trres%x)) then
                res%message = trres%message
                return
            end if
            res%popt = trres%x
            res%fvec = trres%fun
            res%cost = trres%cost
            res%nfev = trres%nfev
            res%njev = trres%njev
            res%status = trres%status
            res%success = trres%success
            res%message = trres%message
            call covariance_from_jacobian(trres%jac, pcov_base, cov_ok)
        end if

        allocate(res%pcov(n, n))
        if (cov_ok) then
            res%pcov = pcov_base
            if (.not. absolute_sigma_use) then
                if (m > n) then
                    scale = 2.0_dp * res%cost / real(m - n, dp)
                    res%pcov = res%pcov * scale
                else
                    res%pcov = ieee_value(1.0_dp, ieee_positive_inf)
                end if
            end if
        else
            res%pcov = ieee_value(1.0_dp, ieee_positive_inf)
        end if

    contains

        subroutine residual_adapter(params, fvec)
            real(dp), intent(in) :: params(:) !! model parameters
            real(dp), intent(out) :: fvec(:) !! weighted model residuals

            call model_ptr(xdata, params, ywork)
            fvec = (ywork - ydata) * weights
        end subroutine residual_adapter

        subroutine complex_residual_adapter(params, fvec)
            complex(dp), intent(in) :: params(:) !! complex-step model parameters
            complex(dp), intent(out) :: fvec(:) !! weighted analytic residual continuation

            complex(dp) :: ymodel(m)

            call complex_model_ptr(xdata, params, ymodel)
            fvec = (ymodel - cmplx(ydata, 0.0_dp, kind=dp)) * weights
        end subroutine complex_residual_adapter

        subroutine jacobian_adapter(params, fjac)
            real(dp), intent(in) :: params(:) !! model parameters
            real(dp), intent(out) :: fjac(:, :) !! weighted model Jacobian

            integer :: i

            call jac_ptr(xdata, params, fjac)
            do i = 1, m
                fjac(i, :) = fjac(i, :) * weights(i)
            end do
        end subroutine jacobian_adapter

    end function curve_fit

    subroutine numerical_jacobian(fun, x, f0, fjac, diff_step, diff_scheme, complex_fun)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! differentiation point
        real(dp), intent(in) :: f0(:) !! residual vector at x
        real(dp), intent(out) :: fjac(:, :) !! numerical Jacobian
        real(dp), intent(in), optional :: diff_step !! user relative difference step
        character(len=*), intent(in) :: diff_scheme !! 2-point, 3-point, or cs scheme
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! complex residual for cs

        real(dp), allocatable :: f1(:), f2(:), x1(:), x2(:)
        real(dp) :: dx, h, rstep, signx
        integer :: j, n

        n = size(x)
        if (trim(diff_scheme) == 'cs') then
            if (present(complex_fun)) call complex_step_jacobian(complex_fun, x, fjac, diff_step)
            return
        end if
        allocate(f1(size(f0)), f2(size(f0)), x1(n), x2(n))
        rstep = finite_difference_rstep(diff_scheme)
        do j = 1, n
            signx = merge(1.0_dp, -1.0_dp, x(j) >= 0.0_dp)
            if (present(diff_step)) then
                h = diff_step * signx * abs(x(j))
                if (same_real((x(j) + h) - x(j), 0.0_dp)) then
                    h = rstep * signx * max(1.0_dp, abs(x(j)))
                end if
            else
                h = rstep * signx * max(1.0_dp, abs(x(j)))
            end if

            if (trim(diff_scheme) == '3-point') then
                h = abs(h)
                x1 = x
                x2 = x
                x1(j) = x(j) - h
                x2(j) = x(j) + h
                call fun(x1, f1)
                call fun(x2, f2)
                dx = x2(j) - x1(j)
                fjac(:, j) = (f2 - f1) / dx
            else
                x1 = x
                x1(j) = x(j) + h
                call fun(x1, f1)
                dx = x1(j) - x(j)
                fjac(:, j) = (f1 - f0) / dx
            end if
        end do
    end subroutine numerical_jacobian

    subroutine numerical_jacobian_bounded(fun, x, f0, fjac, diff_step, lb, ub, diff_scheme, &
            complex_fun, jac_sparsity)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! differentiation point inside bounds
        real(dp), intent(in) :: f0(:) !! residual vector at x
        real(dp), intent(out) :: fjac(:, :) !! bounds-aware numerical Jacobian
        real(dp), intent(in), optional :: diff_step !! user relative difference step
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        character(len=*), intent(in) :: diff_scheme !! 2-point, 3-point, or cs scheme
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! complex residual for cs
        logical, intent(in), optional :: jac_sparsity(:, :) !! structural nonzeros for grouped differences

        real(dp), allocatable :: f1(:), f2(:), x1(:), x2(:)
        real(dp) :: dx, h, lower_dist, rstep, signx, upper_dist
        integer :: j, n
        logical :: one_sided

        n = size(x)
        if (present(jac_sparsity)) then
            call numerical_jacobian_bounded_grouped(fun, x, f0, fjac, diff_step, lb, ub, &
                diff_scheme, jac_sparsity, complex_fun)
            return
        end if
        if (trim(diff_scheme) == 'cs') then
            if (present(complex_fun)) call complex_step_jacobian(complex_fun, x, fjac, diff_step)
            return
        end if
        allocate(f1(size(f0)), f2(size(f0)), x1(n), x2(n))
        rstep = finite_difference_rstep(diff_scheme)
        do j = 1, n
            signx = merge(1.0_dp, -1.0_dp, x(j) >= 0.0_dp)
            if (present(diff_step)) then
                h = diff_step * signx * abs(x(j))
                if (same_real((x(j) + h) - x(j), 0.0_dp)) then
                    h = rstep * signx * max(1.0_dp, abs(x(j)))
                end if
            else
                h = rstep * signx * max(1.0_dp, abs(x(j)))
            end if

            lower_dist = x(j) - lb(j)
            upper_dist = ub(j) - x(j)
            if (trim(diff_scheme) == '3-point') then
                h = abs(h)
                one_sided = .false.
                if (lower_dist < h .or. upper_dist < h) then
                    one_sided = .true.
                    if (upper_dist >= lower_dist) then
                        h = min(h, 0.5_dp * upper_dist)
                    else
                        h = -min(h, 0.5_dp * lower_dist)
                    end if
                    if (abs(h) <= min(lower_dist, upper_dist)) then
                        h = min(lower_dist, upper_dist)
                        one_sided = .false.
                    end if
                end if
                if (same_real(h, 0.0_dp)) then
                    fjac(:, j) = 0.0_dp
                    cycle
                end if
                x1 = x
                x2 = x
                if (one_sided) then
                    x1(j) = x(j) + h
                    x2(j) = x(j) + 2.0_dp * h
                    call fun(x1, f1)
                    call fun(x2, f2)
                    dx = x2(j) - x(j)
                    fjac(:, j) = (-3.0_dp * f0 + 4.0_dp * f1 - f2) / dx
                else
                    x1(j) = x(j) - h
                    x2(j) = x(j) + h
                    call fun(x1, f1)
                    call fun(x2, f2)
                    dx = x2(j) - x1(j)
                    fjac(:, j) = (f2 - f1) / dx
                end if
            else
                if (x(j) + h < lb(j) .or. x(j) + h > ub(j)) then
                    if (abs(h) <= max(lower_dist, upper_dist)) then
                        h = -h
                    else if (upper_dist >= lower_dist) then
                        h = upper_dist
                    else
                        h = -lower_dist
                    end if
                end if
                if (same_real(h, 0.0_dp)) then
                    fjac(:, j) = 0.0_dp
                    cycle
                end if
                x1 = x
                x1(j) = x(j) + h
                call fun(x1, f1)
                dx = x1(j) - x(j)
                fjac(:, j) = (f1 - f0) / dx
            end if
        end do
    end subroutine numerical_jacobian_bounded

    subroutine numerical_jacobian_bounded_grouped(fun, x, f0, fjac, diff_step, lb, ub, &
            diff_scheme, jac_sparsity, complex_fun)
        procedure(least_squares_residual_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! differentiation point inside bounds
        real(dp), intent(in) :: f0(:) !! residual vector at x
        real(dp), intent(out) :: fjac(:, :) !! grouped numerical Jacobian
        real(dp), intent(in), optional :: diff_step !! user relative difference step
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        character(len=*), intent(in) :: diff_scheme !! 2-point, 3-point, or cs scheme
        logical, intent(in) :: jac_sparsity(:, :) !! structural nonzero pattern
        procedure(least_squares_complex_residual_function), optional :: complex_fun !! complex residual for cs

        real(dp), allocatable :: dx(:), f1(:), f2(:), h(:), x1(:), x2(:)
        integer, allocatable :: groups(:)
        logical, allocatable :: one_sided(:)
        real(dp) :: lower_dist, rstep, signx, upper_dist
        integer :: g, j, n, ngroups

        n = size(x)
        allocate(groups(n))
        call group_sparsity_columns(jac_sparsity, groups, ngroups)
        fjac = 0.0_dp
        if (trim(diff_scheme) == 'cs') then
            if (present(complex_fun)) then
                call complex_step_jacobian_grouped(complex_fun, x, fjac, diff_step, &
                    jac_sparsity, groups, ngroups)
            end if
            return
        end if

        allocate(dx(n), f1(size(f0)), f2(size(f0)), h(n), one_sided(n), x1(n), x2(n))
        dx = 0.0_dp
        h = 0.0_dp
        one_sided = .false.
        rstep = finite_difference_rstep(diff_scheme)
        do j = 1, n
            signx = merge(1.0_dp, -1.0_dp, x(j) >= 0.0_dp)
            if (present(diff_step)) then
                h(j) = diff_step * signx * abs(x(j))
                if (same_real((x(j) + h(j)) - x(j), 0.0_dp)) then
                    h(j) = rstep * signx * max(1.0_dp, abs(x(j)))
                end if
            else
                h(j) = rstep * signx * max(1.0_dp, abs(x(j)))
            end if

            lower_dist = x(j) - lb(j)
            upper_dist = ub(j) - x(j)
            if (trim(diff_scheme) == '3-point') then
                h(j) = abs(h(j))
                if (lower_dist < h(j) .or. upper_dist < h(j)) then
                    one_sided(j) = .true.
                    if (upper_dist >= lower_dist) then
                        h(j) = min(h(j), 0.5_dp * upper_dist)
                    else
                        h(j) = -min(h(j), 0.5_dp * lower_dist)
                    end if
                    if (abs(h(j)) <= min(lower_dist, upper_dist)) then
                        h(j) = min(lower_dist, upper_dist)
                        one_sided(j) = .false.
                    end if
                end if
            else
                if (x(j) + h(j) < lb(j) .or. x(j) + h(j) > ub(j)) then
                    if (abs(h(j)) <= max(lower_dist, upper_dist)) then
                        h(j) = -h(j)
                    else if (upper_dist >= lower_dist) then
                        h(j) = upper_dist
                    else
                        h(j) = -lower_dist
                    end if
                end if
            end if
        end do

        do g = 1, ngroups
            x1 = x
            x2 = x
            do j = 1, n
                if (groups(j) /= g .or. same_real(h(j), 0.0_dp)) cycle
                if (trim(diff_scheme) == '3-point') then
                    if (one_sided(j)) then
                        x1(j) = x(j) + h(j)
                        x2(j) = x(j) + 2.0_dp * h(j)
                        dx(j) = x2(j) - x(j)
                    else
                        x1(j) = x(j) - h(j)
                        x2(j) = x(j) + h(j)
                        dx(j) = x2(j) - x1(j)
                    end if
                else
                    x1(j) = x(j) + h(j)
                    dx(j) = x1(j) - x(j)
                end if
            end do

            call fun(x1, f1)
            if (trim(diff_scheme) == '3-point') call fun(x2, f2)
            do j = 1, n
                if (groups(j) /= g .or. same_real(dx(j), 0.0_dp)) cycle
                if (trim(diff_scheme) == '3-point') then
                    if (one_sided(j)) then
                        where (jac_sparsity(:, j))
                            fjac(:, j) = (-3.0_dp * f0 + 4.0_dp * f1 - f2) / dx(j)
                        end where
                    else
                        where (jac_sparsity(:, j))
                            fjac(:, j) = (f2 - f1) / dx(j)
                        end where
                    end if
                else
                    where (jac_sparsity(:, j))
                        fjac(:, j) = (f1 - f0) / dx(j)
                    end where
                end if
            end do
        end do
    end subroutine numerical_jacobian_bounded_grouped

    subroutine complex_step_jacobian_grouped(fun, x, fjac, diff_step, jac_sparsity, &
            groups, ngroups)
        procedure(least_squares_complex_residual_function) :: fun !! analytic complex residual callback
        real(dp), intent(in) :: x(:) !! real differentiation point
        real(dp), intent(out) :: fjac(:, :) !! grouped complex-step Jacobian
        real(dp), intent(in), optional :: diff_step !! user relative difference step
        logical, intent(in) :: jac_sparsity(:, :) !! structural nonzero pattern
        integer, intent(in) :: groups(:) !! one-based finite-difference group per column
        integer, intent(in) :: ngroups !! number of finite-difference groups

        complex(dp), allocatable :: fc(:), xc(:)
        real(dp), allocatable :: h(:)
        real(dp) :: rstep, signx
        integer :: g, j, n

        n = size(x)
        allocate(fc(size(fjac, 1)), h(n), xc(n))
        rstep = sqrt(epsilon(1.0_dp))
        do j = 1, n
            signx = merge(1.0_dp, -1.0_dp, x(j) >= 0.0_dp)
            if (present(diff_step)) then
                h(j) = diff_step * signx * abs(x(j))
                if (same_real((x(j) + h(j)) - x(j), 0.0_dp)) then
                    h(j) = rstep * signx * max(1.0_dp, abs(x(j)))
                end if
            else
                h(j) = rstep * signx * max(1.0_dp, abs(x(j)))
            end if
        end do

        fjac = 0.0_dp
        do g = 1, ngroups
            xc = cmplx(x, 0.0_dp, kind=dp)
            do j = 1, n
                if (groups(j) == g) xc(j) = cmplx(x(j), h(j), kind=dp)
            end do
            call fun(xc, fc)
            do j = 1, n
                if (groups(j) /= g .or. same_real(h(j), 0.0_dp)) cycle
                where (jac_sparsity(:, j))
                    fjac(:, j) = aimag(fc) / h(j)
                end where
            end do
        end do
    end subroutine complex_step_jacobian_grouped

    subroutine group_sparsity_columns(jac_sparsity, groups, ngroups)
        logical, intent(in) :: jac_sparsity(:, :) !! structural nonzero pattern
        integer, intent(out) :: groups(:) !! one-based group assignment for each column
        integer, intent(out) :: ngroups !! number of nonconflicting column groups

        integer :: g, j, k, n
        logical :: conflict

        n = size(jac_sparsity, 2)
        groups = 0
        ngroups = 0
        do j = 1, n
            do g = 1, ngroups
                conflict = .false.
                do k = 1, j - 1
                    if (groups(k) /= g) cycle
                    if (any(jac_sparsity(:, j) .and. jac_sparsity(:, k))) then
                        conflict = .true.
                        exit
                    end if
                end do
                if (.not. conflict) then
                    groups(j) = g
                    exit
                end if
            end do
            if (groups(j) == 0) then
                ngroups = ngroups + 1
                groups(j) = ngroups
            end if
        end do
    end subroutine group_sparsity_columns

    subroutine complex_step_jacobian(fun, x, fjac, diff_step)
        procedure(least_squares_complex_residual_function) :: fun !! analytic complex residual callback
        real(dp), intent(in) :: x(:) !! real differentiation point
        real(dp), intent(out) :: fjac(:, :) !! complex-step numerical Jacobian
        real(dp), intent(in), optional :: diff_step !! user relative difference step

        complex(dp), allocatable :: fc(:), xc(:)
        real(dp) :: h, rstep, signx
        integer :: j, n

        n = size(x)
        allocate(fc(size(fjac, 1)), xc(n))
        rstep = sqrt(epsilon(1.0_dp))
        do j = 1, n
            signx = merge(1.0_dp, -1.0_dp, x(j) >= 0.0_dp)
            if (present(diff_step)) then
                h = diff_step * signx * abs(x(j))
                if (same_real((x(j) + h) - x(j), 0.0_dp)) then
                    h = rstep * signx * max(1.0_dp, abs(x(j)))
                end if
            else
                h = rstep * signx * max(1.0_dp, abs(x(j)))
            end if
            xc = cmplx(x, 0.0_dp, kind=dp)
            xc(j) = xc(j) + cmplx(0.0_dp, h, kind=dp)
            call fun(xc, fc)
            fjac(:, j) = aimag(fc) / h
        end do
    end subroutine complex_step_jacobian

    pure real(dp) function finite_difference_rstep(diff_scheme) result(rstep)
        character(len=*), intent(in) :: diff_scheme !! 2-point, 3-point, or cs scheme

        if (trim(diff_scheme) == '3-point') then
            rstep = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
        else
            rstep = sqrt(epsilon(1.0_dp))
        end if
    end function finite_difference_rstep

    pure logical function valid_trf_tolerances(ftol, xtol, gtol) result(ok)
        real(dp), intent(in) :: ftol !! relative cost tolerance
        real(dp), intent(in) :: xtol !! relative step tolerance
        real(dp), intent(in) :: gtol !! first-order tolerance

        ok = ieee_is_finite(ftol) .and. ieee_is_finite(xtol) .and. &
            ieee_is_finite(gtol) .and. ftol >= 0.0_dp .and. xtol >= 0.0_dp .and. &
            gtol >= 0.0_dp .and. max(ftol, xtol, gtol) > epsilon(1.0_dp)
    end function valid_trf_tolerances

    pure real(dp) function vector_norm2(x) result(value)
        real(dp), intent(in) :: x(:) !! vector whose Euclidean norm is required

        value = sqrt(max(0.0_dp, dot_product(x, x)))
    end function vector_norm2

    subroutine cl_scaling_vector(x, g, lb, ub, v, dv)
        real(dp), intent(in) :: x(:) !! current feasible point
        real(dp), intent(in) :: g(:) !! least-squares gradient
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        real(dp), intent(out) :: v(:) !! Coleman-Li scaling vector
        real(dp), intent(out) :: dv(:) !! diagonal derivative of the scaling vector

        integer :: i

        v = 1.0_dp
        dv = 0.0_dp
        do i = 1, size(x)
            if (g(i) < 0.0_dp .and. ieee_is_finite(ub(i))) then
                v(i) = ub(i) - x(i)
                dv(i) = -1.0_dp
            else if (g(i) > 0.0_dp .and. ieee_is_finite(lb(i))) then
                v(i) = x(i) - lb(i)
                dv(i) = 1.0_dp
            end if
        end do
    end subroutine cl_scaling_vector

    pure logical function point_in_bounds(x, lb, ub) result(ok)
        real(dp), intent(in) :: x(:) !! point to test
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds

        ok = all(x >= lb) .and. all(x <= ub)
    end function point_in_bounds

    subroutine find_active_constraints(x, lb, ub, rtol, active)
        real(dp), intent(in) :: x(:) !! feasible point
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        real(dp), intent(in) :: rtol !! relative activity tolerance
        integer, intent(out) :: active(:) !! -1 lower, 0 inactive, +1 upper

        real(dp) :: lower_dist, lower_threshold, upper_dist, upper_threshold
        integer :: i

        active = 0
        do i = 1, size(x)
            if (same_real(rtol, 0.0_dp)) then
                if (x(i) <= lb(i)) active(i) = -1
                if (x(i) >= ub(i)) active(i) = 1
                cycle
            end if
            lower_dist = x(i) - lb(i)
            upper_dist = ub(i) - x(i)
            lower_threshold = rtol * max(1.0_dp, abs(lb(i)))
            upper_threshold = rtol * max(1.0_dp, abs(ub(i)))
            if (ieee_is_finite(lb(i)) .and. &
                    lower_dist <= min(upper_dist, lower_threshold)) then
                active(i) = -1
            end if
            if (ieee_is_finite(ub(i)) .and. &
                    upper_dist <= min(lower_dist, upper_threshold)) then
                active(i) = 1
            end if
        end do
    end subroutine find_active_constraints

    subroutine make_strictly_feasible(x, lb, ub, rstep, x_new)
        real(dp), intent(in) :: x(:) !! feasible or boundary point
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        real(dp), intent(in) :: rstep !! relative step from active bounds; zero uses nearest
        real(dp), intent(out) :: x_new(:) !! strictly feasible point

        integer, allocatable :: active(:)
        integer :: i, n

        n = size(x)
        allocate(active(n))
        x_new = x
        call find_active_constraints(x, lb, ub, rstep, active)
        do i = 1, n
            if (active(i) == -1) then
                if (same_real(rstep, 0.0_dp)) then
                    x_new(i) = nearest(lb(i), 1.0_dp)
                else
                    x_new(i) = lb(i) + rstep * max(1.0_dp, abs(lb(i)))
                end if
            else if (active(i) == 1) then
                if (same_real(rstep, 0.0_dp)) then
                    x_new(i) = nearest(ub(i), -1.0_dp)
                else
                    x_new(i) = ub(i) - rstep * max(1.0_dp, abs(ub(i)))
                end if
            end if
            if (x_new(i) < lb(i) .or. x_new(i) > ub(i)) then
                x_new(i) = 0.5_dp * (lb(i) + ub(i))
            end if
        end do
    end subroutine make_strictly_feasible

    subroutine step_size_to_bound(x, s, lb, ub, step, hits)
        real(dp), intent(in) :: x(:) !! feasible starting point
        real(dp), intent(in) :: s(:) !! search direction
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        real(dp), intent(out) :: step !! first nonnegative step reaching a bound
        integer, intent(out) :: hits(:) !! -1 lower, 0 none, +1 upper at first step

        real(dp) :: inf, ti
        integer :: i

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        step = inf
        hits = 0
        do i = 1, size(x)
            if (s(i) > 0.0_dp) then
                ti = (ub(i) - x(i)) / s(i)
            else if (s(i) < 0.0_dp) then
                ti = (lb(i) - x(i)) / s(i)
            else
                cycle
            end if
            if (ti < step) step = ti
        end do
        do i = 1, size(x)
            if (s(i) > 0.0_dp) then
                ti = (ub(i) - x(i)) / s(i)
                if (same_real(ti, step)) hits(i) = 1
            else if (s(i) < 0.0_dp) then
                ti = (lb(i) - x(i)) / s(i)
                if (same_real(ti, step)) hits(i) = -1
            end if
        end do
    end subroutine step_size_to_bound

    subroutine intersect_trust_region(x, s, delta, t_neg, t_pos, ok)
        real(dp), intent(in) :: x(:) !! point inside the spherical trust region
        real(dp), intent(in) :: s(:) !! nonzero line direction
        real(dp), intent(in) :: delta !! trust-region radius
        real(dp), intent(out) :: t_neg !! negative intersection parameter
        real(dp), intent(out) :: t_pos !! positive intersection parameter
        logical, intent(out) :: ok !! true if a valid intersection was found

        real(dp) :: a, b, c, d, q, t1, t2

        a = dot_product(s, s)
        c = dot_product(x, x) - delta * delta
        if (same_real(a, 0.0_dp) .or. c > 0.0_dp) then
            ok = .false.
            t_neg = 0.0_dp
            t_pos = 0.0_dp
            return
        end if
        b = dot_product(x, s)
        d = sqrt(max(0.0_dp, b * b - a * c))
        q = -(b + sign(d, b))
        t1 = q / a
        if (same_real(q, 0.0_dp)) then
            t2 = -b / a
        else
            t2 = c / q
        end if
        t_neg = min(t1, t2)
        t_pos = max(t1, t2)
        ok = .true.
    end subroutine intersect_trust_region

    subroutine build_quadratic_1d_operator(matvec, m, g, s, a, b, diag, s0, c)
        procedure(least_squares_linear_operator_function) :: matvec !! compute J*s
        integer, intent(in) :: m !! residual-space dimension
        real(dp), intent(in) :: g(:) !! quadratic linear coefficient
        real(dp), intent(in) :: s(:) !! line direction
        real(dp), intent(out) :: a !! coefficient of t squared
        real(dp), intent(out) :: b !! coefficient of t
        real(dp), intent(in), optional :: diag(:) !! optional nonnegative diagonal term
        real(dp), intent(in), optional :: s0(:) !! optional line origin
        real(dp), intent(out), optional :: c !! constant term when s0 is present

        real(dp), allocatable :: u(:), v(:)

        allocate(v(m))
        call matvec(s, v)
        a = dot_product(v, v)
        if (present(diag)) a = a + dot_product(s * diag, s)
        a = 0.5_dp * a
        b = dot_product(g, s)
        if (present(s0)) then
            allocate(u(m))
            call matvec(s0, u)
            b = b + dot_product(u, v)
            if (present(diag)) b = b + dot_product(s0 * diag, s)
            if (present(c)) then
                c = 0.5_dp * dot_product(u, u) + dot_product(g, s0)
                if (present(diag)) c = c + 0.5_dp * dot_product(s0 * diag, s0)
            end if
        else if (present(c)) then
            c = 0.0_dp
        end if
    end subroutine build_quadratic_1d_operator

    real(dp) function evaluate_quadratic_operator(matvec, m, g, s, diag) result(value)
        procedure(least_squares_linear_operator_function) :: matvec !! compute J*s
        integer, intent(in) :: m !! residual-space dimension
        real(dp), intent(in) :: g(:) !! quadratic linear coefficient
        real(dp), intent(in) :: s(:) !! trial step
        real(dp), intent(in), optional :: diag(:) !! optional diagonal quadratic term

        real(dp), allocatable :: js(:)
        real(dp) :: quadratic

        allocate(js(m))
        call matvec(s, js)
        quadratic = dot_product(js, js)
        if (present(diag)) quadratic = quadratic + dot_product(s * diag, s)
        value = 0.5_dp * quadratic + dot_product(g, s)
    end function evaluate_quadratic_operator

    subroutine build_quadratic_1d(jmat, g, s, a, b, diag, s0, c)
        real(dp), intent(in) :: jmat(:, :) !! dense least-squares Jacobian
        real(dp), intent(in) :: g(:) !! quadratic linear coefficient
        real(dp), intent(in) :: s(:) !! line direction
        real(dp), intent(out) :: a !! coefficient of t squared
        real(dp), intent(out) :: b !! coefficient of t
        real(dp), intent(in), optional :: diag(:) !! optional nonnegative diagonal term
        real(dp), intent(in), optional :: s0(:) !! optional line origin
        real(dp), intent(out), optional :: c !! constant term when s0 is present

        real(dp), allocatable :: u(:), v(:)

        allocate(v(size(jmat, 1)))
        v = matmul(jmat, s)
        a = dot_product(v, v)
        if (present(diag)) a = a + dot_product(s * diag, s)
        a = 0.5_dp * a
        b = dot_product(g, s)
        if (present(s0)) then
            allocate(u(size(jmat, 1)))
            u = matmul(jmat, s0)
            b = b + dot_product(u, v)
            if (present(diag)) b = b + dot_product(s0 * diag, s)
            if (present(c)) then
                c = 0.5_dp * dot_product(u, u) + dot_product(g, s0)
                if (present(diag)) c = c + 0.5_dp * dot_product(s0 * diag, s0)
            end if
        else if (present(c)) then
            c = 0.0_dp
        end if
    end subroutine build_quadratic_1d

    subroutine minimize_quadratic_1d(a, b, lb, ub, t, value, c)
        real(dp), intent(in) :: a !! coefficient of t squared
        real(dp), intent(in) :: b !! coefficient of t
        real(dp), intent(in) :: lb !! finite lower bound on t
        real(dp), intent(in) :: ub !! finite upper bound on t
        real(dp), intent(out) :: t !! minimizing parameter
        real(dp), intent(out) :: value !! quadratic value at t
        real(dp), intent(in), optional :: c !! constant term; default zero

        real(dp) :: cc, extremum, value_lb, value_ub, value_x

        cc = 0.0_dp
        if (present(c)) cc = c
        t = lb
        value_lb = lb * (a * lb + b) + cc
        value = value_lb
        value_ub = ub * (a * ub + b) + cc
        if (value_ub < value) then
            t = ub
            value = value_ub
        end if
        if (.not. same_real(a, 0.0_dp)) then
            extremum = -0.5_dp * b / a
            if (extremum > lb .and. extremum < ub) then
                value_x = extremum * (a * extremum + b) + cc
                if (value_x < value) then
                    t = extremum
                    value = value_x
                end if
            end if
        end if
    end subroutine minimize_quadratic_1d

    real(dp) function evaluate_quadratic(jmat, g, s, diag) result(value)
        real(dp), intent(in) :: jmat(:, :) !! dense least-squares Jacobian
        real(dp), intent(in) :: g(:) !! quadratic linear coefficient
        real(dp), intent(in) :: s(:) !! trial step
        real(dp), intent(in), optional :: diag(:) !! optional diagonal quadratic term

        real(dp), allocatable :: js(:)
        real(dp) :: quadratic

        allocate(js(size(jmat, 1)))
        js = matmul(jmat, s)
        quadratic = dot_product(js, js)
        if (present(diag)) quadratic = quadratic + dot_product(s * diag, s)
        value = 0.5_dp * quadratic + dot_product(g, s)
    end function evaluate_quadratic

    subroutine update_tr_radius(delta, actual_reduction, predicted_reduction, &
            step_norm, bound_hit, delta_new, ratio)
        real(dp), intent(in) :: delta !! current trust-region radius
        real(dp), intent(in) :: actual_reduction !! observed cost reduction
        real(dp), intent(in) :: predicted_reduction !! model-predicted reduction
        real(dp), intent(in) :: step_norm !! trust-region norm of the accepted trial step
        logical, intent(in) :: bound_hit !! true when the trust-region boundary was hit
        real(dp), intent(out) :: delta_new !! updated trust-region radius
        real(dp), intent(out) :: ratio !! actual to predicted reduction ratio

        if (predicted_reduction > 0.0_dp) then
            ratio = actual_reduction / predicted_reduction
        else if (same_real(predicted_reduction, 0.0_dp) .and. &
                same_real(actual_reduction, 0.0_dp)) then
            ratio = 1.0_dp
        else
            ratio = 0.0_dp
        end if
        delta_new = delta
        if (ratio < 0.25_dp) then
            delta_new = 0.25_dp * step_norm
        else if (ratio > 0.75_dp .and. bound_hit) then
            delta_new = 2.0_dp * delta
        end if
    end subroutine update_tr_radius

    subroutine check_ls_termination(df, cost, step_norm, x_norm, ratio, ftol, xtol, status)
        real(dp), intent(in) :: df !! actual cost reduction
        real(dp), intent(in) :: cost !! previous cost
        real(dp), intent(in) :: step_norm !! Euclidean step norm
        real(dp), intent(in) :: x_norm !! Euclidean parameter norm
        real(dp), intent(in) :: ratio !! trust-region agreement ratio
        real(dp), intent(in) :: ftol !! relative cost tolerance
        real(dp), intent(in) :: xtol !! relative step tolerance
        integer, intent(out) :: status !! 2, 3, 4 on convergence; -99 otherwise

        logical :: ftol_satisfied, xtol_satisfied

        ftol_satisfied = df < ftol * cost .and. ratio > 0.25_dp
        xtol_satisfied = step_norm < xtol * (xtol + x_norm)
        if (ftol_satisfied .and. xtol_satisfied) then
            status = 4
        else if (ftol_satisfied) then
            status = 2
        else if (xtol_satisfied) then
            status = 3
        else
            status = -99
        end if
    end subroutine check_ls_termination

    subroutine build_lsmr_subspace_operator(matvec, rmatvec, m, n, diag, g, fvec, &
            delta, regularize, basis, ksub, bsub, gsub, ok)
        procedure(least_squares_linear_operator_function) :: matvec !! scaled Jacobian product
        procedure(least_squares_linear_operator_function) :: rmatvec !! scaled transpose product
        integer, intent(in) :: m !! residual count
        integer, intent(in) :: n !! parameter count
        real(dp), intent(in) :: diag(:) !! nonnegative diagonal Coleman-Li term
        real(dp), intent(in) :: g(:) !! gradient in scaled variables
        real(dp), intent(in) :: fvec(:) !! scaled residual vector
        real(dp), intent(in) :: delta !! current trust-region radius
        logical, intent(in) :: regularize !! add SciPy-style regularization term
        real(dp), intent(out) :: basis(:, :) !! orthonormal one/two-dimensional basis
        integer, intent(out) :: ksub !! active basis dimension
        real(dp), intent(out) :: bsub(:, :) !! reduced quadratic Hessian
        real(dp), intent(out) :: gsub(:) !! reduced gradient
        logical, intent(out) :: ok !! true when a usable subspace was built

        real(dp), allocatable :: baug(:), gn(:), js(:, :), second(:), sqrt_diag(:)
        real(dp), allocatable :: work_m(:), work_n(:)
        real(dp) :: a, ag_value, b, gnorm, reg_term, second_norm, threshold, to_tr, tmin
        integer :: i, j

        basis = 0.0_dp
        bsub = 0.0_dp
        gsub = 0.0_dp
        ksub = 0
        ok = .false.
        if (size(g) /= n .or. size(diag) /= n .or. size(fvec) /= m) return
        if (size(basis, 1) /= n .or. size(basis, 2) < 2) return
        if (size(bsub, 1) < 2 .or. size(bsub, 2) < 2 .or. size(gsub) < 2) return

        gnorm = vector_norm2(g)
        if (same_real(gnorm, 0.0_dp) .or. delta <= 0.0_dp) return

        reg_term = 0.0_dp
        if (regularize) then
            call build_quadratic_1d_operator(matvec, m, g, -g, a, b, diag)
            to_tr = delta / gnorm
            call minimize_quadratic_1d(a, b, 0.0_dp, to_tr, tmin, ag_value)
            reg_term = max(0.0_dp, -ag_value / (delta * delta))
        end if

        allocate(sqrt_diag(n), baug(m + n), gn(n), second(n), js(m, 2))
        allocate(work_m(m), work_n(n))
        sqrt_diag = sqrt(max(0.0_dp, diag + reg_term))
        baug(1:m) = fvec
        baug(m + 1:m + n) = 0.0_dp
        call lsmr_solve(aug_matvec, aug_rmatvec, m + n, n, baug, gn)
        if (.not. all(ieee_is_finite(gn))) return

        basis(:, 1) = g / gnorm
        second = gn - dot_product(basis(:, 1), gn) * basis(:, 1)
        second_norm = vector_norm2(second)
        threshold = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, vector_norm2(gn))
        if (second_norm > threshold) then
            basis(:, 2) = second / second_norm
            ksub = 2
        else
            ksub = 1
        end if

        do i = 1, ksub
            call matvec(basis(:, i), js(:, i))
            gsub(i) = dot_product(basis(:, i), g)
            do j = 1, ksub
                bsub(i, j) = dot_product(js(:, i), js(:, j)) + &
                    dot_product(basis(:, i) * diag, basis(:, j))
            end do
        end do
        ok = all(ieee_is_finite(js(:, 1:ksub))) .and. all(ieee_is_finite(bsub(1:ksub, 1:ksub)))

    contains

        subroutine aug_matvec(xin, yout)
            real(dp), intent(in) :: xin(:) !! subproblem vector
            real(dp), intent(out) :: yout(:) !! augmented Jacobian product

            call matvec(xin, work_m)
            yout(1:m) = work_m
            yout(m + 1:m + n) = sqrt_diag * xin
        end subroutine aug_matvec

        subroutine aug_rmatvec(yin, xout)
            real(dp), intent(in) :: yin(:) !! augmented residual-space vector
            real(dp), intent(out) :: xout(:) !! transposed augmented Jacobian product

            call rmatvec(yin(1:m), work_n)
            xout = work_n + sqrt_diag * yin(m + 1:m + n)
        end subroutine aug_rmatvec

    end subroutine build_lsmr_subspace_operator

    subroutine build_lsmr_subspace(jmat, diag, g, fvec, delta, regularize, &
            basis, ksub, bsub, gsub, ok)
        real(dp), intent(in) :: jmat(:, :) !! Jacobian in trust-region scaled variables
        real(dp), intent(in) :: diag(:) !! nonnegative diagonal Coleman-Li term
        real(dp), intent(in) :: g(:) !! gradient in scaled variables
        real(dp), intent(in) :: fvec(:) !! scaled residual vector
        real(dp), intent(in) :: delta !! current trust-region radius
        logical, intent(in) :: regularize !! add SciPy-style regularization term
        real(dp), intent(out) :: basis(:, :) !! orthonormal one/two-dimensional basis
        integer, intent(out) :: ksub !! active basis dimension
        real(dp), intent(out) :: bsub(:, :) !! reduced quadratic Hessian
        real(dp), intent(out) :: gsub(:) !! reduced gradient
        logical, intent(out) :: ok !! true when a usable subspace was built

        real(dp), allocatable :: baug(:), gn(:), js(:, :), second(:), sqrt_diag(:)
        real(dp) :: a, ag_value, b, gnorm, reg_term, second_norm, threshold, to_tr, tmin
        integer :: i, j, m, n

        m = size(jmat, 1)
        n = size(jmat, 2)
        basis = 0.0_dp
        bsub = 0.0_dp
        gsub = 0.0_dp
        ksub = 0
        ok = .false.
        if (size(g) /= n .or. size(diag) /= n .or. size(fvec) /= m) return
        if (size(basis, 1) /= n .or. size(basis, 2) < 2) return
        if (size(bsub, 1) < 2 .or. size(bsub, 2) < 2 .or. size(gsub) < 2) return

        gnorm = vector_norm2(g)
        if (same_real(gnorm, 0.0_dp) .or. delta <= 0.0_dp) return

        reg_term = 0.0_dp
        if (regularize) then
            call build_quadratic_1d(jmat, g, -g, a, b, diag)
            to_tr = delta / gnorm
            call minimize_quadratic_1d(a, b, 0.0_dp, to_tr, tmin, ag_value)
            reg_term = max(0.0_dp, -ag_value / (delta * delta))
        end if

        allocate(sqrt_diag(n), baug(m + n), gn(n), second(n), js(m, 2))
        sqrt_diag = sqrt(max(0.0_dp, diag + reg_term))
        baug(1:m) = fvec
        baug(m + 1:m + n) = 0.0_dp
        call lsmr_solve(aug_matvec, aug_rmatvec, m + n, n, baug, gn)
        if (.not. all(ieee_is_finite(gn))) return

        basis(:, 1) = g / gnorm
        second = gn - dot_product(basis(:, 1), gn) * basis(:, 1)
        second_norm = vector_norm2(second)
        threshold = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, vector_norm2(gn))
        if (second_norm > threshold) then
            basis(:, 2) = second / second_norm
            ksub = 2
        else
            ksub = 1
        end if

        js(:, 1:ksub) = matmul(jmat, basis(:, 1:ksub))
        do i = 1, ksub
            gsub(i) = dot_product(basis(:, i), g)
            do j = 1, ksub
                bsub(i, j) = dot_product(js(:, i), js(:, j)) + &
                    dot_product(basis(:, i) * diag, basis(:, j))
            end do
        end do
        ok = .true.

    contains

        subroutine aug_matvec(x, y)
            real(dp), intent(in) :: x(:) !! subproblem vector
            real(dp), intent(out) :: y(:) !! augmented Jacobian product

            y(1:m) = matmul(jmat, x)
            y(m + 1:m + n) = sqrt_diag * x
        end subroutine aug_matvec

        subroutine aug_rmatvec(y, x)
            real(dp), intent(in) :: y(:) !! augmented residual-space vector
            real(dp), intent(out) :: x(:) !! transposed augmented Jacobian product

            x = matmul(transpose(jmat), y(1:m)) + sqrt_diag * y(m + 1:m + n)
        end subroutine aug_rmatvec

    end subroutine build_lsmr_subspace

    subroutine solve_trust_region_dense(bmat, g, delta, alpha_initial, p, alpha, ok)
        real(dp), intent(in) :: bmat(:, :) !! symmetric positive-semidefinite model Hessian
        real(dp), intent(in) :: g(:) !! model gradient
        real(dp), intent(in) :: delta !! positive trust-region radius
        real(dp), intent(in) :: alpha_initial !! previous regularization parameter
        real(dp), intent(out) :: p(:) !! trust-region minimizer
        real(dp), intent(out) :: alpha !! final regularization parameter
        logical, intent(out) :: ok !! true on a converged eigensolve/root solve

        real(dp), allocatable :: coeff(:), eig(:), vec(:, :)
        real(dp) :: alpha_hi, alpha_lo, alpha_new, cutoff, deriv, pnorm, scale
        integer :: i, iter, n, status
        logical :: has_hard_component

        n = size(g)
        allocate(coeff(n), eig(n), vec(n, n))
        call symmetric_eigen_jacobi(bmat, eig, vec, status)
        if (status /= linalg_status_success) then
            ok = .false.
            p = 0.0_dp
            alpha = 0.0_dp
            return
        end if
        coeff = matmul(transpose(vec), -g)
        scale = max(1.0_dp, maxval(abs(eig)))
        cutoff = 64.0_dp * epsilon(1.0_dp) * scale
        p = 0.0_dp
        has_hard_component = .false.
        do i = 1, n
            if (eig(i) > cutoff) then
                p = p + vec(:, i) * (coeff(i) / eig(i))
            else if (abs(coeff(i)) > cutoff) then
                has_hard_component = .true.
            end if
        end do
        if (.not. has_hard_component .and. vector_norm2(p) <= delta .and. &
                minval(eig) >= -cutoff) then
            alpha = 0.0_dp
            ok = .true.
            return
        end if

        alpha_lo = max(0.0_dp, -minval(eig) + cutoff)
        alpha_hi = max(1.0_dp, maxval(abs(eig)), vector_norm2(g) / max(delta, tiny(1.0_dp)))
        do
            call trust_step_at_alpha(eig, vec, coeff, alpha_hi, p)
            if (vector_norm2(p) <= delta) exit
            alpha_hi = 2.0_dp * alpha_hi
            if (.not. ieee_is_finite(alpha_hi)) then
                ok = .false.
                alpha = alpha_hi
                return
            end if
        end do
        alpha = min(alpha_hi, max(alpha_lo, alpha_initial))
        if (alpha <= alpha_lo) alpha = 0.5_dp * (alpha_lo + alpha_hi)

        do iter = 1, 80
            call trust_step_at_alpha(eig, vec, coeff, alpha, p)
            pnorm = vector_norm2(p)
            if (abs(pnorm - delta) <= 0.005_dp * delta) then
                ok = .true.
                return
            end if
            if (pnorm > delta) then
                alpha_lo = alpha
            else
                alpha_hi = alpha
            end if
            deriv = trust_norm_derivative(eig, coeff, alpha, pnorm)
            if (deriv < 0.0_dp) then
                alpha_new = alpha - (pnorm - delta) / deriv
            else
                alpha_new = 0.5_dp * (alpha_lo + alpha_hi)
            end if
            if (alpha_new <= alpha_lo .or. alpha_new >= alpha_hi .or. &
                    .not. ieee_is_finite(alpha_new)) then
                alpha_new = 0.5_dp * (alpha_lo + alpha_hi)
            end if
            alpha = alpha_new
        end do
        ok = .true.
    end subroutine solve_trust_region_dense

    subroutine trust_step_at_alpha(eig, vec, coeff, alpha, p)
        real(dp), intent(in) :: eig(:) !! model-Hessian eigenvalues
        real(dp), intent(in) :: vec(:, :) !! model-Hessian eigenvectors
        real(dp), intent(in) :: coeff(:) !! transformed negative gradient
        real(dp), intent(in) :: alpha !! nonnegative regularization parameter
        real(dp), intent(out) :: p(:) !! regularized Newton step

        real(dp) :: denom
        integer :: i

        p = 0.0_dp
        do i = 1, size(eig)
            denom = eig(i) + alpha
            if (denom <= tiny(1.0_dp)) denom = tiny(1.0_dp)
            p = p + vec(:, i) * (coeff(i) / denom)
        end do
    end subroutine trust_step_at_alpha

    pure real(dp) function trust_norm_derivative(eig, coeff, alpha, pnorm) result(value)
        real(dp), intent(in) :: eig(:) !! model-Hessian eigenvalues
        real(dp), intent(in) :: coeff(:) !! transformed negative gradient
        real(dp), intent(in) :: alpha !! regularization parameter
        real(dp), intent(in) :: pnorm !! current step norm

        real(dp) :: denom, total
        integer :: i

        if (pnorm <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        total = 0.0_dp
        do i = 1, size(eig)
            denom = eig(i) + alpha
            if (denom <= tiny(1.0_dp)) then
                value = 0.0_dp
                return
            end if
            total = total + coeff(i) * coeff(i) / (denom * denom * denom)
        end do
        value = -total / pnorm
    end function trust_norm_derivative

    subroutine select_trf_step_operator(x, matvec, m, diag_h, g_h, p_in, p_h_in, &
            d, delta, lb, ub, theta, step, step_h, predicted_reduction)
        real(dp), intent(in) :: x(:) !! current feasible point
        procedure(least_squares_linear_operator_function) :: matvec !! scaled Jacobian product
        integer, intent(in) :: m !! residual-space dimension
        real(dp), intent(in) :: diag_h(:) !! nonnegative Coleman-Li diagonal term
        real(dp), intent(in) :: g_h(:) !! gradient in scaled variables
        real(dp), intent(in) :: p_in(:) !! trust-region step in original variables
        real(dp), intent(in) :: p_h_in(:) !! trust-region step in scaled variables
        real(dp), intent(in) :: d(:) !! original-to-scaled step multipliers
        real(dp), intent(in) :: delta !! trust-region radius
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        real(dp), intent(in) :: theta !! step-back factor from bounds
        real(dp), intent(out) :: step(:) !! selected original-space step
        real(dp), intent(out) :: step_h(:) !! selected scaled-space step
        real(dp), intent(out) :: predicted_reduction !! positive predicted cost reduction

        integer, allocatable :: hits(:)
        real(dp), allocatable :: ag(:), ag_h(:), p(:), p_h(:), r(:), r_h(:), x_on_bound(:)
        real(dp) :: a, ag_limit, ag_stride, ag_value, b, c, p_stride, p_value
        real(dp) :: r_stride, r_stride_l, r_stride_u, r_value, t_neg, to_bound, to_tr
        logical :: ok

        allocate(hits(size(x)), ag(size(x)), ag_h(size(x)), p(size(x)), p_h(size(x)))
        allocate(r(size(x)), r_h(size(x)), x_on_bound(size(x)))
        p = p_in
        p_h = p_h_in
        if (point_in_bounds(x + p, lb, ub)) then
            p_value = evaluate_quadratic_operator(matvec, m, g_h, p_h, diag_h)
            step = p
            step_h = p_h
            predicted_reduction = -p_value
            return
        end if

        call step_size_to_bound(x, p, lb, ub, p_stride, hits)
        r_h = p_h
        where (hits /= 0) r_h = -r_h
        r = d * r_h
        p = p * p_stride
        p_h = p_h * p_stride
        x_on_bound = x + p

        call intersect_trust_region(p_h, r_h, delta, t_neg, to_tr, ok)
        if (.not. ok) to_tr = 0.0_dp
        call step_size_to_bound(x_on_bound, r, lb, ub, to_bound, hits)
        r_stride = min(to_bound, to_tr)
        if (r_stride > 0.0_dp) then
            r_stride_l = (1.0_dp - theta) * p_stride / r_stride
            if (same_real(r_stride, to_bound)) then
                r_stride_u = theta * to_bound
            else
                r_stride_u = to_tr
            end if
        else
            r_stride_l = 0.0_dp
            r_stride_u = -1.0_dp
        end if

        if (r_stride_l <= r_stride_u) then
            call build_quadratic_1d_operator(matvec, m, g_h, r_h, a, b, diag_h, p_h, c)
            call minimize_quadratic_1d(a, b, r_stride_l, r_stride_u, r_stride, r_value, c)
            r_h = p_h + r_stride * r_h
            r = d * r_h
        else
            r_value = ieee_value(1.0_dp, ieee_positive_inf)
        end if

        p = theta * p
        p_h = theta * p_h
        p_value = evaluate_quadratic_operator(matvec, m, g_h, p_h, diag_h)

        ag_h = -g_h
        ag = d * ag_h
        if (vector_norm2(ag_h) > 0.0_dp) then
            to_tr = delta / vector_norm2(ag_h)
            call step_size_to_bound(x, ag, lb, ub, to_bound, hits)
            if (to_bound < to_tr) then
                ag_limit = theta * to_bound
            else
                ag_limit = to_tr
            end if
            call build_quadratic_1d_operator(matvec, m, g_h, ag_h, a, b, diag_h)
            call minimize_quadratic_1d(a, b, 0.0_dp, ag_limit, ag_stride, ag_value)
            ag_h = ag_stride * ag_h
            ag = ag_stride * ag
        else
            ag_value = ieee_value(1.0_dp, ieee_positive_inf)
        end if

        if (p_value < r_value .and. p_value < ag_value) then
            step = p
            step_h = p_h
            predicted_reduction = -p_value
        else if (r_value < p_value .and. r_value < ag_value) then
            step = r
            step_h = r_h
            predicted_reduction = -r_value
        else
            step = ag
            step_h = ag_h
            predicted_reduction = -ag_value
        end if
    end subroutine select_trf_step_operator

    subroutine select_trf_step(x, j_h, diag_h, g_h, p_in, p_h_in, d, delta, lb, ub, &
            theta, step, step_h, predicted_reduction)
        real(dp), intent(in) :: x(:) !! current feasible point
        real(dp), intent(in) :: j_h(:, :) !! Jacobian in scaled variables
        real(dp), intent(in) :: diag_h(:) !! nonnegative Coleman-Li diagonal term
        real(dp), intent(in) :: g_h(:) !! gradient in scaled variables
        real(dp), intent(in) :: p_in(:) !! trust-region step in original variables
        real(dp), intent(in) :: p_h_in(:) !! trust-region step in scaled variables
        real(dp), intent(in) :: d(:) !! original-to-scaled step multipliers
        real(dp), intent(in) :: delta !! trust-region radius
        real(dp), intent(in) :: lb(:) !! lower bounds
        real(dp), intent(in) :: ub(:) !! upper bounds
        real(dp), intent(in) :: theta !! step-back factor from bounds
        real(dp), intent(out) :: step(:) !! selected original-space step
        real(dp), intent(out) :: step_h(:) !! selected scaled-space step
        real(dp), intent(out) :: predicted_reduction !! positive predicted cost reduction

        integer, allocatable :: hits(:)
        real(dp), allocatable :: ag(:), ag_h(:), p(:), p_h(:), r(:), r_h(:), x_on_bound(:)
        real(dp) :: a, ag_limit, ag_stride, ag_value, b, c, p_stride, p_value
        real(dp) :: r_stride, r_stride_l, r_stride_u, r_value, t_neg, to_bound, to_tr
        logical :: ok

        allocate(hits(size(x)), ag(size(x)), ag_h(size(x)), p(size(x)), p_h(size(x)))
        allocate(r(size(x)), r_h(size(x)), x_on_bound(size(x)))
        p = p_in
        p_h = p_h_in
        if (point_in_bounds(x + p, lb, ub)) then
            p_value = evaluate_quadratic(j_h, g_h, p_h, diag_h)
            step = p
            step_h = p_h
            predicted_reduction = -p_value
            return
        end if

        call step_size_to_bound(x, p, lb, ub, p_stride, hits)
        r_h = p_h
        where (hits /= 0) r_h = -r_h
        r = d * r_h
        p = p * p_stride
        p_h = p_h * p_stride
        x_on_bound = x + p

        call intersect_trust_region(p_h, r_h, delta, t_neg, to_tr, ok)
        if (.not. ok) to_tr = 0.0_dp
        call step_size_to_bound(x_on_bound, r, lb, ub, to_bound, hits)
        r_stride = min(to_bound, to_tr)
        if (r_stride > 0.0_dp) then
            r_stride_l = (1.0_dp - theta) * p_stride / r_stride
            if (same_real(r_stride, to_bound)) then
                r_stride_u = theta * to_bound
            else
                r_stride_u = to_tr
            end if
        else
            r_stride_l = 0.0_dp
            r_stride_u = -1.0_dp
        end if

        if (r_stride_l <= r_stride_u) then
            call build_quadratic_1d(j_h, g_h, r_h, a, b, diag_h, p_h, c)
            call minimize_quadratic_1d(a, b, r_stride_l, r_stride_u, r_stride, r_value, c)
            r_h = p_h + r_stride * r_h
            r = d * r_h
        else
            r_value = ieee_value(1.0_dp, ieee_positive_inf)
        end if

        p = theta * p
        p_h = theta * p_h
        p_value = evaluate_quadratic(j_h, g_h, p_h, diag_h)

        ag_h = -g_h
        ag = d * ag_h
        if (vector_norm2(ag_h) > 0.0_dp) then
            to_tr = delta / vector_norm2(ag_h)
            call step_size_to_bound(x, ag, lb, ub, to_bound, hits)
            if (to_bound < to_tr) then
                ag_limit = theta * to_bound
            else
                ag_limit = to_tr
            end if
            call build_quadratic_1d(j_h, g_h, ag_h, a, b, diag_h)
            call minimize_quadratic_1d(a, b, 0.0_dp, ag_limit, ag_stride, ag_value)
            ag_h = ag_stride * ag_h
            ag = ag_stride * ag
        else
            ag_value = ieee_value(1.0_dp, ieee_positive_inf)
        end if

        if (p_value < r_value .and. p_value < ag_value) then
            step = p
            step_h = p_h
            predicted_reduction = -p_value
        else if (r_value < p_value .and. r_value < ag_value) then
            step = r
            step_h = r_h
            predicted_reduction = -r_value
        else
            step = ag
            step_h = ag_h
            predicted_reduction = -ag_value
        end if
    end subroutine select_trf_step

    subroutine solve_linear_least_squares(jmat, rhs, x, ok)
        real(dp), intent(in) :: jmat(:, :) !! dense least-squares coefficient matrix
        real(dp), intent(in) :: rhs(:) !! least-squares right-hand side
        real(dp), intent(out) :: x(:) !! minimum-norm normal-equation solution
        logical, intent(out) :: ok !! true when the eigensolver converged

        real(dp), allocatable :: bmat(:, :), coeff(:), eig(:), vec(:, :)
        real(dp) :: cutoff, scale
        integer :: i, n, status

        n = size(jmat, 2)
        allocate(bmat(n, n), coeff(n), eig(n), vec(n, n))
        bmat = matmul(transpose(jmat), jmat)
        call symmetric_eigen_jacobi(bmat, eig, vec, status)
        if (status /= linalg_status_success) then
            ok = .false.
            x = 0.0_dp
            return
        end if
        coeff = matmul(transpose(vec), matmul(transpose(jmat), rhs))
        scale = maxval(abs(eig))
        x = 0.0_dp
        if (same_real(scale, 0.0_dp)) then
            ok = .true.
            return
        end if
        cutoff = (epsilon(1.0_dp) * real(max(size(jmat, 1), n), dp))**2 * scale
        do i = 1, n
            if (eig(i) > cutoff) x = x + vec(:, i) * (coeff(i) / eig(i))
        end do
        ok = .true.
    end subroutine solve_linear_least_squares

    subroutine dogleg_step_dense(x, newton_step, g, a, b, tr_bounds, lb, ub, &
            step, bound_hits, tr_hit)
        real(dp), intent(in) :: x(:) !! current feasible point
        real(dp), intent(in) :: newton_step(:) !! Gauss-Newton step
        real(dp), intent(in) :: g(:) !! active-set projected gradient
        real(dp), intent(in) :: a !! quadratic coefficient along negative gradient
        real(dp), intent(in) :: b !! linear coefficient along negative gradient
        real(dp), intent(in) :: tr_bounds(:) !! rectangular trust-region half widths
        real(dp), intent(in) :: lb(:) !! lower parameter bounds
        real(dp), intent(in) :: ub(:) !! upper parameter bounds
        real(dp), intent(out) :: step(:) !! dogleg step
        integer, intent(out) :: bound_hits(:) !! original-bound hits by component
        logical, intent(out) :: tr_hit !! true if a trust-region face was hit

        logical, allocatable :: orig_l(:), orig_u(:), tr_l(:), tr_u(:)
        integer, allocatable :: hits(:)
        real(dp), allocatable :: cauchy(:), lb_total(:), step_diff(:), ub_total(:), zero(:)
        real(dp) :: cauchy_scale, dummy_value, step_size, to_bounds
        integer :: i, n

        n = size(x)
        allocate(orig_l(n), orig_u(n), tr_l(n), tr_u(n), hits(n))
        allocate(cauchy(n), lb_total(n), step_diff(n), ub_total(n), zero(n))
        zero = 0.0_dp
        lb_total = max(lb - x, -tr_bounds)
        ub_total = min(ub - x, tr_bounds)
        do i = 1, n
            orig_l(i) = same_real(lb_total(i), lb(i) - x(i))
            orig_u(i) = same_real(ub_total(i), ub(i) - x(i))
            tr_l(i) = same_real(lb_total(i), -tr_bounds(i))
            tr_u(i) = same_real(ub_total(i), tr_bounds(i))
        end do
        bound_hits = 0
        if (point_in_bounds(newton_step, lb_total, ub_total)) then
            step = newton_step
            tr_hit = .false.
            return
        end if

        call step_size_to_bound(zero, -g, lb_total, ub_total, to_bounds, hits)
        call minimize_quadratic_1d(a, b, 0.0_dp, to_bounds, cauchy_scale, dummy_value)
        cauchy = -cauchy_scale * g
        step_diff = newton_step - cauchy
        call step_size_to_bound(cauchy, step_diff, lb_total, ub_total, step_size, hits)
        do i = 1, n
            if (hits(i) < 0 .and. orig_l(i)) bound_hits(i) = -1
            if (hits(i) > 0 .and. orig_u(i)) bound_hits(i) = 1
        end do
        tr_hit = any((hits < 0 .and. tr_l) .or. (hits > 0 .and. tr_u))
        step = cauchy + step_size * step_diff
    end subroutine dogleg_step_dense

    pure subroutine robust_loss_values(z, loss_name, rho0, rho1, rho2)
        real(dp), intent(in) :: z !! squared scaled residual
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(out) :: rho0 !! loss value
        real(dp), intent(out) :: rho1 !! first derivative with respect to z
        real(dp), intent(out) :: rho2 !! second derivative with respect to z

        real(dp) :: t

        select case (trim(loss_name))
        case ('linear')
            rho0 = z
            rho1 = 1.0_dp
            rho2 = 0.0_dp
        case ('soft_l1')
            t = 1.0_dp + z
            rho0 = 2.0_dp * (sqrt(t) - 1.0_dp)
            rho1 = 1.0_dp / sqrt(t)
            rho2 = -0.5_dp / (t * sqrt(t))
        case ('huber')
            if (z <= 1.0_dp) then
                rho0 = z
                rho1 = 1.0_dp
                rho2 = 0.0_dp
            else
                t = sqrt(z)
                rho0 = 2.0_dp * t - 1.0_dp
                rho1 = 1.0_dp / t
                rho2 = -0.5_dp / (z * t)
            end if
        case ('cauchy')
            t = 1.0_dp + z
            rho0 = log(t)
            rho1 = 1.0_dp / t
            rho2 = -1.0_dp / (t * t)
        case ('arctan')
            t = 1.0_dp + z * z
            rho0 = atan(z)
            rho1 = 1.0_dp / t
            rho2 = -2.0_dp * z / (t * t)
        case default
            rho0 = z
            rho1 = 1.0_dp
            rho2 = 0.0_dp
        end select
    end subroutine robust_loss_values

    function robust_cost(f_true, loss_name, f_scale_use, loss_fun) result(cost)
        real(dp), intent(in) :: f_true(:) !! unmodified residual vector
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(in) :: f_scale_use !! positive robust-loss scale
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss
        real(dp) :: cost

        real(dp), allocatable :: rho0(:), rho1(:), rho2(:), z(:)
        integer :: i, m

        if (trim(loss_name) == 'linear') then
            cost = 0.5_dp * dot_product(f_true, f_true)
            return
        end if
        m = size(f_true)
        allocate(z(m), rho0(m), rho1(m), rho2(m))
        z = (f_true / f_scale_use)**2
        if (trim(loss_name) == 'custom') then
            call loss_fun(z, rho0, rho1, rho2)
        else
            do i = 1, m
                call robust_loss_values(z(i), loss_name, rho0(i), rho1(i), rho2(i))
            end do
        end if
        cost = 0.5_dp * f_scale_use * f_scale_use * sum(rho0)
    end function robust_cost

    subroutine update_jac_scale(jmat, scale, scale_inv, update_existing)
        real(dp), intent(in) :: jmat(:, :) !! current residual Jacobian
        real(dp), intent(out) :: scale(:) !! variable scales from inverse column norms
        real(dp), intent(inout) :: scale_inv(:) !! column norms, monotonically nondecreasing
        logical, intent(in) :: update_existing !! retain larger norms from previous iterations

        real(dp) :: column_norm
        integer :: i

        do i = 1, size(jmat, 2)
            column_norm = sqrt(sum(jmat(:, i)**2))
            if (update_existing) then
                scale_inv(i) = max(column_norm, scale_inv(i))
            else if (same_real(column_norm, 0.0_dp)) then
                scale_inv(i) = 1.0_dp
            else
                scale_inv(i) = column_norm
            end if
        end do
        scale = 1.0_dp / scale_inv
    end subroutine update_jac_scale

    subroutine scale_robust_operator_model(f_true, loss_name, f_scale_use, f_model, &
            row_scale, cost, loss_fun)
        real(dp), intent(in) :: f_true(:) !! unmodified residual vector
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(in) :: f_scale_use !! positive robust-loss scale
        real(dp), intent(out) :: f_model(:) !! robust-scaled residual model
        real(dp), intent(out) :: row_scale(:) !! left scaling applied to the Jacobian operator
        real(dp), intent(out) :: cost !! robust objective value
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss

        real(dp), allocatable :: rho0(:), rho1(:), rho2(:), z(:)
        real(dp) :: j_scale_sq
        integer :: i, m

        if (trim(loss_name) == 'linear') then
            f_model = f_true
            row_scale = 1.0_dp
            cost = 0.5_dp * dot_product(f_true, f_true)
            return
        end if
        m = size(f_true)
        allocate(z(m), rho0(m), rho1(m), rho2(m))
        z = (f_true / f_scale_use)**2
        if (trim(loss_name) == 'custom') then
            call loss_fun(z, rho0, rho1, rho2)
        else
            do i = 1, m
                call robust_loss_values(z(i), loss_name, rho0(i), rho1(i), rho2(i))
            end do
        end if
        cost = 0.5_dp * f_scale_use * f_scale_use * sum(rho0)
        do i = 1, m
            rho2(i) = rho2(i) / (f_scale_use * f_scale_use)
            j_scale_sq = rho1(i) + 2.0_dp * rho2(i) * f_true(i) * f_true(i)
            row_scale(i) = sqrt(max(j_scale_sq, epsilon(1.0_dp)))
            f_model(i) = f_true(i) * rho1(i) / row_scale(i)
        end do
    end subroutine scale_robust_operator_model

    subroutine scale_robust_model(f_true, j_true, loss_name, f_scale_use, &
            f_model, j_model, cost, loss_fun)
        real(dp), intent(in) :: f_true(:) !! unmodified residual vector
        real(dp), intent(in) :: j_true(:, :) !! unmodified residual Jacobian
        character(len=*), intent(in) :: loss_name !! validated robust-loss name
        real(dp), intent(in) :: f_scale_use !! positive robust-loss scale
        real(dp), intent(out) :: f_model(:) !! robust-scaled residual model
        real(dp), intent(out) :: j_model(:, :) !! robust-scaled Jacobian model
        real(dp), intent(out) :: cost !! robust objective value
        procedure(robust_loss_function), optional :: loss_fun !! custom robust loss

        real(dp), allocatable :: rho0(:), rho1(:), rho2(:), z(:)
        real(dp) :: j_scale, j_scale_sq
        integer :: i, m

        if (trim(loss_name) == 'linear') then
            f_model = f_true
            j_model = j_true
            cost = 0.5_dp * dot_product(f_true, f_true)
            return
        end if
        m = size(f_true)
        allocate(z(m), rho0(m), rho1(m), rho2(m))
        z = (f_true / f_scale_use)**2
        if (trim(loss_name) == 'custom') then
            call loss_fun(z, rho0, rho1, rho2)
        else
            do i = 1, m
                call robust_loss_values(z(i), loss_name, rho0(i), rho1(i), rho2(i))
            end do
        end if
        cost = 0.5_dp * f_scale_use * f_scale_use * sum(rho0)
        do i = 1, m
            rho2(i) = rho2(i) / (f_scale_use * f_scale_use)
            j_scale_sq = rho1(i) + 2.0_dp * rho2(i) * f_true(i) * f_true(i)
            j_scale = sqrt(max(j_scale_sq, epsilon(1.0_dp)))
            f_model(i) = f_true(i) * rho1(i) / j_scale
            j_model(i, :) = j_true(i, :) * j_scale
        end do
    end subroutine scale_robust_model

    subroutine covariance_from_jacobian(jmat, cov, ok)
        real(dp), intent(in) :: jmat(:, :) !! final dense residual Jacobian
        real(dp), allocatable, intent(out) :: cov(:, :) !! pseudoinverse of J^T J
        logical, intent(out) :: ok !! true when the symmetric eigensolver converged

        real(dp), allocatable :: bmat(:, :), eig(:), vec(:, :)
        real(dp) :: cutoff, scale
        integer :: i, n, status

        n = size(jmat, 2)
        if (n <= 0) then
            ok = .false.
            return
        end if
        allocate(bmat(n, n), eig(n), vec(n, n), cov(n, n))
        bmat = matmul(transpose(jmat), jmat)
        call symmetric_eigen_jacobi(bmat, eig, vec, status)
        if (status /= linalg_status_success) then
            cov = 0.0_dp
            ok = .false.
            return
        end if
        scale = max(1.0_dp, maxval(abs(eig)))
        cutoff = (epsilon(1.0_dp) * real(max(size(jmat, 1), n), dp))**2 * scale
        cov = 0.0_dp
        do i = 1, n
            if (eig(i) > cutoff) then
                cov = cov + spread(vec(:, i), 2, n) * &
                    spread(vec(:, i), 1, n) / eig(i)
            end if
        end do
        ok = .true.
    end subroutine covariance_from_jacobian

    subroutine covariance_from_minpack(fjac_native, ipvt, cov, ok)
        real(dp), intent(in) :: fjac_native(:, :) !! MINPACK m-by-n QR workspace
        integer, intent(in) :: ipvt(:) !! one-based column pivots
        real(dp), allocatable, intent(out) :: cov(:, :) !! inverse normal matrix
        logical, intent(out) :: ok !! true if the triangular factor was invertible

        real(dp), allocatable :: invr(:, :), permuted(:, :), r(:, :)
        integer :: i, j, n

        n = size(ipvt)
        ok = .false.
        if (size(fjac_native, 1) < n .or. size(fjac_native, 2) /= n) return
        if (any(ipvt < 1) .or. any(ipvt > n)) return
        allocate(r(n, n), invr(n, n), permuted(n, n))
        r = 0.0_dp
        do i = 1, n
            do j = i, n
                r(i, j) = fjac_native(i, j)
            end do
        end do
        call invert_upper_triangular(r, invr, ok)
        if (.not. ok) return
        permuted = 0.0_dp
        do i = 1, n
            permuted(ipvt(i), :) = invr(i, :)
        end do
        allocate(cov(n, n))
        cov = matmul(permuted, transpose(permuted))
    end subroutine covariance_from_minpack

    subroutine invert_upper_triangular(r, rinv, ok)
        real(dp), intent(in) :: r(:, :) !! upper-triangular matrix
        real(dp), intent(out) :: rinv(:, :) !! inverse of r
        logical, intent(out) :: ok !! true when r is nonsingular

        real(dp), allocatable :: rhs(:), x(:)
        real(dp) :: scale, tiny_diag
        integer :: i, j, k, n

        n = size(r, 1)
        ok = .false.
        rinv = 0.0_dp
        if (size(r, 2) /= n .or. size(rinv, 1) /= n .or. size(rinv, 2) /= n) return
        scale = max(1.0_dp, maxval(abs(r)))
        tiny_diag = epsilon(1.0_dp) * scale
        if (any(abs([(r(i, i), i=1,n)]) <= tiny_diag)) return
        allocate(rhs(n), x(n))
        do j = 1, n
            rhs = 0.0_dp
            rhs(j) = 1.0_dp
            x = 0.0_dp
            do i = n, 1, -1
                x(i) = rhs(i)
                do k = i + 1, n
                    x(i) = x(i) - r(i, k) * x(k)
                end do
                x(i) = x(i) / r(i, i)
            end do
            rinv(:, j) = x
        end do
        ok = .true.
    end subroutine invert_upper_triangular

    subroutine init_least_squares_result(res, n, m)
        type(least_squares_result), intent(out) :: res !! initialized result
        integer, intent(in) :: n !! parameter count
        integer, intent(in) :: m !! residual count

        allocate(res%active_mask(max(0, n)))
        res%active_mask = 0
        res%method = 'lm'
        if (m < n) res%message = 'method=''lm'' requires residual count >= parameter count'
    end subroutine init_least_squares_result

    subroutine init_leastsq_result(res)
        type(leastsq_result), intent(out) :: res !! initialized result

        res%nfev = 0
        res%njev = -1
        res%status = 0
        res%success = .false.
        res%message = 'Improper input parameters.'
    end subroutine init_leastsq_result

    pure logical function valid_lm_tolerance(value) result(ok)
        real(dp), intent(in) :: value !! tolerance to validate

        ok = ieee_is_finite(value) .and. value > epsilon(1.0_dp)
    end function valid_lm_tolerance

    pure integer function least_squares_status(info) result(status)
        integer, intent(in) :: info !! raw MINPACK lmder status

        select case (info)
        case (1)
            status = 2
        case (2)
            status = 3
        case (3)
            status = 4
        case (4)
            status = 1
        case (5)
            status = 0
        case default
            status = -1
        end select
    end function least_squares_status

    pure function least_squares_message(status) result(message)
        integer, intent(in) :: status !! SciPy least_squares status
        character(len=160) :: message

        select case (status)
        case (-2)
            message = 'Stopped because the least-squares callback requested termination.'
        case (0)
            message = 'The maximum number of function evaluations is exceeded.'
        case (1)
            message = '`gtol` termination condition is satisfied.'
        case (2)
            message = '`ftol` termination condition is satisfied.'
        case (3)
            message = '`xtol` termination condition is satisfied.'
        case (4)
            message = 'Both `ftol` and `xtol` termination conditions are satisfied.'
        case default
            message = 'Improper input parameters status returned from `leastsq`.'
        end select
    end function least_squares_message

    function leastsq_message(info, ftol, xtol, gtol, maxfev) result(message)
        integer, intent(in) :: info !! raw MINPACK status
        real(dp), intent(in) :: ftol !! relative cost tolerance
        real(dp), intent(in) :: xtol !! relative parameter tolerance
        real(dp), intent(in) :: gtol !! orthogonality tolerance
        integer, intent(in) :: maxfev !! maximum residual evaluations
        character(len=160) :: message

        select case (info)
        case (1)
            write(message, '(a,es12.4)') 'Both actual and predicted relative reductions are at most ', ftol
        case (2)
            write(message, '(a,es12.4)') 'The relative error between two consecutive iterates is at most ', xtol
        case (3)
            write(message, '(a,es12.4,a,es12.4)') &
                'Both reduction and iterate tolerances are satisfied: ftol=', ftol, ', xtol=', xtol
        case (4)
            write(message, '(a,es12.4)') 'The residual is orthogonal to the Jacobian columns to gtol=', gtol
        case (5)
            write(message, '(a,i0,a)') 'Number of calls to function has reached maxfev = ', maxfev, '.'
        case (6)
            message = 'ftol is too small; no further reduction is possible.'
        case (7)
            message = 'xtol is too small; no further improvement is possible.'
        case (8)
            message = 'gtol is too small; residual is orthogonal to machine precision.'
        case default
            message = 'Improper input parameters.'
        end select
    end function leastsq_message

    pure logical function same_real(a, b) result(equal)
        real(dp), intent(in) :: a !! first value
        real(dp), intent(in) :: b !! second value

        equal = transfer(a, 0_int64) == transfer(b, 0_int64)
    end function same_real

    pure logical function same_vector(a, b) result(equal)
        real(dp), intent(in) :: a(:) !! first vector
        real(dp), intent(in) :: b(:) !! second vector

        if (size(a) /= size(b)) then
            equal = .false.
        else
            equal = all(transfer(a, 0_int64, size(a)) == &
                transfer(b, 0_int64, size(b)))
        end if
    end function same_vector

    pure function lowercase(text) result(out)
        character(len=*), intent(in) :: text !! source text
        character(len=len(text)) :: out

        integer :: code, i

        out = text
        do i = 1, len(text)
            code = iachar(out(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module scifort_optimize_least_squares
