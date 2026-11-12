! SPDX-License-Identifier: MIT
! Compare LGMRES augmentation policies on changing, ill-conditioned Jacobians.
!
! The test deliberately does not assert which policy is fastest or succeeds on
! the hardest case: such rankings can change with the compiler and solver.
! It does assert honest convergence, accurate evaluation counts, and successful
! convergence for the two less severely conditioned problems.
program test_optimize_root_krylov_pruning_comparison
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: ncases = 3
    integer, parameter :: nmodes = 7
    real(dp), parameter :: conditions(ncases) = [1.0e3_dp, 1.0e4_dp, 1.0e7_dp]
    real(dp), parameter :: exact_tolerance = 1.0e-9_dp
    real(dp), parameter :: numerical_tolerance = 1.0e-8_dp
    real(dp) :: matrix(n, n)
    real(dp) :: left(n, n)
    real(dp) :: right(n, n)
    real(dp) :: left_vector(n)
    real(dp) :: right_vector(n)
    real(dp) :: coupling(n)
    real(dp) :: solution(n)
    real(dp) :: initial(n)
    real(dp) :: y_start(n)
    real(dp) :: y_finish(n)
    real(dp) :: max_residual
    real(dp) :: error_x
    type(vector_root_result) :: result
    type(vector_root_result) :: uncached_result
    integer :: case_index
    integer :: mode
    integer :: jvp_calls
    integer :: residual_calls
    integer :: reference_jvp
    integer :: exact_jvp(ncases, nmodes)
    integer :: finite_evals(nmodes)
    integer :: fails
    logical :: callback_ok

    fails = 0
    call prepare_problem()
    exact_jvp = 0
    finite_evals = 0

    print '(a)', 'condition,product,policy,success,nit,nfev,njvp,max_residual,max_x_error'
    do case_index = 1, ncases
        call build_matrix(conditions(case_index))

        ! A fixed direction sees a different Jacobian at the initial and
        ! final points. This is a *nonlinear* stress test, not a linear solve.
        call analytic_product(initial, coupling, y_start, callback_ok)
        call check(callback_ok, 'initial analytic Jacobian product')
        call analytic_product(solution, coupling, y_finish, callback_ok)
        call check(callback_ok, 'solution analytic Jacobian product')
        call check(norm2(y_start - y_finish) > 1.0e-2_dp * norm2(y_finish), &
            'the Jacobian changes as the nonlinear iterate moves')

        do mode = 1, nmodes
            call run_problem(mode, .true., result)
            exact_jvp(case_index, mode) = jvp_calls
            call verify_result(result, exact_tolerance, .true.)
            if (case_index <= 2) then
                call check(result%success, 'moderately ill-conditioned analytic solve')
                call check(maxval(abs(result%x - solution)) < 1.0e-7_dp, &
                    'moderately ill-conditioned solution accuracy')
            end if
            call report_result(conditions(case_index), 'analytic', mode, result)
        end do
    end do

    call check(exact_jvp(1, 3) > exact_jvp(1, 1), &
        'image filtering incurs additional analytic Jacobian products')
    call check(exact_jvp(2, 3) > exact_jvp(2, 1), &
        'image filtering counts Jacobian products on changing Jacobians')
    call check(exact_jvp(1, 4) > exact_jvp(1, 1), &
        'automatic image filtering incurs extra Jacobian products')
    call check(exact_jvp(1, 5) < exact_jvp(1, 4), &
        'image reuse removes duplicate analytic Jacobian products')
    call check(exact_jvp(2, 5) < exact_jvp(2, 4), &
        'reuse saves products with a changing Jacobian')
    call check(exact_jvp(1, 6) < exact_jvp(1, 4), &
        'candidate image reconstruction saves additional operator products')
    call check(exact_jvp(1, 7) < exact_jvp(1, 5) .and. &
        exact_jvp(1, 7) < exact_jvp(1, 6), &
        'reusing both old and candidate images minimizes redundant Jv calls')

    ! Finite-difference products are actual function evaluations. Verify
    ! that the extra current-Jacobian image checks appear in result%nfev.
    call build_matrix(conditions(1))
    do mode = 1, 5
        call run_problem(mode, .false., result)
        finite_evals(mode) = result%nfev
        call verify_result(result, numerical_tolerance, .false.)
        call check(result%success, 'finite-difference solve converges')
        call report_result(conditions(1), '2-point', mode, result)
    end do
    call check(finite_evals(3) > finite_evals(1), &
        'finite-difference image filtering charges extra residual evaluations')
    call check(finite_evals(4) > finite_evals(1), &
        'automatic image filtering charges extra residual evaluations')
    call check(finite_evals(5) < finite_evals(4), &
        'reused finite-difference images reduce residual evaluations')

    ! Screening the proposed step against the current-Jacobian QR span
    ! must agree with a fresh reconstruction even for changing Jacobians.
    do case_index = 1, ncases
        call build_matrix(conditions(case_index))
        call run_problem(4, .true., uncached_result)
        reference_jvp = jvp_calls
        residual_calls = 0
        jvp_calls = 0
        result = root(residual, initial, method='krylov', &
            inner_method='lgmres', inner_maxiter=12, outer_k=4, &
            jacobian_product=analytic_product, outer_image_prune_auto=.true., &
            outer_image_basis_reuse=.true., &
            inner_orthogonalization='selective', line_search='armijo', &
            fatol=exact_tolerance, maxiter=80)
        call verify_result(result, exact_tolerance, .true.)
        call check(result%success .eqv. uncached_result%success, &
            'cached image basis preserves stress convergence status')
        call check(result%nfev == uncached_result%nfev .and. &
            result%nit == uncached_result%nit .and. jvp_calls == reference_jvp, &
            'cached image basis preserves iteration and callback counts')
        call check(maxval(abs(result%x - uncached_result%x)) < 1.0e-8_dp, &
            'cached image basis matches changing-Jacobian solution')
    end do

    if (fails /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_pruning_comparison: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_pruning_comparison: PASS'

contains

    subroutine prepare_problem()
        integer :: i
        integer :: j

        do i = 1, n
            solution(i) = 0.08_dp * sin(0.3_dp * real(i, dp))
            initial(i) = solution(i) + 0.2_dp + 0.12_dp * sin(1.1_dp * real(i, dp))
            left_vector(i) = sin(0.73_dp * real(i, dp))
            right_vector(i) = cos(0.4_dp * real(i, dp))
            coupling(i) = sin(0.2_dp * real(i, dp) + 0.3_dp)
        end do
        left_vector = left_vector / norm2(left_vector)
        right_vector = right_vector / norm2(right_vector)
        coupling = coupling / norm2(coupling)

        ! Two Householder reflections give dense, orthogonal singular vectors.
        do j = 1, n
            do i = 1, n
                left(i, j) = -2.0_dp * left_vector(i) * left_vector(j)
                right(i, j) = -2.0_dp * right_vector(i) * right_vector(j)
            end do
            left(j, j) = left(j, j) + 1.0_dp
            right(j, j) = right(j, j) + 1.0_dp
        end do
    end subroutine prepare_problem

    subroutine build_matrix(condition_number)
        real(dp), intent(in) :: condition_number !! target condition number of the linear factor
        real(dp) :: diagonal(n, n)
        integer :: i

        diagonal = 0.0_dp
        do i = 1, n
            diagonal(i, i) = condition_number ** &
                (-real(i - 1, dp) / real(n - 1, dp))
        end do
        matrix = matmul(left, matmul(diagonal, transpose(right)))
    end subroutine build_matrix

    pure subroutine model_residual(x, value)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: value(:) !! model residual
        real(dp) :: delta(n)
        real(dp) :: transformed(n)
        real(dp) :: projection

        delta = x - solution
        projection = dot_product(coupling, delta)
        transformed = delta + 0.3_dp * delta * delta + &
            0.03_dp * coupling * projection * projection
        value = matmul(matrix, transformed)
    end subroutine model_residual

    subroutine residual(x, value)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: value(:) !! residual returned to root

        residual_calls = residual_calls + 1
        call model_residual(x, value)
    end subroutine residual

    subroutine analytic_product(x, vec, value, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(in) :: vec(:) !! Jacobian-vector direction
        real(dp), intent(out) :: value(:) !! exact Jacobian-vector product
        logical, intent(out) :: ok !! whether output is finite
        real(dp) :: delta(n)
        real(dp) :: product(n)
        real(dp) :: projection

        jvp_calls = jvp_calls + 1
        delta = x - solution
        projection = dot_product(coupling, delta)
        product = (1.0_dp + 0.6_dp * delta) * vec + &
            0.06_dp * coupling * projection * dot_product(coupling, vec)
        value = matmul(matrix, product)
        ok = all(ieee_is_finite(value))
    end subroutine analytic_product

    subroutine run_problem(policy, use_exact_product, res)
        integer, intent(in) :: policy !! 1 none, 2 solution, 3 image, 4 auto, 5 old reused, 6 candidate, 7 both
        logical, intent(in) :: use_exact_product !! bypass finite differences
        type(vector_root_result), intent(out) :: res !! nonlinear root result

        residual_calls = 0
        jvp_calls = 0
        if (use_exact_product) then
            select case (policy)
            case (1)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            case (2)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, outer_prune_tol=1.0e-7_dp, &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            case (3)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, outer_image_prune_tol=1.0e-7_dp, &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            case (4)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, outer_image_prune_auto=.true., &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            case (5)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, outer_image_prune_auto=.true., &
                    outer_image_reuse=.true., &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            case (6)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, outer_image_prune_auto=.true., &
                    outer_candidate_image_reuse=.true., &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            case (7)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    jacobian_product=analytic_product, outer_image_prune_auto=.true., &
                    outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=exact_tolerance, maxiter=80)
            end select
        else
            select case (policy)
            case (1)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=numerical_tolerance, maxiter=80)
            case (2)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    outer_prune_tol=1.0e-7_dp, &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=numerical_tolerance, maxiter=80)
            case (3)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    outer_image_prune_tol=1.0e-7_dp, &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=numerical_tolerance, maxiter=80)
            case (4)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    outer_image_prune_auto=.true., &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=numerical_tolerance, maxiter=80)
            case (5)
                res = root(residual, initial, method='krylov', &
                    inner_method='lgmres', inner_maxiter=12, outer_k=4, &
                    outer_image_prune_auto=.true., outer_image_reuse=.true., &
                    inner_orthogonalization='selective', line_search='armijo', &
                    fatol=numerical_tolerance, maxiter=80)
            end select
        end if
    end subroutine run_problem

    subroutine verify_result(res, tolerance, with_analytic_product)
        type(vector_root_result), intent(in) :: res !! result to verify independently
        real(dp), intent(in) :: tolerance !! requested maximum residual
        logical, intent(in) :: with_analytic_product !! true for supplied J*v
        real(dp) :: f_check(n)
        real(dp) :: bound

        call check(allocated(res%x) .and. allocated(res%fun), 'result arrays allocated')
        if (.not. allocated(res%x) .or. .not. allocated(res%fun)) return
        call check(all(ieee_is_finite(res%x)) .and. all(ieee_is_finite(res%fun)), &
            'finite final nonlinear iterate and residual')
        call check(res%nfev == residual_calls, 'function evaluation count matches calls')
        call check(res%nit >= 2 .and. res%nit <= 80, 'reasonable nonlinear step count')
        if (with_analytic_product) then
            call check(jvp_calls > 0, 'analytic products were requested')
        else
            call check(jvp_calls == 0, 'finite difference does not call analytic product')
        end if
        call model_residual(res%x, f_check)
        bound = max(1.0e-14_dp, 2.0e-13_dp * maxval(abs(f_check)))
        call check(maxval(abs(f_check - res%fun)) <= bound, &
            'reported residual matches an independent evaluation')
        if (res%success) then
            call check(maxval(abs(f_check)) <= tolerance * (1.0_dp + 1.0e-5_dp), &
                'reported success satisfies the actual residual criterion')
        end if
    end subroutine verify_result

    subroutine report_result(condition_number, product_kind, policy, res)
        real(dp), intent(in) :: condition_number !! nominal linear-factor condition
        character(len=*), intent(in) :: product_kind !! analytic or finite-difference
        integer, intent(in) :: policy !! 1 none, 2 solution, 3 image, 4 auto, 5 old reused, 6 candidate, 7 both
        type(vector_root_result), intent(in) :: res !! result to report
        character(len=16) :: policy_name
        character(len=1) :: comma

        select case (policy)
        case (1)
            policy_name = 'none'
        case (2)
            policy_name = 'solution-space'
        case (3)
            policy_name = 'jacobian-image'
        case (4)
            policy_name = 'image-auto'
        case (5)
            policy_name = 'image-reuse'
        case (6)
            policy_name = 'candidate-reuse'
        case (7)
            policy_name = 'both-reuse'
        end select
        comma = ','
        max_residual = maxval(abs(res%fun))
        error_x = maxval(abs(res%x - solution))
        write (*, '(es10.2, a, a, a, a, a, l1, a, i0, a, i0, a, i0, a, es13.5, a, es13.5)') &
            condition_number, comma, trim(product_kind), comma, trim(policy_name), comma, &
            res%success, comma, res%nit, comma, res%nfev, comma, jvp_calls, comma, &
            max_residual, comma, error_x
    end subroutine report_result

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: label !! human-readable assertion description

        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

end program test_optimize_root_krylov_pruning_comparison
