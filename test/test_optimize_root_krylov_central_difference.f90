! SPDX-License-Identifier: MIT
! SciFort three-point matrix-free Jacobian products and Krylov API coverage.
program test_optimize_root_krylov_central_difference
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    real(dp), parameter :: xstart(2) = [0.8_dp, -0.3_dp]
    real(dp), parameter :: direction(2) = [2.0_dp, -1.0_dp]
    real(dp), parameter :: target(4) = [1.0_dp, -2.0_dp, 3.0_dp, -1.0_dp]
    real(dp) :: a(4, 4), rhs(4)
    real(dp) :: error_forward, error_central, error_default
    real(dp) :: error_central_auto, error_forward_auto
    real(dp) :: measured_error
    type(vector_root_result) :: res
    integer :: failures, probes, i
    character(len=8), parameter :: methods(6) = [character(len=8) :: &
        'gmres', 'lgmres', 'bicgstab', 'tfqmr', 'cgs', 'minres']

    failures = 0
    probes = 0

    ! The custom linear solver exposes the actual J*v approximation without
    ! conflating truncation error with nonlinear/inner solver convergence.
    res = root(curved_residual, xstart, method='krylov', inner_solver=measure_jv, &
        jacobian_diff_scheme='2-point', rdiff=1.0e-4_dp)
    error_forward = measured_error
    call check(.not. res%success .and. res%nfev == 2, &
        'two-point directional product costs one extra evaluation')
    call check(error_forward > 1.0e-6_dp, 'forward truncation error observable')

    res = root(curved_residual, xstart, method='krylov', inner_solver=measure_jv, &
        jacobian_diff_scheme='3-point', rdiff=1.0e-4_dp)
    error_central = measured_error
    call check(.not. res%success .and. res%nfev == 3, &
        'three-point directional product costs two extra evaluations')
    call check(error_central < error_forward / 100.0_dp .and. error_central < 1.0e-8_dp, &
        'centered stencil improves directional derivative accuracy')

    res = root(curved_residual, xstart, method='krylov', inner_solver=measure_jv)
    error_default = measured_error
    call check(res%nfev == 2, 'omitted stencil retains SciPy two-point work count')
    call check(error_default > 0.0_dp, 'default two-point product evaluated')

    res = root(curved_residual, xstart, method='krylov', inner_solver=measure_jv, &
        jacobian_diff_scheme='3-point')
    error_central_auto = measured_error
    call check(res%nfev == 3, 'automatic three-point step counts both samples')
    call check(error_central_auto < 1.0e-8_dp, 'automatic O(eps**(1/3)) centered step')

    res = root(curved_residual, xstart, method='krylov', inner_solver=measure_jv, &
        jacobian_diff_scheme='2-point')
    error_forward_auto = measured_error
    call check(error_central_auto < error_forward_auto, &
        'automatic three-point step improves on automatic forward step')

    res = newton_krylov(curved_residual, xstart, inner_solver=measure_jv, &
        jacobian_diff_scheme='3-point', rdiff=1.0e-4_dp)
    call check(res%nfev == 3 .and. measured_error < 1.0e-8_dp, &
        'newton_krylov forwards three-point stencil and rdiff')
    call check(probes == 6, 'all six probing callbacks invoked')

    ! Public built-in paths use the same stencil as custom inner solvers.
    a = 0.0_dp
    do i = 1, 4
        a(i, i) = real(i + 3, dp)
    end do
    do i = 1, 3
        a(i, i + 1) = 0.1_dp
        a(i + 1, i) = 0.1_dp
    end do
    rhs = matmul(a, target)
    do i = 1, size(methods)
        res = root(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            method='krylov', inner_method=trim(methods(i)), &
            jacobian_diff_scheme='3-point', inner_maxiter=30, &
            line_search='none', fatol=1.0e-8_dp)
        call check(res%success .and. maxval(abs(res%fun)) < 1.0e-8_dp, &
            trim(methods(i))//' converges with centered products')
        call check(maxval(abs(res%x - target)) < 1.0e-7_dp, &
            trim(methods(i))//' solution accuracy with centered products')
    end do

    ! Verify the new stencil also composes with the prior TFQMR
    ! right-preconditioning and true-residual safeguards.
    res = root(linear_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_method='tfqmr', inner_m=inverse_diagonal, &
        inner_tfqmr_right_preconditioned=.true., inner_verify_residual=.true., &
        jacobian_diff_scheme='3-point', inner_maxiter=30, &
        line_search='none', fatol=1.0e-8_dp)
    call check(res%success .and. maxval(abs(res%fun)) < 1.0e-8_dp, &
        'centered TFQMR composes with verified right preconditioning')

    res = root(curved_residual, xstart, method='broyden1', &
        jacobian_diff_scheme='3-point')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') > 0, &
        'non-Krylov methods reject stencil option before calling residual')

    res = newton_krylov(curved_residual, xstart, jacobian_diff_scheme='5-point')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'must be 2-point or 3-point') > 0, &
        'unsupported stencil rejected before evaluating residual')

    ! The negative-direction sample may leave a function's domain. It must
    ! report an unsuccessful product and include that evaluation in nfev.
    res = root(boundary_residual, [1.0e-5_dp, 0.0_dp], &
        method='krylov', inner_solver=probe_boundary, &
        jacobian_diff_scheme='3-point', rdiff=1.0e-4_dp)
    call check(.not. res%success .and. res%nfev == 3, &
        'nonfinite backward stencil evaluation is counted and rejected')

    res = root(curved_residual, xstart, method='krylov', &
        inner_solver=probe_tiny_step, jacobian_diff_scheme='3-point', &
        rdiff=1.0e-30_dp)
    call check(.not. res%success .and. res%nfev == 1, &
        'unrepresentable centered step rejected without evaluating duplicate x')

    res = root(opposite_large_residual, [0.0_dp], method='krylov', &
        inner_solver=probe_large_difference, jacobian_diff_scheme='3-point', &
        rdiff=1.0_dp)
    call check(.not. res%success .and. res%nfev == 3, &
        'centered difference avoids overflow in subtraction of large samples')

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_root_krylov_central_difference: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_central_difference: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! property to verify
        character(len=*), intent(in) :: label !! failure description

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine curved_residual(x, f)
        real(dp), intent(in) :: x(:) !! independent variables
        real(dp), intent(out) :: f(:) !! nonlinear square residual vector

        f(1) = exp(x(1)) + 0.1_dp * x(2)**3 - 1.0_dp
        f(2) = sin(x(2)) + x(1)**2 - 1.0_dp
    end subroutine curved_residual

    subroutine measure_jv(matvec, preconditioner, b, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! numerical Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse-Jacobian action
        real(dp), intent(in) :: b(:) !! right-hand side supplied by outer Newton
        real(dp), intent(in) :: rtol !! relative linear tolerance
        integer, intent(in) :: maxiter !! Krylov iteration budget
        real(dp), intent(out) :: solution(:) !! unused Newton correction
        logical, intent(out) :: ok !! deliberately false after probing

        real(dp) :: measured(2), exact(2), zero_product(2)
        logical :: good

        probes = probes + 1
        call matvec([0.0_dp, 0.0_dp], zero_product, good)
        call check(good .and. all(zero_product == 0.0_dp), &
            'zero direction is exact with no residual calls')
        call matvec(direction, measured, good)
        call check(good, 'finite matrix-free directional product')
        exact(1) = exp(xstart(1)) * direction(1) + &
            0.3_dp * xstart(2)**2 * direction(2)
        exact(2) = cos(xstart(2)) * direction(2) + &
            2.0_dp * xstart(1) * direction(1)
        measured_error = maxval(abs(measured - exact))
        solution = 0.0_dp
        ok = .false.
    end subroutine measure_jv

    subroutine linear_residual(x, f)
        real(dp), intent(in) :: x(:) !! independent variables
        real(dp), intent(out) :: f(:) !! linear square residual vector

        f = matmul(a, x) - rhs
    end subroutine linear_residual

    subroutine inverse_diagonal(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! input direction
        real(dp), intent(out) :: product(:) !! diagonal approximate inverse action
        logical, intent(out) :: ok !! true when application succeeds

        product = vec / [a(1, 1), a(2, 2), a(3, 3), a(4, 4)]
        ok = .true.
    end subroutine inverse_diagonal

    subroutine boundary_residual(x, f)
        real(dp), intent(in) :: x(:) !! evaluation location, first component positive
        real(dp), intent(out) :: f(:) !! residual, NaN outside domain

        if (x(1) <= 0.0_dp) then
            f = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            f = x - [0.5_dp, 0.0_dp]
        end if
    end subroutine boundary_residual

    subroutine probe_boundary(matvec, preconditioner, b, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! numerical Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse-Jacobian action
        real(dp), intent(in) :: b(:) !! right-hand side supplied by outer Newton
        real(dp), intent(in) :: rtol !! relative linear tolerance
        integer, intent(in) :: maxiter !! Krylov iteration budget
        real(dp), intent(out) :: solution(:) !! unused Newton correction
        logical, intent(out) :: ok !! deliberately false after probing

        real(dp) :: measured(2)
        logical :: good

        call matvec([1.0_dp, 0.0_dp], measured, good)
        call check(.not. good, 'invalid negative stencil sample propagates failure')
        solution = 0.0_dp
        ok = .false.
    end subroutine probe_boundary

    subroutine probe_tiny_step(matvec, preconditioner, b, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! numerical Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse-Jacobian action
        real(dp), intent(in) :: b(:) !! nonlinear residual passed as right-hand side
        real(dp), intent(in) :: rtol !! relative linear tolerance
        integer, intent(in) :: maxiter !! Krylov iteration budget
        real(dp), intent(out) :: solution(:) !! unused Newton correction
        logical, intent(out) :: ok !! false to terminate the probing call

        real(dp) :: measured(2)
        logical :: good

        call matvec([1.0_dp, 0.0_dp], measured, good)
        call check(.not. good, 'tiny step is not misinterpreted as a zero Jacobian')
        solution = 0.0_dp
        ok = .false.
    end subroutine probe_tiny_step

    subroutine opposite_large_residual(x, f)
        real(dp), intent(in) :: x(:) !! one-variable evaluation location
        real(dp), intent(out) :: f(:) !! finite, large and opposite-valued samples

        if (x(1) > 0.0_dp) then
            f(1) = 1.0e308_dp
        else if (x(1) < 0.0_dp) then
            f(1) = -1.0e308_dp
        else
            f(1) = 1.0_dp
        end if
    end subroutine opposite_large_residual

    subroutine probe_large_difference(matvec, preconditioner, b, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! numerical Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse-Jacobian action
        real(dp), intent(in) :: b(:) !! nonlinear residual passed as right-hand side
        real(dp), intent(in) :: rtol !! relative linear tolerance
        integer, intent(in) :: maxiter !! Krylov iteration budget
        real(dp), intent(out) :: solution(:) !! unused Newton correction
        logical, intent(out) :: ok !! false to terminate the probing call

        real(dp) :: measured(1)
        logical :: good

        call matvec([1.0_dp], measured, good)
        call check(good, 'large signed sample difference is finite')
        if (good) then
            call check(ieee_is_finite(measured(1)) .and. &
                abs(measured(1) / 1.0e308_dp - 1.0_dp) < 1.0e-14_dp, &
                'halved central difference avoids intermediate overflow')
        end if
        solution = 0.0_dp
        ok = .false.
    end subroutine probe_large_difference
end program test_optimize_root_krylov_central_difference
