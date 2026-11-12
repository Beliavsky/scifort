! SPDX-License-Identifier: MIT
! Regression for memory-bounded restarted GMRES in the Krylov root interface.
program test_optimize_root_krylov_gmres_restart
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    integer, parameter :: n = 32
    real(dp) :: x0(n), target(n)
    type(vector_root_result) :: baseline, equal_restart, restarted, via_newton, invalid
    integer :: i, failures, nfev_calls, jvp_calls, precond_calls
    integer :: baseline_nfev, baseline_jvp, baseline_precond
    integer :: restart_jvp

    failures = 0
    do i = 1, n
        target(i) = 0.75_dp + 0.3_dp * sin(0.37_dp * real(i, dp))
    end do
    x0 = target + 0.6_dp

    call clear_counts()
    baseline = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=30, inner_orthogonalization='selective', &
        jacobian_product=jvp, inner_m=preconditioner, &
        line_search='none', fatol=1.0e-10_dp, maxiter=35)
    call verify(baseline, 'original GMRES')
    baseline_nfev = nfev_calls
    baseline_jvp = jvp_calls
    baseline_precond = precond_calls

    call clear_counts()
    equal_restart = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=30, inner_restart=30, inner_orthogonalization='selective', &
        jacobian_product=jvp, inner_m=preconditioner, &
        line_search='none', fatol=1.0e-10_dp, maxiter=35)
    call verify(equal_restart, 'equal restart uses original path')
    call check(equal_restart%nfev == baseline_nfev, 'equal restart nfev unchanged')
    call check(jvp_calls == baseline_jvp, 'equal restart Jv unchanged')
    call check(precond_calls == baseline_precond, 'equal restart preconditioner unchanged')
    call check(all(equal_restart%x == baseline%x), 'equal restart is bitwise identical')

    call clear_counts()
    restarted = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=30, inner_restart=4, inner_orthogonalization='selective', &
        jacobian_product=jvp, inner_m=preconditioner, &
        inner_check_residual=.true., &
        line_search='none', fatol=1.0e-10_dp, maxiter=35)
    call verify(restarted, 'restarted and preconditioned GMRES')
    call check(maxval(abs(restarted%x - baseline%x)) < 1.0e-8_dp, &
        'restarted and baseline roots agree')
    restart_jvp = jvp_calls

    call clear_counts()
    via_newton = newton_krylov(residual, x0, inner_method='gmres', &
        inner_restart=4, inner_maxiter=30, jacobian_product=jvp, &
        inner_m=preconditioner, inner_orthogonalization='twice', &
        line_search='none', f_tol=1.0e-10_dp, maxiter=35)
    call verify(via_newton, 'newton_krylov forwarding')
    call check(jvp_calls > 0 .and. restart_jvp > 0, 'both paths exercise Jacobian')

    call clear_counts()
    invalid = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_restart=0)
    call rejected(invalid, 'reject zero restart')
    call clear_counts()
    invalid = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_restart=-3)
    call rejected(invalid, 'reject negative restart')
    ! LGMRES now accepts an optional restart width as well.
    call clear_counts()
    invalid = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=30, inner_restart=4, outer_k=6, &
        jacobian_product=jvp, inner_m=preconditioner, &
        line_search='none', fatol=1.0e-10_dp, maxiter=35)
    call verify(invalid, 'accept LGMRES restart')
    call clear_counts()
    invalid = root(residual, x0, method='hybr', inner_restart=4)
    call rejected(invalid, 'reject non-Krylov restart')

    ! A short total budget permits at most seven Arnoldi operator calls
    ! plus one residual verification per (at most three) restart cycles.
    call clear_counts()
    invalid = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=7, inner_restart=3, jacobian_product=jvp, &
        line_search='none', fatol=1.0e-40_dp, maxiter=1)
    call check(jvp_calls >= 1 .and. jvp_calls <= 10, &
        'restart honors total Arnoldi budget')
    call check(invalid%nfev == nfev_calls, 'budget run counts residual evaluations')

    ! Without an exact Jv callback, every additional restart-cycle check
    ! must count the extra finite-difference residual evaluations.
    call clear_counts()
    restarted = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=30, inner_restart=5, line_search='none', &
        fatol=1.0e-9_dp, maxiter=35)
    call check(restarted%success, 'finite-difference restarted solve converges')
    call check(restarted%nfev == nfev_calls, 'finite-difference nfev accounting')
    call check(maxval(abs(restarted%x - target)) < 2.0e-7_dp, &
        'finite-difference restarted root accuracy')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_gmres_restart: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_gmres_restart: PASS'
contains
    subroutine clear_counts()
        nfev_calls = 0
        jvp_calls = 0
        precond_calls = 0
    end subroutine

    subroutine check(condition, label)
        logical, intent(in) :: condition !! condition to verify
        character(len=*), intent(in) :: label !! diagnostic label
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine

    subroutine rejected(r, label)
        type(vector_root_result), intent(in) :: r !! invalid-parameter result
        character(len=*), intent(in) :: label !! diagnostic label
        call check(.not. r%success, trim(label)//' returns failure')
        call check(nfev_calls == 0, trim(label)//' rejects before residual call')
        call check(r%nfev == 0, trim(label)//' reports zero residual evaluations')
    end subroutine

    subroutine verify(r, label)
        type(vector_root_result), intent(in) :: r !! nonlinear root result
        character(len=*), intent(in) :: label !! case name
        call check(r%success, trim(label)//' converges')
        call check(all(ieee_is_finite(r%x)), trim(label)//' has finite solution')
        call check(maxval(abs(r%x - target)) < 1.0e-8_dp, &
            trim(label)//' solves original system')
        call check(r%nfev == nfev_calls, trim(label)//' counts residual calls')
        call check(jvp_calls > 0 .and. precond_calls > 0, &
            trim(label)//' calls Jacobian and preconditioner')
    end subroutine

    subroutine apply_linear(v, product)
        real(dp), intent(in) :: v(:) !! argument
        real(dp), intent(out) :: product(:) !! linear transformation
        integer :: j
        product = 0.0_dp
        do j = 1, n
            product(j) = (2.4_dp + 0.07_dp * real(j, dp)) * v(j)
            if (j > 1) product(j) = product(j) - 0.24_dp * v(j - 1)
            if (j < n) product(j) = product(j) + 0.18_dp * v(j + 1)
        end do
        product = product + 0.04_dp * sum(v) / real(n, dp)
    end subroutine

    subroutine residual(x, values)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: values(:) !! nonlinear residual
        nfev_calls = nfev_calls + 1
        call apply_linear(x - target, values)
        values = values + 0.03_dp * (x*x - target*target)
    end subroutine

    subroutine jvp(x, v, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: v(:) !! direction
        real(dp), intent(out) :: product(:) !! J(x)*v
        logical, intent(out) :: ok !! callback success
        jvp_calls = jvp_calls + 1
        call apply_linear(v, product)
        product = product + 0.06_dp * x * v
        ok = all(ieee_is_finite(product))
    end subroutine

    subroutine preconditioner(v, product, ok)
        real(dp), intent(in) :: v(:) !! residual
        real(dp), intent(out) :: product(:) !! approximate inverse-Jacobian action
        logical, intent(out) :: ok !! callback success
        integer :: j
        precond_calls = precond_calls + 1
        do j = 1, n
            product(j) = v(j) / (2.4_dp + 0.07_dp * real(j, dp) + 0.06_dp * x0(j))
        end do
        ok = all(ieee_is_finite(product))
    end subroutine
end program test_optimize_root_krylov_gmres_restart
