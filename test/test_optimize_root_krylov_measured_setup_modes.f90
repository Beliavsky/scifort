! SPDX-License-Identifier: MIT
! Measured refresh costs include every automatic diagonal/low-rank setup path.
program test_optimize_root_krylov_measured_setup_modes
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result, vector_root_linear_operator
    implicit none

    integer, parameter :: n = 24
    real(dp) :: d(n), u(n,2), v(n,2), w(n), z(n), x0(n)
    integer :: i, mode, failures, forward_calls, solver_calls
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        d(i) = 0.90_dp + 0.013_dp*real(i,dp)
        u(i,1) = 0.10_dp + 0.065_dp*sin(0.63_dp*real(i,dp))
        u(i,2) = -0.07_dp + 0.08_dp*cos(0.38_dp*real(i,dp))
        v(i,1) = 0.13_dp*cos(0.31_dp*real(i,dp)) + 0.03_dp
        v(i,2) = 0.14_dp*sin(0.51_dp*real(i,dp)) - 0.02_dp
        w(i) = 0.20_dp*cos(0.97_dp*real(i,dp))
        z(i) = 0.18_dp*sin(0.79_dp*real(i,dp))
    end do
    x0 = 0.2_dp

    mode = 1
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n, &
        'exact diagonal measures n products and avoids costly refresh')

    mode = 2
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_one=.true., &
        inner_m_diag_rank_one_discover=.true., jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n, &
        'dense discovered rank one reuses its n Jacobian columns')

    mode = 3
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., inner_m_diag_rank_two_stream=.true., &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n, &
        'streamed rank two counts each column once including pilot columns')

    mode = 4
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v, jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n+4, &
        'known-right rank two counts diagonal and four fitting/check products')

    mode = 5
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_one=.true., &
        inner_m_diag_rank_one_right=v(:,1), jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n+3, &
        'known-right rank one counts diagonal and three fitting/check products')

    mode = 6
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v, inner_m_diag_rank_two_sketch=.true., &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == 14, &
        'seven-product rank-two sketch is economical enough to trigger a refresh')

    mode = 7
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., inner_m_diag_rank_two_symmetric_base=d, &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == 8, &
        'symmetric low-product discovery counts four products per refresh')

    mode = 8
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_right=v, inner_m_diag_rank_two_sketch=.true., &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n+9, &
        'rejected rank-two sketch counts probes, failed fit, exact fallback and retry')

    mode = 1
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_probes=2, &
        jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == 4, &
        'two-probe diagonal estimator permits early measured-cost refresh')

    ! Numeric-Jacobian setups count directions, not callback evaluations.
    ! Two-point requires one residual evaluation per Jv; three-point two.
    mode = 1
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_probes=2, &
        jacobian_diff_scheme='2-point')
    call check(fit%nit == 5 .and. fit%nfev == 10, &
        'two-point numerical probes use one residual evaluation per Jv')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_probes=2, &
        jacobian_diff_scheme='3-point')
    call check(fit%nit == 5 .and. fit%nfev == 14, &
        'three-point numerical probes use two residual evaluations per Jv')

    if (failures /= 0) error stop 1
    print '(a)', 'measured Krylov setup modes tests passed'

contains

    subroutine reset()
        forward_calls = 0
        solver_calls = 0
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! true for a passing regression
        character(len=*), intent(in) :: description !! description of expected behavior
        if (.not. condition) then
            print '(a,i0,a,i0,a,i0,a,i0)', 'FAIL: mode=',mode, &
                ' forward=',forward_calls,' nit=',fit%nit,' status=',fit%status
            print '(a)', description
            failures = failures + 1
        end if
    end subroutine check

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! residual vector
        select case (mode)
        case (5)
            fvec = d*x + u(:,1)*dot_product(v(:,1),x)
        case (7)
            fvec = d*x + matmul(u,matmul(transpose(u),x))
        case default
            fvec = d*x + matmul(u,matmul(transpose(v),x))
        end select
        if (mode == 8) fvec = fvec + w*dot_product(z,x)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current linearization point
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector action
        logical, intent(out) :: ok !! true on successful evaluation
        forward_calls = forward_calls + 1
        call residual(vec, product)
        ok = .true.
    end subroutine jvp

    subroutine work_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian operator
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs(:) !! current linear residual
        real(dp), intent(in) :: rtol !! inner relative tolerance
        integer, intent(in) :: maxiter !! inner iteration budget
        real(dp), intent(out) :: solution(:) !! damped step
        logical, intent(out) :: ok !! whether every preconditioning action succeeded
        real(dp) :: tmp(size(rhs))
        integer :: j, count
        solver_calls = solver_calls + 1
        count = 1
        if (solver_calls == 2) count = 12
        do j = 1, count
            call preconditioner(rhs, tmp, ok)
            if (.not. ok) return
        end do
        solution = 0.4_dp*rhs
        ok = .true.
    end subroutine work_solver
end program test_optimize_root_krylov_measured_setup_modes
