! SPDX-License-Identifier: MIT
! Actual rank-two discovery products inform optional work-aware refresh costs.
program test_optimize_root_krylov_measured_rebuild
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 24
    real(dp) :: d(n), u(n,2), v(n,2), w(n), z(n), x0(n)
    integer :: i, failures, forward_calls, transpose_calls, solver_calls
    logical :: incompatible
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        d(i) = 0.9_dp + 0.035_dp*real(i,dp)
        u(i,1) = 0.10_dp + 0.065_dp*sin(0.63_dp*real(i,dp))
        u(i,2) = -0.07_dp + 0.08_dp*cos(0.38_dp*real(i,dp))
        v(i,1) = 0.13_dp*cos(0.31_dp*real(i,dp)) + 0.03_dp
        v(i,2) = 0.14_dp*sin(0.51_dp*real(i,dp)) - 0.02_dp
        w(i) = 0.17_dp*cos(0.97_dp*real(i,dp))
        z(i) = 0.15_dp*sin(0.79_dp*real(i,dp))
    end do
    x0 = 0.2_dp
    incompatible = .false.

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=1.0_dp, &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%nit == 5 .and. forward_calls == 8 .and. transpose_calls == 4, &
        'observed 4 Jv + 2 transpose Jv setup triggers one early refresh')

    call reset()
    fit = newton_krylov(residual, x0, line_search='none', maxiter=5, &
        f_tol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=1.0_dp, &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., &
        inner_m_diag_update_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(fit%nit == 5 .and. forward_calls == 4 .and. transpose_calls == 2, &
        'expensive measured transpose products suppress early refresh')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=1.0_dp, &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., &
        inner_m_diag_update_transpose_cost=huge(1.0_dp), &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(forward_calls == 4 .and. transpose_calls == 2, &
        'extremely expensive but finite transpose cost cannot overflow or trigger refresh')

    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=1.0_dp, &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(forward_calls == 4 .and. transpose_calls == 4, &
        'additional held-out transpose checks count toward measured rebuild cost')

    incompatible = .true.
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, inner_m_diag_update_interval=6, &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_stall_ratio=1.0_dp, &
        inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv)
    call check(forward_calls >= n+3 .and. forward_calls <= n+4 .and. &
        transpose_calls == 2, &
        'rejected low-rank fit counts n exact diagonal probes as rebuild work')
    incompatible = .false.

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_measured_cost=.true.)
    call invalid('measured setup cost requires work growth and nonsymmetric rank two')

    ! Measured setup work also supports ordinary exact-diagonal scaling.
    call reset()
    fit = root(residual, x0, method='krylov', line_search='none', maxiter=5, &
        fatol=1.0e-20_dp, inner_solver=work_solver, inner_m_diag_auto=.true., &
        inner_m_diag_update_interval=6, inner_m_diag_update_adaptive=.true., &
        inner_m_diag_update_stall_ratio=1.0_dp, inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_update_measured_cost=.true., jacobian_product=jvp)
    call check(fit%nit == 5 .and. forward_calls == n, &
        'measured exact diagonal cost suppresses uneconomic refresh')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_transpose_product=jtv, &
        inner_m_diag_update_measured_cost=.true., inner_m_diag_update_rebuild_cost=2)
    call invalid('manual and measured costs are mutually exclusive')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_update_transpose_cost=1.0_dp)
    call invalid('transpose cost requires measured setup mode')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_transpose_product=jtv, &
        inner_m_diag_update_measured_cost=.true., &
        inner_m_diag_update_transpose_cost=0.0_dp)
    call invalid('nonpositive transpose weight rejected')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_update_adaptive=.true., inner_m_diag_update_work_growth=2.0_dp, &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_transpose_product=jtv, &
        inner_m_diag_update_measured_cost=.true., &
        inner_m_diag_update_transpose_cost=ieee_value(0.0_dp,ieee_quiet_nan))
    call invalid('nonfinite transpose weight rejected')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_update_measured_cost=.true.)
    call invalid('non-Krylov methods reject measured setup mode')

    if (failures /= 0) error stop 1
    print '(a)', 'measured Krylov rebuild-cost tests passed'

contains

    subroutine reset()
        forward_calls = 0
        transpose_calls = 0
        solver_calls = 0
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition !! true on success
        character(len=*), intent(in) :: description !! case being tested
        if (.not. condition) then
            print '(a)', 'FAIL: '//description
            failures = failures + 1
        end if
    end subroutine check

    subroutine invalid(description)
        character(len=*), intent(in) :: description !! expected rejection
        call check(fit%status == 0 .and. .not. fit%success .and. fit%nfev == 0, description)
    end subroutine invalid

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        fvec = d*x + matmul(u,matmul(transpose(v),x))
        if (incompatible) fvec = fvec + w*dot_product(z,x)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! linearization point
        real(dp), intent(in) :: vec(:) !! vector for Jacobian product
        real(dp), intent(out) :: product(:) !! Jacobian action
        logical, intent(out) :: ok !! success flag
        forward_calls = forward_calls + 1
        product = d*vec + matmul(u,matmul(transpose(v),vec))
        if (incompatible) product = product + w*dot_product(z,vec)
        ok = .true.
    end subroutine jvp

    subroutine jtv(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! linearization point
        real(dp), intent(in) :: vec(:) !! transpose direction
        real(dp), intent(out) :: product(:) !! transpose action
        logical, intent(out) :: ok !! success flag
        transpose_calls = transpose_calls + 1
        product = d*vec + matmul(v,matmul(transpose(u),vec))
        if (incompatible) product = product + z*dot_product(w,vec)
        ok = .true.
    end subroutine jtv

    subroutine work_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! linear Jacobian product
        procedure(vector_root_linear_operator) :: preconditioner !! approximate inverse
        real(dp), intent(in) :: rhs(:) !! residual vector
        real(dp), intent(in) :: rtol !! inner tolerance
        integer, intent(in) :: maxiter !! inner work budget
        real(dp), intent(out) :: solution(:) !! damped direction
        logical, intent(out) :: ok !! true on success
        real(dp) :: tmp(size(rhs))
        integer :: j, count
        solver_calls = solver_calls + 1
        count = 1
        if (solver_calls == 2) count = 9
        do j = 1, count
            call preconditioner(rhs, tmp, ok)
            if (.not. ok) return
        end do
        solution = 0.4_dp*rhs
        ok = .true.
    end subroutine work_solver
end program test_optimize_root_krylov_measured_rebuild
