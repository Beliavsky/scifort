! SPDX-License-Identifier: MIT
! Cost-aware adjoint verification for nonsymmetric rank-two discovery.
program test_optimize_root_krylov_rank_two_transpose_adaptive
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, &
        vector_root_result, vector_root_linear_operator
    implicit none

    integer, parameter :: n = 24
    real(dp) :: d(n), u(n,2), v(n,2), target(n), x0(n)
    real(dp) :: probe(n), rhs_probe(n), diagonal(n), observed(n)
    real(dp) :: p(n,2), q(n,2), spurious(n), forward_probe(n), theta, norm_value
    integer :: i, j, k, failures, jv_calls, jtv_calls, inner_calls, fail_jtv_at
    logical :: adversarial, borderline_forward
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        theta = real(i,dp)
        d(i) = 1.0_dp + 0.02_dp*theta
        u(i,1) = 0.10_dp + 0.06_dp*sin(0.63_dp*theta)
        u(i,2) = -0.08_dp + 0.09_dp*cos(0.38_dp*theta)
        v(i,1) = 0.13_dp*cos(0.31_dp*theta) + 0.03_dp
        v(i,2) = 0.14_dp*sin(0.51_dp*theta) - 0.02_dp
        target(i) = 0.03_dp*cos(0.46_dp*theta)
        probe(i) = 0.1_dp+sin(0.19_dp*theta)
        p(i,1) = sin(0.53_dp*theta) + cos(0.93_dp*theta)
        p(i,2) = cos(1.31_dp*theta) - sin(0.79_dp*theta)
        spurious(i) = cos(0.27_dp*theta) + sin(1.91_dp*theta)
        q(i,1) = sin(0.67_dp*theta) + cos(1.17_dp*theta)
        q(i,2) = sin(1.41_dp*theta) - cos(0.43_dp*theta)
        forward_probe(i) = sin(0.31_dp*theta) + cos(1.73_dp*theta)
    end do
    do j = 1, 2
        do k = 1, j-1
            p(:,j) = p(:,j) - dot_product(p(:,j),p(:,k))*p(:,k)
        end do
        p(:,j) = p(:,j)/sqrt(sum(p(:,j)**2))
    end do
    do j = 1, 2
        spurious = spurious - dot_product(spurious,p(:,j))*p(:,j)
    end do
    do j = 1, 2
        do k = 1, j-1
            q(:,j) = q(:,j) - dot_product(q(:,j),q(:,k))*q(:,k)
        end do
        q(:,j) = q(:,j)/sqrt(sum(q(:,j)**2))
        forward_probe = forward_probe - dot_product(forward_probe,q(:,j))*q(:,j)
    end do
    forward_probe = forward_probe/sqrt(sum(forward_probe**2))
    norm_value = sqrt(sum(spurious**2))
    if (norm_value <= 0.1_dp) error stop 'invalid test directions'
    spurious = spurious/norm_value
    x0 = 0.0_dp
    diagonal = d + u(:,1)*v(:,1) + u(:,2)*v(:,2)
    call action(probe,rhs_probe)

    adversarial = .false.
    borderline_forward = .false.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(inner_calls == 1 .and. jv_calls == 4 .and. jtv_calls == 2, &
        'default mode retains four forward and two transpose products')
    call check(maxval(abs(observed-probe)) < 2.0e-11_dp, &
        'default verified inverse acts accurately')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(inner_calls == 1 .and. jv_calls == 4 .and. jtv_calls == 4, &
        'opt-in held-out checks add exactly two transpose products')
    call check(maxval(abs(observed-probe)) < 2.0e-11_dp, &
        'held-out checks retain exact Woodbury inverse')

    adversarial = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jtv_calls == 2 .and. maxval(abs(observed-probe)) < 2.0e-11_dp, &
        'adversarial callback is invisible along both fitting directions')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jtv_calls == 3 .and. jv_calls >= n+4 .and. inner_calls == 1, &
        'bad held-out transpose direction triggers exact-diagonal fallback')
    call check(maxval(abs(observed-rhs_probe/diagonal)) < 2.0e-11_dp, &
        'fallback uses true full Jacobian diagonal')

    adversarial = .false.
    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jtv_calls == 4 .and. inner_calls == 1, &
        'newton_krylov forwards transpose-verification option')

    call reset()
    fail_jtv_at = 3
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true., &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(fit%status == 0 .and. .not. fit%success .and. inner_calls == 0, &
        'failed held-out transpose callback halts the linear solve')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true.)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'transpose verification without nonsymmetric base rejected early')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true.)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'transpose verification rejected by non-krylov methods')

    ! An inexpensive transpose callback selects both held-out checks.
    adversarial = .false.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=1.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 4 .and. jtv_calls == 4 .and. inner_calls == 1, &
        'cheap adaptive policy uses two held-out adjoint checks')
    call check(maxval(abs(observed-probe)) < 2.0e-11_dp, &
        'cheap adaptive policy retains inverse accuracy')

    ! A costly transpose callback uses one held-out verification at n=24.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 4 .and. jtv_calls == 3 .and. inner_calls == 1, &
        'expensive adaptive policy uses one held-out adjoint check')
    call check(maxval(abs(observed-probe)) < 2.0e-11_dp, &
        'expensive adaptive policy retains inverse accuracy')

    ! A near-tolerance forward check mandates both adjoint directions even
    ! when the declared adjoint cost would otherwise favor just one.
    borderline_forward = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 4 .and. jtv_calls == 4 .and. inner_calls == 1, &
        'borderline forward consistency forces both held-out adjoint checks')
    borderline_forward = .false.

    ! The first held-out adjoint catches bad callbacks even on the cheap path.
    adversarial = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jtv_calls == 3 .and. jv_calls >= n+4 .and. inner_calls == 1, &
        'adaptive policy rejects inconsistent held-out adjoint')
    call check(maxval(abs(observed-rhs_probe/diagonal)) < 2.0e-11_dp, &
        'adaptive verification failure falls back to exact diagonal')

    adversarial = .false.
    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jtv_calls == 3 .and. inner_calls == 1, &
        'newton_krylov forwards adaptive transpose cost')

    call reset()
    fail_jtv_at = 3
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
        jacobian_product=jvp, jacobian_transpose_product=jtv, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(fit%status == 0 .and. .not. fit%success .and. inner_calls == 0, &
        'adaptive held-out callback failure stops solve')

    call reset()
    fit = root(residual, x0, method='krylov', &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=1.0_dp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'adaptive transpose cost without nonsymmetric base rejected early')

    call reset()
    fit = root(residual, x0, method='hybr', &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=1.0_dp)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'adaptive transpose cost rejected for non-krylov methods')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=0.0_dp, &
        jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'zero callback cost rejected early')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=-1.0_dp, &
        jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'negative callback cost rejected early')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost= &
            ieee_value(0.0_dp,ieee_quiet_nan), &
        jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'nonfinite callback cost rejected early')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost= &
            ieee_value(0.0_dp,ieee_positive_inf), &
        jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'infinite callback cost rejected early')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, &
        inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true., &
        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=1.0_dp, &
        jacobian_transpose_product=jtv)
    call check(fit%status == 0 .and. fit%nfev == 0, &
        'fixed and adaptive transpose policies are mutually exclusive')

    if (failures > 0) error stop 1
    print '(a)', 'adaptive nonsymmetric transpose-verification tests passed'

contains
    subroutine reset()
        jv_calls = 0
        jtv_calls = 0
        inner_calls = 0
        fail_jtv_at = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(ok, label)
        logical, intent(in) :: ok !! expected condition
        character(len=*), intent(in) :: label !! assertion description
        if (.not. ok) then
            print '(a)', 'FAILED: '//label
            failures = failures+1
        end if
    end subroutine check

    subroutine action(vec, product)
        real(dp), intent(in) :: vec(:) !! direction for the true Jacobian
        real(dp), intent(out) :: product(:) !! resulting Jacobian product
        product = d*vec + u(:,1)*dot_product(v(:,1),vec) + &
            u(:,2)*dot_product(v(:,2),vec)
    end subroutine action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current solution estimate
        real(dp), intent(out) :: fvec(:) !! root residual
        call action(x-target,fvec)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current solution estimate
        real(dp), intent(in) :: vec(:) !! forward multiplication direction
        real(dp), intent(out) :: product(:) !! forward Jacobian product
        logical, intent(out) :: ok !! valid Jacobian evaluation
        jv_calls = jv_calls+1
        call action(vec,product)
        if (borderline_forward) product = product + &
            2.0e-6_dp*forward_probe*dot_product(forward_probe,vec)
        ok = .true.
    end subroutine jvp

    subroutine jtv(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current solution estimate
        real(dp), intent(in) :: vec(:) !! transpose multiplication direction
        real(dp), intent(out) :: product(:) !! transpose Jacobian product
        logical, intent(out) :: ok !! valid transpose evaluation
        jtv_calls = jtv_calls+1
        ok = fail_jtv_at == 0 .or. jtv_calls < fail_jtv_at
        product = 0.0_dp
        if (.not. ok) return
        product = d*vec + v(:,1)*dot_product(u(:,1),vec) + &
            v(:,2)*dot_product(u(:,2),vec)
        if (adversarial) product = product + 0.2_dp*dot_product(spurious,vec)*spurious
    end subroutine jtv

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian linear operation
        procedure(vector_root_linear_operator) :: preconditioner !! inverse-Jacobian action
        real(dp), intent(in) :: rhs(:) !! linear system right-hand side
        real(dp), intent(in) :: rtol !! target linear residual tolerance
        integer, intent(in) :: maxiter !! inner iteration budget
        real(dp), intent(out) :: solution(:) !! Newton step
        logical, intent(out) :: ok !! successful preconditioner action
        inner_calls = inner_calls+1
        call preconditioner(rhs_probe,observed,ok)
        solution = rhs
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_transpose_adaptive
