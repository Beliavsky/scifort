! SPDX-License-Identifier: MIT
! Size-independent setup work for a known diagonal plus symmetric rank-two Jacobian.
program test_optimize_root_krylov_rank_two_symmetric_large
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result, vector_root_linear_operator
    implicit none

    integer, parameter :: n = 128
    real(dp) :: d(n), wrong_d(n), u(n,2), probe(n), rhs_probe(n), x0(n), target(n)
    real(dp) :: diagonal(n), observed(n)
    integer :: i, jv_calls, solve_calls, failures
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        d(i) = 10.0_dp**(-3.0_dp+6.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,1) = 0.016_dp*sqrt(d(i))*sin(0.37_dp*real(i,dp))
        u(i,2) = 0.019_dp*sqrt(d(i))*cos(0.17_dp*real(i,dp))
        probe(i) = 0.4_dp*sin(0.21_dp*real(i,dp))+0.13_dp
        target(i) = 0.025_dp*cos(0.13_dp*real(i,dp))
        wrong_d(i) = d(i)*(1.0_dp+0.08_dp*sin(0.53_dp*real(i,dp)))
        diagonal(i) = d(i) + u(i,1)**2 - 0.7_dp*u(i,2)**2
    end do
    x0 = 0.0_dp
    call action(probe,rhs_probe)

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=d, jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 4 .and. solve_calls == 1, &
        '128-variable symmetric discovery needs four Jacobian products')
    call check(maxval(abs(observed-probe)) < 5.0e-10_dp, &
        'accurate rank-two inverse across six decades of base scaling')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_symmetric_base=wrong_d, jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-rhs_probe/diagonal)) < 1.0e-8_dp, &
        'incorrect supplied base diagonal triggers safe exact-diagonal fallback')
    call check(jv_calls >= n+2 .and. jv_calls <= n+4, &
        'fallback performs the full coordinate diagonal evaluation')

    if (failures /= 0) error stop 1
    print '(a)', 'large symmetric rank-two tests passed'

contains
    subroutine reset()
        jv_calls = 0
        solve_calls = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, label)
        logical, intent(in) :: condition !! passing assertion
        character(len=*), intent(in) :: label !! assertion description
        if (.not. condition) then
            print '(a)', 'FAILED: '//label
            failures = failures + 1
        end if
    end subroutine check

    subroutine action(vec, result)
        real(dp), intent(in) :: vec(:) !! Jacobian multiplication direction
        real(dp), intent(out) :: result(:) !! exact Jacobian product
        result = d*vec + u(:,1)*dot_product(u(:,1),vec) - &
            0.7_dp*u(:,2)*dot_product(u(:,2),vec)
    end subroutine action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear evaluation point
        real(dp), intent(out) :: fvec(:) !! linear residual at x
        call action(x-target,fvec)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! Jacobian evaluation point
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! exact Jacobian action
        logical, intent(out) :: ok !! true for a successful action
        jv_calls = jv_calls+1
        call action(vec,product)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian operator
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! nonlinear system right side
        real(dp), intent(in) :: rtol !! requested tolerance
        integer, intent(in) :: maxiter !! linear iteration budget
        real(dp), intent(out) :: solution(:) !! Newton solution estimate
        logical, intent(out) :: ok !! successful application
        solve_calls = solve_calls+1
        call preconditioner(rhs_probe,observed,ok)
        solution = rhs
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_symmetric_large
