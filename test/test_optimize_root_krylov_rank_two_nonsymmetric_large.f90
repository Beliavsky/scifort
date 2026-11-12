! SPDX-License-Identifier: MIT
! Large ill-scaled low-product nonsymmetric rank-two discovery and fallback.
program test_optimize_root_krylov_rank_two_nonsymmetric_large
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result, vector_root_linear_operator
    implicit none

    integer, parameter :: n = 128
    real(dp) :: d(n), wrong_d(n), u(n,2), v(n,2), probe(n), rhs_probe(n)
    real(dp) :: x0(n), target(n), diagonal(n), observed(n)
    integer :: i, jv_calls, transpose_calls, solve_calls, failures
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        d(i) = 10.0_dp**(-3.0_dp+6.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,1) = 0.018_dp*sqrt(d(i))*sin(0.37_dp*real(i,dp))
        u(i,2) = 0.016_dp*sqrt(d(i))*cos(0.17_dp*real(i,dp))
        v(i,1) = 0.020_dp*sqrt(d(i))*cos(0.33_dp*real(i,dp))
        v(i,2) = 0.019_dp*sqrt(d(i))*sin(0.23_dp*real(i,dp))
        probe(i) = 0.4_dp*sin(0.21_dp*real(i,dp))+0.13_dp
        target(i) = 0.025_dp*cos(0.13_dp*real(i,dp))
        wrong_d(i) = d(i)*(1.0_dp+0.08_dp*sin(0.53_dp*real(i,dp)))
        diagonal(i) = d(i) + u(i,1)*v(i,1) + u(i,2)*v(i,2)
    end do
    x0 = 0.0_dp
    call action(probe,rhs_probe)

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(jv_calls == 4 .and. transpose_calls == 2 .and. solve_calls == 1, &
        '128-variable nonsymmetric fit needs four forward and two transpose calls')
    call check(maxval(abs(observed-probe)) < 1.0e-9_dp, &
        'accurate Woodbury inverse across six decades of diagonal scaling')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_nonsymmetric_base=wrong_d, jacobian_product=jvp, &
        jacobian_transpose_product=transpose_jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-rhs_probe/diagonal)) < 1.0e-8_dp, &
        'inaccurate supplied diagonal causes exact-diagonal fallback')
    call check(jv_calls >= n+2 .and. jv_calls <= n+4 .and. transpose_calls == 2, &
        'fallback uses n coordinate products without dense Jacobian storage')

    if (failures /= 0) error stop 1
    print '(a)', 'large nonsymmetric rank-two tests passed'

contains
    subroutine reset()
        jv_calls = 0
        transpose_calls = 0
        solve_calls = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected assertion
        character(len=*), intent(in) :: label !! explanatory check name
        if (.not. condition) then
            print '(a)', 'FAILED: '//label
            failures = failures + 1
        end if
    end subroutine check

    subroutine action(vec, result)
        real(dp), intent(in) :: vec(:) !! Jacobian multiplication direction
        real(dp), intent(out) :: result(:) !! exact Jacobian product
        result = d*vec + u(:,1)*dot_product(v(:,1),vec) + &
            u(:,2)*dot_product(v(:,2),vec)
    end subroutine action

    subroutine transpose_action(vec, result)
        real(dp), intent(in) :: vec(:) !! transpose Jacobian direction
        real(dp), intent(out) :: result(:) !! exact transpose action
        result = d*vec + v(:,1)*dot_product(u(:,1),vec) + &
            v(:,2)*dot_product(u(:,2),vec)
    end subroutine transpose_action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        call action(x-target,fvec)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! evaluation point
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! success flag
        jv_calls = jv_calls+1
        call action(vec,product)
        ok = .true.
    end subroutine jvp

    subroutine transpose_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! evaluation point
        real(dp), intent(in) :: vec(:) !! transpose direction
        real(dp), intent(out) :: product(:) !! Jacobian-transpose product
        logical, intent(out) :: ok !! success flag
        transpose_calls = transpose_calls+1
        call transpose_action(vec,product)
        ok = .true.
    end subroutine transpose_jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! nonlinear equation right-hand side
        real(dp), intent(in) :: rtol !! relative inner tolerance
        integer, intent(in) :: maxiter !! Krylov iteration budget
        real(dp), intent(out) :: solution(:) !! Newton correction
        logical, intent(out) :: ok !! successful preconditioner application
        solve_calls = solve_calls+1
        call preconditioner(rhs_probe,observed,ok)
        solution = rhs
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_nonsymmetric_large
