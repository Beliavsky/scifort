! SPDX-License-Identifier: MIT
! Check O(n)-storage rank-two discovery with more unknown columns than pilots.
program test_optimize_root_krylov_rank_two_stream_large
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 128
    real(dp) :: a(n), u(n,2), v(n,2), target(n), probe(n), x0(n)
    real(dp) :: observed(n), expected(n), diagonal(n), w(n,2)
    real(dp) :: small(2,2), rhs(2), det, coeff(2), t
    integer :: i, failures, jv_calls, inner_calls
    logical :: contaminate, rank_one_only
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        t = real(i,dp)
        a(i) = exp(-5.0_dp + 10.0_dp*(t-1.0_dp)/real(n-1,dp))
        u(i,1) = 0.03_dp*a(i)*sin(0.41_dp*t)
        u(i,2) = 0.02_dp*a(i)*cos(0.29_dp*t)
        v(i,1) = cos(0.17_dp*t)
        v(i,2) = sin(0.23_dp*t)
        target(i) = 0.05_dp*cos(0.11_dp*t)
        probe(i) = sin(0.35_dp*t)
    end do
    x0 = 0.0_dp
    contaminate = .false.
    rank_one_only = .false.
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_stream=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=2, fatol=1.0e-10_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 1.0e-10_dp, &
        'n=128 scaled rank-two nonlinear solve converges')
    call check(maxval(abs(observed-expected)/max(1.0_dp,abs(expected))) < 1.0e-9_dp, &
        'n=128 Woodbury action agrees with analytic inverse')
    call check(jv_calls == n .and. inner_calls == 1, &
        'streamed verified discovery needs n products, even with eight pilot columns')

    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_stream=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=2, f_tol=1.0e-10_dp)
    call check(fit%success .and. jv_calls == n, &
        'newton_krylov forwards streamed discovery option')

    ! Column 33 lies outside the eight evenly spaced pilot columns.
    contaminate = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_stream=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 1.0e-10_dp, &
        'off-pilot rank-three contamination rejects streamed model')
    call check(jv_calls == n, 'rejected stream retains exact measured diagonal')
    contaminate = .false.

    rank_one_only = .true.
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_stream=.true., jacobian_product=jvp, &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 1.0e-10_dp, &
        'rank-one input remains diagonal fallback')
    call check(jv_calls == n, 'rank-deficient pilot still evaluates all n exact diagonals')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_two_stream_large: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_two_stream_large: PASS'

contains
    subroutine reset()
        jv_calls = 0
        inner_calls = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(yes, description)
        logical, intent(in) :: yes !! test property
        character(len=*), intent(in) :: description !! checked condition
        if (.not. yes) then
            failures = failures+1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check

    subroutine reference()
        real(dp) :: u_now(n,2)
        u_now = u
        if (rank_one_only) u_now(:,2) = 0.0_dp
        diagonal = a + u_now(:,1)*v(:,1) + u_now(:,2)*v(:,2)
        w(:,1) = u_now(:,1)/a
        w(:,2) = u_now(:,2)/a
        small = matmul(transpose(v),w)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        det = small(1,1)*small(2,2)-small(1,2)*small(2,1)
        rhs = matmul(transpose(v),probe/a)
        coeff(1) = (small(2,2)*rhs(1)-small(1,2)*rhs(2))/det
        coeff(2) = (small(1,1)*rhs(2)-small(2,1)*rhs(1))/det
        expected = probe/a - matmul(w,coeff)
    end subroutine reference

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: f(:) !! residual
        real(dp) :: delta(n), factor(n,2)
        delta = x-target
        factor = u
        if (rank_one_only) factor(:,2) = 0.0_dp
        f = a*delta + matmul(factor,matmul(transpose(v),delta))
        if (contaminate) f(3) = f(3)+0.17_dp*delta(33)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(in) :: direction(:) !! Jacobian probe
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! true for a finite analytic product
        real(dp) :: factor(n,2)
        jv_calls = jv_calls+1
        factor = u
        if (rank_one_only) factor(:,2) = 0.0_dp
        product = a*direction + matmul(factor,matmul(transpose(v),direction))
        if (contaminate) product(3) = product(3)+0.17_dp*direction(33)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, b, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: b(:) !! Newton right hand side
        real(dp), intent(in) :: rtol !! requested tolerance
        integer, intent(in) :: maxiter !! allowed iteration count
        real(dp), intent(out) :: solution(:) !! preconditioned right hand side
        logical, intent(out) :: ok !! preconditioner validity
        inner_calls = inner_calls+1
        call preconditioner(probe,observed,ok)
        if (.not. ok) return
        call preconditioner(b,solution,ok)
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_stream_large
