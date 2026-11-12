! SPDX-License-Identifier: MIT
! Learn both rank-two coupling factors using exact coordinate Jv products.
program test_optimize_root_krylov_rank_two_discover
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 12
    real(dp) :: a(n), u(n,2), v(n,2), target(n), probe(n), observed(n), x0(n)
    real(dp) :: expected(n), diagonal(n), original_a(n), original_u(n,2)
    real(dp) :: work(n,2), small(2,2), projection(2), coeff(2), expected_short(n)
    integer :: i, failures, f_calls, jv_calls, solver_calls, fail_after, short_n
    logical :: incompatible, rank_one_only, sparse_coupling, singular_model
    type(vector_root_result) :: fit

    failures = 0
    do i = 1, n
        a(i) = 0.8_dp + 0.16_dp*real(i,dp)
        u(i,1) = 0.1_dp + 0.014_dp*real(i,dp)
        u(i,2) = 0.13_dp*sin(0.67_dp*real(i,dp))
        v(i,1) = cos(0.43_dp*real(i,dp))
        v(i,2) = sin(0.59_dp*real(i,dp))
        target(i) = 0.06_dp*sin(0.6_dp*real(i,dp))
        probe(i) = 0.1_dp + 0.039_dp*real(i,dp)
    end do
    v(4,1) = 0.0_dp
    v(9,2) = 0.0_dp
    x0 = 0.0_dp
    incompatible = .false.
    rank_one_only = .false.
    sparse_coupling = .false.
    singular_model = .false.
    call reference()

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, &
        line_search='none', maxiter=2, fatol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 4.0e-12_dp, &
        'learned rank-two inverse solves coupled linear roots')
    call check(maxval(abs(observed-expected)) < 5.0e-12_dp, &
        'learned rank-two inverse action matches analytic Woodbury')
    call check(jv_calls == n .and. solver_calls == 1 .and. fit%nfev == 2, &
        'coordinate products reused without additional verification probes')

    call reset()
    fit = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, &
        line_search='none', maxiter=2, f_tol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(observed-expected)) < 5.0e-12_dp, &
        'newton_krylov forwards unknown rank-two discovery')

    call reset()
    fit = root(residual, x0, method='krylov', inner_method='gmres', &
        inner_maxiter=1, inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., jacobian_product=jvp, &
        line_search='none', maxiter=3, fatol=1.0e-11_dp)
    call check(fit%success .and. maxval(abs(fit%x-target)) < 5.0e-11_dp, &
        'one-dimensional GMRES using the discovered exact inverse')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 5.0e-6_dp, &
        'forward finite differences identify both right factors')
    call check(fit%nfev == n+2 .and. fit%nfev == f_calls, &
        'two-point discovery uses exactly n finite-difference products')

    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_solver=inspect, jacobian_diff_scheme='3-point', &
        line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)) < 5.0e-6_dp, &
        'centered finite differences identify both right factors')
    call check(fit%nfev == 2*n+2 .and. fit%nfev == f_calls, &
        'centered discovery requires two n coordinate perturbations')

    original_a = a
    original_u = u
    do i = 1, n
        a(i) = 10.0_dp**(-4.0_dp + 8.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,1) = 0.04_dp*a(i)
        u(i,2) = 0.05_dp*a(i)*sin(0.37_dp*real(i,dp))
    end do
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-expected)/max(1.0_dp,abs(expected))) < 3.0e-9_dp, &
        'discovery retains unfloored diagonal over eight decades')
    a = original_a
    u = original_u
    call reference()

    incompatible = .true.
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 5.0e-13_dp, &
        'rank-three contamination rejects the rank-two fit')
    call check(jv_calls == n .and. solver_calls == 1, &
        'rejected fit still reuses all coordinate products')
    incompatible = .false.

    rank_one_only = .true.
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 5.0e-13_dp, &
        'unidentifiable rank-one coupling keeps diagonal scaling')
    rank_one_only = .false.
    call reference()

    sparse_coupling = .true.
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 5.0e-13_dp, &
        'sparse rank-two coupling cannot overfit missing diagonals')
    sparse_coupling = .false.
    call reference()

    singular_model = .true.
    call reference()
    call reset()
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect, line_search='none', maxiter=1)
    call check(maxval(abs(observed-probe/diagonal)) < 5.0e-12_dp, &
        'rank-two singular Woodbury factor reverts to diagonal scaling')
    singular_model = .false.
    call reference()

    short_n = 6
    call reference_short(short_n)
    call reset()
    fit = root(residual_short, x0(:short_n), method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., jacobian_product=jvp_short, &
        inner_solver=inspect_short, line_search='none', maxiter=1)
    call check(maxval(abs(observed(:short_n)-expected_short(:short_n))) < 4.0e-12_dp, &
        'smallest identifiable rank-two model uses Woodbury')
    call check(jv_calls == short_n, 'six-variable discovery uses exact coordinate probes')

    short_n = 5
    call reset()
    fit = root(residual_short, x0(:short_n), method='krylov', &
        inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
        inner_m_diag_rank_two_discover=.true., jacobian_product=jvp_short, &
        inner_solver=inspect_short, line_search='none', maxiter=1)
    call check(maxval(abs(observed(:short_n)-probe(:short_n)/diagonal(:short_n))) < 5.0e-13_dp, &
        'nonidentifiable five-variable case reverts to exact diagonal')
    call check(jv_calls == short_n, 'short nonidentifiable case adds no probes')

    fit = root(residual, x0, method='krylov', inner_m_diag_rank_two_discover=.true.)
    call check(.not. fit%success .and. fit%nfev == 0 .and. &
        index(fit%message,'requires inner_m_diag_rank_two') > 0, &
        'discovery requires rank-two option')
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_right=v)
    call check(.not. fit%success .and. fit%nfev == 0 .and. &
        index(fit%message,'conflict') > 0, &
        'known and discovered rank-two right factors conflict')
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_rank_two_sketch=.true.)
    call check(.not. fit%success .and. fit%nfev == 0 .and. &
        index(fit%message,'conflict') > 0, 'discovery excludes the sketch option')
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        inner_m_diag_probes=2)
    call check(.not. fit%success .and. fit%nfev == 0, &
        'discovery rejects inexact diagonal probes')
    fit = root(residual, x0, method='hybr', inner_m_diag_rank_two_discover=.true.)
    call check(.not. fit%success .and. fit%nfev == 0, &
        'non-Krylov methods reject discovery')

    call reset()
    fail_after = 5
    fit = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. fit%success .and. solver_calls == 0 .and. jv_calls == 5, &
        'callback failure stops before the inner Krylov solve')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_two_discover: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_two_discover: PASS'

contains

    subroutine reset()
        f_calls = 0
        jv_calls = 0
        solver_calls = 0
        fail_after = 0
        observed = 0.0_dp
    end subroutine reset

    subroutine check(condition, title)
        logical, intent(in) :: condition !! property verified by this test
        character(len=*), intent(in) :: title !! description of the property
        if (.not. condition) then
            failures = failures+1
            print '(a)', 'FAIL: '//title
        end if
    end subroutine check

    subroutine model_factors(coupling)
        real(dp), intent(out) :: coupling(n,2) !! current left factor
        coupling = u
        if (rank_one_only) coupling(:,2) = 0.0_dp
        if (sparse_coupling) then
            coupling = 0.0_dp
            coupling(1,1) = 0.45_dp
            coupling(5,2) = 0.38_dp
        end if
        if (singular_model) then
            coupling(:,1) = -a*v(:,1)/dot_product(v(:,1),v(:,1))
            coupling(:,2) = 0.1_dp*a*(v(:,2) - &
                dot_product(v(:,1),v(:,2))/dot_product(v(:,1),v(:,1))*v(:,1))
        end if
    end subroutine model_factors

    subroutine reference()
        real(dp) :: factor(n,2), determinant
        call model_factors(factor)
        diagonal = a + factor(:,1)*v(:,1) + factor(:,2)*v(:,2)
        if (singular_model) return
        work(:,1) = factor(:,1)/a
        work(:,2) = factor(:,2)/a
        small = matmul(transpose(v),work)
        small(1,1) = small(1,1)+1.0_dp
        small(2,2) = small(2,2)+1.0_dp
        determinant = small(1,1)*small(2,2)-small(1,2)*small(2,1)
        projection = matmul(transpose(v),probe/a)
        coeff(1) = (small(2,2)*projection(1)-small(1,2)*projection(2))/determinant
        coeff(2) = (small(1,1)*projection(2)-small(2,1)*projection(1))/determinant
        expected = probe/a - matmul(work,coeff)
    end subroutine reference

    subroutine reference_short(m)
        integer, intent(in) :: m !! dimension of reduced rank-two problem
        real(dp) :: w(m,2), h(2,2), q(2), coeffs(2), det
        w(:,1) = u(:m,1)/a(:m)
        w(:,2) = u(:m,2)/a(:m)
        h = matmul(transpose(v(:m,:)),w)
        h(1,1) = h(1,1)+1.0_dp
        h(2,2) = h(2,2)+1.0_dp
        det = h(1,1)*h(2,2)-h(1,2)*h(2,1)
        q = matmul(transpose(v(:m,:)),probe(:m)/a(:m))
        coeffs(1) = (h(2,2)*q(1)-h(1,2)*q(2))/det
        coeffs(2) = (h(1,1)*q(2)-h(2,1)*q(1))/det
        expected_short(:m) = probe(:m)/a(:m)-matmul(w,coeffs)
    end subroutine reference_short

    subroutine residual_short(x, f)
        real(dp), intent(in) :: x(:) !! reduced system iterate
        real(dp), intent(out) :: f(:) !! reduced system residual
        real(dp) :: delta(size(x))
        f_calls = f_calls+1
        delta = x-target(:size(x))
        f = a(:size(x))*delta + matmul(u(:size(x),:), &
            matmul(transpose(v(:size(x),:)),delta))
    end subroutine residual_short

    subroutine jvp_short(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! reduced system iterate
        real(dp), intent(in) :: direction(:) !! reduced Jacobian direction
        real(dp), intent(out) :: product(:) !! reduced Jacobian action
        logical, intent(out) :: ok !! callback success
        jv_calls = jv_calls+1
        product = a(:size(x))*direction + matmul(u(:size(x),:), &
            matmul(transpose(v(:size(x),:)),direction))
        ok = .true.
    end subroutine jvp_short

    subroutine inspect_short(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse action
        real(dp), intent(in) :: rhs(:) !! reduced right hand side
        real(dp), intent(in) :: rtol !! requested relative tolerance
        integer, intent(in) :: maxiter !! available basis count
        real(dp), intent(out) :: solution(:) !! reduced Newton direction
        logical, intent(out) :: ok !! action validity
        solver_calls = solver_calls+1
        call preconditioner(probe(:short_n),observed(:short_n),ok)
        if (.not. ok) return
        call preconditioner(rhs,solution,ok)
    end subroutine inspect_short

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! current nonlinear variable vector
        real(dp), intent(out) :: f(:) !! residuals for coupled linear system
        real(dp) :: delta(n), factor(n,2)
        f_calls = f_calls+1
        call model_factors(factor)
        delta = x-target
        f = a*delta + matmul(factor,matmul(transpose(v),delta))
        if (incompatible) f(2) = f(2) + 0.3_dp*delta(7)
    end subroutine residual

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current point (constant Jacobian)
        real(dp), intent(in) :: direction(:) !! trial Jacobian direction
        real(dp), intent(out) :: product(:) !! applied Jacobian result
        logical, intent(out) :: ok !! false when failure is injected
        real(dp) :: factor(n,2)
        jv_calls = jv_calls+1
        ok = .false.
        if (fail_after > 0 .and. jv_calls >= fail_after) return
        call model_factors(factor)
        product = a*direction + matmul(factor,matmul(transpose(v),direction))
        if (incompatible) product(2) = product(2) + 0.3_dp*direction(7)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! current Jv
        procedure(vector_root_linear_operator) :: preconditioner !! current inverse
        real(dp), intent(in) :: rhs(:) !! right hand side
        real(dp), intent(in) :: rtol !! inner solver relative tolerance
        integer, intent(in) :: maxiter !! iteration budget
        real(dp), intent(out) :: solution(:) !! approximate inverse action
        logical, intent(out) :: ok !! whether inverse application succeeded
        solver_calls = solver_calls+1
        call preconditioner(probe, observed, ok)
        if (.not. ok) return
        call preconditioner(rhs, solution, ok)
    end subroutine inspect
end program test_optimize_root_krylov_rank_two_discover
