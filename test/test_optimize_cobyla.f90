program test_optimize_cobyla
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : cobyla, minimize, multivariate_minimize_result
    use test_optimize_cobyla_reference
    implicit none

    type(multivariate_minimize_result) :: r
    real(dp) :: inf
    real(dp) :: simple_x0(2)
    integer :: failures

    failures = 0
    inf = ieee_value(0.0_dp, ieee_positive_inf)
    simple_x0 = [sqrt(25.0_dp - (2.0_dp / 3.0_dp)**2), &
        2.0_dp / 3.0_dp + 1.0e-4_dp]

    r = cobyla(canonical_objective, [2.0_dp, 0.0_dp], &
        constraints=canonical_constraints, mineq=3)
    call check_result('canonical inequalities', r, canonical_x_ref, canonical_fun_ref, &
        canonical_maxcv_ref, canonical_status_ref, canonical_success_ref, &
        2.0e-5_dp, 2.0e-10_dp, 1.0e-12_dp, canonical_nfev_ref, 2, failures)

    r = minimize(canonical_objective, [2.0_dp, 0.0_dp], method='cobyla', &
        constraints=canonical_constraints, mineq=3)
    call check_result('minimize dispatch', r, canonical_x_ref, canonical_fun_ref, &
        canonical_maxcv_ref, canonical_status_ref, canonical_success_ref, &
        2.0e-5_dp, 2.0e-10_dp, 1.0e-12_dp, canonical_nfev_ref, 2, failures)
    if (trim(r%method) /= 'cobyla') then
        write(*, '(a,a)') 'minimize dispatch method mismatch: ', trim(r%method)
        failures = failures + 1
    end if

    r = cobyla(simple_objective, simple_x0, constraints=circle_pair, mineq=2, &
        rhobeg=1.0_dp, tol=1.0e-6_dp, maxiter=100)
    call check_result('circle pair', r, simple_x_ref, simple_fun_ref, simple_maxcv_ref, &
        simple_status_ref, simple_success_ref, 1.0e-4_dp, 2.0e-9_dp, 1.0e-8_dp, &
        simple_nfev_ref, 3, failures)

    r = cobyla(simple_objective, simple_x0, constraints=sum_equality, meq=1, &
        rhobeg=1.0_dp, tol=1.0e-6_dp, maxiter=100)
    call check_result('equality conversion', r, equality_x_ref, equality_fun_ref, &
        equality_maxcv_ref, equality_status_ref, equality_success_ref, &
        1.0e-4_dp, 2.0e-9_dp, 1.0e-10_dp, equality_nfev_ref, 8, failures)

    r = cobyla(squared_norm, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
        lower=[-1.0_dp, -inf, 1.0_dp, -inf, -0.5_dp], &
        upper=[-0.5_dp, -0.5_dp, inf, inf, -0.5_dp])
    call check_result('mixed and fixed bounds', r, bounds_x_ref, bounds_fun_ref, &
        bounds_maxcv_ref, bounds_status_ref, bounds_success_ref, &
        1.0e-3_dp, 2.0e-8_dp, 1.0e-12_dp, bounds_nfev_ref, 5, failures)

    r = cobyla(squared_norm, [2.0_dp, 3.0_dp], &
        lower=[0.5_dp, 0.5_dp], upper=[0.5_dp, 0.5_dp], &
        constraints=sum_equality, meq=1)
    call check_result('all fixed feasible', r, [0.5_dp, 0.5_dp], 0.5_dp, &
        0.0_dp, 0, .true., 1.0e-14_dp, 1.0e-14_dp, 1.0e-14_dp, 1, 0, failures)

    r = cobyla(squared_norm, [2.0_dp, 3.0_dp], &
        lower=[0.0_dp, 0.0_dp], upper=[0.0_dp, 0.0_dp], &
        constraints=sum_equality, meq=1)
    call check_result('all fixed infeasible', r, [0.0_dp, 0.0_dp], 0.0_dp, &
        1.0_dp, 0, .false., 1.0e-14_dp, 1.0e-14_dp, 1.0e-14_dp, 1, 0, failures)

    r = cobyla(squared_norm, [1.0_dp, 2.0_dp], lower=[1.0_dp, -inf], &
        upper=[inf, inf])
    call check_result('lower-only bounds', r, lower_only_x_ref, lower_only_fun_ref, &
        lower_only_maxcv_ref, lower_only_status_ref, lower_only_success_ref, &
        1.0e-4_dp, 2.0e-9_dp, 1.0e-12_dp, lower_only_nfev_ref, 3, failures)

    r = cobyla(linear_objective, [0.0_dp], constraints=conflicting_constraints, &
        mineq=2, catol=0.6_dp)
    call check_result('conflicting constraints loose catol', r, conflict_loose_x_ref, &
        conflict_loose_fun_ref, conflict_loose_maxcv_ref, conflict_loose_status_ref, &
        conflict_loose_success_ref, 1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp, &
        conflict_loose_nfev_ref, 0, failures)

    r = cobyla(linear_objective, [0.0_dp], constraints=conflicting_constraints, &
        mineq=2, catol=0.4_dp)
    call check_result('conflicting constraints tight catol', r, conflict_tight_x_ref, &
        conflict_tight_fun_ref, conflict_tight_maxcv_ref, conflict_tight_status_ref, &
        conflict_tight_success_ref, 1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp, &
        conflict_tight_nfev_ref, 0, failures)

    r = cobyla(square_objective, [500.0_dp], f_target=250.0_dp)
    call check_result('objective target', r, target_x_ref, target_fun_ref, &
        target_maxcv_ref, target_status_ref, target_success_ref, &
        1.0e-9_dp, 1.0e-8_dp, 1.0e-12_dp, target_nfev_ref, 1, failures)

    r = cobyla(canonical_objective, [2.0_dp, 0.0_dp], &
        constraints=canonical_constraints, mineq=3, maxiter=4)
    call check_result('evaluation limit', r, maxfun_x_ref, maxfun_fun_ref, &
        maxfun_maxcv_ref, maxfun_status_ref, maxfun_success_ref, &
        1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp, maxfun_nfev_ref, 0, failures)

    r = cobyla(squared_norm, [0.0_dp, 0.0_dp], &
        lower=[1.0_dp, 0.0_dp], upper=[0.0_dp, 1.0_dp])
    call check_int('infeasible bounds status', r%status, -1, failures)
    if (r%success) then
        write(*, '(a)') 'infeasible bounds unexpectedly succeeded'
        failures = failures + 1
    end if

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_cobyla: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_cobyla: PASS'

contains

    function canonical_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = (x(1) - 1.0_dp)**2 + (x(2) - 2.5_dp)**2
    end function canonical_objective

    subroutine canonical_constraints(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! feasible when all entries are nonnegative

        c(1) = x(1) - 2.0_dp * x(2) + 2.0_dp
        c(2) = -x(1) - 2.0_dp * x(2) + 6.0_dp
        c(3) = -x(1) + 2.0_dp * x(2) + 2.0_dp
    end subroutine canonical_constraints

    function simple_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = x(1)**2 + abs(x(2))**3
    end function simple_objective

    subroutine circle_pair(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! paired inequalities enforcing a circle
        real(dp) :: h

        h = x(1)**2 + x(2)**2 - 25.0_dp
        c = [h, -h]
    end subroutine circle_pair

    subroutine sum_equality(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! equality constraint value

        c(1) = x(1) + x(2) - 1.0_dp
    end subroutine sum_equality

    function squared_norm(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = dot_product(x, x)
    end function squared_norm

    function linear_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = x(1)
    end function linear_objective

    subroutine conflicting_constraints(x, c)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp), intent(out) :: c(:) !! incompatible inequalities

        c = [4.0_dp - x(1), x(1) - 5.0_dp]
    end subroutine conflicting_constraints

    function square_objective(x) result(fx)
        real(dp), intent(in) :: x(:) !! parameter vector
        real(dp) :: fx

        fx = x(1)**2
    end function square_objective

    subroutine check_result(name, actual, expected_x, expected_fun, expected_maxcv, &
            expected_status, expected_success, xtol, ftol, cvtol, expected_nfev, &
            nfev_slack, nfail)
        character(len=*), intent(in) :: name !! check label
        type(multivariate_minimize_result), intent(in) :: actual !! SciFort result
        real(dp), intent(in) :: expected_x(:) !! SciPy reference point
        real(dp), intent(in) :: expected_fun !! SciPy reference objective
        real(dp), intent(in) :: expected_maxcv !! SciPy reference maximum violation
        integer, intent(in) :: expected_status !! SciPy reference status
        logical, intent(in) :: expected_success !! SciPy reference success flag
        real(dp), intent(in) :: xtol !! point tolerance
        real(dp), intent(in) :: ftol !! objective tolerance
        real(dp), intent(in) :: cvtol !! constraint-violation tolerance
        integer, intent(in) :: expected_nfev !! SciPy reference objective evaluations
        integer, intent(in) :: nfev_slack !! allowed platform work-count difference
        integer, intent(inout) :: nfail !! accumulated failures

        if (size(actual%x) /= size(expected_x) .or. &
                maxval(abs(actual%x - expected_x)) > xtol) then
            write(*, '(a,a)') trim(name), ': x mismatch'
            if (size(actual%x) == size(expected_x)) then
                write(*, '(a,*(es16.8,1x))') '  actual:   ', actual%x
                write(*, '(a,*(es16.8,1x))') '  expected: ', expected_x
            end if
            nfail = nfail + 1
        end if
        if (abs(actual%fun - expected_fun) > ftol) then
            write(*, '(a,a,2es16.8)') trim(name), ': objective mismatch ', &
                actual%fun, expected_fun
            nfail = nfail + 1
        end if
        if (abs(actual%maxcv - expected_maxcv) > cvtol) then
            write(*, '(a,a,2es16.8)') trim(name), ': maxcv mismatch ', &
                actual%maxcv, expected_maxcv
            nfail = nfail + 1
        end if
        call check_int(trim(name) // ' status', actual%status, expected_status, nfail)
        if (actual%success .neqv. expected_success) then
            write(*, '(a,a,2l2)') trim(name), ': success mismatch ', &
                actual%success, expected_success
            nfail = nfail + 1
        end if
        if (abs(actual%nfev - expected_nfev) > nfev_slack) then
            write(*, '(a,a,2i8)') trim(name), ': nfev mismatch ', &
                actual%nfev, expected_nfev
            nfail = nfail + 1
        end if
    end subroutine check_result

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! actual integer
        integer, intent(in) :: expected !! expected integer
        integer, intent(inout) :: nfail !! accumulated failures

        if (actual /= expected) then
            write(*, '(a,a,2i8)') trim(name), ': mismatch ', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

end program test_optimize_cobyla
