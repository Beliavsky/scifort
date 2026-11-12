program test_optimize_bracket
    use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_minimize_scalar, only : bracket, bracket_result, &
        minimize_status_bracket_error, minimize_status_convergence_error, &
        minimize_status_invalid_input, minimize_status_success, &
        minimize_status_value_error
    use test_optimize_bracket_reference
    implicit none

    integer :: failures
    type(bracket_result) :: r

    failures = 0
    r = bracket(quadratic)
    call check_result('default', r, [default_xa_ref, default_xb_ref, default_xc_ref], &
        [default_fa_ref, default_fb_ref, default_fc_ref], default_calls_ref, failures)

    r = bracket(quadratic, xa=3.0_dp, xb=4.0_dp)
    call check_result('reverse', r, [reverse_xa_ref, reverse_xb_ref, reverse_xc_ref], &
        [reverse_fa_ref, reverse_fb_ref, reverse_fc_ref], reverse_calls_ref, failures)
    call check_true('descending bracket', r%xa > r%xb .and. r%xb > r%xc, failures)

    r = bracket(quadratic, xa=-2.0_dp, xb=-1.0_dp)
    call check_result('uphill', r, [uphill_xa_ref, uphill_xb_ref, uphill_xc_ref], &
        [uphill_fa_ref, uphill_fb_ref, uphill_fc_ref], uphill_calls_ref, failures)

    r = bracket(polynomial)
    call check_result('polynomial', r, [polynomial_xa_ref, polynomial_xb_ref, polynomial_xc_ref], &
        [polynomial_fa_ref, polynomial_fb_ref, polynomial_fc_ref], polynomial_calls_ref, failures)

    r = bracket(far_quadratic)
    call check_result('far default', r, &
        [far_default_xa_ref, far_default_xb_ref, far_default_xc_ref], &
        [far_default_fa_ref, far_default_fb_ref, far_default_fc_ref], far_default_calls_ref, failures)

    r = bracket(far_quadratic, grow_limit=2.0_dp)
    call check_result('far limited', r, &
        [far_limited_xa_ref, far_limited_xb_ref, far_limited_xc_ref], &
        [far_limited_fa_ref, far_limited_fb_ref, far_limited_fc_ref], far_limited_calls_ref, failures)
    call check_true('grow_limit changes work', r%function_calls > far_default_calls_ref, failures)

    r = bracket(constant)
    call check_true('flat objective rejects bracket', .not. r%success, failures)
    call check_int('flat status', r%status, minimize_status_bracket_error, failures)
    call check_int('flat evaluations', r%function_calls, 3, failures)

    r = bracket(nan_objective)
    call check_true('NaN objective rejects bracket', .not. r%success, failures)
    call check_int('NaN status', r%status, minimize_status_value_error, failures)
    call check_int('NaN evaluations', r%function_calls, 3, failures)

    r = bracket(linear_decreasing, maxiter=0)
    call check_true('iteration cap rejects monotonic objective', .not. r%success, failures)
    call check_int('iteration cap status', r%status, minimize_status_convergence_error, failures)
    call check_true('iteration cap work', r%function_calls >= 3, failures)

    r = bracket(quadratic, grow_limit=-1.0_dp)
    call check_result('negative grow limit', r, [default_xa_ref, default_xb_ref, default_xc_ref], &
        [default_fa_ref, default_fb_ref, default_fc_ref], default_calls_ref, failures)

    r = bracket(quadratic, grow_limit=ieee_value(0.0_dp, ieee_positive_inf))
    call check_int('invalid growth status', r%status, minimize_status_invalid_input, failures)
    call check_int('invalid growth makes no evaluations', r%function_calls, 0, failures)

    r = bracket(quadratic, maxiter=-1)
    call check_int('invalid iteration cap', r%status, minimize_status_invalid_input, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_bracket: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_bracket: PASS'

contains

    function quadratic(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y
        y = (x - 2.0_dp)**2
    end function quadratic

    function far_quadratic(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y
        y = (x - 50.0_dp)**2
    end function far_quadratic

    function polynomial(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y
        y = (x - 2.0_dp) * x * (x + 2.0_dp)**2
    end function polynomial

    function constant(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y
        y = 1.0_dp + 0.0_dp * x
    end function constant

    function nan_objective(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y
        y = ieee_value(x, ieee_quiet_nan)
    end function nan_objective

    function linear_decreasing(x) result(y)
        real(dp), intent(in) :: x !! input
        real(dp) :: y
        y = -x
    end function linear_decreasing

    subroutine check_result(label, actual, expected_x, expected_f, expected_calls, nfail)
        character(len=*), intent(in) :: label !! case name
        type(bracket_result), intent(in) :: actual !! calculated bracket
        real(dp), intent(in) :: expected_x(3) !! SciPy reference abscissae
        real(dp), intent(in) :: expected_f(3) !! SciPy reference ordinates
        integer, intent(in) :: expected_calls !! SciPy reference function calls
        integer, intent(inout) :: nfail !! failure count
        real(dp), parameter :: xtol = 1.0e-8_dp
        real(dp), parameter :: ftol = 1.0e-8_dp

        call check_true(trim(label)//' success', actual%success, nfail)
        call check_int(trim(label)//' status', actual%status, minimize_status_success, nfail)
        call check_int(trim(label)//' calls', actual%function_calls, expected_calls, nfail)
        call check_vector(trim(label)//' x', [actual%xa, actual%xb, actual%xc], expected_x, xtol, nfail)
        call check_vector(trim(label)//' f', [actual%fa, actual%fb, actual%fc], expected_f, ftol, nfail)
        call check_true(trim(label)//' downhill', &
            (actual%fb < actual%fa .and. actual%fb <= actual%fc) .or. &
            (actual%fb <= actual%fa .and. actual%fb < actual%fc), nfail)
        call check_true(trim(label)//' ordered', &
            (actual%xa < actual%xb .and. actual%xb < actual%xc) .or. &
            (actual%xc < actual%xb .and. actual%xb < actual%xa), nfail)
    end subroutine check_result

    subroutine check_vector(label, actual, expected, tol, nfail)
        character(len=*), intent(in) :: label !! case name
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tol !! tolerance
        integer, intent(inout) :: nfail !! failure count
        if (any(abs(actual - expected) > tol * max(1.0_dp, abs(expected)))) then
            write(*, '(a,3(1x,es15.7))') trim(label)//' actual:', actual
            write(*, '(a,3(1x,es15.7))') trim(label)//' reference:', expected
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_int(label, actual, expected, nfail)
        character(len=*), intent(in) :: label !! case name
        integer, intent(in) :: actual !! computed value
        integer, intent(in) :: expected !! expected value
        integer, intent(inout) :: nfail !! failure count
        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(label)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_true(label, condition, nfail)
        character(len=*), intent(in) :: label !! case name
        logical, intent(in) :: condition !! assertion
        integer, intent(inout) :: nfail !! failure count
        if (.not. condition) then
            write(*, '(a)') trim(label)//': FAIL'
            nfail = nfail + 1
        end if
    end subroutine check_true
end program test_optimize_bracket
