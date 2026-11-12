! SPDX-License-Identifier: MIT
! Regression for SciPy-style custom termination norms in inexact-Newton solvers.
! Numerical work counts come from SciPy 1.17.0 scipy.optimize.root.
program test_optimize_root_tol_norm
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize, only : anderson, broyden1, broyden2, diagbroyden, &
        excitingmixing, linearmixing, newton_krylov, root, vector_root_result
    implicit none

    real(dp), parameter :: x0(2) = [0.0_dp, 0.0_dp]
    integer :: failures, calls
    type(vector_root_result) :: r

    failures = 0
    calls = 0

    r = root(affine, x0, method='broyden1', fatol=0.08_dp)
    call check(r%success .and. r%nit == 3 .and. r%nfev == 8, &
        'default max norm remains SciPy compatible')

    r = root(affine, x0, method='broyden1', fatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted Broyden-1 convergence and work counts')
    call check(abs(r%fun(1) + 0.6975185951049945_dp) < 3.0e-8_dp, &
        'weighted norm changes Broyden stopping location')
    call check(abs(r%fun(2) + 0.06975185951049945_dp) < 3.0e-8_dp, &
        'weighted norm preserves full residual')
    calls = 0
    r = root(affine, x0, method='broyden1', fatol=0.08_dp, tol_norm=tracking_norm)
    call check(r%success .and. calls == 6, &
        'selected norm called for residual, iterate and accepted step')

    r = root(affine, x0, method='broyden2', fatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted Broyden-2 stopping')
    r = root(affine, x0, method='anderson', fatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted Anderson stopping')
    r = root(affine, x0, method='linearmixing', fatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted linearmixing stopping')
    r = root(affine, x0, method='diagbroyden', fatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted diagbroyden stopping')
    r = root(affine, x0, method='excitingmixing', fatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted excitingmixing stopping')

    r = root(affine, [0.2_dp, 0.02_dp], method='broyden1', &
        tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 1 .and. r%nfev == 1, &
        'exact initial root still converges on the first check')

    ! Check each direct public wrapper, not just root(method=...).
    r = broyden1(affine, x0, f_tol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2, 'direct broyden1 norm')
    r = broyden2(affine, x0, f_tol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2, 'direct broyden2 norm')
    r = anderson(affine, x0, f_tol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2, 'direct anderson norm')
    r = linearmixing(affine, x0, f_tol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2, 'direct linearmixing norm')
    r = diagbroyden(affine, x0, f_tol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2, 'direct diagbroyden norm')
    r = excitingmixing(affine, x0, f_tol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2, 'direct excitingmixing norm')

    ! A step tolerance must also use the caller's norm, not just residuals.
    r = root(affine, x0, method='broyden1', fatol=1.0_dp, &
        xatol=0.08_dp, tol_norm=weighted_norm)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 7, &
        'weighted absolute step stopping')
    r = root(affine, x0, method='broyden1', fatol=1.0_dp, xatol=0.08_dp)
    call check(r%success .and. r%nit == 4, &
        'default step stopping parity')

    r = root(affine, x0, method='krylov', fatol=0.001_dp)
    call check(r%success .and. r%nit == 2 .and. r%nfev == 3, &
        'default Krylov tolerance unchanged')
    r = root(affine, x0, method='krylov', fatol=0.001_dp, tol_norm=amplified_norm)
    call check(r%success .and. r%nit == 3 .and. r%nfev == 6, &
        'Krylov custom convergence norm')
    r = newton_krylov(affine, x0, f_tol=0.001_dp, tol_norm=amplified_norm)
    call check(r%success .and. r%nit == 3 .and. r%nfev == 6, &
        'direct newton_krylov norm')

    r = root(affine, x0, method='df-sane', tol_norm=weighted_norm)
    call check(.not. r%success .and. r%status == 0 .and. r%nfev == 0, &
        'custom norm rejected for an unsupported algorithm')
    call check(index(r%message, 'tol_norm') > 0, 'unsupported method diagnostic')

    r = root(affine, x0, method='broyden1', tol_norm=negative_norm)
    call check(.not. r%success .and. r%status == 0 .and. r%nfev == 1, &
        'negative custom norm rejected without extra evaluations')
    r = root(affine, x0, method='krylov', tol_norm=negative_norm)
    call check(.not. r%success .and. r%status == 0 .and. r%nfev == 1, &
        'negative Krylov norm rejected')

    calls = 0
    r = root(affine, x0, method='anderson', tol_norm=bad_later_norm)
    call check(.not. r%success .and. r%status == 0 .and. r%nfev == 1, &
        'nonfinite norm during stopping check rejected')
    call check(calls >= 3, 'callback invoked for residual and iterate norms')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_tol_norm: failures=', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_root_tol_norm: PASS'

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! condition to assert
        character(len=*), intent(in) :: label !! diagnostic label
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', label
        end if
    end subroutine check

    subroutine affine(x, fvec)
        real(dp), intent(in) :: x(:) !! two-dimensional problem point
        real(dp), intent(out) :: fvec(:) !! unscaled affine residual vector
        fvec = x - [0.2_dp, 0.02_dp]
    end subroutine affine

    function weighted_norm(v) result(value)
        real(dp), intent(in) :: v(:) !! residual, iterate, or step
        real(dp) :: value !! anisotropic max norm
        value = maxval(abs(v) * [0.1_dp, 1.0_dp])
    end function weighted_norm

    function tracking_norm(v) result(value)
        real(dp), intent(in) :: v(:) !! residual, iterate, or step
        real(dp) :: value !! weighted max norm for call-count test
        calls = calls + 1
        value = weighted_norm(v)
    end function tracking_norm

    function amplified_norm(v) result(value)
        real(dp), intent(in) :: v(:) !! residual, iterate, or step
        real(dp) :: value !! scaled max norm for tiny Krylov residuals
        value = maxval(abs(v))
        if (value > huge(value) / 1.0e9_dp) then
            value = huge(value)
        else
            value = 1.0e9_dp * value
        end if
    end function amplified_norm

    function negative_norm(v) result(value)
        real(dp), intent(in) :: v(:) !! vector being measured
        real(dp) :: value !! invalid negative response
        value = -real(size(v), dp)
    end function negative_norm

    function bad_later_norm(v) result(value)
        real(dp), intent(in) :: v(:) !! residual, iterate, or step
        real(dp) :: value !! NaN starting at the third invocation
        calls = calls + 1
        value = maxval(abs(v))
        if (calls >= 3) value = ieee_value(0.0_dp, ieee_quiet_nan)
    end function bad_later_norm
end program test_optimize_root_tol_norm
