! SPDX-License-Identifier: MIT
program test_linalg_funm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_funm, only : funm, funm_real
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(3,3), invalid(2,3), real_result(3,3), defect
    real(dp) :: identity(3,3), jordan(3,3), expected(3,3), small(2,2)
    real(dp) :: mixed(7,7)
    complex(dp) :: mixed_f(7,7)
    complex(dp) :: c(3,3), x(3,3), y(2,2)
    integer :: i, status

    identity = 0.0_dp
    do i = 1, 3
        identity(i, i) = 1.0_dp
    end do
    a = reshape([0.4_dp, 0.2_dp, -0.1_dp, &
        0.8_dp, 0.9_dp, -0.3_dp, 0.1_dp, 0.3_dp, 1.5_dp], [3,3])
    call funm(a, square, x, status, commutator_residual=defect)
    if (status /= linalg_status_success) error stop 'real polynomial status'
    if (maxval(abs(x-matmul(a,a))) > 1.0e-10_dp .or. defect > 1.0e-10_dp) &
        error stop 'real polynomial accuracy'
    call funm_real(a, square, real_result, status)
    if (status /= linalg_status_success) error stop 'real output status'
    if (maxval(abs(real_result-matmul(a,a))) > 1.0e-10_dp) error stop 'real output'

    c = cmplx(a, kind=dp)
    c(1,2) = cmplx(0.4_dp,0.1_dp,dp)
    c(2,1) = cmplx(-0.3_dp,0.2_dp,dp)
    call funm(c, square, x, status)
    if (status /= linalg_status_success) error stop 'complex polynomial status'
    if (maxval(abs(x-matmul(c,c))) > 1.0e-10_dp) error stop 'complex polynomial'

    ! Jordan block: the derivative callback supplies the needed Hermite data.
    jordan = 0.0_dp
    do i = 1, 3
        jordan(i,i) = 0.3_dp
    end do
    jordan(1,2) = 2.0_dp
    jordan(2,3) = -0.6_dp
    expected = exp(0.3_dp)*(identity+jordan-0.3_dp*identity+ &
        matmul(jordan-0.3_dp*identity, jordan-0.3_dp*identity)/2.0_dp)
    call funm(jordan, exponential, x, status, derivative=exponential_derivative)
    if (status /= linalg_status_success) error stop 'Jordan block status'
    if (maxval(abs(x-expected)) > 1.0e-10_dp) error stop 'Jordan block values'
    call funm(jordan, exponential, x, status)
    if (status /= linalg_status_no_convergence) error stop 'Jordan block needs derivative'

    ! Repeated eigenvalues in positions 1 and 3 are separated by a
    ! distinct eigenvalue.  The off-diagonal chain creates a Jordan
    ! coupling that the original adjacent-only Schur-Parlett could not
    ! resolve even when analytic derivatives were available.
    a = 0.0_dp
    a(1,1) = 0.5_dp
    a(2,2) = 1.7_dp
    a(3,3) = 0.5_dp
    a(1,2) = 0.2_dp
    a(2,3) = -0.7_dp
    a(1,3) = 0.05_dp
    call funm(a, square, x, status, derivative=square_derivative)
    if (status /= linalg_status_success) error stop 'noncontiguous real status'
    if (maxval(abs(x-matmul(a,a))) > 1.0e-11_dp) &
        error stop 'noncontiguous real values'
    call funm_real(a, square, real_result, status, derivative=square_derivative)
    if (status /= linalg_status_success) error stop 'noncontiguous real-only status'
    if (maxval(abs(real_result-matmul(a,a))) > 1.0e-11_dp) &
        error stop 'noncontiguous real-only values'
    call funm(a, square, x, status)
    if (status /= linalg_status_no_convergence) &
        error stop 'noncontiguous real needs derivative'

    ! The same noncontiguous repeated eigenvalue can be *semisimple*:
    ! choosing the (1,3) entry to cancel the induced Jordan coupling
    ! makes the spectrum diagonalizable.  No derivative is then required.
    a(1,3) = a(1,2)*a(2,3)/(a(2,2)-a(1,1))
    call funm(a, exponential, x, status)
    if (status /= linalg_status_success) error stop 'diagonalizable repeated status'
    ! exp(A) is a similarity transform of exp(diag(A)); verify via
    ! a polynomial identity for this specific three-eigenvalue matrix:
    ! exp(A) = exp(lambda)*I + divided_difference*(A-lambda*I),
    ! since its minimal polynomial is (z-lambda)*(z-mu).
    expected = exp(0.5_dp)*identity + &
        ((exp(1.7_dp)-exp(0.5_dp))/(1.7_dp-0.5_dp))*(a-0.5_dp*identity)
    if (maxval(abs(x-expected)) > 1.0e-11_dp) &
        error stop 'diagonalizable repeated values'

    ! Different *sizes* of repeated eigenvalue blocks require a block
    ! Sylvester solve.  The polynomial reference is exact at all orders.
    mixed = 0.0_dp
    do i = 1, 7
        mixed(i,i) = 0.4_dp
        if (mod(i,2) == 0) mixed(i,i) = 1.3_dp
    end do
    mixed(1,2) = 0.2_dp
    mixed(1,3) = 0.1_dp
    mixed(1,7) = 0.17_dp
    mixed(2,4) = -0.3_dp
    mixed(3,5) = 0.3_dp
    mixed(4,5) = -0.12_dp
    mixed(5,6) = 0.14_dp
    mixed(6,7) = -0.09_dp
    call funm(mixed, square, mixed_f, status, derivative=square_derivative)
    if (status /= linalg_status_success) error stop 'mixed Schur block status'
    if (maxval(abs(mixed_f-matmul(mixed,mixed))) > 2.0e-11_dp) &
        error stop 'mixed Schur block values'

    c = cmplx(a, kind=dp)
    c(1,1) = cmplx(0.5_dp, 0.2_dp, dp)
    c(3,3) = c(1,1)
    c(2,2) = cmplx(1.7_dp, -0.3_dp, dp)
    c(1,2) = cmplx(0.2_dp, 0.1_dp, dp)
    c(2,3) = cmplx(-0.7_dp, -0.15_dp, dp)
    call funm(c, square, x, status, derivative=square_derivative)
    if (status /= linalg_status_success) error stop 'noncontiguous complex status'
    if (maxval(abs(x-matmul(c,c))) > 1.0e-11_dp) &
        error stop 'noncontiguous complex values'

    a = 0.0_dp
    a(1,1) = -1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = -1.0_dp
    call funm(a, square, x, status)
    if (status /= linalg_status_success) error stop 'separate repeated diagonals'
    if (maxval(abs(x-matmul(a,a))) > 1.0e-12_dp) error stop 'repeated diagonal values'

    a = 0.0_dp
    call funm(a, exponential, x, status)
    if (status /= linalg_status_success) error stop 'zero status'
    if (maxval(abs(x-cmplx(identity,kind=dp))) > 1.0e-13_dp) error stop 'zero matrix'
    call funm_real(a, exponential, real_result, status)
    if (status /= linalg_status_success) error stop 'zero real output'
    if (maxval(abs(real_result-identity)) > 1.0e-13_dp) error stop 'zero real values'

    a = 0.0_dp
    a(1,1) = -2.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp
    call funm_real(a, logarithm, real_result, status)
    if (status /= linalg_status_no_convergence) error stop 'non-real output'

    invalid = 0.0_dp
    call funm(invalid, square, x, status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid shape'
    call funm(a, square, x, status, separation_tol=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'invalid gap threshold'
    a(2,1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call funm(a, square, x, status)
    if (status /= linalg_status_invalid_shape) error stop 'nonfinite input'

    small = 0.0_dp
    small(1,1) = 0.5_dp
    small(2,2) = 0.5_dp
    small(1,2) = 1.3_dp
    call funm(small, exponential, y, status, derivative=exponential_derivative)
    if (status /= linalg_status_success) error stop '2x2 Jordan'
    if (abs(y(1,2)-1.3_dp*exp(0.5_dp)) > 1.0e-12_dp) error stop '2x2 derivative'

    ! Unlike a true Jordan block, a clustered but distinct spectrum has a
    ! nonterminating Taylor series.  The historical n-1 truncation must
    ! reject it when its nilpotence screen fails.  The opt-in extended
    ! derivative limit evaluates the complete converged series.
    small = 0.0_dp
    small(1,1) = 1.3_dp
    small(2,2) = 1.303_dp
    small(1,2) = 0.8_dp
    call funm(small, exponential, y, status, separation_tol=0.02_dp, &
        derivative=exponential_derivative)
    if (status /= linalg_status_no_convergence) &
        error stop 'default short Taylor should reject distinct cluster'
    call funm(small, exponential, y, status, separation_tol=0.02_dp, &
        derivative=exponential_derivative, taylor_max_order=96)
    if (status /= linalg_status_success) error stop 'adaptive Taylor exp status'
    if (abs(y(1,1)-exp(small(1,1))) > 1.0e-12_dp) error stop 'adaptive Taylor diagonal'
    if (abs(y(1,2)-0.8_dp*(exp(1.303_dp)-exp(1.3_dp))/0.003_dp) > 5.0e-12_dp) &
        error stop 'adaptive Taylor exp divided difference'
    call funm_real(small, exponential, real_result(1:2,1:2), status, &
        separation_tol=0.02_dp, derivative=exponential_derivative, taylor_max_order=96)
    if (status /= linalg_status_success) error stop 'adaptive Taylor real status'
    if (maxval(abs(real_result(1:2,1:2)-real(y,dp))) > 1.0e-12_dp) &
        error stop 'adaptive Taylor real result'

    ! A low maximum is not silently interpreted as successful truncation.
    call funm(small, exponential, y, status, separation_tol=0.02_dp, &
        derivative=exponential_derivative, taylor_max_order=2)
    if (status /= linalg_status_no_convergence) error stop 'short adaptive series'
    call funm(small, exponential, y, status, taylor_max_order=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid Taylor order low'
    call funm(small, exponential, y, status, taylor_max_order=129)
    if (status /= linalg_status_invalid_shape) error stop 'invalid Taylor order high'

    ! An extended limit must still compute exact nilpotent/Jordan blocks.
    small(2,2) = small(1,1)
    call funm(small, exponential, y, status, derivative=exponential_derivative, &
        taylor_max_order=96)
    if (status /= linalg_status_success) error stop 'adaptive nilpotent status'
    if (abs(y(1,2)-0.8_dp*exp(1.3_dp)) > 1.0e-12_dp) &
        error stop 'adaptive nilpotent offdiagonal'

    print '(a)', 'PASS test_linalg_funm'

contains
    function square(z) result(value)
        complex(dp), intent(in) :: z !! scalar to square
        complex(dp) :: value
        value = z*z
    end function square

    function square_derivative(z, order) result(value)
        complex(dp), intent(in) :: z !! scalar point for the k-th derivative
        integer, intent(in) :: order !! derivative order, k >= 1
        complex(dp) :: value
        value = cmplx(0.0_dp, 0.0_dp, dp)
        if (order == 1) value = 2.0_dp*z
        if (order == 2) value = cmplx(2.0_dp, 0.0_dp, dp)
    end function square_derivative

    function exponential(z) result(value)
        complex(dp), intent(in) :: z !! scalar for exponential
        complex(dp) :: value
        value = exp(z)
    end function exponential

    function logarithm(z) result(value)
        complex(dp), intent(in) :: z !! scalar logarithm argument
        complex(dp) :: value
        value = log(z)
    end function logarithm

    function exponential_derivative(z, order) result(value)
        complex(dp), intent(in) :: z !! derivative evaluation point
        integer, intent(in) :: order !! derivative order, positive
        complex(dp) :: value
        value = exp(z)
    end function exponential_derivative
end program test_linalg_funm
