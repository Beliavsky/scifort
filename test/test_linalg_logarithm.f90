! SPDX-License-Identifier: MIT
program test_linalg_logarithm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence, linalg_status_singular
    use scifort_linalg_logarithm, only : logm, logm_real, fractional_matrix_power, &
        fractional_matrix_power_real
    use scifort_linalg_expm_frechet, only : expm_frechet
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(2,2), err, real_log(2,2), real_pow(2,2)
    complex(dp) :: x(2,2), z(2,2), c(2,2), work(2,2)
    integer :: status, roots, terms
    real(dp), parameter :: pi = acos(-1.0_dp)

    a = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,2) = 4.0_dp
    a(1,2) = 3.0_dp
    call logm(a, x, status, err, roots, terms)
    call assert_ok(status == linalg_status_success, 'real triangular log status')
    call assert_ok(err < 1.0e-11_dp .and. terms >= 1, 'log residual and work')
    call assert_ok(abs(x(1,1)-log(2.0_dp)) < 1.0e-11_dp, 'log diagonal 1')
    call assert_ok(abs(x(2,2)-log(4.0_dp)) < 1.0e-11_dp, 'log diagonal 2')
    call assert_ok(abs(x(1,2)-1.5_dp*log(2.0_dp)) < 1.0e-11_dp, &
        'log divided difference')
    call logm_real(a, real_log, status)
    call assert_ok(status == linalg_status_success, 'real logarithm convenience')
    call assert_ok(maxval(abs(real_log-real(x,dp))) < 1.0e-11_dp, 'real log agreement')

    call fractional_matrix_power(a, 0.5_dp, z, status)
    call assert_ok(status == linalg_status_success, 'sqrt power')
    call assert_ok(maxval(abs(matmul(z,z)-cmplx(a, kind=dp))) < 1.0e-11_dp, &
        'sqrt power reconstruct')
    call fractional_matrix_power(a, 1.5_dp, z, status)
    call assert_ok(status == linalg_status_success, 'fractional 3/2 power')
    call assert_ok(maxval(abs(z-matmul(cmplx(a,kind=dp), &
        matrix_sqrt_reference()))) < 1.0e-10_dp, 'fractional power reference')
    call fractional_matrix_power_real(a, 1.5_dp, real_pow, status)
    call assert_ok(status == linalg_status_success, 'real power convenience')
    call assert_ok(maxval(abs(real_pow-real(z,dp))) < 1.0e-10_dp, 'real power agreement')

    ! Invertible nonnormal Jordan block, log(I+N)=N with N**2=0.
    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 1.0_dp
    a(1,2) = 5.0_dp
    call logm(a, x, status, err)
    call assert_ok(status == linalg_status_success .and. err < 1.0e-10_dp, 'Jordan log')
    call assert_ok(abs(x(1,2)-5.0_dp) < 2.0e-9_dp, 'Jordan log offdiagonal')
    call fractional_matrix_power(a, -2.0_dp, z, status)
    call assert_ok(status == linalg_status_success, 'negative integer power')
    call assert_ok(abs(z(1,2)+10.0_dp) < 1.0e-11_dp, 'negative integer result')

    a = 0.0_dp
    a(1,1) = -1.0_dp
    a(2,2) = 4.0_dp
    call logm(a, x, status, err)
    call assert_ok(status == linalg_status_success, 'real negative logarithm')
    call assert_ok(abs(x(1,1)-cmplx(0.0_dp,pi,dp)) < 1.0e-10_dp, &
        'negative eigenvalue principal phase')
    call assert_ok(abs(x(2,2)-log(4.0_dp)) < 1.0e-10_dp, 'positive eigenvalue log')
    call logm_real(a, real_log, status)
    call assert_ok(status == linalg_status_no_convergence, 'non-real log refused')
    call fractional_matrix_power(a, 0.5_dp, z, status)
    call assert_ok(status == linalg_status_success, 'negative eigenvalue square root')
    call assert_ok(abs(z(1,1)-cmplx(0.0_dp,1.0_dp,dp)) < 1.0e-10_dp, 'principal sqrt phase')
    call fractional_matrix_power_real(a, 0.5_dp, real_pow, status)
    call assert_ok(status == linalg_status_no_convergence, 'non-real power refused')

    c = cmplx(0.0_dp, 0.0_dp, dp)
    c(1,1) = cmplx(1.2_dp,0.3_dp,dp)
    c(2,2) = cmplx(0.7_dp,-0.4_dp,dp)
    c(1,2) = cmplx(0.2_dp,0.15_dp,dp)
    c(2,1) = cmplx(-0.1_dp,0.2_dp,dp)
    call logm(c, x, status, err)
    call assert_ok(status == linalg_status_success .and. err < 1.0e-11_dp, 'complex log')
    call expm_frechet(x, cmplx(0.0_dp, 0.0_dp, dp)*c, z, work, status)
    call assert_ok(status == linalg_status_success, 'complex exp')
    call assert_ok(maxval(abs(z-c)) < 1.0e-10_dp, 'complex exp-log')
    call fractional_matrix_power(c, 0.3_dp, z, status)
    call assert_ok(status == linalg_status_success, 'complex fractional power')
    call fractional_matrix_power(c, -0.3_dp, work, status)
    call assert_ok(status == linalg_status_success, 'complex negative fractional power')
    call assert_ok(maxval(abs(matmul(z,work)-identity2c())) < 1.0e-9_dp, &
        'inverse fractional powers')

    ! Scalar scaling should permit moderate huge and tiny diagonal values.
    a = 0.0_dp
    a(1,1) = 1.0e180_dp
    a(2,2) = 1.0e180_dp
    call logm(a, x, status, err)
    call assert_ok(status == linalg_status_success .and. err < 1.0e-12_dp, 'large scale log')
    call assert_ok(abs(x(1,1)-log(1.0e180_dp)) < 1.0e-10_dp, 'large log value')
    a = 0.0_dp
    a(1,1) = 1.0e-180_dp
    a(2,2) = 1.0e-180_dp
    call logm(a, x, status, err)
    call assert_ok(status == linalg_status_success .and. err < 1.0e-12_dp, 'tiny scale log')
    call assert_ok(abs(x(1,1)-log(1.0e-180_dp)) < 1.0e-10_dp, 'tiny log value')

    ! The half-power path must normalize subnormal-scale matrices before
    ! invoking the existing Schur square-root kernel.
    a = 0.0_dp
    a(1,1) = 1.0e-180_dp
    a(2,2) = 1.0e-180_dp
    call fractional_matrix_power(a, 0.5_dp, z, status)
    call assert_ok(status == linalg_status_success, 'tiny-scale half power')
    call assert_ok(abs(z(1,1)-1.0e-90_dp) < 1.0e-103_dp, &
        'tiny-scale sqrt value')
    a(1,1) = 1.0e180_dp
    a(2,2) = 1.0e180_dp
    call fractional_matrix_power(a, 3.0_dp, z, status)
    call assert_ok(status == linalg_status_no_convergence, 'nonrepresentable power')

    ! Singular input is rejected by logm, but nonnegative integer powers
    ! and the special half-power can still be defined.
    a = 0.0_dp
    a(1,1) = 4.0_dp
    call logm(a, x, status)
    call assert_ok(status == linalg_status_singular, 'singular log')
    call fractional_matrix_power(a, 2.0_dp, z, status)
    call assert_ok(status == linalg_status_success .and. abs(z(1,1)-16.0_dp) < 1.0e-11_dp, &
        'singular positive integer power')
    call fractional_matrix_power(a, 0.5_dp, z, status)
    call assert_ok(status == linalg_status_success .and. abs(z(1,1)-2.0_dp) < 1.0e-11_dp, &
        'singular square root')
    call fractional_matrix_power(a, -1.0_dp, z, status)
    call assert_ok(status == linalg_status_singular, 'singular inverse power')
    call fractional_matrix_power(a, 0.0_dp, z, status)
    call assert_ok(status == linalg_status_success .and. all(z == identity2c()), 'zero power')

    call logm(a(1:2,1:1), x, status)
    call assert_ok(status == linalg_status_invalid_shape, 'invalid log shape')
    call fractional_matrix_power(a, ieee_value(0.0_dp,ieee_quiet_nan), z, status)
    call assert_ok(status == linalg_status_invalid_shape, 'invalid exponent')
    call logm(a, x, status, tol=-1.0_dp)
    call assert_ok(status == linalg_status_invalid_shape, 'invalid log tolerance')
    call logm(c, x, status, max_sqrt=0)
    call assert_ok(status == linalg_status_no_convergence, 'square-root budget')
    call logm(c, x, status, max_terms=1)
    call assert_ok(status == linalg_status_no_convergence, 'series budget')
    a(1,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    call logm(a, x, status)
    call assert_ok(status == linalg_status_invalid_shape, 'nonfinite input')

    print '(a)', 'PASS test_linalg_logarithm'
contains
    subroutine assert_ok(ok, label)
        logical, intent(in) :: ok !! test condition
        character(*), intent(in) :: label !! failed test name
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine assert_ok

    pure function identity2c() result(id)
        complex(dp) :: id(2,2)
        id = cmplx(0.0_dp,0.0_dp,dp)
        id(1,1) = cmplx(1.0_dp,0.0_dp,dp)
        id(2,2) = cmplx(1.0_dp,0.0_dp,dp)
    end function identity2c

    pure function matrix_sqrt_reference() result(id)
        complex(dp) :: id(2,2)
        id = cmplx(0.0_dp,0.0_dp,dp)
        id(1,1) = cmplx(sqrt(2.0_dp),0.0_dp,dp)
        id(2,2) = cmplx(2.0_dp,0.0_dp,dp)
        id(1,2) = cmplx(3.0_dp/(sqrt(2.0_dp)+2.0_dp),0.0_dp,dp)
    end function matrix_sqrt_reference
end program test_linalg_logarithm
