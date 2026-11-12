! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_gamma_family
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : digamma, gamma, gammasgn, loggamma, multigammaln, &
        rgamma
    implicit none

    real(dp) :: nan_value
    real(dp) :: neg_inf
    real(dp) :: neg_zero
    real(dp) :: pos_inf
    real(dp) :: values(5)

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)
    neg_inf = ieee_value(1.0_dp, ieee_negative_inf)
    neg_zero = -0.0_dp

    call check_rel('gamma negative half', gamma(-0.5_dp), &
        -3.5449077018110318_dp, 3.0e-15_dp)
    call check_rel('gamma negative ordinary', gamma(-4.75_dp), &
        -0.05639286923398244_dp, 5.0e-15_dp)
    call check_rel('gamma positive ordinary', gamma(20.0_dp), &
        1.21645100408832e17_dp, 5.0e-15_dp)
    call check_rel('gamma near underflow', gamma(-170.5_dp), &
        -3.3127395215386074e-308_dp, 2.0e-13_dp)
    call check_rel('gamma subnormal 1', gamma(-171.5_dp), &
        1.9316265431711996e-310_dp, 3.0e-13_dp)
    call check_rel('gamma subnormal 2', gamma(-172.5_dp), &
        -1.119783503287652e-312_dp, 3.0e-12_dp)
    call check_true('gamma underflow sign', sign(1.0_dp, gamma(-178.5_dp)) < 0.0_dp)
    call check_close('gamma underflow zero', gamma(-178.5_dp), -0.0_dp, 0.0_dp)
    call check_true('gamma overflow', gamma(172.0_dp) > huge(1.0_dp))

    call check_rel('rgamma negative half', rgamma(-0.5_dp), &
        -0.28209479177387814_dp, 3.0e-15_dp)
    call check_rel('rgamma negative ordinary', rgamma(-4.75_dp), &
        -17.73273843986995_dp, 5.0e-15_dp)
    call check_rel('rgamma positive ordinary', rgamma(20.0_dp), &
        8.22063524662433e-18_dp, 5.0e-15_dp)
    call check_rel('rgamma near overflow', rgamma(-170.5_dp), &
        -3.018649650835054e307_dp, 2.0e-13_dp)
    call check_true('rgamma overflow', rgamma(-171.5_dp) > huge(1.0_dp))
    call check_rel('rgamma positive subnormal', rgamma(172.0_dp), &
        8.057900396443103e-310_dp, 3.0e-13_dp)
    call check_close('rgamma positive underflow', rgamma(179.0_dp), 0.0_dp, 0.0_dp)

    call check_close('gammasgn positive', gammasgn(10.0_dp), 1.0_dp, 0.0_dp)
    call check_close('gammasgn negative interval 1', gammasgn(-0.5_dp), -1.0_dp, 0.0_dp)
    call check_close('gammasgn negative interval 2', gammasgn(-1.5_dp), 1.0_dp, 0.0_dp)
    call check_close('gammasgn negative interval 3', gammasgn(-2.5_dp), -1.0_dp, 0.0_dp)

    call check_rel('loggamma positive', loggamma(20.0_dp), &
        39.339884187199495_dp, 3.0e-15_dp)
    call check_close('loggamma one', loggamma(1.0_dp), 0.0_dp, 0.0_dp)
    call check_true('loggamma negative real branch', ieee_is_nan(loggamma(-0.5_dp)))

    call check_rel('multigammaln d2', multigammaln(2.5_dp, 2), &
        0.8570478133976193_dp, 4.0e-15_dp)
    call check_rel('multigammaln d3', multigammaln(3.5_dp, 3), &
        3.895898482154039_dp, 4.0e-15_dp)
    call check_rel('multigammaln d5', multigammaln(10.0_dp, 5), &
        58.8938418512374_dp, 4.0e-15_dp)
    call check_rel('multigammaln d8', multigammaln(100.0_dp, 8), &
        2824.87514827951_dp, 4.0e-15_dp)

    call check_true('gamma positive zero infinity', gamma(0.0_dp) > huge(1.0_dp))
    call check_true('gamma negative zero infinity', gamma(neg_zero) < -huge(1.0_dp))
    call check_true('gamma negative integer nan', ieee_is_nan(gamma(-2.0_dp)))
    call check_true('gamma positive infinity', gamma(pos_inf) > huge(1.0_dp))
    call check_true('gamma negative infinity nan', ieee_is_nan(gamma(neg_inf)))
    call check_true('gamma nan', ieee_is_nan(gamma(nan_value)))

    call check_true('rgamma positive zero sign', sign(1.0_dp, rgamma(0.0_dp)) > 0.0_dp)
    call check_true('rgamma negative zero sign', sign(1.0_dp, rgamma(neg_zero)) < 0.0_dp)
    call check_close('rgamma negative integer zero', rgamma(-2.0_dp), 0.0_dp, 0.0_dp)
    call check_close('rgamma positive infinity zero', rgamma(pos_inf), 0.0_dp, 0.0_dp)
    call check_close('rgamma negative infinity zero', rgamma(neg_inf), 0.0_dp, 0.0_dp)
    call check_true('rgamma nan', ieee_is_nan(rgamma(nan_value)))

    call check_close('gammasgn positive zero', gammasgn(0.0_dp), 1.0_dp, 0.0_dp)
    call check_close('gammasgn negative zero', gammasgn(neg_zero), -1.0_dp, 0.0_dp)
    call check_true('gammasgn pole nan', ieee_is_nan(gammasgn(-3.0_dp)))
    call check_close('gammasgn positive infinity', gammasgn(pos_inf), 1.0_dp, 0.0_dp)
    call check_true('gammasgn negative infinity nan', ieee_is_nan(gammasgn(neg_inf)))

    call check_true('loggamma zero infinity', loggamma(0.0_dp) > huge(1.0_dp))
    call check_true('loggamma negative zero infinity', loggamma(neg_zero) > huge(1.0_dp))
    call check_true('loggamma positive infinity', loggamma(pos_inf) > huge(1.0_dp))
    call check_true('loggamma negative infinity nan', ieee_is_nan(loggamma(neg_inf)))

    call check_true('multigammaln invalid dimension', ieee_is_nan(multigammaln(2.0_dp, 0)))
    call check_true('multigammaln invalid domain', ieee_is_nan(multigammaln(1.0_dp, 4)))
    call check_true('multigammaln infinity', multigammaln(pos_inf, 4) > huge(1.0_dp))

    call check_true('digamma positive zero', digamma(0.0_dp) < -huge(1.0_dp))
    call check_true('digamma negative zero', digamma(neg_zero) > huge(1.0_dp))

    values = gamma([-2.5_dp, -0.5_dp, 0.5_dp, 1.0_dp, 5.0_dp])
    call check_rel('gamma elemental first', values(1), -0.9453087204829419_dp, 3.0e-15_dp)
    call check_rel('gamma elemental last', values(5), 24.0_dp, 3.0e-15_dp)

    values = rgamma([-2.5_dp, -0.5_dp, 0.5_dp, 1.0_dp, 5.0_dp])
    call check_rel('rgamma elemental first', values(1), -1.057855469152043_dp, 3.0e-15_dp)
    call check_rel('rgamma elemental last', values(5), 1.0_dp / 24.0_dp, 3.0e-15_dp)

    print *, 'Gamma-family special-function tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. &
                abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_gamma_family
