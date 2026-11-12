! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_legendre_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : ellipe, ellipeinc, ellipk, ellipkinc, ellipkm1
    implicit none

    real(dp), parameter :: half_pi = 0.5_dp * scifort_pi
    real(dp) :: m_values(4)
    real(dp) :: expected_e(4)
    real(dp) :: expected_k(4)
    real(dp) :: nan_value
    real(dp) :: pos_inf
    real(dp) :: values(4)

    pos_inf = ieee_value(1.0_dp, ieee_positive_inf)
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)

    m_values = [-2.0_dp, 0.0_dp, 0.5_dp, 0.9_dp]
    expected_k = [1.1714200841467701_dp, 1.5707963267948966_dp, &
        1.8540746773013719_dp, 2.5780921133481733_dp]
    expected_e = [2.1844381427462010_dp, 1.5707963267948966_dp, &
        1.3506438810476755_dp, 1.1047747327040733_dp]
    call check_array_rel('ellipk elemental', ellipk(m_values), expected_k, 2.0e-15_dp)
    call check_array_rel('ellipe elemental', ellipe(m_values), expected_e, 2.0e-15_dp)

    call check_rel('ellipkm1 tiny', ellipkm1(1.0e-300_dp), &
        346.774058310226743221533182445571_dp, 2.0e-15_dp)
    call check_rel('ellipkm1 large', ellipkm1(1.0e100_dp), &
        1.16515549010822174819734036977135e-48_dp, 2.0e-15_dp)
    call check_rel('ellipk extreme negative', ellipk(-1.0e200_dp), &
        2.31644803660524459020633609711353e-98_dp, 2.0e-15_dp)
    call check_rel('ellipe extreme negative', ellipe(-1.0e200_dp), 1.0e100_dp, 2.0e-15_dp)

    call check_rel('ellipkinc ordinary', ellipkinc(1.0_dp, 0.5_dp), &
        1.0832167728451687_dp, 2.0e-15_dp)
    call check_rel('ellipeinc ordinary', ellipeinc(1.0_dp, 0.5_dp), &
        0.9273298836244400_dp, 2.0e-15_dp)
    call check_rel('ellipkinc negative parameter', ellipkinc(1.0_dp, -0.5_dp), &
        0.9414060121877065_dp, 2.0e-15_dp)
    call check_rel('ellipeinc negative parameter', ellipeinc(1.0_dp, -0.5_dp), &
        1.0647286545290244_dp, 2.0e-15_dp)
    call check_rel('ellipkinc periodic', ellipkinc(10.0_dp, 0.5_dp), &
        11.715622315665893_dp, 3.0e-15_dp)
    call check_rel('ellipeinc periodic', ellipeinc(10.0_dp, 0.5_dp), &
        8.663886106525743_dp, 3.0e-15_dp)
    call check_rel('ellipkinc negative periodic', ellipkinc(-10.0_dp, 0.5_dp), &
        -11.715622315665893_dp, 3.0e-15_dp)
    call check_rel('ellipeinc negative periodic', ellipeinc(-10.0_dp, 0.5_dp), &
        -8.663886106525743_dp, 3.0e-15_dp)

    call check_rel('ellipkinc complete identity', ellipkinc(half_pi, 0.9_dp), &
        ellipk(0.9_dp), 2.0e-15_dp)
    call check_rel('ellipeinc complete identity', ellipeinc(half_pi, 0.9_dp), &
        ellipe(0.9_dp), 2.0e-15_dp)
    call check_rel('ellipkinc period identity', ellipkinc(1.0_dp + scifort_pi, 0.5_dp), &
        ellipkinc(1.0_dp, 0.5_dp) + 2.0_dp * ellipk(0.5_dp), 3.0e-15_dp)
    call check_rel('ellipeinc period identity', ellipeinc(1.0_dp + scifort_pi, 0.5_dp), &
        ellipeinc(1.0_dp, 0.5_dp) + 2.0_dp * ellipe(0.5_dp), 3.0e-15_dp)

    call check_close('ellipk zero parameter', ellipk(0.0_dp), half_pi, 0.0_dp)
    call check_close('ellipe zero parameter', ellipe(0.0_dp), half_pi, 0.0_dp)
    call check_close('ellipkm1 unit parameter', ellipkm1(1.0_dp), half_pi, 0.0_dp)
    call check_close('ellipkinc zero parameter', ellipkinc(123.5_dp, 0.0_dp), 123.5_dp, 0.0_dp)
    call check_close('ellipeinc zero parameter', ellipeinc(-123.5_dp, 0.0_dp), -123.5_dp, 0.0_dp)

    call check_true('ellipk m=1 infinity', ellipk(1.0_dp) > huge(1.0_dp))
    call check_close('ellipe m=1', ellipe(1.0_dp), 1.0_dp, 0.0_dp)
    call check_true('ellipkinc m=1 endpoint', ellipkinc(half_pi, 1.0_dp) > huge(1.0_dp))
    call check_true('ellipkinc m=1 negative crossing', &
        ellipkinc(-2.0_dp, 1.0_dp) < -huge(1.0_dp))
    call check_close('ellipeinc m=1 endpoint', ellipeinc(half_pi, 1.0_dp), 1.0_dp, 0.0_dp)
    call check_rel('ellipeinc m=1 periodic', ellipeinc(2.0_dp, 1.0_dp), &
        1.0907025731743185_dp, 2.0e-15_dp)

    call check_true('ellipk invalid m', ieee_is_nan(ellipk(1.1_dp)))
    call check_true('ellipe invalid m', ieee_is_nan(ellipe(1.1_dp)))
    call check_true('ellipkm1 invalid p', ieee_is_nan(ellipkm1(-0.1_dp)))
    call check_true('ellipkinc invalid m', ieee_is_nan(ellipkinc(0.0_dp, 1.1_dp)))
    call check_true('ellipeinc invalid m', ieee_is_nan(ellipeinc(0.0_dp, 1.1_dp)))
    call check_true('ellipk nan', ieee_is_nan(ellipk(nan_value)))
    call check_true('ellipe nan', ieee_is_nan(ellipe(nan_value)))
    call check_true('ellipkinc nan phi', ieee_is_nan(ellipkinc(nan_value, 0.5_dp)))
    call check_true('ellipeinc nan m', ieee_is_nan(ellipeinc(1.0_dp, nan_value)))

    call check_close('ellipk negative infinity', ellipk(-pos_inf), 0.0_dp, 0.0_dp)
    call check_true('ellipe negative infinity', ellipe(-pos_inf) > huge(1.0_dp))
    call check_close('ellipkm1 positive infinity', ellipkm1(pos_inf), 0.0_dp, 0.0_dp)
    call check_true('ellipkm1 zero singular', ellipkm1(0.0_dp) > huge(1.0_dp))
    call check_true('ellipk positive infinity invalid', ieee_is_nan(ellipk(pos_inf)))
    call check_true('ellipe positive infinity invalid', ieee_is_nan(ellipe(pos_inf)))
    call check_true('ellipkinc positive infinity', ellipkinc(pos_inf, 0.5_dp) > huge(1.0_dp))
    call check_true('ellipkinc negative infinity', ellipkinc(-pos_inf, 0.5_dp) < -huge(1.0_dp))
    call check_true('ellipeinc positive infinity', ellipeinc(pos_inf, 0.5_dp) > huge(1.0_dp))
    call check_true('ellipeinc negative infinity', ellipeinc(-pos_inf, 0.5_dp) < -huge(1.0_dp))

    call check_close('ellipkinc m=-infinity finite', ellipkinc(2.0_dp, -pos_inf), 0.0_dp, 0.0_dp)
    call check_true('ellipeinc m=-infinity positive', ellipeinc(2.0_dp, -pos_inf) > huge(1.0_dp))
    call check_true('ellipeinc m=-infinity negative', ellipeinc(-2.0_dp, -pos_inf) < -huge(1.0_dp))
    call check_true('ellipkinc infinity with m=-infinity', &
        ieee_is_nan(ellipkinc(pos_inf, -pos_inf)))

    call check_true('ellipkinc huge overflow guarded', &
        ellipkinc(huge(1.0_dp), 0.9_dp) > huge(1.0_dp))
    call check_true('ellipeinc huge finite', ieee_is_finite(ellipeinc(huge(1.0_dp), 0.9_dp)))
    call check_true('ellipkinc huge negative parameter finite', &
        ieee_is_finite(ellipkinc(huge(1.0_dp), -2.0_dp)))
    call check_true('ellipeinc huge negative parameter overflow guarded', &
        ellipeinc(huge(1.0_dp), -2.0_dp) > huge(1.0_dp))

    values = ellipk([0.0_dp, 0.5_dp, 0.9_dp, 1.0_dp])
    call check_close('ellipk array first', values(1), half_pi, 0.0_dp)
    call check_true('ellipk array last', values(4) > huge(1.0_dp))

    print *, 'Legendre elliptic special-function tests passed'

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

    subroutine check_array_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! nonzero reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise relative error

        if (size(actual) /= size(expected) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance * abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_legendre_elliptic
