program test_special_wright_bessel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_special, only : dp, log_wright_bessel, wright_bessel
    implicit none

    real(dp) :: a_values(3)
    real(dp) :: b_values(3)
    real(dp) :: lhs
    real(dp) :: rhs
    real(dp) :: values(3)
    real(dp) :: x_values(3)

    call check_close(wright_bessel(0.0_dp, 0.0_dp, 1.0_dp), 0.0_dp, 0.0_dp, 'Phi(0,0;1)')
    call require(.not. ieee_is_finite(log_wright_bessel(0.0_dp, 0.0_dp, 1.0_dp)) .and. &
        log_wright_bessel(0.0_dp, 0.0_dp, 1.0_dp) < 0.0_dp, 'log Phi(0,0;1)')
    call check_close(wright_bessel(0.0_dp, 1.0_dp, 1.0_dp), &
        2.7182818284590452353602874713527_dp, 2.0e-15_dp, 'Phi(0,1;1)')
    call check_close(log_wright_bessel(0.0_dp, 2.0_dp, 3.0_dp), 3.0_dp, 0.0_dp, &
        'log Phi(0,2;3)')

    call check_close(wright_bessel(1.0_dp, 0.0_dp, 1.0_dp), &
        1.5906368546373290633822544249997_dp, 3.0e-15_dp, 'Phi(1,0;1)')
    call check_close(wright_bessel(1.0_dp, 1.0_dp, 1.0_dp), &
        2.2795853023360672674372044408115_dp, 3.0e-15_dp, 'Phi(1,1;1)')
    call check_close(wright_bessel(1.5_dp, 1.1_dp, 2.5_dp), &
        3.3027831556598564651664788589047_dp, 4.0e-15_dp, 'Phi(1.5,1.1;2.5)')
    call check_close(wright_bessel(0.5_dp, 0.25_dp, 10.0_dp), &
        2.7159507085794932365379053024222e3_dp, 6.0e-15_dp, 'Phi(0.5,0.25;10)')
    call check_close(log_wright_bessel(0.1_dp, 2.0_dp, 100.0_dp), &
        85.009391477852368246090463829899_dp, 8.0e-15_dp, 'log Phi(0.1,2;100)')
    call check_close(log_wright_bessel(1.0_dp, 1.0_dp, 1.0e6_dp), &
        1995.2806727526574304525311176014_dp, 8.0e-16_dp, 'large log Phi')
    call require(.not. ieee_is_finite(wright_bessel(1.0_dp, 1.0_dp, 1.0e6_dp)) .and. &
        wright_bessel(1.0_dp, 1.0_dp, 1.0e6_dp) > 0.0_dp, 'large Phi overflow')
    call check_close(log_wright_bessel(1.0_dp, 1.0_dp, 1.0e24_dp), &
        1.9999999999849189773185511430e12_dp, 6.0e-15_dp, 'saddle log Phi')

    lhs = wright_bessel(1.5_dp, 0.1_dp, 2.5_dp)
    rhs = 1.5_dp * 2.5_dp * wright_bessel(1.5_dp, 2.6_dp, 2.5_dp) + &
        0.1_dp * wright_bessel(1.5_dp, 1.1_dp, 2.5_dp)
    call check_close(lhs, rhs, 8.0e-15_dp, 'Wright recurrence')

    a_values = [1.0_dp, 1.5_dp, 2.0_dp]
    b_values = [1.0_dp, 1.1_dp, 3.0_dp]
    x_values = [1.0_dp, 2.5_dp, 4.0_dp]
    values = wright_bessel(a_values, b_values, x_values)
    call check_close(values(1), 2.2795853023360672674372044408115_dp, 3.0e-15_dp, &
        'elemental Wright 1')
    call check_close(values(2), 3.3027831556598564651664788589047_dp, 4.0e-15_dp, &
        'elemental Wright 2')
    call check_close(values(3), 0.67804528536995598609896143657115_dp, 4.0e-15_dp, &
        'elemental Wright 3')

    call require(ieee_is_nan(wright_bessel(-1.0_dp, 1.0_dp, 1.0_dp)), 'negative a')
    call require(ieee_is_nan(wright_bessel(1.0_dp, -1.0_dp, 1.0_dp)), 'negative b')
    call require(ieee_is_nan(wright_bessel(1.0_dp, 1.0_dp, -1.0_dp)), 'negative x')
    call require(ieee_is_nan(wright_bessel(positive_infinity(1.0_dp), 1.0_dp, 1.0_dp)), &
        'infinite a')
    call require(ieee_is_nan(log_wright_bessel(1.0_dp, positive_infinity(1.0_dp), 1.0_dp)), &
        'infinite b')
    call require(.not. ieee_is_finite(log_wright_bessel(1.0_dp, 1.0_dp, &
        positive_infinity(1.0_dp))) .and. log_wright_bessel(1.0_dp, 1.0_dp, &
        positive_infinity(1.0_dp)) > 0.0_dp, 'infinite x')
    call require(ieee_is_nan(wright_bessel(quiet_nan(1.0_dp), 1.0_dp, 1.0_dp)), 'NaN a')

    print *, 'test_special_wright_bessel: PASS'

contains

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine check_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: relative_tolerance !! allowed relative error
        character(len=*), intent(in) :: label !! failure label
        real(dp) :: scale

        scale = max(1.0_dp, abs(expected))
        if (abs(actual - expected) > relative_tolerance * scale) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

end program test_special_wright_bessel
