! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_airy_complete
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : ai_zeros, airy_status_invalid_argument, &
        airy_status_success, bi_zeros, dp, itairy
    implicit none

    real(dp), parameter :: x_values(4) = [1.0_dp, -1.0_dp, 8.0_dp, 20.0_dp]
    real(dp), parameter :: apt_ref(4) = [ &
        0.23631734191710978_dp, -0.46567398346706861_dp, &
        0.33333331724248357_dp, 0.33333333333333333_dp]
    real(dp), parameter :: bpt_ref(4) = [ &
        0.87276911673800812_dp, -0.37300500963429492_dp, &
        440065.25804904178_dp, 4.7447423034305476e24_dp]
    real(dp), parameter :: ant_ref(4) = [ &
        0.46567398346706861_dp, -0.23631734191710978_dp, &
        0.78398259657117734_dp, 0.71173925263991846_dp]
    real(dp), parameter :: bnt_ref(4) = [ &
        0.37300500963429492_dp, -0.87276911673800812_dp, &
        -0.014756446293227662_dp, -0.039061733297866577_dp]
    real(dp), parameter :: ai_zero_ref(5) = [ &
        -2.3381074104597670_dp, -4.0879494441309706_dp, &
        -5.5205598280955511_dp, -6.7867080900717589_dp, &
        -7.9441335871208532_dp]
    real(dp), parameter :: aip_zero_ref(5) = [ &
        -1.0187929716474711_dp, -3.2481975821798366_dp, &
        -4.8200992111787356_dp, -6.1633073556394865_dp, &
        -7.3721772550477702_dp]
    real(dp), parameter :: ai_at_aip_ref(5) = [ &
        0.53565665601569980_dp, -0.41901547803256395_dp, &
        0.38040646862815329_dp, -0.35790794371229230_dp, &
        0.34230124441162385_dp]
    real(dp), parameter :: aip_at_ai_ref(5) = [ &
        0.70121082272069136_dp, -0.80311136965486395_dp, &
        0.86520402589415193_dp, -0.91085073704960201_dp, &
        0.94733570944156777_dp]
    real(dp), parameter :: bi_zero_ref(5) = [ &
        -1.1737132227091279_dp, -3.2710933028363527_dp, &
        -4.8307378416620159_dp, -6.1698521283102513_dp, &
        -7.3767620793677637_dp]
    real(dp), parameter :: bip_zero_ref(5) = [ &
        -2.2944396826141232_dp, -4.0731550890718282_dp, &
        -5.5123957296635995_dp, -6.7812944459903054_dp, &
        -7.9401786891685780_dp]
    real(dp), parameter :: bi_at_bip_ref(5) = [ &
        -0.45494438363965740_dp, 0.39652283609446463_dp, &
        -0.36796916148695901_dp, 0.34949911683180427_dp, &
        -0.33602624013365946_dp]
    real(dp), parameter :: bip_at_bi_ref(5) = [ &
        0.60195788797623942_dp, -0.76031014149280110_dp, &
        0.83699101261926122_dp, -0.88947990142653920_dp, &
        0.92998363856802261_dp]

    real(dp), allocatable :: a(:)
    real(dp), allocatable :: ai_values(:)
    real(dp), allocatable :: aip_values(:)
    real(dp), allocatable :: ap(:)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: bi_values(:)
    real(dp), allocatable :: bip_values(:)
    real(dp), allocatable :: bp(:)
    real(dp) :: ant
    real(dp) :: ant_values(4)
    real(dp) :: apt
    real(dp) :: apt_values(4)
    real(dp) :: bnt
    real(dp) :: bnt_values(4)
    real(dp) :: bpt
    real(dp) :: bpt_values(4)
    real(dp) :: nan_value
    real(dp) :: pos_inf
    integer :: status

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call itairy(5.0_dp, apt, bpt, ant, bnt)
    call check_close('itairy Apt(5)', apt, 0.33328759030591788_dp, 4.0e-13_dp)
    call check_rel('itairy Bpt(5)', bpt, 321.47831857046517_dp, 5.0e-15_dp)
    call check_close('itairy Ant(5)', ant, 0.71788220454782765_dp, 4.0e-13_dp)
    call check_close('itairy Bnt(5)', bnt, 0.15873093858143865_dp, 4.0e-13_dp)

    call itairy(10.0_dp, apt, bpt, ant, bnt)
    call check_close('itairy Apt(10)', apt, 0.33333333329916902_dp, 4.0e-13_dp)
    call check_rel('itairy Bpt(10)', bpt, 147809803.10740672_dp, 5.0e-15_dp)
    call check_close('itairy Ant(10)', ant, 0.76569840313421292_dp, 4.0e-13_dp)
    call check_close('itairy Bnt(10)', bnt, 0.015040424806140020_dp, 4.0e-13_dp)

    call itairy(100.0_dp, apt, bpt, ant, bnt)
    call check_close('itairy Apt(100)', apt, one_third_ref(), 4.0e-15_dp)
    call check_rel('itairy Bpt(100)', bpt, 6.0457660973056185e287_dp, 3.0e-13_dp)
    call check_close('itairy Ant(100)', ant, 0.66422602599797836_dp, 2.0e-12_dp)
    call check_close('itairy Bnt(100)', bnt, 0.017673486211800158_dp, 2.0e-12_dp)

    call itairy(x_values, apt_values, bpt_values, ant_values, bnt_values)
    call check_array('itairy elemental Apt', apt_values, apt_ref, 5.0e-13_dp)
    call check_array_rel('itairy elemental Bpt', bpt_values, bpt_ref, 5.0e-13_dp)
    call check_array('itairy elemental Ant', ant_values, ant_ref, 5.0e-13_dp)
    call check_array('itairy elemental Bnt', bnt_values, bnt_ref, 5.0e-13_dp)

    call itairy(pos_inf, apt, bpt, ant, bnt)
    call check_close('itairy +inf Apt', apt, one_third_ref(), 0.0_dp)
    call check_true('itairy +inf Bpt', bpt > huge(1.0_dp))
    call check_close('itairy +inf Ant', ant, 2.0_dp / 3.0_dp, 0.0_dp)
    call check_close('itairy +inf Bnt', bnt, 0.0_dp, 0.0_dp)
    call itairy(-pos_inf, apt, bpt, ant, bnt)
    call check_close('itairy -inf Apt', apt, -2.0_dp / 3.0_dp, 0.0_dp)
    call check_close('itairy -inf Bpt', bpt, -0.0_dp, 0.0_dp)
    call check_close('itairy -inf Ant', ant, -one_third_ref(), 0.0_dp)
    call check_true('itairy -inf Bnt', bnt < -huge(1.0_dp))
    call itairy(nan_value, apt, bpt, ant, bnt)
    call check_true('itairy NaN Apt', ieee_is_nan(apt))
    call check_true('itairy NaN Bpt', ieee_is_nan(bpt))
    call check_true('itairy NaN Ant', ieee_is_nan(ant))
    call check_true('itairy NaN Bnt', ieee_is_nan(bnt))

    call ai_zeros(5, a, ap, ai_values, aip_values, status)
    call check_true('ai_zeros status', status == airy_status_success)
    call check_array('ai_zeros roots', a, ai_zero_ref, 3.0e-12_dp)
    call check_array('ai_zeros derivative roots', ap, aip_zero_ref, 3.0e-12_dp)
    call check_array('ai_zeros Ai(ap)', ai_values, ai_at_aip_ref, 3.0e-12_dp)
    call check_array('ai_zeros Aip(a)', aip_values, aip_at_ai_ref, 3.0e-12_dp)
    call check_true('ai_zeros ordering', all(a(2:) < a(:4)) .and. all(ap(2:) < ap(:4)))

    call bi_zeros(5, b, bp, bi_values, bip_values, status)
    call check_true('bi_zeros status', status == airy_status_success)
    call check_array('bi_zeros roots', b, bi_zero_ref, 3.0e-12_dp)
    call check_array('bi_zeros derivative roots', bp, bip_zero_ref, 3.0e-12_dp)
    call check_array('bi_zeros Bi(bp)', bi_values, bi_at_bip_ref, 3.0e-12_dp)
    call check_array('bi_zeros Bip(b)', bip_values, bip_at_bi_ref, 3.0e-12_dp)
    call check_true('bi_zeros ordering', all(b(2:) < b(:4)) .and. all(bp(2:) < bp(:4)))

    call ai_zeros(0, a, ap, ai_values, aip_values, status)
    call check_true('ai_zeros invalid status', status == airy_status_invalid_argument)
    call check_true('ai_zeros invalid arrays', size(a) == 0 .and. size(ap) == 0 .and. &
        size(ai_values) == 0 .and. size(aip_values) == 0)

    print *, 'Airy integral/zero tests passed'

contains

    pure function one_third_ref() result(value)
        real(dp) :: value
        value = 1.0_dp / 3.0_dp
    end function one_third_ref

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
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum absolute elementwise error

        if (size(actual) /= size(expected) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_array_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! nonzero reference values
        real(dp), intent(in) :: tolerance !! maximum relative elementwise error

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

end program test_special_airy_complete
