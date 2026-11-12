program test_special_parabolic_cylinder
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_special, only : dp, parabolic_cylinder_status_invalid_argument, &
        parabolic_cylinder_status_success, pbdn_seq, pbdv, pbdv_seq, pbvv, pbvv_seq
    implicit none

    real(dp) :: d
    real(dp) :: dpdx
    real(dp) :: d_values(2)
    real(dp) :: dp_values(2)
    real(dp) :: pinf
    real(dp) :: qnan
    real(dp) :: vv
    real(dp) :: vvprime
    real(dp), allocatable :: sequence_values(:)
    real(dp), allocatable :: sequence_derivatives(:)
    complex(dp), allocatable :: complex_values(:)
    complex(dp), allocatable :: complex_derivatives(:)
    integer :: status

    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call pbdv(0.5_dp, 1.0_dp, d, dpdx)
    call check_rel('pbdv ordinary D', d, 0.84220324406983957449_dp, 3.0e-12_dp)
    call check_rel('pbdv ordinary derivative', dpdx, -0.094565608685238832651_dp, 3.0e-12_dp)

    call pbdv(-0.5_dp, -3.0_dp, d, dpdx)
    call check_rel('pbdv negative x D', d, 8.2111204276138111619_dp, 3.0e-12_dp)
    call check_rel('pbdv negative x derivative', dpdx, -10.548825240696614841_dp, 3.0e-12_dp)

    call pbdv(2.3_dp, -5.0_dp, d, dpdx)
    call check_rel('pbdv continued D', d, -6.4418538795883846696_dp, 3.0e-11_dp)
    call check_rel('pbdv continued derivative', dpdx, 10.580992692784509802_dp, 3.0e-11_dp)

    call pbdv(2.3_dp, 10.0_dp, d, dpdx)
    call check_rel('pbdv positive hyperu D', d, 2.7295609304685681295e-9_dp, 2.0e-12_dp)
    call check_rel('pbdv positive hyperu derivative', dpdx, -1.3011711718599137396e-8_dp, 2.0e-12_dp)

    call pbdv(10.2_dp, -10.0_dp, d, dpdx)
    call check_rel('pbdv difficult negative D', d, -2751706.0426883262141_dp, 5.0e-9_dp)
    call check_rel('pbdv difficult negative derivative', dpdx, 10152833.090600011222_dp, 5.0e-9_dp)

    call pbdv(20.2_dp, -10.0_dp, d, dpdx)
    call check_rel('pbdv high-order negative D', d, -2710614206.0014432894_dp, 2.0e-8_dp)
    call check_rel('pbdv high-order negative derivative', dpdx, 5174668166.2047478292_dp, 2.0e-8_dp)

    call pbdv(10.0_dp, -10.0_dp, d, dpdx)
    call check_rel('pbdv integer Hermite D', d, 0.084702174589852322705_dp, 3.0e-12_dp)
    call check_rel('pbdv integer Hermite derivative', dpdx, 0.32955206511453878950_dp, 3.0e-12_dp)

    call pbdv(3.5_dp, 20.0_dp, d, dpdx)
    call check_rel('pbdv large positive D', d, 1.3163845638926592024e-39_dp, 1.0e-12_dp)
    call check_rel('pbdv large positive derivative', dpdx, -1.2932023994352411400e-38_dp, 1.0e-12_dp)

    call pbdv(2.3_dp, 0.0_dp, d, dpdx)
    call check_rel('pbdv origin D', d, -1.0041285564344087235_dp, 3.0e-14_dp)
    call check_rel('pbdv origin derivative', dpdx, -0.86252135743456201982_dp, 3.0e-14_dp)

    call pbdv(0.5_dp, [1.0_dp, -3.0_dp], d_values, dp_values)
    call check_rel('pbdv elemental D first', d_values(1), 0.84220324406983957449_dp, 3.0e-12_dp)
    call check_rel('pbdv elemental derivative first', dp_values(1), -0.094565608685238832651_dp, 3.0e-12_dp)

    call pbdv(1.0_dp, pinf, d, dpdx)
    call check_rel('pbdv positive infinity D', d, 0.0_dp, 0.0_dp)
    call check_rel('pbdv positive infinity derivative', dpdx, -0.0_dp, 0.0_dp)

    call pbdv(qnan, 1.0_dp, d, dpdx)
    call check_nan('pbdv NaN parameter D', d)
    call check_nan('pbdv NaN parameter derivative', dpdx)
    call pbdv(1.0_dp, qnan, d, dpdx)
    call check_nan('pbdv NaN argument D', d)
    call check_nan('pbdv NaN argument derivative', dpdx)

    call pbvv(0.5_dp, 1.0_dp, vv, vvprime)
    call check_rel('pbvv ordinary V', vv, 0.22003508652565433385_dp, 3.0e-12_dp)
    call check_rel('pbvv ordinary derivative', vvprime, 0.92267135561997275506_dp, 3.0e-12_dp)

    call pbvv(2.3_dp, -5.0_dp, vv, vvprime)
    call check_rel('pbvv negative argument V', vv, -1.7780056207885620114_dp, 5.0e-11_dp)
    call check_rel('pbvv negative argument derivative', vvprime, 2.7965831354190559820_dp, 5.0e-11_dp)

    call pbvv(1.0_dp, -8.0_dp, vv, vvprime)
    call check_rel('pbvv integer negative V', vv, 116433.27222791884993_dp, 5.0e-12_dp)
    call check_rel('pbvv integer negative derivative', vvprime, -435082.36885445760913_dp, 5.0e-12_dp)

    call pbvv(8.2_dp, 15.0_dp, vv, vvprime)
    call check_rel('pbvv asymptotic V', vv, 4.0386386194065131020e13_dp, 2.0e-12_dp)
    call check_rel('pbvv asymptotic derivative', vvprime, 2.7688242356280221328e14_dp, 2.0e-12_dp)

    call pbvv(2.3_dp, 0.0_dp, vv, vvprime)
    call check_rel('pbvv origin V', vv, 0.19066181961531129939_dp, 3.0e-14_dp)
    call check_rel('pbvv origin derivative', vvprime, -0.63083025104535043250_dp, 3.0e-14_dp)

    call pbvv(qnan, 1.0_dp, vv, vvprime)
    call check_nan('pbvv NaN parameter V', vv)
    call check_nan('pbvv NaN parameter derivative', vvprime)

    call pbdv_seq(2.3_dp, 0.7_dp, sequence_values, sequence_derivatives, status)
    call require('pbdv_seq status', status == parabolic_cylinder_status_success)
    call require('pbdv_seq size', size(sequence_values) == 3 .and. size(sequence_derivatives) == 3)
    call check_rel('pbdv_seq D0', sequence_values(1), 0.87586419772839181410_dp, 4.0e-12_dp)
    call check_rel('pbdv_seq D1', sequence_values(2), 0.38300855638087207671_dp, 4.0e-12_dp)
    call check_rel('pbdv_seq D2', sequence_values(3), -0.87051746758029890462_dp, 4.0e-12_dp)
    call check_rel('pbdv_seq derivative D2', sequence_derivatives(3), 1.1856007933291103931_dp, 4.0e-12_dp)

    call pbdv_seq(-2.3_dp, 0.7_dp, sequence_values, sequence_derivatives, status)
    call require('pbdv_seq negative compatibility size', size(sequence_values) == 2)
    call check_rel('pbdv_seq negative compatibility D0', sequence_values(1), &
        0.77616256997626906052_dp, 4.0e-12_dp)
    call check_rel('pbdv_seq negative compatibility derivative D1', sequence_derivatives(2), &
        1.3381239450829721407_dp, 4.0e-12_dp)

    call pbvv_seq(2.3_dp, 0.7_dp, sequence_values, sequence_derivatives, status)
    call require('pbvv_seq status', status == parabolic_cylinder_status_success)
    call require('pbvv_seq size', size(sequence_values) == 3 .and. size(sequence_derivatives) == 3)
    call check_rel('pbvv_seq V0', sequence_values(1), 0.18159183536645170018_dp, 4.0e-12_dp)
    call check_rel('pbvv_seq V1', sequence_values(2), -0.63966127658882544861_dp, 4.0e-12_dp)
    call check_rel('pbvv_seq V2', sequence_values(3), -0.27363249086027370183_dp, 4.0e-12_dp)
    call check_rel('pbvv_seq derivative V2', sequence_derivatives(3), -0.54388990478772965297_dp, 4.0e-12_dp)

    call pbdn_seq(4, cmplx(1.0_dp, 0.5_dp, dp), complex_values, complex_derivatives, status)
    call require('pbdn_seq status', status == parabolic_cylinder_status_success)
    call require('pbdn_seq size', size(complex_values) == 5 .and. size(complex_derivatives) == 5)
    call check_complex_rel('pbdn_seq D0', complex_values(1), &
        cmplx(0.80325661056481202984_dp, -0.20510508617511684095_dp, dp), 3.0e-14_dp)
    call check_complex_rel('pbdn_seq D4', complex_values(5), &
        cmplx(-2.4792825707573490921_dp, -3.2172636430773652549_dp, dp), 3.0e-14_dp)
    call check_complex_rel('pbdn_seq derivative D4', complex_derivatives(5), &
        cmplx(-8.5030498846911575164_dp, 4.0829801068718991360_dp, dp), 3.0e-14_dp)

    call pbdn_seq(0, cmplx(1.0_dp, 0.5_dp, dp), complex_values, complex_derivatives, status)
    call require('pbdn_seq zero-order compatibility size', size(complex_values) == 2)

    call pbdn_seq(-1, cmplx(1.0_dp, 0.5_dp, dp), complex_values, complex_derivatives, status)
    call require('pbdn_seq negative status', status == parabolic_cylinder_status_invalid_argument)
    call require('pbdn_seq negative size', size(complex_values) == 0 .and. size(complex_derivatives) == 0)

    print *, 'parabolic-cylinder function and sequence tests passed'

contains

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed finite value
        real(dp), intent(in) :: expected !! finite reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance

        if (.not. ieee_is_finite(actual) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_complex_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed finite value
        complex(dp), intent(in) :: expected !! finite reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance

        if (.not. ieee_is_finite(real(actual, dp)) .or. .not. ieee_is_finite(aimag(actual)) .or. &
                abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_rel

    subroutine require(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine require

    subroutine check_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_nan

end program test_special_parabolic_cylinder
