! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_mathieu
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_math, only : positive_infinity
    use scifort_special, only : dp, mathieu_a, mathieu_b, mathieu_cem, mathieu_even_coef, &
        mathieu_odd_coef, mathieu_sem, mathieu_status_invalid_argument, mathieu_status_success
    implicit none

    real(dp), allocatable :: coeff(:)
    real(dp) :: derivative
    real(dp) :: values(3)
    real(dp) :: derivatives(3)
    real(dp) :: value
    integer :: status

    call check_close(mathieu_a(0, 1.0_dp), -0.45513860410741364_dp, 2.0e-13_dp, 'a0 q1')
    call check_close(mathieu_a(1, 1.0_dp), 1.8591080725143634_dp, 2.0e-13_dp, 'a1 q1')
    call check_close(mathieu_a(2, 1.0_dp), 4.371300982735086_dp, 2.0e-13_dp, 'a2 q1')
    call check_close(mathieu_a(3, 10.0_dp), 15.502784369732575_dp, 2.0e-13_dp, 'a3 q10')
    call check_close(mathieu_a(8, 100.0_dp), 98.15039353420741_dp, 2.0e-13_dp, 'a8 q100')
    call check_close(mathieu_b(1, 1.0_dp), -0.11024881699209521_dp, 2.0e-13_dp, 'b1 q1')
    call check_close(mathieu_b(2, 1.0_dp), 3.917024772998471_dp, 2.0e-13_dp, 'b2 q1')
    call check_close(mathieu_b(3, 10.0_dp), 7.986069144681660_dp, 2.0e-13_dp, 'b3 q10')
    call check_close(mathieu_b(8, 100.0_dp), 68.10919350015939_dp, 2.0e-13_dp, 'b8 q100')

    call check_close(mathieu_a(7, 0.0_dp), 49.0_dp, 0.0_dp, 'a q0')
    call check_close(mathieu_b(7, 0.0_dp), 49.0_dp, 0.0_dp, 'b q0')
    call check_close(mathieu_a(1, -10.0_dp), mathieu_b(1, 10.0_dp), 0.0_dp, &
        'odd a q reflection')
    call check_close(mathieu_b(1, -10.0_dp), mathieu_a(1, 10.0_dp), 0.0_dp, &
        'odd b q reflection')
    call check_close(mathieu_a(6, -10.0_dp), mathieu_a(6, 10.0_dp), 0.0_dp, &
        'even a q reflection')
    call check_close(mathieu_b(6, -10.0_dp), mathieu_b(6, 10.0_dp), 0.0_dp, &
        'even b q reflection')

    call mathieu_cem(1, 1.0_dp, 30.0_dp, value, derivative)
    call check_close(value, 0.8524251395426153_dp, 2.0e-12_dp, 'cem 1 1 30 value')
    call check_close(derivative, -0.09208605025523009_dp, 8.0e-12_dp, &
        'cem 1 1 30 derivative')
    call mathieu_cem(2, 10.0_dp, 45.0_dp, value, derivative)
    call check_close(value, 1.0634128143472312_dp, 2.0e-12_dp, 'cem 2 10 45 value')
    call check_close(derivative, 0.9059099757066005_dp, 8.0e-12_dp, &
        'cem 2 10 45 derivative')
    call mathieu_cem(8, 100.0_dp, 12.3_dp, value, derivative)
    call check_close(value, 0.08228439355879553_dp, 3.0e-12_dp, 'cem 8 100 12.3 value')
    call check_close(derivative, 0.7587954696590246_dp, 8.0e-12_dp, &
        'cem 8 100 12.3 derivative')

    call mathieu_sem(1, 1.0_dp, 30.0_dp, value, derivative)
    call check_close(value, 0.38962739991246825_dp, 2.0e-12_dp, 'sem 1 1 30 value')
    call check_close(derivative, 0.8424282238210571_dp, 8.0e-12_dp, &
        'sem 1 1 30 derivative')
    call mathieu_sem(2, 10.0_dp, 45.0_dp, value, derivative)
    call check_close(value, 0.6916564591555207_dp, 2.0e-12_dp, 'sem 2 10 45 value')
    call check_close(derivative, 1.9539757458185834_dp, 8.0e-12_dp, &
        'sem 2 10 45 derivative')
    call mathieu_sem(8, 100.0_dp, 12.3_dp, value, derivative)
    call check_close(value, 0.028988782904512324_dp, 3.0e-12_dp, 'sem 8 100 12.3 value')
    call check_close(derivative, 0.32236796739741547_dp, 8.0e-12_dp, &
        'sem 8 100 12.3 derivative')

    call mathieu_cem(3, -10.0_dp, 30.0_dp, value, derivative)
    call check_close(value, -0.9227428982227223_dp, 2.0e-12_dp, 'cem negative q value')
    call check_close(derivative, -2.6391164400054654_dp, 8.0e-12_dp, &
        'cem negative q derivative')
    call mathieu_sem(3, -10.0_dp, 30.0_dp, value, derivative)
    call check_close(value, 0.11886541099608393_dp, 2.0e-12_dp, 'sem negative q value')
    call check_close(derivative, -3.9550530722656836_dp, 8.0e-12_dp, &
        'sem negative q derivative')

    call mathieu_cem(2, 0.0_dp, 30.0_dp, value, derivative)
    call check_close(value, 0.5_dp, 2.0e-15_dp, 'cem degree input q0')
    call check_close(derivative, -1.7320508075688772935_dp, 2.0e-15_dp, &
        'cem radian derivative q0')
    call mathieu_sem(2, 0.0_dp, 30.0_dp, value, derivative)
    call check_close(value, 0.8660254037844386468_dp, 2.0e-15_dp, 'sem degree input q0')
    call check_close(derivative, 1.0_dp, 2.0e-15_dp, 'sem radian derivative q0')

    call mathieu_even_coef(0, 1.0_dp, coeff, status)
    call require(status == mathieu_status_success .and. size(coeff) >= 6, 'even coef status')
    call check_close(coeff(1), 0.6729896723164999_dp, 2.0e-12_dp, 'even coef 0')
    call check_close(coeff(2), -0.30630358003683744_dp, 2.0e-12_dp, 'even coef 2')
    call check_close(coeff(3), 0.018645559365419891_dp, 2.0e-12_dp, 'even coef 4')
    call check_close(2.0_dp * coeff(1)**2 + sum(coeff(2:)**2), 1.0_dp, &
        3.0e-13_dp, 'even coefficient normalization')
    call check_close(mathieu_a(0, 1.0_dp) * coeff(1) - coeff(2), 0.0_dp, &
        4.0e-12_dp, 'even coefficient recurrence')

    call mathieu_odd_coef(1, 1.0_dp, coeff, status)
    call require(status == mathieu_status_success .and. size(coeff) >= 6, 'odd coef status')
    call check_close(coeff(1), 0.9939679613989341_dp, 2.0e-12_dp, 'odd coef 1')
    call check_close(coeff(2), -0.10958379187227703_dp, 2.0e-12_dp, 'odd coef 3')
    call check_close(sum(coeff**2), 1.0_dp, 3.0e-13_dp, 'odd coefficient normalization')

    values = mathieu_a([0, 1, 2], [1.0_dp, 1.0_dp, 1.0_dp])
    call check_close(values(1), -0.45513860410741364_dp, 2.0e-13_dp, 'a elemental 1')
    call check_close(values(3), 4.371300982735086_dp, 2.0e-13_dp, 'a elemental 3')
    call mathieu_cem([0, 1, 2], [0.0_dp, 1.0_dp, 1.0_dp], &
        [30.0_dp, 30.0_dp, 30.0_dp], values, derivatives)
    call check_close(values(1), 0.70710678118654752440_dp, 2.0e-15_dp, 'cem elemental 1')
    call check_close(values(2), 0.8524251395426153_dp, 2.0e-12_dp, 'cem elemental 2')

    call require(ieee_is_nan(mathieu_a(-1, 1.0_dp)), 'negative a order')
    call require(ieee_is_nan(mathieu_b(0, 1.0_dp)), 'zero b order')
    call require(ieee_is_nan(mathieu_a(0, positive_infinity(1.0_dp))), 'infinite q')
    call mathieu_even_coef(0, -1.0_dp, coeff, status)
    call require(status == mathieu_status_invalid_argument .and. size(coeff) == 0, &
        'negative coefficient q')
    call mathieu_odd_coef(0, 1.0_dp, coeff, status)
    call require(status == mathieu_status_invalid_argument .and. size(coeff) == 0, &
        'zero odd coefficient order')

    ! At q=300 SciPy 1.17 duplicates a_19 and a_21. The DLMF Fourier
    ! eigenproblem keeps the same-parity Sturm-Liouville ordering distinct.
    call require(mathieu_a(19, 300.0_dp) < mathieu_a(21, 300.0_dp), &
        'high-q same-parity ordering')

    print *, 'test_special_mathieu: PASS'

contains

    subroutine check_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! independent reference value
        real(dp), intent(in) :: relative_tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > relative_tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_special_mathieu
