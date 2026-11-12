! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_spheroidal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_math, only : positive_infinity
    use scifort_special, only : dp, obl_cv, obl_cv_seq, pro_cv, pro_cv_seq, &
        spheroidal_status_invalid_argument, spheroidal_status_success
    implicit none

    real(dp), allocatable :: sequence(:)
    real(dp) :: values(4)
    integer :: status

    call check_close(pro_cv(0, 0, 0.0_dp), 0.0_dp, 0.0_dp, 'pro c0')
    call check_close(pro_cv(0, 0, 1.0_dp), 0.31900005514689334_dp, 3.0e-13_dp, &
        'pro 0 0 1')
    call check_close(pro_cv(1, 1, 1.0_dp), 2.1955483554130018_dp, 3.0e-13_dp, &
        'pro 1 1 1')
    call check_close(pro_cv(0, 3, 2.0_dp), 14.100203876205317_dp, 3.0e-13_dp, &
        'pro odd family')
    call check_close(pro_cv(2, 4, sqrt(10.0_dp)), 23.979073449847263_dp, &
        4.0e-13_dp, 'pro DLMF example')
    call check_close(pro_cv(3, 7, 5.0_dp), 66.5235949585232_dp, 4.0e-13_dp, &
        'pro higher mode')
    call check_close(pro_cv(1, 1, 50.0_dp), 50.25646486503018_dp, 8.0e-13_dp, &
        'pro large c')

    call check_close(obl_cv(0, 0, 1.0_dp), -0.34860239947026983_dp, 3.0e-13_dp, &
        'obl 0 0 1')
    call check_close(obl_cv(1, 1, 1.0_dp), 1.7953045872818243_dp, 3.0e-13_dp, &
        'obl 1 1 1')
    call check_close(obl_cv(0, 3, 2.0_dp), 10.003863804730251_dp, 3.0e-13_dp, &
        'obl odd family')
    call check_close(obl_cv(2, 4, sqrt(10.0_dp)), 15.990084610139247_dp, &
        4.0e-13_dp, 'obl DLMF scale')
    call check_close(obl_cv(5, 5, 20.0_dp), -166.55392690574655_dp, 8.0e-13_dp, &
        'obl large c')
    call check_close(obl_cv(5, 8, 20.0_dp), -102.71411983969165_dp, 8.0e-13_dp, &
        'obl large higher mode')

    call check_close(pro_cv(2, 4, -3.0_dp), pro_cv(2, 4, 3.0_dp), 0.0_dp, &
        'pro even in c')
    call check_close(obl_cv(2, 4, -3.0_dp), obl_cv(2, 4, 3.0_dp), 0.0_dp, &
        'obl even in c')
    call require(pro_cv(2, 4, 3.0_dp) < pro_cv(2, 5, 3.0_dp), &
        'pro degree ordering')
    call require(obl_cv(2, 4, 3.0_dp) < obl_cv(2, 5, 3.0_dp), &
        'obl degree ordering')

    values = pro_cv([0, 1, 2, 3], [0, 1, 4, 7], [1.0_dp, 1.0_dp, 2.0_dp, 5.0_dp])
    call check_close(values(1), 0.31900005514689334_dp, 3.0e-13_dp, &
        'pro elemental 1')
    call check_close(values(4), 66.5235949585232_dp, 4.0e-13_dp, &
        'pro elemental 4')

    call pro_cv_seq(0, 3, 1.0_dp, sequence, status)
    call require(status == spheroidal_status_success .and. size(sequence) == 4, &
        'pro sequence status')
    call check_close(sequence(1), 0.31900005514689334_dp, 3.0e-13_dp, &
        'pro sequence first')
    call check_close(sequence(4), 12.51446214509403_dp, 4.0e-13_dp, &
        'pro sequence last')
    call obl_cv_seq(2, 4, 3.0_dp, sequence, status)
    call require(status == spheroidal_status_success .and. size(sequence) == 3, &
        'obl sequence status')
    call check_close(sequence(1), 4.526460462194931_dp, 4.0e-13_dp, &
        'obl sequence first')
    call check_close(sequence(3), 16.386106019099607_dp, 4.0e-13_dp, &
        'obl sequence last')
    call pro_cv_seq(2, 1, 3.0_dp, sequence, status)
    call require(status == spheroidal_status_success .and. size(sequence) == 0, &
        'empty pro sequence')
    call obl_cv_seq(-1, 2, 1.0_dp, sequence, status)
    call require(status == spheroidal_status_invalid_argument .and. size(sequence) == 0, &
        'invalid obl sequence')

    call require(ieee_is_nan(pro_cv(-1, 2, 1.0_dp)), 'pro negative m')
    call require(ieee_is_nan(pro_cv(2, 1, 1.0_dp)), 'pro n below m')
    call require(ieee_is_nan(obl_cv(-1, 2, 1.0_dp)), 'obl negative m')
    call require(ieee_is_nan(obl_cv(2, 1, 1.0_dp)), 'obl n below m')
    call require(ieee_is_nan(pro_cv(2, 4, positive_infinity(1.0_dp))), &
        'pro infinite c')
    call require(ieee_is_nan(obl_cv(2, 4, positive_infinity(1.0_dp))), &
        'obl infinite c')
    call require(ieee_is_nan(pro_cv(0, 0, 2.0_dp * sqrt(sqrt(huge(1.0_dp))))), &
        'pro extreme finite c')

    print *, 'test_special_spheroidal: PASS'

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

end program test_special_spheroidal
