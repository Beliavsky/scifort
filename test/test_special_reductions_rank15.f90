! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank15
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : log_softmax, logsumexp, &
        logsumexp_keepdims, logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    real(dp) :: x7(2, 1, 1, 1, 1, 1, 2)
    real(dp) :: x8(2, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x9(2, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x10(2, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x11(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x12(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x13(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x14(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x15(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3)
    complex(dp) :: z15(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3)
    real(dp), allocatable :: fourteen(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: fifteen(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: signs(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    complex(dp), allocatable :: cfourteen(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    complex(dp), allocatable :: cfifteen(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    complex(dp), allocatable :: phases(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: empty15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    complex(dp) :: cvalue
    complex(dp) :: phase
    integer :: i

    x7 = reshape([(real(i, dp), i = 1, size(x7))], shape(x7))
    x8 = reshape([(real(i, dp), i = 1, size(x8))], shape(x8))
    x9 = reshape([(real(i, dp), i = 1, size(x9))], shape(x9))
    x10 = reshape([(real(i, dp), i = 1, size(x10))], shape(x10))
    x11 = reshape([(real(i, dp), i = 1, size(x11))], shape(x11))
    x12 = reshape([(real(i, dp), i = 1, size(x12))], shape(x12))
    x13 = reshape([(real(i, dp), i = 1, size(x13))], shape(x13))
    x14 = reshape([(real(i, dp), i = 1, size(x14))], shape(x14))
    x15 = reshape([(real(i, dp), i = 1, size(x15))], shape(x15))
    z15 = reshape([(cmplx(real(i, dp), 0.07_dp * real(i, dp) * &
        merge(-1.0_dp, 1.0_dp, mod(i, 2) == 1), kind=dp), i = 1, size(z15))], shape(z15))

    call check_close('rank7 all', logsumexp(x7), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank7 softmax', sum(softmax(x7)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank8 all', logsumexp(x8), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank8 softmax', sum(softmax(x8)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank9 all', logsumexp(x9), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank9 softmax', sum(softmax(x9)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank10 all', logsumexp(x10), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank10 softmax', sum(softmax(x10)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank11 all', logsumexp(x11), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank11 softmax', sum(softmax(x11)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank12 all', logsumexp(x12), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank12 softmax', sum(softmax(x12)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank13 all', logsumexp(x13), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank13 softmax', sum(softmax(x13)), 1.0_dp, 4.0e-15_dp)
    call check_close('rank14 all', logsumexp(x14), 4.440189698561196_dp, 5.0e-14_dp)
    call check_close('rank14 softmax', sum(softmax(x14)), 1.0_dp, 4.0e-15_dp)

    call check_close('rank15 all', logsumexp(x15), 24.458675145349332_dp, 8.0e-14_dp)

    fourteen = logsumexp(x15, 1)
    call check_shape('rank15 dim1 shape', shape(fourteen), &
        [2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3])
    call check_close('rank15 dim1 first', fourteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        2.313261687518223_dp, 3.0e-14_dp)
    call check_close('rank15 dim1 last', fourteen(2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3), &
        24.313261687518224_dp, 8.0e-14_dp)

    fourteen = logsumexp(x15, 15)
    call check_close('rank15 dim15 first', fourteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        17.000335518870326_dp, 8.0e-14_dp)
    call check_close('rank15 dim15 last', fourteen(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        24.000335518870326_dp, 8.0e-14_dp)

    fifteen = logsumexp_keepdims(x15, [1, 4, 15])
    call check_shape('rank15 tuple shape', shape(fifteen), &
        [1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    call check_close('rank15 tuple first', fifteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        22.331747134306358_dp, 8.0e-14_dp)
    call check_close('rank15 tuple second', fifteen(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        24.331747134306358_dp, 8.0e-14_dp)
    call check_rank15_real('rank15 reversed tuple', &
        logsumexp_keepdims(x15, [15, 4, 1]), fifteen, 4.0e-15_dp)

    call logsumexp_signed(x15, logabs, sign_value)
    call check_close('rank15 signed all logabs', logabs, 24.458675145349332_dp, 8.0e-14_dp)
    call check_close('rank15 signed all sign', sign_value, 1.0_dp, 0.0_dp)

    call logsumexp_signed_keepdims(x15, [1, 4, 15], fifteen, signs)
    call check_close('rank15 signed tuple logabs', &
        fifteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        22.331747134306358_dp, 8.0e-14_dp)
    call check_true('rank15 signed tuple signs', all(signs > 0.0_dp))

    fifteen = softmax(x15, [1, 4, 15])
    call check_close('rank15 tuple softmax last', &
        fifteen(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3), &
        0.7176687737243994_dp, 6.0e-15_dp)
    call check_close('rank15 tuple softmax group1', &
        sum(fifteen(:, 1, :, :, :, :, :, :, :, :, :, :, :, :, :)), 1.0_dp, 5.0e-15_dp)
    call check_close('rank15 tuple softmax group2', &
        sum(fifteen(:, 2, :, :, :, :, :, :, :, :, :, :, :, :, :)), 1.0_dp, 5.0e-15_dp)

    fifteen = log_softmax(x15, [2, 4, 15])
    call check_close('rank15 tuple log_softmax first', &
        fifteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        -22.145413457831108_dp, 8.0e-14_dp)
    call check_close('rank15 tuple log_softmax last', &
        fifteen(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3), &
        -0.1454134578311077_dp, 8.0e-14_dp)

    fifteen = softmax(x15, [integer ::])
    call check_close('rank15 empty softmax', maxval(abs(fifteen - 1.0_dp)), 0.0_dp, 0.0_dp)
    fifteen = log_softmax(x15, [integer ::])
    call check_close('rank15 empty log_softmax', maxval(abs(fifteen)), 0.0_dp, 0.0_dp)

    allocate(empty15(2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
    fifteen = logsumexp_keepdims(empty15, [2])
    call check_shape('rank15 zero-extent shape', shape(fifteen), &
        [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    call check_true('rank15 zero-extent value', all(fifteen < -huge(1.0_dp)))
    fifteen = logsumexp_keepdims(x15, [1, 1])
    call check_shape('rank15 duplicate axes rejected', shape(fifteen), &
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    cvalue = logsumexp(z15)
    call check_complex_close('complex rank15 all', cvalue, &
        cmplx(23.689996230271078_dp, 1.7188068674562635_dp, dp), 1.2e-13_dp)
    call logsumexp_signed(z15, logabs, phase)
    call check_close('complex rank15 signed logabs', logabs, 23.689996230271078_dp, 1.2e-13_dp)
    call check_complex_close('complex rank15 signed phase', phase, &
        exp(cmplx(0.0_dp, 1.7188068674562635_dp, dp)), 1.0e-14_dp)

    call logsumexp_signed_keepdims(z15, [1, 4, 15], fifteen, phases)
    call check_close('complex rank15 signed tuple logabs', &
        fifteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        21.56823936459147_dp, 1.2e-13_dp)
    call check_complex_close('complex rank15 signed tuple phase', &
        phases(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        exp(cmplx(0.0_dp, 1.453204197176173_dp, dp)), 1.0e-14_dp)

    cfourteen = logsumexp(z15, 15)
    call check_complex_close('complex rank15 dim15 first', &
        cfourteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        cmplx(17.000284246943586_dp, -1.189821756260675_dp, dp), 1.2e-13_dp)
    call check_complex_close('complex rank15 dim15 last', &
        cfourteen(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        cmplx(24.000284246943586_dp, 1.679821756260675_dp, dp), 1.2e-13_dp)

    cfifteen = logsumexp_keepdims(z15, [1, 4, 15])
    call check_complex_close('complex rank15 tuple first', &
        cfifteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        cmplx(21.56823936459147_dp, 1.453204197176173_dp, dp), 1.2e-13_dp)
    call check_complex_close('complex rank15 tuple second', &
        cfifteen(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        cmplx(23.567760510611553_dp, 1.7543556432066514_dp, dp), 1.2e-13_dp)

    cfifteen = softmax(z15, [1, 4, 15])
    call check_complex_close('complex rank15 softmax last', &
        cfifteen(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3), &
        cmplx(1.5364469224209218_dp, -0.11445450752107039_dp, dp), 2.0e-14_dp)

    cfifteen = log_softmax(z15, [2, 4, 15])
    call check_complex_close('complex rank15 log_softmax first', &
        cfifteen(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        cmplx(-22.143645531615036_dp, 1.5181958921239103_dp, dp), 1.2e-13_dp)
    call check_complex_close('complex rank15 log_softmax last', &
        cfifteen(2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3), &
        cmplx(-0.14364553161503452_dp, 0.021804107876089714_dp, dp), 1.2e-13_dp)

    print *, 'rank-7 through rank-15 special reduction tests passed'

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

    subroutine check_complex_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum complex absolute error

        if (actual == expected) return
        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_shape(name, actual, expected)
        character(len=*), intent(in) :: name !! check label printed on failure
        integer, intent(in) :: actual(:) !! computed extents
        integer, intent(in) :: expected(:) !! expected extents

        if (size(actual) /= size(expected) .or. any(actual /= expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_shape

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition that must hold

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

    subroutine check_rank15_real(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! result
        real(dp), intent(in) :: expected(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! reference
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_rank15_real

    pure function complex_is_nan(value) result(is_nan)
        complex(dp), intent(in) :: value !! value inspected for a NaN component
        logical :: is_nan

        is_nan = ieee_is_nan(real(value, dp)) .or. ieee_is_nan(aimag(value))
    end function complex_is_nan

end program test_special_reductions_rank15
