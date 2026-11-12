! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank15_weighted
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims
    implicit none

    real(dp) :: x6(2, 1, 1, 1, 1, 2)
    real(dp) :: x7(2, 1, 1, 1, 1, 1, 2)
    real(dp) :: x8(2, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x9(2, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x10(2, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x11(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x12(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x13(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x14(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    complex(dp) :: z15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    real(dp) :: b6(3, 2)
    real(dp) :: b15(3, 1, 2)
    complex(dp) :: cb15(3, 1, 2)
    real(dp) :: bneg(2)
    real(dp) :: bad_b(3)
    real(dp) :: too_high_b(1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w6(1, 1, 1, 1, 1, 1)
    real(dp) :: w7(1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w8(1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w9(1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w10(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w11(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w12(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w13(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w14(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: w15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    complex(dp) :: cw15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: r14(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: r15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: s15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    complex(dp), allocatable :: c14(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    complex(dp), allocatable :: c15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: empty15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    complex(dp) :: cvalue
    complex(dp) :: phase

    x6 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x6))
    x7 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x7))
    x8 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x8))
    x9 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x9))
    x10 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x10))
    x11 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x11))
    x12 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x12))
    x13 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x13))
    x14 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x14))
    x15 = reshape([-0.4_dp, 0.2_dp, 0.7_dp, -0.1_dp], shape(x15))
    z15 = reshape([cmplx(-0.4_dp, 0.2_dp, dp), cmplx(0.2_dp, -0.7_dp, dp), &
        cmplx(0.7_dp, 1.1_dp, dp), cmplx(-0.1_dp, -0.3_dp, dp)], shape(z15))
    b6 = reshape([0.5_dp, 1.5_dp, 2.0_dp, 0.75_dp, 1.25_dp, 0.25_dp], shape(b6))
    b15 = reshape([0.5_dp, 1.5_dp, 2.0_dp, 0.75_dp, 1.25_dp, 0.25_dp], shape(b15))
    cb15 = reshape([cmplx(0.5_dp, 0.25_dp, dp), cmplx(1.5_dp, -0.5_dp, dp), &
        cmplx(2.0_dp, 0.75_dp, dp), cmplx(0.75_dp, -0.25_dp, dp), &
        cmplx(1.25_dp, 0.4_dp, dp), cmplx(0.25_dp, -0.1_dp, dp)], shape(cb15))
    bneg = [1.0_dp, -2.0_dp]
    bad_b = [1.0_dp, 2.0_dp, 3.0_dp]
    too_high_b = 1.0_dp

    w6 = 2.0_dp
    w7 = 2.0_dp
    w8 = 2.0_dp
    w9 = 2.0_dp
    w10 = 2.0_dp
    w11 = 2.0_dp
    w12 = 2.0_dp
    w13 = 2.0_dp
    w14 = 2.0_dp
    w15 = 2.0_dp
    cw15 = cmplx(2.0_dp, 0.5_dp, dp)

    call check_close('rank6 weighted broadcast all', logsumexp(x6, b6), &
        2.6485633608748156_dp, 6.0e-14_dp)
    call check_true('rank6 weight rank too high rejected', ieee_is_nan(logsumexp(x6, too_high_b)))

    call check_close('rank6 full-rank weight', logsumexp(x6, w6), &
        logsumexp(x6) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank7 full-rank weight', logsumexp(x7, w7), &
        logsumexp(x7) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank8 full-rank weight', logsumexp(x8, w8), &
        logsumexp(x8) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank9 full-rank weight', logsumexp(x9, w9), &
        logsumexp(x9) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank10 full-rank weight', logsumexp(x10, w10), &
        logsumexp(x10) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank11 full-rank weight', logsumexp(x11, w11), &
        logsumexp(x11) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank12 full-rank weight', logsumexp(x12, w12), &
        logsumexp(x12) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank13 full-rank weight', logsumexp(x13, w13), &
        logsumexp(x13) + log(2.0_dp), 4.0e-14_dp)
    call check_close('rank14 full-rank weight', logsumexp(x14, w14), &
        logsumexp(x14) + log(2.0_dp), 4.0e-14_dp)

    call check_close('rank15 weighted broadcast all', logsumexp(x15, b15), &
        2.6485633608748156_dp, 6.0e-14_dp)

    call check_close('rank15 full-rank real weight', logsumexp(x15, w15), &
        logsumexp(x15) + log(2.0_dp), 4.0e-14_dp)
    call check_complex_close('rank15 full-rank complex weight', logsumexp(x15, cw15), &
        cmplx(logsumexp(x15), 0.0_dp, dp) + log(cmplx(2.0_dp, 0.5_dp, dp)), 5.0e-14_dp)
    call check_complex_close('rank15 complex x full-rank real weight', logsumexp(z15, w15), &
        logsumexp(z15) + log(cmplx(2.0_dp, 0.0_dp, dp)), 7.0e-14_dp)
    call check_complex_close('rank15 complex x full-rank complex weight', logsumexp(z15, cw15), &
        logsumexp(z15) + log(cmplx(2.0_dp, 0.5_dp, dp)), 7.0e-14_dp)

    r14 = logsumexp(x15, 13, b15)
    call check_shape('rank15 weighted dim shape', shape(r14), &
        [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2])
    call check_real_array('rank15 weighted dim values', reshape(r14, [size(r14)]), &
        [0.9862943611198906_dp, 1.5862943611198906_dp, &
         1.5109302162163288_dp, 0.7109302162163288_dp], 7.0e-14_dp)

    r15 = logsumexp_keepdims(x15, [1, 13], b15)
    call check_shape('rank15 weighted tuple shape', shape(r15), &
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2])
    call check_real_array('rank15 weighted tuple values', reshape(r15, [size(r15)]), &
        [2.0237823116057765_dp, 1.8820308821641065_dp], 7.0e-14_dp)

    call logsumexp_signed(x15, logabs, sign_value, bneg)
    call check_close('rank15 signed weighted logabs', logabs, 1.3725649036760983_dp, 6.0e-14_dp)
    call check_close('rank15 signed weighted sign', sign_value, -1.0_dp, 0.0_dp)
    call check_true('rank15 unsigned negative weighted sum', ieee_is_nan(logsumexp(x15, bneg)))

    cvalue = logsumexp(x15, cb15)
    call check_complex_close('rank15 real x complex b', cvalue, &
        cmplx(2.6515380631720893_dp, 0.07709414827988947_dp, dp), 7.0e-14_dp)
    call logsumexp_signed(x15, logabs, phase, cb15)
    call check_close('rank15 real x complex b logabs', logabs, 2.6515380631720893_dp, 7.0e-14_dp)
    call check_complex_close('rank15 real x complex b phase', phase, &
        cmplx(0.9970297177457482_dp, 0.07701780269414214_dp, dp), 2.0e-14_dp)

    c14 = logsumexp(x15, 13, cb15)
    call check_complex_array('rank15 real x complex b dim', reshape(c14, [size(c14)]), &
        [cmplx(0.9940464543878732_dp, 0.12435499454676144_dp, dp), &
         cmplx(1.5940464543878732_dp, 0.12435499454676144_dp, dp), &
         cmplx(1.511177068850323_dp, 0.022218565326719057_dp, dp), &
         cmplx(0.7111770688503231_dp, 0.022218565326719057_dp, dp)], 8.0e-14_dp)

    c15 = logsumexp_keepdims(x15, [1, 13], cb15)
    call check_complex_array('rank15 real x complex b tuple', reshape(c15, [size(c15)]), &
        [cmplx(2.0315344048737587_dp, 0.12435499454676144_dp, dp), &
         cmplx(1.882277734798101_dp, 0.022218565326719057_dp, dp)], 8.0e-14_dp)

    cvalue = logsumexp(z15, b15)
    call check_complex_close('rank15 complex x real b', cvalue, &
        cmplx(2.341540821799777_dp, 0.0791088543187568_dp, dp), 8.0e-14_dp)
    call logsumexp_signed(z15, logabs, phase, b15)
    call check_close('rank15 complex x real b logabs', logabs, 2.341540821799777_dp, 8.0e-14_dp)
    call check_complex_close('rank15 complex x real b phase', phase, &
        cmplx(0.996872526127253_dp, 0.07902636681936677_dp, dp), 2.0e-14_dp)

    cvalue = logsumexp(z15, cb15)
    call check_complex_close('rank15 complex x complex b', cvalue, &
        cmplx(2.3750392636185547_dp, 0.1593611769741965_dp, dp), 8.0e-14_dp)
    call logsumexp_signed(z15, logabs, phase, cb15)
    call check_close('rank15 complex x complex b logabs', logabs, 2.3750392636185547_dp, 8.0e-14_dp)
    call check_complex_close('rank15 complex x complex b phase', phase, &
        cmplx(0.9873288580663374_dp, 0.15868751062835_dp, dp), 2.0e-14_dp)

    c14 = logsumexp(z15, 13, b15)
    call check_complex_array('rank15 complex x real b dim', reshape(c14, [size(c14)]), &
        [cmplx(0.9862943611198906_dp, 0.2_dp, dp), &
         cmplx(1.5862943611198906_dp, -0.7_dp, dp), &
         cmplx(1.5109302162163288_dp, 1.1_dp, dp), &
         cmplx(0.7109302162163288_dp, -0.3_dp, dp)], 8.0e-14_dp)

    c15 = logsumexp_keepdims(z15, [1, 13], cb15)
    call check_complex_array('rank15 complex x complex b tuple', reshape(c15, [size(c15)]), &
        [cmplx(1.9364748542772017_dp, -0.265447068770062_dp, dp), &
         cmplx(1.662944508326411_dp, 0.7319456352966522_dp, dp)], 9.0e-14_dp)

    call check_true('rank15 incompatible broadcast rejected', ieee_is_nan(logsumexp(x15, bad_b)))
    r15 = logsumexp_keepdims(x15, [1], bad_b)
    call check_true('rank15 incompatible keepdims is NaN', all(ieee_is_nan(r15)))

    allocate(empty15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 2))
    call logsumexp_signed_keepdims(empty15, [13], r15, s15, 1.0_dp)
    call check_shape('rank15 weighted zero-extent shape', shape(r15), &
        [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2])
    call check_true('rank15 weighted zero-extent logabs', all(r15 < -huge(1.0_dp)))
    call check_true('rank15 weighted zero-extent sign', all(s15 == 0.0_dp))

    print *, 'rank-6 through rank-15 weighted special reduction tests passed'

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

    subroutine check_real_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (size(actual) /= size(expected) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_real_array

    subroutine check_complex_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:) !! computed complex values
        complex(dp), intent(in) :: expected(:) !! reference complex values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_complex_array

    subroutine check_shape(name, actual, expected)
        character(len=*), intent(in) :: name !! check label printed on failure
        integer, intent(in) :: actual(:) !! computed extents
        integer, intent(in) :: expected(:) !! reference extents

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

    pure function complex_is_nan(value) result(is_nan)
        complex(dp), intent(in) :: value !! value inspected for a NaN component
        logical :: is_nan

        is_nan = ieee_is_nan(real(value, dp)) .or. ieee_is_nan(aimag(value))
    end function complex_is_nan

end program test_special_reductions_rank15_weighted
