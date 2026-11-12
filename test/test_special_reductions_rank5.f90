! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank5
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    real(dp), parameter :: x(2, 2, 2, 2, 3) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, &
        9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp, 15.0_dp, 16.0_dp, &
        17.0_dp, 18.0_dp, 19.0_dp, 20.0_dp, 21.0_dp, 22.0_dp, 23.0_dp, 24.0_dp, &
        25.0_dp, 26.0_dp, 27.0_dp, 28.0_dp, 29.0_dp, 30.0_dp, 31.0_dp, 32.0_dp, &
        33.0_dp, 34.0_dp, 35.0_dp, 36.0_dp, 37.0_dp, 38.0_dp, 39.0_dp, 40.0_dp, &
        41.0_dp, 42.0_dp, 43.0_dp, 44.0_dp, 45.0_dp, 46.0_dp, 47.0_dp, 48.0_dp], &
        [2, 2, 2, 2, 3])
    complex(dp), parameter :: z(2, 2, 2, 2, 3) = reshape([ &
        cmplx(1.0_dp, -0.1_dp, dp), cmplx(2.0_dp, 0.2_dp, dp), &
        cmplx(3.0_dp, -0.3_dp, dp), cmplx(4.0_dp, 0.4_dp, dp), &
        cmplx(5.0_dp, -0.5_dp, dp), cmplx(6.0_dp, 0.6_dp, dp), &
        cmplx(7.0_dp, -0.7_dp, dp), cmplx(8.0_dp, 0.8_dp, dp), &
        cmplx(9.0_dp, -0.9_dp, dp), cmplx(10.0_dp, 1.0_dp, dp), &
        cmplx(11.0_dp, -1.1_dp, dp), cmplx(12.0_dp, 1.2_dp, dp), &
        cmplx(13.0_dp, -1.3_dp, dp), cmplx(14.0_dp, 1.4_dp, dp), &
        cmplx(15.0_dp, -1.5_dp, dp), cmplx(16.0_dp, 1.6_dp, dp), &
        cmplx(17.0_dp, -1.7_dp, dp), cmplx(18.0_dp, 1.8_dp, dp), &
        cmplx(19.0_dp, -1.9_dp, dp), cmplx(20.0_dp, 2.0_dp, dp), &
        cmplx(21.0_dp, -2.1_dp, dp), cmplx(22.0_dp, 2.2_dp, dp), &
        cmplx(23.0_dp, -2.3_dp, dp), cmplx(24.0_dp, 2.4_dp, dp), &
        cmplx(25.0_dp, -2.5_dp, dp), cmplx(26.0_dp, 2.6_dp, dp), &
        cmplx(27.0_dp, -2.7_dp, dp), cmplx(28.0_dp, 2.8_dp, dp), &
        cmplx(29.0_dp, -2.9_dp, dp), cmplx(30.0_dp, 3.0_dp, dp), &
        cmplx(31.0_dp, -3.1_dp, dp), cmplx(32.0_dp, 3.2_dp, dp), &
        cmplx(33.0_dp, -3.3_dp, dp), cmplx(34.0_dp, 3.4_dp, dp), &
        cmplx(35.0_dp, -3.5_dp, dp), cmplx(36.0_dp, 3.6_dp, dp), &
        cmplx(37.0_dp, -3.7_dp, dp), cmplx(38.0_dp, 3.8_dp, dp), &
        cmplx(39.0_dp, -3.9_dp, dp), cmplx(40.0_dp, 4.0_dp, dp), &
        cmplx(41.0_dp, -4.1_dp, dp), cmplx(42.0_dp, 4.2_dp, dp), &
        cmplx(43.0_dp, -4.3_dp, dp), cmplx(44.0_dp, 4.4_dp, dp), &
        cmplx(45.0_dp, -4.5_dp, dp), cmplx(46.0_dp, 4.6_dp, dp), &
        cmplx(47.0_dp, -4.7_dp, dp), cmplx(48.0_dp, 4.8_dp, dp)], &
        [2, 2, 2, 2, 3])

    real(dp), allocatable :: four(:, :, :, :)
    real(dp), allocatable :: five(:, :, :, :, :)
    real(dp), allocatable :: empty5(:, :, :, :, :)
    real(dp), allocatable :: signs(:, :, :, :, :)
    real(dp), allocatable :: complex_logabs(:, :, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    complex(dp), allocatable :: cfour(:, :, :, :)
    complex(dp), allocatable :: cfive(:, :, :, :, :)
    complex(dp), allocatable :: phases(:, :, :, :, :)
    complex(dp) :: phase
    complex(dp) :: cvalue
    real(dp) :: edge(1, 1, 1, 1, 2)
    integer :: i1
    integer :: i2
    integer :: i4

    call check_close('rank5 all', logsumexp(x), 48.45867514538708_dp, 4.0e-14_dp)

    four = logsumexp(x, 1)
    call check_shape('rank5 dim1 shape', shape(four), [2, 2, 2, 3])
    call check_close('rank5 dim1 first', four(1, 1, 1, 1), &
        2.313261687518223_dp, 2.0e-14_dp)
    call check_close('rank5 dim1 last', four(2, 2, 2, 3), &
        48.31326168751822_dp, 4.0e-14_dp)

    four = logsumexp(x, 5)
    call check_shape('rank5 dim5 shape', shape(four), [2, 2, 2, 2])
    call check_close('rank5 dim5 first', four(1, 1, 1, 1), &
        33.00000011253518_dp, 4.0e-14_dp)
    call check_close('rank5 dim5 last', four(2, 2, 2, 2), &
        48.00000011253518_dp, 4.0e-14_dp)

    five = logsumexp_keepdims(x, [1, 3, 5])
    call check_shape('rank5 tuple shape', shape(five), [1, 2, 1, 2, 1])
    call check_close('rank5 tuple first', five(1, 1, 1, 1, 1), &
        38.331411727971215_dp, 4.0e-14_dp)
    call check_close('rank5 tuple last', five(1, 2, 1, 2, 1), &
        48.331411727971215_dp, 4.0e-14_dp)
    call check_five('rank5 reversed tuple', logsumexp_keepdims(x, [5, 3, 1]), &
        five, 2.0e-15_dp)

    five = logsumexp_keepdims(x, [2, 3, 4, 5])
    call check_shape('rank5 four-axis shape', shape(five), [2, 1, 1, 1, 1])
    call check_close('rank5 four-axis first', five(1, 1, 1, 1, 1), &
        47.14541345786886_dp, 4.0e-14_dp)
    call check_close('rank5 four-axis second', five(2, 1, 1, 1, 1), &
        48.14541345786886_dp, 4.0e-14_dp)

    five = logsumexp_keepdims(x, [integer ::])
    call check_five('rank5 empty axes', five, x, 0.0_dp)

    allocate(empty5(2, 0, 2, 1, 1))
    five = logsumexp_keepdims(empty5, [2])
    call check_shape('rank5 zero-extent keepdims shape', shape(five), [2, 1, 2, 1, 1])
    call check_true('rank5 zero-extent logsumexp', all(five < -huge(1.0_dp)))

    call logsumexp_signed(x, logabs, sign_value)
    call check_close('rank5 signed all logabs', logabs, 48.45867514538708_dp, 4.0e-14_dp)
    call check_close('rank5 signed all sign', sign_value, 1.0_dp, 0.0_dp)

    edge(1, 1, 1, 1, 1) = ieee_value(0.0_dp, ieee_negative_inf)
    edge(1, 1, 1, 1, 2) = 1.0_dp
    call logsumexp_signed_keepdims(edge, [integer ::], five, signs)
    call check_close('rank5 empty-axis minus-inf sign', signs(1, 1, 1, 1, 1), 0.0_dp, 0.0_dp)
    call check_close('rank5 empty-axis finite sign', signs(1, 1, 1, 1, 2), 1.0_dp, 0.0_dp)

    five = softmax(x)
    call check_close('rank5 all softmax normalization', sum(five), 1.0_dp, 2.0e-15_dp)
    call check_close('rank5 all softmax first', five(1, 1, 1, 1, 1), &
        2.448833545946252e-21_dp, 2.0e-35_dp)
    call check_close('rank5 all softmax last', five(2, 2, 2, 2, 3), &
        0.6321205588285577_dp, 3.0e-15_dp)

    five = log_softmax(x)
    call check_close('rank5 all log_softmax normalization', sum(exp(five)), &
        1.0_dp, 4.0e-15_dp)
    call check_close('rank5 all log_softmax first', five(1, 1, 1, 1, 1), &
        -47.45867514538708_dp, 5.0e-14_dp)
    call check_close('rank5 all log_softmax last', five(2, 2, 2, 2, 3), &
        -0.45867514538708193_dp, 5.0e-14_dp)

    five = softmax(x, [1, 3, 5])
    do i4 = 1, size(x, 4)
        do i2 = 1, size(x, 2)
            call check_close('rank5 tuple softmax normalization', &
                sum(five(:, i2, :, i4, :)), 1.0_dp, 2.0e-15_dp)
        end do
    end do
    call check_close('rank5 tuple softmax first', five(1, 1, 1, 1, 1), &
        6.125956165667844e-17_dp, 3.0e-31_dp)
    call check_close('rank5 tuple softmax last', five(2, 2, 2, 2, 3), &
        0.7179095247500945_dp, 3.0e-15_dp)

    five = log_softmax(x, [2, 3, 4, 5])
    do i1 = 1, size(x, 1)
        call check_close('rank5 tuple log_softmax normalization', &
            sum(exp(five(i1, :, :, :, :))), 1.0_dp, 5.0e-15_dp)
    end do
    call check_close('rank5 tuple log_softmax first', five(1, 1, 1, 1, 1), &
        -46.14541345786886_dp, 4.0e-14_dp)
    call check_close('rank5 tuple log_softmax last', five(2, 2, 2, 2, 3), &
        -0.14541345786885912_dp, 4.0e-14_dp)

    five = softmax(x, [integer ::])
    call check_close('rank5 empty softmax', maxval(abs(five - 1.0_dp)), 0.0_dp, 0.0_dp)
    five = log_softmax(x, [integer ::])
    call check_close('rank5 empty log_softmax', maxval(abs(five)), 0.0_dp, 0.0_dp)

    five = logsumexp_keepdims(x, [1, 1])
    call check_shape('rank5 duplicate keepdims rejected', shape(five), [0, 0, 0, 0, 0])
    five = softmax(x, [1, 1])
    call check_true('rank5 duplicate softmax rejected', all(ieee_is_nan(five)))

    cvalue = logsumexp(z)
    call check_complex_close('complex rank5 all', cvalue, &
        cmplx(47.68322375978535_dp, -1.5064665313240377_dp, dp), 5.0e-14_dp)
    call logsumexp_signed(z, logabs, phase)
    call check_close('complex rank5 signed logabs', logabs, 47.68322375978535_dp, 5.0e-14_dp)
    call check_complex_close('complex rank5 signed phase', phase, &
        cmplx(0.06428543507648364_dp, -0.9979315521803223_dp, dp), 4.0e-15_dp)

    cfour = logsumexp(z, 5)
    call check_shape('complex rank5 dim5 shape', shape(cfour), [2, 2, 2, 2])
    call check_complex_close('complex rank5 dim5 first', cfour(1, 1, 1, 1), &
        cmplx(32.99999999671402_dp, 2.9831854196667766_dp, dp), 5.0e-14_dp)
    call check_complex_close('complex rank5 dim5 last', cfour(2, 2, 2, 2), &
        cmplx(47.99999999671402_dp, -1.4831854196667758_dp, dp), 5.0e-14_dp)

    cfive = logsumexp_keepdims(z, [1, 3, 5])
    call check_shape('complex rank5 tuple shape', shape(cfive), [1, 2, 1, 2, 1])
    call check_complex_close('complex rank5 tuple first', cfive(1, 1, 1, 1, 1), &
        cmplx(38.18500101699504_dp, -2.7845248714768083_dp, dp), 5.0e-14_dp)
    call check_complex_close('complex rank5 tuple last', cfive(1, 2, 1, 2, 1), &
        cmplx(47.55979886931292_dp, -1.4546619716751952_dp, dp), 5.0e-14_dp)

    call logsumexp_signed_keepdims(z, [1, 3, 5], complex_logabs, phases)
    call check_complex_close('complex rank5 tuple phase', phases(1, 1, 1, 1, 1), &
        cmplx(-0.9369257436341061_dp, -0.3495284693922329_dp, dp), 4.0e-15_dp)

    cfive = softmax(z, [1, 3, 5])
    do i4 = 1, size(z, 4)
        do i2 = 1, size(z, 2)
            call check_complex_close('complex rank5 softmax normalization', &
                sum(cfive(:, i2, :, i4, :)), cmplx(1.0_dp, 0.0_dp, dp), 5.0e-15_dp)
        end do
    end do
    call check_complex_close('complex rank5 softmax first', cfive(1, 1, 1, 1, 1), &
        cmplx(-6.363870120969607e-17_dp, 3.129764489586486e-17_dp, dp), 4.0e-31_dp)
    call check_complex_close('complex rank5 softmax last', cfive(2, 2, 2, 2, 3), &
        cmplx(1.552387836479639_dp, -0.044291291257165084_dp, dp), 5.0e-15_dp)

    cfive = log_softmax(z, [2, 3, 4, 5])
    do i1 = 1, size(z, 1)
        call check_complex_close('complex rank5 log_softmax normalization', &
            sum(exp(cfive(i1, :, :, :, :))), cmplx(1.0_dp, 0.0_dp, dp), 5.0e-15_dp)
    end do
    call check_complex_close('complex rank5 log_softmax first', cfive(1, 1, 1, 1, 1), &
        cmplx(-46.141818158396_dp, 4.569011381219808_dp, dp), 5.0e-14_dp)
    call check_complex_close('complex rank5 log_softmax last', cfive(2, 2, 2, 2, 3), &
        cmplx(-0.14181815839600226_dp, 0.030988618780192413_dp, dp), 5.0e-14_dp)

    print *, 'rank-5 real/complex special reduction tests passed'

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

    subroutine check_five(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :, :, :, :) !! computed rank-5 array
        real(dp), intent(in) :: expected(:, :, :, :, :) !! reference rank-5 array
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_five

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

    pure function complex_is_nan(value) result(has_nan)
        complex(dp), intent(in) :: value !! complex value inspected componentwise
        logical :: has_nan

        has_nan = ieee_is_nan(real(value, dp)) .or. ieee_is_nan(aimag(value))
    end function complex_is_nan

end program test_special_reductions_rank5
