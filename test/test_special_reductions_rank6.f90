! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank6
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    real(dp) :: x(2, 2, 2, 2, 2, 3)
    complex(dp) :: z(2, 2, 2, 2, 2, 3)
    real(dp), allocatable :: five(:, :, :, :, :)
    real(dp), allocatable :: six(:, :, :, :, :, :)
    real(dp), allocatable :: empty6(:, :, :, :, :, :)
    real(dp), allocatable :: signs(:, :, :, :, :, :)
    real(dp), allocatable :: complex_logabs(:, :, :, :, :, :)
    complex(dp), allocatable :: cfive(:, :, :, :, :)
    complex(dp), allocatable :: csix(:, :, :, :, :, :)
    complex(dp), allocatable :: phases(:, :, :, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    complex(dp) :: phase
    complex(dp) :: cvalue
    real(dp) :: edge(1, 1, 1, 1, 1, 2)
    integer :: i
    integer :: i1
    integer :: i2
    integer :: i4
    integer :: i5

    x = reshape([(real(i, dp), i = 1, size(x))], shape(x))
    z = reshape([(cmplx(real(i, dp), 0.07_dp * real(i, dp) * &
        merge(-1.0_dp, 1.0_dp, mod(i, 2) == 1), kind=dp), i = 1, size(z))], shape(z))

    call check_close('rank6 all', logsumexp(x), 96.45867514538708_dp, 8.0e-14_dp)

    five = logsumexp(x, 1)
    call check_shape('rank6 dim1 shape', shape(five), [2, 2, 2, 2, 3])
    call check_close('rank6 dim1 first', five(1, 1, 1, 1, 1), &
        2.313261687518223_dp, 3.0e-14_dp)
    call check_close('rank6 dim1 last', five(2, 2, 2, 2, 3), &
        96.31326168751822_dp, 8.0e-14_dp)

    five = logsumexp(x, 6)
    call check_shape('rank6 dim6 shape', shape(five), [2, 2, 2, 2, 2])
    call check_close('rank6 dim6 first', five(1, 1, 1, 1, 1), &
        65.00000000000001_dp, 8.0e-14_dp)
    call check_close('rank6 dim6 last', five(2, 2, 2, 2, 2), &
        96.00000000000001_dp, 8.0e-14_dp)

    six = logsumexp_keepdims(x, [1, 3, 6])
    call check_shape('rank6 tuple shape', shape(six), [1, 2, 1, 2, 2, 1])
    call check_close('rank6 tuple first', six(1, 1, 1, 1, 1, 1), &
        70.33141161543604_dp, 8.0e-14_dp)
    call check_close('rank6 tuple last', six(1, 2, 1, 2, 2, 1), &
        96.33141161543604_dp, 8.0e-14_dp)
    call check_six('rank6 reversed tuple', logsumexp_keepdims(x, [6, 3, 1]), &
        six, 3.0e-15_dp)

    six = logsumexp_keepdims(x, [2, 3, 4, 5, 6])
    call check_shape('rank6 five-axis shape', shape(six), [2, 1, 1, 1, 1, 1])
    call check_close('rank6 five-axis first', six(1, 1, 1, 1, 1, 1), &
        95.14541345786886_dp, 8.0e-14_dp)
    call check_close('rank6 five-axis second', six(2, 1, 1, 1, 1, 1), &
        96.14541345786886_dp, 8.0e-14_dp)

    six = logsumexp_keepdims(x, [integer ::])
    call check_six('rank6 empty axes', six, x, 0.0_dp)

    allocate(empty6(2, 0, 2, 1, 1, 1))
    six = logsumexp_keepdims(empty6, [2])
    call check_shape('rank6 zero-extent keepdims shape', shape(six), [2, 1, 2, 1, 1, 1])
    call check_true('rank6 zero-extent logsumexp', all(six < -huge(1.0_dp)))

    call logsumexp_signed(x, logabs, sign_value)
    call check_close('rank6 signed all logabs', logabs, 96.45867514538708_dp, 8.0e-14_dp)
    call check_close('rank6 signed all sign', sign_value, 1.0_dp, 0.0_dp)

    edge(1, 1, 1, 1, 1, 1) = ieee_value(0.0_dp, ieee_negative_inf)
    edge(1, 1, 1, 1, 1, 2) = 1.0_dp
    call logsumexp_signed_keepdims(edge, [integer ::], six, signs)
    call check_close('rank6 empty-axis minus-inf sign', signs(1, 1, 1, 1, 1, 1), 0.0_dp, 0.0_dp)
    call check_close('rank6 empty-axis finite sign', signs(1, 1, 1, 1, 1, 2), 1.0_dp, 0.0_dp)

    six = softmax(x)
    call check_close('rank6 all softmax normalization', sum(six), 1.0_dp, 3.0e-15_dp)
    call check_close('rank6 all softmax last', six(2, 2, 2, 2, 2, 3), &
        0.6321205588285577_dp, 3.0e-15_dp)

    six = log_softmax(x)
    call check_close('rank6 all log_softmax normalization', sum(exp(six)), &
        1.0_dp, 5.0e-15_dp)
    call check_close('rank6 all log_softmax first', six(1, 1, 1, 1, 1, 1), &
        -95.45867514538708_dp, 8.0e-14_dp)
    call check_close('rank6 all log_softmax last', six(2, 2, 2, 2, 2, 3), &
        -0.45867514538708193_dp, 8.0e-14_dp)

    six = softmax(x, [1, 3, 6])
    do i5 = 1, size(x, 5)
        do i4 = 1, size(x, 4)
            do i2 = 1, size(x, 2)
                call check_close('rank6 tuple softmax normalization', &
                    sum(six(:, i2, :, i4, i5, :)), 1.0_dp, 3.0e-15_dp)
            end do
        end do
    end do
    call check_close('rank6 tuple softmax last', six(2, 2, 2, 2, 2, 3), &
        0.7179096055401682_dp, 4.0e-15_dp)

    six = log_softmax(x, [2, 3, 4, 5, 6])
    do i1 = 1, size(x, 1)
        call check_close('rank6 tuple log_softmax normalization', &
            sum(exp(six(i1, :, :, :, :, :))), 1.0_dp, 6.0e-15_dp)
    end do
    call check_close('rank6 tuple log_softmax first', six(1, 1, 1, 1, 1, 1), &
        -94.14541345786886_dp, 8.0e-14_dp)
    call check_close('rank6 tuple log_softmax last', six(2, 2, 2, 2, 2, 3), &
        -0.14541345786885912_dp, 8.0e-14_dp)

    six = softmax(x, [integer ::])
    call check_close('rank6 empty softmax', maxval(abs(six - 1.0_dp)), 0.0_dp, 0.0_dp)
    six = log_softmax(x, [integer ::])
    call check_close('rank6 empty log_softmax', maxval(abs(six)), 0.0_dp, 0.0_dp)

    six = logsumexp_keepdims(x, [1, 1])
    call check_shape('rank6 duplicate keepdims rejected', shape(six), [0, 0, 0, 0, 0, 0])
    six = softmax(x, [1, 1])
    call check_true('rank6 duplicate softmax rejected', all(ieee_is_nan(six)))

    cvalue = logsumexp(z)
    call check_complex_close('complex rank6 all', cvalue, &
        cmplx(96.3996467228793_dp, 0.2175267401320209_dp, dp), 1.0e-13_dp)
    call logsumexp_signed(z, logabs, phase)
    call check_close('complex rank6 signed logabs', logabs, 96.3996467228793_dp, 1.0e-13_dp)
    call check_complex_close('complex rank6 signed phase', phase, &
        cmplx(0.9764342025612696_dp, 0.21581531008836494_dp, dp), 8.0e-15_dp)

    cfive = logsumexp(z, 6)
    call check_shape('complex rank6 dim6 shape', shape(cfive), [2, 2, 2, 2, 2])
    call check_complex_close('complex rank6 dim6 first', cfive(1, 1, 1, 1, 1), &
        cmplx(65.0_dp, 1.7331853071795953_dp, dp), 1.0e-13_dp)
    call check_complex_close('complex rank6 dim6 last', cfive(2, 2, 2, 2, 2), &
        cmplx(96.0_dp, 0.43681469282040464_dp, dp), 1.0e-13_dp)

    csix = logsumexp_keepdims(z, [1, 3, 6])
    call check_shape('complex rank6 tuple shape', shape(csix), [1, 2, 1, 2, 2, 1])
    call check_complex_close('complex rank6 tuple first', csix(1, 1, 1, 1, 1, 1), &
        cmplx(69.59712295693967_dp, -1.2244422472967251_dp, dp), 1.0e-13_dp)
    call check_complex_close('complex rank6 tuple last', csix(1, 2, 1, 2, 2, 1), &
        cmplx(96.26823338286135_dp, 0.22629331837296274_dp, dp), 1.0e-13_dp)

    call logsumexp_signed_keepdims(z, [1, 3, 6], complex_logabs, phases)
    call check_complex_close('complex rank6 tuple phase', phases(1, 1, 1, 1, 1, 1), &
        cmplx(0.3394706577879423_dp, -0.9406166448139337_dp, dp), 8.0e-15_dp)

    csix = softmax(z, [1, 3, 6])
    do i5 = 1, size(z, 5)
        do i4 = 1, size(z, 4)
            do i2 = 1, size(z, 2)
                call check_complex_close('complex rank6 softmax normalization', &
                    sum(csix(:, i2, :, i4, i5, :)), cmplx(1.0_dp, 0.0_dp, dp), 8.0e-15_dp)
            end do
        end do
    end do
    call check_complex_close('complex rank6 softmax last', csix(2, 2, 2, 2, 2, 3), &
        cmplx(0.747845665915911_dp, 0.15980531960168323_dp, dp), 1.0e-14_dp)

    csix = log_softmax(z, [2, 3, 4, 5, 6])
    do i1 = 1, size(z, 1)
        call check_complex_close('complex rank6 log_softmax normalization', &
            sum(exp(csix(i1, :, :, :, :, :))), cmplx(1.0_dp, 0.0_dp, dp), 1.0e-14_dp)
    end do
    call check_complex_close('complex rank6 log_softmax first', csix(1, 1, 1, 1, 1, 1), &
        cmplx(-94.14364553161091_dp, 6.558195892086384_dp, dp), 1.0e-13_dp)
    call check_complex_close('complex rank6 log_softmax last', csix(2, 2, 2, 2, 2, 3), &
        cmplx(-0.14364553161092009_dp, 0.021804107913616134_dp, dp), 1.0e-13_dp)

    print *, 'rank-6 real/complex special reduction tests passed'

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

    subroutine check_six(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :, :, :, :, :) !! computed rank-6 array
        real(dp), intent(in) :: expected(:, :, :, :, :, :) !! reference rank-6 array
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_six

    subroutine check_shape(name, actual, expected)
        character(len=*), intent(in) :: name !! check label printed on failure
        integer, intent(in) :: actual(:) !! computed extents
        integer, intent(in) :: expected(:) !! expected extents

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', name, ': rank mismatch'
            error stop 1
        end if
        if (any(actual /= expected)) then
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

end program test_special_reductions_rank6
