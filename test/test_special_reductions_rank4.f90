! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank4
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, ieee_value
    use scifort_special, only : dp, log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    real(dp), parameter :: x(2, 2, 2, 3) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, &
        10.0_dp, 11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp, 15.0_dp, 16.0_dp, 17.0_dp, &
        18.0_dp, 19.0_dp, 20.0_dp, 21.0_dp, 22.0_dp, 23.0_dp, 24.0_dp], [2, 2, 2, 3])
    real(dp), allocatable :: cube(:, :, :)
    real(dp) :: edge(1, 1, 1, 2)
    real(dp), allocatable :: hyper(:, :, :, :)
    real(dp), allocatable :: signs(:, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    integer :: i
    integer :: j
    integer :: l

    call check_close('rank4 all', logsumexp(x), 24.458675145349332_dp, 2.0e-14_dp)

    cube = logsumexp(x, 1)
    call check_shape('rank4 dim1 shape', shape(cube), [2, 2, 3])
    call check_close('rank4 dim1 first', cube(1, 1, 1), 2.313261687518223_dp, 2.0e-14_dp)
    call check_close('rank4 dim1 last', cube(2, 2, 3), 24.313261687518224_dp, 2.0e-14_dp)

    cube = logsumexp(x, 4)
    call check_shape('rank4 dim4 shape', shape(cube), [2, 2, 2])
    call check_close('rank4 dim4 first', cube(1, 1, 1), 17.000335518870326_dp, 2.0e-14_dp)
    call check_close('rank4 dim4 last', cube(2, 2, 2), 24.000335518870326_dp, 2.0e-14_dp)

    hyper = logsumexp_keepdims(x, [1, 3])
    call check_shape('rank4 pair keepdims shape', shape(hyper), [1, 2, 1, 3])
    call check_close('rank4 pair keepdims first', hyper(1, 1, 1, 1), &
        6.331411615436033_dp, 2.0e-14_dp)
    call check_close('rank4 pair keepdims last', hyper(1, 2, 1, 3), &
        24.331411615436032_dp, 2.0e-14_dp)
    call check_hyper('rank4 reversed pair keepdims', logsumexp_keepdims(x, [3, 1]), &
        hyper, 1.0e-15_dp)

    hyper = logsumexp_keepdims(x, [2, 3, 4])
    call check_shape('rank4 triple keepdims shape', shape(hyper), [2, 1, 1, 1])
    call check_close('rank4 triple keepdims first', hyper(1, 1, 1, 1), &
        23.145413457831108_dp, 2.0e-14_dp)
    call check_close('rank4 triple keepdims second', hyper(2, 1, 1, 1), &
        24.145413457831108_dp, 2.0e-14_dp)

    hyper = logsumexp_keepdims(x, [4, 1])
    call check_shape('rank4 nonadjacent keepdims shape', shape(hyper), [1, 2, 2, 1])
    call check_close('rank4 nonadjacent keepdims first', hyper(1, 1, 1, 1), &
        18.313597206388547_dp, 2.0e-14_dp)
    call check_close('rank4 nonadjacent keepdims last', hyper(1, 2, 2, 1), &
        24.313597206388547_dp, 2.0e-14_dp)

    hyper = logsumexp_keepdims(x, [1, 2, 3, 4])
    call check_shape('rank4 all keepdims shape', shape(hyper), [1, 1, 1, 1])
    call check_close('rank4 all keepdims', hyper(1, 1, 1, 1), &
        24.458675145349332_dp, 2.0e-14_dp)

    hyper = logsumexp_keepdims(x, [integer ::])
    call check_hyper('rank4 empty keepdims', hyper, x, 0.0_dp)

    call logsumexp_signed(x, logabs, sign_value)
    call check_close('rank4 signed all logabs', logabs, 24.458675145349332_dp, 2.0e-14_dp)
    call check_close('rank4 signed all sign', sign_value, 1.0_dp, 0.0_dp)

    edge(1, 1, 1, 1) = ieee_value(0.0_dp, ieee_negative_inf)
    edge(1, 1, 1, 2) = 1.0_dp
    call logsumexp_signed_keepdims(edge, [integer ::], hyper, signs)
    call check_close('rank4 empty-axis minus-inf sign', signs(1, 1, 1, 1), 0.0_dp, 0.0_dp)
    call check_close('rank4 empty-axis finite sign', signs(1, 1, 1, 2), 1.0_dp, 0.0_dp)

    call logsumexp_signed_keepdims(x, [2, 3, 4], hyper, signs)
    call check_close('rank4 signed keepdims first logabs', hyper(1, 1, 1, 1), &
        23.145413457831108_dp, 2.0e-14_dp)
    call check_close('rank4 signed keepdims sign', signs(2, 1, 1, 1), 1.0_dp, 0.0_dp)

    hyper = softmax(x, [1, 3])
    do l = 1, size(x, 4)
        do j = 1, size(x, 2)
            call check_close('rank4 tuple softmax normalization', sum(hyper(:, j, :, l)), &
                1.0_dp, 8.0e-16_dp)
        end do
    end do
    call check_close('rank4 tuple softmax reference', hyper(1, 1, 1, 1), &
        0.00483723687226407_dp, 2.0e-16_dp)
    call check_close('rank4 tuple softmax last reference', hyper(2, 2, 2, 3), &
        0.7179096055401775_dp, 2.0e-15_dp)

    hyper = log_softmax(x, [2, 3, 4])
    do i = 1, size(x, 1)
        call check_close('rank4 tuple log_softmax normalization', sum(exp(hyper(i, :, :, :))), &
            1.0_dp, 2.0e-15_dp)
    end do
    call check_close('rank4 tuple log_softmax first', hyper(1, 1, 1, 1), &
        -22.145413457831108_dp, 2.0e-14_dp)
    call check_close('rank4 tuple log_softmax last', hyper(2, 2, 2, 3), &
        -0.1454134578311077_dp, 2.0e-14_dp)

    hyper = softmax(x, [integer ::])
    call check_close('rank4 empty softmax', maxval(abs(hyper - 1.0_dp)), 0.0_dp, 0.0_dp)
    hyper = log_softmax(x, [integer ::])
    call check_close('rank4 empty log_softmax', maxval(abs(hyper)), 0.0_dp, 0.0_dp)

    hyper = logsumexp_keepdims(x, [1, 1])
    call check_shape('rank4 duplicate keepdims rejected', shape(hyper), [0, 0, 0, 0])
    hyper = softmax(x, [1, 1])
    call check_true('rank4 duplicate softmax rejected', all(ieee_is_nan(hyper)))

    print *, 'rank-4 real special reduction tests passed'

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

    subroutine check_hyper(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :, :, :) !! computed rank-4 array
        real(dp), intent(in) :: expected(:, :, :, :) !! reference rank-4 array
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_hyper

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

end program test_special_reductions_rank4
