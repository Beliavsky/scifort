! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank4_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    complex(dp), parameter :: z(2, 2, 2, 3) = reshape([ &
        cmplx(1.0_dp, 0.1_dp, dp), cmplx(2.0_dp, -0.2_dp, dp), &
        cmplx(3.0_dp, 0.3_dp, dp), cmplx(4.0_dp, -0.4_dp, dp), &
        cmplx(5.0_dp, 0.5_dp, dp), cmplx(6.0_dp, -0.6_dp, dp), &
        cmplx(7.0_dp, 0.7_dp, dp), cmplx(8.0_dp, -0.8_dp, dp), &
        cmplx(9.0_dp, 0.9_dp, dp), cmplx(10.0_dp, -1.0_dp, dp), &
        cmplx(11.0_dp, 1.1_dp, dp), cmplx(12.0_dp, -1.2_dp, dp), &
        cmplx(13.0_dp, 1.3_dp, dp), cmplx(14.0_dp, -1.4_dp, dp), &
        cmplx(15.0_dp, 1.5_dp, dp), cmplx(16.0_dp, -1.6_dp, dp), &
        cmplx(17.0_dp, 1.7_dp, dp), cmplx(18.0_dp, -1.8_dp, dp), &
        cmplx(19.0_dp, 1.9_dp, dp), cmplx(20.0_dp, -2.0_dp, dp), &
        cmplx(21.0_dp, 2.1_dp, dp), cmplx(22.0_dp, -2.2_dp, dp), &
        cmplx(23.0_dp, 2.3_dp, dp), cmplx(24.0_dp, -2.4_dp, dp)], [2, 2, 2, 3])

    complex(dp), allocatable :: cube(:, :, :)
    complex(dp), allocatable :: hyper(:, :, :, :)
    complex(dp), allocatable :: phases(:, :, :, :)
    complex(dp) :: phase
    complex(dp) :: value
    real(dp), allocatable :: logabs_hyper(:, :, :, :)
    real(dp) :: logabs
    integer :: i
    integer :: j
    integer :: l

    value = logsumexp(z)
    call check_complex_close('complex rank4 all', value, &
        cmplx(24.180608812664094_dp, -2.7296846851578334_dp, dp), 3.0e-14_dp)

    call logsumexp_signed(z, logabs, phase)
    call check_close('complex rank4 signed logabs', logabs, 24.180608812664094_dp, 3.0e-14_dp)
    call check_complex_close('complex rank4 signed phase', phase, &
        cmplx(-0.9163586199466364_dp, -0.4003584389637564_dp, dp), 3.0e-15_dp)

    cube = logsumexp(z, 4)
    call check_shape('complex rank4 dim4 shape', shape(cube), [2, 2, 2])
    call check_complex_close('complex rank4 dim4 first', cube(1, 1, 1), &
        cmplx(17.000233717439162_dp, 1.6997592976140101_dp, dp), 3.0e-14_dp)
    call check_complex_close('complex rank4 dim4 last', cube(2, 2, 2), &
        cmplx(24.000233717439162_dp, -2.3997592976140103_dp, dp), 3.0e-14_dp)

    hyper = logsumexp_keepdims(z, [1, 3])
    call check_shape('complex rank4 pair keepdims shape', shape(hyper), [1, 2, 1, 3])
    call check_complex_close('complex rank4 pair keepdims first', hyper(1, 1, 1, 1), &
        cmplx(6.212179081325873_dp, -0.32197518416586307_dp, dp), 3.0e-14_dp)
    call check_complex_close('complex rank4 pair keepdims last', hyper(1, 2, 1, 3), &
        cmplx(24.0715847835967_dp, -2.7485613021267463_dp, dp), 3.0e-14_dp)

    hyper = logsumexp_keepdims(z, [2, 3, 4])
    call check_complex_close('complex rank4 triple keepdims first', hyper(1, 1, 1, 1), &
        cmplx(23.14181815842384_dp, 2.2690113812453077_dp, dp), 3.0e-14_dp)
    call check_complex_close('complex rank4 triple keepdims second', hyper(2, 1, 1, 1), &
        cmplx(24.14181815842384_dp, -2.3690113812453077_dp, dp), 3.0e-14_dp)

    call logsumexp_signed_keepdims(z, [1, 3], logabs_hyper, phases)
    call check_close('complex rank4 signed keepdims logabs', logabs_hyper(1, 1, 1, 1), &
        6.212179081325873_dp, 3.0e-14_dp)
    call check_complex_close('complex rank4 signed keepdims phase', phases(1, 1, 1, 1), &
        exp(cmplx(0.0_dp, -0.32197518416586307_dp, dp)), 3.0e-15_dp)

    hyper = softmax(z, [1, 3])
    do l = 1, size(z, 4)
        do j = 1, size(z, 2)
            call check_complex_close('complex rank4 softmax normalization', &
                sum(hyper(:, j, :, l)), cmplx(1.0_dp, 0.0_dp, dp), 3.0e-15_dp)
        end do
    end do
    call check_complex_close('complex rank4 softmax reference', hyper(1, 1, 1, 1), &
        cmplx(0.004971739634308217_dp, 0.0022320313369974904_dp, dp), 3.0e-16_dp)
    call check_complex_close('complex rank4 softmax last reference', hyper(2, 2, 2, 3), &
        cmplx(0.8749366952473199_dp, 0.31795107738543293_dp, dp), 3.0e-15_dp)

    hyper = log_softmax(z, [2, 3, 4])
    do i = 1, size(z, 1)
        call check_complex_close('complex rank4 log_softmax normalization', &
            sum(exp(hyper(i, :, :, :))), cmplx(1.0_dp, 0.0_dp, dp), 3.0e-15_dp)
    end do
    call check_complex_close('complex rank4 log_softmax first', hyper(1, 1, 1, 1), &
        cmplx(-22.14181815842384_dp, -2.1690113812453076_dp, dp), 3.0e-14_dp)
    call check_complex_close('complex rank4 log_softmax last', hyper(2, 2, 2, 3), &
        cmplx(-0.14181815842383977_dp, -0.030988618754692786_dp, dp), 3.0e-14_dp)

    hyper = logsumexp_keepdims(z, [integer ::])
    call check_complex_close('complex rank4 empty-axis principal branch', hyper(1, 1, 1, 1), &
        z(1, 1, 1, 1), 2.0e-15_dp)
    hyper = softmax(z, [integer ::])
    call check_complex_close('complex rank4 empty softmax', hyper(2, 2, 2, 3), &
        cmplx(1.0_dp, 0.0_dp, dp), 0.0_dp)
    hyper = log_softmax(z, [integer ::])
    call check_complex_close('complex rank4 empty log_softmax', hyper(2, 2, 2, 3), &
        cmplx(0.0_dp, 0.0_dp, dp), 0.0_dp)

    print *, 'rank-4 complex special reduction tests passed'

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

    pure function complex_is_nan(value) result(has_nan)
        complex(dp), intent(in) :: value !! complex value inspected componentwise
        logical :: has_nan

        has_nan = ieee_is_nan(real(value, dp)) .or. ieee_is_nan(aimag(value))
    end function complex_is_nan

end program test_special_reductions_rank4_complex
