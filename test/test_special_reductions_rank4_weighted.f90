! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank4_weighted
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims
    implicit none

    real(dp), parameter :: x(2, 1, 3, 1) = reshape([ &
        0.2_dp, -0.3_dp, 0.7_dp, -1.1_dp, 1.2_dp, 0.4_dp], [2, 1, 3, 1])
    real(dp), parameter :: zimag(2, 1, 3, 1) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp], [2, 1, 3, 1])
    real(dp) :: b1(5)
    real(dp) :: b2(3, 5)
    real(dp) :: b3(4, 3, 5)
    real(dp) :: b4(2, 4, 3, 5)
    real(dp) :: bneg(4, 3, 5)
    real(dp) :: bad(2, 5)
    real(dp) :: bzero(0, 3, 5)
    real(dp) :: xfull(2, 2, 2, 3)
    complex(dp) :: bc(4, 3, 5)
    complex(dp) :: z(2, 1, 3, 1)
    complex(dp) :: zfull(2, 2, 2, 3)
    real(dp), allocatable :: cube(:, :, :)
    real(dp), allocatable :: hyper(:, :, :, :)
    real(dp), allocatable :: base_hyper(:, :, :, :)
    real(dp), allocatable :: signs(:, :, :, :)
    real(dp), allocatable :: base_signs(:, :, :, :)
    complex(dp), allocatable :: ccube(:, :, :)
    complex(dp), allocatable :: chyper(:, :, :, :)
    complex(dp), allocatable :: base_chyper(:, :, :, :)
    complex(dp), allocatable :: phases(:, :, :, :)
    complex(dp), allocatable :: base_phases(:, :, :, :)
    real(dp), allocatable :: logabs_hyper(:, :, :, :)
    real(dp), allocatable :: base_logabs(:, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    complex(dp) :: phase
    complex(dp) :: cvalue
    integer :: d
    integer :: dims(4)
    integer :: i
    integer :: mask
    integer :: nd

    b1 = [(0.4_dp + 0.2_dp * real(i, dp), i = 1, size(b1))]
    b2 = reshape([(0.5_dp + 0.07_dp * real(i - 1, dp), i = 1, size(b2))], shape(b2))
    b3 = reshape([(0.4_dp + 0.03_dp * real(i, dp), i = 1, size(b3))], shape(b3))
    b4 = reshape([(0.3_dp + 0.01_dp * real(i, dp), i = 1, size(b4))], shape(b4))
    bneg = b3
    bneg(2, :, :) = -0.7_dp * b3(2, :, :)
    bneg(4, :, :) = -0.7_dp * b3(4, :, :)
    bc = reshape([(cmplx(0.4_dp + 0.03_dp * real(i, dp), &
        0.2_dp * (0.4_dp + 0.03_dp * real(i, dp)) * sin(real(i, dp)), kind=dp), &
        i = 1, size(bc))], shape(bc))
    z = cmplx(x, zimag, kind=dp)
    xfull = reshape([(real(i, dp), i = 1, size(xfull))], shape(xfull))
    zfull = cmplx(xfull, 0.02_dp * xfull, kind=dp)
    bad = 1.0_dp

    call check_close('rank4 scalar weight', logsumexp(x, 1.7_dp), &
        2.7411842784473315_dp, 3.0e-14_dp)
    call check_close('rank4 rank1 weight', logsumexp(x, b1), &
        3.8199939398192617_dp, 3.0e-14_dp)
    call check_close('rank4 rank2 weight', logsumexp(x, b2), &
        3.831795002489291_dp, 3.0e-14_dp)
    call check_close('rank4 rank3 expanding weight', logsumexp(x, b=b3), &
        5.508237924578352_dp, 3.0e-14_dp)
    call check_close('rank4 rank4 expanding weight', logsumexp(x, b4), &
        5.131358068801565_dp, 3.0e-14_dp)

    cube = logsumexp(x, dim=2, b=b3)
    call check_shape3('rank4 weighted dim2 shape', shape(cube), [2, 3, 5])
    call check_close('rank4 weighted dim2 first', cube(1, 1, 1), &
        0.8418538861723948_dp, 3.0e-14_dp)
    call check_close('rank4 weighted dim2 last', cube(2, 3, 5), &
        2.5540850846756014_dp, 3.0e-14_dp)

    hyper = logsumexp_keepdims(x, dims=[2, 4], b=b3)
    call check_shape4('rank4 weighted keepdims shape', shape(hyper), [2, 1, 3, 1])
    call check_close('rank4 weighted keepdims first', hyper(1, 1, 1, 1), &
        3.3738784589374653_dp, 3.0e-14_dp)
    call check_close('rank4 weighted keepdims last', hyper(2, 1, 3, 1), &
        3.7568971227655754_dp, 3.0e-14_dp)

    call logsumexp_signed(x, logabs, sign_value, bneg)
    call check_close('rank4 weighted signed logabs', logabs, 3.546209321881256_dp, 3.0e-14_dp)
    call check_close('rank4 weighted signed sign', sign_value, 1.0_dp, 0.0_dp)
    call logsumexp_signed_keepdims(x, [2, 4], hyper, signs, bneg)
    call check_close('rank4 weighted signed keepdims first', hyper(1, 1, 1, 1), &
        1.402972303992352_dp, 3.0e-14_dp)
    call check_close('rank4 weighted signed keepdims last', hyper(2, 1, 3, 1), &
        1.7987168811184482_dp, 3.0e-14_dp)
    call check_true('rank4 weighted signed keepdims signs', &
        all(abs(signs - 1.0_dp) <= epsilon(1.0_dp)))

    cvalue = logsumexp(x, bc)
    call check_complex_close('rank4 real input complex weights', cvalue, &
        cmplx(5.508271445921257_dp, 0.008187914062446185_dp, dp), 4.0e-14_dp)
    ccube = logsumexp(x, 4, bc)
    call check_complex_close('rank4 real input complex weights dim4 first', ccube(1, 1, 1), &
        cmplx(1.9521554682019777_dp, -0.07684663502030213_dp, dp), 4.0e-14_dp)
    call check_complex_close('rank4 real input complex weights dim4 last', ccube(2, 4, 3), &
        cmplx(2.4104355309214713_dp, -0.13363254701122768_dp, dp), 4.0e-14_dp)
    chyper = logsumexp_keepdims(x, [1, 3], bc)
    call check_shape4('rank4 complex weight keepdims shape', shape(chyper), [1, 4, 1, 5])
    call check_complex_close('rank4 complex weight keepdims first', chyper(1, 1, 1, 1), &
        cmplx(1.6791119912393941_dp, 0.02988945725172095_dp, dp), 4.0e-14_dp)
    call check_complex_close('rank4 complex weight keepdims last', chyper(1, 4, 1, 5), &
        cmplx(2.960998645261328_dp, -0.020480153179275574_dp, dp), 4.0e-14_dp)
    call logsumexp_signed(x, logabs, phase, bc)
    call check_close('rank4 complex weight signed logabs', logabs, &
        5.508271445921257_dp, 4.0e-14_dp)
    call check_complex_close('rank4 complex weight signed phase', phase, &
        cmplx(0.9999664792189286_dp, 0.008187822573816811_dp, dp), 4.0e-15_dp)

    cvalue = logsumexp(z, b3)
    call check_complex_close('rank4 complex input real weights', cvalue, &
        cmplx(5.42255308083794_dp, 0.14689861736024773_dp, dp), 4.0e-14_dp)
    ccube = logsumexp(z, 3, b3)
    call check_complex_close('rank4 complex real-weight dim3 first', ccube(1, 1, 1), &
        cmplx(1.3395674492858647_dp, 0.3885514157394176_dp, dp), 4.0e-14_dp)
    call check_complex_close('rank4 complex real-weight dim3 last', ccube(2, 4, 5), &
        cmplx(1.6760502853994144_dp, -0.46808016299668237_dp, dp), 4.0e-14_dp)
    chyper = logsumexp_keepdims(z, [2, 4], b3)
    call check_complex_close('rank4 complex real-weight keepdims first', chyper(1, 1, 1, 1), &
        cmplx(3.3738784589374653_dp, 0.1_dp, dp), 4.0e-14_dp)
    call check_complex_close('rank4 complex real-weight keepdims last', chyper(2, 1, 3, 1), &
        cmplx(3.7568971227655754_dp, -0.6_dp, dp), 4.0e-14_dp)
    call logsumexp_signed(z, logabs, phase, b3)
    call check_close('rank4 complex real-weight signed logabs', logabs, &
        5.42255308083794_dp, 4.0e-14_dp)
    call check_complex_close('rank4 complex real-weight signed phase', phase, &
        cmplx(0.9892297867426179_dp, 0.14637086124346943_dp, dp), 4.0e-15_dp)

    cvalue = logsumexp(z, bc)
    call check_complex_close('rank4 complex input complex weights', cvalue, &
        cmplx(5.425192736859498_dp, 0.154262922590703_dp, dp), 4.0e-14_dp)
    call logsumexp_signed(z, logabs, phase, bc)
    call check_close('rank4 complex-complex signed logabs', logabs, &
        5.425192736859498_dp, 4.0e-14_dp)
    call check_complex_close('rank4 complex-complex signed phase', phase, &
        cmplx(0.988125052462594_dp, 0.15365181644157683_dp, dp), 4.0e-15_dp)
    call logsumexp_signed_keepdims(z, [1, 3], logabs_hyper, phases, bc)
    call check_complex_close('rank4 complex-complex signed keepdims phase', phases(1, 1, 1, 1), &
        exp(cmplx(0.0_dp, 0.1798330473200907_dp, dp)), 5.0e-15_dp)


    do mask = 0, 15
        nd = 0
        do d = 1, 4
            if (btest(mask, d - 1)) then
                nd = nd + 1
                dims(nd) = d
            end if
        end do

        base_hyper = logsumexp_keepdims(xfull, dims(:nd))
        hyper = logsumexp_keepdims(xfull, dims(:nd), 2.0_dp)
        call check_true('rank4 scalar-weight all real keepdims modes shape', &
            all(shape(hyper) == shape(base_hyper)))
        call check_true('rank4 scalar-weight all real keepdims modes value', &
            maxval(abs(hyper - base_hyper - log(2.0_dp))) <= 8.0e-14_dp)
        call logsumexp_signed_keepdims(xfull, dims(:nd), base_logabs, base_signs)
        call logsumexp_signed_keepdims(xfull, dims(:nd), logabs_hyper, signs, 2.0_dp)
        call check_true('rank4 scalar-weight all real signed modes logabs', &
            maxval(abs(logabs_hyper - base_logabs - log(2.0_dp))) <= 8.0e-14_dp)
        call check_true('rank4 scalar-weight all real signed modes sign', &
            maxval(abs(signs - base_signs)) <= 2.0e-15_dp)

        base_chyper = logsumexp_keepdims(zfull, dims(:nd))
        chyper = logsumexp_keepdims(zfull, dims(:nd), cmplx(2.0_dp, 0.0_dp, kind=dp))
        call check_true('rank4 scalar-weight all complex keepdims modes shape', &
            all(shape(chyper) == shape(base_chyper)))
        call check_true('rank4 scalar-weight all complex keepdims modes value', &
            maxval(abs(chyper - base_chyper - log(2.0_dp))) <= 1.0e-13_dp)
        call logsumexp_signed_keepdims(zfull, dims(:nd), base_logabs, base_phases)
        call logsumexp_signed_keepdims(zfull, dims(:nd), logabs_hyper, phases, &
            cmplx(2.0_dp, 0.0_dp, kind=dp))
        call check_true('rank4 scalar-weight all complex signed modes logabs', &
            maxval(abs(logabs_hyper - base_logabs - log(2.0_dp))) <= 1.0e-13_dp)
        call check_true('rank4 scalar-weight all complex signed modes phase', &
            maxval(abs(phases - base_phases)) <= 8.0e-15_dp)
    end do

    call check_true('rank4 incompatible weights return nan', ieee_is_nan(logsumexp(x, bad)))
    call check_true('rank4 empty broadcast returns minus infinity', &
        .not. ieee_is_nan(logsumexp(x, bzero)) .and. logsumexp(x, bzero) < -huge(0.0_dp))

    print *, 'rank-4 weighted special reduction tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (.not. ieee_is_nan(actual) .and. abs(actual - expected) <= tolerance) return
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

        if (.not. complex_is_nan(actual) .and. abs(actual - expected) <= tolerance) return
        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_shape3(name, actual, expected)
        character(len=*), intent(in) :: name !! check label printed on failure
        integer, intent(in) :: actual(3) !! computed extents
        integer, intent(in) :: expected(3) !! expected extents

        if (any(actual /= expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_shape3

    subroutine check_shape4(name, actual, expected)
        character(len=*), intent(in) :: name !! check label printed on failure
        integer, intent(in) :: actual(4) !! computed extents
        integer, intent(in) :: expected(4) !! expected extents

        if (any(actual /= expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_shape4

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

end program test_special_reductions_rank4_weighted
