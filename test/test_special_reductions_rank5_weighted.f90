! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank5_weighted
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims
    implicit none

    real(dp), parameter :: x(2, 1, 2, 1, 3) = reshape([ &
        0.2_dp, -0.3_dp, 0.7_dp, -1.1_dp, 1.2_dp, 0.4_dp, &
        -0.8_dp, 0.5_dp, 1.5_dp, -0.2_dp, 0.9_dp, -1.4_dp], [2, 1, 2, 1, 3])
    real(dp), parameter :: zimag(2, 1, 2, 1, 3) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 1.1_dp, -1.2_dp], [2, 1, 2, 1, 3])
    real(dp), parameter :: b1(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    real(dp), parameter :: b2(4, 3) = reshape([ &
        1.0_dp, -0.5_dp, 2.0_dp, 0.25_dp, &
        0.75_dp, 1.5_dp, -1.25_dp, 2.5_dp, &
        -0.4_dp, 1.2_dp, 0.8_dp, 1.1_dp], [4, 3])
    real(dp), parameter :: b3(1, 4, 3) = reshape([ &
        1.0_dp, 0.5_dp, -1.5_dp, 2.0_dp, &
        0.75_dp, -0.25_dp, 1.25_dp, 0.8_dp, &
        -0.6_dp, 1.4_dp, 0.9_dp, 1.1_dp], [1, 4, 3])
    real(dp) :: b4(3, 1, 4, 3)
    real(dp) :: b5(1, 3, 1, 4, 3)
    real(dp) :: b6(1, 1, 1, 1, 1, 1)
    complex(dp) :: cb2(4, 3)
    complex(dp) :: cb5(1, 3, 1, 4, 3)
    complex(dp) :: z(2, 1, 2, 1, 3)
    real(dp) :: xfull(2, 2, 2, 2, 2)
    complex(dp) :: zfull(2, 2, 2, 2, 2)
    real(dp), allocatable :: xempty(:, :, :, :, :)
    real(dp), allocatable :: bempty(:, :, :)
    real(dp), allocatable :: four(:, :, :, :)
    real(dp), allocatable :: kept(:, :, :, :, :)
    real(dp), allocatable :: base_kept(:, :, :, :, :)
    real(dp), allocatable :: signs(:, :, :, :, :)
    real(dp), allocatable :: base_signs(:, :, :, :, :)
    real(dp), allocatable :: logabs_kept(:, :, :, :, :)
    real(dp), allocatable :: base_logabs(:, :, :, :, :)
    complex(dp), allocatable :: cfour(:, :, :, :)
    complex(dp), allocatable :: ckept(:, :, :, :, :)
    complex(dp), allocatable :: base_ckept(:, :, :, :, :)
    complex(dp), allocatable :: phases(:, :, :, :, :)
    complex(dp), allocatable :: base_phases(:, :, :, :, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    complex(dp) :: phase
    complex(dp) :: cvalue
    integer :: d
    integer :: dims(5)
    integer :: i
    integer :: mask
    integer :: nd

    b4 = reshape([(0.1_dp * real(i, dp), i = 1, size(b4))], shape(b4))
    b4(2, 1, 3, 2) = -b4(2, 1, 3, 2)
    b5 = reshape([(real(i, dp) / 13.0_dp, i = 1, size(b5))], shape(b5))
    b5(1, 3, 1, 2, 3) = -b5(1, 3, 1, 2, 3)
    b6 = 1.0_dp
    cb2 = cmplx(b2, 0.2_dp * b2, kind=dp)
    cb2(2, 2) = cb2(2, 2) * cmplx(0.5_dp, -1.1_dp, kind=dp)
    cb5 = b5 * cmplx(0.7_dp, -0.3_dp, kind=dp)
    cb5(1, 1, 1, 1, 1) = cmplx(1.0_dp, 0.5_dp, kind=dp)
    z = cmplx(x, zimag, kind=dp)
    xfull = reshape([(0.05_dp * real(i, dp), i = 1, size(xfull))], shape(xfull))
    zfull = cmplx(xfull, 0.01_dp * xfull, kind=dp)

    call check_close('rank5 scalar weight', logsumexp(x, 0.75_dp), &
        2.668553189015345_dp, 4.0e-14_dp)
    call logsumexp_signed(x, logabs, sign_value, b1)
    call check_close('rank5 rank1 signed logabs', logabs, &
        1.7061705330210652_dp, 4.0e-14_dp)
    call check_close('rank5 rank1 signed sign', sign_value, -1.0_dp, 0.0_dp)
    call check_true('rank5 rank1 negative plain is nan', ieee_is_nan(logsumexp(x, b1)))
    call check_close('rank5 rank2 expanding weight', logsumexp(x, b2), &
        4.054427858092769_dp, 4.0e-14_dp)
    call check_close('rank5 rank3 expanding weight', logsumexp(x, b3), &
        3.8847896097559023_dp, 4.0e-14_dp)
    call check_close('rank5 rank4 expanding weight', logsumexp(x, b4), &
        6.11459576516539_dp, 5.0e-14_dp)
    call check_close('rank5 rank5 expanding weight', logsumexp(x, b5), &
        5.806092289236909_dp, 5.0e-14_dp)

    four = logsumexp(x, dim=4, b=b2)
    call check_shape4('rank5 weighted dim4 shape', shape(four), [2, 1, 2, 3])
    call check_close('rank5 weighted dim4 first', four(1, 1, 1, 1), &
        1.2116009116784798_dp, 4.0e-14_dp)
    call check_close('rank5 weighted dim4 last', four(2, 1, 2, 3), &
        -0.40674822698971647_dp, 4.0e-14_dp)

    kept = logsumexp_keepdims(x, [1, 3, 5], b2)
    call check_shape5('rank5 weighted keepdims shape', shape(kept), [1, 1, 1, 4, 1])
    call check_close('rank5 weighted keepdims first', kept(1, 1, 1, 1, 1), &
        1.8387528935602484_dp, 4.0e-14_dp)
    call check_close('rank5 weighted keepdims last', kept(1, 1, 1, 4, 1), &
        3.3017258886160565_dp, 4.0e-14_dp)
    call logsumexp_signed_keepdims(x, [1, 3, 5], logabs_kept, signs, b2)
    call check_close('rank5 signed keepdims middle', logabs_kept(1, 1, 1, 2, 1), &
        2.880237763998256_dp, 4.0e-14_dp)
    call check_true('rank5 signed keepdims signs', &
        maxval(abs(signs - 1.0_dp)) <= epsilon(1.0_dp))

    cvalue = logsumexp(x, cb5)
    call check_complex_close('rank5 real input rank5 complex weights', cvalue, &
        cmplx(5.545066558970053_dp, -0.39052421562009276_dp, dp), 8.0e-14_dp)

    cvalue = logsumexp(x, cb2)
    call check_complex_close('rank5 real input complex weights', cvalue, &
        cmplx(4.002914220755019_dp, -0.016574776805369146_dp, dp), 6.0e-14_dp)
    cfour = logsumexp(x, 4, cb2)
    call check_complex_close('rank5 real complex-weight dim4 first', cfour(1, 1, 1, 1), &
        cmplx(1.2312112682551206_dp, 0.19739555984988078_dp, dp), 6.0e-14_dp)
    call check_complex_close('rank5 real complex-weight dim4 last', cfour(2, 1, 2, 3), &
        cmplx(-0.3871378704130757_dp, 0.19739555984988075_dp, dp), 6.0e-14_dp)
    call logsumexp_signed(x, logabs, phase, cb2)
    call check_close('rank5 real complex-weight signed logabs', logabs, &
        4.002914220755019_dp, 6.0e-14_dp)
    call check_complex_close('rank5 real complex-weight signed phase', phase, &
        cmplx(0.9998626415315994_dp, -0.016574017903100854_dp, dp), 8.0e-15_dp)

    ckept = logsumexp_keepdims(x, [1, 3, 5], cb2)
    call check_complex_close('rank5 real complex-weight keepdims first', &
        ckept(1, 1, 1, 1, 1), &
        cmplx(1.858363250136889_dp, 0.19739555984988044_dp, dp), 8.0e-14_dp)
    call check_complex_close('rank5 real complex-weight keepdims last', &
        ckept(1, 1, 1, 4, 1), &
        cmplx(3.321336245192697_dp, 0.19739555984988078_dp, dp), 8.0e-14_dp)
    call logsumexp_signed_keepdims(x, [1, 3, 5], logabs_kept, phases, cb2)
    call check_close('rank5 real complex-weight keepdims logabs', &
        logabs_kept(1, 1, 1, 2, 1), 2.853969067009297_dp, 8.0e-14_dp)
    call check_complex_close('rank5 real complex-weight keepdims phase', &
        phases(1, 1, 1, 2, 1), &
        cmplx(0.859406158981976_dp, -0.5112935105238929_dp, dp), 8.0e-15_dp)

    cvalue = logsumexp(z, 2.0_dp)
    call check_complex_close('rank5 complex input scalar real weight', cvalue, &
        cmplx(3.409246634565723_dp, 0.3426330572047245_dp, dp), 7.0e-14_dp)
    cvalue = logsumexp(z, b5)
    call check_complex_close('rank5 complex input rank5 real weights', cvalue, &
        cmplx(5.530493210605167_dp, 0.46352320023748_dp, dp), 8.0e-14_dp)

    cvalue = logsumexp(z, b2)
    call check_complex_close('rank5 complex input real weights', cvalue, &
        cmplx(3.812668917271383_dp, 0.3056718017047606_dp, dp), 6.0e-14_dp)
    cfour = logsumexp(z, 4, b2)
    call check_complex_close('rank5 complex real-weight dim4 first', cfour(1, 1, 1, 1), &
        cmplx(1.21160091167848_dp, 0.1_dp, dp), 6.0e-14_dp)
    call check_complex_close('rank5 complex real-weight dim4 last', cfour(2, 1, 2, 3), &
        cmplx(-0.40674822698971647_dp, -1.2_dp, dp), 6.0e-14_dp)

    ckept = logsumexp_keepdims(z, [1, 3, 5], b2)
    call check_complex_close('rank5 complex real-weight keepdims first', &
        ckept(1, 1, 1, 1, 1), &
        cmplx(1.9187527245319118_dp, -0.2342404020200135_dp, dp), 8.0e-14_dp)
    call check_complex_close('rank5 complex real-weight keepdims last', &
        ckept(1, 1, 1, 4, 1), &
        cmplx(3.0267910465958767_dp, 0.24542570447596912_dp, dp), 8.0e-14_dp)
    call logsumexp_signed_keepdims(z, [1, 3, 5], logabs_kept, phases, b2)
    call check_close('rank5 complex real-weight keepdims logabs', &
        logabs_kept(1, 1, 1, 3, 1), 1.928743757773339_dp, 8.0e-14_dp)
    call check_complex_close('rank5 complex real-weight keepdims phase', &
        phases(1, 1, 1, 3, 1), &
        cmplx(0.7039122086743662_dp, 0.7102869859987407_dp, dp), 8.0e-15_dp)

    cvalue = logsumexp(z, cb5)
    call check_complex_close('rank5 complex input complex weights', cvalue, &
        cmplx(5.278109348878839_dp, 0.07075635795506374_dp, dp), 7.0e-14_dp)
    call logsumexp_signed(z, logabs, phase, cb5)
    call check_close('rank5 complex-complex signed logabs', logabs, &
        5.278109348878839_dp, 7.0e-14_dp)
    call check_complex_close('rank5 complex-complex signed phase', phase, &
        cmplx(0.9974978130911839_dp, 0.07069733289385814_dp, dp), 8.0e-15_dp)

    ckept = logsumexp_keepdims(z, [1, 3, 5], cb2)
    call check_shape5('rank5 complex keepdims shape', shape(ckept), [1, 1, 1, 4, 1])
    call check_complex_close('rank5 complex keepdims first', ckept(1, 1, 1, 1, 1), &
        cmplx(1.9383630811085526_dp, -0.03684484217013384_dp, dp), 8.0e-14_dp)
    call check_complex_close('rank5 complex keepdims last', ckept(1, 1, 1, 4, 1), &
        cmplx(3.046401403172517_dp, 0.4428212643258498_dp, dp), 8.0e-14_dp)

    do mask = 0, 31
        nd = 0
        do d = 1, 5
            if (btest(mask, d - 1)) then
                nd = nd + 1
                dims(nd) = d
            end if
        end do

        base_kept = logsumexp_keepdims(xfull, dims(:nd))
        kept = logsumexp_keepdims(xfull, dims(:nd), 2.0_dp)
        call check_true('rank5 scalar-weight all real keepdims modes shape', &
            all(shape(kept) == shape(base_kept)))
        call check_true('rank5 scalar-weight all real keepdims modes value', &
            maxval(abs(kept - base_kept - log(2.0_dp))) <= 1.2e-13_dp)
        call logsumexp_signed_keepdims(xfull, dims(:nd), base_logabs, base_signs)
        call logsumexp_signed_keepdims(xfull, dims(:nd), logabs_kept, signs, 2.0_dp)
        call check_true('rank5 scalar-weight all real signed modes logabs', &
            maxval(abs(logabs_kept - base_logabs - log(2.0_dp))) <= 1.2e-13_dp)
        call check_true('rank5 scalar-weight all real signed modes sign', &
            maxval(abs(signs - base_signs)) <= 2.0e-15_dp)

        base_ckept = logsumexp_keepdims(zfull, dims(:nd))
        ckept = logsumexp_keepdims(zfull, dims(:nd), cmplx(2.0_dp, 0.0_dp, kind=dp))
        call check_true('rank5 scalar-weight all complex keepdims modes shape', &
            all(shape(ckept) == shape(base_ckept)))
        call check_true('rank5 scalar-weight all complex keepdims modes value', &
            maxval(abs(ckept - base_ckept - log(2.0_dp))) <= 1.5e-13_dp)
        call logsumexp_signed_keepdims(zfull, dims(:nd), base_logabs, base_phases)
        call logsumexp_signed_keepdims(zfull, dims(:nd), logabs_kept, phases, &
            cmplx(2.0_dp, 0.0_dp, kind=dp))
        call check_true('rank5 scalar-weight all complex signed modes logabs', &
            maxval(abs(logabs_kept - base_logabs - log(2.0_dp))) <= 1.5e-13_dp)
        call check_true('rank5 scalar-weight all complex signed modes phase', &
            maxval(abs(phases - base_phases)) <= 1.0e-14_dp)
    end do

    call check_true('rank5 rank6 weights return nan', ieee_is_nan(logsumexp(x, b6)))
    cvalue = logsumexp(z, b6)
    call check_true('rank5 complex rank6 real weights return nan', complex_is_nan(cvalue))
    ckept = logsumexp_keepdims(z, [1, 3], b6)
    call check_shape5('rank5 complex rank6 keepdims shape', shape(ckept), [1, 1, 1, 1, 3])
    call check_true('rank5 complex rank6 keepdims nan', complex_is_nan(ckept(1, 1, 1, 1, 1)))
    allocate(xempty(2, 1, 0, 1, 3), bempty(1, 4, 3))
    bempty = 1.0_dp
    call logsumexp_signed(xempty, logabs, sign_value, bempty)
    call check_true('rank5 zero-extent weighted logabs is minus infinity', &
        .not. ieee_is_nan(logabs) .and. logabs < -huge(0.0_dp))
    call check_close('rank5 zero-extent weighted sign', sign_value, 0.0_dp, 0.0_dp)

    print *, 'rank-5 weighted special reduction tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (.not. ieee_is_nan(actual) .and. abs(actual - expected) <= tolerance) return
        print *, 'FAIL ', name
        print *, ' actual  = ', actual
        print *, ' expected= ', expected
        print *, ' tolerance=', tolerance
        error stop 1
    end subroutine check_close

    subroutine check_complex_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum complex absolute error

        if (.not. complex_is_nan(actual) .and. abs(actual - expected) <= tolerance) return
        print *, 'FAIL ', name
        print *, ' actual  = ', actual
        print *, ' expected= ', expected
        print *, ' tolerance=', tolerance
        error stop 1
    end subroutine check_complex_close

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

    subroutine check_shape5(name, actual, expected)
        character(len=*), intent(in) :: name !! check label printed on failure
        integer, intent(in) :: actual(5) !! computed extents
        integer, intent(in) :: expected(5) !! expected extents

        if (any(actual /= expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_shape5

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

end program test_special_reductions_rank5_weighted
