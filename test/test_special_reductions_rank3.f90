! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_rank3
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    real(dp), parameter :: x(2, 2, 3) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
        7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp], [2, 2, 3])
    real(dp), parameter :: bvec(3) = [1.0_dp, -0.5_dp, 2.0_dp]
    real(dp), parameter :: bmid(2, 1) = reshape([1.0_dp, 2.0_dp], [2, 1])
    real(dp), parameter :: b3(1, 2, 1) = reshape([1.0_dp, -2.0_dp], [1, 2, 1])
    complex(dp), parameter :: z(2, 2, 3) = reshape([ &
        cmplx(1.0_dp, 0.1_dp, kind=dp), cmplx(2.0_dp, -0.2_dp, kind=dp), &
        cmplx(3.0_dp, 0.3_dp, kind=dp), cmplx(4.0_dp, -0.4_dp, kind=dp), &
        cmplx(5.0_dp, 0.5_dp, kind=dp), cmplx(6.0_dp, -0.6_dp, kind=dp), &
        cmplx(7.0_dp, 0.7_dp, kind=dp), cmplx(8.0_dp, -0.8_dp, kind=dp), &
        cmplx(9.0_dp, 0.9_dp, kind=dp), cmplx(10.0_dp, -1.0_dp, kind=dp), &
        cmplx(11.0_dp, 1.1_dp, kind=dp), cmplx(12.0_dp, -1.2_dp, kind=dp)], [2, 2, 3])
    complex(dp), parameter :: cbvec(3) = [ &
        cmplx(1.0_dp, 0.5_dp, kind=dp), cmplx(-0.5_dp, 0.25_dp, kind=dp), &
        cmplx(2.0_dp, -0.75_dp, kind=dp)]

    complex(dp), allocatable :: ccube(:, :, :)
    complex(dp), allocatable :: cmatrix(:, :)
    complex(dp) :: phase
    complex(dp) :: zvalue
    real(dp), allocatable :: cube(:, :, :)
    real(dp), allocatable :: logabs_cube(:, :, :)
    real(dp), allocatable :: logabs_matrix(:, :)
    real(dp), allocatable :: matrix(:, :)
    real(dp), allocatable :: signs(:, :, :)
    real(dp), allocatable :: sign_matrix(:, :)
    real(dp) :: logabs
    real(dp) :: sign_value
    integer :: i
    integer :: j

    call check_close('rank3 weighted all', logsumexp(x, bvec), &
        13.128915942758809_dp, 1.5e-14_dp)
    call check_close('rank3 scalar broadcast', logsumexp(x, 2.0_dp), &
        13.151816181715798_dp, 1.5e-14_dp)
    call logsumexp_signed(x, logabs, sign_value, b3)
    call check_close('rank3 singleton broadcast logabs', logabs, &
        12.954822250512544_dp, 1.5e-14_dp)
    call check_close('rank3 singleton broadcast sign', sign_value, -1.0_dp, 0.0_dp)

    call logsumexp_signed(x, 1, logabs_matrix, sign_matrix, bvec)
    call check_shape('rank3 dim1 shape', shape(logabs_matrix), [2, 3])
    call check_matrix('rank3 dim1 signed logabs', logabs_matrix, reshape([ &
        2.313261687518223_dp, 4.313261687518223_dp, &
        5.620114506958277_dp, 7.620114506958277_dp, &
        11.006408868078168_dp, 13.006408868078168_dp], [2, 3]), 1.5e-14_dp)
    call check_matrix('rank3 dim1 signs', sign_matrix, reshape([ &
        1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], [2, 3]), 0.0_dp)

    matrix = logsumexp(x, 3, bvec)
    call check_matrix('rank3 dim3 reduction', matrix, reshape([ &
        9.688726244197614_dp, 10.688726244197614_dp, &
        11.688726244197614_dp, 12.688726244197614_dp], [2, 2]), 1.5e-14_dp)

    cube = logsumexp_keepdims(x, [1, 3], bmid)
    call check_shape('rank3 keepdims pair shape', shape(cube), [1, 2, 1])
    call check_close('rank3 keepdims pair first', cube(1, 1, 1), &
        10.33174099011288_dp, 1.5e-14_dp)
    call check_close('rank3 keepdims pair second', cube(1, 2, 1), &
        13.024888170672826_dp, 1.5e-14_dp)
    call check_cube('rank3 reversed keepdims axes', logsumexp_keepdims(x, [3, 1], bmid), &
        cube, 1.0e-15_dp)

    call logsumexp_signed_keepdims(x, [integer ::], logabs_cube, signs, b3)
    call check_shape('rank3 empty keepdims shape', shape(logabs_cube), shape(x))
    call check_close('rank3 empty keepdims logabs', logabs_cube(2, 2, 3), &
        12.693147180559945_dp, 1.0e-14_dp)
    call check_close('rank3 empty keepdims negative sign', signs(1, 2, 2), -1.0_dp, 0.0_dp)

    cube = logsumexp_keepdims(x, [1, 1])
    call check_true('rank3 duplicate axes invalid', size(cube) == 0)
    matrix = logsumexp(x, 1, reshape([1.0_dp, 2.0_dp, 3.0_dp], [3, 1]))
    call check_true('rank3 broadcast mismatch', all(ieee_is_nan(matrix)))

    cube = softmax(x, [1, 3])
    do j = 1, size(x, 2)
        call check_close('rank3 tuple softmax normalization', sum(cube(:, j, :)), &
            1.0_dp, 5.0e-16_dp)
    end do
    cube = softmax(x, [integer ::])
    call check_close('rank3 empty softmax', maxval(abs(cube - 1.0_dp)), 0.0_dp, 0.0_dp)

    cube = log_softmax(x, [2, 3])
    do i = 1, size(x, 1)
        call check_close('rank3 tuple log_softmax normalization', sum(exp(cube(i, :, :))), &
            1.0_dp, 5.0e-16_dp)
    end do
    call check_close('rank3 log_softmax reference', cube(1, 1, 1), &
        -10.14540731363763_dp, 1.5e-14_dp)

    zvalue = logsumexp(x, cbvec)
    call check_complex_close('real rank3 complex weights', zvalue, &
        cmplx(13.194466820674249_dp, -0.35813725758065457_dp, kind=dp), 2.0e-14_dp)
    call logsumexp_signed(x, logabs, phase, cbvec)
    call check_close('real rank3 complex signed logabs', logabs, &
        13.194466820674249_dp, 2.0e-14_dp)
    call check_complex_close('real rank3 complex signed phase', phase, &
        cmplx(0.936551395764728_dp, -0.3505302884076351_dp, kind=dp), 2.0e-15_dp)

    zvalue = logsumexp(z, bvec)
    call check_complex_close('complex rank3 real weights', zvalue, &
        cmplx(12.613480507321034_dp, -0.8224100584555297_dp, kind=dp), 2.0e-14_dp)

    zvalue = logsumexp(z, cbvec)
    call check_complex_close('complex rank3 weighted all', zvalue, &
        cmplx(12.678594480020225_dp, -1.1805204345092404_dp, kind=dp), 2.0e-14_dp)
    call logsumexp_signed(z, logabs, phase, cbvec)
    call check_close('complex rank3 signed logabs', logabs, 12.678594480020225_dp, 2.0e-14_dp)
    call check_complex_close('complex rank3 signed phase', phase, &
        cmplx(0.38044357592515676_dp, -0.9248041336073708_dp, kind=dp), 1.0e-15_dp)

    cmatrix = logsumexp(z, 2, cbvec)
    call check_complex_matrix('complex rank3 dim2', cmatrix, reshape([ &
        cmplx(3.2364025128240996_dp, 0.7399136964176286_dp, kind=dp), &
        cmplx(4.2364025128241_dp, 0.08738152158398355_dp, kind=dp), &
        cmplx(6.543255332264154_dp, -2.928974175173777_dp, kind=dp), &
        cmplx(7.543255332264154_dp, 1.9016789571721644_dp, kind=dp), &
        cmplx(11.8837660966213_dp, 0.7174954171462504_dp, kind=dp), &
        cmplx(12.8837660966213_dp, -1.5350367576873947_dp, kind=dp)], [2, 3]), 2.0e-14_dp)

    ccube = logsumexp_keepdims(z, [1, 3], cbvec)
    call check_shape('complex rank3 keepdims shape', shape(ccube), [1, 2, 1])
    call check_complex_close('complex rank3 keepdims first', ccube(1, 1, 1), &
        cmplx(10.69887295900714_dp, -0.9834180942723016_dp, kind=dp), 2.0e-14_dp)
    call check_complex_close('complex rank3 keepdims second', ccube(1, 2, 1), &
        cmplx(12.533556274712879_dp, -1.2117923201380876_dp, kind=dp), 2.0e-14_dp)

    ccube = logsumexp_keepdims(x, [1, 3], cbvec)
    call check_complex_close('real rank3 complex keepdims first', ccube(1, 1, 1), &
        cmplx(11.067538809631277_dp, -0.35813725758065457_dp, kind=dp), 2.0e-14_dp)
    call check_complex_close('real rank3 complex keepdims second', ccube(1, 2, 1), &
        cmplx(13.067538809631277_dp, -0.35813725758065457_dp, kind=dp), 2.0e-14_dp)

    call logsumexp_signed_keepdims(z, [1, 3], logabs_cube, ccube, bvec)
    call check_close('complex rank3 real keepdims logabs first', logabs_cube(1, 1, 1), &
        10.633656712638292_dp, 2.0e-14_dp)
    call check_complex_close('complex rank3 real keepdims phase first', ccube(1, 1, 1), &
        cmplx(0.8107749022709215_dp, -0.5853580595221847_dp, kind=dp), 2.0e-15_dp)

    ccube = softmax(z, [1, 3])
    do j = 1, size(z, 2)
        call check_complex_close('complex rank3 softmax normalization', sum(ccube(:, j, :)), &
            cmplx(1.0_dp, 0.0_dp, kind=dp), 2.0e-15_dp)
    end do
    ccube = log_softmax(z, [2, 3])
    do i = 1, size(z, 1)
        call check_complex_close('complex rank3 log_softmax normalization', &
            sum(exp(ccube(i, :, :))), cmplx(1.0_dp, 0.0_dp, kind=dp), 2.0e-15_dp)
    end do

    print *, 'rank-3 special reduction tests passed'

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

    subroutine check_matrix(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :) !! computed matrix
        real(dp), intent(in) :: expected(:, :) !! reference matrix
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_matrix

    subroutine check_cube(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :, :) !! computed rank-3 array
        real(dp), intent(in) :: expected(:, :, :) !! reference rank-3 array
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_cube

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

    subroutine check_complex_matrix(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:, :) !! computed matrix
        complex(dp), intent(in) :: expected(:, :) !! reference matrix
        real(dp), intent(in) :: tolerance !! maximum elementwise complex absolute error

        if (any(shape(actual) /= shape(expected)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_complex_matrix

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

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition that must hold

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_reductions_rank3
