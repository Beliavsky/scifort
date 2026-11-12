! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_broadcast
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, ieee_value
    use scifort_special, only : dp, log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    implicit none

    real(dp), parameter :: x(2, 3) = reshape([ &
        1.0_dp, 4.0_dp, 2.0_dp, 5.0_dp, 3.0_dp, 6.0_dp], [2, 3])
    real(dp), parameter :: bvec(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    real(dp), parameter :: bcol(2, 1) = reshape([1.0_dp, 2.0_dp], [2, 1])
    real(dp), parameter :: brow(1, 3) = reshape([1.0_dp, -2.0_dp, 0.5_dp], [1, 3])
    complex(dp), parameter :: z(2, 3) = reshape([ &
        cmplx(1.0_dp, 0.2_dp, kind=dp), cmplx(4.0_dp, -1.0_dp, kind=dp), &
        cmplx(2.0_dp, -0.5_dp, kind=dp), cmplx(5.0_dp, 0.3_dp, kind=dp), &
        cmplx(3.0_dp, 1.0_dp, kind=dp), cmplx(6.0_dp, -0.7_dp, kind=dp)], [2, 3])
    complex(dp), parameter :: cbvec(3) = [ &
        cmplx(1.0_dp, 0.5_dp, kind=dp), cmplx(-2.0_dp, 0.25_dp, kind=dp), &
        cmplx(0.5_dp, -0.75_dp, kind=dp)]

    complex(dp), allocatable :: cmat(:, :)
    complex(dp), allocatable :: cvector(:)
    complex(dp), allocatable :: phases(:, :)
    complex(dp), allocatable :: phase_vector(:)
    complex(dp) :: phase
    complex(dp) :: zvalue
    real(dp), allocatable :: logabs_matrix(:, :)
    real(dp), allocatable :: logabs_vector(:)
    real(dp), allocatable :: matrix(:, :)
    real(dp), allocatable :: sign_vector(:)
    real(dp), allocatable :: signs(:, :)
    real(dp), allocatable :: vector(:)
    real(dp) :: logabs

    call check_close('vector scalar broadcast', logsumexp([1.0_dp, 2.0_dp], 2.0_dp), &
        3.006408868078168_dp, 4.0e-15_dp)
    call check_close('vector singleton broadcast', logsumexp([1.0_dp, 2.0_dp], [2.0_dp]), &
        3.006408868078168_dp, 4.0e-15_dp)

    call logsumexp_signed(x, 1, logabs_vector, sign_vector, bvec)
    call check_array('trailing-vector broadcast dim1 logabs', logabs_vector, [ &
        4.048587351573742_dp, 5.741734532133687_dp, 5.355440171013797_dp], 6.0e-15_dp)
    call check_array('trailing-vector broadcast dim1 sign', sign_vector, &
        [1.0_dp, -1.0_dp, 1.0_dp], 0.0_dp)

    call logsumexp_signed(x, 1, logabs_vector, sign_vector, brow)
    call check_array('row broadcast matches vector', logabs_vector, [ &
        4.048587351573742_dp, 5.741734532133687_dp, 5.355440171013797_dp], 6.0e-15_dp)
    call check_array('row broadcast signs', sign_vector, [1.0_dp, -1.0_dp, 1.0_dp], 0.0_dp)

    call check_close('column broadcast all', logsumexp(x, bcol), &
        7.125341883111083_dp, 6.0e-15_dp)
    call check_true('matrix broadcast mismatch', &
        ieee_is_nan(logsumexp(x, reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2]))))
    call check_true('matrix rank3 weight invalid', ieee_is_nan(logsumexp(x, &
        reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1, 2, 2]))))

    matrix = logsumexp_keepdims(x, [1], 2.0_dp)
    call check_shape('keepdims dim1 shape', shape(matrix), [1, 3])
    call check_matrix('keepdims dim1 values', matrix, reshape([ &
        4.741734532133687_dp, 5.741734532133687_dp, 6.741734532133687_dp], [1, 3]), &
        6.0e-15_dp)

    matrix = logsumexp_keepdims(x, [2], 2.0_dp)
    call check_shape('keepdims dim2 shape', shape(matrix), [2, 1])
    call check_matrix('keepdims dim2 values', matrix, reshape([ &
        4.100753145004326_dp, 7.100753145004326_dp], [2, 1]), 6.0e-15_dp)

    matrix = logsumexp_keepdims(x, [1, 2], 2.0_dp)
    call check_shape('keepdims both shape', shape(matrix), [1, 1])
    call check_close('keepdims both value', matrix(1, 1), 7.149340496578068_dp, 6.0e-15_dp)
    matrix = logsumexp_keepdims(x, [2, 1], 2.0_dp)
    call check_close('keepdims reversed axes', matrix(1, 1), 7.149340496578068_dp, 6.0e-15_dp)

    call logsumexp_signed_keepdims(x, [integer ::], logabs_matrix, signs, bvec)
    call check_shape('keepdims empty shape', shape(logabs_matrix), shape(x))
    call check_matrix('keepdims empty logabs', logabs_matrix, reshape([ &
        1.0_dp, 4.0_dp, 2.6931471805599453_dp, 5.6931471805599453_dp, &
        2.3068528194400546_dp, 5.306852819440055_dp], [2, 3]), 6.0e-15_dp)
    call check_matrix('keepdims empty signs', signs, reshape([ &
        1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], [2, 3]), 0.0_dp)

    matrix = logsumexp_keepdims(x, [1, 1])
    call check_true('duplicate keepdims axes invalid', size(matrix) == 0)

    matrix = softmax(x, [1, 2])
    call check_matrix('softmax tuple matches all', matrix, softmax(x), 3.0e-16_dp)
    call check_close('softmax tuple normalization', sum(matrix), 1.0_dp, 3.0e-16_dp)
    matrix = softmax(x, [integer ::])
    call check_matrix('softmax empty axes', matrix, spread([1.0_dp, 1.0_dp], 2, 3), 0.0_dp)
    matrix = log_softmax(x, [integer ::])
    call check_close('log_softmax empty axes', maxval(abs(matrix)), 0.0_dp, 0.0_dp)
    matrix = softmax(reshape([1.0_dp, ieee_value(0.0_dp, ieee_positive_inf)], [1, 2]), &
        [integer ::])
    call check_close('softmax empty finite element', matrix(1, 1), 1.0_dp, 0.0_dp)
    call check_true('softmax empty infinite element', ieee_is_nan(matrix(1, 2)))
    call logsumexp_signed_keepdims(reshape([ieee_value(0.0_dp, ieee_positive_inf), &
        -ieee_value(0.0_dp, ieee_positive_inf)], [1, 2]), [integer ::], &
        logabs_matrix, signs)
    call check_close('empty axes positive infinity sign', signs(1, 1), 1.0_dp, 0.0_dp)
    call check_close('empty axes negative infinity sign', signs(1, 2), 0.0_dp, 0.0_dp)
    matrix = log_softmax(x, [2, 1])
    call check_matrix('log_softmax tuple matches all', matrix, log_softmax(x), 3.0e-15_dp)
    matrix = softmax(x, [1, 1])
    call check_true('softmax duplicate axes invalid', all(ieee_is_nan(matrix)))

    zvalue = logsumexp([1.0_dp, 2.0_dp], [ &
        cmplx(1.0_dp, 1.0_dp, kind=dp), cmplx(2.0_dp, -0.5_dp, kind=dp)])
    call check_complex_close('real input complex weights', zvalue, &
        cmplx(2.8635490382699356_dp, -0.05573919920320018_dp, kind=dp), 6.0e-15_dp)
    call logsumexp_signed([1.0_dp, 2.0_dp], logabs, phase, [ &
        cmplx(1.0_dp, 1.0_dp, kind=dp), cmplx(2.0_dp, -0.5_dp, kind=dp)])
    call check_close('real input complex weights logabs', logabs, 2.8635490382699356_dp, 6.0e-15_dp)
    call check_complex_close('real input complex weights phase', phase, &
        cmplx(0.9984469729848049_dp, -0.05571034138721724_dp, kind=dp), 7.0e-16_dp)

    call logsumexp_signed(z, logabs, phase, cbvec)
    call check_close('complex vector broadcast all logabs', logabs, 6.2409852731097555_dp, 7.0e-15_dp)
    call check_complex_close('complex vector broadcast all phase', phase, &
        cmplx(-0.5343037075379818_dp, -0.8452925813652729_dp, kind=dp), 8.0e-16_dp)

    call logsumexp_signed(z, logabs, phase, bcol)
    call check_close('complex real-column broadcast logabs', logabs, 7.007159956976995_dp, 8.0e-15_dp)
    call check_complex_close('complex real-column broadcast phase', phase, &
        cmplx(0.8869858952624203_dp, -0.4617965153674532_dp, kind=dp), 8.0e-16_dp)

    cmat = logsumexp_keepdims(z, [1, 2], cbvec)
    call check_shape('complex keepdims shape', shape(cmat), [1, 1])
    call check_complex_close('complex keepdims value', cmat(1, 1), &
        cmplx(6.2409852731097555_dp, -2.134480129745498_dp, kind=dp), 8.0e-15_dp)

    call logsumexp_signed_keepdims(z, [integer ::], logabs_matrix, phases, bcol)
    call check_shape('complex empty keepdims shape', shape(logabs_matrix), shape(z))
    call check_close('complex empty keepdims logabs 21', logabs_matrix(2, 1), &
        4.693147180559945_dp, 5.0e-15_dp)
    call check_complex_close('complex empty keepdims phase 11', phases(1, 1), &
        cmplx(0.9800665778412416_dp, 0.19866933079506122_dp, kind=dp), 7.0e-16_dp)
    call check_complex_close('complex empty keepdims phase 23', phases(2, 3), &
        cmplx(0.7648421872844885_dp, -0.6442176872376910_dp, kind=dp), 8.0e-16_dp)

    cmat = logsumexp_keepdims(x, [1, 2], cmplx(2.0_dp, -1.0_dp, kind=dp))
    call check_complex_close('real matrix complex scalar keepdims', cmat(1, 1), &
        cmplx(7.260912272235172_dp, -0.4636476090008061_dp, kind=dp), 8.0e-15_dp)

    cmat = softmax(z, [1, 2])
    call check_complex_close('complex tuple softmax normalization', sum(cmat), &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 1.0e-15_dp)
    cmat = softmax(z, [integer ::])
    call check_complex_matrix_constant('complex empty softmax', cmat, &
        cmplx(1.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    cmat = log_softmax(z, [integer ::])
    call check_complex_matrix_constant('complex empty log_softmax', cmat, &
        cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)

    cvector = logsumexp(z, 1, cbvec)
    call check_complex_close('complex broadcast dim1 first', cvector(1), &
        cmplx(4.130489447242830_dp, -0.49080274286065895_dp, kind=dp), 1.0e-14_dp)
    call logsumexp_signed(z, 1, logabs_vector, phase_vector, cbvec)
    call check_array('complex broadcast dim1 logabs', logabs_vector, [ &
        4.130489447242830_dp, 5.735593607017857_dp, 5.890977932869935_dp], 1.0e-14_dp)

    print *, 'special reduction broadcasting/keepdims tests passed'

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

    subroutine check_array(name, actual, expected, tolerance)
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
    end subroutine check_array

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

    subroutine check_complex_matrix_constant(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:, :) !! computed matrix
        complex(dp), intent(in) :: expected !! expected value for every element
        real(dp), intent(in) :: tolerance !! maximum complex absolute error

        if (any_complex_nan(reshape(actual, [size(actual)])) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_complex_matrix_constant

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

    pure function any_complex_nan(values) result(has_nan)
        complex(dp), intent(in) :: values(:) !! complex values inspected componentwise
        logical :: has_nan
        integer :: i

        has_nan = .false.
        do i = 1, size(values)
            if (complex_is_nan(values(i))) then
                has_nan = .true.
                return
            end if
        end do
    end function any_complex_nan

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition that must hold

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_reductions_broadcast
