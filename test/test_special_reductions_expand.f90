! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_expand
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_special_reductions, only : logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims
    implicit none

    real(dp), parameter :: xcol(2, 1) = reshape([1.0_dp, 2.0_dp], [2, 1])
    real(dp), parameter :: brow(3) = [1.0_dp, 2.0_dp, 4.0_dp]
    real(dp), parameter :: bsigned(3) = [1.0_dp, -2.0_dp, 4.0_dp]
    real(dp), parameter :: xrow(1, 3) = reshape([1.0_dp, 2.0_dp, 3.0_dp], [1, 3])
    real(dp), parameter :: bcol(2, 1) = reshape([1.0_dp, 3.0_dp], [2, 1])
    complex(dp), parameter :: cbrow(3) = [ &
        cmplx(1.0_dp, 0.5_dp, kind=dp), &
        cmplx(-2.0_dp, 0.25_dp, kind=dp), &
        cmplx(3.0_dp, -1.0_dp, kind=dp)]
    complex(dp), parameter :: zcol(2, 1) = reshape([ &
        cmplx(1.0_dp, 0.2_dp, kind=dp), &
        cmplx(2.0_dp, -0.3_dp, kind=dp)], [2, 1])
    real(dp), parameter :: xcube(2, 1, 1) = reshape([1.0_dp, 2.0_dp], [2, 1, 1])
    real(dp), parameter :: bcube(3, 4) = reshape([ &
        1.0_dp, 5.0_dp, 9.0_dp, 2.0_dp, 6.0_dp, 10.0_dp, &
        3.0_dp, 7.0_dp, 11.0_dp, 4.0_dp, 8.0_dp, 12.0_dp], [3, 4])
    real(dp), parameter :: bcube_signed(3, 4) = reshape([ &
        1.0_dp, 5.0_dp, 9.0_dp, -2.0_dp, -6.0_dp, -10.0_dp, &
        3.0_dp, 7.0_dp, 11.0_dp, -4.0_dp, -8.0_dp, -12.0_dp], [3, 4])
    complex(dp), parameter :: zcube(2, 1, 1) = reshape([ &
        cmplx(1.0_dp, 0.2_dp, kind=dp), &
        cmplx(2.0_dp, -0.3_dp, kind=dp)], [2, 1, 1])

    complex(dp), allocatable :: ccube(:, :, :)
    complex(dp), allocatable :: cmat(:, :)
    complex(dp), allocatable :: cphase(:, :)
    complex(dp), allocatable :: cvec(:)
    complex(dp) :: phase
    complex(dp) :: zvalue
    real(dp), allocatable :: cube(:, :, :)
    real(dp), allocatable :: logabs_mat(:, :)
    real(dp), allocatable :: logabs_vec(:)
    real(dp), allocatable :: mat(:, :)
    real(dp), allocatable :: signs(:, :)
    real(dp), allocatable :: sign_vec(:)
    real(dp), allocatable :: vec(:)
    real(dp) :: logabs
    complex(dp) :: cb_cube(3, 4)

    cb_cube = cmplx(bcube, 0.25_dp * bcube, kind=dp)

    call check_close('vector x expansion', logsumexp([2.0_dp], [1.0_dp, 2.0_dp, 3.0_dp]), &
        3.791759469228055_dp, 1.0e-14_dp)

    call check_close('rank2 trailing expansion all', logsumexp(xcol, brow), &
        4.259171836573536_dp, 1.0e-14_dp)
    vec = logsumexp(xcol, 1, brow)
    call check_array('rank2 trailing expansion dim1', vec, [ &
        2.313261687518223_dp, 3.006408868078168_dp, 3.6995560486381134_dp], 1.0e-14_dp)
    vec = logsumexp(xcol, 2, brow)
    call check_array('rank2 trailing expansion dim2', vec, [ &
        2.9459101490553135_dp, 3.9459101490553135_dp], 1.0e-14_dp)

    mat = logsumexp_keepdims(xcol, [1], brow)
    call check_shape('rank2 expanded keepdims shape', shape(mat), [1, 3])
    call check_matrix('rank2 expanded keepdims values', mat, reshape([ &
        2.313261687518223_dp, 3.006408868078168_dp, 3.6995560486381134_dp], [1, 3]), &
        1.0e-14_dp)
    mat = logsumexp_keepdims(xcol, [integer ::], brow)
    call check_shape('rank2 expanded empty axes shape', shape(mat), [2, 3])
    call check_matrix('rank2 expanded empty axes values', mat, reshape([ &
        1.0_dp, 2.0_dp, 1.6931471805599454_dp, 2.6931471805599454_dp, &
        2.386294361119891_dp, 3.386294361119891_dp], [2, 3]), 1.0e-14_dp)

    call logsumexp_signed(xcol, 1, logabs_vec, sign_vec, bsigned)
    call check_array('rank2 expanded signed logabs', logabs_vec, [ &
        2.313261687518223_dp, 3.006408868078168_dp, 3.6995560486381134_dp], 1.0e-14_dp)
    call check_array('rank2 expanded signed sign', sign_vec, [1.0_dp, -1.0_dp, 1.0_dp], 0.0_dp)
    call logsumexp_signed_keepdims(xcol, [integer ::], logabs_mat, signs, bsigned)
    call check_shape('rank2 signed empty axes expanded shape', shape(logabs_mat), [2, 3])
    call check_matrix('rank2 signed empty axes signs', signs, reshape([ &
        1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], [2, 3]), 0.0_dp)

    call check_close('rank2 leading expansion all', logsumexp(xrow, bcol), &
        4.793900325564271_dp, 1.0e-14_dp)
    vec = logsumexp(xrow, 1, bcol)
    call check_array('rank2 leading expansion dim1', vec, [ &
        2.386294361119891_dp, 3.386294361119891_dp, 4.386294361119891_dp], 1.0e-14_dp)
    vec = logsumexp(xrow, 2, bcol)
    call check_array('rank2 leading expansion dim2', vec, [ &
        3.40760596444438_dp, 4.50621825311249_dp], 1.0e-14_dp)

    zvalue = logsumexp(xcol, cbrow)
    call check_complex_close('real x complex b expanded all', zvalue, &
        cmplx(3.0141609613461506_dp, -0.12435499454676144_dp, kind=dp), 1.0e-14_dp)
    cvec = logsumexp(xcol, 1, cbrow)
    call check_complex_array('real x complex b expanded dim1', cvec, [ &
        cmplx(2.424833463175328_dp, 0.4636476090008061_dp, kind=dp), &
        cmplx(3.0141609613461506_dp, 3.017237659043032_dp, kind=dp), &
        cmplx(3.464554234015246_dp, -0.3217505543966422_dp, kind=dp)], 1.5e-14_dp)
    call logsumexp_signed(xcol, 1, logabs_vec, cvec, cbrow)
    call check_array('real x complex b expanded signed logabs', logabs_vec, [ &
        2.424833463175328_dp, 3.0141609613461506_dp, 3.464554234015246_dp], 1.0e-14_dp)
    call check_complex_array('real x complex b expanded phase', cvec, [ &
        cmplx(0.8944271909999159_dp, 0.4472135954999579_dp, kind=dp), &
        cmplx(-0.9922778767136677_dp, 0.12403473458920847_dp, kind=dp), &
        cmplx(0.9486832980505138_dp, -0.31622776601683794_dp, kind=dp)], 2.0e-15_dp)

    zvalue = logsumexp(zcol, brow)
    call check_complex_close('complex x real b expanded all', zvalue, &
        cmplx(4.234504514760136_dp, -0.16745476500169593_dp, kind=dp), 1.5e-14_dp)
    cvec = logsumexp(zcol, 1, brow)
    call check_complex_array('complex x real b expanded dim1', cvec, [ &
        cmplx(2.288594365704823_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(2.981741546264768_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(3.6748887268247135_dp, -0.16745476500169595_dp, kind=dp)], 1.5e-14_dp)
    cmat = logsumexp_keepdims(zcol, [1], brow)
    call check_shape('complex x real b expanded keepdims shape', shape(cmat), [1, 3])
    call check_complex_array('complex x real b expanded keepdims values', reshape(cmat, [3]), [ &
        cmplx(2.288594365704823_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(2.981741546264768_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(3.6748887268247135_dp, -0.16745476500169595_dp, kind=dp)], 1.5e-14_dp)

    zvalue = logsumexp([cmplx(1.0_dp, 0.2_dp, kind=dp)], [ &
        cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(0.0_dp, 2.0_dp, kind=dp)])
    call check_complex_close('complex vector expansion', zvalue, &
        cmplx(1.8047189562170503_dp, 1.3071487177940904_dp, kind=dp), 1.0e-14_dp)
    call logsumexp_signed([cmplx(1.0_dp, 0.2_dp, kind=dp)], logabs, phase, [ &
        cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(0.0_dp, 2.0_dp, kind=dp)])
    call check_close('complex vector expanded signed logabs', logabs, 1.8047189562170503_dp, 1.0e-14_dp)
    call check_complex_close('complex vector expanded signed phase', phase, &
        cmplx(0.26060384662486136_dp, 0.9654458219518719_dp, kind=dp), 2.0e-15_dp)

    call logsumexp_signed_keepdims(zcol, [1], logabs_mat, cphase, brow)
    call check_shape('complex signed expanded keepdims shape', shape(logabs_mat), [1, 3])
    call check_complex_array('complex signed expanded keepdims phase', reshape(cphase, [3]), [ &
        cmplx(0.9860121828642311_dp, -0.16667325893290225_dp, kind=dp), &
        cmplx(0.9860121828642311_dp, -0.16667325893290225_dp, kind=dp), &
        cmplx(0.9860121828642311_dp, -0.16667325893290225_dp, kind=dp)], 2.0e-15_dp)

    call check_close('rank3 trailing expansion all', logsumexp(xcube, bcube), &
        6.669970514207814_dp, 1.0e-14_dp)
    mat = logsumexp(xcube, 1, bcube)
    call check_shape('rank3 trailing expansion dim1 shape', shape(mat), [3, 4])
    call check_matrix('rank3 trailing expansion dim1', mat, reshape([ &
        2.313261687518223_dp, 3.922699599952323_dp, 4.510486264854443_dp, &
        3.006408868078168_dp, 4.105021156746278_dp, 4.615846780512269_dp, &
        3.4118739761863326_dp, 4.259171836573536_dp, 4.711156960316593_dp, &
        3.6995560486381134_dp, 4.392703229198059_dp, 4.798168337306223_dp], [3, 4]), &
        1.5e-14_dp)
    cube = logsumexp_keepdims(xcube, [1, 3], bcube)
    call check_shape('rank3 expanded keepdims shape', shape(cube), [1, 3, 1])
    call check_array('rank3 expanded keepdims values', reshape(cube, [size(cube)]), [ &
        4.615846780512269_dp, 5.571358225539705_dp, 6.050931305801591_dp], 1.5e-14_dp)
    cube = logsumexp_keepdims(xcube, [integer ::], bcube)
    call check_shape('rank3 expanded empty axes shape', shape(cube), [2, 3, 4])
    call check_close('rank3 expanded empty axes first', cube(1, 1, 1), 1.0_dp, 1.0e-14_dp)
    call check_close('rank3 expanded empty axes last', cube(2, 3, 4), &
        4.484906649788_dp, 1.0e-14_dp)

    call logsumexp_signed(xcube, 1, logabs_mat, signs, bcube_signed)
    call check_matrix('rank3 expanded signed logabs', logabs_mat, reshape([ &
        2.313261687518223_dp, 3.922699599952323_dp, 4.510486264854443_dp, &
        3.006408868078168_dp, 4.105021156746278_dp, 4.615846780512269_dp, &
        3.4118739761863326_dp, 4.259171836573536_dp, 4.711156960316593_dp, &
        3.6995560486381134_dp, 4.392703229198059_dp, 4.798168337306223_dp], [3, 4]), &
        1.5e-14_dp)
    call check_matrix('rank3 expanded signed signs', signs, reshape([ &
        1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp], [3, 4]), 0.0_dp)

    zvalue = logsumexp(zcube, bcube)
    call check_complex_close('rank3 complex x real b expanded all', zvalue, &
        cmplx(6.645303192394414_dp, -0.16745476500169593_dp, kind=dp), 2.0e-14_dp)
    cmat = logsumexp(zcube, 1, bcube)
    call check_shape('rank3 complex x real b dim1 shape', shape(cmat), [3, 4])
    call check_complex_array('rank3 complex x real b dim1', reshape(cmat, [size(cmat)]), [ &
        cmplx(2.288594365704823_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(3.898032278138923_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.485818943041043_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(2.981741546264768_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.080353834932878_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.5911794586988695_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(3.3872066543729327_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.234504514760136_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.686489638503193_dp, -0.16745476500169593_dp, kind=dp), &
        cmplx(3.6748887268247135_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.368035907384659_dp, -0.16745476500169595_dp, kind=dp), &
        cmplx(4.7735010154928235_dp, -0.16745476500169595_dp, kind=dp)], 2.0e-14_dp)

    cmat = logsumexp(xcube, 3, cb_cube)
    call check_shape('rank3 real x complex b dim3 shape', shape(cmat), [2, 3])
    call check_complex_array('rank3 real x complex b dim3', reshape(cmat, [size(cmat)]), [ &
        cmplx(3.332897403902263_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(4.332897403902264_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(4.288408848929699_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(5.288408848929699_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(4.767981929191586_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(5.767981929191586_dp, 0.24497866312686414_dp, kind=dp)], 2.0e-14_dp)
    ccube = logsumexp_keepdims(xcube, [1, 3], cb_cube)
    call check_shape('rank3 real x complex b keepdims shape', shape(ccube), [1, 3, 1])
    call check_complex_array('rank3 real x complex b keepdims', reshape(ccube, [size(ccube)]), [ &
        cmplx(4.646159091420486_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(5.601670536447923_dp, 0.24497866312686414_dp, kind=dp), &
        cmplx(6.081243616709808_dp, 0.24497866312686414_dp, kind=dp)], 2.0e-14_dp)

    print *, 'special reduction expanding-broadcast tests passed'

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
        complex(dp), intent(in) :: actual !! computed value
        complex(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum complex absolute error

        if (actual == expected) return
        if (complex_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_complex_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:) !! computed values
        complex(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise complex absolute error

        if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
        if (any_complex_nan(actual)) then
            print *, 'FAIL ', name, ' contains NaN'
            error stop 1
        end if
    end subroutine check_complex_array

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

end program test_special_reductions_expand
