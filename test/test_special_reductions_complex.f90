! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_reductions_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, log_softmax, logsumexp, logsumexp_signed, softmax
    implicit none

    complex(dp), parameter :: x(3) = [ &
        cmplx(1000.0_dp, 1.0_dp, kind=dp), &
        cmplx(999.0_dp, 2.0_dp, kind=dp), &
        cmplx(-1000.0_dp, -3.0_dp, kind=dp)]
    real(dp), parameter :: real_weights(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    complex(dp), parameter :: complex_weights(3) = [ &
        cmplx(1.0_dp, 1.0_dp, kind=dp), &
        cmplx(-2.0_dp, 0.5_dp, kind=dp), &
        cmplx(0.3_dp, -0.2_dp, kind=dp)]
    complex(dp), parameter :: matrix(2, 3) = reshape([ &
        cmplx(1.0_dp, 0.2_dp, kind=dp), cmplx(2.0_dp, -0.5_dp, kind=dp), &
        cmplx(0.5_dp, -1.0_dp, kind=dp), cmplx(4.0_dp, 0.7_dp, kind=dp), &
        cmplx(3.0_dp, 2.0_dp, kind=dp), cmplx(-1.0_dp, 1.5_dp, kind=dp)], [2, 3])
    real(dp), parameter :: matrix_real_weights(2, 3) = reshape([ &
        1.0_dp, 0.2_dp, -2.0_dp, 1.0_dp, 0.5_dp, -1.0_dp], [2, 3])
    complex(dp), parameter :: matrix_complex_weights(2, 3) = reshape([ &
        cmplx(1.0_dp, 0.1_dp, kind=dp), cmplx(0.2_dp, -0.4_dp, kind=dp), &
        cmplx(-2.0_dp, 0.3_dp, kind=dp), cmplx(1.0_dp, 0.5_dp, kind=dp), &
        cmplx(0.5_dp, -0.2_dp, kind=dp), cmplx(-1.0_dp, 0.7_dp, kind=dp)], [2, 3])

    complex(dp), allocatable :: cmat(:, :)
    complex(dp), allocatable :: cvalues(:)
    complex(dp), allocatable :: phases(:)
    real(dp), allocatable :: logabs_values(:)
    complex(dp) :: phase
    real(dp) :: logabs

    call check_complex_close('complex logsumexp vector', logsumexp(x), &
        cmplx(1000.2135700834424_dp, 1.2527113782185881_dp, kind=dp), 4.0e-13_dp)
    call logsumexp_signed(x, logabs, phase)
    call check_close('complex signed logabs', logabs, 1000.2135700834424_dp, 4.0e-13_dp)
    call check_complex_close('complex signed phase', phase, &
        cmplx(0.3127481502644492_dp, 0.9498360882311039_dp, kind=dp), 5.0e-16_dp)

    call check_complex_close('complex logsumexp real weights', logsumexp(x, real_weights), &
        cmplx(999.8536705763354_dp, 0.20097126138952082_dp, kind=dp), 4.0e-13_dp)
    call check_complex_close('complex logsumexp complex weights', logsumexp(x, complex_weights), &
        cmplx(999.5792660082635_dp, 1.820488470763454_dp, kind=dp), 4.0e-13_dp)
    call logsumexp_signed(x, logabs, phase, complex_weights)
    call check_close('complex weighted signed logabs', logabs, 999.5792660082635_dp, 4.0e-13_dp)
    call check_complex_close('complex weighted signed phase', phase, &
        cmplx(-0.24710566200237485_dp, 0.968988540596001_dp, kind=dp), 7.0e-16_dp)

    cvalues = softmax(x)
    call check_complex_array('complex softmax vector', cvalues, [ &
        cmplx(0.7820415706337492_dp, -0.2019482276580129_dp, kind=dp), &
        cmplx(0.2179584293662508_dp, 0.20194822765801285_dp, kind=dp), &
        cmplx(0.0_dp, 0.0_dp, kind=dp)], 8.0e-14_dp)
    cvalues = log_softmax(x)
    call check_complex_array('complex log_softmax vector', cvalues, [ &
        cmplx(-0.21357008344233991_dp, -0.252711378218588_dp, kind=dp), &
        cmplx(-1.21357008344234_dp, 0.747288621781412_dp, kind=dp), &
        cmplx(-2000.2135700834424_dp, -4.252711378218588_dp, kind=dp)], 5.0e-13_dp)

    cvalues = log_softmax([cmplx(1.0_dp, -0.1_dp, kind=dp), &
        cmplx(2.0_dp, 4.8_dp, kind=dp)])
    call check_complex_array('complex log_softmax shifted branch', cvalues, [ &
        cmplx(-1.1205166518631033_dp, -5.226139962459079_dp, kind=dp), &
        cmplx(-0.1205166518631033_dp, -0.3261399624590792_dp, kind=dp)], 8.0e-15_dp)

    call check_complex_close('complex logsumexp matrix all', logsumexp(matrix), &
        cmplx(4.186779388741846_dp, 0.8493515801940533_dp, kind=dp), 5.0e-15_dp)
    cvalues = logsumexp(matrix, 1)
    call check_complex_array('complex logsumexp matrix dim1', cvalues, [ &
        cmplx(2.2647475446742416_dp, -0.3171127207551251_dp, kind=dp), &
        cmplx(3.9965533220475407_dp, 0.669946402335191_dp, kind=dp), &
        cmplx(3.0159830159414067_dp, 1.9913581383906815_dp, kind=dp)], 6.0e-15_dp)
    cvalues = logsumexp(matrix, 2)
    call check_complex_array('complex logsumexp matrix dim2', cvalues, [ &
        cmplx(2.8940715463868414_dp, 1.8399157583099874_dp, kind=dp), &
        cmplx(4.058922810031681_dp, 0.5853860494108294_dp, kind=dp)], 6.0e-15_dp)

    cvalues = logsumexp(matrix, 1, matrix_real_weights)
    call check_complex_array('complex matrix real-weight values dim1', cvalues, [ &
        cmplx(1.3773997292410338_dp, -0.0425043701429223_dp, kind=dp), &
        cmplx(4.009514200225349_dp, 0.7593591012040041_dp, kind=dp), &
        cmplx(2.2743423856852196_dp, 2.018143294995201_dp, kind=dp)], 8.0e-15_dp)
    call logsumexp_signed(matrix, 1, logabs_values, phases, matrix_real_weights)
    call check_array('complex matrix real-weight logabs dim1', logabs_values, [ &
        1.3773997292410338_dp, 4.009514200225349_dp, 2.2743423856852196_dp], 6.0e-15_dp)
    call check_complex_array('complex matrix real-weight phase dim1', phases, [ &
        cmplx(0.9990968252462362_dp, -0.0424915730809248_dp, kind=dp), &
        cmplx(0.7252773907717361_dp, 0.6884567571281018_dp, kind=dp), &
        cmplx(-0.43257509134870786_dp, 0.9015979094611173_dp, kind=dp)], 8.0e-16_dp)

    cvalues = logsumexp(matrix, 2, matrix_complex_weights)
    call check_complex_array('complex matrix complex-weight values dim2', cvalues, [ &
        cmplx(2.685774140775134_dp, 1.5218208046840007_dp, kind=dp), &
        cmplx(4.0524886352729075_dp, 1.1449511975724893_dp, kind=dp)], 8.0e-15_dp)
    call logsumexp_signed(matrix, logabs, phase, matrix_complex_weights)
    call check_close('complex matrix complex-weight logabs all', logabs, &
        4.268087300252341_dp, 7.0e-15_dp)
    call check_complex_close('complex matrix complex-weight phase all', phase, &
        cmplx(0.34303567834140347_dp, 0.9393223745790651_dp, kind=dp), 8.0e-16_dp)

    cmat = softmax(matrix, 1)
    call check_complex_array('complex softmax dim1 normalization', sum(cmat, dim=1), [ &
        cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(1.0_dp, 0.0_dp, kind=dp), &
        cmplx(1.0_dp, 0.0_dp, kind=dp)], 8.0e-16_dp)
    call check_complex_close('complex softmax dim1 reference', cmat(1, 1), &
        cmplx(0.24539850674586_dp, 0.1395665576571547_dp, kind=dp), 8.0e-16_dp)
    cmat = log_softmax(matrix, 2)
    call check_complex_close('complex log_softmax dim2 reference', cmat(1, 3), &
        cmplx(0.10592845361315877_dp, 0.16008424169001248_dp, kind=dp), 8.0e-16_dp)
    call check_complex_array('complex log_softmax dim2 normalization', sum(exp(cmat), dim=2), [ &
        cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(1.0_dp, 0.0_dp, kind=dp)], 9.0e-16_dp)

    call logsumexp_signed([cmplx(0.0_dp, 0.0_dp, kind=dp), &
        cmplx(0.0_dp, 0.0_dp, kind=dp)], logabs, phase, [ &
        cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(-1.0_dp, 0.0_dp, kind=dp)])
    call check_true('complex exact cancellation logabs', logabs < -huge(1.0_dp))
    call check_complex_close('complex exact cancellation phase', phase, &
        cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex_close('complex zero weight suppresses huge term', &
        logsumexp([cmplx(1000.0_dp, 1.0_dp, kind=dp), cmplx(1.0_dp, 2.0_dp, kind=dp)], &
        [0.0_dp, 1.0_dp]), cmplx(1.0_dp, 2.0_dp, kind=dp), 2.0e-16_dp)
    call check_true('complex weight size mismatch', complex_is_nan(logsumexp(x, [1.0_dp, 2.0_dp])))
    cvalues = logsumexp(matrix, 3)
    call check_true('complex invalid dim empty', size(cvalues) == 0)

    print *, 'complex reduction tests passed'

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
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_array

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

    subroutine check_complex_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual(:) !! computed complex values
        complex(dp), intent(in) :: expected(:) !! reference complex values
        real(dp), intent(in) :: tolerance !! maximum elementwise complex absolute error

        if (size(actual) /= size(expected) .or. any_complex_nan(actual) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_array

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

end program test_special_reductions_complex
