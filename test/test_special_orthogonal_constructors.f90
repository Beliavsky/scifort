! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_orthogonal_constructors
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : chebyc, chebys, chebyt, chebyu, dp, eval_genlaguerre, &
        eval_jacobi, eval_sh_chebyu, gegenbauer, genlaguerre, hermite, hermitenorm, &
        jacobi, laguerre, legendre, sh_chebyt, sh_chebyu, sh_jacobi, sh_legendre
    implicit none

    real(dp), allocatable :: coefficients(:)

    call check_array('legendre', legendre(4), &
        [4.375_dp, 0.0_dp, -3.75_dp, 0.0_dp, 0.375_dp], 2.0e-15_dp)
    call check_array('jacobi', jacobi(4, 0.2_dp, 0.4_dp), &
        [6.2909_dp, -0.5852_dp, -4.9434_dp, 0.2772_dp, 0.4389_dp], 2.0e-14_dp)
    call check_array('gegenbauer', gegenbauer(4, 0.7_dp), &
        [7.9254_dp, 0.0_dp, -6.426_dp, 0.0_dp, 0.595_dp], 2.0e-14_dp)
    call check_array('laguerre', laguerre(4), &
        [1.0_dp / 24.0_dp, -2.0_dp / 3.0_dp, 3.0_dp, -4.0_dp, 1.0_dp], 2.0e-15_dp)
    call check_array('genlaguerre', genlaguerre(4, 0.7_dp), &
        [1.0_dp / 24.0_dp, -0.7833333333333333_dp, 4.3475_dp, &
        -7.8255_dp, 3.3258375_dp], 2.0e-14_dp)
    call check_array('hermite', hermite(4), &
        [16.0_dp, 0.0_dp, -48.0_dp, 0.0_dp, 12.0_dp], 0.0_dp)
    call check_array('hermitenorm', hermitenorm(4), &
        [1.0_dp, 0.0_dp, -6.0_dp, 0.0_dp, 3.0_dp], 0.0_dp)
    call check_array('chebyt', chebyt(4), &
        [8.0_dp, 0.0_dp, -8.0_dp, 0.0_dp, 1.0_dp], 0.0_dp)
    call check_array('chebyu', chebyu(4), &
        [16.0_dp, 0.0_dp, -12.0_dp, 0.0_dp, 1.0_dp], 0.0_dp)
    call check_array('chebyc', chebyc(4), &
        [1.0_dp, 0.0_dp, -4.0_dp, 0.0_dp, 2.0_dp], 0.0_dp)
    call check_array('chebys', chebys(4), &
        [1.0_dp, 0.0_dp, -3.0_dp, 0.0_dp, 1.0_dp], 0.0_dp)
    call check_array('shifted legendre', sh_legendre(4), &
        [70.0_dp, -140.0_dp, 90.0_dp, -20.0_dp, 1.0_dp], 2.0e-14_dp)
    call check_array('shifted jacobi', sh_jacobi(4, 0.2_dp, 0.4_dp), &
        [1.0_dp, -1.8888888888888889_dp, 1.0967741935483871_dp, &
        -0.19685690653432575_dp, 0.004687069203198217_dp], 3.0e-15_dp)
    call check_array('shifted chebyt', sh_chebyt(4), &
        [128.0_dp, -256.0_dp, 160.0_dp, -32.0_dp, 1.0_dp], 0.0_dp)
    call check_array('shifted chebyu', sh_chebyu(4), &
        [256.0_dp, -512.0_dp, 336.0_dp, -80.0_dp, 5.0_dp], 0.0_dp)

    call check_array('monic legendre', legendre(4, .true.), &
        [1.0_dp, 0.0_dp, -6.0_dp / 7.0_dp, 0.0_dp, 3.0_dp / 35.0_dp], 2.0e-16_dp)
    call check_array('monic laguerre', laguerre(4, .true.), &
        [1.0_dp, -16.0_dp, 72.0_dp, -96.0_dp, 24.0_dp], 2.0e-14_dp)
    call check_array('monic chebyc zero', chebyc(0, .true.), [1.0_dp], 0.0_dp)
    call check_array('gegenbauer alpha zero monic', gegenbauer(6, 0.0_dp, .true.), &
        [1.0_dp, 0.0_dp, -1.5_dp, 0.0_dp, 0.5625_dp, 0.0_dp, -0.03125_dp], 0.0_dp)

    coefficients = jacobi(10, 0.2_dp, 1.3_dp)
    call check_close('jacobi coefficient evaluation', horner(coefficients, -0.45_dp), &
        eval_jacobi(10, 0.2_dp, 1.3_dp, -0.45_dp), 3.0e-13_dp)
    coefficients = genlaguerre(10, 1.2_dp)
    call check_close('genlaguerre coefficient evaluation', horner(coefficients, 0.4_dp), &
        eval_genlaguerre(10, 1.2_dp, 0.4_dp), 3.0e-13_dp)
    coefficients = sh_chebyu(10)
    call check_close('shifted chebyu coefficient evaluation', horner(coefficients, 0.37_dp), &
        eval_sh_chebyu(10, 0.37_dp), 2.0e-11_dp)

    call check_true('negative degree empty', size(legendre(-1)) == 0)
    call check_true('jacobi alpha invalid', size(jacobi(3, -1.0_dp, 0.2_dp)) == 0)
    call check_true('jacobi beta invalid', size(jacobi(3, 0.2_dp, -1.0_dp)) == 0)
    call check_true('gegenbauer invalid', size(gegenbauer(3, -0.5_dp)) == 0)
    call check_true('genlaguerre invalid', size(genlaguerre(3, -1.0_dp)) == 0)
    call check_true('shifted jacobi q invalid', size(sh_jacobi(3, 1.0_dp, 0.0_dp)) == 0)
    call check_true('shifted jacobi p invalid', size(sh_jacobi(3, -0.8_dp, 0.3_dp)) == 0)

    print *, 'special orthogonal-polynomial constructor tests passed'

contains

    pure function horner(coefficients, x) result(value)
        real(dp), intent(in) :: coefficients(:) !! descending polynomial coefficients
        real(dp), intent(in) :: x !! evaluation argument
        real(dp) :: value

        integer :: j

        if (size(coefficients) == 0) then
            value = 0.0_dp
            return
        end if
        value = coefficients(1)
        do j = 2, size(coefficients)
            value = value * x + coefficients(j)
        end do
    end function horner

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed coefficient array
        real(dp), intent(in) :: expected(:) !! reference coefficient array
        real(dp), intent(in) :: tolerance !! maximum absolute coefficient error

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', name, ' size ', size(actual), size(expected)
            error stop 1
        end if
        if (any(ieee_is_nan(actual)) .or. maxval(abs(actual - expected)) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_array

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

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_orthogonal_constructors
