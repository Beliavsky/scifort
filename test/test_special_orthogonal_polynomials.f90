! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_orthogonal_polynomials
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : dp, eval_chebyc, eval_chebys, eval_chebyt, eval_chebyu, &
        eval_genlaguerre, eval_gegenbauer, eval_hermite, eval_hermitenorm, eval_jacobi, &
        eval_laguerre, eval_legendre, eval_sh_chebyt, eval_sh_chebyu, eval_sh_jacobi, &
        eval_sh_legendre
    implicit none

    complex(dp) :: actual_complex
    real(dp) :: actual_array(3)
    real(dp) :: x_array(3)

    call check_close('legendre', eval_legendre(8, 0.37_dp), &
        -0.28260834215867053_dp, 2.0e-15_dp)
    call check_close('jacobi', eval_jacobi(7, 0.2_dp, 1.3_dp, -0.45_dp), &
        -0.6928812821187839_dp, 3.0e-15_dp)
    call check_close('gegenbauer', eval_gegenbauer(9, 0.7_dp, 0.33_dp), &
        -0.061896572471990874_dp, 3.0e-15_dp)
    call check_close('laguerre', eval_laguerre(12, 2.5_dp), &
        -0.5458447451798171_dp, 3.0e-15_dp)
    call check_close('genlaguerre', eval_genlaguerre(10, 1.2_dp, 3.4_dp), &
        -2.581676422800673_dp, 6.0e-15_dp)
    call check_close('hermite', eval_hermite(9, 0.7_dp), &
        902.062697984_dp, 2.0e-12_dp)
    call check_close('hermitenorm', eval_hermitenorm(9, 0.7_dp), &
        289.926058807_dp, 8.0e-13_dp)
    call check_close('chebyt', eval_chebyt(17, 0.123_dp), &
        0.865067030840482_dp, 4.0e-15_dp)
    call check_close('chebyu', eval_chebyu(17, 0.123_dp), &
        0.8028911856420371_dp, 4.0e-15_dp)
    call check_close('chebyc', eval_chebyc(12, 1.7_dp), &
        1.8613468309609997_dp, 7.0e-15_dp)
    call check_close('chebys', eval_chebys(12, -1.3_dp), &
        -1.2823874614190005_dp, 7.0e-15_dp)
    call check_close('shifted legendre', eval_sh_legendre(8, 0.37_dp), &
        -0.17253030691607302_dp, 3.0e-15_dp)
    call check_close('shifted jacobi', eval_sh_jacobi(7, 2.5_dp, 1.2_dp, 0.37_dp), &
        4.6960261682069505e-5_dp, 2.0e-18_dp)
    call check_close('shifted chebyt', eval_sh_chebyt(17, 0.37_dp), &
        0.9710970343977504_dp, 5.0e-15_dp)
    call check_close('shifted chebyu', eval_sh_chebyu(17, 0.37_dp), &
        1.0353654244531119_dp, 6.0e-15_dp)

    actual_complex = eval_legendre(8, cmplx(0.6_dp, 0.2_dp, dp))
    call check_complex('complex legendre', actual_complex, &
        cmplx(0.57785109999999978_dp, 1.0237752000000002_dp, dp), 3.0e-15_dp)
    actual_complex = eval_jacobi(7, 0.2_dp, 1.3_dp, cmplx(-0.3_dp, 0.4_dp, dp))
    call check_complex('complex jacobi', actual_complex, &
        cmplx(2.6009591899893962_dp, 6.4487135212600717_dp, dp), 2.0e-14_dp)
    actual_complex = eval_gegenbauer(9, 0.7_dp, cmplx(0.3_dp, -0.4_dp, dp))
    call check_complex('complex gegenbauer', actual_complex, &
        cmplx(5.6272956442434070_dp, 10.442504228733856_dp, dp), 3.0e-14_dp)
    actual_complex = eval_laguerre(12, cmplx(1.2_dp, 0.8_dp, dp))
    call check_complex('complex laguerre', actual_complex, &
        cmplx(2.0738028557854647_dp, -1.9532862016106308_dp, dp), 8.0e-15_dp)
    actual_complex = eval_genlaguerre(10, 1.2_dp, cmplx(1.2_dp, -0.8_dp, dp))
    call check_complex('complex genlaguerre', actual_complex, &
        cmplx(3.4629945308169031_dp, -7.6843030177636692_dp, dp), 1.5e-14_dp)
    actual_complex = eval_chebyt(17, cmplx(0.2_dp, 0.5_dp, dp))
    call check_complex('complex chebyt', actual_complex, &
        cmplx(190.00578434971827_dp, -2011.4647930552747_dp, dp), 3.0e-12_dp)
    actual_complex = eval_chebyu(17, cmplx(0.2_dp, 0.5_dp, dp))
    call check_complex('complex chebyu', actual_complex, &
        cmplx(-9.4361832159030769_dp, -2973.4830619878228_dp, dp), 4.0e-12_dp)
    actual_complex = eval_sh_jacobi(7, 2.5_dp, 1.2_dp, cmplx(0.37_dp, 0.2_dp, dp))
    call check_complex('complex shifted jacobi', actual_complex, &
        cmplx(5.7354675934890e-4_dp, -1.6963981410779e-4_dp, dp), 3.0e-17_dp)

    x_array = [-0.5_dp, 0.0_dp, 0.5_dp]
    actual_array = eval_legendre(4, x_array)
    call check_close('elemental 1', actual_array(1), -0.2890625_dp, 2.0e-16_dp)
    call check_close('elemental 2', actual_array(2), 0.375_dp, 0.0_dp)
    call check_close('elemental 3', actual_array(3), -0.2890625_dp, 2.0e-16_dp)

    call check_true('negative degree', ieee_is_nan(eval_legendre(-1, 0.2_dp)))
    call check_true('jacobi alpha domain', ieee_is_nan(eval_jacobi(2, -1.0_dp, 0.2_dp, 0.3_dp)))
    call check_true('gegenbauer alpha domain', ieee_is_nan(eval_gegenbauer(2, -0.5_dp, 0.3_dp)))
    call check_true('genlaguerre alpha domain', ieee_is_nan(eval_genlaguerre(2, -1.0_dp, 0.3_dp)))
    call check_true('shifted jacobi q domain', ieee_is_nan(eval_sh_jacobi(2, 1.0_dp, 0.0_dp, 0.3_dp)))

    print *, 'special orthogonal-polynomial evaluator tests passed'

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

    subroutine check_complex(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! high-precision reference value
        real(dp), intent(in) :: tolerance !! maximum complex absolute error

        if (actual == expected) return
        if (ieee_is_nan(real(actual, dp)) .or. ieee_is_nan(aimag(actual)) .or. &
                abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_orthogonal_polynomials
