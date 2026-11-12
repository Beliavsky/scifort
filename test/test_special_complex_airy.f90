! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_airy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_special, only : airy, dp
    implicit none

    complex(dp) :: ai
    complex(dp) :: ai_array(3)
    complex(dp) :: aip
    complex(dp) :: aip_array(3)
    complex(dp) :: bi
    complex(dp) :: bi_array(3)
    complex(dp) :: bip
    complex(dp) :: bip_array(3)
    complex(dp) :: z
    complex(dp) :: z_array(3)
    real(dp) :: ai_real
    real(dp) :: aip_real
    real(dp) :: bi_real
    real(dp) :: bip_real
    real(dp) :: nan_value
    real(dp) :: pinf
    integer :: k

    z = cmplx(0.5_dp, 0.25_dp, kind=dp)
    call airy(z, ai, aip, bi, bip)
    call check_complex('complex airy Ai ordinary', ai, &
        cmplx(0.2280095141493811_dp, -0.05653481492820948_dp, kind=dp), 5.0e-14_dp)
    call check_complex('complex airy Aip ordinary', aip, &
        cmplx(-0.22856986552474257_dp, 0.02998452304318803_dp, kind=dp), 5.0e-14_dp)
    call check_complex('complex airy Bi ordinary', bi, &
        cmplx(0.8411392509244162_dp, 0.13322430492795645_dp, kind=dp), 5.0e-14_dp)
    call check_complex('complex airy Bip ordinary', bip, &
        cmplx(0.5096651975476902_dp, 0.10343395316881238_dp, kind=dp), 5.0e-14_dp)

    z = cmplx(3.0_dp, -2.0_dp, kind=dp)
    call airy(z, ai, aip, bi, bip)
    call check_complex('complex airy Ai outer', ai, &
        cmplx(-0.009677201058610238_dp, -0.0055246891117327074_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airy Aip outer', aip, &
        cmplx(0.020990085245160245_dp, 0.005347465695574647_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airy Bi outer', bi, &
        cmplx(-7.337352267867277_dp, 1.661523038578517_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airy Bip outer', bip, &
        cmplx(-11.918089566947568_dp, 7.254625440784052_dp, kind=dp), 8.0e-12_dp)

    z = cmplx(-4.0_dp, 1.0_dp, kind=dp)
    call airy(z, ai, aip, bi, bip)
    call check_complex('complex airy Ai negative sector', ai, &
        cmplx(-0.360008730636869_dp, -1.408384507108835_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airy Aip negative sector', aip, &
        cmplx(-2.8609722702644484_dp, 0.9599498413090052_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airy Bi negative sector', bi, &
        cmplx(1.4580011545176899_dp, -0.3412741000410008_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airy Bip negative sector', bip, &
        cmplx(-0.9829525270582962_dp, -2.7544015435341587_dp, kind=dp), 8.0e-12_dp)

    z = cmplx(-4.0_dp, 6.928203230275509_dp, kind=dp)
    call airy(z, ai, aip, bi, bip)
    call check_complex('complex airy Ai Stokes sector', ai, &
        cmplx(519435.9767980277_dp, -299896.50103109365_dp, kind=dp), 3.0e-11_dp)
    call check_complex('complex airy Bi Stokes sector', bi, &
        cmplx(299896.50103116524_dp, 519435.97679798084_dp, kind=dp), 3.0e-11_dp)

    z = cmplx(5.0_dp, 5.0_dp, kind=dp)
    call airy(z, ai, aip, bi, bip)
    call check_complex('complex airy Wronskian', ai * bip - aip * bi, &
        cmplx(1.0_dp / scifort_pi, 0.0_dp, kind=dp), 5.0e-11_dp)

    call airy(conjg(z), ai_array(1), aip_array(1), bi_array(1), bip_array(1))
    call airy(z, ai, aip, bi, bip)
    call check_complex('complex airy Ai conjugate', ai_array(1), conjg(ai), 5.0e-11_dp)
    call check_complex('complex airy Aip conjugate', aip_array(1), conjg(aip), 5.0e-11_dp)
    call check_complex('complex airy Bi conjugate', bi_array(1), conjg(bi), 5.0e-11_dp)
    call check_complex('complex airy Bip conjugate', bip_array(1), conjg(bip), 5.0e-11_dp)

    call airy(2.0_dp, ai_real, aip_real, bi_real, bip_real)
    call airy(cmplx(2.0_dp, 0.0_dp, kind=dp), ai, aip, bi, bip)
    call check_complex('complex airy exact real Ai', ai, cmplx(ai_real, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('complex airy exact real Aip', aip, cmplx(aip_real, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('complex airy exact real Bi', bi, cmplx(bi_real, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('complex airy exact real Bip', bip, cmplx(bip_real, 0.0_dp, kind=dp), 0.0_dp)

    z_array = [cmplx(0.1_dp, 0.2_dp, kind=dp), &
        cmplx(-1.0_dp, 0.4_dp, kind=dp), cmplx(8.0_dp, -0.5_dp, kind=dp)]
    call airy(z_array, ai_array, aip_array, bi_array, bip_array)
    do k = 1, 3
        call check_elemental(k, z_array(k), ai_array(k), aip_array(k), bi_array(k), bip_array(k))
    end do

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    call airy(cmplx(nan_value, 0.0_dp, kind=dp), ai, aip, bi, bip)
    call check_nan_complex('complex airy nan Ai', ai)
    call check_nan_complex('complex airy nan Aip', aip)
    call check_nan_complex('complex airy nan Bi', bi)
    call check_nan_complex('complex airy nan Bip', bip)

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    call airy(cmplx(pinf, 1.0_dp, kind=dp), ai, aip, bi, bip)
    call check_nan_complex('complex airy infinite Ai', ai)
    call check_nan_complex('complex airy infinite Aip', aip)
    call check_nan_complex('complex airy infinite Bi', bi)
    call check_nan_complex('complex airy infinite Bip', bip)

    print *, 'special complex Airy tests passed'

contains

    subroutine check_elemental(k, input, ai_value, aip_value, bi_value, bip_value)
        integer, intent(in) :: k !! array position used in the check label
        complex(dp), intent(in) :: input !! complex Airy argument at this position
        complex(dp), intent(in) :: ai_value !! elemental-array Ai result
        complex(dp), intent(in) :: aip_value !! elemental-array Ai-prime result
        complex(dp), intent(in) :: bi_value !! elemental-array Bi result
        complex(dp), intent(in) :: bip_value !! elemental-array Bi-prime result
        complex(dp) :: ai_scalar
        complex(dp) :: aip_scalar
        complex(dp) :: bi_scalar
        complex(dp) :: bip_scalar
        character(len=40) :: label

        call airy(input, ai_scalar, aip_scalar, bi_scalar, bip_scalar)
        write(label, '(a,i0)') 'complex airy elemental Ai ', k
        call check_complex(trim(label), ai_value, ai_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex airy elemental Aip ', k
        call check_complex(trim(label), aip_value, aip_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex airy elemental Bi ', k
        call check_complex(trim(label), bi_value, bi_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex airy elemental Bip ', k
        call check_complex(trim(label), bip_value, bip_scalar, 0.0_dp)
    end subroutine check_elemental

    subroutine check_complex(name, actual_value, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum scaled absolute error

        if (actual_value == expected) return
        if (ieee_is_nan(real(actual_value, dp)) .or. ieee_is_nan(aimag(actual_value)) .or. &
                abs(actual_value - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual_value
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex

    subroutine check_nan_complex(name, actual_value)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual_value !! result expected to have NaN components

        if (.not. ieee_is_nan(real(actual_value, dp)) .or. .not. ieee_is_nan(aimag(actual_value))) then
            print *, 'FAIL ', name
            print *, ' actual = ', actual_value
            error stop 1
        end if
    end subroutine check_nan_complex

end program test_special_complex_airy
