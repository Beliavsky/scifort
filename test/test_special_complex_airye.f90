! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_complex_airye
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_special, only : airye, dp
    implicit none

    complex(dp) :: eai
    complex(dp) :: eai_array(3)
    complex(dp) :: eaip
    complex(dp) :: eaip_array(3)
    complex(dp) :: ebi
    complex(dp) :: ebi_array(3)
    complex(dp) :: ebip
    complex(dp) :: ebip_array(3)
    complex(dp) :: z
    complex(dp) :: z_array(3)
    real(dp) :: eai_real
    real(dp) :: eaip_real
    real(dp) :: ebi_real
    real(dp) :: ebip_real
    real(dp) :: nan_value
    real(dp) :: pinf
    integer :: k

    z = cmplx(0.5_dp, 0.25_dp, kind=dp)
    call airye(z, eai, eaip, ebi, ebip)
    call check_complex('complex airye eAi ordinary', eai, &
        cmplx(0.29034344097400183_dp, -0.018755704972517259_dp, kind=dp), 8.0e-14_dp)
    call check_complex('complex airye eAip ordinary', eaip, &
        cmplx(-0.28518659109146494_dp, -0.013728254276261051_dp, kind=dp), 8.0e-14_dp)
    call check_complex('complex airye eBi ordinary', ebi, &
        cmplx(0.67914150291153674_dp, 0.10756620211657136_dp, kind=dp), 8.0e-14_dp)
    call check_complex('complex airye eBip ordinary', ebip, &
        cmplx(0.41150711711983451_dp, 0.083513271232970906_dp, kind=dp), 8.0e-14_dp)

    z = cmplx(5.0_dp, 5.0_dp, kind=dp)
    call airye(z, eai, eaip, ebi, ebip)
    call check_complex('complex airye direct-sector eAi', eai, &
        cmplx(0.16944898843759476_dp, -0.032833146846162210_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airye direct-sector eAip', eaip, &
        cmplx(-0.45321418699590027_dp, -0.086884543712037632_dp, kind=dp), 8.0e-12_dp)
    call check_complex('complex airye direct-sector eBi', ebi, &
        cmplx(0.12984696753176367_dp, -0.32140807879690747_dp, kind=dp), 3.0e-11_dp)
    call check_complex('complex airye direct-sector eBip', ebip, &
        cmplx(0.65133108206565860_dp, -0.64625796948702074_dp, kind=dp), 3.0e-11_dp)

    z = cmplx(-4.0_dp, 6.928203230275509_dp, kind=dp)
    call airye(z, eai, eaip, ebi, ebip)
    call check_complex('complex airye Stokes eAi', eai, &
        cmplx(0.14595662018722996_dp, -0.084268093955097573_dp, kind=dp), 4.0e-11_dp)
    call check_complex('complex airye Stokes eAip', eaip, &
        cmplx(-0.40813119295813849_dp, -0.23563465411904264_dp, kind=dp), 4.0e-11_dp)
    call check_complex('complex airye Stokes eBi', ebi, &
        cmplx(0.084268093955118126_dp, 0.14595662018721670_dp, kind=dp), 5.0e-11_dp)
    call check_complex('complex airye Stokes eBip', ebip, &
        cmplx(0.23563465411911422_dp, -0.40813119295810113_dp, kind=dp), 5.0e-11_dp)

    z = cmplx(-20.0_dp, 0.0_dp, kind=dp)
    call airye(z, eai, eaip, ebi, ebip)
    call check_complex('complex airye negative-real eAi', eai, &
        cmplx(0.17606957238525048_dp, 0.010891618370520799_dp, kind=dp), 4.0e-11_dp)
    call check_complex('complex airye negative-real eAip', eaip, &
        cmplx(-0.89115941712595448_dp, -0.055126891870865689_dp, kind=dp), 4.0e-11_dp)
    call check_complex('complex airye negative-real eBi', ebi, &
        cmplx(-0.20013930932265164_dp, 0.0_dp, kind=dp), 4.0e-11_dp)
    call check_complex('complex airye negative-real eBip', ebip, &
        cmplx(-0.79142903383953511_dp, 0.0_dp, kind=dp), 4.0e-11_dp)

    call airye(cmplx(-20.0_dp, -0.0_dp, kind=dp), &
        eai_array(1), eaip_array(1), ebi_array(1), ebip_array(1))
    call check_complex('complex airye lower-cut eAi', eai_array(1), conjg(eai), 4.0e-11_dp)
    call check_complex('complex airye lower-cut eAip', eaip_array(1), conjg(eaip), 4.0e-11_dp)
    call check_complex('complex airye lower-cut eBi', ebi_array(1), conjg(ebi), 4.0e-11_dp)
    call check_complex('complex airye lower-cut eBip', ebip_array(1), conjg(ebip), 4.0e-11_dp)

    z = cmplx(-200.0_dp, 300.0_dp, kind=dp)
    call airye(z, eai, eaip, ebi, ebip)
    call check_complex('complex airye scaled outer eAi', eai, &
        cmplx(0.055536216875284185_dp, -0.033267500839148118_dp, kind=dp), 2.0e-10_dp)
    call check_complex('complex airye scaled outer eAip', eaip, &
        cmplx(-1.0544990460549861_dp, -0.63167058538333876_dp, kind=dp), 2.0e-10_dp)
    call check_complex('complex airye scaled outer eBi', ebi, &
        cmplx(0.0076634174204135708_dp, 0.064282735087594273_dp, kind=dp), 2.0e-10_dp)
    call check_complex('complex airye scaled outer eBip', ebip, &
        cmplx(1.0074922776693838_dp, -0.70422672272923825_dp, kind=dp), 2.0e-10_dp)

    call airye(20.0_dp, eai_real, eaip_real, ebi_real, ebip_real)
    call airye(cmplx(20.0_dp, 0.0_dp, kind=dp), eai, eaip, ebi, ebip)
    call check_complex('complex airye positive-real eAi', eai, &
        cmplx(eai_real, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('complex airye positive-real eAip', eaip, &
        cmplx(eaip_real, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('complex airye positive-real eBi', ebi, &
        cmplx(ebi_real, 0.0_dp, kind=dp), 0.0_dp)
    call check_complex('complex airye positive-real eBip', ebip, &
        cmplx(ebip_real, 0.0_dp, kind=dp), 0.0_dp)

    z = cmplx(50.0_dp, 10.0_dp, kind=dp)
    call airye(conjg(z), eai_array(1), eaip_array(1), ebi_array(1), ebip_array(1))
    call airye(z, eai, eaip, ebi, ebip)
    call check_complex('complex airye eAi conjugate', eai_array(1), conjg(eai), 7.0e-11_dp)
    call check_complex('complex airye eAip conjugate', eaip_array(1), conjg(eaip), 7.0e-11_dp)
    call check_complex('complex airye eBi conjugate', ebi_array(1), conjg(ebi), 7.0e-11_dp)
    call check_complex('complex airye eBip conjugate', ebip_array(1), conjg(ebip), 7.0e-11_dp)

    z_array = [cmplx(0.1_dp, 0.2_dp, kind=dp), &
        cmplx(-3.0_dp, 0.4_dp, kind=dp), cmplx(50.0_dp, -10.0_dp, kind=dp)]
    call airye(z_array, eai_array, eaip_array, ebi_array, ebip_array)
    do k = 1, 3
        call check_elemental(k, z_array(k), eai_array(k), eaip_array(k), &
            ebi_array(k), ebip_array(k))
    end do

    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    call airye(cmplx(nan_value, 0.0_dp, kind=dp), eai, eaip, ebi, ebip)
    call check_nan_complex('complex airye nan eAi', eai)
    call check_nan_complex('complex airye nan eAip', eaip)
    call check_nan_complex('complex airye nan eBi', ebi)
    call check_nan_complex('complex airye nan eBip', ebip)

    pinf = ieee_value(1.0_dp, ieee_positive_inf)
    call airye(cmplx(pinf, 1.0_dp, kind=dp), eai, eaip, ebi, ebip)
    call check_nan_complex('complex airye infinite eAi', eai)
    call check_nan_complex('complex airye infinite eAip', eaip)
    call check_nan_complex('complex airye infinite eBi', ebi)
    call check_nan_complex('complex airye infinite eBip', ebip)

    print *, 'special complex scaled Airy tests passed'

contains

    subroutine check_elemental(k, input, eai_value, eaip_value, ebi_value, ebip_value)
        integer, intent(in) :: k !! array position used in the check label
        complex(dp), intent(in) :: input !! complex scaled-Airy argument at this position
        complex(dp), intent(in) :: eai_value !! elemental-array scaled Ai result
        complex(dp), intent(in) :: eaip_value !! elemental-array scaled Ai-prime result
        complex(dp), intent(in) :: ebi_value !! elemental-array scaled Bi result
        complex(dp), intent(in) :: ebip_value !! elemental-array scaled Bi-prime result
        complex(dp) :: eai_scalar
        complex(dp) :: eaip_scalar
        complex(dp) :: ebi_scalar
        complex(dp) :: ebip_scalar
        character(len=40) :: label

        call airye(input, eai_scalar, eaip_scalar, ebi_scalar, ebip_scalar)
        write(label, '(a,i0)') 'complex airye elemental eAi ', k
        call check_complex(trim(label), eai_value, eai_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex airye elemental eAip ', k
        call check_complex(trim(label), eaip_value, eaip_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex airye elemental eBi ', k
        call check_complex(trim(label), ebi_value, ebi_scalar, 0.0_dp)
        write(label, '(a,i0)') 'complex airye elemental eBip ', k
        call check_complex(trim(label), ebip_value, ebip_scalar, 0.0_dp)
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

        if (.not. ieee_is_nan(real(actual_value, dp)) .or. &
                .not. ieee_is_nan(aimag(actual_value))) then
            print *, 'FAIL ', name
            print *, ' actual = ', actual_value
            error stop 1
        end if
    end subroutine check_nan_complex

end program test_special_complex_airye
