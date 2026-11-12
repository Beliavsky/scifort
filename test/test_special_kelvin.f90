! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_kelvin
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_special, only : bei, bei_zeros, beip, beip_zeros, ber, ber_zeros, berp, &
        berp_zeros, dp, kei, kei_zeros, keip, keip_zeros, kelvin, kelvin_max_zeros, &
        kelvin_status_invalid_argument, kelvin_status_success, kelvin_zeros, ker, ker_zeros, &
        kerp, kerp_zeros
    implicit none

    complex(dp) :: be
    complex(dp) :: bep
    complex(dp) :: ke
    complex(dp) :: kep
    complex(dp) :: be_array(3)
    complex(dp) :: bep_array(3)
    complex(dp) :: ke_array(3)
    complex(dp) :: kep_array(3)
    real(dp), allocatable :: z1(:)
    real(dp), allocatable :: z2(:)
    real(dp), allocatable :: z3(:)
    real(dp), allocatable :: z4(:)
    real(dp), allocatable :: z5(:)
    real(dp), allocatable :: z6(:)
    real(dp), allocatable :: z7(:)
    real(dp), allocatable :: z8(:)
    real(dp) :: scalar_array(3)
    integer :: status

    call check_close(ber(0.1_dp), 0.999998437500067816839806827499_dp, 2.0e-14_dp, 'ber 0.1')
    call check_close(bei(0.1_dp), 0.002499999565972229003906225972_dp, 2.0e-14_dp, 'bei 0.1')
    call check_close(ker(0.1_dp), 2.42047398103817142116927510343_dp, 3.0e-14_dp, 'ker 0.1')
    call check_close(kei(0.1_dp), -0.776850646536661068870242896722_dp, 3.0e-14_dp, 'kei 0.1')
    call check_close(berp(0.1_dp), -0.0000624999945746528342918111924_dp, &
        2.0e-14_dp, 'berp 0.1')
    call check_close(beip(0.1_dp), 0.049999973958334011501732747181_dp, 2.0e-14_dp, 'beip 0.1')
    call check_close(kerp(0.1_dp), -9.96095939454749907192444700719_dp, 3.0e-14_dp, 'kerp 0.1')
    call check_close(keip(0.1_dp), 0.145974811428596182188045831651_dp, 3.0e-14_dp, 'keip 0.1')

    call check_close(ber(10.0_dp), 138.840465941632647206379898774_dp, 3.0e-14_dp, 'ber 10')
    call check_close(bei(10.0_dp), 56.3704585539066382327749170916_dp, 3.0e-14_dp, 'bei 10')
    call check_close(ker(10.0_dp), 0.000129466330214806122154239374410_dp, &
        2.0e-13_dp, 'ker 10')
    call check_close(kei(10.0_dp), -0.000307524569088144199016772703501_dp, &
        2.0e-13_dp, 'kei 10')
    call check_close(berp(10.0_dp), 51.1952583935996113822673273510_dp, 3.0e-14_dp, 'berp 10')
    call check_close(beip(10.0_dp), 135.309301715964613584056848162_dp, 3.0e-14_dp, 'beip 10')
    call check_close(kerp(10.0_dp), -0.000315596934420897502199359967003_dp, &
        2.0e-13_dp, 'kerp 10')
    call check_close(keip(10.0_dp), 0.000140913837612616130626223640454_dp, &
        2.0e-13_dp, 'keip 10')

    call check_close(ber(25.0_dp), 9797.71694973549028112928403612_dp, 1.0e-12_dp, 'ber 25')
    call check_close(bei(25.0_dp), -3808789.91144403568678926720351_dp, 5.0e-13_dp, 'bei 25')
    call check_close(ker(25.0_dp), 3.72329131364329485162252098589e-9_dp, 3.0e-14_dp, 'ker 25')
    call check_close(kei(25.0_dp), 3.70270353525262407342589523518e-9_dp, 3.0e-14_dp, 'kei 25')

    call kelvin(1.0_dp, be, ke, bep, kep)
    call check_complex(be, cmplx(0.984381781213086883965555659140_dp, &
        0.249566040036659721419359848310_dp, kind=dp), 3.0e-14_dp, 'kelvin Be')
    call check_complex(ke, cmplx(0.286706208728316045954197314872_dp, &
        -0.494994636518719900347050652586_dp, kind=dp), 3.0e-14_dp, 'kelvin Ke')
    call check_complex(bep, cmplx(-0.062445752179030960235859443122_dp, &
        0.497396511468097326857087607634_dp, kind=dp), 3.0e-14_dp, 'kelvin Bep')
    call check_complex(kep, cmplx(-0.694603891100690521230750673513_dp, &
        0.352369913336170534389139824763_dp, kind=dp), 3.0e-14_dp, 'kelvin Kep')

    call require(ber(-2.0_dp) == ber(2.0_dp), 'ber even parity')
    call require(bei(-2.0_dp) == bei(2.0_dp), 'bei even parity')
    call require(berp(-2.0_dp) == -berp(2.0_dp), 'berp odd parity')
    call require(beip(-2.0_dp) == -beip(2.0_dp), 'beip odd parity')
    call require(ieee_is_nan(ker(-1.0_dp)) .and. ieee_is_nan(kei(-1.0_dp)), &
        'negative K-family domain')
    call require(.not. ieee_is_finite(ker(0.0_dp)) .and. ker(0.0_dp) > 0.0_dp, 'ker zero')
    call check_close(kei(0.0_dp), -0.25_dp * scifort_pi, 0.0_dp, 'kei zero')
    call require(.not. ieee_is_finite(kerp(0.0_dp)) .and. kerp(0.0_dp) < 0.0_dp, 'kerp zero')
    call require(keip(0.0_dp) == 0.0_dp, 'keip zero')

    scalar_array = ber([0.1_dp, 1.0_dp, 5.0_dp])
    call check_close(scalar_array(1), 0.999998437500067816839806827499_dp, &
        2.0e-14_dp, 'ber elemental 1')
    call check_close(scalar_array(2), 0.984381781213086883965555659140_dp, &
        2.0e-14_dp, 'ber elemental 2')
    call check_close(scalar_array(3), -6.23008247866635773318507867901_dp, &
        3.0e-14_dp, 'ber elemental 3')
    call kelvin([0.1_dp, 1.0_dp, 5.0_dp], be_array, ke_array, bep_array, kep_array)
    call check_close(real(be_array(2), dp), ber(1.0_dp), 0.0_dp, 'kelvin elemental Be')
    call check_close(aimag(ke_array(2)), kei(1.0_dp), 0.0_dp, 'kelvin elemental Ke')
    call check_close(real(bep_array(3), dp), berp(5.0_dp), 0.0_dp, 'kelvin elemental Bep')
    call check_close(aimag(kep_array(1)), keip(0.1_dp), 0.0_dp, 'kelvin elemental Kep')

    call ber_zeros(4, z1, status)
    call require(status == kelvin_status_success, 'ber zeros status')
    call check_close(z1(1), 2.8489178207951396_dp, 3.0e-14_dp, 'ber zero 1')
    call check_close(z1(4), 16.113563827387901_dp, 3.0e-14_dp, 'ber zero 4')
    call bei_zeros(4, z2, status)
    call check_close(z2(1), 5.0262239519531510_dp, 3.0e-14_dp, 'bei zero 1')
    call check_close(z2(4), 18.333983455774560_dp, 3.0e-14_dp, 'bei zero 4')
    call ker_zeros(4, z3, status)
    call check_close(z3(1), 1.7185429596232309_dp, 3.0e-14_dp, 'ker zero 1')
    call check_close(z3(4), 15.002688121534660_dp, 3.0e-14_dp, 'ker zero 4')
    call kei_zeros(4, z4, status)
    call check_close(z4(1), 3.9146676068432563_dp, 3.0e-14_dp, 'kei zero 1')
    call check_close(z4(4), 17.223143723434830_dp, 3.0e-14_dp, 'kei zero 4')
    call berp_zeros(4, z5, status)
    call check_close(z5(1), 6.0387108067212780_dp, 3.0e-14_dp, 'berp zero 1')
    call check_close(z5(4), 19.417574926000770_dp, 3.0e-14_dp, 'berp zero 4')
    call beip_zeros(4, z6, status)
    call check_close(z6(1), 3.7726733049349530_dp, 3.0e-14_dp, 'beip zero 1')
    call check_close(z6(4), 17.193431752512570_dp, 3.0e-14_dp, 'beip zero 4')
    call kerp_zeros(4, z7, status)
    call check_close(z7(1), 2.6658397930175576_dp, 3.0e-14_dp, 'kerp zero 1')
    call check_close(z7(4), 16.083120249405760_dp, 3.0e-14_dp, 'kerp zero 4')
    call keip_zeros(4, z8, status)
    call check_close(z8(1), 4.9318119411525510_dp, 3.0e-14_dp, 'keip zero 1')
    call check_close(z8(4), 18.307172935599080_dp, 3.0e-14_dp, 'keip zero 4')

    call kelvin_zeros(3, z1, z2, z3, z4, z5, z6, z7, z8, status)
    call require(status == kelvin_status_success, 'kelvin zeros status')
    call check_close(z1(3), 11.673963549642490_dp, 4.0e-14_dp, 'kelvin zeros ber')
    call check_close(z8(3), 13.858269159613410_dp, 4.0e-14_dp, 'kelvin zeros keip')

    call ber_zeros(kelvin_max_zeros, z1, status)
    call require(status == kelvin_status_success, 'maximum zero table status')
    call check_close(z1(kelvin_max_zeros), 886.9106475809855_dp, 2.0e-14_dp, &
        'ber maximum zero')

    call ber_zeros(0, z1, status)
    call require(status == kelvin_status_invalid_argument .and. size(z1) == 0, &
        'zero-count invalid status')
    call ber_zeros(kelvin_max_zeros + 1, z1, status)
    call require(status == kelvin_status_invalid_argument .and. size(z1) == 0, &
        'maximum-count invalid status')

    print *, 'test_special_kelvin: PASS'

contains

    subroutine check_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! independent reference value
        real(dp), intent(in) :: relative_tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > relative_tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex(actual, expected, relative_tolerance, label)
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! independent complex reference value
        real(dp), intent(in) :: relative_tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > relative_tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_complex

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_special_kelvin
