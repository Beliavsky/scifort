! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_jacobi_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, &
        ieee_positive_inf, ieee_quiet_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_special, only : ellipj, ellipk, ellipkinc
    implicit none

    real(dp) :: cn
    real(dp) :: cn_values(3)
    real(dp) :: dn
    real(dp) :: dn_values(3)
    real(dp) :: k
    real(dp) :: nan_value
    real(dp) :: ph
    real(dp) :: ph_values(3)
    real(dp) :: pos_inf
    real(dp) :: sn
    real(dp) :: sn_values(3)

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    pos_inf = ieee_value(0.0_dp, ieee_positive_inf)

    call ellipj(0.5_dp, 0.5_dp, sn, cn, dn, ph)
    call check_rel('ellipj sn ordinary', sn, 0.47075047365565736_dp, 3.0e-15_dp)
    call check_rel('ellipj cn ordinary', cn, 0.8822663948904402_dp, 3.0e-15_dp)
    call check_rel('ellipj dn ordinary', dn, 0.9429724257773857_dp, 3.0e-15_dp)
    call check_rel('ellipj ph ordinary', ph, 0.4901412054142550_dp, 3.0e-15_dp)

    call ellipj(2.0_dp, 0.5_dp, sn, cn, dn, ph)
    call check_rel('ellipj sn beyond K', sn, 0.9946623253580177_dp, 3.0e-15_dp)
    call check_rel('ellipj cn beyond K', cn, -0.10318361552776205_dp, 3.0e-14_dp)
    call check_rel('ellipj dn beyond K', dn, 0.7108610477840870_dp, 3.0e-15_dp)
    call check_rel('ellipj ph beyond K', ph, 1.6741639220482394_dp, 3.0e-15_dp)

    call ellipj(10.0_dp, 0.5_dp, sn, cn, dn, ph)
    call check_rel('ellipj sn many periods', sn, 0.8588125059527789_dp, 2.0e-14_dp)
    call check_rel('ellipj cn many periods', cn, -0.5122900346669922_dp, 2.0e-14_dp)
    call check_rel('ellipj dn many periods', dn, 0.7944938890951594_dp, 2.0e-14_dp)
    call check_rel('ellipj ph many periods', ph, 8.39183082303414_dp, 2.0e-14_dp)

    call ellipj(-2.0_dp, 0.5_dp, sn, cn, dn, ph)
    call check_rel('ellipj odd sn', sn, -0.9946623253580177_dp, 3.0e-15_dp)
    call check_rel('ellipj even cn', cn, -0.10318361552776205_dp, 3.0e-14_dp)
    call check_rel('ellipj even dn', dn, 0.7108610477840870_dp, 3.0e-15_dp)
    call check_rel('ellipj odd ph', ph, -1.6741639220482394_dp, 3.0e-15_dp)

    call ellipj(1.0_dp, 0.2_dp, sn, cn, dn, ph)
    call check_rel('ellipj m=0.2 sn', sn, 0.8264684982946885_dp, 3.0e-15_dp)
    call check_rel('ellipj m=0.2 cn', cn, 0.5629829671726513_dp, 3.0e-15_dp)
    call check_rel('ellipj m=0.2 dn', dn, 0.9291877981685427_dp, 3.0e-15_dp)
    call check_rel('ellipj m=0.2 ph', ph, 0.9728056518120449_dp, 3.0e-15_dp)

    call ellipj(3.0_dp, 0.9_dp, sn, cn, dn, ph)
    call check_rel('ellipj m=0.9 sn', sn, 0.9906305999378325_dp, 3.0e-14_dp)
    call check_rel('ellipj m=0.9 cn', cn, -0.1365687170138534_dp, 3.0e-14_dp)
    call check_rel('ellipj m=0.9 dn', dn, 0.3417395397376911_dp, 3.0e-14_dp)
    call check_rel('ellipj m=0.9 ph', ph, 1.7077931710728238_dp, 3.0e-14_dp)

    call ellipj(10.0_dp, 1.0e-12_dp, sn, cn, dn, ph)
    call check_rel('ellipj small m sn', sn, -0.5440211108873678_dp, 3.0e-14_dp)
    call check_rel('ellipj small m cn', cn, -0.8390715290777504_dp, 3.0e-14_dp)
    call check_rel('ellipj small m dn', dn, 0.9999999999998520_dp, 3.0e-15_dp)
    call check_rel('ellipj small m ph', ph, 9.999999999997614_dp, 3.0e-15_dp)

    call ellipj(10.0_dp, 1.0_dp - 1.0e-12_dp, sn, cn, dn, ph)
    call check_rel('ellipj near-one sn', sn, 0.9999999958779427_dp, 3.0e-14_dp)
    call check_rel('ellipj near-one cn', cn, 9.079710609074482e-5_dp, 3.0e-11_dp)
    call check_rel('ellipj near-one dn', dn, 9.080261258534365e-5_dp, 3.0e-9_dp)
    call check_rel('ellipj near-one ph', ph, 1.5707055296886814_dp, 3.0e-14_dp)

    ! Independent 100-digit mpmath references cover regimes where SciPy's
    ! Cephes-backed near-m=1 approximation is no longer reliable.
    call ellipj(1000.0_dp, 1.0_dp - 1.0e-12_dp, sn, cn, dn, ph)
    call check_close('ellipj near-one large-u sn', sn, &
        0.9973885998408305_dp, 5.0e-14_dp)
    call check_close('ellipj near-one large-u cn', cn, &
        -0.07222174816180854_dp, 5.0e-14_dp)
    call check_close('ellipj near-one large-u dn', dn, &
        0.07222174816869540_dp, 5.0e-14_dp)
    call check_close('ellipj near-one large-u ph', ph, &
        102.17404592219985_dp, 8.0e-14_dp)

    call ellipj(1.0e6_dp, 0.5_dp, sn, cn, dn, ph)
    call check_close('ellipj high-precision large-u sn', sn, &
        0.8547923424610302_dp, 5.0e-11_dp)
    call check_close('ellipj high-precision large-u cn', cn, &
        0.5189701834113255_dp, 8.0e-11_dp)
    call check_close('ellipj high-precision large-u dn', dn, &
        0.7966586631895698_dp, 3.0e-11_dp)
    call check_close('ellipj high-precision large-u ph', ph, &
        847213.1656000549_dp, 3.0e-10_dp)

    call ellipj(1.0e6_dp, 0.9_dp, sn, cn, dn, ph)
    call check_close('ellipj high-precision large-u m=.9 sn', sn, &
        -0.5961125458548071_dp, 2.0e-10_dp)
    call check_close('ellipj high-precision large-u m=.9 cn', cn, &
        0.8029008859594692_dp, 1.0e-10_dp)
    call check_close('ellipj high-precision large-u m=.9 dn', dn, &
        0.8247331989237795_dp, 1.0e-10_dp)
    call check_close('ellipj high-precision large-u m=.9 ph', ph, &
        609286.1237719168_dp, 3.0e-10_dp)

    ! Exact endpoints use the analytic trigonometric/hyperbolic forms.
    call ellipj(2.0_dp, 0.0_dp, sn, cn, dn, ph)
    call check_close('ellipj m=0 sn', sn, sin(2.0_dp), 0.0_dp)
    call check_close('ellipj m=0 cn', cn, cos(2.0_dp), 0.0_dp)
    call check_close('ellipj m=0 dn', dn, 1.0_dp, 0.0_dp)
    call check_close('ellipj m=0 ph', ph, 2.0_dp, 0.0_dp)

    call ellipj(1000.0_dp, 1.0_dp, sn, cn, dn, ph)
    call check_close('ellipj m=1 large sn', sn, 1.0_dp, 0.0_dp)
    call check_close('ellipj m=1 large cn', cn, 0.0_dp, 0.0_dp)
    call check_close('ellipj m=1 large dn', dn, 0.0_dp, 0.0_dp)
    call check_close('ellipj m=1 large ph', ph, 0.5_dp * scifort_pi, 0.0_dp)

    ! Quarter-period and incomplete-integral consistency.
    k = ellipk(0.5_dp)
    call ellipj(k, 0.5_dp, sn, cn, dn, ph)
    call check_close('ellipj quarter sn', sn, 1.0_dp, 4.0e-15_dp)
    call check_close('ellipj quarter cn', cn, 0.0_dp, 5.0e-15_dp)
    call check_rel('ellipj quarter dn', dn, sqrt(0.5_dp), 4.0e-15_dp)
    call check_close('ellipj quarter ph', ph, 0.5_dp * scifort_pi, 5.0e-15_dp)
    call check_rel('ellipj inverse F', ellipkinc(ph, 0.5_dp), k, 5.0e-15_dp)

    call ellipj(1.0e6_dp, 0.5_dp, sn, cn, dn, ph)
    call check_true('ellipj large u finite sn', ieee_is_finite(sn))
    call check_true('ellipj large u finite cn', ieee_is_finite(cn))
    call check_true('ellipj large u finite dn', ieee_is_finite(dn))
    call check_true('ellipj large u finite ph', ieee_is_finite(ph))
    call check_close('ellipj identity sn/cn', sn * sn + cn * cn, 1.0_dp, 4.0e-15_dp)
    call check_close('ellipj identity dn', dn * dn + 0.5_dp * sn * sn, 1.0_dp, 4.0e-15_dp)

    call ellipj([0.0_dp, 0.5_dp, 2.0_dp], 0.5_dp, &
        sn_values, cn_values, dn_values, ph_values)
    call check_close('ellipj array sn0', sn_values(1), 0.0_dp, 0.0_dp)
    call check_rel('ellipj array sn1', sn_values(2), 0.47075047365565736_dp, 3.0e-15_dp)
    call check_rel('ellipj array sn2', sn_values(3), 0.9946623253580177_dp, 3.0e-15_dp)
    call check_close('ellipj array identity', &
        maxval(abs(sn_values * sn_values + cn_values * cn_values - 1.0_dp)), 0.0_dp, 4.0e-15_dp)

    call ellipj(0.0_dp, -0.1_dp, sn, cn, dn, ph)
    call check_true('ellipj negative m sn nan', ieee_is_nan(sn))
    call check_true('ellipj negative m cn nan', ieee_is_nan(cn))
    call check_true('ellipj negative m dn nan', ieee_is_nan(dn))
    call check_true('ellipj negative m ph nan', ieee_is_nan(ph))

    call ellipj(0.0_dp, 1.1_dp, sn, cn, dn, ph)
    call check_true('ellipj m>1 nan', ieee_is_nan(sn) .and. ieee_is_nan(cn) .and. &
        ieee_is_nan(dn) .and. ieee_is_nan(ph))

    call ellipj(nan_value, 0.5_dp, sn, cn, dn, ph)
    call check_true('ellipj u nan', ieee_is_nan(sn) .and. ieee_is_nan(cn) .and. &
        ieee_is_nan(dn) .and. ieee_is_nan(ph))

    call ellipj(1.0_dp, nan_value, sn, cn, dn, ph)
    call check_true('ellipj m nan', ieee_is_nan(sn) .and. ieee_is_nan(cn) .and. &
        ieee_is_nan(dn) .and. ieee_is_nan(ph))

    call ellipj(pos_inf, 0.5_dp, sn, cn, dn, ph)
    call check_true('ellipj u infinity nan', ieee_is_nan(sn) .and. ieee_is_nan(cn) .and. &
        ieee_is_nan(dn) .and. ieee_is_nan(ph))

    print *, 'Jacobi elliptic special-function tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_rel(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error

        if (actual == expected) return
        if (ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
        if (abs(actual - expected) > tolerance * abs(expected)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' relative tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_rel

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_jacobi_elliptic
