program test_special_pbwa
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_special, only : pbwa
    implicit none

    real(dp) :: a(3)
    real(dp) :: inf_value
    real(dp) :: nan_value
    real(dp) :: w
    real(dp) :: w_array(3)
    real(dp) :: w_minus
    real(dp) :: w_plus
    real(dp) :: wp
    real(dp) :: wp_array(3)
    real(dp) :: wp_minus
    real(dp) :: wp_plus
    real(dp) :: x(3)

    inf_value = ieee_value(0.0_dp, ieee_positive_inf)
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)

    ! Independent 80-digit references from DLMF 12.14.1-2 and 12.14.15-16.
    call pbwa(0.0_dp, 0.0_dp, w, wp)
    call check_close(w, 1.0227656721131687_dp, 3.0e-14_dp, 'origin value')
    call check_close(wp, -0.48887053372346190_dp, 3.0e-14_dp, 'origin derivative')

    call pbwa(0.5_dp, 1.0_dp, w, wp)
    call check_close(w, 0.47413985718436124_dp, 4.0e-14_dp, 'ordinary positive value')
    call check_close(wp, -0.29012438869130212_dp, 4.0e-14_dp, &
        'ordinary positive derivative')

    call pbwa(-0.5_dp, 2.0_dp, w, wp)
    call check_close(w, -0.64186315714489415_dp, 5.0e-14_dp, 'negative parameter value')
    call check_close(wp, -0.50836211359434181_dp, 5.0e-14_dp, &
        'negative parameter derivative')

    call pbwa(1.2_dp, -1.5_dp, w, wp)
    call check_close(w, 3.2862061535012684_dp, 7.0e-14_dp, 'negative argument value')
    call check_close(wp, -3.0990596142323124_dp, 7.0e-14_dp, &
        'negative argument derivative')

    call pbwa(3.0_dp, 4.0_dp, w, wp)
    call check_close(w, 0.0030231060194659366_dp, 8.0e-14_dp, 'near-domain-edge value')
    call check_close(wp, -0.0067993156933778040_dp, 8.0e-14_dp, &
        'near-domain-edge derivative')

    ! SciPy's documented legacy domain includes the endpoint abs(a),abs(x)=5.
    call pbwa(5.0_dp, 0.0_dp, w, wp)
    call check_close(w, 0.47347857648660443_dp, 8.0e-14_dp, 'parameter boundary value')
    call check_close(wp, -1.0560139884473652_dp, 8.0e-14_dp, &
        'parameter boundary derivative')

    call pbwa(0.0_dp, 5.0_dp, w, wp)
    call check_close(w, 0.30095490623907518_dp, 2.0e-13_dp, 'argument boundary value')
    call check_close(wp, -0.71482331488183804_dp, 2.0e-13_dp, &
        'argument boundary derivative')

    call pbwa(0.0_dp, -5.0_dp, w, wp)
    call check_close(w, 0.65996337123562829_dp, 2.0e-13_dp, 'negative boundary value')
    call check_close(wp, -1.7552224081475017_dp, 2.0e-13_dp, &
        'negative boundary derivative')

    ! DLMF 12.14.3 Wronskian for W(a,x) and W(a,-x).
    call pbwa(2.25_dp, 2.75_dp, w_plus, wp_plus)
    call pbwa(2.25_dp, -2.75_dp, w_minus, wp_minus)
    call check_close(-(w_plus * wp_minus + w_minus * wp_plus), 1.0_dp, 2.0e-13_dp, &
        'Wronskian')

    ! Elemental array agreement.
    a = [-1.0_dp, 0.0_dp, 1.0_dp]
    x = [-1.0_dp, 0.5_dp, 2.0_dp]
    call pbwa(a, x, w_array, wp_array)
    call pbwa(a(1), x(1), w, wp)
    call check_close(w_array(1), w, 0.0_dp, 'elemental value 1')
    call check_close(wp_array(1), wp, 0.0_dp, 'elemental derivative 1')
    call pbwa(a(2), x(2), w, wp)
    call check_close(w_array(2), w, 0.0_dp, 'elemental value 2')
    call check_close(wp_array(2), wp, 0.0_dp, 'elemental derivative 2')
    call pbwa(a(3), x(3), w, wp)
    call check_close(w_array(3), w, 0.0_dp, 'elemental value 3')
    call check_close(wp_array(3), wp, 0.0_dp, 'elemental derivative 3')

    call pbwa(5.01_dp, 0.0_dp, w, wp)
    call require(ieee_is_nan(w) .and. ieee_is_nan(wp), 'parameter outside legacy domain')
    call pbwa(0.0_dp, -5.01_dp, w, wp)
    call require(ieee_is_nan(w) .and. ieee_is_nan(wp), 'argument outside legacy domain')
    call pbwa(nan_value, 0.0_dp, w, wp)
    call require(ieee_is_nan(w) .and. ieee_is_nan(wp), 'NaN parameter')
    call pbwa(inf_value, 0.0_dp, w, wp)
    call require(ieee_is_nan(w) .and. ieee_is_nan(wp), 'infinite parameter')
    call pbwa(0.0_dp, -inf_value, w, wp)
    call require(ieee_is_nan(w) .and. ieee_is_nan(wp), 'infinite argument')

    print *, 'test_special_pbwa: PASS'

contains

    subroutine check_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! independent reference value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

end program test_special_pbwa
