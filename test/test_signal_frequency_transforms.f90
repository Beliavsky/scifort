! SPDX-License-Identifier: MIT
! Regression checks for analog lowpass frequency transformations.
program test_signal_frequency_transforms
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_frequency_transforms, only : &
        lp2lp, lp2hp, lp2bp, lp2bs, lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, &
        analog_tf_transform_result, analog_zpk_transform_result, &
        frequency_transform_success, frequency_transform_invalid_input, frequency_transform_too_expensive
    implicit none
    type(analog_tf_transform_result) :: ba
    type(analog_zpk_transform_result) :: zpk
    complex(dp) :: zeros(1), poles(2), gain, s
    complex(dp) :: h1, h2
    integer :: i, k

    ba = lp2lp([2.0_dp,3.0_dp],[1.0_dp,4.0_dp,5.0_dp],wo=2.0_dp)
    call check_ba(ba,[4.0_dp,12.0_dp],[1.0_dp,8.0_dp,20.0_dp])
    ba = lp2hp([2.0_dp,3.0_dp],[1.0_dp,4.0_dp,5.0_dp],wo=2.0_dp)
    call check_ba(ba,[0.6_dp,0.8_dp,0.0_dp],[1.0_dp,1.6_dp,0.8_dp])
    ba = lp2bp([2.0_dp,3.0_dp],[1.0_dp,4.0_dp,5.0_dp],wo=2.0_dp,bw=0.75_dp)
    call check_ba(ba,[1.5_dp,1.6875_dp,6.0_dp,0.0_dp], &
        [1.0_dp,3.0_dp,10.8125_dp,12.0_dp,16.0_dp])
    ba = lp2bs([2.0_dp,3.0_dp],[1.0_dp,4.0_dp,5.0_dp],wo=2.0_dp,bw=0.75_dp)
    call check_ba(ba,[0.6_dp,0.3_dp,4.8_dp,1.2_dp,9.6_dp], &
        [1.0_dp,0.6_dp,8.1125_dp,2.4_dp,16.0_dp])

    zeros = [cmplx(-2.0_dp,0.0_dp,dp)]
    poles = [cmplx(-1.0_dp,2.0_dp,dp),cmplx(-1.0_dp,-2.0_dp,dp)]
    gain = cmplx(3.0_dp,0.0_dp,dp)
    zpk = lp2lp_zpk(zeros,poles,gain,2.0_dp)
    if (zpk%status /= frequency_transform_success) error stop 'lp2lp ZPK status'
    if (abs(zpk%gain-6.0_dp) > 1.0e-12_dp) error stop 'lp2lp ZPK gain'
    zpk = lp2hp_zpk(zeros,poles,gain,2.0_dp)
    if (zpk%status /= frequency_transform_success) error stop 'lp2hp ZPK status'
    if (abs(zpk%gain-1.2_dp) > 1.0e-12_dp) error stop 'lp2hp ZPK gain'
    zpk = lp2bp_zpk(zeros,poles,gain,2.0_dp,0.75_dp)
    if (zpk%status /= frequency_transform_success) error stop 'lp2bp ZPK status'
    if (size(zpk%zeros) /= 3 .or. size(zpk%poles) /= 4) error stop 'lp2bp ZPK shapes'
    if (abs(zpk%gain-2.25_dp) > 1.0e-12_dp) error stop 'lp2bp ZPK gain'
    zpk = lp2bs_zpk(zeros,poles,gain,2.0_dp,0.75_dp)
    if (zpk%status /= frequency_transform_success) error stop 'lp2bs ZPK status'
    if (size(zpk%zeros) /= 4 .or. size(zpk%poles) /= 4) error stop 'lp2bs ZPK shapes'
    if (abs(zpk%gain-1.2_dp) > 1.0e-12_dp) error stop 'lp2bs ZPK gain'

    ! Check the transformed systems by evaluating their rational functions.
    do k = 1, 4
        select case(k)
        case(1)
            zpk = lp2lp_zpk(zeros,poles,gain,2.0_dp)
        case(2)
            zpk = lp2hp_zpk(zeros,poles,gain,2.0_dp)
        case(3)
            zpk = lp2bp_zpk(zeros,poles,gain,2.0_dp,0.75_dp)
        case(4)
            zpk = lp2bs_zpk(zeros,poles,gain,2.0_dp,0.75_dp)
        end select
        do i = 1, 12
            s = cmplx(0.3_dp+0.08_dp*i,0.2_dp*i,dp)
            select case(k)
            case(1)
                h1 = zpk_eval(zeros,poles,gain,s/2.0_dp)
            case(2)
                h1 = zpk_eval(zeros,poles,gain,2.0_dp/s)
            case(3)
                h1 = zpk_eval(zeros,poles,gain,(s*s+4.0_dp)/(0.75_dp*s))
            case(4)
                h1 = zpk_eval(zeros,poles,gain,0.75_dp*s/(s*s+4.0_dp))
            end select
            h2 = zpk_eval(zpk%zeros,zpk%poles,zpk%gain,s)
            if (abs(h1-h2) > 3.0e-12_dp*max(1.0_dp,abs(h1))) error stop 'ZPK functional check'
        end do
    end do

    ba = lp2hp([0.0_dp,1.0_dp],[1.0_dp,2.0_dp,3.0_dp],2.0_dp)
    call check_ba(ba,[1.0_dp/3.0_dp,0.0_dp,0.0_dp], &
        [1.0_dp,4.0_dp/3.0_dp,4.0_dp/3.0_dp])
    ba = lp2bp([1.0_dp],[1.0_dp,1.0_dp],wo=2.0_dp,bw=0.75_dp)
    call check_ba(ba,[0.75_dp,0.0_dp],[1.0_dp,0.75_dp,4.0_dp])
    ba = lp2lp([0.0_dp],[1.0_dp,1.0_dp])
    if (ba%status /= frequency_transform_success) error stop 'zero numerator status'
    if (size(ba%b) /= 1 .or. abs(ba%b(1)) /= 0.0_dp) error stop 'zero numerator result'
    ba = lp2lp([1.0_dp],[1.0_dp,2.0_dp],wo=0.0_dp)
    if (ba%status /= frequency_transform_invalid_input) error stop 'invalid frequency'
    ba = lp2bp([1.0_dp],[1.0_dp],bw=-1.0_dp)
    if (ba%status /= frequency_transform_invalid_input) error stop 'invalid bandwidth'
    ba = lp2lp([1.0_dp],[real(dp) ::])
    if (ba%status /= frequency_transform_invalid_input) error stop 'empty denominator'
    ba = lp2lp([1.0_dp],[0.0_dp])
    if (ba%status /= frequency_transform_invalid_input) error stop 'zero denominator'
    ba = lp2lp([1.0_dp],[1.0_dp,2.0_dp],wo=ieee_value(0.0_dp,ieee_quiet_nan))
    if (ba%status /= frequency_transform_invalid_input) error stop 'nan frequency'
    ba = lp2lp([1.0_dp],[(1.0_dp, i=1,35)])
    if (ba%status /= frequency_transform_too_expensive) error stop 'large order'
    zpk = lp2hp_zpk([cmplx(0.0_dp,0.0_dp,dp)],poles,3.0_dp)
    if (zpk%status /= frequency_transform_invalid_input) error stop 'zero root highpass'
    zpk = lp2bs_zpk(zeros,[poles,cmplx(0.0_dp,0.0_dp,dp)],3.0_dp)
    if (zpk%status /= frequency_transform_invalid_input) error stop 'zero root bandstop'
    ! Non-conjugate complex prototypes also retain their exact complex gain.
    zpk = lp2hp_zpk([cmplx(-0.5_dp,0.8_dp,dp)], &
        [cmplx(-1.0_dp,0.4_dp,dp)],cmplx(2.0_dp,0.7_dp,dp),1.2_dp)
    if (zpk%status /= frequency_transform_success) error stop "complex-root status"
    do i = 1, 5
        s = cmplx(0.2_dp*i,0.7_dp,dp)
        h1 = zpk_eval([cmplx(-0.5_dp,0.8_dp,dp)], &
            [cmplx(-1.0_dp,0.4_dp,dp)],cmplx(2.0_dp,0.7_dp,dp),1.2_dp/s)
        h2 = zpk_eval(zpk%zeros,zpk%poles,zpk%gain,s)
        if (abs(h1-h2) > 1.0e-12_dp) error stop "complex-root functional check"
    end do
    zpk = lp2lp_zpk([poles,zeros],poles,3.0_dp)
    if (zpk%status /= frequency_transform_invalid_input) error stop 'improper ZPK'
    print '(a)', 'test_signal_frequency_transforms passed'
contains
    subroutine check_ba(res, b, a)
        type(analog_tf_transform_result), intent(in) :: res !! computed polynomial coefficients
        real(dp), intent(in) :: b(:) !! expected numerator
        real(dp), intent(in) :: a(:) !! expected denominator
        if (res%status /= frequency_transform_success) error stop 'BA transformation status'
        if (size(res%b) /= size(b) .or. size(res%a) /= size(a)) then
            print *, 'shapes', size(res%b),size(b),size(res%a),size(a)
            error stop 'BA transformation shape'
        end if
        if (any(abs(res%b-cmplx(b,0.0_dp,dp)) > 1.0e-11_dp)) error stop 'BA numerator coefficients'
        if (any(abs(res%a-cmplx(a,0.0_dp,dp)) > 1.0e-11_dp)) error stop 'BA denominator coefficients'
    end subroutine check_ba

    pure function zpk_eval(z, p, k, s) result(h)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: k !! gain
        complex(dp), intent(in) :: s !! evaluation point
        complex(dp) :: h
        integer :: j
        h = k
        do j = 1, size(z)
            h = h*(s-z(j))
        end do
        do j = 1, size(p)
            h = h/(s-p(j))
        end do
    end function zpk_eval
end program test_signal_frequency_transforms
