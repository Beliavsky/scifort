! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_group_delay_zpk
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : group_delay, group_delay_zpk, group_delay_result, &
        design_status_success, design_status_invalid_input, design_status_singular
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: w(5) = [0.13_dp,0.54_dp,1.17_dp,2.12_dp,2.90_dp]
    type(group_delay_result) :: gd, reference
    real(dp) :: zr(2), pr(2)
    complex(dp) :: z(2), p(2), zz(2), pp(2), gain

    ! A constant ZPK response has zero delay, including complex gain.
    gd = group_delay_zpk([complex(dp) ::], [complex(dp) ::], (1.0_dp, 0.7_dp), w=w)
    call check(gd%status == design_status_success, 'constant response status')
    call close(gd%delay,spread(0.0_dp,1,size(w)), 'constant response delay')

    ! Unequal root counts include a net advance or pure delay.  In particular,
    ! H(z)=z is a one-sample advance, while H(z)=1/z is one sample of delay.
    gd = group_delay_zpk([0.0_dp], [real(dp) ::], 1.0_dp, w=w)
    call check(gd%status == design_status_success, 'advance status')
    call close(gd%delay,spread(-1.0_dp,1,size(w)), 'advance')
    gd = group_delay_zpk([real(dp) ::], [0.0_dp], 1.0_dp, w=w)
    call close(gd%delay,spread(1.0_dp,1,size(w)), 'delay')

    gd = group_delay_zpk([0.3_dp], [0.8_dp], 1.0_dp, w=w)
    reference = group_delay([1.0_dp,-0.3_dp], [1.0_dp,-0.8_dp],w=w)
    call check(gd%status == design_status_success, 'real ZPK status')
    call close(gd%frequencies, w, 'real ZPK frequency grid')
    call close(gd%delay,reference%delay, 'real ZPK and transfer-function equivalence')

    ! Complex pole/zero pairs, with complex gain, match the analytic TF result.
    z = [cmplx(0.4_dp,0.25_dp,dp),cmplx(0.4_dp,-0.25_dp,dp)]
    p = [cmplx(0.74_dp,0.12_dp,dp),cmplx(0.74_dp,-0.12_dp,dp)]
    reference = group_delay([cmplx(1.0_dp,0.0_dp,dp),-sum(z),product(z)], &
        [cmplx(1.0_dp,0.0_dp,dp),-sum(p),product(p)],w=w)
    gain = cmplx(0.4_dp,-0.3_dp,dp)
    gd = group_delay(z,p,gain,w=w)
    call check(gd%status == design_status_success, 'complex ZPK generic status')
    call close(gd%delay,reference%delay, 'complex ZPK/TF equivalence')
    gd = group_delay_zpk(z,p,4.0_dp,w=w)
    call close(gd%delay,reference%delay, 'gain invariance')

    ! Real zero/complex pole and complex zero/real pole overloads.
    zr = [0.2_dp,-0.6_dp]
    pr = [0.5_dp,-0.3_dp]
    gd = group_delay_zpk(zr,p,gain,w=w)
    call check(gd%status == design_status_success, 'real-zero complex-pole gain')
    gd = group_delay_zpk(zr,p,1.0_dp,w=w)
    call check(gd%status == design_status_success, 'real-zero real-gain')
    gd = group_delay_zpk(z,pr,gain,w=w)
    call check(gd%status == design_status_success, 'complex-zero real-pole gain')
    gd = group_delay_zpk(z,pr,1.0_dp,w=w)
    call check(gd%status == design_status_success, 'complex-zero real-pole real-gain')
    gd = group_delay(zr,pr,gain,w=w)
    call check(gd%status == design_status_success, 'real-roots complex-gain generic')
    gd = group_delay(zr,pr,1.0_dp,w=w)
    call check(gd%status == design_status_success, 'real-roots real-gain generic')

    gd = group_delay_zpk(z,p,gain,worn=24,whole=.true.,fs=2000.0_dp)
    call check(gd%status == design_status_success, 'whole-circle status')
    call check(size(gd%delay) == 24, 'whole-circle grid count')
    call check(abs(gd%frequencies(24)-2000.0_dp*23.0_dp/24.0_dp) < 1.0e-10_dp, &
        'whole-circle endpoint convention')
    gd = group_delay_zpk(z,p,gain,worn=5,fs=2000.0_dp)
    call check(abs(gd%frequencies(5)-800.0_dp) < 1.0e-11_dp, 'default half-circle grid')

    ! Undefined at exact unit-circle zeros or poles, but other bins survive.
    gd = group_delay_zpk([1.0_dp], [0.4_dp], 1.0_dp,w=[0.0_dp,0.6_dp])
    call check(gd%status == design_status_singular, 'unit-circle zero status')
    call check(ieee_is_nan(gd%delay(1)), 'zero NaN')
    call check(.not. ieee_is_nan(gd%delay(2)), 'other bin from unit-circle zero')
    gd = group_delay_zpk([0.4_dp], [1.0_dp], 1.0_dp,w=[0.0_dp,0.6_dp])
    call check(gd%status == design_status_singular, 'unit-circle pole status')
    call check(ieee_is_nan(gd%delay(1)), 'pole NaN')
    gd = group_delay_zpk(z, [cmplx(cos(0.5_dp),sin(0.5_dp),dp)], &
        cmplx(1.0_dp,0.0_dp,dp),w=[0.5_dp,0.8_dp])
    call check(gd%status == design_status_singular, 'complex unit-circle pole')
    call check(ieee_is_nan(gd%delay(1)), 'complex pole NaN')
    call check(.not. ieee_is_nan(gd%delay(2)), 'complex pole other bin')

    gd = group_delay_zpk(z,p,(0.0_dp,0.0_dp))
    call invalid(gd,'zero gain')
    gd = group_delay_zpk(z,p,(1.0_dp,0.0_dp),fs=0.0_dp)
    call invalid(gd,'bad fs')
    gd = group_delay_zpk(z,p,(1.0_dp,0.0_dp),worn=0)
    call invalid(gd,'bad bin count')
    gd = group_delay_zpk(z,p,(1.0_dp,0.0_dp),w=[ieee_value(0.0_dp,ieee_quiet_nan)])
    call invalid(gd,'NaN grid')
    zz = z
    zz(1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    gd = group_delay_zpk(zz,p,(1.0_dp,0.0_dp))
    call invalid(gd,'NaN zero')
    pp = p
    pp(1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    gd = group_delay_zpk(z,pp,(1.0_dp,0.0_dp))
    call invalid(gd,'NaN pole')

    print '(a)', 'test_signal_group_delay_zpk passed'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! condition to verify
        character(*), intent(in) :: label !! check description
        if (.not. ok) then
            print '(a)', 'failure: '//label
            error stop 'group_delay_zpk check failed'
        end if
    end subroutine check

    subroutine close(actual, expected, label)
        real(dp), intent(in) :: actual(:) !! actual group delay or frequencies
        real(dp), intent(in) :: expected(:) !! reference values
        character(*), intent(in) :: label !! failing comparison
        call check(size(actual)==size(expected),label//' shape')
        call check(all(abs(actual-expected)<3.0e-11_dp),label//' values')
    end subroutine close

    subroutine invalid(actual, label)
        type(group_delay_result), intent(in) :: actual !! invalid request response
        character(*), intent(in) :: label !! check description
        call check(actual%status==design_status_invalid_input,label//' status')
        call check(.not. allocated(actual%delay),label//' no output')
    end subroutine invalid
end program test_signal_group_delay_zpk
