! SPDX-License-Identifier: MIT
program test_signal_state_space
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, tf2ss, ss2tf, zpk2ss, ss2zpk, &
        state_space_result, state_space_tf_result, filter_zpk_result, &
        state_space_success, state_space_invalid_input, zpk_conversion_success
    implicit none
    type(state_space_result) :: ss
    type(state_space_tf_result) :: tf
    type(filter_zpk_result) :: zp
    real(dp), parameter :: tol = 3.0e-9_dp
    real(dp) :: nanvalue
    real(dp) :: am(2,2), bm(2,2), cm(2,2), dm(2,2)
    real(dp) :: empty_a(0,0), empty_b(0,1), empty_c(1,0), static_d(1,1)
    complex(dp) :: num(3), den(3)

    ! Controller canonical form, as in scipy.signal.tf2ss documentation.
    ss = tf2ss([1.0_dp,3.0_dp,3.0_dp],[1.0_dp,2.0_dp,1.0_dp])
    call check(ss%status == state_space_success,'tf2ss status')
    call check(maxval(abs(ss%a-reshape([(-2.0_dp,0.0_dp),(1.0_dp,0.0_dp), &
        (-1.0_dp,0.0_dp),(0.0_dp,0.0_dp)],[2,2]))) < tol,'controller A')
    call check(maxval(abs(ss%b(:,1)-[cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(0.0_dp,0.0_dp,dp)])) < tol,'controller B')
    call check(maxval(abs(ss%c(1,:)-cmplx([1.0_dp,2.0_dp],0.0_dp,dp))) < tol,'controller C')
    call check(abs(ss%d(1,1)-1.0_dp) < tol,'controller D')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(tf%status == state_space_success,'ss2tf status')
    call check(maxval(abs(tf%b(1,:)-cmplx([1.0_dp,3.0_dp,3.0_dp],0.0_dp,dp))) < tol, &
        'ss2tf numerator')
    call check(maxval(abs(tf%a-cmplx([1.0_dp,2.0_dp,1.0_dp],0.0_dp,dp))) < tol, &
        'ss2tf denominator')

    ! Normalization and numerator padding.
    ss = tf2ss([0.0_dp,2.0_dp],[0.0_dp,2.0_dp,4.0_dp,2.0_dp])
    call check(ss%status == state_space_success,'leading zero normalization')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(maxval(abs(tf%b(1,:)-cmplx([0.0_dp,0.0_dp,1.0_dp],0.0_dp,dp))) < tol, &
        'normalized numerator')
    call check(maxval(abs(tf%a-cmplx([1.0_dp,2.0_dp,1.0_dp],0.0_dp,dp))) < tol, &
        'normalized denominator')

    ss = tf2ss([5.0_dp],[2.0_dp])
    call check(ss%status == state_space_success,'static gain')
    call check(size(ss%a,1)==1 .and. abs(ss%d(1,1)-2.5_dp)<tol,'static gain shape')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(maxval(abs(tf%b(1,:)-cmplx([2.5_dp,0.0_dp],0.0_dp,dp))) < tol, &
        'static gain numerator')
    ss = tf2ss([0.0_dp],[1.0_dp,2.0_dp])
    call check(ss%status == state_space_success,'zero numerator')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(maxval(abs(tf%b))==0.0_dp,'zero numerator roundtrip')

    ! Multi-output transfer function (one input).
    ss = tf2ss(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], &
        [2,3]),[1.0_dp,2.0_dp,3.0_dp])
    call check(ss%status == state_space_success,'SIMO tf2ss')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(tf%status == state_space_success,'SIMO ss2tf')
    call check(maxval(abs(tf%b-cmplx(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
        5.0_dp,6.0_dp],[2,3]),0.0_dp,dp))) < tol,'SIMO coefficients')

    ! MIMO selected input (one-based, unlike Python).
    am = reshape([-1.0_dp,0.0_dp,0.0_dp,-2.0_dp],[2,2])
    bm = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
    cm = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
    dm = reshape([0.0_dp,0.0_dp,0.0_dp,0.5_dp],[2,2])
    tf = ss2tf(am,bm,cm,dm,input=2)
    call check(tf%status == state_space_success,'MIMO input selection')
    call check(maxval(abs(tf%a-cmplx([1.0_dp,3.0_dp,2.0_dp],0.0_dp,dp))) < tol, &
        'MIMO denominator')
    call check(maxval(abs(tf%b(1,:))) < tol,'MIMO disconnected first output')
    call check(maxval(abs(tf%b(2,:)-cmplx([0.5_dp,2.5_dp,2.0_dp],0.0_dp,dp))) < tol, &
        'MIMO selected numerator')
    zp = ss2zpk(am,bm,cm,dm,input=2,output=2)
    call check(zp%status == zpk_conversion_success,'MIMO selected ZPK')
    call check(size(zp%poles)==2,'MIMO ZPK number of poles')
    zp = ss2zpk(am,bm,cm,dm,input=3)
    call check(zp%status /= zpk_conversion_success,'invalid MIMO input')

    ss = zpk2ss([-0.5_dp],[-1.0_dp,-2.0_dp],2.0_dp)
    call check(ss%status == state_space_success,'zpk2ss real')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(maxval(abs(tf%b(1,:)-cmplx([0.0_dp,2.0_dp,1.0_dp],0.0_dp,dp))) < tol, &
        'ZPK numerator reconstruction')
    zp = ss2zpk(ss%a,ss%b,ss%c,ss%d)
    call check(zp%status == zpk_conversion_success,'ss2zpk')
    call check(abs(zp%gain-2.0_dp) < tol,'ss2zpk gain')
    call check(minval(abs(zp%poles+1.0_dp)) < tol .and. &
        minval(abs(zp%poles+2.0_dp)) < tol,'ss2zpk poles')
    call check(abs(zp%zeros(1)+0.5_dp) < tol,'ss2zpk zeros')

    num = [cmplx(1.0_dp,0.0_dp,dp),cmplx(0.5_dp,-0.2_dp,dp), &
        cmplx(0.1_dp,0.3_dp,dp)]
    den = [cmplx(1.0_dp,0.0_dp,dp),cmplx(1.2_dp,0.25_dp,dp), &
        cmplx(0.5_dp,-0.1_dp,dp)]
    ss = tf2ss(num,den)
    call check(ss%status==state_space_success,'complex tf2ss')
    tf = ss2tf(ss%a,ss%b,ss%c,ss%d)
    call check(maxval(abs(tf%b(1,:)-num))<tol,'complex numerator')
    call check(maxval(abs(tf%a-den))<tol,'complex denominator')

    ! Explicit zero-dimensional state-space realization.
    static_d(1,1) = 3.0_dp
    tf = ss2tf(empty_a,empty_b,empty_c,static_d)
    call check(tf%status==state_space_success,'zero-dimensional state')
    call check(size(tf%a)==1 .and. size(tf%b,2)==1,'zero-state coefficient shape')
    call check(abs(tf%a(1)-1.0_dp)<tol .and. abs(tf%b(1,1)-3.0_dp)<tol, &
        'zero-state feedthrough')

    ss = zpk2ss([0.0_dp,1.0_dp],[2.0_dp],1.0_dp)
    call check(ss%status==state_space_invalid_input,'improper ZPK rejection')

    ss = tf2ss([1.0_dp,2.0_dp,3.0_dp],[1.0_dp,2.0_dp])
    call check(ss%status==state_space_invalid_input,'improper rejection')
    ss = tf2ss([1.0_dp],[0.0_dp,0.0_dp])
    call check(ss%status==state_space_invalid_input,'zero denominator rejection')
    ss = tf2ss([real(dp) :: ],[1.0_dp])
    call check(ss%status==state_space_invalid_input,'empty numerator rejection')
    nanvalue = ieee_value(0.0_dp,ieee_quiet_nan)
    ss = tf2ss([nanvalue],[1.0_dp])
    call check(ss%status==state_space_invalid_input,'NaN rejection')
    tf = ss2tf(am,bm,cm,dm,input=0)
    call check(tf%status==state_space_invalid_input,'invalid input selection')
    tf = ss2tf(am,bm(1:1,:),cm,dm)
    call check(tf%status==state_space_invalid_input,'shape mismatch rejection')

    print '(a)', 'PASS: state-space conversions'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: label !! diagnostic name
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_state_space
