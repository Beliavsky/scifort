! SPDX-License-Identifier: MIT
program test_signal_zpk_conversion
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_signal, only : dp, tf2zpk, zpk2tf, zpk2sos, sos2zpk, &
        filter_zpk_result, filter_tf_result, filter_sos_result, &
        zpk_conversion_success, zpk_conversion_invalid_input, &
        freqz_zpk, sosfreqz, frequency_response
    implicit none
    type(filter_zpk_result) :: zp
    type(filter_tf_result) :: tf
    type(filter_sos_result) :: ss
    type(frequency_response) :: from_roots, from_sos
    complex(dp), parameter :: cc = cmplx(0.4_dp,0.3_dp,dp)
    complex(dp) :: zeros(4), poles(4), z1(3), p1(3)
    complex(dp) :: q, href, hsections, numerator, denominator
    complex(dp) :: one_sos(1,6)
    real(dp) :: frequencies(5), inf
    integer :: i, j

    zp = tf2zpk([3.0_dp,9.0_dp,6.0_dp],[2.0_dp,6.0_dp,4.0_dp])
    call check(zp%status == zpk_conversion_success,'tf2zpk status')
    call check(abs(zp%gain - 1.5_dp) < 1.0e-12_dp,'tf2zpk gain')
    call check(size(zp%zeros) == 2 .and. size(zp%poles) == 2,'tf2zpk degrees')
    call check(maxval(abs(zp%zeros-zp%poles)) < 1.0e-9_dp,'tf2zpk roots')
    tf = zpk2tf(zp%zeros,zp%poles,zp%gain)
    call check(tf%status == zpk_conversion_success,'zpk2tf status')
    call check(maxval(abs(tf%b-cmplx([1.5_dp,4.5_dp,3.0_dp],0.0_dp,dp))) < &
        1.0e-9_dp,'zpk2tf numerator')
    call check(maxval(abs(tf%a-cmplx([1.0_dp,3.0_dp,2.0_dp],0.0_dp,dp))) < &
        1.0e-9_dp,'zpk2tf denominator')

    zp = tf2zpk([0.0_dp,0.0_dp],[0.0_dp,2.0_dp,1.0_dp])
    call check(zp%status == zpk_conversion_success,'zero numerator')
    call check(size(zp%zeros) == 0 .and. size(zp%poles) == 1,'zero roots')
    call check(abs(zp%gain) == 0.0_dp,'zero gain')
    zp = tf2zpk([0.0_dp,3.0_dp],[1.0_dp,1.0_dp])
    call check(size(zp%zeros) == 0,'trim numerator leading zero')
    call check(abs(zp%gain-3.0_dp) < 1.0e-12_dp,'trim numerator gain')
    zp = tf2zpk([1.0_dp],[0.0_dp,0.0_dp])
    call check(zp%status == zpk_conversion_invalid_input,'zero denominator')

    zeros = [cc,conjg(cc),cmplx(-0.7_dp,0.0_dp,dp),cmplx(0.2_dp,0.0_dp,dp)]
    poles = [cmplx(0.7_dp,0.2_dp,dp),cmplx(0.7_dp,-0.2_dp,dp), &
        cmplx(-0.4_dp,0.0_dp,dp),cmplx(-0.2_dp,0.0_dp,dp)]
    tf = zpk2tf(zeros,poles,1.2_dp)
    call check(tf%status == zpk_conversion_success,'real-conjugate transfer')
    call check(maxval(abs(aimag(tf%b))) < 1.0e-13_dp,'real b')
    call check(maxval(abs(aimag(tf%a))) < 1.0e-13_dp,'real a')
    ss = zpk2sos(zeros,poles,1.2_dp)
    call check(ss%status == zpk_conversion_success,'digital SOS status')
    call check(size(ss%sos,1) == 2 .and. size(ss%sos,2) == 6,'SOS dimensions')
    call check(maxval(abs(aimag(ss%sos))) < 1.0e-12_dp,'conjugate SOS real')
    frequencies = [0.0_dp,0.17_dp,0.8_dp,1.8_dp,3.0_dp]
    from_roots = freqz_zpk(zeros,poles,1.2_dp,w=frequencies)
    from_sos = sosfreqz(ss%sos,w=frequencies)
    call check(maxval(abs(from_roots%response-from_sos%response)) < 1.0e-10_dp, &
        'matched SOS frequency response')
    ss = zpk2sos(zeros(1:3),poles(1:3),1.2_dp,pairing='keep_odd')
    call check(ss%status == zpk_conversion_success,'keep_odd SOS status')
    call check(any(abs(ss%sos(:,3)) < 1.0e-12_dp .and. &
        abs(ss%sos(:,6)) < 1.0e-12_dp),'keep_odd preserves first-order section')
    from_roots = freqz_zpk(zeros(1:3),poles(1:3),1.2_dp,w=frequencies)
    from_sos = sosfreqz(ss%sos,w=frequencies)
    call check(maxval(abs(from_roots%response-from_sos%response)) < 1.0e-10_dp, &
        'odd SOS frequency response')

    tf = zpk2tf([0.25_dp],[0.5_dp],cmplx(1.0_dp,-0.5_dp,dp))
    call check(tf%status == zpk_conversion_success,'complex gain with real roots')
    ss = zpk2sos([0.25_dp],[0.5_dp],cmplx(1.0_dp,-0.5_dp,dp))
    call check(ss%status == zpk_conversion_success,'complex gain SOS real roots')
    ss = zpk2sos([0.5_dp],[0.7_dp,0.8_dp],1.0_dp)
    call check(ss%status == zpk_conversion_success,'unequal order SOS status')
    call check(size(ss%sos,1) == 1,'one-section unequal order')
    ss = zpk2sos([real(dp) :: ],[real(dp) :: ],2.0_dp)
    call check(ss%status == zpk_conversion_success,'constant SOS')
    call check(abs(ss%sos(1,1)-2.0_dp) < 1.0e-14_dp,'constant gain')
    ss = zpk2sos(zeros,poles,1.0_dp,pairing='nonsense')
    call check(ss%status == zpk_conversion_invalid_input,'invalid pairing')
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    tf = zpk2tf([inf],[0.5_dp],1.0_dp)
    call check(tf%status == zpk_conversion_invalid_input,'nonfinite roots')
    ! Minimal pairing uses no artificial roots at the origin; a0 can be
    ! zero and such sections are not necessarily valid for sosfilt.
    z1 = [cmplx(-1.0_dp,0.0_dp,dp),cmplx(-0.5_dp,0.5_dp,dp), &
        cmplx(-0.5_dp,-0.5_dp,dp)]
    p1 = [cmplx(0.75_dp,0.0_dp,dp),cmplx(0.8_dp,0.1_dp,dp), &
        cmplx(0.8_dp,-0.1_dp,dp)]
    ss = zpk2sos(z1,p1,1.0_dp,pairing='minimal')
    call check(ss%status == zpk_conversion_success,'minimal SOS status')
    call check(size(ss%sos,1) == 2,'minimal section count')
    call check(any(abs(ss%sos(:,4)) < 1.0e-15_dp),'minimal first-order denominator')
    call check(any(abs(ss%sos(:,1)) < 1.0e-15_dp),'minimal first-order numerator')
    do i = 1, size(frequencies)
        q = exp(cmplx(0.0_dp,frequencies(i),dp))
        href = product(q-z1)/product(q-p1)
        hsections = cmplx(1.0_dp,0.0_dp,dp)
        do j = 1, size(ss%sos,1)
            numerator = (ss%sos(j,1)*q+ss%sos(j,2))*q+ss%sos(j,3)
            denominator = (ss%sos(j,4)*q+ss%sos(j,5))*q+ss%sos(j,6)
            hsections = hsections*numerator/denominator
        end do
        call check(abs(hsections-href) < 1.0e-11_dp,'minimal positive-power response')
    end do
    zp = sos2zpk(ss%sos)
    call check(zp%status == zpk_conversion_success,'minimal SOS inverse status')
    call check(size(zp%zeros) == 4 .and. size(zp%poles) == 4,'sos2zpk two roots per section')
    call check(abs(zp%gain-1.0_dp) < 1.0e-13_dp,'minimal inverse gain')
    tf = zpk2tf(zp%zeros,zp%poles,zp%gain)
    call check(tf%status == zpk_conversion_success,'SOS inverse polynomial status')
    block
        type(filter_tf_result) :: expected
        expected = zpk2tf([z1,cmplx(0.0_dp,0.0_dp,dp)], &
            [p1,cmplx(0.0_dp,0.0_dp,dp)],1.0_dp)
        call check(maxval(abs(tf%b-expected%b)) < 1.0e-10_dp,'inverse zeros')
        call check(maxval(abs(tf%a-expected%a)) < 1.0e-10_dp,'inverse poles')
    end block

    ! A proper system with an empty numerator-root list requires
    ! leading polynomial zeros, not artificial transfer-function zeros.
    ss = zpk2sos([real(dp) :: ],[0.2_dp,0.4_dp,0.6_dp],2.0_dp,pairing='minimal')
    call check(ss%status == zpk_conversion_success,'minimal unequal degrees')
    call check(size(ss%sos,1) == 2,'minimal unequal section count')
    call check(count(abs(ss%sos(:,1)) < 1.0e-15_dp) >= 1,'minimal numerator padding')
    zp = sos2zpk(ss%sos)
    call check(zp%status == zpk_conversion_success,'minimal unequal inverse')
    call check(abs(zp%gain-2.0_dp) < 1.0e-12_dp,'minimal unequal gain')
    call check(count(abs(zp%zeros) < 1.0e-12_dp) == 4,'minimal padded zeros')
    ss = zpk2sos([0.1_dp,0.2_dp],[0.3_dp],1.0_dp,pairing='minimal')
    call check(ss%status == zpk_conversion_invalid_input,'minimal improper rejection')
    ss = zpk2sos([real(dp) :: ],[real(dp) :: ],2.0_dp,pairing='minimal')
    call check(ss%status == zpk_conversion_success,'minimal constant gain')

    ! Complex-gain SOS, section normalization, and invalid inputs.
    one_sos(1,:) = cmplx([1.0_dp,-0.4_dp,0.0_dp,2.0_dp,-0.6_dp,0.0_dp],0.0_dp,dp)
    one_sos(1,1:3) = cmplx(1.0_dp,0.2_dp,dp)*one_sos(1,1:3)
    zp = sos2zpk(one_sos)
    call check(zp%status == zpk_conversion_success,'complex gain SOS inverse')
    call check(abs(zp%gain-cmplx(0.5_dp,0.1_dp,dp)) < 1.0e-13_dp,'complex SOS gain')
    call check(size(zp%zeros) == 2 .and. size(zp%poles) == 2,'complex SOS sizes')
    zp = sos2zpk(reshape([real(dp) :: 1.0_dp,2.0_dp,3.0_dp,1.0_dp,2.0_dp],[1,5]))
    call check(zp%status == zpk_conversion_invalid_input,'bad SOS shape')
    one_sos(1,4:6) = cmplx(0.0_dp,0.0_dp,dp)
    zp = sos2zpk(one_sos)
    call check(zp%status == zpk_conversion_invalid_input,'zero denominator SOS')
    one_sos(1,4) = cmplx(inf,0.0_dp,dp)
    zp = sos2zpk(one_sos)
    call check(zp%status == zpk_conversion_invalid_input,'nonfinite SOS')
    print '(a)', 'PASS: signal ZPK, TF, and digital SOS conversions (including minimal and inverse)'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: label !! message on failure
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_zpk_conversion
