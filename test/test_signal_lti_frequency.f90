! SPDX-License-Identifier: MIT
! Tests for analog/digital frequency responses and Bode magnitude/phase.
program test_signal_lti_frequency
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, freqresp, bode, dfreqresp, dbode, &
        freqresp_zpk, bode_zpk, dfreqresp_zpk, dbode_zpk, &
        frequency_response, bode_result, design_status_success, design_status_invalid_input
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tol = 2.0e-10_dp
    real(dp) :: w(4), aa(1,1), bb(1,1), cc(1,1), dd(1,1)
    real(dp) :: dm(1,2)
    complex(dp) :: ac(1,1), bc(1,1), dc(1,1), ccx(1,1)
    type(frequency_response) :: resp, other
    type(bode_result) :: plot
    integer :: i

    w = [0.0_dp,0.5_dp,1.0_dp,1.5_dp]
    resp = freqresp([1.0_dp],[1.0_dp,1.0_dp],w=w)
    call check(resp%status == design_status_success,'analog TF status')
    do i=1,size(w)
        call check(abs(resp%response(i)-1.0_dp/cmplx(1.0_dp,w(i),dp)) < tol,'analog TF values')
    end do
    plot = bode([1.0_dp],[1.0_dp,1.0_dp],w=w)
    call check(plot%status == design_status_success,'analog Bode status')
    call check(maxval(abs(plot%frequencies-w)) < tol,'analog frequency units')
    call check(maxval(abs(plot%magnitude+10.0_dp*log10(1.0_dp+w*w))) < tol,'analog dB')
    call check(maxval(abs(plot%phase+atan(w)*180.0_dp/pi)) < tol,'analog phase')
    resp = freqresp_zpk([real(dp) ::],[-1.0_dp],1.0_dp,w=w)
    call check(maxval(abs(resp%response-1.0_dp/cmplx(1.0_dp,w,dp))) < tol,'analog ZPK')
    plot = bode_zpk([real(dp) ::],[-1.0_dp],1.0_dp,w=w)
    call check(abs(plot%phase(4)+atan(w(4))*180.0_dp/pi) < tol,'ZPK Bode')

    aa = -1.0_dp
    bb = 1.0_dp
    cc = 1.0_dp
    dd = 0.0_dp
    other = freqresp(aa,bb,cc,dd,w=w)
    call check(other%status == design_status_success,'state space analog status')
    call check(maxval(abs(other%response-resp%response)) < tol,'state space analog agreement')
    plot = bode(aa,bb,cc,dd,w=w)
    call check(abs(plot%magnitude(3)+10.0_dp*log10(2.0_dp)) < tol,'state space Bode')
    ac = cmplx(-1.0_dp,0.3_dp,dp)
    bc = (1.0_dp,0.0_dp)
    ccx = (1.0_dp,0.0_dp)
    dc = (0.0_dp,0.0_dp)
    other = freqresp(ac,bc,ccx,dc,w=w)
    call check(other%status == design_status_success,'complex state space analog')
    call check(maxval(abs(other%response-1.0_dp/cmplx(1.0_dp,w-0.3_dp,dp))) < tol, &
        'complex analog state space values')

    resp = dfreqresp([1.0_dp],[1.0_dp,-0.5_dp],w=w)
    call check(resp%status == design_status_success,'digital TF status')
    do i=1,size(w)
        call check(abs(resp%response(i)-1.0_dp/(exp(cmplx(0.0_dp,w(i),dp))-0.5_dp))<tol, &
            'descending z polynomials')
    end do
    other = dfreqresp_zpk([real(dp) ::],[0.5_dp],1.0_dp,w=w)
    call check(maxval(abs(other%response-resp%response))<tol,'digital ZPK')
    plot = dbode([1.0_dp],[1.0_dp,-0.5_dp],w=w,dt=0.25_dp)
    call check(plot%status == design_status_success,'digital Bode status')
    call check(maxval(abs(plot%frequencies-4.0_dp*w))<tol,'digital frequency scale')
    call check(abs(plot%magnitude(1)-20.0_dp*log10(2.0_dp))<tol,'digital Bode dB')
    plot = dbode_zpk([real(dp) ::],[0.5_dp],1.0_dp,w=w,dt=0.25_dp)
    call check(abs(plot%magnitude(1)-20.0_dp*log10(2.0_dp))<tol,'digital ZPK Bode')
    aa = 0.5_dp
    other = dfreqresp(aa,bb,cc,dd,w=w)
    call check(other%status == design_status_success,'discrete state space status')
    call check(maxval(abs(other%response-resp%response)) < tol,'digital state space values')
    plot = dbode(aa,bb,cc,dd,w=w,dt=0.5_dp)
    call check(plot%status == design_status_success,'digital state space Bode status')
    call check(maxval(abs(plot%frequencies-2.0_dp*w))<tol,'digital state-space units')

    ! Unwrapping produces continuous phase when sample spacing is less than pi.
    plot = dbode([1.0_dp],[1.0_dp,0.0_dp,0.0_dp,0.0_dp],w=w)
    call check(plot%status == design_status_success,'pure digital delay status')
    call check(maxval(abs(plot%phase+3.0_dp*w*180.0_dp/pi))<tol,'phase unwrap')

    resp = freqresp([1.0_dp],[1.0_dp,1.0_dp],n=29)
    call check(resp%status == design_status_success,'automatic analog grid')
    call check(size(resp%frequencies)==29,'analog grid count')
    call check(all(resp%frequencies(2:)>resp%frequencies(:28)),'analog log grid')
    resp = dfreqresp([1.0_dp],[1.0_dp,-0.5_dp],n=32,whole=.true.)
    call check(resp%status == design_status_success,'automatic full digital grid')
    call check(size(resp%frequencies)==32,'digital grid count')
    call check(abs(resp%frequencies(32)-2.0_dp*pi*31.0_dp/32.0_dp)<tol,'digital whole grid')
    plot = bode([1.0_dp],[1.0_dp,1.0_dp],n=7)
    call check(size(plot%frequencies)==7,'Bode auto grid count')
    plot = dbode([1.0_dp],[1.0_dp,-0.5_dp],n=9,dt=0.2_dp)
    call check(size(plot%frequencies)==9,'dbode auto grid count')

    resp = dfreqresp([1.0_dp],[1.0_dp,-0.5_dp],dt=-1.0_dp,w=w)
    call check(resp%status == design_status_invalid_input,'negative sample interval')
    plot = dbode([1.0_dp],[1.0_dp,-0.5_dp],dt=0.0_dp,w=w)
    call check(plot%status == design_status_invalid_input,'zero sample interval')
    dm = 0.0_dp
    resp = freqresp(aa,bb,dm,dd,w=w)
    call check(resp%status == design_status_invalid_input,'non-SISO state space')
    resp = dfreqresp([1.0_dp],[0.0_dp],w=w)
    call check(resp%status == design_status_invalid_input,'zero denominator')
    ! Complex TF and ZPK gain cases are checked directly because SciPy's
    ! transfer-function wrapper can discard imaginary numerator components.
    resp = freqresp([cmplx(1.0_dp,2.0_dp,dp)], &
        [cmplx(1.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp)],w=w)
    call check(resp%status == design_status_success,'complex TF analog status')
    call check(maxval(abs(resp%response-cmplx(1.0_dp,2.0_dp,dp)/ &
        cmplx(2.0_dp,w,dp))) < tol,'complex TF analog values')
    resp = dfreqresp([cmplx(1.0_dp,2.0_dp,dp)], &
        [cmplx(1.0_dp,0.0_dp,dp),cmplx(-0.5_dp,0.0_dp,dp)],w=w)
    call check(resp%status == design_status_success,'complex digital TF')
    call check(maxval(abs(resp%response-cmplx(1.0_dp,2.0_dp,dp)/ &
        (exp(cmplx(0.0_dp,w,dp))-0.5_dp))) < tol,'complex digital TF values')
    resp = freqresp_zpk([complex(dp) ::],[-(1.0_dp,0.0_dp)], &
        cmplx(0.4_dp,0.2_dp,dp),w=w)
    call check(resp%status == design_status_success,'complex analog ZPK gain')
    call check(maxval(abs(resp%response-cmplx(0.4_dp,0.2_dp,dp)/ &
        cmplx(1.0_dp,w,dp))) < tol,'complex ZPK frequency response')

    ! Zero gain yields -infinity dB without triggering an invalid division.
    plot = bode([0.0_dp],[1.0_dp,1.0_dp],w=w)
    call check(plot%status == design_status_success,'zero-gain Bode status')
    call check(all(plot%magnitude < 0.0_dp) .and. &
        .not. any(ieee_is_finite(plot%magnitude)),'zero gain is minus infinity dB')
    resp = dfreqresp([1.0_dp],[1.0_dp,-0.5_dp],n=0)
    call check(resp%status == design_status_invalid_input,'invalid auto grid count')

    print '(a)', 'PASS: continuous/discrete LTI frequency response and Bode'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! logical test expression
        character(*), intent(in) :: label !! assertion name
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            stop 1
        end if
    end subroutine check
end program test_signal_lti_frequency
