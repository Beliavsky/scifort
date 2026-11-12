! SPDX-License-Identifier: MIT
! Complex and real multiple-input/multiple-output LTI frequency responses.
program test_signal_lti_mimo
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, &
        ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, freqresp_mimo, dfreqresp_mimo, &
        bode_mimo, dbode_mimo, freqresp, dfreqresp, &
        mimo_frequency_response, mimo_bode_result, frequency_response, &
        design_status_success, design_status_invalid_input, design_status_singular
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tol = 2.0e-12_dp
    real(dp) :: a(2,2), b(2,3), c(2,2), d(2,3), w(5)
    real(dp) :: ah(60,60), bh(60,3), ch(2,60), dh(2,3)
    real(dp) :: a0(0,0), b0(0,2), c0(3,0), d0(3,2)
    real(dp) :: aa(1,1), bb(1,1), cc(1,1), dd(1,1)
    real(dp) :: ai(2,2), bi(2,2), ci(2,2), di(2,2)
    complex(dp) :: ac(1,1), bc(1,2), cx(2,1), dc(2,2)
    complex(dp) :: expected, z, x1, x2, tmp(2,3)
    real(dp) :: wh(3)
    type(mimo_frequency_response) :: r
    type(mimo_bode_result) :: plot
    type(frequency_response) :: scalar
    integer :: i, j, k

    w = [0.0_dp, 0.3_dp, 0.8_dp, 1.7_dp, 3.0_dp]
    a = reshape([-2.0_dp,0.0_dp,1.5_dp,-3.0_dp],[2,2])
    b = reshape([1.0_dp,2.0_dp,3.0_dp,-1.0_dp,-2.0_dp,0.3_dp],[2,3])
    c = reshape([1.0_dp,0.2_dp,0.5_dp,-1.0_dp],[2,2])
    d = reshape([0.1_dp,0.4_dp,-0.3_dp,0.2_dp,0.0_dp,0.5_dp],[2,3])
    r = freqresp_mimo(a,b,c,d,w=w)
    call check(r%status == design_status_success,'analog status')
    call check(all(shape(r%response) == [2,3,5]),'analog response shape')
    do k = 1,size(w)
        z = cmplx(0.0_dp,w(k),dp)
        do j = 1,3
            x2 = b(2,j)/(z+3.0_dp)
            x1 = (b(1,j)+1.5_dp*x2)/(z+2.0_dp)
            call check(maxval(abs(r%response(:,j,k)- &
                [c(1,1)*x1+c(1,2)*x2+d(1,j), &
                 c(2,1)*x1+c(2,2)*x2+d(2,j)])) < tol,'analog channels')
        end do
    end do
    plot = bode_mimo(a,b,c,d,w=w)
    call check(plot%status == design_status_success,'analog bode status')
    call check(all(shape(plot%magnitude) == [2,3,5]),'analog bode shape')
    call check(maxval(abs(plot%magnitude-20.0_dp*log10(abs(r%response)))) < tol, &
        'analog bode magnitude')
    r = dfreqresp_mimo(a,b,c,d,w=w)
    call check(r%status == design_status_success,'digital status')
    do k = 1,size(w)
        z = exp(cmplx(0.0_dp,w(k),dp))
        do j = 1,3
            x2 = b(2,j)/(z+3.0_dp)
            x1 = (b(1,j)+1.5_dp*x2)/(z+2.0_dp)
            call check(maxval(abs(r%response(:,j,k)- &
                [c(1,1)*x1+c(1,2)*x2+d(1,j), &
                 c(2,1)*x1+c(2,2)*x2+d(2,j)])) < tol,'digital channels')
        end do
    end do
    plot = dbode_mimo(a,b,c,d,w=w,dt=0.25_dp)
    call check(plot%status == design_status_success,'digital bode status')
    call check(maxval(abs(plot%frequencies-w/0.25_dp)) < tol,'digital bode units')
    r = dfreqresp_mimo(a,b,c,d,n=8,whole=.true.)
    call check(r%status == design_status_success,'digital automatic grid')
    call check(all(shape(r%response) == [2,3,8]),'digital automatic shape')
    call check(abs(r%frequencies(8)-7.0_dp*pi/4.0_dp) < tol,'digital full-circle')
    r = freqresp_mimo(a,b,c,d,n=11)
    call check(size(r%frequencies) == 11,'analog automatic count')
    call check(all(r%frequencies(2:)>r%frequencies(:10)),'analog automatic log grid')
    plot = bode_mimo(a,b,c,d,n=12)
    call check(size(plot%frequencies) == 12,'analog bode automatic grid')
    plot = dbode_mimo(a,b,c,d,n=12,dt=0.2_dp)
    call check(size(plot%frequencies) == 12,'digital bode automatic grid')

    ! High-order resolvent: characteristic polynomial would overflow real(dp).
    ah = 0.0_dp
    bh = 0.0_dp
    ch = 0.0_dp
    do i = 1,60
        ah(i,i) = -1.0e10_dp
        bh(i,:) = [real(mod(i,7),dp),real(mod(i,5),dp),1.0_dp]
        ch(:,i) = [real(mod(i,4),dp),real(mod(i,3),dp)]
    end do
    dh = reshape([0.1_dp,0.2_dp,0.3_dp,0.4_dp,0.5_dp,0.6_dp],[2,3])
    wh = [0.0_dp,1.0e10_dp,1.0e12_dp]
    r = freqresp_mimo(ah,bh,ch,dh,w=wh)
    call check(r%status == design_status_success,'high-order status')
    tmp = matmul(cmplx(ch,0.0_dp,dp),cmplx(bh,0.0_dp,dp))
    do k = 1,3
        expected = cmplx(1.0e10_dp,wh(k),dp)
        call check(maxval(abs(r%response(:,:,k)-dh-tmp/expected)) < 2.0e-11_dp, &
            'high-order response')
    end do

    ! Complex coefficients and two inputs, two outputs.
    ac(1,1) = cmplx(-1.0_dp,0.4_dp,dp)
    bc(1,:) = [cmplx(1.0_dp,0.2_dp,dp),cmplx(-0.3_dp,1.0_dp,dp)]
    cx(:,1) = [cmplx(0.8_dp,-0.1_dp,dp),cmplx(-0.2_dp,0.5_dp,dp)]
    dc = cmplx(0.0_dp,0.0_dp,dp)
    dc(1,2) = cmplx(0.1_dp,-0.3_dp,dp)
    r = freqresp_mimo(ac,bc,cx,dc,w=w)
    call check(r%status == design_status_success,'complex analog status')
    do k = 1,5
        z = cmplx(0.0_dp,w(k),dp)
        do j = 1,2
            do i = 1,2
                expected = cx(i,1)*bc(1,j)/(z-ac(1,1))+dc(i,j)
                call check(abs(r%response(i,j,k)-expected) < tol,'complex analog value')
            end do
        end do
    end do
    r = dfreqresp_mimo(ac,bc,cx,dc,w=w)
    do k = 1,5
        z = exp(cmplx(0.0_dp,w(k),dp))
        do j = 1,2
            do i = 1,2
                expected = cx(i,1)*bc(1,j)/(z-ac(1,1))+dc(i,j)
                call check(abs(r%response(i,j,k)-expected) < tol,'complex digital value')
            end do
        end do
    end do

    ! Compare 1x1 MIMO with the backwards-compatible SISO API.
    aa = -2.0_dp
    bb = 0.3_dp
    cc = 1.7_dp
    dd = 0.4_dp
    r = freqresp_mimo(aa,bb,cc,dd,w=w)
    scalar = freqresp(aa,bb,cc,dd,w=w)
    call check(maxval(abs(r%response(1,1,:)-scalar%response)) < tol,'SISO analog parity')
    r = dfreqresp_mimo(aa,bb,cc,dd,w=w)
    scalar = dfreqresp(aa,bb,cc,dd,w=w)
    call check(maxval(abs(r%response(1,1,:)-scalar%response)) < tol,'SISO digital parity')

    ! Zero-state feedthrough response, including multiple input/output channels.
    d0 = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[3,2])
    r = freqresp_mimo(a0,b0,c0,d0,w=w)
    call check(r%status == design_status_success,'static analog')
    do k = 1,5
        call check(all(abs(r%response(:,:,k)-d0) < tol),'static feedthrough')
    end do
    r = dfreqresp_mimo(a0,b0,c0,d0,n=7)
    call check(size(r%frequencies)==7,'static digital frequency count')
    call check(all(abs(r%response(:,:,7)-d0) < tol),'static digital feedthrough')

    ! For a unit delay in two uncoupled channels, phase unwrap is independent.
    ai = 0.0_dp
    bi = 0.0_dp
    ci = 0.0_dp
    di = 0.0_dp
    bi(1,1) = 1.0_dp
    bi(2,2) = 1.0_dp
    ci(1,1) = 1.0_dp
    ci(2,2) = 1.0_dp
    w = [0.0_dp,1.5_dp,3.0_dp,4.5_dp,6.0_dp]
    plot = dbode_mimo(ai,bi,ci,di,w=w,dt=0.1_dp)
    call check(plot%status == design_status_success,'delay bode status')
    call check(maxval(abs(plot%phase(1,1,:)+w*180.0_dp/pi)) < tol,'delay phase unwrap')
    call check(maxval(abs(plot%phase(2,2,:)+w*180.0_dp/pi)) < tol,'independent delay phase')
    call check(all(plot%magnitude(1,2,:) < 0.0_dp),'zero channel -infinity dB')
    call check(.not. any(ieee_is_finite(plot%magnitude(1,2,:))),'zero channel not finite dB')

    ! Failure at a frequency does not discard values at other frequencies.
    ai = 0.0_dp
    bi = 1.0_dp
    ci = 1.0_dp
    di = 0.0_dp
    r = freqresp_mimo(ai,bi,ci,di,w=[0.0_dp,1.0_dp])
    call check(r%status == design_status_singular,'analog singular status')
    call check(all(ieee_is_nan(real(r%response(:,:,1),dp))),'singular nan slice')
    call check(all(ieee_is_finite(real(r%response(:,:,2),dp))),'nonsingular slice kept')
    plot = bode_mimo(ai,bi,ci,di,w=[0.0_dp,1.0_dp])
    call check(plot%status == design_status_singular,'singular bode status')
    call check(all(ieee_is_nan(plot%magnitude(:,:,1))),'singular bode nan slice')
    call check(all(ieee_is_finite(plot%magnitude(:,:,2))),'singular bode surviving slice')
    ai(1,1) = 1.0_dp
    ai(2,2) = 0.5_dp
    r = dfreqresp_mimo(ai,bi,ci,di,w=[0.0_dp,0.8_dp])
    call check(r%status == design_status_singular,'digital singular status')
    call check(all(ieee_is_nan(real(r%response(:,:,1),dp))),'digital singular slice')

    r = freqresp_mimo(a,b,c,d,n=0)
    call check(r%status == design_status_invalid_input,'invalid frequency count')
    r = dfreqresp_mimo(a,b,c,d,n=6,dt=-0.1_dp)
    call check(r%status == design_status_invalid_input,'invalid sample interval')
    r = freqresp_mimo(a,b,c,d(:1,:),w=w)
    call check(r%status == design_status_invalid_input,'invalid D dimensions')
    r = freqresp_mimo(a,b,c(:,:1),d,w=w)
    call check(r%status == design_status_invalid_input,'invalid C dimensions')
    ai(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = freqresp_mimo(ai,bi,ci,di,w=w)
    call check(r%status == design_status_invalid_input,'nonfinite matrix')
    r = freqresp_mimo(a,b,c,d,w=[0.1_dp,ieee_value(0.0_dp,ieee_quiet_nan)])
    call check(r%status == design_status_invalid_input,'nonfinite frequencies')

    print '(a)', 'PASS: MIMO continuous/discrete resolvent frequency responses'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion expression
        character(*), intent(in) :: label !! test label
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            stop 1
        end if
    end subroutine check
end program test_signal_lti_mimo
