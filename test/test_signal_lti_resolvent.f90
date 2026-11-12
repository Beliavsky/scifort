! SPDX-License-Identifier: MIT
! Direct state-space frequency responses, including high orders and singularities.
program test_signal_lti_resolvent
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, freqresp, dfreqresp, bode, dbode, &
        frequency_response, bode_result, design_status_success, &
        design_status_invalid_input, design_status_singular
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: a(60,60), b(60,1), c(1,60), d(1,1), w(4)
    real(dp) :: at(2,2), bt(2,1), ct(1,2), dt(1,1)
    real(dp) :: empty_a(0,0), empty_b(0,1), empty_c(1,0)
    complex(dp) :: ac(1,1), bc(1,1), cc(1,1), dc(1,1)
    complex(dp) :: expected, z
    type(frequency_response) :: r
    type(bode_result) :: p
    integer :: i, j

    w = [0.0_dp,0.3_dp,0.8_dp,2.4_dp]
    at = reshape([-2.0_dp,0.0_dp,10.0_dp,-3.0_dp],[2,2])
    bt(:,1) = [0.0_dp,1.0_dp]
    ct(1,:) = [1.0_dp,0.0_dp]
    dt = 0.0_dp
    r = freqresp(at,bt,ct,dt,w=w)
    call check(r%status == design_status_success,'triangular analog')
    do i = 1,size(w)
        z = cmplx(0.0_dp,w(i),dp)
        expected = 10.0_dp/((z+2.0_dp)*(z+3.0_dp))
        call check(abs(r%response(i)-expected) < 1.0e-12_dp,'triangular analog value')
    end do
    at = reshape([0.2_dp,0.0_dp,0.4_dp,0.5_dp],[2,2])
    r = dfreqresp(at,bt,ct,dt,w=w)
    call check(r%status == design_status_success,'triangular digital')
    do i = 1,size(w)
        z = exp(cmplx(0.0_dp,w(i),dp))
        expected = 0.4_dp/((z-0.2_dp)*(z-0.5_dp))
        call check(abs(r%response(i)-expected) < 1.0e-12_dp,'triangular digital value')
    end do

    ! A 60-state diagonal system with characteristic polynomial overflow:
    ! prod(s+1e10) is not representable in real(dp), but the resolvent is.
    a = 0.0_dp
    b = 1.0_dp
    c = 1.0_dp
    d = 0.25_dp
    do j = 1,60
        a(j,j) = -1.0e10_dp
    end do
    w = [0.0_dp,1.0e8_dp,1.0e10_dp,1.0e12_dp]
    r = freqresp(a,b,c,d,w=w)
    call check(r%status == design_status_success,'high-order analog status')
    do i = 1,size(w)
        expected = d(1,1)+60.0_dp/cmplx(1.0e10_dp,w(i),dp)
        call check(abs(r%response(i)-expected) <= 1.0e-12_dp*max(1.0_dp,abs(expected)), &
            'high-order analog resolvent')
    end do
    p = bode(a,b,c,d,w=w)
    call check(p%status == design_status_success,'high-order analog bode')
    call check(all(abs(p%magnitude-20.0_dp*log10(abs(r%response))) < 1.0e-11_dp), &
        'high-order analog bode values')
    r = freqresp(a,b,c,d,n=9)
    call check(r%status == design_status_success,'high-order analog automatic grid')
    call check(size(r%frequencies)==9,'high-order analog automatic size')
    call check(all(r%frequencies(2:)>r%frequencies(:8)), 'high-order analog log grid')

    a = 0.0_dp
    do j = 1,60
        a(j,j) = 0.92_dp
    end do
    w = [0.0_dp,0.2_dp,0.8_dp,pi]
    r = dfreqresp(a,b,c,d,w=w)
    call check(r%status == design_status_success,'high-order digital status')
    do i = 1,size(w)
        expected = d(1,1)+60.0_dp/(exp(cmplx(0.0_dp,w(i),dp))-0.92_dp)
        call check(abs(r%response(i)-expected) <= 2.0e-12_dp*max(1.0_dp,abs(expected)), &
            'high-order digital resolvent')
    end do
    p = dbode(a,b,c,d,w=w,dt=0.2_dp)
    call check(p%status == design_status_success,'high-order dbode')
    call check(maxval(abs(p%frequencies-w/0.2_dp)) < 1.0e-12_dp,'high-order dbode time units')
    r = dfreqresp(a,b,c,d,n=8,whole=.true.)
    call check(r%status == design_status_success,'digital auto grid')
    call check(abs(r%frequencies(8)-7.0_dp*pi/4.0_dp) < 1.0e-12_dp,'full-circle grid')

    ac(1,1) = cmplx(-2.0_dp,0.8_dp,dp)
    bc(1,1) = cmplx(1.0_dp,0.4_dp,dp)
    cc(1,1) = cmplx(0.3_dp,-0.2_dp,dp)
    dc(1,1) = cmplx(0.1_dp,0.3_dp,dp)
    w = [0.0_dp,0.3_dp,0.8_dp,2.4_dp]
    r = freqresp(ac,bc,cc,dc,w=w)
    call check(r%status == design_status_success,'complex analog status')
    do i = 1,size(w)
        expected = dc(1,1)+cc(1,1)*bc(1,1)/ &
            (cmplx(0.0_dp,w(i),dp)-ac(1,1))
        call check(abs(r%response(i)-expected) < 1.0e-12_dp,'complex analog resolvent')
    end do
    r = dfreqresp(ac,bc,cc,dc,w=w)
    call check(r%status == design_status_success,'complex digital status')
    do i = 1,size(w)
        expected = dc(1,1)+cc(1,1)*bc(1,1)/ &
            (exp(cmplx(0.0_dp,w(i),dp))-ac(1,1))
        call check(abs(r%response(i)-expected) < 1.0e-12_dp,'complex digital resolvent')
    end do

    ! No states: pure feedthrough, including automatic response generation.
    dt = 3.0_dp
    r = freqresp(empty_a,empty_b,empty_c,dt,w=w)
    call check(r%status == design_status_success,'zero states analog status')
    call check(all(abs(r%response-3.0_dp)<1.0e-14_dp),'zero states analog gain')
    r = dfreqresp(empty_a,empty_b,empty_c,dt,n=11)
    call check(r%status == design_status_success,'zero states digital status')
    call check(size(r%frequencies)==11,'zero states digital samples')
    call check(all(abs(r%response-3.0_dp)<1.0e-14_dp),'zero states digital gain')

    ! A pole on the evaluation contour yields a singular status, not a crash.
    at = 0.0_dp
    bt = 1.0_dp
    ct = 1.0_dp
    dt = 0.0_dp
    r = freqresp(at,bt,ct,dt,w=w)
    call check(r%status == design_status_singular,'analog singular status')
    at = 0.0_dp
    at(1,1) = 1.0_dp
    at(2,2) = 0.4_dp
    r = dfreqresp(at,bt,ct,dt,w=w)
    call check(r%status == design_status_singular,'digital singular status')
    r = dfreqresp(at,bt,ct,dt,w=w,dt=-0.01_dp)
    call check(r%status == design_status_invalid_input,'invalid sampling interval')
    r = freqresp(at,bt,ct,dt,n=0)
    call check(r%status == design_status_invalid_input,'invalid frequency count')
    at(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = freqresp(at,bt,ct,dt,w=w)
    call check(r%status == design_status_invalid_input,'nonfinite state matrix')

    print '(a)', 'PASS: direct state-space resolvents'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion expression
        character(*), intent(in) :: label !! assertion identifier
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            stop 1
        end if
    end subroutine check
end program test_signal_lti_resolvent
