! SPDX-License-Identifier: MIT
! Regression for dependency-free complex Schur response plans.
program test_signal_lti_schur
    use scifort_signal, only : dp, freqresp_mimo, dfreqresp_mimo, &
        bode_mimo, dbode_mimo, prepare_mimo_frequency, mimo_frequency_plan, &
        mimo_frequency_response, mimo_bode_result, design_status_success, &
        design_status_invalid_input, design_status_singular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    implicit none
    integer, parameter :: ns = 20, ni = 3, no = 2, nw = 33
    real(dp), parameter :: tol = 3.0e-10_dp
    real(dp) :: a(ns,ns), b(ns,ni), c(no,ns), d(no,ni), w(nw)
    complex(dp) :: ac(ns,ns), bc(ns,ni), cc(no,ns), dc(no,ni)
    real(dp) :: a1(1,1), b1(1,2), c1(2,1), d1(2,2)
    real(dp) :: a0(0,0), b0(0,2), c0(2,0), d0(2,2)
    type(mimo_frequency_plan) :: plan, cp, bad, static_plan
    type(mimo_frequency_response) :: base, hs, ss, prepared
    type(mimo_bode_result) :: bd, pb
    integer :: i, j, k

    do j = 1,ns
        do i = 1,ns
            a(i,j) = 0.03_dp*sin(real(5*i+2*j,dp))
            if (i == j) a(i,j) = a(i,j)-0.25_dp-real(i,dp)*0.05_dp
            ac(i,j) = cmplx(a(i,j),0.02_dp*cos(real(i+4*j,dp)),dp)
        end do
    end do
    do j = 1,ni
        do i = 1,ns
            b(i,j) = cos(real(i+3*j,dp))
            bc(i,j) = cmplx(b(i,j),0.03_dp*sin(real(2*i+j,dp)),dp)
        end do
    end do
    do j = 1,ns
        do i = 1,no
            c(i,j) = sin(real(i+5*j,dp))
            cc(i,j) = cmplx(c(i,j),0.01_dp*cos(real(i+j,dp)),dp)
        end do
    end do
    do j = 1,ni
        do i = 1,no
            d(i,j) = 0.01_dp*real(i+j,dp)
            dc(i,j) = cmplx(d(i,j),0.03_dp*sin(real(i+j,dp)),dp)
        end do
    end do
    w = [(0.04_dp*real(k-1,dp),k=1,nw)]

    plan = prepare_mimo_frequency(a,b,c,d,solver='schur')
    call check(plan%status == design_status_success,'prepared real valid')
    call check(plan%method == 'schur','real Schur converged')
    call check(all(abs([(plan%h(i,i-1),i=2,ns)]) <= 1.0e-12_dp), 'triangular Schur form')
    base = freqresp_mimo(a,b,c,d,w=w,solver='direct')
    hs = freqresp_mimo(a,b,c,d,w=w,solver='hessenberg')
    ss = freqresp_mimo(a,b,c,d,w=w,solver='schur')
    prepared = freqresp_mimo(plan,w=w)
    call check(base%status == design_status_success,'direct status')
    call close_response(base,hs,'hessenberg response')
    call close_response(base,ss,'Schur response')
    call close_response(base,prepared,'prepared Schur response')
    prepared = freqresp_mimo(plan,w=w,inputs=[3,1],outputs=[2])
    call check(maxval(abs(prepared%response-base%response([2],[3,1],:))) < tol, &
        'Schur selected channels')

    base = dfreqresp_mimo(a,b,c,d,w=w,solver='direct')
    prepared = dfreqresp_mimo(plan,w=w)
    call close_response(base,prepared,'discrete prepared Schur')
    ss = dfreqresp_mimo(a,b,c,d,w=w,solver='schur')
    call close_response(base,ss,'discrete direct Schur')
    bd = bode_mimo(a,b,c,d,w=w,solver='direct')
    pb = bode_mimo(plan,w=w)
    call check(maxval(abs(bd%magnitude-pb%magnitude)) < 1.0e-8_dp,'Schur Bode magnitude')
    call check(maxval(abs(modulo(bd%phase-pb%phase+180.0_dp,360.0_dp)-180.0_dp)) < 1.0e-7_dp,'Schur Bode phase')
    bd = dbode_mimo(a,b,c,d,w=w,dt=0.25_dp,solver='direct')
    pb = dbode_mimo(plan,w=w,dt=0.25_dp)
    call check(maxval(abs(bd%magnitude-pb%magnitude)) < 1.0e-8_dp,'digital Schur Bode')

    cp = prepare_mimo_frequency(ac,bc,cc,dc,solver='schur')
    call check(cp%status == design_status_success,'complex plan')
    call check(cp%method == 'schur','complex Schur converged')
    base = freqresp_mimo(ac,bc,cc,dc,w=w,solver='direct')
    prepared = freqresp_mimo(cp,w=w)
    call close_response(base,prepared,'complex analog Schur')
    base = dfreqresp_mimo(ac,bc,cc,dc,w=w,solver='direct')
    prepared = dfreqresp_mimo(cp,w=w)
    call close_response(base,prepared,'complex digital Schur')

    bad = prepare_mimo_frequency(a,b,c,d,solver='unsupported')
    call check(bad%status == design_status_invalid_input,'invalid solver option')
    plan = prepare_mimo_frequency(a,b,c,d)
    call check(plan%method == 'hessenberg','default Hessenberg plan')

    a1 = 0.0_dp
    b1 = 1.0_dp
    c1 = 1.0_dp
    d1 = 0.0_dp
    plan = prepare_mimo_frequency(a1,b1,c1,d1,solver='schur')
    call check(plan%method == 'schur','scalar Schur')
    ss = freqresp_mimo(plan,w=[0.0_dp,1.0_dp])
    call check(ss%status == design_status_singular,'scalar singular frequency')
    call check(ieee_is_nan(real(ss%response(1,1,1),dp)),'scalar nan on singular')

    d0 = 0.2_dp
    static_plan = prepare_mimo_frequency(a0,b0,c0,d0,solver='schur')
    call check(static_plan%method == 'schur','static Schur')
    ss = freqresp_mimo(static_plan,w=w)
    call check(all(abs(ss%response-cmplx(0.2_dp,0.0_dp,dp)) < tol), 'static response')

    print '(a)', 'PASS: optional complex Schur MIMO frequency plans'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! expected true condition
        character(len=*), intent(in) :: message !! failure diagnostic
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check

    subroutine close_response(x,y,message)
        type(mimo_frequency_response), intent(in) :: x !! reference response
        type(mimo_frequency_response), intent(in) :: y !! comparison response
        character(len=*), intent(in) :: message !! failure diagnostic
        call check(y%status == design_status_success,message//' status')
        call check(maxval(abs(x%response-y%response)) < tol,message)
    end subroutine close_response
end program test_signal_lti_schur
