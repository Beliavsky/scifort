! SPDX-License-Identifier: MIT
! Prepared Hessenberg sweeps and one-based MIMO channel selection.
program test_signal_lti_mimo_plan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, prepare_mimo_frequency, mimo_frequency_plan, &
        freqresp_mimo, dfreqresp_mimo, bode_mimo, dbode_mimo, &
        mimo_frequency_response, mimo_bode_result, &
        design_status_success, design_status_invalid_input, design_status_singular
    implicit none
    integer, parameter :: ns = 12, ni = 4, no = 3, nw = 18
    integer, parameter :: chosen_inputs(3) = [4,2,4], chosen_outputs(2) = [3,1]
    real(dp) :: a(ns,ns), b(ns,ni), c(no,ns), d(no,ni), w(nw)
    real(dp) :: aa(1,1), bb(1,2), cc(2,1), dd(2,2)
    real(dp) :: a0(0,0), b0(0,2), c0(3,0), d0(3,2)
    complex(dp) :: ac(ns,ns), bc(ns,ni), dc(no,ni), xc(no,ns)
    type(mimo_frequency_response) :: complete, selected, reused, digital, target
    type(mimo_bode_result) :: plot1, plot2
    type(mimo_frequency_plan) :: prepared, complex_plan, invalid_plan, static_plan
    integer :: i, j, k
    real(dp), parameter :: tolerance = 3.0e-11_dp

    do j = 1,ns
        do i = 1,ns
            a(i,j) = 0.02_dp*cos(real(3*i+j,dp))
            if (i == j) a(i,j) = a(i,j)-0.6_dp-0.15_dp*real(i,dp)
            ac(i,j) = cmplx(a(i,j),0.01_dp*sin(real(i+3*j,dp)),dp)
        end do
    end do
    do j = 1,ni
        do i = 1,ns
            b(i,j) = cos(real(i+2*j,dp))
            bc(i,j) = cmplx(b(i,j),0.1_dp*sin(real(i+4*j,dp)),dp)
        end do
    end do
    do j = 1,ns
        do i = 1,no
            c(i,j) = sin(real(2*i+3*j,dp))
            xc(i,j) = cmplx(c(i,j),0.09_dp*cos(real(4*i+j,dp)),dp)
        end do
    end do
    do j = 1,ni
        do i = 1,no
            d(i,j) = 0.05_dp*real(i+j,dp)
            dc(i,j) = cmplx(d(i,j),0.01_dp*real(j-i,dp),dp)
        end do
    end do
    w = [(0.05_dp*real(k-1,dp),k=1,nw)]

    prepared = prepare_mimo_frequency(a,b,c,d)
    call check(prepared%status == design_status_success,'real plan construction')
    complete = freqresp_mimo(a,b,c,d,w=w,solver='direct')
    selected = freqresp_mimo(a,b,c,d,w=w,inputs=chosen_inputs, &
        outputs=chosen_outputs,solver='hessenberg')
    reused = freqresp_mimo(prepared,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    call check(selected%status == design_status_success,'real subset status')
    call check(reused%status == design_status_success,'real prepared subset status')
    call check(all(shape(selected%response) == [2,3,nw]),'subset dimensions')
    target%response = complete%response(chosen_outputs,chosen_inputs,:)
    call check(maxval(abs(selected%response-target%response)) < tolerance,'channel order and repetition')
    call check(maxval(abs(reused%response-target%response)) < tolerance,'prepared analog channels')
    call check(maxval(abs(reused%response-selected%response)) < tolerance,'prepared vs direct API')

    digital = dfreqresp_mimo(a,b,c,d,w=w,solver='direct')
    selected = dfreqresp_mimo(a,b,c,d,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    reused = dfreqresp_mimo(prepared,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    target%response = digital%response(chosen_outputs,chosen_inputs,:)
    call check(maxval(abs(selected%response-target%response)) < tolerance,'selected digital channels')
    call check(maxval(abs(reused%response-target%response)) < tolerance,'prepared digital channels')
    plot1 = bode_mimo(a,b,c,d,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    plot2 = bode_mimo(prepared,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    call check(maxval(abs(plot1%magnitude-plot2%magnitude)) < 1.0e-9_dp,'analog Bode magnitude')
    call check(maxval(abs(plot1%phase-plot2%phase)) < 1.0e-9_dp,'analog Bode phase')
    plot1 = dbode_mimo(a,b,c,d,w=w,dt=0.2_dp,inputs=chosen_inputs,outputs=chosen_outputs)
    plot2 = dbode_mimo(prepared,w=w,dt=0.2_dp,inputs=chosen_inputs,outputs=chosen_outputs)
    call check(maxval(abs(plot1%magnitude-plot2%magnitude)) < 1.0e-9_dp,'digital Bode magnitude')
    call check(maxval(abs(plot1%phase-plot2%phase)) < 1.0e-9_dp,'digital Bode phase')
    call check(maxval(abs(plot2%frequencies-w/0.2_dp)) < tolerance,'dt frequency conversion')

    complex_plan = prepare_mimo_frequency(ac,bc,xc,dc)
    call check(complex_plan%status == design_status_success,'complex plan construction')
    complete = freqresp_mimo(ac,bc,xc,dc,w=w,solver='direct')
    reused = freqresp_mimo(complex_plan,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    call check(maxval(abs(reused%response-complete%response(chosen_outputs,chosen_inputs,:))) &
        < tolerance,'complex analog prepared')
    complete = dfreqresp_mimo(ac,bc,xc,dc,w=w,solver='direct')
    reused = dfreqresp_mimo(complex_plan,w=w,inputs=chosen_inputs,outputs=chosen_outputs)
    call check(maxval(abs(reused%response-complete%response(chosen_outputs,chosen_inputs,:))) &
        < tolerance,'complex digital prepared')

    ! One preparation can serve disjoint frequency sweeps and channel sets.
    selected = freqresp_mimo(prepared,w=w(:6),inputs=[1])
    complete = freqresp_mimo(a,b,c,d,w=w(:6),solver='direct',inputs=[1])
    call check(maxval(abs(selected%response-complete%response)) < tolerance,'first sweep')
    selected = freqresp_mimo(prepared,w=w(7:),outputs=[2])
    complete = freqresp_mimo(a,b,c,d,w=w(7:),solver='direct',outputs=[2])
    call check(maxval(abs(selected%response-complete%response)) < tolerance,'second sweep')

    complete = freqresp_mimo(a,b,c,d,n=19,solver='hessenberg')
    reused = freqresp_mimo(prepared,n=19)
    call check(all(complete%frequencies == reused%frequencies),'same automatic analog grid')
    call check(maxval(abs(complete%response-reused%response)) < tolerance, &
        'same automatic analog response')

    selected = freqresp_mimo(a,b,c,d,w=w,inputs=[0])
    call check(selected%status == design_status_invalid_input,'zero input index rejected')
    selected = freqresp_mimo(a,b,c,d,w=w,outputs=[no+1])
    call check(selected%status == design_status_invalid_input,'large output index rejected')
    selected = freqresp_mimo(prepared,w=w,inputs=[0])
    call check(selected%status == design_status_invalid_input,'prepared invalid selector')
    selected = freqresp_mimo(prepared,w=w,outputs=[integer ::])
    call check(selected%status == design_status_invalid_input,'empty output indices')
    selected = dfreqresp_mimo(prepared,w=w,dt=-1.0_dp)
    call check(selected%status == design_status_invalid_input,'prepared bad dt')
    selected = freqresp_mimo(prepared,n=0)
    call check(selected%status == design_status_invalid_input,'prepared bad frequency count')
    invalid_plan = prepare_mimo_frequency(a,b,c,d(:1,:))
    call check(invalid_plan%status == design_status_invalid_input,'invalid model rejected by plan')
    selected = freqresp_mimo(invalid_plan,w=w)
    call check(selected%status == design_status_invalid_input,'unprepared evaluation rejected')
    call check(.not. allocated(selected%response),'unprepared evaluation not allocated')

    d0 = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[3,2])
    static_plan = prepare_mimo_frequency(a0,b0,c0,d0)
    call check(static_plan%status == design_status_success,'static plan')
    selected = freqresp_mimo(static_plan,w=w(:3),inputs=[2],outputs=[3,1])
    call check(all(abs(selected%response(1,1,:)-d0(3,2)) < tolerance),'static first channel')
    call check(all(abs(selected%response(2,1,:)-d0(1,2)) < tolerance),'static second channel')
    aa = 0.0_dp
    bb = 1.0_dp
    cc = 1.0_dp
    dd = 0.0_dp
    prepared = prepare_mimo_frequency(aa,bb,cc,dd)
    selected = freqresp_mimo(prepared,w=[0.0_dp,1.0_dp],outputs=[2])
    call check(selected%status == design_status_singular,'prepared singular status')
    call check(ieee_is_nan(real(selected%response(1,1,1),dp)),'prepared singular frequency')
    call check(abs(selected%response(1,1,2)-cmplx(0.0_dp,-1.0_dp,dp)) < tolerance, &
        'prepared finite frequency')

    print '(a)', 'PASS: reusable Hessenberg MIMO plans and channel selection'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! predicate to test
        character(len=*), intent(in) :: message !! diagnostic
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check
end program test_signal_lti_mimo_plan
