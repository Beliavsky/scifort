! SPDX-License-Identifier: MIT
! Analog object-based transfer, response, interconnection and discretization.
program test_signal_clti
    use scifort_signal, only : dp, clti_system, clti_tf, clti_zpk, clti_ss, &
        clti_success, clti_invalid_input, dlti_system, dlti_success, &
        frequency_response, mimo_frequency_response, &
        continuous_simulation_result, continuous_response_result
    implicit none
    type(clti_system) :: first,zpk,ss,parallel,series,feedback,complex_sys,bad,static
    type(dlti_system) :: disc
    type(frequency_response) :: response
    type(mimo_frequency_response) :: multi
    type(continuous_simulation_result) :: output
    type(continuous_response_result) :: steps,impulses
    complex(dp) :: a(2,2),b(2,2),c(2,2),d(2,2), &
        complex_expected
    real(dp) :: w(4),t(31),u(31)
    integer :: i,j

    first = clti_tf([2.0_dp],[1.0_dp,3.0_dp])
    call check(first%status() == clti_success,'tf valid')
    call check(first%nstates() == 1 .and. first%ninputs() == 1,'tf dimensions')
    call check(first%stable(),'one-pole continuous stable')
    w = [0.0_dp,0.3_dp,1.2_dp,10.0_dp]
    response = first%freqresp(w=w)
    call check(response%status == 0,'frequency status')
    do i=1,size(w)
        call check(abs(response%response(i)-2.0_dp/cmplx(3.0_dp,w(i),dp)) < 1.0e-13_dp, &
            'analog frequency response')
    end do
    zpk = first%to_zpk()
    call check(zpk%status() == clti_success,'ZPK conversion')
    response = zpk%freqresp(w=w)
    call check(maxval(abs(response%response-[ &
        (2.0_dp/cmplx(3.0_dp,w(i),dp),i=1,4)])) < 1.0e-11_dp,'ZPK spectrum')
    ss = first%to_ss()
    call check(ss%status() == clti_success,'SS conversion')
    do i=1,size(t)
        t(i) = real(i-1,dp)/10.0_dp
    end do
    u = 1.0_dp
    output = ss%output(u,t)
    call check(output%status == 0,'continuous lsim')
    do i=1,size(t)
        call check(abs(output%y(i,1)-(2.0_dp/3.0_dp)*(1.0_dp-exp(-3.0_dp*t(i)))) < &
            5.0e-13_dp,'constant input output')
    end do
    steps = first%step(t=t)
    impulses = first%impulse(t=t)
    call check(steps%status == 0 .and. impulses%status == 0,'transients')
    do i=1,size(t)
        call check(abs(steps%y(i,1,1)-output%y(i,1)) < 2.0e-12_dp,'unit step')
        call check(abs(impulses%y(i,1,1)-2.0_dp*exp(-3.0_dp*t(i))) < &
            2.0e-12_dp,'impulse')
    end do
    disc = first%to_discrete(0.2_dp)
    call check(disc%status == dlti_success,'ZOH discretize')
    response = disc%freqresp(w=[0.0_dp])
    call check(abs(response%response(1)-2.0_dp/3.0_dp) < 1.0e-12_dp,'ZOH dc gain')
    a = (0.0_dp,0.0_dp)
    b = (0.0_dp,0.0_dp)
    c = (0.0_dp,0.0_dp)
    d = (0.0_dp,0.0_dp)
    a(1,1) = cmplx(-1.0_dp,0.25_dp,dp)
    a(2,2) = cmplx(-3.0_dp,-0.4_dp,dp)
    b(1,1) = (1.0_dp,0.0_dp)
    b(2,2) = (1.0_dp,0.0_dp)
    c(1,1) = cmplx(0.7_dp,0.2_dp,dp)
    c(2,2) = cmplx(0.4_dp,-0.1_dp,dp)
    d(1,1) = cmplx(0.2_dp,0.05_dp,dp)
    complex_sys = clti_ss(a,b,c,d)
    call check(complex_sys%status() == clti_success,'complex MIMO')
    call check(complex_sys%stable(),'complex MIMO continuous stability')
    multi = complex_sys%freqresp_mimo(w=w,solver='direct')
    call check(multi%status == 0,'MIMO response')
    do i=1,size(w)
        do j=1,2
            complex_expected = d(j,j)+c(j,j)/(cmplx(0.0_dp,w(i),dp)-a(j,j))
            call check(abs(multi%response(j,j,i)-complex_expected) < 2.0e-13_dp, &
                'complex MIMO diagonal response')
        end do
    end do
    parallel = complex_sys%parallel(complex_sys)
    multi = parallel%freqresp_mimo(w=w,solver='direct')
    complex_expected = 2.0_dp*(d(1,1)+c(1,1)/(cmplx(0.0_dp,w(2),dp)-a(1,1)))
    call check(abs(multi%response(1,1,2)-complex_expected) < 5.0e-13_dp, &
        'parallel continuous systems')
    series = complex_sys%series(complex_sys)
    multi = series%freqresp_mimo(w=w,solver='direct')
    complex_expected = (d(1,1)+c(1,1)/(cmplx(0.0_dp,w(2),dp)-a(1,1)))**2
    call check(abs(multi%response(1,1,2)-complex_expected) < 5.0e-13_dp, &
        'series continuous systems')
    static = clti_ss(a(:0,:0),b(:0,:),c(:,:0),d)
    feedback = static%feedback(static)
    call check(feedback%status() == clti_success,'continuous feedback')
    multi = feedback%freqresp_mimo(w=w,solver='direct')
    call check(abs(multi%response(1,1,1)-d(1,1)/(1.0_dp+d(1,1)**2)) < &
        1.0e-13_dp,'static negative feedback')
    bad = clti_tf([1.0_dp],[1.0_dp,-2.0_dp])
    call check(.not. bad%stable(),'unstable right-half-plane pole')
    bad = clti_tf([1.0_dp],[0.0_dp,0.0_dp])
    call check(bad%status() == clti_invalid_input,'invalid analog denominator')
    print '(a)', 'PASS: continuous-time LTI conversions, simulation, MIMO and interconnection'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! asserted condition
        character(*), intent(in) :: label !! diagnostic
        if (.not. ok) then
            print '(a,1x,a)', 'FAIL:',label
            error stop 1
        end if
    end subroutine check
end program test_signal_clti
