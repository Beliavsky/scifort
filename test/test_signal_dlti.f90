! SPDX-License-Identifier: MIT
! Discrete-time LTI system constructors, conversions and inherited kernels.
program test_signal_dlti
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti, dlti_tf, dlti_zpk, dlti_ss, dlti_system, &
        dlti_success, dlti_invalid_input, dlti_transfer_function, &
        dlti_zeros_poles_gain, dlti_state_space
    use scifort_signal_discrete_simulation, only : dlsim, dstep, dimpulse, &
        discrete_simulation_result, discrete_response_result, discrete_sim_success
    use scifort_signal_lti_frequency, only : dfreqresp, dfreqresp_mimo
    use scifort_signal_design, only : frequency_response
    implicit none
    type(dlti_system) :: tf, ss, zp, direct, converted, invalid, mimo, unstable, &
        cascade, parallel_system
    type(discrete_simulation_result) :: observed, reference, upstream, downstream
    type(discrete_response_result) :: response, expected
    type(frequency_response) :: frequency, freq_reference
    complex(dp) :: a(2,2), b(2,2), c(2,2), d(2,2), ucomplex(8,2)
    real(dp) :: u(20), w(5), ur(8,2), t(20)
    integer :: i

    tf = dlti([0.3_dp,0.1_dp],[1.0_dp,-0.25_dp],dt=0.125_dp)
    call check(tf%status == dlti_success,'transfer constructor')
    call check(tf%representation == dlti_transfer_function,'transfer representation')
    call check(tf%ninputs() == 1 .and. tf%noutputs() == 1,'transfer shape')
    call check(tf%nstates() == 1,'transfer state count')
    call check(tf%stable(),'transfer stability')
    ss = tf%to_ss()
    call check(ss%status == dlti_success .and. ss%representation == dlti_state_space,'to ss')
    zp = tf%to_zpk()
    call check(zp%status == dlti_success .and. zp%representation == dlti_zeros_poles_gain,'to zpk')
    converted = ss%to_tf()
    call check(converted%status == dlti_success,'ss to tf')
    call check(maxval(abs(converted%numerator-[cmplx(0.3_dp,0.0_dp,dp), &
        cmplx(0.1_dp,0.0_dp,dp)])) < 1.0e-12_dp,'ss to tf numerator')
    converted = zp%to_tf()
    call check(converted%status == dlti_success,'zpk to tf')
    call check(maxval(abs(converted%denominator-tf%denominator)) < 1.0e-11_dp,'zpk tf denominator')
    do i = 1, size(u)
        u(i) = sin(0.3_dp*i)
        t(i) = real(i-1,dp)*tf%dt
    end do
    observed = tf%output(u)
    reference = dlsim(tf%numerator,tf%denominator,cmplx(u,0.0_dp,dp),tf%dt)
    call check(observed%status == discrete_sim_success,'output status')
    call check(maxval(abs(observed%y-reference%y)) < 1.0e-12_dp,'output samples')
    observed = ss%output(u,t=t)
    call check(maxval(abs(observed%y-reference%y)) < 1.0e-12_dp,'ss output')
    observed = zp%output(cmplx(u,0.0_dp,dp))
    call check(maxval(abs(observed%y-reference%y)) < 1.0e-12_dp,'zpk output')
    response = tf%step(n=20)
    expected = dstep(tf%numerator,tf%denominator,tf%dt,n=20)
    call check(response%status == discrete_sim_success,'step status')
    call check(maxval(abs(response%y-expected%y)) < 1.0e-12_dp,'step samples')
    response = zp%impulse(n=20)
    expected = dimpulse(tf%numerator,tf%denominator,tf%dt,n=20)
    call check(maxval(abs(response%y-expected%y)) < 1.0e-12_dp,'impulse samples')
    w = [0.05_dp,0.3_dp,0.6_dp,1.1_dp,2.5_dp]
    frequency = tf%freqresp(w=w)
    freq_reference = dfreqresp(tf%numerator,tf%denominator,w=w,dt=tf%dt)
    call check(frequency%status == 0,'frequency status')
    call check(maxval(abs(frequency%response-freq_reference%response)) < 1.0e-12_dp,'frequency')

    a = cmplx(0.0_dp,0.0_dp,dp)
    a(1,1) = cmplx(0.7_dp,0.05_dp,dp)
    a(2,2) = cmplx(-0.1_dp,0.1_dp,dp)
    b(1,:) = [cmplx(0.2_dp,0.1_dp,dp),cmplx(0.5_dp,0.0_dp,dp)]
    b(2,:) = [cmplx(-0.3_dp,0.0_dp,dp),cmplx(0.4_dp,0.1_dp,dp)]
    c(1,:) = [cmplx(1.0_dp,0.0_dp,dp),cmplx(0.1_dp,0.0_dp,dp)]
    c(2,:) = [cmplx(0.0_dp,0.2_dp,dp),cmplx(0.9_dp,0.0_dp,dp)]
    d = cmplx(0.0_dp,0.0_dp,dp)
    mimo = dlti_ss(a,b,c,d,dt=0.25_dp)
    call check(mimo%status == dlti_success,'MIMO constructor')
    call check(mimo%ninputs() == 2 .and. mimo%noutputs() == 2,'MIMO shape')
    call check(mimo%stable(),'MIMO internal stability')
    converted = mimo%to_tf()
    call check(converted%status == dlti_invalid_input,'reject MIMO to_tf')
    do i = 1,8
        ur(i,:) = [sin(0.2_dp*i),cos(0.15_dp*i)]
    end do
    ucomplex = cmplx(ur,0.0_dp,dp)
    observed = mimo%output(ucomplex)
    reference = dlsim(a,b,c,d,ucomplex,0.25_dp)
    call check(maxval(abs(observed%y-reference%y)) < 1.0e-12_dp,'MIMO output')
    call check(maxval(abs(observed%x-reference%x)) < 1.0e-12_dp,'MIMO state history')
    cascade = mimo%series(mimo)
    call check(cascade%status == dlti_success .and. cascade%nstates() == 4,'MIMO cascade')
    upstream = mimo%output(ucomplex)
    downstream = mimo%output(upstream%y)
    observed = cascade%output(ucomplex)
    call check(maxval(abs(observed%y-downstream%y)) < 2.0e-12_dp,'MIMO cascade output')
    parallel_system = mimo%parallel(mimo)
    call check(parallel_system%status == dlti_success,'MIMO parallel')
    observed = parallel_system%output(ucomplex)
    call check(maxval(abs(observed%y-2.0_dp*upstream%y)) < 2.0e-12_dp,'MIMO parallel output')
    converted = mimo%series(tf)
    call check(converted%status == dlti_invalid_input,'cascade channel mismatch')
    converted = mimo%parallel(tf)
    call check(converted%status == dlti_invalid_input,'parallel channel mismatch')
    converted = mimo%parallel(dlti_ss(a,b,c,d,dt=2.0_dp))
    call check(converted%status == dlti_invalid_input,'parallel dt mismatch')
    converted = tf%series(zp)
    call check(converted%status == dlti_success,'SISO cascade')
    observed = converted%output(cmplx(u,0.0_dp,dp))
    upstream = zp%output(cmplx(u,0.0_dp,dp))
    downstream = tf%output(upstream%y(:,1))
    call check(maxval(abs(observed%y-downstream%y)) < 2.0e-12_dp,'SISO cascade output')
    response = mimo%step(n=8)
    call check(size(response%y,2) == 2 .and. size(response%y,3) == 2,'MIMO step shape')
    block
        use scifort_signal_lti_frequency, only : mimo_frequency_response
        type(mimo_frequency_response) :: mr, mf
        mr = mimo%freqresp_mimo(w=w,solver='direct')
        mf = dfreqresp_mimo(a,b,c,d,w=w,dt=0.25_dp,solver='direct')
        call check(mr%status == 0,'MIMO freq status')
        call check(maxval(abs(mr%response-mf%response)) < 1.0e-12_dp,'MIMO freq')
    end block
    a(1,1) = cmplx(1.05_dp,0.0_dp,dp)
    unstable = dlti_ss(a,b,c,d)
    call check(.not. unstable%stable(),'unstable MIMO eigenvalue')
    invalid = dlti_tf([1.0_dp],[0.0_dp,0.0_dp])
    call check(invalid%status == dlti_invalid_input,'zero denominator')
    invalid = dlti_tf([1.0_dp],[1.0_dp,0.1_dp],dt=-1.0_dp)
    call check(invalid%status == dlti_invalid_input,'bad dt')
    invalid = dlti_zpk([0.5_dp],[0.4_dp],1.0_dp,dt=0.0_dp)
    call check(invalid%status == dlti_invalid_input,'bad ZPK dt')
    invalid = dlti_ss(a,b,c,d(:,1:1))
    call check(invalid%status == dlti_invalid_input,'bad MIMO dimensions')
    print '(a)', 'PASS: dlti constructors, conversions, MIMO, simulation and analysis'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! passed assertion
        character(*), intent(in) :: label !! context
        if (.not. ok) then
            print '(a,1x,a)', 'FAIL:', label
            error stop 1
        end if
    end subroutine check
end program test_signal_dlti
