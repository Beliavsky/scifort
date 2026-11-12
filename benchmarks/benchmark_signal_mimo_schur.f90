! SPDX-License-Identifier: MIT
! Compare reusable Hessenberg and complex Schur MIMO frequency sweeps.
program benchmark_signal_mimo_schur
    use scifort_signal, only : dp, freqresp_mimo, prepare_mimo_frequency, &
        mimo_frequency_plan, mimo_frequency_response, design_status_success
    implicit none
    integer, parameter :: ns = 48, ni = 3, no = 2, nw = 1500, repeats = 4
    real(dp) :: a(ns,ns), b(ns,ni), c(no,ns), d(no,ni), w(nw)
    real(dp) :: t0, t1, hsetup, ssetup, htime, stime, difference
    type(mimo_frequency_plan) :: ph, ps
    type(mimo_frequency_response) :: rh, rs
    integer :: i, j, k

    do j = 1,ns
        do i = 1,ns
            a(i,j) = 0.015_dp*cos(real(i*7+j*3,dp))
            if (i == j) a(i,j) = a(i,j)-0.2_dp-0.05_dp*real(i,dp)
        end do
    end do
    do j = 1,ni
        do i = 1,ns
            b(i,j) = cos(real(i+3*j,dp))
        end do
    end do
    do j = 1,ns
        do i = 1,no
            c(i,j) = sin(real(2*i+j,dp))
        end do
    end do
    d = 0.05_dp
    w = [(3.0_dp*real(k-1,dp)/real(nw,dp),k=1,nw)]

    call cpu_time(t0)
    ph = prepare_mimo_frequency(a,b,c,d)
    call cpu_time(t1)
    hsetup = t1-t0
    call cpu_time(t0)
    ps = prepare_mimo_frequency(a,b,c,d,solver='schur')
    call cpu_time(t1)
    ssetup = t1-t0
    if (ph%status /= design_status_success) error stop 'Hessenberg setup failed'
    if (ps%status /= design_status_success) error stop 'Schur setup failed'
    if (ps%method /= 'schur') error stop 'Schur did not converge'

    call cpu_time(t0)
    do k = 1,repeats
        rh = freqresp_mimo(ph,w=w)
    end do
    call cpu_time(t1)
    htime = t1-t0
    call cpu_time(t0)
    do k = 1,repeats
        rs = freqresp_mimo(ps,w=w)
    end do
    call cpu_time(t1)
    stime = t1-t0
    difference = maxval(abs(rs%response-rh%response)/(1.0_dp+abs(rh%response)))
    if (difference > 1.0e-9_dp) error stop 'Schur and Hessenberg disagree'
    print '(a,i0,a,i0,a,i0)', 'states=',ns,' frequencies=',nw,' sweeps=',repeats
    print '(a,f11.6)', 'Hessenberg setup seconds: ',hsetup
    print '(a,f11.6)', 'Schur setup seconds: ',ssetup
    print '(a,f11.6)', 'Hessenberg sweeps seconds: ',htime
    print '(a,f11.6)', 'Schur sweeps seconds: ',stime
    print '(a,f8.3)', 'Sweep speedup (Schur vs Hessenberg): ',htime/stime
    print '(a,es12.4)', 'Maximum scaled response difference: ',difference
end program benchmark_signal_mimo_schur
