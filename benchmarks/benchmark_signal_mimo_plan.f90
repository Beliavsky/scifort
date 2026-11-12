! SPDX-License-Identifier: MIT
! Time separately reduced MIMO sweeps versus one reusable reduction.
program benchmark_signal_mimo_plan
    use scifort_signal, only : dp, freqresp_mimo, prepare_mimo_frequency, &
        mimo_frequency_plan, mimo_frequency_response, design_status_success
    implicit none
    integer, parameter :: ns = 80, ni = 6, no = 5, nw = 20, repeats = 25
    real(dp) :: a(ns,ns), b(ns,ni), c(no,ns), d(no,ni), w(nw)
    real(dp) :: t0, t1, dense_time, prepared_time, error
    type(mimo_frequency_plan) :: plan
    type(mimo_frequency_response) :: direct, prepared
    integer :: i,j,k

    do j = 1,ns
        do i = 1,ns
            a(i,j) = 0.01_dp*sin(real(i+3*j,dp))
            if (i == j) a(i,j) = a(i,j)-0.2_dp*real(i,dp)
        end do
    end do
    do j = 1,ni
        do i = 1,ns
            b(i,j) = cos(real(i+2*j,dp))
        end do
    end do
    do j = 1,ns
        do i = 1,no
            c(i,j) = sin(real(2*i+j,dp))
        end do
    end do
    d = 0.1_dp
    w = [(0.04_dp*real(k,dp),k=1,nw)]

    call cpu_time(t0)
    do k = 1,repeats
        direct = freqresp_mimo(a,b,c,d,w=w,solver='hessenberg', &
            inputs=[2,4],outputs=[1,5])
    end do
    call cpu_time(t1)
    dense_time = t1-t0
    call cpu_time(t0)
    plan = prepare_mimo_frequency(a,b,c,d)
    if (plan%status /= design_status_success) error stop 'Plan creation failed'
    do k = 1,repeats
        prepared = freqresp_mimo(plan,w=w,inputs=[2,4],outputs=[1,5])
    end do
    call cpu_time(t1)
    prepared_time = t1-t0
    if (direct%status /= design_status_success .or. &
        prepared%status /= design_status_success) error stop 'Response failure'
    error = maxval(abs(direct%response-prepared%response)/ &
        (1.0_dp+abs(direct%response)))
    print '(a,i0,a,i0,a,i0)', 'States: ', ns, ', frequencies/sweep: ', nw, ', repetitions: ', repeats
    print '(a,f11.6)', 'Repeated independent reduction CPU seconds: ', dense_time
    print '(a,f11.6)', 'Reusable preparation and sweeps CPU seconds: ', prepared_time
    print '(a,f10.3)', 'Speedup: ', dense_time/max(prepared_time,tiny(prepared_time))
    print '(a,es12.4)', 'Maximum scaled difference: ', error
    if (error > 1.0e-10_dp) error stop 'Inconsistent responses'
end program benchmark_signal_mimo_plan
