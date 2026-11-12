! SPDX-License-Identifier: MIT
! Compare original dense shifted solves with reusable Hessenberg reduction.
program benchmark_signal_mimo_hessenberg
    use scifort_signal, only : dp, freqresp_mimo, mimo_frequency_response, &
        design_status_success
    implicit none
    integer, parameter :: ns = 80, ni = 3, no = 2, nf = 600
    real(dp) :: a(ns,ns), b(ns,ni), c(no,ns), d(no,ni), w(nf)
    real(dp) :: t0, t1, t2, t3, scaled
    type(mimo_frequency_response) :: dense, reduced
    integer :: i,j

    do j = 1,ns
        do i = 1,ns
            a(i,j) = 0.006_dp*sin(real(i+3*j,dp))
            if (i == j) a(i,j) = a(i,j)-0.3_dp-0.015_dp*real(i,dp)
        end do
    end do
    do j = 1,ni
        do i = 1,ns
            b(i,j) = sin(real(i+2*j,dp))
        end do
    end do
    do j = 1,ns
        do i = 1,no
            c(i,j) = cos(real(i+2*j,dp))
        end do
    end do
    d = 0.01_dp
    do i = 1,nf
        w(i) = real(i-1,dp)*20.0_dp/real(nf-1,dp)
    end do

    call cpu_time(t0)
    dense = freqresp_mimo(a,b,c,d,w=w,solver='direct')
    call cpu_time(t1)
    reduced = freqresp_mimo(a,b,c,d,w=w,solver='hessenberg')
    call cpu_time(t2)
    if (dense%status /= design_status_success .or. &
        reduced%status /= design_status_success) error stop 'solver failure'
    scaled = maxval(abs(reduced%response-dense%response)/(1.0_dp+abs(dense%response)))
    if (scaled > 2.0e-11_dp) error stop 'Hessenberg response mismatch'
    print '(a,i0,a,i0)', 'states=',ns,' frequencies=',nf
    print '(a,f10.5)', 'direct CPU seconds:     ', t1-t0
    print '(a,f10.5)', 'Hessenberg CPU seconds: ', t2-t1
    print '(a,f10.4)', 'speedup:                ', (t1-t0)/max(1.0e-9_dp,t2-t1)
    print '(a,es12.4)', 'max scaled difference:  ', scaled
end program benchmark_signal_mimo_hessenberg
