! SPDX-License-Identifier: MIT
! End-to-end cost of direct and FFT spectral-concentration calculations.
! Build with: gfortran -std=f2018 -O3 src/kinds.f90 src/signal_dpss.f90 &
!             benchmarks/bench_signal_dpss_fft_ratios.f90 -o bench_ratios
program bench_signal_dpss_fft_ratios
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result, dpss_status_success
    implicit none
    integer, parameter :: lengths(*) = [512,1024,2048,4096,8192]
    integer :: j, repeat, method
    real(dp) :: start, elapsed, best(2), diff, checksum
    type(dpss_result) :: direct, fourier

    print '(a)', 'n direct_s fft_s speedup max_ratio_diff'
    checksum = 0.0_dp
    do j = 1,size(lengths)
        do method = 1,2
            best(method) = huge(1.0_dp)
            do repeat = 1,2
                call cpu_time(start)
                if (method == 1) then
                    direct = dpss(lengths(j),3.5_dp,kmax=6, &
                        solver='partial',ratio_method='direct')
                    if (direct%status /= dpss_status_success) error stop 'direct failed'
                else
                    fourier = dpss(lengths(j),3.5_dp,kmax=6, &
                        solver='partial',ratio_method='fft')
                    if (fourier%status /= dpss_status_success) error stop 'fft failed'
                end if
                call cpu_time(elapsed)
                best(method) = min(best(method),elapsed-start)
            end do
        end do
        diff = maxval(abs(direct%ratios-fourier%ratios))
        if (diff > 1.0e-10_dp) error stop 'ratio mismatch'
        checksum = checksum+sum(fourier%ratios)
        print '(i7,2(1x,f12.6),1x,f9.2,1x,es12.4)', lengths(j), &
            best(1),best(2),best(1)/max(best(2),tiny(1.0_dp)),diff
    end do
    print '(a,es14.6)', 'checksum: ',checksum
end program bench_signal_dpss_fft_ratios
