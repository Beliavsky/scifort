! SPDX-License-Identifier: MIT
! Run separately with 'ql' and 'partial' to compare CPU time and RSS.
program bench_signal_dpss_partial
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result
    implicit none
    character(16) :: solver, arg
    type(dpss_result) :: res
    real(dp) :: start, finish
    integer :: n, ios

    n = 2048
    solver = 'partial'
    call get_command_argument(1, solver)
    call get_command_argument(2, arg)
    if (len_trim(arg) /= 0) then
        read(arg, *, iostat=ios) n
        if (ios /= 0 .or. n < 32) error stop 'invalid window length'
    end if
    call cpu_time(start)
    res = dpss(n, 3.5_dp, kmax=6, solver=trim(solver))
    call cpu_time(finish)
    if (res%status /= 0) error stop 'DPSS benchmark failure'
    print '(a,1x,a,1x,a,i0,1x,a,f9.5,1x,a,es12.4)', &
        'solver',trim(solver),'n=',n,'cpu_seconds=',finish-start, &
        'checksum=',sum(res%ratios)
end program bench_signal_dpss_partial
