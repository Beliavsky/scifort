! SPDX-License-Identifier: MIT
! Compare indexed label-count lookup with a direct per-label scan.
program bench_ndimage_label_lookup
    use scifort_ndimage, only : count_labels_nd
    implicit none
    integer, parameter :: n=100000,k=1200
    integer :: labels(n),requested(k),fast(k),slow(k)
    integer :: i,j,status
    real :: start_time,end_time,fast_seconds,slow_seconds

    do i=1,n
        labels(i)=modulo(197*i+73*(i/13),1009)-504
    end do
    do j=1,k
        requested(j)=modulo(83*j,1009)-504
    end do
    call cpu_time(start_time)
    call count_labels_nd(labels,requested,fast,status)
    call cpu_time(end_time)
    if (status/=0) error stop 'indexed count_labels_nd failed'
    fast_seconds=end_time-start_time

    slow=0
    call cpu_time(start_time)
    do i=1,n
        do j=1,k
            if (labels(i)==requested(j)) slow(j)=slow(j)+1
        end do
    end do
    call cpu_time(end_time)
    slow_seconds=end_time-start_time
    if (any(fast/=slow)) error stop 'indexed/direct label-count disagreement'
    print '(a,i0,a,i0)', 'pixels=',n,', requested IDs=',k
    print '(a,f9.5)', 'indexed seconds=',fast_seconds
    print '(a,f9.5)', 'naive seconds  =',slow_seconds
    if (fast_seconds>0.0) print '(a,f9.2)', 'speedup=',slow_seconds/fast_seconds
end program bench_ndimage_label_lookup
