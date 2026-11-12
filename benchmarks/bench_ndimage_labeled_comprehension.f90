! SPDX-License-Identifier: MIT
! Compare sorted-label bucket grouping with independent repeated image scans.
program bench_ndimage_labeled_comprehension
    use scifort_kinds, only : dp
    use scifort_ndimage_measurements, only : labeled_comprehension_nd
    implicit none
    integer, parameter :: n1=320,n2=250,nrequests=1200,trials=3
    real(dp) :: x(n1,n2),got(nrequests),expected(nrequests)
    integer :: labels(n1,n2),ids(nrequests),i,j,k,iter,code
    real(dp) :: t0,t1,bucket_time,scan_time,delta
    do j=1,n2
        do i=1,n1
            x(i,j)=sin(real(i+j,dp)*0.017_dp)
            labels(i,j)=modulo(37*i+103*j,1601)-800
        end do
    end do
    do k=1,nrequests
        ids(k)=modulo(17*k,1400)-700
    end do
    ids(2)=ids(1)
    call labeled_comprehension_nd(x,labels,ids,reduce,0.0_dp,got,code)
    if (code/=0) error stop 'indexed reduction failed'
    call direct_scan()
    delta=maxval(abs(got-expected))
    if (delta>1.0e-11_dp) error stop 'inconsistent reductions'

    call cpu_time(t0)
    do iter=1,trials
        call labeled_comprehension_nd(x,labels,ids,reduce,0.0_dp,got,code)
        if (code/=0) error stop 'indexed reduction failed'
    end do
    call cpu_time(t1)
    bucket_time=(t1-t0)/real(trials,dp)
    call cpu_time(t0)
    do iter=1,trials
        call direct_scan()
    end do
    call cpu_time(t1)
    scan_time=(t1-t0)/real(trials,dp)
    print '(a,f10.6)', 'Indexed bucket (s): ',bucket_time
    print '(a,f10.6)', 'Repeated scans (s): ',scan_time
    if (bucket_time>0.0_dp) then
        print '(a,f9.2)', 'Relative speedup: ',scan_time/bucket_time
    end if
    print '(a,es12.4)', 'Maximum difference: ',delta
contains
    function reduce(values) result(v)
        real(dp), intent(in) :: values(:) !! all selected samples
        real(dp) :: v !! region sum
        v=sum(values)
    end function reduce

    subroutine direct_scan()
        integer :: i,j,k
        expected=0.0_dp
        do k=1,nrequests
            do j=1,n2
                do i=1,n1
                    if (labels(i,j)==ids(k)) expected(k)=expected(k)+x(i,j)
                end do
            end do
        end do
    end subroutine direct_scan
end program bench_ndimage_labeled_comprehension
