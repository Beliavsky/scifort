! SPDX-License-Identifier: MIT
! Compare blocked overlap-add implementations, built separately against each.
program bench_signal_oaconvolve_tiling
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : oaconvolve_nd
    implicit none
    real(dp), allocatable :: a(:,:), b(:,:), y(:,:)
    real(dp) :: t0, t1, checks
    integer :: rep, status, i, j, n1, n2

    do j=1,2
        select case(j)
        case(1)
            n1=513
            n2=257
        case(2)
            n1=1025
            n2=513
        end select
        allocate(a(n1,n2),b(5,5),y(n1+4,n2+4))
        a = reshape([(sin(0.01_dp*i)+cos(0.013_dp*i),i=1,size(a))],shape(a))
        b = reshape([(cos(0.17_dp*i),i=1,size(b))],shape(b))
        call cpu_time(t0)
        do rep=1,2
            call oaconvolve_nd(a,b,y,status=status)
            if (status /= 0) error stop 'oaconvolve error'
        end do
        call cpu_time(t1)
        checks = sum(abs(y))
        print '(a,i0,a,i0,a,es13.5,a,f9.4)', 'shape=',n1,'x',n2, &
            ' checksum=',checks,' elapsed_s=',t1-t0
        deallocate(a,b,y)
    end do
end program bench_signal_oaconvolve_tiling
