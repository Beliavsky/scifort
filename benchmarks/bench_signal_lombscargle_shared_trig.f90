! SPDX-License-Identifier: MIT
! Bench shared-frequency trigonometry versus independent channel evaluation.
program bench_signal_lombscargle_shared_trig
    use scifort_kinds, only : dp
    use scifort_signal_irregular_nd, only : lombscargle_into
    implicit none
    real(dp), allocatable :: times(:), y(:,:), omega(:), weights(:), power(:,:)
    real(dp) :: t0, t1, checksum
    integer :: n, nf, nc, repeats, i, j, k, st

    read(*,*) n, nf, nc, repeats
    allocate(times(n),y(n,nc),omega(nf),weights(n),power(nf,nc))
    do i=1,n
        times(i)=0.31_dp*(i-1)+0.0001_dp*i*i
        weights(i)=0.15_dp+0.12_dp*mod(i*13,11)
    end do
    do j=1,nf
        omega(j)=0.04_dp+0.016_dp*j
    end do
    do k=1,nc
        do i=1,n
            y(i,k)=cos((0.2_dp+0.003_dp*k)*times(i)) &
                +0.35_dp*sin(0.36_dp*times(i)+0.11_dp*k)+0.04_dp*k
        end do
    end do
    call cpu_time(t0)
    do i=1,repeats
        call lombscargle_into(times,y,omega,power,axis=1,weights=weights, &
            floating_mean=.true.,normalize='normalize',status=st)
        if(st/=0) error stop 'Lomb-Scargle benchmark status'
    end do
    call cpu_time(t1)
    checksum=sum(power)
    print '(a,1x,es18.10,1x,es18.10)', 'BENCH',t1-t0,checksum
end program bench_signal_lombscargle_shared_trig
