! SPDX-License-Identifier: MIT
! Deterministic driver for SciPy comparisons of multichannel Lomb-Scargle.
program signal_irregular_nd_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_irregular_nd, only : lombscargle_into, lombscargle_amplitude_into
    implicit none
    integer :: n, nf, nc, variant, ios, st, i, j, k, channel
    real(dp), allocatable :: x(:), freqs(:), w(:), y(:,:,:), power(:,:,:)
    complex(dp), allocatable :: amp(:,:,:)

    do
        read(*,*,iostat=ios) n,nf,nc,variant
        if (ios < 0) exit
        if (ios /= 0) error stop 'reference driver input'
        if (n < 3 .or. nf < 1 .or. nc < 1) error stop 'reference driver dimensions'
        allocate(x(n),freqs(nf),w(n),y(2,n,nc),power(2,nf,nc),amp(2,nf,nc))
        do i = 1,n
            x(i) = 0.37_dp*real(i-1,dp) + 0.009_dp*real(i*i,dp) + 0.13_dp*sin(0.17_dp*i)
            w(i) = 0.2_dp + 0.15_dp*real(mod(i*7,13),dp)
        end do
        do j = 1,nf
            freqs(j) = 0.13_dp + 0.31_dp*j + 0.004_dp*j*j
        end do
        do channel = 1,nc
            do i = 1,n
                do k = 1,2
                    y(k,i,channel) = cos(0.85_dp*x(i)+0.23_dp*k*channel) + &
                        0.37_dp*sin(1.74_dp*x(i)-0.13_dp*channel) + &
                        0.05_dp*cos(0.19_dp*i*k) + 0.4_dp*channel + 0.13_dp*k
                end do
            end do
        end do
        select case(variant)
        case (0)
            call lombscargle_into(x,y,freqs,power,axis=2,status=st)
        case (1)
            call lombscargle_into(x,y,freqs,power,axis=2,weights=w, &
                floating_mean=.true.,normalize='normalize',status=st)
        case (2)
            call lombscargle_amplitude_into(x,y,freqs,amp,axis=2,weights=w, &
                floating_mean=.true.,status=st)
        case (3)
            call lombscargle_into(x,y,freqs,power,axis=2,weights=w, &
                precenter=.true.,status=st)
        case (4)
            call lombscargle_into(x,y,freqs,power,axis=2,weights=w, &
                floating_mean=.true.,status=st)
        case (5)
            call lombscargle_amplitude_into(x,y,freqs,amp,axis=2,status=st)
        case default
            error stop 'unknown reference variant'
        end select
        write(*,'(a,1x,5(i0,1x))') 'CASE',n,nf,nc,variant,st
        do channel = 1,nc
            do j = 1,nf
                do k = 1,2
                    if (variant == 2 .or. variant == 5) then
                        write(*,'(2(es27.17e3,1x))') real(amp(k,j,channel),dp), &
                            aimag(amp(k,j,channel))
                    else
                        write(*,'(es27.17e3)') power(k,j,channel)
                    end if
                end do
            end do
        end do
        deallocate(x,freqs,w,y,power,amp)
    end do
end program signal_irregular_nd_reference_driver
