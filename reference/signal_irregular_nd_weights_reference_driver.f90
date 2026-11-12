! SPDX-License-Identifier: MIT
! Deterministic multichannel driver for SciPy weighted/omission comparisons.
program signal_irregular_nd_weights_reference_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_irregular_nd, only : lombscargle_into, lombscargle_amplitude_into
    implicit none
    integer :: n, nf, nc, variant, ios, st, i, j, k, channel
    real(dp), allocatable :: x(:), freqs(:), w(:), wn(:,:,:), y(:,:,:), power(:,:,:)
    complex(dp), allocatable :: amp(:,:,:)

    do
        read(*,*,iostat=ios) n,nf,nc,variant
        if (ios<0) exit
        if (ios/=0) error stop 'input read'
        allocate(x(n),freqs(nf),w(n),wn(2,n,nc),y(2,n,nc),power(2,nf,nc),amp(2,nf,nc))
        do i=1,n
            x(i)=0.31_dp*(i-1)+0.006_dp*i*i+0.17_dp*sin(0.23_dp*i)
            w(i)=0.20_dp+0.11_dp*mod(i*7,11)
        end do
        do j=1,nf
            freqs(j)=0.27_dp+0.31_dp*j+0.0007_dp*j*j
        end do
        do channel=1,nc
            do i=1,n
                do k=1,2
                    y(k,i,channel)=cos(0.85_dp*x(i)+0.23_dp*k*channel) &
                        +0.37_dp*sin(1.74_dp*x(i)-0.13_dp*channel) &
                        +0.05_dp*cos(0.19_dp*i*k)+0.4_dp*channel+0.13_dp*k
                    wn(k,i,channel)=0.18_dp+0.10_dp*mod(i*5+k*3+channel*7,13)
                end do
            end do
        end do
        if (variant==6 .or. variant==7 .or. variant==8 .or. variant==9) then
            do channel=1,nc
                do i=1,n
                    do k=1,2
                        if (mod(i+3*k+channel,7)==0) y(k,i,channel)= &
                            ieee_value(0.0_dp,ieee_quiet_nan)
                        if (mod(i+2*k+channel,6)==0) wn(k,i,channel)=0.0_dp
                    end do
                end do
            end do
        end if
        if (variant==8 .or. variant==9) then
            do i=1,n
                if (mod(i,5)==0) w(i)=0.0_dp
            end do
        end if
        select case(variant)
        case(0)
            call lombscargle_into(x,y,freqs,power,axis=2,weights=w,status=st)
        case(1)
            call lombscargle_into(x,y,freqs,power,axis=2,weights=w, &
                floating_mean=.true.,normalize='normalize',status=st)
        case(2)
            call lombscargle_amplitude_into(x,y,freqs,amp,axis=2,weights=w, &
                floating_mean=.true.,status=st)
        case(3)
            call lombscargle_into(x,y,freqs,power,axis=2,weights_nd=wn,status=st)
        case(4)
            call lombscargle_into(x,y,freqs,power,axis=2,weights_nd=wn, &
                floating_mean=.true.,normalize='normalize',status=st)
        case(5)
            call lombscargle_amplitude_into(x,y,freqs,amp,axis=2,weights_nd=wn, &
                floating_mean=.true.,status=st)
        case(6)
            call lombscargle_into(x,y,freqs,power,axis=2,weights_nd=wn, &
                nan_policy='omit',status=st)
        case(7)
            call lombscargle_into(x,y,freqs,power,axis=2,weights_nd=wn, &
                nan_policy='omit',normalize='normalize',floating_mean=.true.,status=st)
        case(8)
            call lombscargle_into(x,y,freqs,power,axis=2,weights=w, &
                nan_policy='omit',floating_mean=.true.,status=st)
        case(9)
            call lombscargle_amplitude_into(x,y,freqs,amp,axis=2,weights=w, &
                nan_policy='omit',floating_mean=.true.,status=st)
        case(10)
            call lombscargle_into(x,y,freqs,power,axis=2,weights_nd=wn, &
                precenter=.true.,status=st)
        case(11)
            call lombscargle_amplitude_into(x,y,freqs,amp,axis=2,weights=w, &
                precenter=.true.,status=st)
        case(12)
            call lombscargle_into(x,y,freqs,power,axis=2,weights_nd=wn, &
                precenter=.true.,floating_mean=.true.,normalize='normalize',status=st)
        case(13)
            call lombscargle_into(x,y,freqs,power,axis=2,status=st)
        case default
            error stop 'unknown variant'
        end select
        write(*,'(a,1x,5(i0,1x))') 'CASE',n,nf,nc,variant,st
        do channel=1,nc
            do j=1,nf
                do k=1,2
                    if(variant==2 .or. variant==5 .or. variant==9 .or. variant==11) then
                        write(*,'(2(es27.17e3,1x))') real(amp(k,j,channel),dp), &
                            aimag(amp(k,j,channel))
                    else
                        write(*,'(es27.17e3)') power(k,j,channel)
                    end if
                end do
            end do
        end do
        deallocate(x,freqs,w,wn,y,power,amp)
    end do
end program signal_irregular_nd_weights_reference_driver
