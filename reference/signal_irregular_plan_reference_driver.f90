! SPDX-License-Identifier: MIT
! Deterministic repeated-grid Lomb-Scargle plan reference driver.
program signal_irregular_plan_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : irregular_status_success
    use scifort_signal_irregular_nd, only : lombscargle_plan, &
        lombscargle_plan_into, lombscargle_plan_amplitude_into
    implicit none
    type(lombscargle_plan) :: plan
    real(dp), allocatable :: x(:), w(:), f(:), y(:,:,:), power(:,:,:)
    complex(dp), allocatable :: amp(:,:,:)
    integer :: n, nf, nc, variant, ios, st, i, j, channel, k
    logical :: floating, center, unweighted, amplitude, normalized

    do
        read(*,*,iostat=ios) n,nf,nc,variant
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid case header'
        allocate(x(n),w(n),f(nf),y(2,n,nc),power(2,nf,nc),amp(2,nf,nc))
        do i=1,n
            x(i) = 0.17_dp*(i-1)+0.009_dp*i*i+0.07_dp*sin(0.21_dp*i)
            w(i) = 0.17_dp+0.14_dp*mod(i*5,11)
        end do
        do j=1,nf
            f(j) = 0.15_dp+0.19_dp*j+0.003_dp*j*j
        end do
        do channel=1,nc
            do k=1,2
                do i=1,n
                    y(k,i,channel)=0.27_dp*cos(0.7_dp*x(i)+0.14_dp*k*channel) &
                        +sin(1.8_dp*x(i)-0.3_dp*channel)+0.05_dp*i/n+0.27_dp*k
                end do
            end do
        end do
        floating = (variant == 1 .or. variant == 2 .or. variant == 6)
        center = (variant == 4 .or. variant == 5 .or. variant == 7)
        unweighted = (variant == 3 .or. variant == 7)
        amplitude = (variant == 2 .or. variant == 5 .or. variant == 7)
        normalized = (variant == 1 .or. variant == 6)
        if (variant == 6) then
            do i=1,n
                if (mod(i,4)==0) w(i)=0.0_dp
            end do
        end if
        if (unweighted) then
            call plan%initialize(x,f,floating_mean=floating,status=st)
        else
            call plan%initialize(x,f,weights=w,floating_mean=floating,status=st)
        end if
        if (st /= irregular_status_success) error stop 'plan initialization'
        if (amplitude) then
            call lombscargle_plan_amplitude_into(plan,y,amp,axis=2, &
                precenter=center,status=st)
        else if (normalized) then
            call lombscargle_plan_into(plan,y,power,axis=2, &
                precenter=center,normalize='normalize',status=st)
        else
            call lombscargle_plan_into(plan,y,power,axis=2, &
                precenter=center,status=st)
        end if
        if (st /= irregular_status_success) error stop 'planned transform'
        write(*,'(a,1x,5(i0,1x))') 'CASE',n,nf,nc,variant,st
        do channel=1,nc
            do j=1,nf
                do k=1,2
                    if (amplitude) then
                        write(*,'(2(es27.17e3,1x))') real(amp(k,j,channel),dp), &
                            aimag(amp(k,j,channel))
                    else
                        write(*,'(es27.17e3)') power(k,j,channel)
                    end if
                end do
            end do
        end do
        deallocate(x,w,f,y,power,amp)
    end do
end program signal_irregular_plan_reference_driver
