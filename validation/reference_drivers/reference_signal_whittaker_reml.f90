! SPDX-License-Identifier: MIT
! Numerical test driver for automatic Whittaker-Henderson REML selection.
program reference_signal_whittaker_reml
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_whittaker, only : whittaker_henderson, whittaker_result, &
        whittaker_status_success
    implicit none
    type(whittaker_result) :: res
    real(dp), allocatable :: y(:), w(:)
    integer :: ic, nc, n, p, mode, i

    open(unit=11,file='whittaker_reml_cases.txt',status='old',action='read')
    open(unit=12,file='whittaker_reml_results.txt',status='replace',action='write')
    read(11,*) nc
    do ic=1,nc
        read(11,*) n,p,mode
        allocate(y(n),w(n))
        do i=1,n
            y(i) = 0.03_dp*(i-1)+sin(0.12_dp*(i-1))+0.23_dp*cos(0.51_dp*(i-1)) + &
                0.2_dp*sin(2.1_dp*(i-1)+0.15_dp*ic)
            w(i) = 1.0_dp
            if (mode == 1 .and. mod(i+ic,7) == 0) w(i) = 0.0_dp
            if (mode == 2 .and. mod(i+ic,7) == 0) w(i) = 0.3_dp
            if (mode == 2 .and. mod(i+ic,11) == 0) w(i) = 0.0_dp
            if (w(i) == 0.0_dp) y(i) = ieee_value(0.0_dp,ieee_quiet_nan)
        end do
        res = whittaker_henderson(y,order=p,weights=w)
        if (res%status /= whittaker_status_success) then
            print *, 'Failed REML case:',ic,n,p,mode,res%status
            error stop 'REML solver failed'
        end if
        write(12,*) ic,n,p,mode
        write(12,'(es26.17)') res%lamb
        do i=1,n
            write(12,'(es26.17)') res%x(i)
        end do
        deallocate(y,w)
    end do
    close(11)
    close(12)
    print *, 'Whittaker-Henderson REML cases:',nc
end program reference_signal_whittaker_reml
