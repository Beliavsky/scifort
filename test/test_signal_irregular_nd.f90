! SPDX-License-Identifier: MIT
! Independent native regressions for multichannel irregularly sampled spectra.
program test_signal_irregular_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : lombscargle, lombscargle_amplitude, &
        irregular_status_success, irregular_status_invalid_input, irregular_status_singular, &
        lombscargle_into, lombscargle_amplitude_into
    implicit none
    integer, parameter :: nt = 11, nf = 7
    real(dp) :: times(nt), omega(nf), w(nt), y(nt,3), p(nf,3)
    real(dp) :: cube(2,nt,3), pcube(2,nf,3)
    real(dp) :: high(2,1,1,1,1,1,nt,1,1,1,1,1,1,1,2)
    real(dp) :: hp(2,1,1,1,1,1,nf,1,1,1,1,1,1,1,2)
    real(dp) :: parent(2*nt,4), output(2*nf,4), pbad(nf,3)
    complex(dp) :: amps(nf,3), acube(2,nf,3), ahigh(2,1,1,1,1,1,nf,1,1,1,1,1,1,1,2)
    complex(dp) :: abad(nf,3)
    real(dp), allocatable :: reference(:)
    complex(dp), allocatable :: areference(:)
    integer :: i, j, k, st

    do i = 1,nt
        times(i) = 0.32_dp * real(i-1,dp) + 0.018_dp * real(i*i,dp)
        w(i) = 0.5_dp + 0.13_dp * real(mod(i*7,11),dp)
    end do
    do i = 1,nf
        omega(i) = 0.39_dp + 0.23_dp * real(i,dp)
    end do
    do j = 1,3
        do i = 1,nt
            y(i,j) = cos(0.45_dp*times(i)+0.23_dp*j) + &
                0.31_dp*sin(1.7_dp*times(i)-0.27_dp*j) + 0.4_dp*j
        end do
    end do

    call lombscargle_into(times,y,omega,p,axis=1,weights=w, &
        floating_mean=.true.,normalize='normalize',status=st)
    if (st /= irregular_status_success) error stop 'rank-two power status'
    call lombscargle_amplitude_into(times,y,omega,amps,axis=1,weights=w, &
        floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'rank-two amplitude status'
    do j = 1,3
        reference = lombscargle(times,y(:,j),omega,normalize='normalize', &
            weights=w,floating_mean=.true.)
        areference = lombscargle_amplitude(times,y(:,j),omega,weights=w,floating_mean=.true.)
        call check(p(:,j),reference)
        if (maxval(abs(amps(:,j)-areference)) > 2.0e-13_dp) error stop 'rank-two amplitude mismatch'
    end do

    do k = 1,3
        do j = 1,2
            cube(j,:,k) = y(:,k) + 0.2_dp*real(j,dp)*sin(0.8_dp*times)
        end do
    end do
    call lombscargle_into(times,cube,omega,pcube,axis=2,weights=w, &
        floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'rank-three power status'
    call lombscargle_amplitude_into(times,cube,omega,acube,axis=-2,weights=w, &
        floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'negative-axis amplitude status'
    do k = 1,3
        do j = 1,2
            call check(pcube(j,:,k),lombscargle(times,cube(j,:,k),omega, &
                weights=w,floating_mean=.true.))
            areference = lombscargle_amplitude(times,cube(j,:,k),omega, &
                weights=w,floating_mean=.true.)
            if (maxval(abs(acube(j,:,k)-areference)) > 2.0e-13_dp) error stop 'cube amplitude'
        end do
    end do

    do k = 1,2
        do j = 1,2
            high(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k) = y(:,k) + 0.31_dp*j
        end do
    end do
    call lombscargle_into(times,high,omega,hp,axis=7,weights=w, &
        floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'rank-fifteen power status'
    call lombscargle_amplitude_into(times,high,omega,ahigh,axis=-9,weights=w, &
        floating_mean=.true.,status=st)
    if (st /= irregular_status_success) error stop 'rank-fifteen amplitude status'
    do k = 1,2
        do j = 1,2
            call check(hp(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k), &
                lombscargle(times,high(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k),omega, &
                weights=w,floating_mean=.true.))
            areference = lombscargle_amplitude(times,high(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k),omega, &
                weights=w,floating_mean=.true.)
            if (maxval(abs(ahigh(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k)-areference)) > 2.0e-13_dp) &
                error stop 'rank-fifteen amplitude'
        end do
    end do

    parent = -1.0_dp
    output = 817.0_dp
    parent(1:2*nt:2,1) = y(:,1)
    parent(1:2*nt:2,3) = y(:,2)
    call lombscargle_into(times,parent(1:2*nt:2,1:4:2),omega, &
        output(1:2*nf:2,1:4:2),axis=1,weights=w,status=st)
    if (st /= irregular_status_success) error stop 'strided output status'
    call check(output(1:2*nf:2,1),lombscargle(times,y(:,1),omega,weights=w))
    call check(output(1:2*nf:2,3),lombscargle(times,y(:,2),omega,weights=w))
    if (any(output(2:2*nf:2,:) /= 817.0_dp) .or. &
        any(output(:,2) /= 817.0_dp)) error stop 'strided output touched unrelated elements'

    pbad = 817.0_dp
    call lombscargle_into(times,y,omega,pbad,axis=1,weights=-w,status=st)
    call unchanged_real(st,irregular_status_invalid_input)
    call lombscargle_into(times,y,omega,pbad,axis=1,weights=w,normalize='broken',status=st)
    call unchanged_real(st,irregular_status_invalid_input)
    call lombscargle_into(times,y,omega,pbad,axis=3,weights=w,status=st)
    call unchanged_real(st,irregular_status_invalid_input)
    call lombscargle_into(times,y,omega,pbad(:6,:),axis=1,status=st)
    call unchanged_real(st,irregular_status_invalid_input)
    call lombscargle_into(times,y,omega,pbad,axis=1,weights=w(:nt-1),status=st)
    call unchanged_real(st,irregular_status_invalid_input)
    abad = cmplx(817.0_dp,-28.0_dp,dp)
    call lombscargle_amplitude_into(times,y,omega,abad,axis=1,weights=w(:nt-1),status=st)
    if (st /= irregular_status_invalid_input .or. &
        any(abad /= cmplx(817.0_dp,-28.0_dp,dp))) error stop 'amplitude invalid preservation'

    ! A singular last channel must not commit successful earlier channels.
    pbad = 817.0_dp
    cube(1,:,3) = 5.0_dp
    call lombscargle_into(times,cube,omega,pcube,axis=2, &
        floating_mean=.true.,normalize='normalize',status=st)
    if (st /= irregular_status_singular) error stop 'expected singular channel'
    call lombscargle_into(times,reshape([y(:,1),y(:,2),spread(2.0_dp,1,nt)], &
        [nt,3]),omega,pbad,axis=1,floating_mean=.true.,status=st)
    if (st /= irregular_status_singular .or. any(pbad /= 817.0_dp)) &
        error stop 'atomic power commit on singular channel'
    call lombscargle_amplitude_into(times,reshape([y(:,1),y(:,2),spread(2.0_dp,1,nt)], &
        [nt,3]),omega,abad,axis=1,floating_mean=.true.,status=st)
    if (st /= irregular_status_singular .or. any(abad /= cmplx(817.0_dp,-28.0_dp,dp))) &
        error stop 'atomic amplitude commit on singular channel'

    print '(a)', 'test_signal_irregular_nd passed'
contains
    subroutine check(actual, expected)
        real(dp), intent(in) :: actual(:) !! candidate power series
        real(dp), intent(in) :: expected(:) !! independent scalar power series
        if (size(actual) /= size(expected)) error stop 'shape mismatch'
        if (any(abs(actual-expected) > 2.0e-13_dp)) error stop 'numerical mismatch'
    end subroutine check

    subroutine unchanged_real(code,expected)
        integer, intent(in) :: code !! reported error status
        integer, intent(in) :: expected !! expected error status
        if (code /= expected .or. any(pbad /= 817.0_dp)) &
            error stop 'error handling changed caller output'
    end subroutine unchanged_real
end program test_signal_irregular_nd
