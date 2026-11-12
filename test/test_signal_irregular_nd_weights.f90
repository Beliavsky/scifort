! SPDX-License-Identifier: MIT
! Channel-specific Lomb-Scargle weighting, omitted samples and atomic writes.
program test_signal_irregular_nd_weights
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : lombscargle, lombscargle_amplitude, &
        irregular_status_success, irregular_status_invalid_input, irregular_status_singular, &
        lombscargle_into, lombscargle_amplitude_into
    implicit none
    integer, parameter :: nt=19, nf=9
    real(dp) :: t(nt), freq(nf), y(2,nt,3), weight(2,nt,3)
    real(dp) :: p(2,nf,3), pkeep(2,nf,3), wshared(nt)
    real(dp) :: high(2,1,1,1,1,1,nt,1,1,1,1,1,1,1,2)
    real(dp) :: highw(2,1,1,1,1,1,nt,1,1,1,1,1,1,1,2)
    real(dp) :: highp(2,1,1,1,1,1,nf,1,1,1,1,1,1,1,2)
    real(dp) :: long_data(2*nt,4), long_w(2*nt,4), long_p(2*nf,4)
    real(dp), allocatable :: ref(:)
    complex(dp) :: amp(2,nf,3), ampkeep(2,nf,3)
    complex(dp), allocatable :: aref(:)
    logical :: keep(nt)
    integer :: i, j, k, st

    do i=1,nt
        t(i)=0.28_dp*(i-1)+0.013_dp*i*i
        wshared(i)=0.2_dp+0.04_dp*mod(i*7,11)
    end do
    do i=1,nf
        freq(i)=0.18_dp+0.29_dp*i
    end do
    do k=1,3
        do j=1,2
            do i=1,nt
                y(j,i,k)=cos((0.8_dp+0.13_dp*j)*t(i)+0.21_dp*k) &
                    +0.18_dp*sin(2.12_dp*t(i)+0.11_dp*j*k)+0.35_dp*k
                weight(j,i,k)=0.1_dp+0.09_dp*mod(i*7+j*3+k*5,13)
            end do
        end do
    end do

    ! Distinct per-channel weights and the old shared-weight fast path.
    call lombscargle_into(t,y,freq,p,axis=2,weights_nd=weight, &
        floating_mean=.true.,normalize='normalize',status=st)
    if(st/=irregular_status_success) error stop 'channel-weight power status'
    call lombscargle_amplitude_into(t,y,freq,amp,axis=-2,weights_nd=weight, &
        floating_mean=.true.,status=st)
    if(st/=irregular_status_success) error stop 'channel-weight amplitude status'
    do k=1,3
        do j=1,2
            ref=lombscargle(t,y(j,:,k),freq,weights=weight(j,:,k), &
                floating_mean=.true.,normalize='normalize')
            call check_real(p(j,:,k),ref,5.0e-13_dp)
            aref=lombscargle_amplitude(t,y(j,:,k),freq,weights=weight(j,:,k), &
                floating_mean=.true.)
            call check_complex(amp(j,:,k),aref,5.0e-13_dp)
        end do
    end do
    call lombscargle_into(t,y,freq,p,axis=2,weights=wshared, &
        floating_mean=.true.,normalize='normalize',status=st)
    if(st/=irregular_status_success) error stop 'shared-weight fast path'
    do k=1,3
        do j=1,2
            call check_real(p(j,:,k),lombscargle(t,y(j,:,k),freq, &
                weights=wshared,floating_mean=.true.,normalize='normalize'),7.0e-13_dp)
        end do
    end do

    ! Omission excludes nonfinite observations and zero-weight samples in each lane.
    weight(1,3,1)=0.0_dp
    y(1,3,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    y(2,7,2)=ieee_value(0.0_dp,ieee_quiet_nan)
    weight(1,11,3)=0.0_dp
    call lombscargle_into(t,y,freq,p,axis=2,weights_nd=weight, &
        nan_policy='omit',floating_mean=.true.,normalize='normalize',status=st)
    if(st/=irregular_status_success) error stop 'omitted power status'
    call lombscargle_amplitude_into(t,y,freq,amp,axis=2,weights_nd=weight, &
        nan_policy='omit',floating_mean=.true.,status=st)
    if(st/=irregular_status_success) error stop 'omitted amplitude status'
    do k=1,3
        do j=1,2
            keep=(y(j,:,k)==y(j,:,k)).and.(weight(j,:,k)>0.0_dp)
            ref=lombscargle(pack(t,keep),pack(y(j,:,k),keep),freq, &
                weights=pack(weight(j,:,k),keep), &
                floating_mean=.true.,normalize='normalize')
            call check_real(p(j,:,k),ref,5.0e-13_dp)
            aref=lombscargle_amplitude(pack(t,keep),pack(y(j,:,k),keep), &
                freq,weights=pack(weight(j,:,k),keep),floating_mean=.true.)
            call check_complex(amp(j,:,k),aref,5.0e-13_dp)
        end do
    end do

    ! Test shared-weight omissions without changing the scalar reference API.
    call lombscargle_into(t,y,freq,p,axis=2,weights=wshared, &
        nan_policy='omit',status=st)
    if(st/=irregular_status_success) error stop 'shared omission status'
    do k=1,3
        do j=1,2
            keep=(y(j,:,k)==y(j,:,k)).and.(wshared>0.0_dp)
            call check_real(p(j,:,k),lombscargle(pack(t,keep), &
                pack(y(j,:,k),keep),freq,weights=pack(wshared,keep)),5.0e-13_dp)
        end do
    end do

    ! Rank 15 and negative axes preserve independent per-channel weighting.
    y(1,3,1)=0.0_dp
    y(2,7,2)=0.0_dp
    do k=1,2
        do j=1,2
            high(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k)=y(j,:,k)
            highw(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k)=weight(j,:,k)
        end do
    end do
    call lombscargle_into(t,high,freq,highp,axis=-9, &
        weights_nd=highw,normalize='normalize',status=st)
    if(st/=irregular_status_success) error stop 'rank 15 weights status'
    do k=1,2
        do j=1,2
            call check_real(highp(j,1,1,1,1,1,:,1,1,1,1,1,1,1,k), &
                lombscargle(t,y(j,:,k),freq,weights=weight(j,:,k), &
                normalize='normalize'),5.0e-13_dp)
        end do
    end do

    long_data=-1.0_dp
    long_w=0.0_dp
    long_p=456.0_dp
    long_data(1:2*nt:2,1)=y(1,:,1)
    long_data(1:2*nt:2,3)=y(2,:,2)
    long_w(1:2*nt:2,1)=weight(1,:,1)
    long_w(1:2*nt:2,3)=weight(2,:,2)
    call lombscargle_into(t,long_data(1:2*nt:2,1:4:2),freq, &
        long_p(1:2*nf:2,1:4:2),weights_nd=long_w(1:2*nt:2,1:4:2), &
        axis=1,status=st)
    if(st/=irregular_status_success) error stop 'strided per-channel weights'
    call check_real(long_p(1:2*nf:2,1), &
        lombscargle(t,y(1,:,1),freq,weights=weight(1,:,1)),5.0e-13_dp)
    if(any(long_p(2:2*nf:2,:)/=456.0_dp) .or. &
       any(long_p(:,2)/=456.0_dp)) error stop 'strided output overwritten'

    ! Invalid arguments and later-channel errors preserve all caller outputs.
    pkeep=456.0_dp
    call lombscargle_into(t,y,freq,pkeep,axis=2,weights_nd=weight, &
        weights=wshared,status=st)
    call reject(st)
    call lombscargle_into(t,y,freq,pkeep,axis=2,weights_nd=weight, &
        nan_policy='unknown',status=st)
    call reject(st)
    call lombscargle_into(t,y,freq,pkeep,axis=2,weights_nd=weight(:,1:nt-1,:),status=st)
    call reject(st)
    y(1,6,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call lombscargle_into(t,y,freq,pkeep,axis=2,weights_nd=weight,status=st)
    call reject(st)
    y(1,6,1)=0.0_dp
    weight(2,5,3)=-1.0_dp
    call lombscargle_into(t,y,freq,pkeep,axis=2,weights_nd=weight,status=st)
    call reject(st)
    weight(2,5,3)=0.0_dp
    ampkeep=cmplx(456.0_dp,-456.0_dp,dp)
    call lombscargle_amplitude_into(t,y,freq,ampkeep,axis=2, &
        weights_nd=weight,nan_policy='bad',status=st)
    if(st/=irregular_status_invalid_input .or. &
       any(ampkeep/=cmplx(456.0_dp,-456.0_dp,dp))) error stop 'invalid amplitude commit'
    ! A channel with zero positive weights is rejected after the previous channels.
    weight(2,:,3)=0.0_dp
    call lombscargle_into(t,y,freq,pkeep,axis=2,weights_nd=weight,status=st)
    call reject(st)
    call lombscargle_amplitude_into(t,y,freq,ampkeep,axis=2,weights_nd=weight,status=st)
    if(st/=irregular_status_invalid_input .or. &
       any(ampkeep/=cmplx(456.0_dp,-456.0_dp,dp))) error stop 'late amplitude commit'
    ! A constant last channel is singular with the floating intercept enabled.
    weight(2,:,3)=1.0_dp
    y(2,:,3)=2.0_dp
    call lombscargle_into(t,y,freq,pkeep,axis=2, &
        weights_nd=weight,floating_mean=.true.,status=st)
    if(st/=irregular_status_singular .or. any(pkeep/=456.0_dp)) &
        error stop 'late singular power commit'

    print '(a)', 'test_signal_irregular_nd_weights passed'
contains
    subroutine check_real(a,b,tol)
        real(dp), intent(in) :: a(:) !! actual spectrum
        real(dp), intent(in) :: b(:) !! reference spectrum
        real(dp), intent(in) :: tol !! maximum coefficient difference
        if(size(a)/=size(b)) error stop 'real shape mismatch'
        if(any(abs(a-b)>tol)) error stop 'real spectral mismatch'
    end subroutine check_real
    subroutine check_complex(a,b,tol)
        complex(dp), intent(in) :: a(:) !! actual complex coefficient array
        complex(dp), intent(in) :: b(:) !! reference complex coefficients
        real(dp), intent(in) :: tol !! maximum complex error
        if(size(a)/=size(b)) error stop 'complex shape mismatch'
        if(any(abs(a-b)>tol)) error stop 'complex spectral mismatch'
    end subroutine check_complex
    subroutine reject(code)
        integer, intent(in) :: code !! error code from failed call
        if(code/=irregular_status_invalid_input .or. any(pkeep/=456.0_dp)) &
            error stop 'invalid call committed partial output'
    end subroutine reject
end program test_signal_irregular_nd_weights
