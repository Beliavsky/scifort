! SPDX-License-Identifier: MIT
! Native regressions for rank-independent real analytic signals.
program test_signal_hilbert_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_hilbert_nd, only : hilbert_into, hilbert2_into
    use scifort_signal_hilbert, only : hilbert, hilbert2, &
        hilbert_status_success, hilbert_status_invalid_input
    implicit none
    real(dp) :: x(5,6,4), nan, x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,6)
    real(dp) :: x4(4,5,6), empty(0,6)
    complex(dp) :: y(5,8,4), z(5,6,7), z4(4,5,6)
    complex(dp) :: scratch(5,8,8), before(5,8,8), yempty(0,6)
    complex(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,9)
    complex(dp), allocatable :: ref(:), ref2(:,:)
    integer :: i,j,k, st
    real(dp), parameter :: tol = 6.0e-12_dp

    do k=1,4
        do j=1,6
            do i=1,5
                x(i,j,k) = sin(real(3*i+2*j+7*k,dp)) + 0.1_dp*i*j
            end do
        end do
    end do
    call hilbert_into(x,y,axis=2,n=8,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert axis status'
    do k=1,4
        do i=1,5
            ref = hilbert(x(i,:,k),n=8)
            if (maxval(abs(y(i,:,k)-ref)) > tol) error stop 'hilbert axis values'
        end do
    end do

    call hilbert2_into(x,z,axes=[1,3],n=[5,7],status=st)
    if (st /= hilbert_status_success) error stop 'hilbert2 axis status'
    do j=1,6
        ref2 = hilbert2(x(:,j,:),n=[5,7])
        if (maxval(abs(z(:,j,:)-ref2)) > tol) error stop 'hilbert2 axis values'
    end do

    do k=1,6
        do j=1,5
            do i=1,4
                x4(i,j,k) = cos(real(i+2*j+3*k,dp))
            end do
        end do
    end do
    call hilbert2_into(x4,z4,axes=[-3,-1],status=st)
    if (st /= hilbert_status_success) error stop 'hilbert2 even status'
    do j=1,5
        ref2 = hilbert2(x4(:,j,:))
        if (maxval(abs(z4(:,j,:)-ref2)) > tol) error stop 'hilbert2 even parity'
    end do

    do k=1,6
        x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,k) = real(k*k-3*k,dp)
    end do
    call hilbert_into(x15,y15,axis=-1,n=9,status=st)
    if (st /= hilbert_status_success) error stop 'rank 15 status'
    ref = hilbert(x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:),n=9)
    if (maxval(abs(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-ref)) > tol) &
        error stop 'rank 15 values'

    scratch = cmplx(-81.0_dp,29.0_dp,dp)
    call hilbert_into(x(:,:,1:4:2),scratch(:,:,2:8:4),axis=2,n=8,status=st)
    if (st /= hilbert_status_success) error stop 'strided destination status'
    do k=1,2
        do i=1,5
            ref = hilbert(x(i,:,2*k-1),n=8)
            if (maxval(abs(scratch(i,:,4*k-2)-ref)) > tol) &
                error stop 'strided destination values'
        end do
    end do
    if (any(scratch(:,:,1) /= cmplx(-81.0_dp,29.0_dp,dp))) &
        error stop 'strided destination guard'

    before = scratch
    call hilbert_into(x(:,:,1:4:2),scratch(:,:,2:8:4),axis=2,n=0,status=st)
    if (st /= hilbert_status_invalid_input .or. any(scratch /= before)) &
        error stop 'invalid fft length modified destination'
    call hilbert_into(x(:,:,1:4:2),scratch(:,:,2:8:4),axis=4,n=8,status=st)
    if (st /= hilbert_status_invalid_input .or. any(scratch /= before)) &
        error stop 'invalid axis modified destination'
    call hilbert2_into(x(:,:,1:4:2),scratch(:,:,2:8:4),axes=[2,2],status=st)
    if (st /= hilbert_status_invalid_input .or. any(scratch /= before)) &
        error stop 'duplicate axes modified destination'
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    x(2,1,3)=nan
    call hilbert_into(x(:,:,1:4:2),scratch(:,:,2:8:4),axis=2,n=8,status=st)
    if (st /= hilbert_status_invalid_input .or. any(scratch /= before)) &
        error stop 'nonfinite later lane modified destination'

    call hilbert_into(empty,yempty,axis=2,status=st)
    if (st /= hilbert_status_success) error stop 'empty batch status'

    print '(a)', 'test_signal_hilbert_nd passed'
end program test_signal_hilbert_nd
