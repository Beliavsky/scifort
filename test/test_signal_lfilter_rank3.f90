! SPDX-License-Identifier: MIT
! Tests for rank-three signal filtering with channel-specific initial states.
program test_signal_lfilter_rank3
    use scifort_signal, only : dp, lfilter, lfilter_status_success, lfilter_status_invalid_input
    implicit none
    real(dp) :: x(3,4,5), b(3), a(2), zi(2,4,5)
    real(dp), allocatable :: y(:,:,:), zf(:,:,:), line(:), state(:)
    complex(dp) :: cx(3,4,5), cb(3), ca(2)
    complex(dp), allocatable :: cy(:,:,:), czf(:,:,:), cline(:), cstate(:)
    integer :: i,j,k,st

    do k = 1,5
        do j = 1,4
            do i = 1,3
                x(i,j,k) = sin(real(i+3*j+2*k,dp))
            end do
        end do
    end do
    b = [0.1_dp,0.3_dp,0.2_dp]
    a = [1.0_dp,-0.4_dp]
    zi = 0.2_dp
    y = lfilter(b,a,x,axis=1,zi=zi,zf=zf,status=st)
    if (st /= lfilter_status_success) error stop 'rank3 axis1 status'
    if (any(shape(zf) /= [2,4,5])) error stop 'rank3 axis1 state shape'
    do k = 1,5
        do j = 1,4
            line = lfilter(b,a,x(:,j,k),zi=zi(:,j,k),zf=state)
            if (maxval(abs(line-y(:,j,k))) > 1.0e-13_dp) error stop 'rank3 axis1 result'
            if (maxval(abs(state-zf(:,j,k))) > 1.0e-13_dp) error stop 'rank3 axis1 state'
        end do
    end do
    y = lfilter(b,a,x,axis=2,zf=zf,status=st)
    if (st /= lfilter_status_success) error stop 'rank3 axis2 status'
    do k = 1,5
        do i = 1,3
            line = lfilter(b,a,x(i,:,k),zf=state)
            if (maxval(abs(line-y(i,:,k))) > 1.0e-13_dp) error stop 'rank3 axis2 result'
            if (maxval(abs(state-zf(i,:,k))) > 1.0e-13_dp) error stop 'rank3 axis2 state'
        end do
    end do
    cx = cmplx(x,-0.2_dp*x,dp)
    cb = cmplx(b,0.1_dp*b,dp)
    ca = cmplx(a,0.05_dp*a,dp)
    cy = lfilter(cb,ca,cx,zf=czf,status=st)
    if (st /= lfilter_status_success) error stop 'complex rank3 status'
    do j = 1,4
        do i = 1,3
            cline = lfilter(cb,ca,cx(i,j,:),zf=cstate)
            if (maxval(abs(cline-cy(i,j,:))) > 1.0e-13_dp) error stop 'complex rank3 result'
            if (maxval(abs(cstate-czf(i,j,:))) > 1.0e-13_dp) error stop 'complex rank3 state'
        end do
    end do
    y = lfilter(b,a,x,axis=4,status=st)
    if (st /= lfilter_status_invalid_input) error stop 'invalid rank3 axis'
    print '(a)', 'test_signal_lfilter_rank3 passed'
end program test_signal_lfilter_rank3
