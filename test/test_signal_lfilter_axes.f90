! SPDX-License-Identifier: MIT
! Tests for rank-two filter axes and initial/final states.
program test_signal_lfilter_axes
    use scifort_signal, only : dp, lfilter, lfilter_status_success, lfilter_status_invalid_input
    implicit none
    real(dp) :: b(2), a(2), x(4,3), states(4,1), states1(1,3)
    real(dp), allocatable :: y(:,:), zf(:,:), row(:), endstate(:)
    complex(dp), allocatable :: cy(:,:), czf(:,:), crow(:), cend(:)
    complex(dp) :: cb(2), ca(2), cx(4,3)
    integer :: i,j,code

    b = [0.6_dp,0.2_dp]
    a = [1.0_dp,-0.3_dp]
    do j = 1,3
        do i = 1,4
            x(i,j) = real(2*i-j,dp)
        end do
    end do
    states = 0.1_dp
    y = lfilter(b,a,x,zi=states,zf=zf,status=code)
    if (code /= lfilter_status_success) error stop 'rank2 lfilter status'
    if (any(shape(zf) /= [4,1])) error stop 'rank2 final states shape'
    do i = 1,4
        row = lfilter(b,a,x(i,:),zi=states(i,:),zf=endstate)
        if (maxval(abs(row-y(i,:))) > 1.0e-13_dp) error stop 'row filter mismatch'
        if (maxval(abs(endstate-zf(i,:))) > 1.0e-13_dp) error stop 'row state mismatch'
    end do
    states1 = 0.25_dp
    y = lfilter(b,a,x,zi=states1,zf=zf,status=code,axis=1)
    if (code /= lfilter_status_success) error stop 'column filter status'
    if (any(shape(zf) /= [1,3])) error stop 'column final states shape'
    do j = 1,3
        row = lfilter(b,a,x(:,j),zi=states1(:,j),zf=endstate)
        if (maxval(abs(row-y(:,j))) > 1.0e-13_dp) error stop 'column filter mismatch'
        if (maxval(abs(endstate-zf(:,j))) > 1.0e-13_dp) error stop 'column state mismatch'
    end do
    cb = cmplx(b,-0.2_dp*b,dp)
    ca = cmplx(a,0.1_dp*a,dp)
    cx = cmplx(x,x*0.1_dp,dp)
    cy = lfilter(cb,ca,cx,zf=czf,status=code,axis=2)
    if (code /= lfilter_status_success) error stop 'complex matrix lfilter'
    do i = 1,4
        crow = lfilter(cb,ca,cx(i,:),zf=cend)
        if (maxval(abs(crow-cy(i,:))) > 1.0e-13_dp) error stop 'complex row mismatch'
        if (maxval(abs(cend-czf(i,:))) > 1.0e-13_dp) error stop 'complex row state mismatch'
    end do
    y = lfilter(b,a,x,axis=0,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid filter axis'
    y = lfilter(b,a,x,zi=states1,status=code)
    if (code /= lfilter_status_invalid_input) error stop 'invalid initial state shape'
    print '(a)', 'test_signal_lfilter_axes passed'
end program test_signal_lfilter_axes
