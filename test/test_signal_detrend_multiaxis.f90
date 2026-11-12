! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_detrend_multiaxis
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_detrend, only : detrend, detrend_status_success, detrend_status_invalid_input
    implicit none
    real(dp) :: r2(17,9), r3(7,11,13)
    real(dp), allocatable :: d2(:,:), d3(:,:,:), original(:)
    complex(dp) :: c2(17,9), c3(7,11,13)
    complex(dp), allocatable :: z2(:,:), z3(:,:,:), z(:)
    integer :: i,j,k,st,axis

    do j=1,9
        do i=1,17
            r2(i,j) = 3.0_dp+2.0_dp*i+real(j,dp)*0.2_dp
        end do
    end do
    c2 = cmplx(r2, -2.0_dp*r2+3.0_dp,dp)
    d2 = detrend(r2,axis=1,status=st)
    if (st /= detrend_status_success) error stop 'real axis1 status'
    if (maxval(abs(d2)) > 1.0e-12_dp) error stop 'real axis1 linear'
    z2 = detrend(c2,axis=1,status=st)
    if (st /= detrend_status_success) error stop 'complex axis1 status'
    if (maxval(abs(z2)) > 1.0e-12_dp) error stop 'complex axis1 linear'
    d2 = detrend(r2,kind='constant',status=st)
    if (st /= detrend_status_success) error stop 'constant status'
    if (maxval(abs(sum(d2,dim=2))) > 1.0e-12_dp) error stop 'constant mean'
    d2 = detrend(r2,axis=2,status=st)
    if (st /= detrend_status_success) error stop 'real axis2 status'
    if (maxval(abs(d2)) > 1.0e-12_dp) error stop 'real axis2 linear'
    z2 = detrend(c2,axis=2)
    if (maxval(abs(z2)) > 1.0e-12_dp) error stop 'complex axis2 linear'
    d2 = detrend(r2,axis=0,status=st)
    if (st /= detrend_status_invalid_input) error stop 'bad axis'
    if (.not. all(ieee_is_nan(d2))) error stop 'bad axis output'
    d2 = detrend(r2,axis=1,breaks=[18,6,6,1,11],status=st)
    if (st /= detrend_status_success) error stop 'broken segments'
    if (maxval(abs(d2)) > 1.0e-12_dp) error stop 'broken segments residual'
    do k=1,13
        do j=1,11
            do i=1,7
                r3(i,j,k) = real(i,dp)+2.0_dp*real(j,dp)-real(k,dp)
            end do
        end do
    end do
    c3 = cmplx(r3,3.0_dp*r3+1.0_dp,dp)
    do axis=1,3
        d3 = detrend(r3,axis=axis,status=st)
        if (st /= detrend_status_success) error stop 'real rank3 status'
        if (maxval(abs(d3)) > 2.0e-12_dp) error stop 'real rank3'
        z3 = detrend(c3,axis=axis,status=st)
        if (st /= detrend_status_success) error stop 'complex rank3 status'
        if (maxval(abs(z3)) > 2.0e-12_dp) error stop 'complex rank3'
    end do
    z3 = detrend(c3(::2,::3,::4),axis=1)
    if (maxval(abs(z3)) > 2.0e-12_dp) error stop 'strided complex'
    z = detrend(c2(:,3),kind='constant')
    if (abs(sum(z)) > 1.0e-12_dp) error stop 'complex constant'
    d3 = detrend(r3,axis=3,kind='invalid',status=st)
    if (st /= detrend_status_invalid_input) error stop 'invalid trend'
    if (.not. all(ieee_is_nan(d3))) error stop 'invalid trend output'
    z2 = detrend(c2,axis=3,status=st)
    if (st /= detrend_status_invalid_input) error stop 'invalid complex axis'
    original = detrend([1.0_dp,3.0_dp,4.0_dp,8.0_dp])
    if (maxval(abs(original-[0.3_dp,0.1_dp,-1.1_dp,0.7_dp])) > 1.0e-12_dp) &
        error stop 'backward compatibility'
    print '(a)', 'test_signal_detrend_multiaxis passed'
end program test_signal_detrend_multiaxis
