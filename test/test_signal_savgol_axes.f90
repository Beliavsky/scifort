! SPDX-License-Identifier: MIT
! Rank-2/rank-3 Savitzky-Golay filters and complex-valued signals.
program test_signal_savgol_axes
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_filter, &
        savgol_status_success, savgol_status_invalid_input
    implicit none
    real(dp) :: x2(13,7), x3(9,11,6), xempty(0,9)
    complex(dp) :: z2(13,7), z3(9,11,6), shift
    real(dp), allocatable :: y2(:,:), y3(:,:,:), ref(:)
    complex(dp), allocatable :: cz2(:,:), cz3(:,:,:), cz1(:)
    integer :: i, j, k, ax, code

    do j = 1, size(x2,2)
        do i = 1, size(x2,1)
            x2(i,j) = cos(real(2*i-j,dp)*0.19_dp) + 0.03_dp*i*j
        end do
    end do
    do k = 1, size(x3,3)
        do j = 1, size(x3,2)
            do i = 1, size(x3,1)
                x3(i,j,k) = sin(real(i+2*j-k,dp)*0.22_dp) + 0.01_dp*j*k
            end do
        end do
    end do
    z2 = cmplx(x2,0.5_dp*x2,dp)
    z3 = cmplx(x3,-0.3_dp*x3,dp)
    shift = cmplx(4.0_dp,-2.0_dp,dp)

    do ax = 1, 2
        y2 = savgol_filter(x2,5,2,axis=ax,status=code)
        if (code /= savgol_status_success) error stop 'rank2 interpolation status'
        if (any(shape(y2) /= shape(x2))) error stop 'rank2 output shape'
        if (ax == 1) then
            do j = 1, size(x2,2)
                ref = savgol_filter(x2(:,j),5,2)
                if (any(abs(y2(:,j)-ref)>3.0e-14_dp)) error stop 'rank2 first axis'
            end do
        else
            do i = 1, size(x2,1)
                ref = savgol_filter(x2(i,:),5,2)
                if (any(abs(y2(i,:)-ref)>3.0e-14_dp)) error stop 'rank2 second axis'
            end do
        end if
    end do
    y2 = savgol_filter(x2,5,2)
    if (any(abs(y2-savgol_filter(x2,5,2,axis=2))>3.0e-14_dp)) error stop 'rank2 default axis'

    do ax = 1, 3
        y3 = savgol_filter(x3,5,3,axis=ax,deriv=1,delta=0.3_dp,mode='mirror',status=code)
        if (code /= savgol_status_success) error stop 'rank3 derivative status'
        select case (ax)
        case (1)
            do k = 1, size(x3,3)
                do j = 1, size(x3,2)
                    ref = savgol_filter(x3(:,j,k),5,3,deriv=1,delta=0.3_dp,mode='mirror')
                    if (any(abs(y3(:,j,k)-ref)>3.0e-13_dp)) error stop 'rank3 first axis'
                end do
            end do
        case (2)
            do k = 1, size(x3,3)
                do i = 1, size(x3,1)
                    ref = savgol_filter(x3(i,:,k),5,3,deriv=1,delta=0.3_dp,mode='mirror')
                    if (any(abs(y3(i,:,k)-ref)>3.0e-13_dp)) error stop 'rank3 second axis'
                end do
            end do
        case (3)
            do j = 1, size(x3,2)
                do i = 1, size(x3,1)
                    ref = savgol_filter(x3(i,j,:),5,3,deriv=1,delta=0.3_dp,mode='mirror')
                    if (any(abs(y3(i,j,:)-ref)>3.0e-13_dp)) error stop 'rank3 third axis'
                end do
            end do
        end select
    end do

    ! Complex samples are filtered componentwise, including constant padding.
    do ax = 1, 2
        cz2 = savgol_filter(z2,4,2,axis=ax,mode='constant',cval=shift,status=code)
        if (code /= savgol_status_success) error stop 'rank2 complex status'
        y2 = savgol_filter(real(z2,dp),4,2,axis=ax,mode='constant',cval=real(shift,dp))
        if (any(abs(real(cz2,dp)-y2)>3.0e-14_dp)) error stop 'rank2 complex real'
        y2 = savgol_filter(aimag(z2),4,2,axis=ax,mode='constant',cval=aimag(shift))
        if (any(abs(aimag(cz2)-y2)>3.0e-14_dp)) error stop 'rank2 complex imaginary'
    end do
    do ax = 1, 3
        cz3 = savgol_filter(z3,5,2,axis=ax,mode='wrap',status=code)
        if (code /= savgol_status_success) error stop 'rank3 complex status'
        y3 = savgol_filter(real(z3,dp),5,2,axis=ax,mode='wrap')
        if (any(abs(real(cz3,dp)-y3)>3.0e-14_dp)) error stop 'rank3 complex real'
        y3 = savgol_filter(aimag(z3),5,2,axis=ax,mode='wrap')
        if (any(abs(aimag(cz3)-y3)>3.0e-14_dp)) error stop 'rank3 complex imaginary'
    end do
    cz1 = savgol_filter(z2(:,1),5,2,mode='nearest',status=code)
    if (code /= savgol_status_success) error stop 'rank1 complex status'
    if (any(abs(real(cz1,dp)-savgol_filter(real(z2(:,1),dp),5,2,mode='nearest'))>3.0e-14_dp)) &
        error stop 'rank1 complex real'

    y3 = savgol_filter(x3,5,2,axis=0,status=code)
    if (code /= savgol_status_invalid_input) error stop 'bad rank3 axis status'
    if (.not. all(ieee_is_nan(y3))) error stop 'bad rank3 axis NaNs'
    y2 = savgol_filter(x2,9,2,axis=2,status=code)
    if (code /= savgol_status_invalid_input) error stop 'short rank2 interpolation'
    y2 = savgol_filter(x2,9,2,axis=2,mode='nearest',status=code)
    if (code /= savgol_status_success) error stop 'short rank2 boundary mode'
    y2 = savgol_filter(xempty,5,2,axis=2,status=code)
    if (code /= savgol_status_success .or. size(y2)/=0) error stop 'empty nonfiltered axis'
    y2 = savgol_filter(xempty,5,2,axis=1,status=code)
    if (code /= savgol_status_invalid_input) error stop 'empty filtered axis'

    print '(a)', 'test_signal_savgol_axes passed'
end program test_signal_savgol_axes
