! SPDX-License-Identifier: MIT
! Assumed-rank Savitzky-Golay filtering at ranks 5 and 15.
program test_signal_savgol_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, savgol_filter, savgol_filter_nd, &
        savgol_status_success, savgol_status_invalid_input
    implicit none
    real(dp) :: x5(6,4,3,2,7), y5(6,4,3,2,7), bad5(6,4,3,2,6)
    real(dp) :: x15(7,2,3,2,5,1,1,1,1,1,1,1,1,1,7)
    real(dp) :: y15(7,2,3,2,5,1,1,1,1,1,1,1,1,1,7)
    real(dp) :: bad2(2,5), x2(3,5), out2(3,5)
    complex(dp) :: c5(6,4,3,2,7), d5(6,4,3,2,7)
    complex(dp) :: c15(7,2,3,2,5,1,1,1,1,1,1,1,1,1,7)
    complex(dp) :: d15(7,2,3,2,5,1,1,1,1,1,1,1,1,1,7)
    real(dp), allocatable :: reference(:)
    character(len=8), parameter :: modes(5) = [character(len=8) :: &
        'interp', 'mirror', 'nearest', 'wrap', 'constant']
    integer :: i, ax, m, code

    x5 = reshape([(sin(0.035_dp*real(i,dp)) + 0.003_dp*i, i=1,size(x5))], shape(x5))
    c5 = cmplx(x5, 0.4_dp - 0.6_dp*x5, dp)
    do m = 1, size(modes)
        do ax = 1, 5
            call savgol_filter_nd(x5, y5, 2, 1, axis=ax, mode=trim(modes(m)), &
                cval=-1.2_dp, status=code)
            if (code /= savgol_status_success) error stop 'rank5 status'
            select case (ax)
            case (1)
                reference = savgol_filter(x5(:,2,2,2,2), 2, 1, mode=trim(modes(m)), cval=-1.2_dp)
                if (any(abs(y5(:,2,2,2,2)-reference)>1.0e-12_dp)) error stop 'rank5 axis1'
            case (2)
                reference = savgol_filter(x5(2,:,2,2,2), 2, 1, mode=trim(modes(m)), cval=-1.2_dp)
                if (any(abs(y5(2,:,2,2,2)-reference)>1.0e-12_dp)) error stop 'rank5 axis2'
            case (3)
                reference = savgol_filter(x5(2,2,:,2,2), 2, 1, mode=trim(modes(m)), cval=-1.2_dp)
                if (any(abs(y5(2,2,:,2,2)-reference)>1.0e-12_dp)) error stop 'rank5 axis3'
            case (4)
                reference = savgol_filter(x5(2,2,2,:,2), 2, 1, mode=trim(modes(m)), cval=-1.2_dp)
                if (any(abs(y5(2,2,2,:,2)-reference)>1.0e-12_dp)) error stop 'rank5 axis4'
            case (5)
                reference = savgol_filter(x5(2,2,2,2,:), 2, 1, mode=trim(modes(m)), cval=-1.2_dp)
                if (any(abs(y5(2,2,2,2,:)-reference)>1.0e-12_dp)) error stop 'rank5 axis5'
            end select
            call savgol_filter_nd(c5, d5, 2, 1, axis=ax, mode=trim(modes(m)), &
                cval=cmplx(-1.2_dp, 0.2_dp,dp), status=code)
            if (code /= savgol_status_success) error stop 'rank5 complex status'
        end do
    end do

    ! Rank 15 is the maximum standard rank. Check non-unit strides and
    ! filtering at the beginning, middle, and end of the shape.
    x15 = reshape([(sin(0.047_dp*real(i,dp)) + 0.003_dp*i, i=1,size(x15))], shape(x15))
    c15 = cmplx(x15, 0.5_dp - 0.3_dp*x15, dp)
    do m = 1, size(modes)
        call savgol_filter_nd(x15, y15, 5, 2, axis=15, mode=trim(modes(m)), status=code)
        if (code /= savgol_status_success) error stop 'rank15 last axis status'
        reference = savgol_filter(x15(2,1,2,1,3,1,1,1,1,1,1,1,1,1,:), &
            5, 2, mode=trim(modes(m)))
        if (any(abs(y15(2,1,2,1,3,1,1,1,1,1,1,1,1,1,:)-reference)>1.0e-12_dp)) &
            error stop 'rank15 last axis'
        call savgol_filter_nd(x15, y15, 3, 2, axis=3, mode=trim(modes(m)), status=code)
        if (code /= savgol_status_success) error stop 'rank15 middle status'
        reference = savgol_filter(x15(2,1,:,1,3,1,1,1,1,1,1,1,1,1,2), &
            3, 2, mode=trim(modes(m)))
        if (any(abs(y15(2,1,:,1,3,1,1,1,1,1,1,1,1,1,2)-reference)>1.0e-12_dp)) &
            error stop 'rank15 middle axis'
        call savgol_filter_nd(x15, y15, 5, 3, axis=1, mode=trim(modes(m)), status=code)
        if (code /= savgol_status_success) error stop 'rank15 first status'
        reference = savgol_filter(x15(:,1,2,1,3,1,1,1,1,1,1,1,1,1,2), &
            5, 3, mode=trim(modes(m)))
        if (any(abs(y15(:,1,2,1,3,1,1,1,1,1,1,1,1,1,2)-reference)>1.0e-12_dp)) &
            error stop 'rank15 first axis'
        call savgol_filter_nd(c15, d15, 5, 3, axis=1, mode=trim(modes(m)), status=code)
        if (code /= savgol_status_success) error stop 'rank15 complex status'
    end do

    call savgol_filter_nd(x15, y15, 5, 2, status=code)
    if (code /= savgol_status_success) error stop 'rank15 default axis'
    reference = savgol_filter(x15(2,1,2,1,3,1,1,1,1,1,1,1,1,1,:), 5, 2)
    if (any(abs(y15(2,1,2,1,3,1,1,1,1,1,1,1,1,1,:)-reference)>1.0e-12_dp)) &
        error stop 'rank15 default last axis'

    ! Noncontiguous actual arguments must copy in and out correctly.
    y5 = -999.0_dp
    call savgol_filter_nd(x5(1:6:2,:,:,:,:), y5(1:6:2,:,:,:,:), &
        3, 1, axis=2, mode='wrap', status=code)
    if (code /= savgol_status_success) error stop 'noncontiguous status'
    reference = savgol_filter(x5(1,:,2,2,2), 3, 1, mode='wrap')
    if (any(abs(y5(1,:,2,2,2)-reference)>1.0e-12_dp)) error stop 'noncontiguous result'
    if (any(abs(y5(2,:,:,:,:)+999.0_dp)>0.0_dp)) error stop 'noncontiguous overwrite'

    call savgol_filter_nd(x5, bad5, 3, 1, status=code)
    if (code /= savgol_status_invalid_input) error stop 'shape mismatch status'
    if (.not. all(ieee_is_nan(bad5))) error stop 'shape mismatch output'
    bad2 = 0.0_dp
    call savgol_filter_nd(x5, bad2, 3, 1, status=code)
    if (code /= savgol_status_invalid_input) error stop 'rank mismatch status'
    if (.not. all(ieee_is_nan(bad2))) error stop 'rank mismatch output'
    call savgol_filter_nd(x15, y15, 5, 2, axis=16, status=code)
    if (code /= savgol_status_invalid_input) error stop 'bad rank15 axis'
    if (.not. all(ieee_is_nan(y15))) error stop 'bad rank15 output'

    x2 = 0.0_dp
    call savgol_filter_nd(x2(:,1:0), out2(:,1:0), 1, 0, status=code)
    if (code /= savgol_status_invalid_input) error stop 'empty default axis status'
    call savgol_filter_nd(x2(:,1:0), out2(:,1:0), 1, 0, axis=1, status=code)
    if (code /= savgol_status_success) error stop 'empty other axis status'

    print '(a)', 'test_signal_savgol_nd passed'
end program test_signal_savgol_nd
