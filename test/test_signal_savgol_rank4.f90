! SPDX-License-Identifier: MIT
! Rank-4 Savitzky-Golay filtering with shared axis coefficients.
program test_signal_savgol_rank4
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_filter, &
        savgol_status_success, savgol_status_invalid_input
    implicit none

    real(dp) :: x(9,10,7,8), empty(0,8,7,6)
    complex(dp) :: z(9,10,7,8), border
    real(dp), allocatable :: y(:,:,:,:), ref(:), yr(:,:,:,:), yi(:,:,:,:)
    complex(dp), allocatable :: yc(:,:,:,:)
    character(len=8), parameter :: modes(5) = [character(len=8) :: &
        'interp', 'mirror', 'nearest', 'wrap', 'constant']
    integer :: i, j, k, l, ax, m, code

    do l = 1, size(x,4)
        do k = 1, size(x,3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    x(i,j,k,l) = sin(real(i + 2*j + k - l, dp)*0.11_dp) &
                        + 0.004_dp * real(i*j + k*l,dp)
                end do
            end do
        end do
    end do
    z = cmplx(x, -0.5_dp*x + 0.03_dp, dp)
    border = cmplx(-1.25_dp, 2.5_dp, dp)

    do m = 1, size(modes)
        do ax = 1, 4
            y = savgol_filter(x, 5, 3, deriv=1, delta=0.4_dp, &
                mode=trim(modes(m)), cval=-1.25_dp, axis=ax, status=code)
            if (code /= savgol_status_success) error stop 'rank4 real status'
            if (any(shape(y) /= shape(x))) error stop 'rank4 real shape'
            select case (ax)
            case (1)
                ref = savgol_filter(x(:,2,3,4), 5, 3, deriv=1, delta=0.4_dp, &
                    mode=trim(modes(m)), cval=-1.25_dp)
                if (any(abs(y(:,2,3,4)-ref)>5.0e-13_dp)) error stop 'rank4 axis 1'
            case (2)
                ref = savgol_filter(x(2,:,3,4), 5, 3, deriv=1, delta=0.4_dp, &
                    mode=trim(modes(m)), cval=-1.25_dp)
                if (any(abs(y(2,:,3,4)-ref)>5.0e-13_dp)) error stop 'rank4 axis 2'
            case (3)
                ref = savgol_filter(x(2,3,:,4), 5, 3, deriv=1, delta=0.4_dp, &
                    mode=trim(modes(m)), cval=-1.25_dp)
                if (any(abs(y(2,3,:,4)-ref)>5.0e-13_dp)) error stop 'rank4 axis 3'
            case (4)
                ref = savgol_filter(x(2,3,4,:), 5, 3, deriv=1, delta=0.4_dp, &
                    mode=trim(modes(m)), cval=-1.25_dp)
                if (any(abs(y(2,3,4,:)-ref)>5.0e-13_dp)) error stop 'rank4 axis 4'
            end select

            yc = savgol_filter(z, 5, 3, deriv=1, delta=0.4_dp, &
                mode=trim(modes(m)), cval=border, axis=ax, status=code)
            if (code /= savgol_status_success) error stop 'rank4 complex status'
            yr = savgol_filter(real(z,dp), 5, 3, deriv=1, delta=0.4_dp, &
                mode=trim(modes(m)), cval=real(border,dp), axis=ax)
            yi = savgol_filter(aimag(z), 5, 3, deriv=1, delta=0.4_dp, &
                mode=trim(modes(m)), cval=aimag(border), axis=ax)
            if (any(abs(real(yc,dp)-yr)>5.0e-13_dp)) error stop 'rank4 complex real'
            if (any(abs(aimag(yc)-yi)>5.0e-13_dp)) error stop 'rank4 complex imaginary'
        end do
    end do

    y = savgol_filter(x, 1, 0, axis=1, status=code)
    if (code /= savgol_status_success) error stop 'rank4 one-tap status'
    if (any(abs(y-x)>1.0e-13_dp)) error stop 'rank4 one-tap real'
    yc = savgol_filter(z, 1, 0, axis=4, status=code)
    if (code /= savgol_status_success) error stop 'rank4 one-tap complex status'
    if (any(abs(yc-z)>1.0e-13_dp)) error stop 'rank4 one-tap complex'
    y = savgol_filter(x, 5, 2, deriv=4, axis=2, status=code)
    if (code /= savgol_status_success) error stop 'rank4 high derivative status'
    if (any(abs(y)>1.0e-13_dp)) error stop 'rank4 high derivative nonzero'

    y = savgol_filter(x, 4, 2)
    yr = savgol_filter(x, 4, 2, axis=4)
    if (any(abs(y-yr)>5.0e-14_dp)) error stop 'rank4 default axis'
    y = savgol_filter(x, 5, 2, axis=5, status=code)
    if (code /= savgol_status_invalid_input) error stop 'rank4 invalid axis status'
    if (.not. all(ieee_is_nan(y))) error stop 'rank4 invalid axis NaN'
    yc = savgol_filter(z, 5, 2, axis=0, status=code)
    if (code /= savgol_status_invalid_input) error stop 'rank4 complex invalid axis'
    if (.not. all(ieee_is_nan(real(yc,dp)))) error stop 'rank4 complex invalid NaN'
    y = savgol_filter(empty, 5, 2, axis=2, status=code)
    if (code /= savgol_status_success .or. size(y) /= 0) error stop 'rank4 empty other axis'
    y = savgol_filter(empty, 5, 2, axis=1, status=code)
    if (code /= savgol_status_invalid_input) error stop 'rank4 empty interp axis'
    y = savgol_filter(empty, 5, 2, axis=1, mode='nearest', status=code)
    if (code /= savgol_status_success .or. size(y) /= 0) error stop 'rank4 empty nearest axis'
    y = savgol_filter(x, 12, 2, axis=1, status=code)
    if (code /= savgol_status_invalid_input) error stop 'rank4 short interp window'
    y = savgol_filter(x, 12, 2, axis=1, mode='mirror', status=code)
    if (code /= savgol_status_success) error stop 'rank4 short mirror window'
    print '(a)', 'test_signal_savgol_rank4 passed'
end program test_signal_savgol_rank4
