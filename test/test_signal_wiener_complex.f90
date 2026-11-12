! SPDX-License-Identifier: MIT
program test_signal_wiener_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_denoise, only : wiener, denoise_status_success, denoise_status_invalid_input
    implicit none

    complex(dp) :: x(13), whole(26, 6), cube(7, 4, 5)
    complex(dp), allocatable :: y(:), y2(:,:), y3(:,:,:), ref2(:,:), ref3(:,:,:)
    real(dp), allocatable :: real_result(:)
    integer :: i, j, k, status

    x = [(cmplx(real(mod(17*i,29)-14,dp)/7.0_dp, &
        real(mod(19*i+3,23)-11,dp)/9.0_dp,dp), i=1,13)]
    y = wiener(x, kernel_size=3, status=status)
    if (status /= denoise_status_success) error stop '1D implicit complex noise'
    if (size(y) /= 13) error stop 'complex result shape'
    if (any(ieee_is_nan(real(y,dp)))) error stop 'unexpected complex NaNs'
    if (maxval(abs(y - wiener(x, kernel_shape=[3]))) > 1.0e-14_dp) &
        error stop 'scalar/window-shape disagreement'

    ! Directly verify the SciPy pseudo-variance convention for a 3-sample
    ! fully populated neighborhood, with known non-Hermitian second moment.
    y = wiener([cmplx(1.0_dp,0.1_dp,dp), cmplx(2.0_dp,0.2_dp,dp), &
        cmplx(3.0_dp,0.3_dp,dp)], kernel_size=3, &
        noise=cmplx(0.0_dp,0.0_dp,dp), status=status)
    if (status /= denoise_status_success) error stop 'complex noise overload'
    if (abs(y(2) - cmplx(2.0_dp,0.2_dp,dp)) > 1.0e-14_dp) &
        error stop 'zero-noise complex identity interior'

    do j = 1, size(whole,2)
        do i = 1, size(whole,1)
            whole(i,j) = cmplx(real(mod(7*i+11*j,37)-18,dp)/8.0_dp, &
                real(mod(13*i+26*j,31)-15,dp)/10.0_dp,dp)
        end do
    end do
    y2 = wiener(whole(2:26:2,2:6:2),kernel_shape=[4,3], &
        noise=cmplx(0.3_dp,-0.14_dp,dp),status=status)
    if (status /= denoise_status_success) error stop '2D noncontiguous input'
    ref2 = wiener(reshape(whole(2:26:2,2:6:2),[13,3]),kernel_shape=[4,3], &
        noise=cmplx(0.3_dp,-0.14_dp,dp))
    if (any(shape(y2) /= [13,3])) error stop 'strided shape mismatch'
    if (maxval(abs(y2 - ref2)) > 1.0e-14_dp) error stop 'strided input mismatch'

    do k = 1, size(cube,3)
        do j = 1, size(cube,2)
            do i = 1, size(cube,1)
                cube(i,j,k) = cmplx(real(mod(11*i+7*j+5*k,41)-20,dp)/11.0_dp, &
                    real(mod(17*i+13*j+19*k,43)-21,dp)/13.0_dp,dp)
            end do
        end do
    end do
    y3 = wiener(cube(1:7:2,:,:),kernel_shape=[2,3,4],noise=0.15_dp,status=status)
    if (status /= denoise_status_success) error stop '3D real noise overload'
    ref3 = wiener(reshape(cube(1:7:2,:,:),[4,4,5]),kernel_shape=[2,3,4],noise=0.15_dp)
    if (any(shape(y3) /= [4,4,5])) error stop '3D result shape'
    if (maxval(abs(y3 - ref3)) > 1.0e-14_dp) error stop '3D strided input mismatch'

    y = wiener(x,kernel_size=0,status=status)
    if (status /= denoise_status_invalid_input) error stop 'zero window accepted'
    if (.not. all(ieee_is_nan(real(y,dp)))) error stop 'invalid window not NaN'
    y = wiener(x,kernel_shape=[2,2],status=status)
    if (status /= denoise_status_invalid_input) error stop 'wrong kernel rank accepted'
    y = wiener(x,kernel_size=2,kernel_shape=[2],status=status)
    if (status /= denoise_status_invalid_input) error stop 'conflicting kernel options accepted'
    y = wiener(x, noise=ieee_value(0.0_dp,ieee_quiet_nan),status=status)
    if (status /= denoise_status_invalid_input) error stop 'nonfinite real noise accepted'
    y = wiener(x, noise=cmplx(0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),dp),status=status)
    if (status /= denoise_status_invalid_input) error stop 'nonfinite complex noise accepted'

    ! Regression check: existing real and integer Wiener overloads remain available.
    real_result = wiener([1.0_dp,2.0_dp,3.0_dp,2.0_dp,1.0_dp],kernel_size=3,status=status)
    if (status /= denoise_status_success) error stop 'legacy real Wiener overload'
    if (any(ieee_is_nan(real_result))) error stop 'legacy real Wiener NaN'
    real_result = wiener([1,2,3,2,1],kernel_size=3,status=status)
    if (status /= denoise_status_success) error stop 'legacy integer Wiener overload'

    print '(a)', 'PASS: complex Wiener native tests (1D/2D/3D, strides, noise, validation, legacy real)'
end program test_signal_wiener_complex
