! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_denoise
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_denoise, only : medfilt, medfilt2d, wiener, &
        denoise_status_success, denoise_status_invalid_input
    implicit none
    real(dp) :: image(5,5), expected(5,5), cube(3,4,2)
    real(dp) :: strip(5), frame(4,6), expected_frame(4,6)
    real(dp), allocatable :: a(:), z2(:,:), z3(:,:,:)
    integer, allocatable :: zi2(:,:), zi3(:,:,:)
    integer :: ints(5,5), j, i, k, status

    do j = 1, 5
        do i = 1, 5
            image(i,j) = real(5*(i-1) + j-1, dp)
        end do
    end do
    ints = int(image)
    expected(1,:) = [0,1,2,3,0]
    expected(2,:) = [1,6,7,8,4]
    expected(3,:) = [6,11,12,13,9]
    expected(4,:) = [11,16,17,18,14]
    expected(5,:) = [0,16,17,18,0]
    z2 = medfilt2d(image, status=status)
    if (status /= denoise_status_success) error stop 'medfilt2d status'
    if (any(z2 /= expected)) error stop 'medfilt2d default example'
    zi2 = medfilt2d(ints, kernel_shape=[5,3], status=status)
    if (status /= denoise_status_success) error stop 'medfilt2d integer status'
    expected(1,:) = [0,1,2,3,0]
    expected(2,:) = [0,6,7,8,3]
    expected(3,:) = [5,11,12,13,8]
    expected(4,:) = [5,11,12,13,8]
    expected(5,:) = [0,11,12,13,0]
    if (any(zi2 /= int(expected))) error stop 'medfilt2d integer rectangular'

    strip = [1.0_dp,2.0_dp,50.0_dp,3.0_dp,4.0_dp]
    a = medfilt(strip, status=status)
    if (status /= denoise_status_success) error stop 'median rank1 status'
    if (any(a /= [1.0_dp,2.0_dp,3.0_dp,4.0_dp,3.0_dp])) error stop 'median rank1'
    a = medfilt(strip, kernel_size=0, status=status)
    if (status /= denoise_status_invalid_input .or. any(a /= 0.0_dp)) error stop 'median zero window'
    z2 = medfilt(image, kernel_shape=[3,2], status=status)
    if (status /= denoise_status_invalid_input .or. any(z2 /= 0.0_dp)) error stop 'median even window'
    z2 = medfilt(image, kernel_shape=[3], status=status)
    if (status /= denoise_status_invalid_input) error stop 'median shape mismatch'
    z2 = medfilt(image, kernel_size=3, kernel_shape=[3,3], status=status)
    if (status /= denoise_status_invalid_input) error stop 'median conflicting arguments'
    a = medfilt([1,2,3], kernel_size=7, status=status)
    if (status /= denoise_status_success) error stop 'median kernel exceeds array'
    if (any(a /= 0.0_dp)) error stop 'median zero padding oversized kernel'

    do k = 1, 2
        do j = 1, 4
            do i = 1, 3
                cube(i,j,k) = sin(real(i+7*j+11*k, dp)) + 0.1_dp * real(i*j*k,dp)
            end do
        end do
    end do
    z3 = medfilt(cube, kernel_shape=[3,1,3], status=status)
    if (status /= denoise_status_success) error stop 'median rank3 status'
    if (any(shape(z3) /= shape(cube))) error stop 'median rank3 shape'
    zi3 = medfilt(int(cube * 1000.0_dp), kernel_shape=[1,3,3], status=status)
    if (status /= denoise_status_success) error stop 'median integer rank3 status'
    if (any(shape(zi3) /= shape(cube))) error stop 'median integer rank3 shape'

    ! Noise power zero should preserve nonconstant observations (local variance > 0).
    a = wiener(strip, kernel_size=3, noise=0.0_dp, status=status)
    if (status /= denoise_status_success) error stop 'wiener rank1 status'
    if (maxval(abs(a - strip)) > 2.0e-13_dp) error stop 'wiener noiseless'
    a = wiener(strip, kernel_size=1, noise=0.5_dp, status=status)
    if (status /= denoise_status_success) error stop 'wiener singleton status'
    if (maxval(abs(a - strip)) > 2.0e-13_dp) error stop 'wiener scalar noise'
    a = wiener(strip, kernel_size=1, noise=0.0_dp, status=status)
    if (status /= denoise_status_success) error stop 'wiener zero local variance status'
    if (.not. all(ieee_is_nan(a))) error stop 'wiener zero variance follows SciPy'

    z2 = wiener(image, kernel_shape=[2,5], noise=1.0_dp, status=status)
    if (status /= denoise_status_success) error stop 'wiener even rectangular'
    if (any(shape(z2) /= shape(image))) error stop 'wiener rank2 shape'
    z2 = wiener(image, kernel_shape=[-1,2], status=status)
    if (status /= denoise_status_invalid_input) error stop 'wiener invalid window'
    z3 = wiener(cube, kernel_shape=[2,3,4], status=status)
    if (status /= denoise_status_success) error stop 'wiener rank3 status'
    if (any(shape(z3) /= shape(cube))) error stop 'wiener rank3 shape'
    z3 = wiener(int(cube), kernel_shape=[1,2,3], noise=1.2_dp, status=status)
    if (status /= denoise_status_success) error stop 'wiener integer rank3 status'
    if (any(shape(z3) /= shape(cube))) error stop 'wiener integer rank3 shape'

    frame = 1.0_dp
    expected_frame = 1.0_dp
    z2 = wiener(frame, kernel_size=1, noise=0.1_dp, status=status)
    if (any(z2 /= expected_frame)) error stop 'wiener flat image supplied noise'
    z2 = medfilt(frame, kernel_size=1, status=status)
    if (any(z2 /= frame)) error stop 'median identity window'
    z2 = medfilt(frame(:0,:), status=status)
    if (size(z2) /= 0 .or. status /= denoise_status_success) error stop 'median empty input'
    z2 = wiener(frame(:0,:), status=status)
    if (size(z2) /= 0 .or. status /= denoise_status_success) error stop 'wiener empty input'

    print *, 'PASS test_signal_denoise'
end program test_signal_denoise
