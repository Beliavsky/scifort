! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_denoise_public
    use scifort_signal, only : dp, medfilt, medfilt2d, wiener, &
        denoise_status_success, denoise_status_invalid_input
    implicit none
    integer :: status
    integer :: x(5) = [1, 2, 90, 3, 4]
    integer, allocatable :: med(:)
    real(dp), allocatable :: estimate(:), image(:,:)
    complex(dp) :: complex_signal(5)
    complex(dp), allocatable :: complex_result(:)

    med = medfilt(x, status=status)
    if (status /= denoise_status_success) error stop 'public median status'
    if (any(med /= [1, 2, 3, 4, 3])) error stop 'public median result'
    estimate = wiener(x, noise=0.0_dp, status=status)
    if (status /= denoise_status_success) error stop 'public wiener status'
    if (maxval(abs(estimate - real(x,dp))) > 1.0e-12_dp) error stop 'public wiener result'
    image = medfilt2d(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp], [2,2]), &
        kernel_shape=[1,3], status=status)
    if (status /= denoise_status_success .or. any(shape(image) /= [2,2])) error stop 'public medfilt2d'
    image = medfilt2d(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp], [2,2]), &
        kernel_shape=[2,3], status=status)
    if (status /= denoise_status_invalid_input) error stop 'public invalid medfilt2d'
    complex_signal = cmplx(real(x, dp), real([2,1,0,-1,-2],dp), dp)
    complex_result = wiener(complex_signal, noise=0.1_dp, status=status)
    if (status /= denoise_status_success) error stop 'public real-noise complex Wiener'
    if (size(complex_result) /= size(complex_signal)) error stop 'public complex Wiener shape'
    complex_result = wiener(complex_signal, noise=cmplx(0.1_dp,-0.2_dp,dp), status=status)
    if (status /= denoise_status_success) error stop 'public complex-noise Wiener'
    print *, 'PASS test_signal_denoise_public' 
end program test_signal_denoise_public
