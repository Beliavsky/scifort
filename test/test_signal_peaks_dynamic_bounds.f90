! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style position-dependent bounds are sampled at peak centers.
program test_signal_peaks_dynamic_bounds
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_signal, only : dp, find_peaks, peak_result, &
        peak_status_success, peak_status_invalid_input
    implicit none

    real(dp) :: x(14), bound(14), nanvalue
    real(dp) :: strided_x(28), strided_bound(28)
    integer :: tested
    type(peak_result) :: r

    tested = 0
    x = [0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 4.0_dp, &
        0.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp]
    r = find_peaks(x)
    call check(r, [2, 4, 6, 9, 12])

    bound = -100.0_dp
    bound(4) = 3.01_dp
    r = find_peaks(x, height_min_samples=bound)
    call check(r, [2, 6, 9, 12])
    if (size(r%peak_heights) /= 4) error stop 'height properties missing'

    bound = 100.0_dp
    bound(12) = 4.99_dp
    r = find_peaks(x, height_max_samples=bound)
    call check(r, [2, 4, 6, 9])

    bound = 0.0_dp
    bound(9) = 0.01_dp
    r = find_peaks(x, threshold_min_samples=bound)
    call check(r, [2, 4, 6, 12])
    if (size(r%left_thresholds) /= 4) error stop 'threshold properties missing'

    bound = 100.0_dp
    bound(6) = 3.99_dp
    r = find_peaks(x, threshold_max_samples=bound)
    call check(r, [2, 4, 9, 12])

    bound = 0.0_dp
    bound(6) = 4.01_dp
    r = find_peaks(x, prominence_min_samples=bound)
    call check(r, [2, 4, 9, 12])
    if (size(r%prominences) /= 4) error stop 'prominence properties missing'

    bound = 100.0_dp
    bound(12) = 4.99_dp
    r = find_peaks(x, prominence_max_samples=bound)
    call check(r, [2, 4, 6, 9])

    bound = 0.0_dp
    bound(9) = 3.01_dp
    r = find_peaks(x, width_min_samples=bound)
    call check(r, [2, 4, 6, 12])
    if (size(r%widths) /= 4) error stop 'width properties missing'

    bound = 100.0_dp
    bound(9) = 2.99_dp
    r = find_peaks(x, width_max_samples=bound)
    call check(r, [2, 4, 6, 12])

    bound = 0.0_dp
    bound(9) = 3.01_dp
    r = find_peaks(x, plateau_size_min_samples=bound)
    call check(r, [2, 4, 6, 12])
    if (size(r%plateau_sizes) /= 4) error stop 'plateau properties missing'

    bound = 100.0_dp
    bound(9) = 2.99_dp
    r = find_peaks(x, plateau_size_max_samples=bound)
    call check(r, [2, 4, 6, 12])

    ! Scalar and sample-wise bounds are intersected, without changing
    ! the behavior of callers that specify only the older scalar bounds.
    bound = 1.0_dp
    bound(6) = 4.0_dp
    r = find_peaks(x, height_min=3.0_dp, height_min_samples=bound)
    call check(r, [4, 6, 9, 12])

    bound = 4.0_dp
    r = find_peaks(x, height_min_samples=bound, width_min=0.0_dp, &
        plateau_size_min=1, prominence_min=0.0_dp)
    call check(r, [6, 12])
    if (size(r%widths) /= 2 .or. size(r%left_edges) /= 2) then
        error stop 'mixed scalar/array properties missing'
    end if

    ! Noncontiguous Fortran sections must be valid, including the bounds.
    strided_x = -999.0_dp
    strided_x(1::2) = x
    strided_bound = -999.0_dp
    strided_bound(1::2) = bound
    r = find_peaks(strided_x(1::2), height_min_samples=strided_bound(1::2))
    call check(r, [6, 12])

    r = find_peaks(x, height_min_samples=bound(1:13))
    if (r%status /= peak_status_invalid_input) error stop 'short bound accepted'
    if (size(r%peaks) /= 0) error stop 'invalid bound produced peaks'
    r = find_peaks(x, height_min_samples=[bound, 0.0_dp])
    if (r%status /= peak_status_invalid_input) error stop 'long bound accepted'

    nanvalue = ieee_value(1.0_dp, ieee_quiet_nan)
    bound(2) = nanvalue
    r = find_peaks(x, width_max_samples=bound)
    if (r%status /= peak_status_invalid_input) error stop 'nonfinite bound accepted'
    bound = 1.0_dp
    bound(9) = -0.01_dp
    r = find_peaks(x, plateau_size_min_samples=bound)
    if (r%status /= peak_status_invalid_input) error stop 'negative plateau bound accepted'

    r = find_peaks([real(dp) :: ], width_min_samples=[real(dp) :: ])
    call check(r, [integer :: ])
    r = find_peaks([real(dp) :: ], threshold_min_samples=[0.0_dp])
    if (r%status /= peak_status_invalid_input) error stop 'empty signal bound length'

    print '(a,i0,a)', 'test_signal_peaks_dynamic_bounds passed (', tested, ' cases)'

contains
    subroutine check(res, expected)
        type(peak_result), intent(in) :: res
        integer, intent(in) :: expected(:)

        tested = tested + 1
        if (res%status /= peak_status_success) error stop 'unexpected peak status'
        if (size(res%peaks) /= size(expected)) error stop 'peak count mismatch'
        if (any(res%peaks /= expected)) then
            print *, 'got:', res%peaks
            print *, 'expected:', expected
            error stop 'peak locations mismatch'
        end if
    end subroutine check
end program test_signal_peaks_dynamic_bounds
