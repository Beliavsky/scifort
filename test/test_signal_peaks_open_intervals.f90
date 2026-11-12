! SPDX-License-Identifier: MIT
! Open-interval properties are requested independently of filtering.
program test_signal_peaks_open_intervals
    use scifort_signal, only : dp, find_peaks, peak_result, peak_status_success
    implicit none
    real(dp) :: x(11), spaced(22)
    real(dp) :: limits(11)
    type(peak_result) :: res, ref
    integer :: mask, trials

    trials = 0
    x = [0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, &
         0.0_dp, 4.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    ref = find_peaks(x)
    call check_peaks(ref, [2, 5, 8, 10])
    call check_property_counts(ref, 0, 0, 0, 0, 0)
    do mask = 0, 31
        res = find_peaks(x, height_open=btest(mask, 0), &
            threshold_open=btest(mask, 1), prominence_open=btest(mask, 2), &
            width_open=btest(mask, 3), plateau_size_open=btest(mask, 4))
        call check_peaks(res, ref%peaks)
        call check_property_counts(res, merge(4, 0, btest(mask, 0)), &
            merge(4, 0, btest(mask, 1)), &
            merge(4, 0, btest(mask, 2) .or. btest(mask, 3)), &
            merge(4, 0, btest(mask, 3)), &
            merge(4, 0, btest(mask, 4)))
    end do
    res = find_peaks(x, height_open=.true., threshold_open=.true., &
        prominence_open=.true., width_open=.true., plateau_size_open=.true.)
    if (any(res%plateau_sizes /= [1, 3, 1, 1])) error stop 'plateau lengths'
    if (any(res%peak_heights /= [2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp])) error stop 'peak heights'
    if (any(res%left_thresholds /= [2.0_dp, 0.0_dp, 4.0_dp, 1.0_dp])) error stop 'left threshold'
    if (any(res%right_thresholds /= [2.0_dp, 0.0_dp, 4.0_dp, 1.0_dp])) error stop 'right threshold'
    if (any(abs(res%prominences - [2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp]) > 1.e-12_dp)) &
        error stop 'prominence values'
    if (any(res%left_edges /= [2, 4, 8, 10])) error stop 'left edge'
    if (any(res%right_edges /= [2, 6, 8, 10])) error stop 'right edge'
    if (any(abs(res%widths - [1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp]) > 1.e-12_dp)) &
        error stop 'width values'

    ! The open request does not disable restrictive bounds for the same property.
    limits = 0.0_dp
    limits(8) = 4.5_dp
    res = find_peaks(x, height_open=.true., height_min_samples=limits, &
        width_open=.true., plateau_size_open=.true.)
    call check_peaks(res, [2, 5, 10])
    call check_property_counts(res, 3, 0, 3, 3, 3)
    res = find_peaks(x, height_min=3.0_dp, prominence_open=.true., &
        threshold_open=.true.)
    call check_peaks(res, [5, 8])
    call check_property_counts(res, 2, 2, 2, 0, 0)
    ! Detrending, wlen and distance must not change in the property-only path.
    res = find_peaks(x, width_open=.true., distance=3, wlen=3, rel_height=1.0_dp)
    call check_peaks(res, [2, 5, 8])
    call check_property_counts(res, 0, 0, 3, 3, 0)

    ! A strided actual argument and empty input must both be supported.
    spaced = -10.0_dp
    spaced(1::2) = x
    res = find_peaks(spaced(1::2), width_open=.true.)
    call check_peaks(res, ref%peaks)
    res = find_peaks([real(dp) :: ], width_open=.true., plateau_size_open=.true.)
    call check_peaks(res, [integer :: ])
    call check_property_counts(res, 0, 0, 0, 0, 0)

    print '(a,i0)', 'test_signal_peaks_open_intervals passed: ', trials
contains
    subroutine check_peaks(r, expected)
        type(peak_result), intent(in) :: r
        integer, intent(in) :: expected(:)
        trials = trials + 1
        if (r%status /= peak_status_success) error stop 'bad status'
        if (size(r%peaks) /= size(expected)) error stop 'wrong peak count'
        if (any(r%peaks /= expected)) error stop 'wrong peak positions'
    end subroutine check_peaks

    subroutine check_property_counts(r, nh, nt, np, nw, ns)
        type(peak_result), intent(in) :: r
        integer, intent(in) :: nh, nt, np, nw, ns
        trials = trials + 1
        if (size(r%peak_heights) /= nh) error stop 'height count'
        if (size(r%left_thresholds) /= nt) error stop 'threshold count'
        if (size(r%right_thresholds) /= nt) error stop 'threshold right count'
        if (size(r%prominences) /= np) error stop 'prominence count'
        if (size(r%left_bases) /= np) error stop 'left base count'
        if (size(r%right_bases) /= np) error stop 'right base count'
        if (size(r%widths) /= nw) error stop 'width count'
        if (size(r%width_heights) /= nw) error stop 'width height count'
        if (size(r%left_ips) /= nw) error stop 'left ip count'
        if (size(r%right_ips) /= nw) error stop 'right ip count'
        if (size(r%plateau_sizes) /= ns) error stop 'plateau size count'
        if (size(r%left_edges) /= ns) error stop 'plateau left edge count'
        if (size(r%right_edges) /= ns) error stop 'plateau right edge count'
    end subroutine check_property_counts
end program test_signal_peaks_open_intervals
