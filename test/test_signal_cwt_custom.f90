! SPDX-License-Identifier: MIT
program test_signal_cwt_custom
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : find_peaks_cwt, cwt_ricker, cwt_wavelet, &
        ricker_wavelet, cwt_peak_success, cwt_peak_invalid_input
    implicit none
    real(dp) :: x(127), widths(9), t
    real(dp), allocatable :: baseline(:,:), custom(:,:), shifted(:,:), strided(:,:)
    real(dp), allocatable :: xcopy(:), wcopy(:)
    integer, allocatable :: defaults(:), peaks(:), strided_peaks(:), explicit_dist(:)
    integer :: i, status

    widths = [0.75_dp, 1.25_dp, 2.0_dp, 3.5_dp, 5.0_dp, 6.5_dp, 8.0_dp, 10.0_dp, 12.0_dp]
    do i = 1, size(x)
        t = real(i, dp)
        x(i) = 3.0_dp*exp(-0.5_dp*((t-39.0_dp)/4.0_dp)**2) + &
            4.0_dp*exp(-0.5_dp*((t-91.0_dp)/6.0_dp)**2) + &
            0.02_dp*sin(0.9_dp*t)
    end do
    baseline = cwt_ricker(x, widths, status)
    if (status /= cwt_peak_success) error stop 'baseline CWT status'
    custom = cwt_wavelet(x, widths, ricker_wavelet, status)
    if (status /= cwt_peak_success) error stop 'callback CWT status'
    if (any(shape(custom) /= shape(baseline))) error stop 'callback CWT shape'
    if (maxval(abs(custom-baseline)) > 1.0e-13_dp) error stop 'Ricker callback CWT'
    defaults = find_peaks_cwt(x, widths, status=status)
    if (status /= cwt_peak_success) error stop 'baseline peaks status'
    peaks = find_peaks_cwt(x, widths, status=status, wavelet=ricker_wavelet)
    if (status /= cwt_peak_success) error stop 'callback peaks status'
    if (size(peaks) /= size(defaults)) error stop 'callback peaks shape'
    if (any(peaks /= defaults)) error stop 'callback Ricker peaks'

    custom = cwt_wavelet(x, widths, asymmetric_wavelet, status)
    if (status /= cwt_peak_success) error stop 'asymmetric CWT status'
    shifted = cwt_wavelet(x, widths, gaussian_wavelet, status)
    if (status /= cwt_peak_success) error stop 'Gaussian CWT status'
    if (maxval(abs(shifted-custom)) < 1.0e-6_dp) error stop 'wavelets ignored'
    peaks = find_peaks_cwt(x, widths, wavelet=gaussian_wavelet, status=status)
    if (status /= cwt_peak_success) error stop 'Gaussian peaks status'
    if (any(peaks < 1) .or. any(peaks > size(x))) error stop 'Gaussian peak range'
    explicit_dist = find_peaks_cwt(x, widths, max_distances=widths/4.0_dp, &
        wavelet=gaussian_wavelet, status=status)
    if (status /= cwt_peak_success) error stop 'Gaussian distance status'
    if (size(explicit_dist) /= size(peaks)) error stop 'explicit distance shape'
    if (any(explicit_dist /= peaks)) error stop 'explicit default distance'

    ! Noncontiguous data and widths must match contiguous operation.
    xcopy = x(1::2)
    wcopy = widths(1::2)
    strided = cwt_wavelet(x(1::2), widths(1::2), asymmetric_wavelet, status)
    custom = cwt_wavelet(xcopy, wcopy, asymmetric_wavelet, status)
    if (status /= cwt_peak_success) error stop 'strided CWT status'
    if (maxval(abs(strided-custom)) > 1.0e-13_dp) error stop 'strided CWT'
    strided_peaks = find_peaks_cwt(x(1::2), widths(1::2), wavelet=asymmetric_wavelet, &
        status=status)
    if (status /= cwt_peak_success) error stop 'strided peak status'
    if (any(strided_peaks < 1) .or. any(strided_peaks > size(x(1::2)))) &
        error stop 'strided peak range'

    custom = cwt_wavelet(x, widths, short_wavelet, status)
    if (status /= cwt_peak_invalid_input) error stop 'short kernel status'
    if (any(shape(custom) /= 0)) error stop 'short kernel output'
    peaks = find_peaks_cwt(x, widths, wavelet=short_wavelet, status=status)
    if (status /= cwt_peak_invalid_input .or. size(peaks) /= 0) &
        error stop 'short peak result'
    custom = cwt_wavelet(x, widths, nonfinite_wavelet, status)
    if (status /= cwt_peak_invalid_input) error stop 'nan kernel status'
    if (any(shape(custom) /= 0)) error stop 'nan kernel output'
    peaks = find_peaks_cwt(x, widths, wavelet=nonfinite_wavelet, status=status)
    if (status /= cwt_peak_invalid_input .or. size(peaks) /= 0) &
        error stop 'nan peak result'
    custom = cwt_wavelet(x, widths, late_bad_wavelet, status)
    if (status /= cwt_peak_invalid_input) error stop 'late kernel status'
    if (any(shape(custom) /= 0)) error stop 'late kernel output'
    peaks = find_peaks_cwt(x, widths, wavelet=late_bad_wavelet, status=status)
    if (status /= cwt_peak_invalid_input .or. size(peaks) /= 0) &
        error stop 'late peak result'
    custom = cwt_wavelet(x, [-2.0_dp], gaussian_wavelet, status)
    if (status /= cwt_peak_invalid_input) error stop 'invalid widths status'
    if (any(shape(custom) /= 0)) error stop 'invalid widths output'
    print *, 'signal custom CWT wavelets passed'
contains
    function gaussian_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested wavelet sample count
        real(dp), intent(in) :: width !! positive scale
        real(dp), allocatable :: w(:) !! sampled Gaussian wavelet
        real(dp) :: u
        integer :: k
        allocate(w(ceiling(points)))
        do k = 1, size(w)
            u = (real(k-1,dp) - (points-1.0_dp)/2.0_dp)/width
            w(k) = exp(-0.5_dp*u*u)/(sqrt(2.0_dp*acos(-1.0_dp))*width)
        end do
    end function gaussian_wavelet

    function asymmetric_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested wavelet sample count
        real(dp), intent(in) :: width !! positive scale
        real(dp), allocatable :: w(:) !! sampled shifted Gaussian wavelet
        real(dp) :: u
        integer :: k
        allocate(w(ceiling(points)))
        do k = 1, size(w)
            u = (real(k-1,dp) - (points-1.0_dp)/2.0_dp)/width - 0.35_dp
            w(k) = exp(-0.5_dp*u*u)/(sqrt(2.0_dp*acos(-1.0_dp))*width)
        end do
    end function asymmetric_wavelet

    function short_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested wavelet sample count
        real(dp), intent(in) :: width !! positive scale
        real(dp), allocatable :: w(:) !! deliberately malformed length
        allocate(w(max(0,ceiling(points)-1)))
        w = width
    end function short_wavelet

    function late_bad_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested wavelet sample count
        real(dp), intent(in) :: width !! positive scale
        real(dp), allocatable :: w(:) !! invalid only after several valid widths
        if (width > 4.0_dp) then
            allocate(w(0))
        else
            w = gaussian_wavelet(points, width)
        end if
    end function late_bad_wavelet

    function nonfinite_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! requested wavelet sample count
        real(dp), intent(in) :: width !! positive scale
        real(dp), allocatable :: w(:) !! deliberately nonfinite samples
        allocate(w(ceiling(points)))
        w = width
        w(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nonfinite_wavelet
end program test_signal_cwt_custom
