! SPDX-License-Identifier: MIT
! Analog IIR coefficient references from SciPy 1.17.0, independent formulas.
program test_signal_analog_design
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, butter, cheby1, cheby2, iirfilter, &
        butter_result, freqs, frequency_response, design_status_success, &
        design_status_invalid_input
    implicit none
    type(butter_result) :: filt, wrapper
    type(frequency_response) :: response, expected
    integer :: i
    real(dp), parameter :: w(6) = [0.0_dp,0.2_dp,0.5_dp,1.0_dp,3.0_dp,7.0_dp]

    filt = butter(2, 2.0_dp, btype='lowpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'butter lowpass analog status'
    if (allocated(filt%sos)) error stop 'butter analog SOS should be unallocated'
    call compare(filt%b, [ &
            4.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 2.82842712474619029e+00_dp, 4.00000000000000000e+00_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, 2.0_dp, btype='lowpass', ftype='butter', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter butter'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = butter(2, 2.0_dp, btype='highpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'butter highpass analog status'
    if (allocated(filt%sos)) error stop 'butter analog SOS should be unallocated'
    call compare(filt%b, [ &
            1.00000000000000000e+00_dp, 0.00000000000000000e+00_dp, 0.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 2.82842712474619029e+00_dp, 4.00000000000000000e+00_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, 2.0_dp, btype='highpass', ftype='butter', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter butter'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = butter(2, [2.0_dp,5.0_dp], btype='bandpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'butter bandpass analog status'
    if (allocated(filt%sos)) error stop 'butter analog SOS should be unallocated'
    call compare(filt%b, [ &
            9.00000000000000000e+00_dp, 0.00000000000000000e+00_dp, 0.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 4.24264068711928566e+00_dp, 2.90000000000000071e+01_dp, &
            4.24264068711928672e+01_dp, 1.00000000000000043e+02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, [2.0_dp,5.0_dp], btype='bandpass', ftype='butter', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter butter'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = butter(2, [2.0_dp,5.0_dp], btype='bandstop', analog=.true.)
    if (filt%status /= design_status_success) error stop 'butter bandstop analog status'
    if (allocated(filt%sos)) error stop 'butter analog SOS should be unallocated'
    call compare(filt%b, [ &
            1.00000000000000000e+00_dp, 0.00000000000000000e+00_dp, 2.00000000000000036e+01_dp, &
            0.00000000000000000e+00_dp, 1.00000000000000028e+02_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 4.24264068711928566e+00_dp, 2.90000000000000071e+01_dp, &
            4.24264068711928672e+01_dp, 1.00000000000000043e+02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, [2.0_dp,5.0_dp], btype='bandstop', ftype='butter', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter butter'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby1(2, 0.5_dp, 2.0_dp, btype='lowpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby1 lowpass analog status'
    if (allocated(filt%sos)) error stop 'cheby1 analog SOS should be unallocated'
    call compare(filt%b, [ &
            5.72555032248637996e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 2.85124902728040430e+00_dp, 6.06481050778372577e+00_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, 2.0_dp, btype='lowpass', ftype='cheby1', &
        rp=0.5_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby1'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby1(2, 0.5_dp, 2.0_dp, btype='highpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby1 highpass analog status'
    if (allocated(filt%sos)) error stop 'cheby1 analog SOS should be unallocated'
    call compare(filt%b, [ &
            9.44060876285923278e-01_dp, 0.00000000000000000e+00_dp, 0.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 1.88051977790306357e+00_dp, 2.63816981247232896e+00_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, 2.0_dp, btype='highpass', ftype='cheby1', &
        rp=0.5_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby1'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby1(2, 0.5_dp, [2.0_dp,5.0_dp], btype='bandpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby1 bandpass analog status'
    if (allocated(filt%sos)) error stop 'cheby1 analog SOS should be unallocated'
    call compare(filt%b, [ &
            1.28824882255943542e+01_dp, 0.00000000000000000e+00_dp, 0.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 4.27687354092060623e+00_dp, 3.36458236425133848e+01_dp, &
            4.27687354092060588e+01_dp, 1.00000000000000014e+02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, [2.0_dp,5.0_dp], btype='bandpass', ftype='cheby1', &
        rp=0.5_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby1'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby1(2, 0.5_dp, [2.0_dp,5.0_dp], btype='bandstop', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby1 bandstop analog status'
    if (allocated(filt%sos)) error stop 'cheby1 analog SOS should be unallocated'
    call compare(filt%b, [ &
            9.44060876285923278e-01_dp, 0.00000000000000000e+00_dp, 1.88812175257184691e+01_dp, &
            0.00000000000000000e+00_dp, 9.44060876285923598e+01_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 2.82077966685459502e+00_dp, 2.59358820780627433e+01_dp, &
            2.82077966685459494e+01_dp, 1.00000000000000028e+02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, [2.0_dp,5.0_dp], btype='bandstop', ftype='cheby1', &
        rp=0.5_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby1'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby2(2, 40.0_dp, 2.0_dp, btype='lowpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby2 lowpass analog status'
    if (allocated(filt%sos)) error stop 'cheby2 analog SOS should be unallocated'
    call compare(filt%b, [ &
            9.99999999999999327e-03_dp, 0.00000000000000000e+00_dp, 7.99999999999999600e-02_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 3.97994974842647953e-01_dp, 7.99999999999999600e-02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, 2.0_dp, btype='lowpass', ftype='cheby2', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby2'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby2(2, 40.0_dp, 2.0_dp, btype='highpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby2 highpass analog status'
    if (allocated(filt%sos)) error stop 'cheby2 analog SOS should be unallocated'
    call compare(filt%b, [ &
            1.00000000000000022e+00_dp, 0.00000000000000000e+00_dp, 2.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 1.98997487421324060e+01_dp, 2.00000000000000057e+02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, 2.0_dp, btype='highpass', ftype='cheby2', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby2'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby2(2, 40.0_dp, [2.0_dp,5.0_dp], btype='bandpass', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby2 bandpass analog status'
    if (allocated(filt%sos)) error stop 'cheby2 analog SOS should be unallocated'
    call compare(filt%b, [ &
            9.99999999999999327e-03_dp, 0.00000000000000000e+00_dp, 3.79999999999999893e-01_dp, &
            0.00000000000000000e+00_dp, 1.00000000000000000e+00_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 5.96992462263971957e-01_dp, 2.01800000000000033e+01_dp, &
            5.96992462263972090e+00_dp, 1.00000000000000043e+02_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, [2.0_dp,5.0_dp], btype='bandpass', ftype='cheby2', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby2'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    filt = cheby2(2, 40.0_dp, [2.0_dp,5.0_dp], btype='bandstop', analog=.true.)
    if (filt%status /= design_status_success) error stop 'cheby2 bandstop analog status'
    if (allocated(filt%sos)) error stop 'cheby2 analog SOS should be unallocated'
    call compare(filt%b, [ &
            1.00000000000000022e+00_dp, 0.00000000000000000e+00_dp, 2.45000000000000071e+01_dp, &
            0.00000000000000000e+00_dp, 1.00000000000000071e+02_dp ], 6.0e-12_dp)
    call compare(filt%a, [ &
            1.00000000000000000e+00_dp, 2.98496231131986107e+01_dp, 4.70000000000000114e+02_dp, &
            2.98496231131985894e+02_dp, 9.99999999999998295e+01_dp ], 6.0e-12_dp)
    wrapper = iirfilter(2, [2.0_dp,5.0_dp], btype='bandstop', ftype='cheby2', &
        rp=40.0_dp, rs=40.0_dp, analog=.true.)
    if (wrapper%status /= design_status_success) error stop 'iirfilter cheby2'
    call compare(filt%b, wrapper%b, 1.0e-14_dp)
    call compare(filt%a, wrapper%a, 1.0e-14_dp)

    ! An odd-order inverse Chebyshev prototype has one zero at infinity.
    filt = cheby2(3, 40.0_dp, 5.0_dp, analog=.true.)
    if (filt%status /= design_status_success) error stop 'odd Chebyshev II'
    if (size(filt%b) /= 3 .or. size(filt%a) /= 4) error stop 'Cheby2 degrees'
    response = freqs(filt%b, filt%a, w=w)
    if (response%status /= design_status_success) error stop 'analog freqs interoperability'
    if (abs(response%response(1) - cmplx(1.0_dp, 0.0_dp, dp)) > 1.0e-12_dp) &
        error stop 'Cheby2 unity DC gain'

    ! Butterworth cutoff response is 1/sqrt(2) of DC for all orders.
    do i = 1, 8
        filt = butter(i, 3.0_dp, analog=.true.)
        if (filt%status /= design_status_success) error stop 'Butterworth order'
        response = freqs(filt%b, filt%a, w=[0.0_dp,3.0_dp])
        if (response%status /= design_status_success) error stop 'Butterworth freqs'
        if (abs(abs(response%response(2)) - sqrt(0.5_dp)) > 1.0e-9_dp) &
            error stop 'Butterworth -3 dB cutoff'
    end do

    ! Previously supported digital API is unaffected by the analog option.
    filt = butter(3, 0.3_dp)
    wrapper = butter(3, 0.3_dp, analog=.false.)
    if (filt%status /= design_status_success) error stop 'digital default'
    if (.not. allocated(filt%sos)) error stop 'digital SOS missing'
    call compare(filt%b, wrapper%b, 0.0_dp)
    call compare(filt%a, wrapper%a, 0.0_dp)
    filt = cheby2(3, 35.0_dp, 0.4_dp)
    wrapper = cheby2(3, 35.0_dp, 0.4_dp, analog=.false.)
    if (filt%status /= design_status_success) error stop 'digital Cheby2 default'
    call compare(filt%a, wrapper%a, 0.0_dp)

    ! Reject invalid input without allocating misleading polynomial outputs.
    filt = butter(3, -1.0_dp, analog=.true.)
    call invalid(filt)
    filt = butter(3, 4.0_dp, analog=.true., fs=500.0_dp)
    call invalid(filt)
    filt = cheby1(3, 0.0_dp, 5.0_dp, analog=.true.)
    call invalid(filt)
    filt = cheby2(3, 40.0_dp, [2.0_dp,2.0_dp], btype='bandpass', analog=.true.)
    call invalid(filt)
    filt = butter(21, 2.0_dp, analog=.true.)
    call invalid(filt)
    filt = iirfilter(3, 2.0_dp, ftype='unknown', analog=.true.)
    call invalid(filt)
    filt = butter(3, ieee_value(1.0_dp,ieee_quiet_nan), analog=.true.)
    call invalid(filt)

    print '(a)', 'test_signal_analog_design passed'
contains
    subroutine compare(actual, expected, tol)
        real(dp), intent(in) :: actual(:) !! actual analog coefficients
        real(dp), intent(in) :: expected(:) !! SciPy reference coefficients
        real(dp), intent(in) :: tol !! relative and absolute tolerance
        if (size(actual) /= size(expected)) error stop 'analog coefficient count'
        if (any(abs(actual-expected) > tol * max(1.0_dp,abs(expected)))) then
            print *, 'actual:',actual
            print *, 'expected:',expected
            error stop 'analog coefficient mismatch'
        end if
    end subroutine compare
    subroutine invalid(f)
        type(butter_result), intent(in) :: f !! filter expected to fail validation
        if (f%status /= design_status_invalid_input) error stop 'invalid analog status'
        if (allocated(f%b) .or. allocated(f%a)) error stop 'invalid analog allocations'
    end subroutine invalid
end program test_signal_analog_design
