! SPDX-License-Identifier: MIT
! Test harmonic amplitudes, independent frequency failure, and prediction.
program test_signal_harmonic_fit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : harmonic_fit_result, lombscargle_harmonic_fit, &
        lombscargle_harmonic_predict, lombscargle_multiharmonic, irregular_status_success, &
        irregular_status_singular, irregular_status_invalid_input
    implicit none
    integer, parameter :: n = 79
    real(dp), parameter :: omega = 1.17_dp
    real(dp) :: x(n), y(n), w(n), future(7), expected(7), shifted(2*n)
    real(dp), allocatable :: prediction(:), reference(:)
    type(harmonic_fit_result) :: fit, shifted_fit
    integer :: i, st, h

    do i = 1, n
        x(i) = -2.0_dp + 0.14_dp*real(i,dp) + 0.019_dp*sin(real(i*i,dp))
        w(i) = 0.01_dp+real(mod(i,9),dp)/7.0_dp
    end do
    w(1:5) = 0.0_dp
    y = 2.3_dp + 1.4_dp*cos(omega*x) - 0.35_dp*sin(omega*x) &
        + 0.7_dp*cos(2.0_dp*omega*x) + 0.9_dp*sin(2.0_dp*omega*x) &
        - 0.42_dp*cos(3.0_dp*omega*x) + 0.18_dp*sin(3.0_dp*omega*x)
    fit = lombscargle_harmonic_fit(x, y, [0.0_dp,omega,1.61_dp], 3, &
        floating_mean=.true., weights=w)
    if (fit%status /= irregular_status_singular) error stop 'partially singular status'
    if (any(fit%valid .neqv. [.false.,.true.,.true.])) error stop 'partial validity'
    if (.not. ieee_is_nan(fit%cosine(1,1))) error stop 'singular coefficient NaN'
    if (abs(fit%offset(2)-2.3_dp) > 3.0e-12_dp) error stop 'intercept'
    if (maxval(abs(fit%cosine(:,2)-[1.4_dp,0.7_dp,-0.42_dp])) > 3.0e-12_dp) &
        error stop 'cosine coefficients'
    if (maxval(abs(fit%sine(:,2)-[-0.35_dp,0.9_dp,0.18_dp])) > 3.0e-12_dp) &
        error stop 'sine coefficients'
    if (fit%rms(2) > 3.0e-12_dp) error stop 'residual RMS'
    if (abs(fit%normalized_power(2)-1.0_dp) > 3.0e-12_dp) error stop 'normalized power'
    reference = lombscargle_multiharmonic(x, y, [omega,1.61_dp], 3, &
        floating_mean=.true., weights=w)
    if (maxval(abs(fit%power(2:)-reference)) > 5.0e-12_dp) error stop 'periodogram comparison'
    future = [-6.0_dp, -1.0_dp, 0.0_dp, 0.2_dp, 9.0_dp, 12.0_dp, 20.0_dp]
    expected = 2.3_dp + 1.4_dp*cos(omega*future) - 0.35_dp*sin(omega*future) &
        + 0.7_dp*cos(2.0_dp*omega*future) + 0.9_dp*sin(2.0_dp*omega*future) &
        - 0.42_dp*cos(3.0_dp*omega*future) + 0.18_dp*sin(3.0_dp*omega*future)
    prediction = lombscargle_harmonic_predict(fit, future, 2, st)
    if (st /= irregular_status_success) error stop 'prediction status'
    if (maxval(abs(prediction-expected)) > 3.0e-12_dp) error stop 'extrapolation'
    prediction = lombscargle_harmonic_predict(fit, future, 1, st)
    if (st /= irregular_status_invalid_input .or. .not. all(ieee_is_nan(prediction))) &
        error stop 'singular prediction rejected'
    prediction = lombscargle_harmonic_predict(fit, future, 7, st)
    if (st /= irregular_status_invalid_input .or. .not. all(ieee_is_nan(prediction))) &
        error stop 'invalid frequency index'

    ! Precentering must NOT change the returned model on the original data scale.
    shifted_fit = lombscargle_harmonic_fit(x, y, [omega], 3, precenter=.true., &
        floating_mean=.true., weights=w)
    if (shifted_fit%status /= irregular_status_success) error stop 'precentered fit'
    if (abs(shifted_fit%offset(1)-fit%offset(2)) > 3.0e-12_dp) error stop 'precentered offset'
    if (maxval(abs(shifted_fit%cosine(:,1)-fit%cosine(:,2))) > 3.0e-12_dp) &
        error stop 'precentered coefficients'
    if (maxval(abs(shifted_fit%sine(:,1)-fit%sine(:,2))) > 3.0e-12_dp) &
        error stop 'precentered sine'
    shifted(1:2*n:2) = x
    shifted(2:2*n:2) = y
    shifted_fit = lombscargle_harmonic_fit(shifted(1:2*n:2), shifted(2:2*n:2), [omega], &
        3, floating_mean=.true., weights=w)
    if (maxval(abs(shifted_fit%cosine(:,1)-fit%cosine(:,2))) > 3.0e-12_dp) &
        error stop 'strided inputs'

    ! A one-term model is also supported; a weighted intercept is optional.
    shifted_fit = lombscargle_harmonic_fit(x, y, [omega], 1, floating_mean=.false.)
    if (shifted_fit%status /= irregular_status_success) error stop 'one-harmonic fit'
    reference = lombscargle_multiharmonic(x, y, [omega], 1)
    if (abs(shifted_fit%power(1)-reference(1)) > 5.0e-12_dp) error stop 'one harmonic power'

    fit = lombscargle_harmonic_fit(x, y, [omega], huge(1), weights=w)
    if (fit%status /= irregular_status_invalid_input) error stop 'huge harmonic count'
    if (size(fit%cosine,1) /= 0) error stop 'huge count must not allocate'
    fit = lombscargle_harmonic_fit(x, y, [omega], 0)
    if (fit%status /= irregular_status_invalid_input) error stop 'zero harmonic count'
    fit = lombscargle_harmonic_fit(x, y, [omega], 3, weights=-w)
    if (fit%status /= irregular_status_invalid_input) error stop 'negative weight'
    fit = lombscargle_harmonic_fit(x, y, [omega], 3, weights=w(:n-1))
    if (fit%status /= irregular_status_invalid_input) error stop 'weight dimensions'
    fit = lombscargle_harmonic_fit(x, y(:10), [omega], 3)
    if (fit%status /= irregular_status_invalid_input) error stop 'signal dimensions'
    fit = lombscargle_harmonic_fit(x, 1.0_dp+0.0_dp*y, [omega], 3, &
        floating_mean=.true.)
    if (fit%status /= irregular_status_singular) error stop 'constant signal'
    if (any(fit%valid)) error stop 'constant signal valid'
    print '(a)', 'test_signal_harmonic_fit passed'
end program test_signal_harmonic_fit
