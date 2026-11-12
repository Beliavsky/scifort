! SPDX-License-Identifier: MIT
! Verify that harmonic fitting and forecasting are exposed by scifort_signal.
program test_signal_harmonic_fit_public
    use scifort_signal, only : dp, harmonic_fit_result, lombscargle_harmonic_fit, &
        lombscargle_harmonic_predict, irregular_status_success
    implicit none
    real(dp) :: x(45), y(45)
    real(dp), allocatable :: predicted(:)
    type(harmonic_fit_result) :: fit
    integer :: i, status

    do i = 1, size(x)
        x(i) = 0.19_dp*real(i,dp)+0.03_dp*sin(real(i*i,dp))
    end do
    y = 1.8_dp + 0.5_dp*cos(1.24_dp*x)-0.3_dp*sin(1.24_dp*x) &
        + 0.27_dp*cos(2.48_dp*x)+0.12_dp*sin(2.48_dp*x)
    fit = lombscargle_harmonic_fit(x,y,[1.24_dp],2,floating_mean=.true.)
    if (fit%status /= irregular_status_success) error stop 'public fit status'
    predicted = lombscargle_harmonic_predict(fit,x,1,status)
    if (status /= irregular_status_success) error stop 'public prediction status'
    if (maxval(abs(predicted-y)) > 2.0e-12_dp) error stop 'public reconstruction'
    print '(a)', 'test_signal_harmonic_fit_public passed'
end program test_signal_harmonic_fit_public
