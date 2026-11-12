! SPDX-License-Identifier: MIT
! Verify the new uncertainty API through the public scifort_signal module.
program test_signal_harmonic_uncertainty_public
    use scifort_signal, only : dp, harmonic_fit_result, lombscargle_harmonic_fit, &
        lombscargle_harmonic_predict_std, irregular_status_success
    implicit none

    real(dp) :: t(57), y(57)
    real(dp), allocatable :: se(:)
    type(harmonic_fit_result) :: fit
    integer :: j, status

    do j = 1, size(t)
        t(j) = 0.17_dp*real(j,dp)+0.008_dp*sin(real(j*j,dp))
    end do
    y = 0.4_dp+1.1_dp*cos(0.94_dp*t)-0.2_dp*sin(0.94_dp*t) + &
        0.3_dp*cos(1.88_dp*t)+0.25_dp*sin(1.88_dp*t) + 0.04_dp*sin(3.7_dp*t)
    fit = lombscargle_harmonic_fit(t,y,[0.94_dp],2, &
        floating_mean=.true.,estimate_covariance=.true.)
    if (fit%status /= irregular_status_success) error stop 'public covariance fit'
    if (any(shape(fit%covariance) /= [5,5,1])) error stop 'public covariance shape'
    se = lombscargle_harmonic_predict_std(fit,[0.0_dp,0.5_dp,2.7_dp],1,status)
    if (status /= irregular_status_success) error stop 'public prediction SE status'
    if (any(se <= 0.0_dp)) error stop 'public prediction uncertainty'
    print '(a)', 'test_signal_harmonic_uncertainty_public passed'
end program test_signal_harmonic_uncertainty_public
