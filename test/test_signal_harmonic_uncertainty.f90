! SPDX-License-Identifier: MIT
! Weighted coefficient covariance, prediction SE, and input validation.
program test_signal_harmonic_uncertainty
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : harmonic_fit_result, lombscargle_harmonic_fit, &
        lombscargle_harmonic_predict_std, irregular_status_success, irregular_status_singular, &
        irregular_status_invalid_input
    implicit none

    integer, parameter :: n = 75
    real(dp), parameter :: omega = 0.83_dp
    real(dp) :: time(n), observation(n), weight(n), strided(2*n)
    real(dp) :: predict_times(4), nan_value
    real(dp), allocatable :: se(:), other_se(:)
    type(harmonic_fit_result) :: fit, another, no_intercept, without_covariance
    integer :: j, status

    do j = 1, n
        time(j) = -2.0_dp + 0.173_dp*real(j,dp) + 0.013_dp*cos(real(j*j,dp))
        weight(j) = 0.15_dp + real(mod(j,11),dp)/7.0_dp
        observation(j) = 1.3_dp + 0.7_dp*cos(omega*time(j)) - 0.4_dp*sin(omega*time(j)) + &
            0.31_dp*cos(2.0_dp*omega*time(j)) + 0.23_dp*sin(2.0_dp*omega*time(j)) + &
            0.05_dp*cos(real(4*j*j,dp))
    end do
    weight(3:6) = 0.0_dp
    predict_times = [-5.0_dp,0.0_dp,3.0_dp,17.0_dp]
    fit = lombscargle_harmonic_fit(time,observation,[0.0_dp,omega,1.39_dp],2, &
        floating_mean=.true., weights=weight, estimate_covariance=.true.)
    if (fit%status /= irregular_status_singular) error stop 'partially singular'
    if (any(fit%valid .neqv. [.false.,.true.,.true.])) error stop 'valid bins'
    if (.not. allocated(fit%covariance)) error stop 'covariance not allocated'
    if (any(shape(fit%covariance) /= [5,5,3])) error stop 'covariance shape'
    if (any(shape(fit%standard_error) /= [5,3])) error stop 'SE shape'
    if (fit%degrees_of_freedom /= n-4-5) error stop 'degrees of freedom'
    if (.not. ieee_is_nan(fit%covariance(1,1,1))) error stop 'invalid covariance'
    if (any(fit%standard_error(:,2) <= 0.0_dp)) error stop 'positive SE'
    if (abs(fit%covariance(1,3,2)-fit%covariance(3,1,2)) > 1.0e-14_dp) &
        error stop 'covariance symmetry'
    if (abs(fit%standard_error(2,2)**2-fit%covariance(2,2,2)) > 1.0e-13_dp) &
        error stop 'variance / standard error'
    if (fit%residual_variance(2) <= 0.0_dp) error stop 'residual noise scale'
    se = lombscargle_harmonic_predict_std(fit,predict_times,2,status)
    if (status /= irregular_status_success) error stop 'prediction SE status'
    if (any(se <= 0.0_dp)) error stop 'prediction SE positive'
    other_se = lombscargle_harmonic_predict_std(fit,predict_times,1,status)
    if (status /= irregular_status_invalid_input) error stop 'singular SE status'
    if (.not. all(ieee_is_nan(other_se))) error stop 'singular SE values'
    other_se = lombscargle_harmonic_predict_std(fit,predict_times,8,status)
    if (status /= irregular_status_invalid_input) error stop 'invalid index'

    another = lombscargle_harmonic_fit(time,observation,[omega],2, &
        floating_mean=.true., weights=1.0e7_dp*weight, estimate_covariance=.true.)
    if (another%status /= irregular_status_success) error stop 'scaled weight status'
    if (maxval(abs(another%covariance(:,:,1)-fit%covariance(:,:,2))) > 2.0e-13_dp) &
        error stop 'scale invariant covariance'
    other_se = lombscargle_harmonic_predict_std(another,predict_times,1,status)
    if (maxval(abs(other_se-se)) > 2.0e-13_dp) error stop 'scale invariant prediction SE'

    strided(1:2*n:2) = time
    strided(2:2*n:2) = observation
    another = lombscargle_harmonic_fit(strided(1:2*n:2),strided(2:2*n:2), &
        [omega],2, floating_mean=.true., weights=weight, estimate_covariance=.true.)
    if (maxval(abs(another%covariance(:,:,1)-fit%covariance(:,:,2))) > 2.0e-13_dp) &
        error stop 'strided data covariance'

    no_intercept = lombscargle_harmonic_fit(time,observation,[omega],2, &
        precenter=.true.,floating_mean=.false.,estimate_covariance=.true.)
    if (no_intercept%status /= irregular_status_success) error stop 'fixed intercept status'
    if (.not. ieee_is_nan(no_intercept%standard_error(1,1))) &
        error stop 'fixed offset must have undefined standard error'
    if (.not. all(ieee_is_nan(no_intercept%covariance(1,:,1)))) &
        error stop 'fixed offset covariance row'
    if (any(no_intercept%standard_error(2:,1) <= 0.0_dp)) &
        error stop 'harmonic SE without floating mean'
    other_se = lombscargle_harmonic_predict_std(no_intercept,predict_times,1,status)
    if (status /= irregular_status_success .or. any(other_se <= 0.0_dp)) &
        error stop 'fixed-offset prediction uncertainty'

    without_covariance = lombscargle_harmonic_fit(time,observation,[omega],2)
    if (allocated(without_covariance%covariance)) error stop 'default must be allocation-free'
    other_se = lombscargle_harmonic_predict_std(without_covariance,predict_times,1,status)
    if (status /= irregular_status_invalid_input) error stop 'missing covariance'
    if (.not. all(ieee_is_nan(other_se))) error stop 'missing covariance values'
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    other_se = lombscargle_harmonic_predict_std(fit,[0.0_dp,nan_value],2,status)
    if (status /= irregular_status_invalid_input) error stop 'nonfinite time status'
    if (.not. all(ieee_is_nan(other_se))) error stop 'nonfinite time values'
    other_se = lombscargle_harmonic_predict_std(fit,[huge(1.0_dp)],3,status)
    if (status /= irregular_status_invalid_input) error stop 'overflow time status'

    another = lombscargle_harmonic_fit(time,observation,[omega],2, &
        weights=-weight,estimate_covariance=.true.)
    if (another%status /= irregular_status_invalid_input) error stop 'negative weights'
    if (.not. all(ieee_is_nan(another%covariance))) error stop 'invalid covariance NaNs'
    another = lombscargle_harmonic_fit(time,observation,[omega],huge(1), &
        estimate_covariance=.true.)
    if (another%status /= irregular_status_invalid_input) error stop 'huge harmonic count'
    if (size(another%covariance,1) /= 1) error stop 'huge coefficient allocation'

    print '(a)', 'test_signal_harmonic_uncertainty passed'
end program test_signal_harmonic_uncertainty
