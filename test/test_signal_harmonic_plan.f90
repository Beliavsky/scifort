! SPDX-License-Identifier: MIT
! Prepared multiharmonic weighted QR reused across independent channels.
program test_signal_harmonic_plan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : harmonic_fit_plan, harmonic_fit_result, &
        lombscargle_harmonic_prepare, lombscargle_harmonic_fit_prepared, &
        lombscargle_harmonic_fit_batch, lombscargle_harmonic_fit, &
        lombscargle_harmonic_predict_std, irregular_status_success, &
        irregular_status_singular, irregular_status_invalid_input
    implicit none

    integer, parameter :: n = 83, nch = 4, nh = 3, nf = 4
    real(dp) :: t(n), y(n,nch), w(n), input(2*n,2*nch), freqs(nf)
    real(dp) :: future(5), nan, maxdiff
    real(dp), allocatable :: prediction_std(:), direct_std(:)
    type(harmonic_fit_plan) :: plan, invalid, other
    type(harmonic_fit_result) :: ref, single
    type(harmonic_fit_result), allocatable :: fits(:), strided_fits(:), invalid_fits(:)
    integer :: i, j, k, st

    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    do i = 1, n
        t(i) = -3.0_dp+0.139_dp*real(i,dp)+0.027_dp*sin(real(i*i,dp))
        w(i) = 0.2_dp+real(mod(i,13),dp)/6.0_dp
        do j = 1, nch
            y(i,j) = 0.4_dp*real(j,dp)+0.32_dp*sin(real(i*j,dp))
            do k = 1, nh
                y(i,j) = y(i,j)+0.21_dp*real(k*j,dp)*cos(real(k,dp)*0.93_dp*t(i))+ &
                    0.13_dp*real(j-k,dp)*sin(real(k,dp)*0.93_dp*t(i))
            end do
        end do
    end do
    w(2:6) = 0.0_dp
    freqs = [0.0_dp,0.93_dp,1.41_dp,2.13_dp]
    future = [-10.0_dp,-1.0_dp,0.0_dp,3.0_dp,14.0_dp]
    plan = lombscargle_harmonic_prepare(t,freqs,nh,weights=w,floating_mean=.true.)
    if (plan%status /= irregular_status_singular) error stop 'partly singular plan'
    if (any(plan%valid .neqv. [.false.,.true.,.true.,.true.])) error stop 'plan valid bins'
    fits = lombscargle_harmonic_fit_batch(plan,y,precenter=.true.,estimate_covariance=.true.)
    if (size(fits) /= nch) error stop 'batch channel count'
    do j = 1, nch
        ref = lombscargle_harmonic_fit(t,y(:,j),freqs,nh,precenter=.true., &
            weights=w,floating_mean=.true.,estimate_covariance=.true.)
        single = lombscargle_harmonic_fit_prepared(plan,y(:,j),precenter=.true., &
            estimate_covariance=.true.)
        if (fits(j)%status /= irregular_status_singular) error stop 'batch singular bin status'
        if (any(fits(j)%valid .neqv. ref%valid)) error stop 'batch valid bins'
        if (.not. ieee_is_nan(fits(j)%cosine(1,1))) error stop 'bad bin NaN'
        if (.not. ieee_is_nan(fits(j)%covariance(1,1,1))) error stop 'bad covariance NaN'
        if (fits(j)%degrees_of_freedom /= count(w>0.0_dp)-2*nh-1) &
            error stop 'batch dof'
        do k = 2, nf
            if (maxval(abs(fits(j)%cosine(:,k)-ref%cosine(:,k))) > 1.0e-10_dp) &
                error stop 'cosine versus direct'
            if (maxval(abs(fits(j)%sine(:,k)-ref%sine(:,k))) > 1.0e-10_dp) &
                error stop 'sine versus direct'
            if (abs(fits(j)%offset(k)-ref%offset(k)) > 1.0e-10_dp) error stop 'offset'
            if (abs(fits(j)%power(k)-ref%power(k)) > 1.0e-10_dp) error stop 'power'
            if (abs(fits(j)%normalized_power(k)-ref%normalized_power(k)) > 1.0e-10_dp) &
                error stop 'normalized power'
            if (abs(fits(j)%rms(k)-ref%rms(k)) > 1.0e-10_dp) error stop 'rms'
            maxdiff = maxval(abs(fits(j)%covariance(:,:,k)-ref%covariance(:,:,k)))
            if (maxdiff > 1.0e-10_dp) error stop 'covariance'
            if (maxval(abs(single%standard_error(:,k)-fits(j)%standard_error(:,k))) > 1.0e-12_dp) &
                error stop 'single/batch'
            prediction_std = lombscargle_harmonic_predict_std(fits(j),future,k,st)
            direct_std = lombscargle_harmonic_predict_std(ref,future,k,st)
            if (st /= irregular_status_success) error stop 'prediction std status'
            if (maxval(abs(prediction_std-direct_std)) > 1.0e-10_dp) &
                error stop 'prediction std'
        end do
    end do
    ! A reused plan does not accumulate information from previous observations.
    single = lombscargle_harmonic_fit_prepared(plan,y(:,1),estimate_covariance=.false.)
    if (allocated(single%covariance)) error stop 'optional covariance allocation'
    ref = lombscargle_harmonic_fit(t,y(:,1),freqs,nh,weights=w,floating_mean=.true.)
    if (maxval(abs(single%cosine(:,2:)-ref%cosine(:,2:))) > 1.0e-10_dp) &
        error stop 'reused plan without precentering'
    ! Strided first and second dimensions and zero channels.
    input = 0.0_dp
    input(1:2*n:2,1:2*nch:2) = y
    strided_fits = lombscargle_harmonic_fit_batch(plan,input(1:2*n:2,1:2*nch:2), &
        precenter=.true.,estimate_covariance=.true.)
    do j = 1, nch
        if (maxval(abs(strided_fits(j)%cosine(:,2:)-fits(j)%cosine(:,2:))) > 1.0e-11_dp) &
            error stop 'strided batch'
    end do
    strided_fits = lombscargle_harmonic_fit_batch(plan,y(:,:0))
    if (size(strided_fits) /= 0) error stop 'empty channels'
    ! A bad channel must not spoil results from other channels.
    input(1,3) = nan
    invalid_fits = lombscargle_harmonic_fit_batch(plan,input(1:2*n:2,1:2*nch:2))
    if (invalid_fits(2)%status /= irregular_status_invalid_input) error stop 'NaN channel status'
    if (.not. all(ieee_is_nan(invalid_fits(2)%power))) error stop 'NaN channel power'
    if (.not. invalid_fits(1)%valid(2)) error stop 'other channels still valid'
    single = lombscargle_harmonic_fit_prepared(plan,y(:n-1,1))
    if (single%status /= irregular_status_invalid_input) error stop 'sample mismatch'
    if (.not. all(ieee_is_nan(single%power))) error stop 'bad sample dimension NaNs'
    invalid = lombscargle_harmonic_prepare(t,freqs,nh,weights=-w)
    if (invalid%status /= irregular_status_invalid_input) error stop 'bad weights'
    single = lombscargle_harmonic_fit_prepared(invalid,y(:,1),estimate_covariance=.true.)
    if (single%status /= irregular_status_invalid_input) error stop 'bad plan fit'
    if (.not. all(ieee_is_nan(single%covariance))) error stop 'bad plan covariance'
    invalid = lombscargle_harmonic_prepare(t,freqs,huge(1))
    if (invalid%status /= irregular_status_invalid_input) error stop 'huge harmonics'
    if (invalid%nterms /= 0) error stop 'huge allocation prevention'
    other = lombscargle_harmonic_prepare(t,freqs,nh,weights=1.0e8_dp*w,floating_mean=.true.)
    if (any(other%valid .neqv. plan%valid)) error stop 'rescaled weights valid'
    single = lombscargle_harmonic_fit_prepared(other,y(:,1),precenter=.true.)
    if (maxval(abs(single%cosine(:,2:)-fits(1)%cosine(:,2:))) > 1.0e-10_dp) &
        error stop 'weight scale invariance'
    print '(a)', 'test_signal_harmonic_plan passed'
end program test_signal_harmonic_plan
