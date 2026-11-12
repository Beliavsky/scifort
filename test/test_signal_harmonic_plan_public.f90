! SPDX-License-Identifier: MIT
! Test the prepared multiharmonic API through the public SciFort module.
program test_signal_harmonic_plan_public
    use scifort_signal, only : dp, harmonic_fit_result, harmonic_fit_plan, &
        lombscargle_harmonic_prepare, lombscargle_harmonic_fit_prepared, &
        lombscargle_harmonic_fit_batch, lombscargle_harmonic_predict_std, &
        irregular_status_success
    implicit none
    integer, parameter :: n = 67
    real(dp) :: t(n), w(n), y(n,3), omega
    real(dp), allocatable :: se(:)
    type(harmonic_fit_plan) :: plan
    type(harmonic_fit_result) :: single
    type(harmonic_fit_result), allocatable :: fitted(:)
    integer :: i, j, status

    omega = 0.83_dp
    do i = 1, n
        t(i) = 0.13_dp*real(i,dp)+0.014_dp*cos(real(i*i,dp))
        w(i) = 0.3_dp+real(mod(i,7),dp)
    end do
    do j = 1, 3
        y(:,j) = real(j,dp)+cos(omega*t)*real(j,dp) &
            -0.3_dp*sin(omega*t)+0.15_dp*cos(2.0_dp*omega*t) &
            +0.02_dp*sin(real(j+2,dp)*t)
    end do
    plan = lombscargle_harmonic_prepare(t,[omega,1.71_dp],2, &
        weights=w,floating_mean=.true.)
    if (plan%status /= irregular_status_success) error stop 'public plan status'
    fitted = lombscargle_harmonic_fit_batch(plan,y,estimate_covariance=.true.)
    if (size(fitted) /= 3) error stop 'public batch count'
    do j = 1, 3
        if (fitted(j)%status /= irregular_status_success) error stop 'public batch fit'
        single = lombscargle_harmonic_fit_prepared(plan,y(:,j),estimate_covariance=.true.)
        if (maxval(abs(single%cosine-fitted(j)%cosine)) > 1.0e-12_dp) &
            error stop 'public single/batch agreement'
        se = lombscargle_harmonic_predict_std(fitted(j),[-2.0_dp,0.0_dp,4.0_dp],1,status)
        if (status /= irregular_status_success .or. any(se <= 0.0_dp)) &
            error stop 'public uncertainty'
    end do
    print '(a)', 'test_signal_harmonic_plan_public passed'
end program test_signal_harmonic_plan_public
