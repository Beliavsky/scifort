! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_halfnorm
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        halfnorm_censored_loglikelihood, halfnorm_censored_nnlf, halfnorm_fit, &
        halfnorm_logpdf, halfnorm_logcdf, halfnorm_logsf, halfnorm_loglikelihood, &
        normal_logpdf, normal_logsf, normal_sf, &
        fit_status_invalid_input, fit_status_max_iter, fit_status_no_finite_objective
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty
    type(fit_result) :: fit, parent
    real(dp) :: x, inf, nan, actual, expected, a, b, c0, base
    integer :: i, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    c0 = log(2.0_dp/acos(-1.0_dp))/2.0_dp
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    actual = halfnorm_censored_loglikelihood(data,scale=1.3_dp)
    call close(actual,-19.203965730287618_dp,7.0e-14_dp,'mixed mpmath reference')
    call close(halfnorm_censored_nnlf(data,scale=1.3_dp),-actual,0.0_dp,'negative likelihood')
    do i = 1, 12
        x = 0.3_dp*real(i,dp)
        single = make_censored_data(uncensored=[x])
        call close(halfnorm_censored_loglikelihood(single),halfnorm_logpdf(x), &
            2.0e-14_dp,'scalar density and defaults')
        single = make_censored_data(left=[x])
        call close(halfnorm_censored_loglikelihood(single),halfnorm_logcdf(x), &
            2.0e-14_dp,'scalar CDF')
        single = make_censored_data(right=[x])
        call close(halfnorm_censored_loglikelihood(single),halfnorm_logsf(x), &
            2.0e-14_dp,'scalar SF')
    end do
    single = make_censored_data(uncensored=[1.0_dp,2.0_dp,3.0_dp])
    call close(halfnorm_censored_loglikelihood(single,scale=1.3_dp), &
        halfnorm_loglikelihood(single%uncensored,scale=1.3_dp),2.0e-14_dp,'array agreement')
    call halfnorm_fit(single,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%success,'exact scale MLE converges')
    call close(fit%params(2),sqrt(14.0_dp/3.0_dp),4.0e-6_dp,'analytic RMS scale MLE')
    call halfnorm_fit(single%uncensored,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],parent)
    call halfnorm_fit(single,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],fit)
    call check(fit%success .and. parent%success,'array and censored fixed fits')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'array fit unchanged')
    call halfnorm_fit(data,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],fit)
    call check(fit%success,'both parameters fixed')
    base = fit%nllf
    call halfnorm_fit(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%success,'mixed fixed-location fit')
    call close(fit%params(2),2.7028686618572966_dp,4.0e-4_dp,'SciPy mixed fitted scale')
    call halfnorm_fit(data,[-2.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'both parameters free')
    call close(fit%nllf,-halfnorm_censored_loglikelihood(data,fit%params(1),fit%params(2)), &
        0.0_dp,'mixed objective agreement')
    call halfnorm_fit(data,[-2.0_dp,1.3_dp],[0.29_dp,1.3_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'fixed scale and free location')
    call halfnorm_fit(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit, &
        guess=[0.0_dp,1.3_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls forwarded')
    call halfnorm_fit(empty,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call halfnorm_fit(data,[0.0_dp,1.0_dp],[0.0_dp,0.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'reversed bounds rejected')
    call close(halfnorm_censored_loglikelihood(empty),0.0_dp,0.0_dp,'valid empty likelihood')
    call check(ieee_is_nan(halfnorm_censored_loglikelihood(empty,loc=inf)),'infinite location')
    call check(ieee_is_nan(halfnorm_censored_loglikelihood(empty,scale=nan)),'NaN scale')
    call check(ieee_is_nan(halfnorm_censored_loglikelihood(empty,scale=inf)),'infinite scale')
    call check(ieee_is_nan(halfnorm_censored_loglikelihood(empty,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(halfnorm_censored_loglikelihood(empty,scale=-1.0_dp)),'negative scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(halfnorm_censored_loglikelihood(single)),'invalid observation')
    single = make_censored_data(uncensored=[0.0_dp])
    call close(halfnorm_censored_loglikelihood(single),c0,2.0e-16_dp,'finite endpoint density')
    call halfnorm_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%success,'endpoint included in fit')
    single = make_censored_data(uncensored=[-1.0_dp])
    call check(halfnorm_censored_loglikelihood(single) == -inf,'below-support density')
    call halfnorm_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'impossible exact fit')
    single = make_censored_data(left=[0.0_dp])
    call check(halfnorm_censored_loglikelihood(single) == -inf,'zero endpoint CDF')
    single = make_censored_data(right=[-1.0_dp])
    call close(halfnorm_censored_loglikelihood(single),0.0_dp,0.0_dp,'unit SF below support')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(halfnorm_censored_loglikelihood(single),halfnorm_logcdf(2.0_dp), &
        2.0e-14_dp,'infinite lower endpoint')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(halfnorm_censored_loglikelihood(single),halfnorm_logsf(2.0_dp), &
        2.0e-14_dp,'infinite upper endpoint')
    single = make_censored_data(interval=reshape([-2.0_dp,1.0_dp],[1,2]))
    call close(halfnorm_censored_loglikelihood(single),halfnorm_logcdf(1.0_dp), &
        2.0e-14_dp,'clipped interval')
    single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
    call close(halfnorm_censored_loglikelihood(single),halfnorm_logpdf(1.0_dp), &
        2.0e-14_dp,'equal endpoints mean exact')
    single = make_censored_data(right=[40.0_dp])
    call close(halfnorm_censored_loglikelihood(single),normal_logsf(40.0_dp)+log(2.0_dp), &
        3.0e-13_dp,'extreme SF')
    single = make_censored_data(left=[10.0_dp])
    expected = -2.0_dp*normal_sf(10.0_dp)
    call close(halfnorm_censored_loglikelihood(single)/expected,1.0_dp,5.0e-14_dp, &
        'near-one CDF retains excluded mass')
    single = make_censored_data(right=[1.0e-300_dp])
    expected = -sqrt(2.0_dp/acos(-1.0_dp))*1.0e-300_dp
    call close(halfnorm_censored_loglikelihood(single)/expected,1.0_dp,3.0e-13_dp, &
        'near-one SF retains tiny CDF')
    single = make_censored_data(interval=reshape([1.0e-20_dp,10.0_dp],[1,2]))
    expected = -sqrt(2.0_dp/acos(-1.0_dp))*1.0e-20_dp-2.0_dp*normal_sf(10.0_dp)
    call close(halfnorm_censored_loglikelihood(single)/expected,1.0_dp,3.0e-13_dp, &
        'almost-full-support interval retains both excluded tails')
    a = 40.0_dp
    b = nearest(a,1.0_dp)
    single = make_censored_data(interval=reshape([a,b],[1,2]))
    call close(halfnorm_censored_loglikelihood(single), &
        halfnorm_logpdf(a)+log(b-a),2.0e-12_dp,'adjacent far-tail interval')
    a = nearest(0.0_dp,1.0_dp)
    single = make_censored_data(interval=reshape([0.0_dp,a],[1,2]))
    call close(halfnorm_censored_loglikelihood(single),c0+log(a),3.0e-13_dp, &
        'subnormal endpoint interval')
    single = make_censored_data(interval=reshape([1.0e-300_dp,2.0e-300_dp],[1,2]))
    call close(halfnorm_censored_loglikelihood(single,scale=1.0e300_dp), &
        c0+log(1.0e-300_dp)-log(1.0e300_dp),5.0e-13_dp,'underflowed scaled interval')
    single = make_censored_data(left=[1.0e-300_dp])
    call close(halfnorm_censored_loglikelihood(single,scale=1.0e300_dp), &
        c0+log(1.0e-300_dp)-log(1.0e300_dp),5.0e-13_dp,'underflowed scaled CDF')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(halfnorm_censored_loglikelihood(single,-1.0e308_dp,1.0e308_dp), &
        c0-2.0_dp-log(1.0e308_dp),3.0e-13_dp,'overflowing raw displacement')
    single = make_censored_data(uncensored=[huge(1.0_dp)])
    call check(halfnorm_censored_loglikelihood(single,scale=tiny(1.0_dp)) == -inf, &
        'unrepresentable log density')
    write(*,'(a,i0,a)') 'Half-normal censored: ',checks,' checks'
    if (failures /= 0) error stop 'Half-normal censoring failures'
contains
    subroutine check(condition,label)
        implicit none
        logical, intent(in) :: condition !! true when expected behavior holds
        character(*), intent(in) :: label !! failure diagnostic
        checks = checks+1
        if (condition) return
        failures = failures+1
        write(*,'(a)') 'FAIL: '//label
    end subroutine check
    subroutine close(value,reference,tolerance,label)
        implicit none
        real(dp), intent(in) :: value !! computed finite result
        real(dp), intent(in) :: reference !! expected finite value
        real(dp), intent(in) :: tolerance !! absolute tolerance for this regime
        character(*), intent(in) :: label !! failure diagnostic
        call check(ieee_is_finite(value),label//' finite')
        if (.not. ieee_is_finite(value)) return
        call check(abs(value-reference) <= tolerance,label)
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_halfnorm
