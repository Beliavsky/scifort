! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_gibrat
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        gibrat_censored_loglikelihood, gibrat_censored_nnlf, gibrat_fit, lognormal_fit, &
        lognormal_censored_loglikelihood, normal_logsf, normal_logcdf, &
        gibrat_logpdf, gibrat_logcdf, gibrat_logsf, gibrat_loglikelihood, &
        fit_status_invalid_input, fit_status_max_iter, fit_status_no_finite_objective
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty
    type(fit_result) :: fit, parent
    real(dp) :: x, inf, nan, actual, expected, a, b
    integer :: i, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    actual = gibrat_censored_loglikelihood(data,scale=1.3_dp)
    call close(actual,-14.6038170449131_dp,5.0e-14_dp,'mixed mpmath reference')
    call close(gibrat_censored_nnlf(data,scale=1.3_dp),-actual,0.0_dp,'negative likelihood')
    do i = 1, 12
        x = exp(real(i-6,dp))
        single = make_censored_data(uncensored=[x])
        call close(gibrat_censored_loglikelihood(single),gibrat_logpdf(x), &
            2.0e-14_dp,'scalar density and defaults')
        single = make_censored_data(left=[x])
        call close(gibrat_censored_loglikelihood(single),gibrat_logcdf(x), &
            2.0e-14_dp,'scalar left tail')
        single = make_censored_data(right=[x])
        call close(gibrat_censored_loglikelihood(single),gibrat_logsf(x), &
            2.0e-14_dp,'scalar right tail')
        call close(gibrat_censored_loglikelihood(data,loc=-x,scale=x), &
            lognormal_censored_loglikelihood(data,1.0_dp,-x,x),0.0_dp,'mixed parent identity')
    end do
    single = make_censored_data(uncensored=[1.0_dp,2.0_dp,4.0_dp])
    call close(gibrat_censored_loglikelihood(single,scale=1.3_dp), &
        gibrat_loglikelihood(single%uncensored,scale=1.3_dp),2.0e-14_dp,'array agreement')
    call gibrat_fit(single,[0.0_dp,0.1_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%success,'exact-data scale fit converges')
    call close(fit%params(2),2.0_dp,3.0e-6_dp,'analytic geometric-mean MLE')
    call gibrat_fit(single%uncensored,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],parent)
    call gibrat_fit(single,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],fit)
    call check(fit%success .and. parent%success,'fixed fits converge')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'original array fit unchanged')
    call gibrat_fit(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%success,'mixed fit converges')
    call close(fit%params(2),1.3704101562500006_dp,2.0e-4_dp,'SciPy fitted scale')
    call close(fit%nllf,-gibrat_censored_loglikelihood(data,fit%params(1),fit%params(2)), &
        0.0_dp,'mixed objective agreement')
    call lognormal_fit(data,[1.0_dp,0.0_dp,0.05_dp],[1.0_dp,0.0_dp,10.0_dp],parent)
    call close(fit%params(2),parent%params(3),0.0_dp,'parent fit scale')
    call close(fit%nllf,parent%nllf,0.0_dp,'parent fit objective')
    call gibrat_fit(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit, &
        guess=[0.0_dp,1.3_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls forwarded')
    call gibrat_fit(data,[-0.1_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'free location rejected')
    call check(size(fit%params) == 2 .and. all(ieee_is_nan(fit%params)), &
        'invalid fit preserves two-parameter result')
    call gibrat_fit(data,[nan,0.05_dp],[nan,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'nonfinite location bounds rejected')
    call gibrat_fit(data,[0.0_dp,0.0_dp],[0.0_dp,0.0_dp],fit)
    call check(.not. fit%success,'zero fixed scale rejected')
    call gibrat_fit(empty,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call close(gibrat_censored_loglikelihood(empty),0.0_dp,0.0_dp,'valid empty likelihood')
    call check(ieee_is_nan(gibrat_censored_loglikelihood(empty,loc=inf)),'infinite location')
    call check(ieee_is_nan(gibrat_censored_loglikelihood(empty,scale=nan)),'NaN scale')
    call check(ieee_is_nan(gibrat_censored_loglikelihood(empty,scale=inf)),'infinite scale')
    call check(ieee_is_nan(gibrat_censored_loglikelihood(empty,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(gibrat_censored_loglikelihood(empty,scale=-1.0_dp)),'negative scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(gibrat_censored_loglikelihood(single)),'invalid observation')
    single = make_censored_data(uncensored=[0.0_dp])
    call check(gibrat_censored_loglikelihood(single) == -inf,'zero endpoint density')
    call gibrat_fit(single,[0.0_dp,1.0_dp],[0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'impossible exact fit')
    single = make_censored_data(left=[0.0_dp])
    call check(gibrat_censored_loglikelihood(single) == -inf,'zero endpoint CDF')
    single = make_censored_data(right=[-1.0_dp])
    call close(gibrat_censored_loglikelihood(single),0.0_dp,0.0_dp,'unit SF below support')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(gibrat_censored_loglikelihood(single),gibrat_logcdf(2.0_dp), &
        2.0e-14_dp,'infinite lower endpoint')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(gibrat_censored_loglikelihood(single),gibrat_logsf(2.0_dp), &
        2.0e-14_dp,'infinite upper endpoint')
    single = make_censored_data(interval=reshape([-2.0_dp,1.0_dp],[1,2]))
    call close(gibrat_censored_loglikelihood(single),-log(2.0_dp),2.0e-14_dp,'clipped interval')
    single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
    call close(gibrat_censored_loglikelihood(single),gibrat_logpdf(1.0_dp), &
        2.0e-14_dp,'equal endpoints mean exact observation')
    single = make_censored_data(right=[exp(40.0_dp)])
    call close(gibrat_censored_loglikelihood(single),normal_logsf(40.0_dp), &
        3.0e-13_dp,'extreme right log tail')
    single = make_censored_data(left=[exp(-40.0_dp)])
    call close(gibrat_censored_loglikelihood(single),normal_logcdf(-40.0_dp), &
        3.0e-13_dp,'extreme left log tail')
    a = 1.0_dp
    b = nearest(a,1.0_dp)
    single = make_censored_data(interval=reshape([a,b],[1,2]))
    call close(gibrat_censored_loglikelihood(single), &
        gibrat_logpdf(a)+log(b-a),3.0e-14_dp,'adjacent interval')
    single = make_censored_data(interval=reshape([1.0e-300_dp,2.0e-300_dp],[1,2]))
    call close(gibrat_censored_loglikelihood(single,scale=1.0e300_dp), &
        lognormal_censored_loglikelihood(single,1.0_dp,scale=1.0e300_dp), &
        0.0_dp,'underflowed standardized interval')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(gibrat_censored_loglikelihood(single,-1.0e308_dp,1.0e308_dp), &
        -log(2.0_dp)-log(1.0e308_dp)-log(2.0_dp)**2/2.0_dp &
        -log(2.0_dp*acos(-1.0_dp))/2.0_dp,3.0e-13_dp,'overflowing raw displacement')
    write(*,'(a,i0,a)') 'Gibrat censored: ',checks,' checks'
    if (failures /= 0) error stop 'Gibrat censoring failures'
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
end program test_censored_gibrat
