! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_pareto
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        pareto_censored_loglikelihood, pareto_censored_nnlf, pareto_fit, &
        pareto_logpdf, pareto_logcdf, pareto_logsf, pareto_loglikelihood, &
        exponential_censored_loglikelihood, loglaplace_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_max_iter, fit_status_no_finite_objective
    use scifort_math, only : quiet_nan, positive_infinity, log1p_safe
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, transformed
    type(fit_result) :: fit, parent
    real(dp) :: inf, nan, actual, expected, x, shape, a, b, small
    integer :: i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    small = nearest(0.0_dp,1.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
        left=[0.3_dp],right=[3.5_dp], &
        interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    actual = pareto_censored_loglikelihood(data,1.7_dp,scale=0.2_dp)
    call close(actual,-24.25414363441046_dp,6.0e-14_dp,'mixed mpmath reference')
    call close(pareto_censored_nnlf(data,1.7_dp,scale=0.2_dp),-actual,0.0_dp,'NNLF')
    do j = 1, 3
        shape = 0.3_dp
        if (j == 2) shape = 1.0_dp
        if (j == 3) shape = 3.0_dp
        do i = 1, 12
            x = exp(real(i,dp)/2.0_dp)
            single = make_censored_data(uncensored=[x])
            call close(pareto_censored_loglikelihood(single,shape),pareto_logpdf(x,shape), &
                3.0e-14_dp,'scalar density and defaults')
            ! Above z=1, Pareto density and survival are twice log-Laplace's.
            call close(pareto_censored_loglikelihood(single,shape), &
                loglaplace_censored_loglikelihood(single,shape)+log(2.0_dp), &
                3.0e-14_dp,'log-Laplace density identity')
            single = make_censored_data(left=[x])
            call close(pareto_censored_loglikelihood(single,shape),pareto_logcdf(x,shape), &
                3.0e-14_dp,'scalar CDF')
            single = make_censored_data(right=[x])
            call close(pareto_censored_loglikelihood(single,shape),pareto_logsf(x,shape), &
                3.0e-14_dp,'scalar SF')
        end do
        transformed = make_censored_data(uncensored=log(data%uncensored/0.2_dp), &
            left=log(data%left/0.2_dp),right=log(data%right/0.2_dp), &
            interval=log(data%interval/0.2_dp))
        expected = exponential_censored_loglikelihood(transformed,scale=1.0_dp/shape) &
            -sum(log(data%uncensored))
        call close(pareto_censored_loglikelihood(data,shape,scale=0.2_dp),expected, &
            4.0e-14_dp,'exponential transform with exact-density Jacobian')
    end do
    single = make_censored_data(uncensored=[1.0_dp,2.0_dp,4.0_dp])
    call close(pareto_censored_loglikelihood(single,1.7_dp), &
        pareto_loglikelihood(single%uncensored,1.7_dp),2.0e-14_dp,'array agreement')
    call pareto_fit(single%uncensored,[1.7_dp,0.0_dp,1.0_dp],[1.7_dp,0.0_dp,1.0_dp],parent)
    call pareto_fit(single,[1.7_dp,0.0_dp,1.0_dp],[1.7_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success .and. parent%success,'fixed fits converge')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'array fit unchanged')
    call pareto_fit(single,[0.1_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success,'analytic shape MLE converges')
    call close(fit%params(1),1.0_dp/log(2.0_dp),4.0e-6_dp,'analytic shape MLE')
    call pareto_fit(single,[2.0_dp,0.0_dp,0.1_dp],[2.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success,'scale MLE converges at exact endpoint')
    call close(fit%params(3),1.0_dp,3.0e-6_dp,'scale MLE is sample minimum')
    call pareto_fit(data,[0.05_dp,0.0_dp,0.2_dp],[10.0_dp,0.0_dp,0.2_dp],fit)
    call check(fit%success,'mixed shape fit converges')
    call close(fit%params(1),0.45595703124999953_dp,2.0e-4_dp,'SciPy fixed-support shape')
    call close(fit%nllf,pareto_censored_nnlf(data,fit%params(1), &
        fit%params(2),fit%params(3)),0.0_dp,'fit objective')
    call pareto_fit(data,[0.1_dp,-0.1_dp,0.05_dp],[3.0_dp,0.1_dp,0.25_dp],fit)
    call check(fit%success,'all parameters free')
    call pareto_fit(data,[0.1_dp,0.0_dp,0.2_dp],[3.0_dp,0.0_dp,0.2_dp],fit, &
        guess=[1.0_dp,0.0_dp,0.2_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls forwarded')
    call pareto_fit(empty,[1.0_dp,0.0_dp,0.1_dp],[2.0_dp,0.0_dp,2.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call close(pareto_censored_loglikelihood(empty,1.0_dp),0.0_dp,0.0_dp,'valid empty sample')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,0.0_dp)),'zero shape')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,-1.0_dp)),'negative shape')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,nan)),'NaN shape')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,inf)),'infinite shape')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,1.0_dp,loc=inf)),'infinite loc')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,1.0_dp,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(pareto_censored_loglikelihood(empty,1.0_dp,scale=nan)),'NaN scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(pareto_censored_loglikelihood(single,1.0_dp)),'invalid data')
    single = make_censored_data(uncensored=[0.5_dp])
    call check(pareto_censored_loglikelihood(single,1.0_dp) == -inf,'below support density')
    call pareto_fit(single,[1.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'impossible fit objective')
    single = make_censored_data(uncensored=[1.0_dp])
    call close(pareto_censored_loglikelihood(single,2.0_dp),log(2.0_dp),0.0_dp,'endpoint density')
    single = make_censored_data(left=[1.0_dp])
    call check(pareto_censored_loglikelihood(single,1.0_dp) == -inf,'endpoint CDF')
    single = make_censored_data(right=[0.5_dp])
    call close(pareto_censored_loglikelihood(single,1.0_dp),0.0_dp,0.0_dp,'unit SF')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(pareto_censored_loglikelihood(single,1.7_dp),pareto_logcdf(2.0_dp,1.7_dp), &
        2.0e-14_dp,'infinite lower interval endpoint')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(pareto_censored_loglikelihood(single,1.7_dp),pareto_logsf(2.0_dp,1.7_dp), &
        2.0e-14_dp,'infinite upper interval endpoint')
    single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
    call close(pareto_censored_loglikelihood(single,1.7_dp),log(1.7_dp),0.0_dp,'equal endpoints')
    single = make_censored_data(interval=reshape([-2.0_dp,2.0_dp],[1,2]))
    call close(pareto_censored_loglikelihood(single,1.7_dp),pareto_logcdf(2.0_dp,1.7_dp), &
        2.0e-14_dp,'support-clipped interval')
    single = make_censored_data(interval=reshape([-2.0_dp,1.0_dp],[1,2]))
    call check(pareto_censored_loglikelihood(single,1.7_dp) == -inf,'impossible interval')
    a = 2.0_dp
    b = nearest(a,1.0_dp)
    single = make_censored_data(interval=reshape([a,b],[1,2]))
    call close(pareto_censored_loglikelihood(single,1.7_dp), &
        pareto_logpdf(a,1.7_dp)+log(b-a),3.0e-14_dp,'adjacent tail interval')
    single = make_censored_data(left=[nearest(1.0_dp,1.0_dp)])
    call close(pareto_censored_loglikelihood(single,small), &
        log(small)+log(log1p_safe(epsilon(1.0_dp))),2.0e-13_dp,'underflowed cumulative hazard')
    single = make_censored_data(interval=reshape([1.0_dp,1.0e100_dp],[1,2]))
    actual = pareto_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-1.0e-100_dp),1.0_dp,3.0e-14_dp,'near-unit interval relative error')
    single = make_censored_data(right=[1.0e300_dp])
    call close(pareto_censored_loglikelihood(single,2.0_dp,scale=1.0e-300_dp), &
        -2.0_dp*(log(1.0e300_dp)-log(1.0e-300_dp)),4.0e-13_dp,'overflowing ratio tail')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(pareto_censored_loglikelihood(single,2.0_dp,-1.0e308_dp,1.0e308_dp), &
        -log(1.0e308_dp)-2.0_dp*log(2.0_dp),3.0e-13_dp,'overflowing raw displacement')
    single = make_censored_data(right=[2.0_dp])
    call check(pareto_censored_loglikelihood(single,huge(1.0_dp),scale=0.1_dp) == -inf, &
        'unrepresentable log tail')
    write(*,'(a,i0,a)') 'Pareto censored: ',checks,' checks'
    if (failures /= 0) error stop 'Pareto censoring failures'
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
end program test_censored_pareto
