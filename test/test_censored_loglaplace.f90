! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_loglaplace
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        loglaplace_censored_loglikelihood, loglaplace_censored_nnlf, loglaplace_fit, &
        loglaplace_logpdf, loglaplace_logcdf, loglaplace_logsf, loglaplace_loglikelihood, &
        laplace_censored_loglikelihood, fit_status_invalid_input, fit_status_max_iter, &
        fit_status_no_finite_objective
    use scifort_math, only : quiet_nan, positive_infinity, log1p_safe
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, transformed
    type(fit_result) :: fit, parent
    real(dp) :: inf, nan, actual, expected, x, c, a, b, small
    integer :: i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
        left=[0.3_dp],right=[3.5_dp], &
        interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    actual = loglaplace_censored_loglikelihood(data,1.7_dp,scale=1.3_dp)
    call close(actual,-15.530991142547814_dp,5.0e-14_dp,'mixed mpmath reference')
    call close(loglaplace_censored_nnlf(data,1.7_dp,scale=1.3_dp),-actual,0.0_dp,'NNLF')
    do j = 1, 3
        c = 0.3_dp
        if (j == 2) c = 1.0_dp
        if (j == 3) c = 3.0_dp
        do i = 1, 12
            x = exp(real(i-6,dp))
            single = make_censored_data(uncensored=[x])
            call close(loglaplace_censored_loglikelihood(single,c),loglaplace_logpdf(x,c), &
                3.0e-14_dp,'scalar density and defaults')
            single = make_censored_data(left=[x])
            call close(loglaplace_censored_loglikelihood(single,c),loglaplace_logcdf(x,c), &
                3.0e-14_dp,'scalar CDF')
            single = make_censored_data(right=[x])
            call close(loglaplace_censored_loglikelihood(single,c),loglaplace_logsf(x,c), &
                3.0e-14_dp,'scalar SF')
        end do
        transformed = make_censored_data(uncensored=log(data%uncensored/1.3_dp), &
            left=log(data%left/1.3_dp),right=log(data%right/1.3_dp), &
            interval=log(data%interval/1.3_dp))
        expected = laplace_censored_loglikelihood(transformed,scale=1.0_dp/c) &
            -sum(log(data%uncensored))
        call close(loglaplace_censored_loglikelihood(data,c,scale=1.3_dp),expected, &
            3.0e-14_dp,'Laplace transformation with exact-density Jacobian')
    end do
    single = make_censored_data(uncensored=[1.0_dp,2.0_dp,4.0_dp])
    call close(loglaplace_censored_loglikelihood(single,1.7_dp,scale=1.3_dp), &
        loglaplace_loglikelihood(single%uncensored,1.7_dp,scale=1.3_dp), &
        2.0e-14_dp,'original array likelihood')
    call loglaplace_fit(single%uncensored,[1.7_dp,0.0_dp,1.3_dp], &
        [1.7_dp,0.0_dp,1.3_dp],parent)
    call loglaplace_fit(single,[1.7_dp,0.0_dp,1.3_dp],[1.7_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success .and. parent%success,'fixed fit convergence')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'array fit compatibility')
    ! At fixed location and scale, c = n / sum(abs(log(x/scale))).
    call loglaplace_fit(single,[0.1_dp,0.0_dp,2.0_dp],[10.0_dp,0.0_dp,2.0_dp],fit)
    call check(fit%success,'analytic shape MLE converges')
    call close(fit%params(1),3.0_dp/(2.0_dp*log(2.0_dp)),4.0e-6_dp,'analytic shape MLE')
    call loglaplace_fit(data,[0.05_dp,0.0_dp,0.05_dp],[10.0_dp,0.0_dp,10.0_dp],fit)
    call check(fit%success,'mixed fit converges')
    call close(fit%params(1),0.9300659117552635_dp,4.0e-4_dp,'SciPy fitted shape')
    call close(fit%params(3),1.240969675551833_dp,4.0e-4_dp,'SciPy fitted scale')
    call close(fit%nllf,loglaplace_censored_nnlf(data,fit%params(1), &
        fit%params(2),fit%params(3)),0.0_dp,'fit objective')
    call loglaplace_fit(data,[1.7_dp,-0.1_dp,1.3_dp],[1.7_dp,0.1_dp,1.3_dp],fit)
    call check(fit%success,'free location supported')
    call loglaplace_fit(data,[1.7_dp,0.0_dp,0.05_dp],[1.7_dp,0.0_dp,10.0_dp],fit)
    call check(fit%success,'fixed shape scale fit converges')
    call close(fit%params(3),1.2410156250000006_dp,4.0e-4_dp,'SciPy fixed-shape scale')
    call loglaplace_fit(data,[0.1_dp,-0.1_dp,0.1_dp],[3.0_dp,0.1_dp,3.0_dp],fit)
    call check(fit%success,'all parameters free')
    call loglaplace_fit(data,[1.0_dp,0.0_dp,0.0_dp],[1.0_dp,0.0_dp,0.0_dp],fit)
    call check(.not. fit%success,'invalid fixed scale fit')
    call loglaplace_fit(data,[0.05_dp,0.0_dp,0.05_dp],[10.0_dp,0.0_dp,10.0_dp],fit, &
        guess=[1.0_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls')
    call loglaplace_fit(empty,[1.0_dp,0.0_dp,0.1_dp],[2.0_dp,0.0_dp,2.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call close(loglaplace_censored_loglikelihood(empty,1.0_dp),0.0_dp,0.0_dp,'empty likelihood')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,0.0_dp)),'zero shape')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,-1.0_dp)),'negative shape')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,nan)),'NaN shape')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,inf)),'infinite shape')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,1.0_dp,loc=inf)),'bad loc')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,1.0_dp,scale=0.0_dp)), &
        'zero scale')
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(empty,1.0_dp,scale=nan)), &
        'NaN scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(loglaplace_censored_loglikelihood(single,1.0_dp)),'NaN data')
    single = make_censored_data(uncensored=[0.0_dp])
    call check(loglaplace_censored_loglikelihood(single,0.5_dp) == inf,'singular endpoint')
    call close(loglaplace_censored_loglikelihood(single,1.0_dp),-log(2.0_dp), &
        0.0_dp,'finite endpoint')
    call check(loglaplace_censored_loglikelihood(single,2.0_dp) == -inf,'zero endpoint')
    call loglaplace_fit(single,[0.5_dp,0.0_dp,1.0_dp],[0.5_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'singular fit rejected')
    single = make_censored_data(left=[0.0_dp])
    call check(loglaplace_censored_loglikelihood(single,1.0_dp) == -inf,'zero endpoint CDF')
    single = make_censored_data(right=[-1.0_dp])
    call close(loglaplace_censored_loglikelihood(single,1.0_dp),0.0_dp,0.0_dp,'unit SF')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(loglaplace_censored_loglikelihood(single,1.7_dp),loglaplace_logcdf(2.0_dp,1.7_dp), &
        2.0e-14_dp,'infinite lower interval')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(loglaplace_censored_loglikelihood(single,1.7_dp),loglaplace_logsf(2.0_dp,1.7_dp), &
        2.0e-14_dp,'infinite upper interval')
    single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
    call close(loglaplace_censored_loglikelihood(single,1.7_dp),log(1.7_dp/2.0_dp), &
        2.0e-14_dp,'equal interval endpoints')
    a = 1.0_dp
    b = nearest(a,1.0_dp)
    single = make_censored_data(interval=reshape([a,b],[1,2]))
    call close(loglaplace_censored_loglikelihood(single,1.7_dp), &
        log(1.7_dp/2.0_dp)+log(b-a),3.0e-14_dp,'adjacent interval at join')
    small = nearest(0.0_dp,1.0_dp)
    single = make_censored_data(interval=reshape([0.5_dp,2.0_dp],[1,2]))
    call close(loglaplace_censored_loglikelihood(single,small), &
        log(small)+log(log(2.0_dp)),2.0e-13_dp,'smallest shape interval')
    single = make_censored_data(interval=reshape([1.0e-100_dp,1.0e100_dp],[1,2]))
    actual = loglaplace_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-1.0e-100_dp),1.0_dp,3.0e-14_dp,'near-unit interval relative error')
    single = make_censored_data(right=[1.0e300_dp])
    call close(loglaplace_censored_loglikelihood(single,2.0_dp), &
        -log(2.0_dp)-2.0_dp*log(1.0e300_dp),3.0e-13_dp,'underflowed right probability')
    single = make_censored_data(left=[1.0e-300_dp])
    call close(loglaplace_censored_loglikelihood(single,2.0_dp), &
        -log(2.0_dp)+2.0_dp*log(1.0e-300_dp),3.0e-13_dp,'underflowed left probability')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(loglaplace_censored_loglikelihood(single,2.0_dp,-1.0e308_dp,1.0e308_dp), &
        -log(1.0e308_dp)-3.0_dp*log(2.0_dp),3.0e-13_dp,'overflowing raw displacement')
    single = make_censored_data(right=[2.0_dp])
    call check(loglaplace_censored_loglikelihood(single,huge(1.0_dp),scale=0.1_dp) == -inf, &
        'unrepresentable log tail')
    write(*,'(a,i0,a)') 'Log-Laplace censored: ',checks,' checks'
    if (failures /= 0) error stop 'Log-Laplace censoring failures'
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
end program test_censored_loglaplace
