! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_rayleigh_max
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        rayleigh_censored_loglikelihood, rayleigh_censored_nnlf, rayleigh_fit, &
        weibull_max_censored_loglikelihood, weibull_max_censored_nnlf, weibull_max_fit, &
        rayleigh_logpdf, rayleigh_logcdf, rayleigh_logsf, rayleigh_loglikelihood, &
        weibull_max_logpdf, weibull_max_logcdf, weibull_max_logsf, weibull_max_loglikelihood, &
        weibull_censored_loglikelihood, exponential_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    implicit none
    type(censored_data) :: data, reflected, empty, single
    type(fit_result) :: fit
    real(dp) :: inf, nan, logtwo, tiny_width, base, expected, x, a, b, c
    integer :: checks, failures, i
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    logtwo = log(2.0_dp)
    tiny_width = nearest(0.0_dp,1.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    reflected = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp],right=[-0.3_dp], &
        left=[-3.5_dp],interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
    call close(rll(data,0.0_dp,1.3_dp),-17.949503441709364_dp,6.0e-14_dp,'Rayleigh mixed mpmath')
    call close(rayleigh_censored_nnlf(data,0.0_dp,1.3_dp),17.949503441709364_dp, &
        6.0e-14_dp,'Rayleigh negative likelihood')
    call rayleigh_fit(data,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%success,'Rayleigh fixed-location fit')
    call close(fit%params(2),1.8108567877921697_dp,4.0e-4_dp,'Rayleigh SciPy fitted scale')
    call close(fit%nllf,-rll(data,fit%params(1),fit%params(2)),0.0_dp,'Rayleigh fit objective')
    call rayleigh_fit(data,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],fit)
    call check(fit%success,'Rayleigh both fixed')
    base = fit%nllf
    call rayleigh_fit(data,[-3.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'Rayleigh both free')
    call rayleigh_fit(data,[-3.0_dp,1.3_dp],[0.29_dp,1.3_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'Rayleigh fixed scale')
    call rayleigh_fit(data,[-3.0_dp,0.05_dp],[0.29_dp,10.0_dp],fit, &
        guess=[0.0_dp,1.3_dp],max_iter=1)
    call check(fit%status == fit_status_max_iter,'Rayleigh iteration limit')
    call close(rll(data,0.0_dp,1.3_dp), &
        weibull_censored_loglikelihood(data,2.0_dp,0.0_dp,sqrt(2.0_dp)*1.3_dp), &
        6.0e-14_dp,'Rayleigh Weibull shape-two identity')
    call close(wll(reflected,2.0_dp,0.0_dp,sqrt(2.0_dp)*1.3_dp), &
        rll(data,0.0_dp,1.3_dp),6.0e-14_dp,'reflected Rayleigh identity')
    call close(wll(reflected,1.7_dp,0.0_dp,1.3_dp),-20.5129071280305_dp, &
        7.0e-14_dp,'Weibull maximum mixed mpmath')
    call close(weibull_max_censored_nnlf(reflected,1.7_dp,0.0_dp,1.3_dp), &
        20.5129071280305_dp,7.0e-14_dp,'Weibull maximum negative likelihood')
    do i = 1, 5
        c = 0.25_dp*real(i*i,dp)
        call close(wll(reflected,c,0.0_dp,1.3_dp), &
            weibull_censored_loglikelihood(data,c,0.0_dp,1.3_dp), &
            2.0e-12_dp,'mixed reflection for multiple shapes')
    end do
    call close(wll(reflected,1.0_dp,0.0_dp,1.3_dp), &
        exponential_censored_loglikelihood(data,0.0_dp,1.3_dp), &
        6.0e-14_dp,'shape-one reflected exponential')
    call weibull_max_fit(reflected,[0.1_dp,0.0_dp,0.05_dp],[6.0_dp,0.0_dp,10.0_dp],fit)
    call check(fit%success,'Weibull maximum fixed-location fit')
    call close(fit%params(1),1.0192387961160727_dp,4.0e-4_dp,'Weibull maximum SciPy shape')
    call close(fit%params(3),2.22130765498472_dp,4.0e-4_dp,'Weibull maximum SciPy scale')
    call close(fit%nllf,-wll(reflected,fit%params(1),fit%params(2),fit%params(3)), &
        0.0_dp,'Weibull maximum fit objective')
    call weibull_max_fit(reflected,[1.7_dp,0.0_dp,1.3_dp],[1.7_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success,'Weibull maximum all parameters fixed')
    base = fit%nllf
    call weibull_max_fit(reflected,[1.7_dp,-0.29_dp,0.05_dp],[1.7_dp,3.0_dp,10.0_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'Weibull maximum free location and scale')
    call weibull_max_fit(reflected,[0.1_dp,0.0_dp,0.05_dp],[6.0_dp,0.0_dp,10.0_dp],fit, &
        guess=[1.7_dp,0.0_dp,1.3_dp],max_iter=1)
    call check(fit%status == fit_status_max_iter,'Weibull maximum iteration limit')
    do i = 1, 10
        x = 0.5_dp*real(i,dp)
        single = make_censored_data(uncensored=[x])
        call close(rll(single,0.1_dp,1.3_dp),rayleigh_logpdf(x,0.1_dp,1.3_dp), &
            2.0e-14_dp,'Rayleigh scalar density')
        single = make_censored_data(left=[x])
        call close(rll(single,0.1_dp,1.3_dp),rayleigh_logcdf(x,0.1_dp,1.3_dp), &
            2.0e-14_dp,'Rayleigh scalar CDF')
        single = make_censored_data(right=[x])
        call close(rll(single,0.1_dp,1.3_dp),rayleigh_logsf(x,0.1_dp,1.3_dp), &
            2.0e-14_dp,'Rayleigh scalar SF')
        single = make_censored_data(uncensored=[-x])
        call close(wll(single,1.7_dp,0.1_dp,1.3_dp),weibull_max_logpdf(-x,1.7_dp,0.1_dp,1.3_dp), &
            3.0e-14_dp,'Weibull maximum scalar density')
        single = make_censored_data(left=[-x])
        call close(wll(single,1.7_dp,0.1_dp,1.3_dp),weibull_max_logcdf(-x,1.7_dp,0.1_dp,1.3_dp), &
            3.0e-14_dp,'Weibull maximum scalar CDF')
        single = make_censored_data(right=[-x])
        call close(wll(single,1.7_dp,0.1_dp,1.3_dp),weibull_max_logsf(-x,1.7_dp,0.1_dp,1.3_dp), &
            3.0e-14_dp,'Weibull maximum scalar SF')
    end do
    single = make_censored_data(uncensored=[1.0_dp,2.0_dp,3.0_dp],right=[4.0_dp])
    call rayleigh_fit(single,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call close(fit%params(2),sqrt(5.0_dp),3.0e-6_dp,'analytic exact/right Rayleigh scale')
    single = make_censored_data(uncensored=[0.5_dp,1.0_dp,3.0_dp])
    call close(rll(single,0.0_dp,1.3_dp),rayleigh_loglikelihood(single%uncensored,0.0_dp,1.3_dp), &
        2.0e-14_dp,'Rayleigh array likelihood')
    call rayleigh_fit(single%uncensored,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    base = fit%nllf
    call rayleigh_fit(single,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call close(fit%nllf,base,2.0e-12_dp,'Rayleigh array fit')
    single = make_censored_data(uncensored=[-0.5_dp,-1.0_dp,-3.0_dp])
    call close(wll(single,1.7_dp,0.0_dp,1.3_dp), &
        weibull_max_loglikelihood(single%uncensored,1.7_dp,0.0_dp,1.3_dp), &
        3.0e-14_dp,'Weibull maximum array likelihood')
    call weibull_max_fit(single%uncensored,[1.7_dp,0.0_dp,0.05_dp],[1.7_dp,0.0_dp,10.0_dp],fit)
    base = fit%nllf
    call weibull_max_fit(single,[1.7_dp,0.0_dp,0.05_dp],[1.7_dp,0.0_dp,10.0_dp],fit)
    call close(fit%nllf,base,2.0e-12_dp,'Weibull maximum array fit')

    single = make_censored_data(uncensored=[0.0_dp])
    call check(rll(single) == -inf,'Rayleigh zero endpoint density')
    call check(wll(single,0.5_dp) == inf,'Weibull maximum infinite endpoint density')
    call close(wll(single,1.0_dp,scale=1.3_dp),-log(1.3_dp),0.0_dp,'shape-one endpoint density')
    call check(wll(single,2.0_dp) == -inf,'Weibull maximum zero endpoint density')
    single = make_censored_data(uncensored=[0.0_dp,1.0_dp])
    call check(ieee_is_nan(wll(single,0.5_dp)),'mixed infinite endpoint and impossible density')
    single = make_censored_data(left=[0.0_dp])
    call check(rll(single) == -inf,'Rayleigh zero endpoint CDF')
    call close(wll(single,2.0_dp),0.0_dp,0.0_dp,'Weibull maximum unit endpoint CDF')
    single = make_censored_data(right=[0.0_dp])
    call close(rll(single),0.0_dp,0.0_dp,'Rayleigh unit endpoint SF')
    call check(wll(single,2.0_dp) == -inf,'Weibull maximum zero endpoint SF')
    call interval(-1.0_dp,1.0_dp)
    base = rll(single)
    call interval(0.0_dp,1.0_dp)
    call close(rll(single),base,0.0_dp,'Rayleigh interval clipping')
    call interval(-1.0_dp,0.0_dp)
    call check(rll(single) == -inf,'Rayleigh impossible interval')
    base = wll(single,2.0_dp)
    call interval(-1.0_dp,1.0_dp)
    call close(wll(single,2.0_dp),base,0.0_dp,'Weibull maximum interval clipping')
    call interval(0.0_dp,tiny_width)
    call close(rll(single),2.0_dp*log(tiny_width)-logtwo,4.0e-13_dp,'smallest Rayleigh interval')
    call interval(-tiny_width,0.0_dp)
    call close(wll(single,0.3_dp),0.3_dp*log(tiny_width),2.0e-13_dp,'smallest Weibull interval')
    call interval(1.0e-300_dp,2.0e-300_dp)
    expected = 2.0_dp*log(1.0e-300_dp)-2.0_dp*log(1.0e300_dp)+log(1.5_dp)
    call close(rll(single,scale=1.0e300_dp),expected,8.0e-13_dp,'Rayleigh underflowed standardized width')
    call interval(-2.0e-300_dp,-1.0e-300_dp)
    expected = 2.0_dp*(log(1.0e-300_dp)-log(1.0e300_dp))+log(3.0_dp)
    call close(wll(single,2.0_dp,scale=1.0e300_dp),expected,8.0e-13_dp,'Weibull underflowed width')
    a = 40.0_dp
    b = nearest(a,1.0_dp)
    call interval(a,b)
    call close(rll(single),log(a)-0.5_dp*a*a+log(b-a),1.0e-12_dp,'Rayleigh adjacent tail')
    call interval(-b,-a)
    call close(wll(single,2.0_dp),logtwo+log(a)-a*a+log(b-a), &
        2.0e-12_dp,'Weibull maximum adjacent tail')
    b = sqrt(1400.0_dp)
    single = make_censored_data(left=[b])
    expected = -exp(-0.5_dp*b*b)
    call close(rll(single)/expected,1.0_dp,3.0e-12_dp,'Rayleigh tiny near-one CDF log')
    single = make_censored_data(right=[1.0e-150_dp])
    call close(rll(single)/(-0.5e-300_dp),1.0_dp,5.0e-13_dp,'Rayleigh tiny near-one SF log')
    call interval(1.0e-150_dp,b)
    expected = -0.5e-300_dp-exp(-0.5_dp*b*b)
    call close(rll(single)/expected,1.0_dp,3.0e-12_dp,'Rayleigh almost whole support')
    call interval(-b,-1.0e-150_dp)
    expected = -1.0e-300_dp-exp(-b*b)
    call close(wll(single,2.0_dp)/expected,1.0_dp,3.0e-12_dp,'Weibull almost whole support')
    call interval(-1.0e308_dp,1.0e308_dp)
    call close(rll(single,scale=1.5e308_dp),log(1.0_dp-exp(-2.0_dp/9.0_dp)), &
        5.0e-14_dp,'Rayleigh scale conversion would overflow')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(rll(single,-1.0e308_dp,1.0e308_dp),logtwo-log(1.0e308_dp)-2.0_dp, &
        3.0e-13_dp,'Rayleigh raw displacement overflow')
    call check(rll(single,scale=1.0e-300_dp) == -inf,'Rayleigh unrepresentable log density')
    single = make_censored_data(uncensored=[-1.0e308_dp])
    call close(wll(single,2.0_dp,1.0e308_dp,1.0e308_dp),2.0_dp*logtwo-log(1.0e308_dp)-4.0_dp, &
        5.0e-13_dp,'Weibull raw displacement overflow')
    call check(wll(single,2.0_dp,scale=1.0e-300_dp) == -inf,'Weibull unrepresentable log density')
    single = make_censored_data(uncensored=[-1.0_dp])
    call close(wll(single,huge(1.0_dp)),log(huge(1.0_dp))-1.0_dp, &
        2.0e-13_dp,'largest shape at unit displacement')
    call close(wll(single,tiny_width),log(tiny_width)-1.0_dp, &
        2.0e-13_dp,'smallest shape at unit displacement')
    call interval(-inf,inf)
    call check(ieee_is_nan(rll(single)),'Rayleigh invalid infinite interval')
    call check(ieee_is_nan(wll(single,2.0_dp)),'Weibull invalid infinite interval')
    call interval(-inf,0.0_dp)
    call close(wll(single,2.0_dp),0.0_dp,0.0_dp,'Weibull full support')
    call interval(0.0_dp,inf)
    call close(rll(single),0.0_dp,0.0_dp,'Rayleigh full support')
    call close(rll(empty),0.0_dp,0.0_dp,'Rayleigh empty sample')
    call close(wll(empty,2.0_dp),0.0_dp,0.0_dp,'Weibull empty sample')
    call check(ieee_is_nan(rll(empty,scale=0.0_dp)),'Rayleigh empty invalid parameter')
    call check(ieee_is_nan(wll(empty,0.0_dp)),'Weibull empty invalid shape')
    call check(ieee_is_nan(rll(data,loc=inf)),'Rayleigh infinite location')
    call check(ieee_is_nan(rll(data,scale=nan)),'Rayleigh NaN scale')
    call check(ieee_is_nan(wll(reflected,inf)),'Weibull infinite shape')
    call check(ieee_is_nan(wll(reflected,nan)),'Weibull NaN shape')
    call check(ieee_is_nan(wll(reflected,2.0_dp,scale=-1.0_dp)),'Weibull negative scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(rll(single)),'Rayleigh NaN observation')
    call check(ieee_is_nan(wll(single,2.0_dp)),'Weibull NaN observation')
    call rayleigh_fit(empty,[0.0_dp,0.1_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'Rayleigh empty fit')
    call weibull_max_fit(empty,[0.1_dp,0.0_dp,0.1_dp],[3.0_dp,0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'Weibull empty fit')
    single = make_censored_data(uncensored=[0.0_dp])
    call rayleigh_fit(single,[0.0_dp,0.1_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'Rayleigh impossible endpoint fit')
    call weibull_max_fit(single,[0.5_dp,0.0_dp,1.0_dp],[0.5_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'Weibull infinite endpoint objective')
    write(*,'(a,i0,a,i0)') 'Rayleigh/Weibull max checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    pure function rll(sample,loc,scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = rayleigh_censored_loglikelihood(sample,loc,scale)
    end function rll
    pure function wll(sample,shape,loc,scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        real(dp), intent(in) :: shape !! finite positive Weibull shape
        real(dp), intent(in), optional :: loc !! finite upper support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = weibull_max_censored_loglikelihood(sample,shape,loc,scale)
    end function wll
    subroutine interval(lower,upper)
        implicit none
        real(dp), intent(in) :: lower !! lower interval endpoint
        real(dp), intent(in) :: upper !! upper interval endpoint
        single = make_censored_data(interval=reshape([lower,upper],[1,2]))
    end subroutine interval
    subroutine check(condition,label)
        implicit none
        logical, intent(in) :: condition !! true when expected behavior holds
        character(*), intent(in) :: label !! test description for failure diagnostics
        checks = checks+1
        if (condition) return
        failures = failures+1
        write(*,'(a)') 'FAIL: '//label
    end subroutine check
    subroutine close(value,reference,tolerance,label)
        implicit none
        real(dp), intent(in) :: value !! computed finite value
        real(dp), intent(in) :: reference !! expected finite result
        real(dp), intent(in) :: tolerance !! absolute error allowed in this regime
        character(*), intent(in) :: label !! test description for failure diagnostics
        call check(ieee_is_finite(value),label//' is finite')
        if (.not. ieee_is_finite(value)) return
        call check(abs(value-reference) <= tolerance,label)
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_rayleigh_max
