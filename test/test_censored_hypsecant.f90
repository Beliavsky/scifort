! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_hypsecant
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        hypsecant_censored_loglikelihood, hypsecant_censored_nnlf, hypsecant_fit, &
        hypsecant_logpdf, hypsecant_logcdf, hypsecant_logsf, hypsecant_loglikelihood, &
        fit_status_invalid_input, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    type(censored_data) :: data, single, reflected, empty
    type(fit_result) :: fit, parent
    real(dp) :: x, inf, nan, actual, expected, a, b, base, logpi, logtwo
    integer :: i, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    logpi = log(acos(-1.0_dp))
    logtwo = log(2.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    actual = hypsecant_censored_loglikelihood(data,0.4_dp,1.3_dp)
    call close(actual,-16.31028768324201_dp,7.0e-14_dp,'mixed mpmath reference')
    call close(hypsecant_censored_nnlf(data,0.4_dp,1.3_dp),-actual,0.0_dp,'negative likelihood')
    reflected = make_censored_data(uncensored=-data%uncensored,left=-data%right, &
        right=-data%left,interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
    call close(hypsecant_censored_loglikelihood(reflected,-0.4_dp,1.3_dp),actual, &
        2.0e-14_dp,'mixed reflection')
    do i = -12, 12
        x = 0.5_dp*real(i,dp)
        single = make_censored_data(uncensored=[x])
        call close(hypsecant_censored_loglikelihood(single),hypsecant_logpdf(x), &
            2.0e-14_dp,'scalar density and defaults')
        single = make_censored_data(left=[x])
        call close(hypsecant_censored_loglikelihood(single),hypsecant_logcdf(x), &
            2.0e-14_dp,'scalar CDF')
        single = make_censored_data(right=[x])
        call close(hypsecant_censored_loglikelihood(single),hypsecant_logsf(x), &
            2.0e-14_dp,'scalar SF')
    end do
    single = make_censored_data(uncensored=[-1.0_dp,0.0_dp,1.0_dp])
    call close(hypsecant_censored_loglikelihood(single,scale=1.3_dp), &
        hypsecant_loglikelihood(single%uncensored,scale=1.3_dp),2.0e-14_dp,'array agreement')
    call hypsecant_fit(single,[-2.0_dp,1.0_dp],[2.0_dp,1.0_dp],fit)
    call check(fit%success,'fixed-scale location MLE converges')
    call close(fit%params(1),0.0_dp,3.0e-6_dp,'symmetric sample analytic location MLE')
    call hypsecant_fit(single%uncensored,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],parent)
    call hypsecant_fit(single,[0.0_dp,1.3_dp],[0.0_dp,1.3_dp],fit)
    call check(fit%success .and. parent%success,'array and censored fixed fits')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'array fit unchanged')
    call hypsecant_fit(data,[0.4_dp,1.3_dp],[0.4_dp,1.3_dp],fit)
    call check(fit%success,'both parameters fixed')
    base = fit%nllf
    call hypsecant_fit(data,[-5.0_dp,0.05_dp],[6.0_dp,10.0_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'both parameters free')
    call close(fit%params(1),1.7600566713787145_dp,4.0e-4_dp,'SciPy mixed fitted location')
    call close(fit%params(2),1.3906362832127406_dp,4.0e-4_dp,'SciPy mixed fitted scale')
    call close(fit%nllf,-hypsecant_censored_loglikelihood(data,fit%params(1),fit%params(2)), &
        0.0_dp,'mixed objective agreement')
    call hypsecant_fit(data,[0.4_dp,0.05_dp],[0.4_dp,10.0_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'fixed location and free scale')
    call hypsecant_fit(data,[-5.0_dp,0.05_dp],[6.0_dp,10.0_dp],fit, &
        guess=[0.4_dp,1.3_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls forwarded')
    call hypsecant_fit(empty,[0.0_dp,0.05_dp],[0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call hypsecant_fit(data,[0.0_dp,1.0_dp],[0.0_dp,0.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'reversed bounds rejected')
    call close(hypsecant_censored_loglikelihood(empty),0.0_dp,0.0_dp,'valid empty likelihood')
    call check(ieee_is_nan(hypsecant_censored_loglikelihood(empty,loc=inf)),'infinite location')
    call check(ieee_is_nan(hypsecant_censored_loglikelihood(empty,scale=nan)),'NaN scale')
    call check(ieee_is_nan(hypsecant_censored_loglikelihood(empty,scale=inf)),'infinite scale')
    call check(ieee_is_nan(hypsecant_censored_loglikelihood(empty,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(hypsecant_censored_loglikelihood(empty,scale=-1.0_dp)),'negative scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(hypsecant_censored_loglikelihood(single)),'invalid observation')
    single = make_censored_data(uncensored=[0.0_dp])
    call close(hypsecant_censored_loglikelihood(single),-logpi,3.0e-16_dp,'center density')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(hypsecant_censored_loglikelihood(single),hypsecant_logcdf(2.0_dp), &
        2.0e-14_dp,'infinite lower endpoint')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(hypsecant_censored_loglikelihood(single),hypsecant_logsf(2.0_dp), &
        2.0e-14_dp,'infinite upper endpoint')
    single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
    call close(hypsecant_censored_loglikelihood(single),hypsecant_logpdf(1.0_dp), &
        2.0e-14_dp,'equal endpoints mean exact')
    single = make_censored_data(right=[1000.0_dp])
    call close(hypsecant_censored_loglikelihood(single),logtwo-logpi-1000.0_dp, &
        2.0e-13_dp,'extreme SF')
    single = make_censored_data(left=[-1000.0_dp])
    call close(hypsecant_censored_loglikelihood(single),logtwo-logpi-1000.0_dp, &
        2.0e-13_dp,'extreme CDF')
    single = make_censored_data(left=[40.0_dp])
    expected = -2.0_dp*exp(-40.0_dp)/acos(-1.0_dp)
    call close(hypsecant_censored_loglikelihood(single)/expected,1.0_dp,3.0e-14_dp, &
        'near-one CDF retains excluded tail')
    single = make_censored_data(interval=reshape([-40.0_dp,41.0_dp],[1,2]))
    expected = -2.0_dp*(exp(-40.0_dp)+exp(-41.0_dp))/acos(-1.0_dp)
    call close(hypsecant_censored_loglikelihood(single)/expected,1.0_dp,3.0e-14_dp, &
        'almost-full-support interval retains both excluded tails')
    a = 1000.0_dp
    b = nearest(a,1.0_dp)
    single = make_censored_data(interval=reshape([a,b],[1,2]))
    expected = hypsecant_logpdf(a)+log(b-a)
    call close(hypsecant_censored_loglikelihood(single),expected,3.0e-12_dp, &
        'adjacent far-tail interval')
    single = make_censored_data(interval=reshape([-b,-a],[1,2]))
    call close(hypsecant_censored_loglikelihood(single),expected,3.0e-12_dp, &
        'negative adjacent far-tail interval')
    a = nearest(0.0_dp,1.0_dp)
    single = make_censored_data(interval=reshape([-a,a],[1,2]))
    call close(hypsecant_censored_loglikelihood(single),logtwo+log(a)-logpi, &
        3.0e-13_dp,'subnormal crossing interval')
    single = make_censored_data(interval=reshape([1.0e-300_dp,2.0e-300_dp],[1,2]))
    expected = log(1.0e-300_dp)-log(1.0e300_dp)-logpi
    call close(hypsecant_censored_loglikelihood(single,scale=1.0e300_dp),expected, &
        5.0e-13_dp,'underflowed scaled positive interval')
    single = make_censored_data(interval=reshape([-2.0e-300_dp,-1.0e-300_dp],[1,2]))
    call close(hypsecant_censored_loglikelihood(single,scale=1.0e300_dp),expected, &
        5.0e-13_dp,'underflowed scaled negative interval')
    single = make_censored_data(interval=reshape([-1.0e-300_dp,1.0e-300_dp],[1,2]))
    call close(hypsecant_censored_loglikelihood(single,scale=1.0e300_dp), &
        expected+logtwo,5.0e-13_dp,'underflowed scaled crossing interval')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(hypsecant_censored_loglikelihood(single,-1.0e308_dp,1.0e308_dp), &
        -logpi-log(cosh(2.0_dp))-log(1.0e308_dp),3.0e-13_dp,'overflowing raw displacement')
    single = make_censored_data(uncensored=[huge(1.0_dp)])
    call check(hypsecant_censored_loglikelihood(single,scale=tiny(1.0_dp)) == -inf, &
        'unrepresentable log density')
    single = make_censored_data(interval=reshape([1.0_dp,2.0_dp],[1,2]))
    call check(hypsecant_censored_loglikelihood(single,scale=nearest(0.0_dp,1.0_dp)) == -inf, &
        'unrepresentable interval log mass')
    write(*,'(a,i0,a)') 'Hypsecant censored: ',checks,' checks'
    if (failures /= 0) error stop 'Hypsecant censoring failures'
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
end program test_censored_hypsecant
