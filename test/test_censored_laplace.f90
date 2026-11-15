! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_laplace
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        laplace_censored_loglikelihood, laplace_censored_nnlf, laplace_fit, &
        laplace_loglikelihood, laplace_logpdf, laplace_logcdf, laplace_logsf, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    implicit none
    type(censored_data) :: data, other, empty
    type(fit_result) :: fit
    real(dp) :: inf, nan, base, expected, x, a, b, tiny_width, logtwo
    integer :: checks, failures, i

    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    logtwo = log(2.0_dp)
    tiny_width = nearest(0.0_dp,1.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], left=[0.3_dp], &
        right=[3.5_dp], interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    call close(ll(data,0.4_dp,1.3_dp),-16.62410935200545_dp,4.0e-14_dp,'mixed mpmath reference')
    call close(laplace_censored_nnlf(data,0.4_dp,1.3_dp),16.62410935200545_dp, &
        4.0e-14_dp,'negative log-likelihood')
    call laplace_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call check(fit%success,'mixed fit converges')
    call close(fit%params(1),1.2500073446335311_dp,4.0e-4_dp,'SciPy fitted center')
    call close(fit%params(2),1.6532044386643419_dp,4.0e-4_dp,'SciPy fitted scale')
    call close(fit%nllf,-ll(data,fit%params(1),fit%params(2)),0.0_dp,'fit objective')
    call laplace_fit(data,[0.4_dp,1.3_dp],[0.4_dp,1.3_dp],fit)
    call check(fit%success,'both parameters fixed')
    call close(fit%nllf,16.62410935200545_dp,4.0e-14_dp,'fixed objective')
    base = fit%nllf
    call laplace_fit(data,[0.4_dp,0.05_dp],[0.4_dp,10.0_dp],fit)
    call check(fit%success .and. fit%params(1) == 0.4_dp,'fixed center')
    call check(fit%nllf <= base,'free scale improves objective')
    call laplace_fit(data,[-5.0_dp,1.3_dp],[8.0_dp,1.3_dp],fit)
    call check(fit%success .and. fit%params(2) == 1.3_dp,'fixed scale')
    call check(fit%nllf <= base,'free center improves objective')
    other = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp], right=[-0.3_dp], &
        left=[-3.5_dp],interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
    call close(ll(data,0.4_dp,1.3_dp),ll(other,-0.4_dp,1.3_dp),2.0e-14_dp,'reflection')
    call laplace_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit, &
        guess=[0.4_dp,1.3_dp], max_iter=1)
    call check(fit%status == fit_status_max_iter,'iteration limit forwarded')

    data = make_censored_data(uncensored=[-3.0_dp,-1.0_dp,0.0_dp,2.0_dp,7.0_dp])
    call close(ll(data,0.4_dp,1.3_dp),laplace_loglikelihood(data%uncensored,0.4_dp,1.3_dp), &
        2.0e-14_dp,'array likelihood unchanged')
    call laplace_fit(data%uncensored,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    base = fit%nllf
    call laplace_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call check(fit%success,'exact fit converges')
    call close(fit%nllf,base,2.0e-13_dp,'array fit agreement')
    call close(fit%params(1),0.0_dp,3.0e-6_dp,'sample median MLE')
    call close(fit%params(2),2.6_dp,3.0e-6_dp,'sample mean absolute deviation MLE')

    do i = -10, 10
        x = real(i,dp)
        data = make_censored_data(uncensored=[x])
        call close(ll(data,0.4_dp,1.3_dp),laplace_logpdf(x,0.4_dp,1.3_dp), &
            4.0e-15_dp,'scalar density agreement')
        data = make_censored_data(left=[x])
        call close(ll(data,0.4_dp,1.3_dp),laplace_logcdf(x,0.4_dp,1.3_dp), &
            4.0e-15_dp,'scalar CDF agreement')
        data = make_censored_data(right=[x])
        call close(ll(data,0.4_dp,1.3_dp),laplace_logsf(x,0.4_dp,1.3_dp), &
            4.0e-15_dp,'scalar SF agreement')
    end do
    do i = -3, 3
        a = real(i,dp)
        b = a + 0.75_dp
        call interval(a,b)
        base = ll(data)
        if (b <= 0.0_dp) then
            expected = -logtwo+b+log(1.0_dp-exp(a-b))
        else if (a >= 0.0_dp) then
            expected = -logtwo-a+log(1.0_dp-exp(a-b))
        else
            expected = log(1.0_dp-0.5_dp*exp(a)-0.5_dp*exp(-b))
        end if
        call close(base,expected,8.0e-16_dp,'ordinary interval formula')
        call interval(a-0.1_dp,b+0.1_dp)
        call check(ll(data) >= base .and. ll(data) <= 0.0_dp,'nested masses')
    end do
    call interval(-0.25_dp,0.75_dp)
    call close(ll(data),log(1.0_dp-0.5_dp*exp(-0.25_dp)-0.5_dp*exp(-0.75_dp)), &
        5.0e-16_dp,'cusp-crossing interval')
    call interval(0.0_dp,tiny_width)
    call close(ll(data),log(tiny_width)-logtwo,2.0e-13_dp,'subnormal width')
    call interval(-tiny_width,tiny_width)
    call close(ll(data),log(tiny_width),2.0e-13_dp,'subnormal cusp-crossing width')
    call interval(1.0e-300_dp,2.0e-300_dp)
    call close(ll(data,scale=1.0e300_dp),log(1.0e-300_dp)-log(1.0e300_dp)-logtwo, &
        3.0e-13_dp,'standardized width underflows')
    call interval(-1.0e-300_dp,1.0e-300_dp)
    call close(ll(data,scale=1.0e300_dp),log(1.0e-300_dp)-log(1.0e300_dp), &
        3.0e-13_dp,'crossing standardized widths underflow')
    a = 1000.0_dp
    b = nearest(a,1.0_dp)
    call interval(a,b)
    expected = -a-logtwo+log(b-a)
    call close(ll(data),expected,3.0e-13_dp,'adjacent far right tail')
    call interval(-b,-a)
    call close(ll(data),expected,3.0e-13_dp,'adjacent far left tail')
    call interval(-1.0e308_dp,1.0e308_dp)
    call close(ll(data),0.0_dp,0.0_dp,'overflowing raw width')
    call close(ll(data,scale=1.0e308_dp),log(1.0_dp-exp(-1.0_dp)), &
        4.0e-14_dp,'overflowing raw width finite scaled width')
    data = make_censored_data(uncensored=[1.0e308_dp])
    call close(ll(data,-1.0e308_dp,1.0e308_dp),-2.0_dp-logtwo-log(1.0e308_dp), &
        2.0e-13_dp,'overflowing raw displacement finite standardized value')
    data = make_censored_data(uncensored=[0.0_dp])
    call close(ll(data,scale=tiny_width),-logtwo-log(tiny_width),2.0e-13_dp,'smallest scale at cusp')
    data = make_censored_data(left=[-1.0e308_dp])
    call check(ll(data,scale=1.0e-300_dp) == -inf,'unrepresentable log CDF')
    data = make_censored_data(right=[1.0e308_dp])
    call check(ll(data,scale=1.0e-300_dp) == -inf,'unrepresentable log SF')
    data = make_censored_data(uncensored=[1.0e308_dp])
    call check(ll(data,scale=1.0e-300_dp) == -inf,'unrepresentable log density')
    call interval(-inf,0.4_dp)
    call close(ll(data,0.4_dp,1.3_dp),-logtwo,0.0_dp,'left infinite endpoint')
    call interval(0.4_dp,inf)
    call close(ll(data,0.4_dp,1.3_dp),-logtwo,0.0_dp,'right infinite endpoint')
    call interval(0.4_dp,0.4_dp)
    call close(ll(data,0.4_dp,1.3_dp),laplace_logpdf(0.4_dp,0.4_dp,1.3_dp), &
        0.0_dp,'equal endpoints mean exact observation')
    call close(ll(empty),0.0_dp,0.0_dp,'empty valid sample')
    call check(ieee_is_nan(ll(empty,scale=0.0_dp)),'invalid scale even for empty sample')
    call check(ieee_is_nan(ll(data,scale=-1.0_dp)),'negative scale')
    call check(ieee_is_nan(ll(data,scale=inf)),'infinite scale')
    call check(ieee_is_nan(ll(data,loc=inf)),'infinite center')
    call check(ieee_is_nan(ll(data,scale=nan)),'NaN scale')
    call check(ieee_is_nan(ll(data,loc=nan)),'NaN center')
    data = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(ll(data)),'NaN exact observation')
    data = make_censored_data(uncensored=[inf])
    call check(ieee_is_nan(ll(data)),'infinite exact observation')
    call interval(1.0_dp,0.0_dp)
    call check(ieee_is_nan(ll(data)),'reversed interval')
    call interval(-inf,inf)
    call check(ieee_is_nan(ll(data)),'two infinite endpoints rejected')
    call laplace_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'invalid fit data')
    call laplace_fit(empty,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit sample')
    data = make_censored_data(uncensored=[0.0_dp])
    call laplace_fit(data,[0.0_dp,0.0_dp],[0.0_dp,0.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'invalid fixed scale objective')
    call laplace_fit(data,[1.0_dp,0.1_dp],[-1.0_dp,2.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'reversed fit bounds')
    call laplace_fit(data,[-1.0_dp,0.1_dp],[1.0_dp,2.0_dp],fit,tolerance=0.0_dp)
    call check(fit%status == fit_status_invalid_input,'invalid fit tolerance')
    data = make_censored_data(right=[3.0_dp])
    call laplace_fit(data,[-5.0_dp,1.0_dp],[8.0_dp,1.0_dp],fit)
    call check(fit%success,'all-right-censored bounded fit')
    call close(fit%params(1),8.0_dp,0.0_dp,'all-right fit reaches upper center bound')
    write(*,'(a,i0,a,i0)') 'Laplace censored checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    pure function ll(sample,loc,scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = laplace_censored_loglikelihood(sample,loc,scale)
    end function ll
    subroutine interval(lower,upper)
        implicit none
        real(dp), intent(in) :: lower !! lower endpoint for the test observation
        real(dp), intent(in) :: upper !! upper endpoint for the test observation
        data = make_censored_data(interval=reshape([lower,upper],[1,2]))
    end subroutine interval
    subroutine check(condition,label)
        implicit none
        logical, intent(in) :: condition !! true when the expected behavior holds
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
        real(dp), intent(in) :: tolerance !! absolute error allowed for this regime
        character(*), intent(in) :: label !! test description for failure diagnostics
        call check(ieee_is_finite(value),label//' is finite')
        if (.not. ieee_is_finite(value)) return
        call check(abs(value-reference) <= tolerance,label)
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_laplace
