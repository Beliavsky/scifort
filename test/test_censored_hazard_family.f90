! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_hazard_family
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        genlogistic_censored_loglikelihood, genlogistic_censored_nnlf, genlogistic_fit, &
        genlogistic_logpdf, genlogistic_logcdf, genlogistic_logsf, genlogistic_ppf, &
        gompertz_censored_loglikelihood, gompertz_censored_nnlf, gompertz_fit, &
        gompertz_logpdf, gompertz_logcdf, gompertz_logsf, gompertz_ppf, &
        invweibull_censored_loglikelihood, invweibull_censored_nnlf, invweibull_fit, &
        invweibull_logpdf, invweibull_logcdf, invweibull_logsf, invweibull_ppf, &
        logistic_censored_loglikelihood, weibull_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, other
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, actual, expected, x, c, p, previous, reference(3), lower(3), upper(3)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    reference = [-15.095845974501499_dp,-71.584738650955799_dp,-24.434523984882752_dp]
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    do family = 1, 3
        actual = likelihood(data,family,1.7_dp,scale=1.3_dp)
        call close(actual,reference(family),3.0e-13_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = genlogistic_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (2)
            expected = gompertz_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (3)
            expected = invweibull_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = 0.4_dp
            if (j == 2) c = 1.0_dp
            if (j == 3) c = 3.0_dp
            previous = -inf
            do i = 1, 12
                x = real(i,dp)/3.0_dp
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = genlogistic_logpdf(x,c)
                case (2)
                    expected = gompertz_logpdf(x,c)
                case (3)
                    expected = invweibull_logpdf(x,c)
                end select
                call close(likelihood(single,family,c),expected,6.0e-14_dp,'scalar density')
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = genlogistic_logcdf(x,c)
                case (2)
                    expected = gompertz_logcdf(x,c)
                case (3)
                    expected = invweibull_logcdf(x,c)
                end select
                actual = likelihood(single,family,c)
                call close(actual,expected,6.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity and range')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = genlogistic_logsf(x,c)
                case (2)
                    expected = gompertz_logsf(x,c)
                case (3)
                    expected = invweibull_logsf(x,c)
                end select
                call close(likelihood(single,family,c),expected,6.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = genlogistic_ppf(p,c)
                case (2)
                    x = gompertz_ppf(p,c)
                case (3)
                    x = invweibull_ppf(p,c)
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,c)),p,5.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        single = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
        lower = [1.7_dp,0.0_dp,1.3_dp]
        upper = lower
        call fitting(single,family,lower,upper,fit)
        select case (family)
        case (1)
            call genlogistic_fit(single%uncensored,lower,upper,old)
        case (2)
            call gompertz_fit(single%uncensored,lower,upper,old)
        case (3)
            call invweibull_fit(single%uncensored,lower,upper,old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,2.0e-13_dp,'array fit compatibility')
        call close(fit%nllf,-likelihood(single,family,1.7_dp,scale=1.3_dp), &
            0.0_dp,'fit objective agreement')
        call close(likelihood(empty,family,1.0_dp),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
        call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'negative shape')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=nan)),'NaN scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=inf)),'infinite scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=-1.0_dp)),'negative scale')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'NaN exact value')
        single = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'infinite exact value')
        single = make_censored_data(interval=reshape([1.0_dp,0.0_dp],[1,2]))
        call check(ieee_is_nan(likelihood(single,family,1.0_dp)),'reversed interval')
        single = make_censored_data(interval=reshape([-inf,1.0_dp],[1,2]))
        other = make_censored_data(left=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite lower interval endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,inf],[1,2]))
        other = make_censored_data(right=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite upper interval endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
        other = make_censored_data(uncensored=[1.0_dp])
        call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'equal endpoints represent exact data')
        if (family > 1) then
            single = make_censored_data(left=[0.0_dp])
            call check(likelihood(single,family,1.0_dp) == -inf,'zero CDF at support')
            single = make_censored_data(right=[-1.0_dp])
            call close(likelihood(single,family,1.0_dp),0.0_dp,0.0_dp,'SF below support')
            single = make_censored_data(interval=reshape([-1.0_dp,0.0_dp],[1,2]))
            call check(likelihood(single,family,1.0_dp) == -inf,'interval outside support')
            single = make_censored_data(interval=reshape([-1.0_dp,1.0_dp],[1,2]))
            other = make_censored_data(left=[1.0_dp])
            call close(likelihood(single,family,1.0_dp),likelihood(other,family,1.0_dp), &
                0.0_dp,'interval clipped to support')
            single = make_censored_data(uncensored=[0.0_dp])
            if (family == 2) then
                call close(likelihood(single,family,0.4_dp),log(0.4_dp),0.0_dp,'Gompertz endpoint')
            else
                call check(likelihood(single,family,0.4_dp) == -inf,'inverse Weibull endpoint')
            end if
        end if
        lower = [0.3_dp,-0.1_dp,0.8_dp]
        upper = [3.0_dp,0.1_dp,2.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'all-free bounded fit')
        if (allocated(fit%params)) then
            call check(all(fit%params >= lower) .and. all(fit%params <= upper),'fit within bounds')
            call close(fit%nllf,-likelihood(data,family,fit%params(1), &
                fit%params(2),fit%params(3)),0.0_dp,'free fit objective')
        end if
        lower = [1.7_dp,0.0_dp,1.3_dp]
        upper = [1.7_dp,0.0_dp,3.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'partly fixed fit')
        lower = [0.01_dp,0.0_dp,1.3_dp]
        upper = [10.0_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'fixed support shape fit')
        call check(fit%nllf <= -reference(family),'shape fit improves objective')
        select case (family)
        case (1)
            expected = 3.0464843750000044_dp
        case (2)
            expected = 0.12675781249999923_dp
        case (3)
            expected = 0.7249999999999989_dp
        end select
        call close(fit%params(1),expected,2.0e-4_dp,'SciPy fixed-support shape fit')
        call fitting(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        upper(3) = 0.1_dp
        call fitting(data,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'reversed bounds rejected')
        single = make_censored_data(uncensored=[-1.0_dp])
        if (family == 1) single = make_censored_data(uncensored=[-huge(1.0_dp)])
        lower = [2.0_dp,0.0_dp,0.01_dp]
        call fitting(single,family,lower,lower,fit)
        call check(fit%status == fit_status_no_finite_objective,'impossible fixed objective')
    end do
    call close(genlogistic_censored_loglikelihood(data,1.0_dp), &
        logistic_censored_loglikelihood(data),3.0e-14_dp,'logistic shape-one identity')
    single = make_censored_data(uncensored=[2.0_dp])
    other = make_censored_data(uncensored=[0.5_dp])
    call close(invweibull_censored_loglikelihood(single,1.7_dp), &
        weibull_censored_loglikelihood(other,1.7_dp)-2.0_dp*log(2.0_dp), &
        3.0e-14_dp,'reciprocal Weibull Jacobian')
    single = make_censored_data(interval=reshape([2.0_dp,4.0_dp],[1,2]))
    other = make_censored_data(interval=reshape([0.25_dp,0.5_dp],[1,2]))
    call close(invweibull_censored_loglikelihood(single,1.7_dp), &
        weibull_censored_loglikelihood(other,1.7_dp),3.0e-14_dp,'reciprocal interval identity')
    single = make_censored_data(left=[1000.0_dp])
    expected = -1.0e300_dp*exp(-1000.0_dp/2.0_dp)*exp(-1000.0_dp/2.0_dp)
    actual = genlogistic_censored_loglikelihood(single,1.0e300_dp)
    call close(actual/expected,1.0_dp,2.0e-13_dp,'large shape restores tiny softplus')
    single = make_censored_data(left=[-1.0e300_dp])
    actual = genlogistic_censored_loglikelihood(single,1.0e-300_dp,scale=1.0e-300_dp)
    call close(actual/(-1.0e300_dp),1.0_dp,2.0e-13_dp,'small shape restores huge coordinate')
    single = make_censored_data(right=[1.0e-300_dp])
    actual = gompertz_censored_loglikelihood(single,1.0e300_dp,scale=1.0e300_dp)
    call close(actual/(-1.0e-300_dp),1.0_dp,2.0e-13_dp,'underflowed Gompertz distance')
    single = make_censored_data(right=[1000.0_dp])
    expected = -exp(1000.0_dp+log(1.0e-300_dp))
    actual = gompertz_censored_loglikelihood(single,1.0e-300_dp)
    call close(actual/expected,1.0_dp,2.0e-13_dp,'small shape restores huge exponential')
    single = make_censored_data(right=[1.0e300_dp])
    call close(invweibull_censored_loglikelihood(single,2.0_dp), &
        -600.0_dp*log(10.0_dp),3.0e-13_dp,'finite log SF after power underflows')
    single = make_censored_data(left=[1.0e300_dp])
    actual = invweibull_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-1.0e-300_dp),1.0_dp,2.0e-13_dp,'near-unit inverse Weibull CDF')
    single = make_censored_data(interval=reshape([-100.0_dp,100.0_dp],[1,2]))
    actual = genlogistic_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-2.0_dp*exp(-100.0_dp)),1.0_dp,2.0e-13_dp, &
        'nearly full generalized logistic interval ratio')
    single = make_censored_data(interval=reshape([1.0e-100_dp,1000.0_dp],[1,2]))
    actual = gompertz_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-1.0e-100_dp),1.0_dp,2.0e-13_dp,'nearly full Gompertz interval ratio')
    single = make_censored_data(interval=reshape([1.0e-100_dp,1.0e100_dp],[1,2]))
    actual = invweibull_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-1.0e-100_dp),1.0_dp,2.0e-13_dp,'nearly full inverse Weibull interval ratio')
    single = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    call close(genlogistic_censored_loglikelihood(single,1.0_dp), &
        log(nearest(0.0_dp,1.0_dp))-log(4.0_dp),3.0e-13_dp,'smallest logistic width')
    single = make_censored_data(right=[1.0_dp])
    call close(genlogistic_censored_loglikelihood(single,nearest(0.0_dp,1.0_dp)), &
        log(nearest(0.0_dp,1.0_dp))+log(log(1.0_dp+exp(-1.0_dp))), &
        3.0e-13_dp,'smallest shape generalized logistic SF')
    single = make_censored_data(uncensored=[1.0e300_dp])
    call check(invweibull_censored_loglikelihood(single,huge(1.0_dp)) == -inf, &
        'inverse Weibull unrepresentable log power')
    single = make_censored_data(uncensored=[1.0e-300_dp])
    call check(gompertz_censored_loglikelihood(single,1.0_dp,scale=nearest(0.0_dp,1.0_dp)) &
        == -inf,'Gompertz unrepresentable standardized coordinate')
    write(*,'(a,i0,a)') 'Hazard family censored: ',checks,' checks'
    if (failures /= 0) error stop 'Hazard family censoring failures'
contains
    pure function likelihood(data, family, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! 1 generalized logistic, 2 Gompertz, 3 inverse Weibull
        real(dp), intent(in) :: c !! positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = genlogistic_censored_loglikelihood(data,c,loc,scale)
        case (2)
            value = gompertz_censored_loglikelihood(data,c,loc,scale)
        case (3)
            value = invweibull_censored_loglikelihood(data,c,loc,scale)
        end select
    end function likelihood
    subroutine fitting(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! distribution identifier used by this test
        real(dp), intent(in) :: lower(3) !! shape, location, scale lower bounds
        real(dp), intent(in) :: upper(3) !! shape, location, scale upper bounds
        type(fit_result), intent(out) :: result !! bounded fit diagnostics
        select case (family)
        case (1)
            call genlogistic_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (2)
            call gompertz_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (3)
            call invweibull_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        end select
    end subroutine fitting
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
end program test_censored_hazard_family
