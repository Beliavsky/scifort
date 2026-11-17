! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_lifetime
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        exponpow_censored_loglikelihood, exponpow_censored_nnlf, exponpow_fit, &
        exponpow_logpdf, exponpow_logcdf, exponpow_logsf, exponpow_ppf, &
        exponweib_censored_loglikelihood, exponweib_censored_nnlf, exponweib_fit, &
        exponweib_logpdf, exponweib_logcdf, exponweib_logsf, exponweib_ppf, &
        mielke_censored_loglikelihood, mielke_censored_nnlf, mielke_fit, &
        mielke_logpdf, mielke_logcdf, mielke_logsf, mielke_ppf, &
        gompertz_censored_loglikelihood, weibull_censored_loglikelihood, &
        fisk_censored_loglikelihood, fit_status_invalid_input, fit_status_max_iter, &
        fit_status_no_finite_objective
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty, other
    type(fit_result) :: fit, old
    real(dp) :: actual, expected, c, x, p, previous, lower(3), upper(3), inf, nan
    real(dp), parameter :: references(3) = &
        [-1097.7750366596047_dp,-15.776394603077286_dp,-14.709016801278233_dp]
    real(dp), parameter :: scipy_shapes(3) = &
        [0.3581054687499994_dp,1.3926757812500012_dp,1.302636718750001_dp]
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    do family = 1, 3
        c = 1.7_dp
        actual = likelihood(data,family,c,scale=1.3_dp)
        call close(actual,references(family),3.0e-13_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = exponpow_censored_nnlf(data,c,scale=1.3_dp)
        case (2)
            expected = exponweib_censored_nnlf(data,c,1.2_dp,scale=1.3_dp)
        case (3)
            expected = mielke_censored_nnlf(data,c,1.2_dp,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = real(j,dp)/2.0_dp
            previous = -inf
            do i = 1, 15
                x = exp(real(i-8,dp)/2.0_dp)
                one = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = exponpow_logpdf(x,c)
                case (2)
                    expected = exponweib_logpdf(x,c,1.2_dp)
                case (3)
                    expected = mielke_logpdf(x,c,1.2_dp)
                end select
                call close(likelihood(one,family,c),expected,7.0e-14_dp,'scalar density')
                one = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = exponpow_logcdf(x,c)
                case (2)
                    expected = exponweib_logcdf(x,c,1.2_dp)
                case (3)
                    expected = mielke_logcdf(x,c,1.2_dp)
                end select
                actual = likelihood(one,family,c)
                call close(actual,expected,7.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = exponpow_logsf(x,c)
                case (2)
                    expected = exponweib_logsf(x,c,1.2_dp)
                case (3)
                    expected = mielke_logsf(x,c,1.2_dp)
                end select
                call close(likelihood(one,family,c),expected,7.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = exponpow_ppf(p,c)
                case (2)
                    x = exponweib_ppf(p,c,1.2_dp)
                case (3)
                    x = mielke_ppf(p,c,1.2_dp)
                end select
                one = make_censored_data(left=[x])
                call close(exp(likelihood(one,family,c)),p,8.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        c = 1.7_dp
        call close(likelihood(empty,family,c),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'invalid shape')
        call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero power')
        call check(ieee_is_nan(likelihood(empty,family,c,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,c,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=-1.0_dp)),'negative scale')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=inf)),'infinite scale')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=nan)),'NaN scale')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(one,family,c)),'NaN exact observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(one,family,c)),'infinite exact observation')
        one = make_censored_data(interval=reshape([1.0_dp,0.0_dp],[1,2]))
        call check(ieee_is_nan(likelihood(one,family,c)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,1.0_dp],[1,2]))
        other = make_censored_data(left=[1.0_dp])
        call close(likelihood(one,family,c),likelihood(other,family,c),0.0_dp,'left interval')
        one = make_censored_data(interval=reshape([1.0_dp,inf],[1,2]))
        other = make_censored_data(right=[1.0_dp])
        call close(likelihood(one,family,c),likelihood(other,family,c),0.0_dp,'right interval')
        one = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
        other = make_censored_data(uncensored=[1.0_dp])
        call close(likelihood(one,family,c),likelihood(other,family,c),0.0_dp,'equal endpoints')
        lower = [c,0.0_dp,1.3_dp]
        upper = lower
        one = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
        call fitting(one,family,lower,upper,fit)
        select case (family)
        case (1)
            call exponpow_fit(one%uncensored,lower,upper,old)
        case (2)
            call exponweib_fit(one%uncensored,[c,1.2_dp,0.0_dp,1.3_dp], &
                [c,1.2_dp,0.0_dp,1.3_dp],old)
        case (3)
            call mielke_fit(one%uncensored,[c,1.2_dp,0.0_dp,1.3_dp], &
                [c,1.2_dp,0.0_dp,1.3_dp],old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,3.0e-13_dp,'original array fit compatibility')
        call close(fit%nllf,-likelihood(one,family,c,scale=1.3_dp),0.0_dp,'objective')
        lower(1) = 0.05_dp
        upper(1) = 3.0_dp
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'shape fit converges')
        call close(fit%params(1),scipy_shapes(family),2.0e-4_dp,'SciPy constrained shape fit')
        call check(fit%nllf <= -likelihood(data,family,c,scale=1.3_dp),'fit improves objective')
        lower(2:3) = [-0.2_dp,1.2_dp]
        upper(2:3) = [0.1_dp,2.0_dp]
        call fitting(data,family,lower,upper,fit,max_iter=10000)
        call check(fit%success,'free location and scale fit')
        call fitting(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        call fitting(data,family,upper,lower,fit)
        call check(fit%status == fit_status_invalid_input,'reversed bounds rejected')
        call fitting(data,family,lower,upper,fit,max_iter=1)
        call check(fit%status == fit_status_max_iter,'iteration limit')
    end do
    do i = 1, 6
        x = exp(real(i-3,dp))
        one = make_censored_data(left=[x],right=[x+0.25_dp], &
            interval=reshape([x,x+0.125_dp],[1,2]))
        call close(likelihood(one,1,1.0_dp),gompertz_censored_loglikelihood(one,1.0_dp), &
            1.0e-13_dp,'unit power Gompertz identity')
        call close(likelihood(one,2,1.0_dp),weibull_censored_loglikelihood(one,1.2_dp), &
            1.0e-13_dp,'unit exponent Weibull identity')
        call close(likelihood(one,3,1.2_dp),fisk_censored_loglikelihood(one,1.2_dp), &
            1.0e-13_dp,'equal Mielke shapes Fisk identity')
    end do
    call check(ieee_is_nan(exponweib_censored_loglikelihood(empty,1.0_dp,0.0_dp)), &
        'invalid second shape')
    call check(ieee_is_nan(mielke_censored_loglikelihood(empty,1.0_dp,nan)), &
        'NaN second shape')
    call check(ieee_is_nan(mielke_censored_loglikelihood(empty,1.0_dp,inf)), &
        'infinite second shape')
    do family = 1, 3
        one = make_censored_data(uncensored=[0.0_dp])
        call check(likelihood(one,family,0.5_dp) == inf,'singular endpoint density')
        call check(likelihood(one,family,2.0_dp) == -inf,'zero endpoint density')
        one = make_censored_data(left=[0.0_dp])
        call check(likelihood(one,family,1.0_dp) == -inf,'zero endpoint CDF')
        one = make_censored_data(right=[-1.0_dp])
        call close(likelihood(one,family,1.0_dp),0.0_dp,0.0_dp,'below support SF')
        one = make_censored_data(interval=reshape([-1.0_dp,0.0_dp],[1,2]))
        call check(likelihood(one,family,1.0_dp) == -inf,'below support interval')
        one = make_censored_data(uncensored=[0.0_dp])
        lower = [2.0_dp,0.0_dp,1.0_dp]
        call fitting(one,family,lower,lower,fit)
        call check(fit%status == fit_status_no_finite_objective,'impossible fixed objective')
        lower(1) = 0.5_dp
        call fitting(one,family,lower,lower,fit)
        call check(fit%status == fit_status_no_finite_objective,'singular fixed objective')
    end do
    one = make_censored_data(left=[1.0e-300_dp])
    actual = exponweib_censored_loglikelihood(one,1.0e-300_dp,1.0e308_dp)
    expected = 1.0e8_dp*log(1.0e-300_dp)
    call close(actual/expected,1.0_dp,3.0e-13_dp,'opposing shapes rescue overflowing log power')
    actual = mielke_censored_loglikelihood(one,1.0e-300_dp,1.0e308_dp)
    expected = 1.0e-300_dp*log(1.0e-300_dp)
    call close(actual/expected,1.0_dp,3.0e-13_dp,'large Mielke shape cancellation')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(exponweib_censored_loglikelihood(one,0.5_dp,2.0_dp),0.0_dp,0.0_dp, &
        'finite unit product endpoint')
    call close(mielke_censored_loglikelihood(one,1.0_dp,2.0_dp),0.0_dp,0.0_dp, &
        'finite unit first shape endpoint')
    call close(exponweib_logsf(80.0_dp,1.0_dp,1.0_dp),-80.0_dp, &
        3.0e-14_dp,'scalar unit-shape upper tail')
    call check(exponweib_censored_loglikelihood(one,1.0e300_dp,1.0e300_dp) == -inf, &
        'overflowing shape product endpoint')
    one = make_censored_data(right=[1.0e-300_dp])
    actual = exponpow_censored_loglikelihood(one,0.4_dp,scale=1.0e300_dp)
    call close(actual/(-1.0e-240_dp),1.0_dp,3.0e-13_dp,'tiny negative log survival')
    one = make_censored_data(right=[1.0_dp])
    actual = mielke_censored_loglikelihood(one,1.0e-300_dp,1.0e300_dp)
    expected = log(1.0e-300_dp)-log(1.0e300_dp)+log(log(2.0_dp))
    call close(actual,expected,3.0e-13_dp,'underflowing Mielke shape ratio survival')
    one = make_censored_data(right=[1000.0_dp])
    actual = exponweib_censored_loglikelihood(one,1.0e300_dp,1.0_dp)
    call close(actual,log(1.0e300_dp)-1000.0_dp,3.0e-13_dp,'large shape far survival')
    call exponweib_fit(data,[0.1_dp,0.5_dp,-0.2_dp,0.8_dp], &
        [3.0_dp,2.0_dp,0.1_dp,2.0_dp],fit,max_iter=10000,tolerance=1.0e-6_dp)
    call check(fit%success,'four free parameter Weibull fit')
    call close(fit%nllf,exponweib_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4)),0.0_dp,'four parameter fit objective')
    call mielke_fit(data,[0.1_dp,0.5_dp,-0.2_dp,0.8_dp], &
        [3.0_dp,2.0_dp,0.1_dp,2.0_dp],fit,max_iter=10000,tolerance=1.0e-6_dp)
    call check(fit%success,'four free parameter Mielke fit')
    print '(a,i0,a)', 'Censored lifetime checks: ', checks, ' passed'
    if (failures /= 0) error stop 'censored lifetime failures'
contains
    function likelihood(sample, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 exponpow, 2 exponweib, 3 mielke
        real(dp), intent(in) :: shape !! positive first shape parameter
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = exponpow_censored_loglikelihood(sample,shape,loc,scale)
        case (2)
            value = exponweib_censored_loglikelihood(sample,shape,1.2_dp,loc,scale)
        case (3)
            value = mielke_censored_loglikelihood(sample,shape,1.2_dp,loc,scale)
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function likelihood
    subroutine fitting(sample, family, lower, upper, result, max_iter)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! family selector as in likelihood
        real(dp), intent(in) :: lower(3) !! lower bounds for shape, location, scale
        real(dp), intent(in) :: upper(3) !! upper bounds; equality fixes a parameter
        type(fit_result), intent(out) :: result !! fit parameters and diagnostics
        integer, intent(in), optional :: max_iter !! positive sweep limit, default 1000
        select case (family)
        case (1)
            call exponpow_fit(sample,lower,upper,result,max_iter=max_iter)
        case (2)
            call exponweib_fit(sample,[lower(1),1.2_dp,lower(2:3)], &
                [upper(1),1.2_dp,upper(2:3)],result,max_iter=max_iter)
        case (3)
            call mielke_fit(sample,[lower(1),1.2_dp,lower(2:3)], &
                [upper(1),1.2_dp,upper(2:3)],result,max_iter=max_iter)
        end select
    end subroutine fitting
    subroutine check(condition, label)
        implicit none
        logical, intent(in) :: condition !! assertion that must hold
        character(*), intent(in) :: label !! diagnostic description
        checks = checks+1
        if (condition) return
        failures = failures+1
        print '(a)', 'FAIL: '//label
    end subroutine check
    subroutine close(actual, expected, tolerance, label)
        implicit none
        real(dp), intent(in) :: actual !! computed finite value
        real(dp), intent(in) :: expected !! independently expected finite value
        real(dp), intent(in) :: tolerance !! scaled absolute tolerance
        character(*), intent(in) :: label !! diagnostic description
        call check(ieee_is_finite(actual),label//' finite')
        if (.not. ieee_is_finite(actual)) return
        call check(abs(actual-expected) <= tolerance*max(1.0_dp,abs(expected)),label)
        if (abs(actual-expected) > tolerance*max(1.0_dp,abs(expected))) &
            print '(2es25.16)', actual, expected
    end subroutine close
end program test_censored_lifetime
