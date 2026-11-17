! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_transform_hazard
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fatiguelife_censored_loglikelihood, fatiguelife_censored_nnlf, fatiguelife_fit, &
        fatiguelife_logpdf, fatiguelife_logcdf, fatiguelife_logsf, fatiguelife_ppf, &
        johnsonsu_censored_loglikelihood, johnsonsu_censored_nnlf, johnsonsu_fit, &
        johnsonsu_logpdf, johnsonsu_logcdf, johnsonsu_logsf, johnsonsu_ppf, &
        genexpon_censored_loglikelihood, genexpon_censored_nnlf, genexpon_fit, &
        genexpon_logpdf, genexpon_logcdf, genexpon_logsf, genexpon_ppf, &
        exponential_censored_loglikelihood, fit_status_invalid_input, fit_status_max_iter, &
        fit_status_no_finite_objective
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty, other
    type(fit_result) :: fit, old
    real(dp) :: actual, expected, c, x, p, previous, lower(3), upper(3), inf, nan
    real(dp), parameter :: references(3) = &
        [-14.670431355890935_dp,-42.437977289016104_dp,-23.52587064286926_dp]
    real(dp), parameter :: scipy_shapes(3) = &
        [1.2929687500000007_dp,-1.2559570312500017_dp,0.15458984374999926_dp]
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
            expected = fatiguelife_censored_nnlf(data,c,scale=1.3_dp)
        case (2)
            expected = johnsonsu_censored_nnlf(data,c,1.2_dp,scale=1.3_dp)
        case (3)
            expected = genexpon_censored_nnlf(data,c,1.2_dp,0.8_dp,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = real(j,dp)/2.0_dp
            if (family == 2) c = real(j-2,dp)
            previous = -inf
            do i = 1, 15
                x = exp(real(i-8,dp)/2.0_dp)
                if (family == 2) x = real(i-8,dp)/2.0_dp
                one = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = fatiguelife_logpdf(x,c)
                case (2)
                    expected = johnsonsu_logpdf(x,c,1.2_dp)
                case (3)
                    expected = genexpon_logpdf(x,c,1.2_dp,0.8_dp)
                end select
                call close(likelihood(one,family,c),expected,7.0e-14_dp,'scalar density')
                one = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = fatiguelife_logcdf(x,c)
                case (2)
                    expected = johnsonsu_logcdf(x,c,1.2_dp)
                case (3)
                    expected = genexpon_logcdf(x,c,1.2_dp,0.8_dp)
                end select
                actual = likelihood(one,family,c)
                call close(actual,expected,7.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = fatiguelife_logsf(x,c)
                case (2)
                    expected = johnsonsu_logsf(x,c,1.2_dp)
                case (3)
                    expected = genexpon_logsf(x,c,1.2_dp,0.8_dp)
                end select
                call close(likelihood(one,family,c),expected,7.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = fatiguelife_ppf(p,c)
                case (2)
                    x = johnsonsu_ppf(p,c,1.2_dp)
                case (3)
                    x = genexpon_ppf(p,c,1.2_dp,0.8_dp)
                end select
                one = make_censored_data(left=[x])
                call close(exp(likelihood(one,family,c)),p,8.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        c = 1.7_dp
        call close(likelihood(empty,family,c),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        if (family /= 2) call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'invalid shape')
        if (family /= 2) call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
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
            call fatiguelife_fit(one%uncensored,lower,upper,old)
        case (2)
            call johnsonsu_fit(one%uncensored,[c,1.2_dp,0.0_dp,1.3_dp], &
                [c,1.2_dp,0.0_dp,1.3_dp],old)
        case (3)
            call genexpon_fit(one%uncensored,[c,1.2_dp,0.8_dp,0.0_dp,1.3_dp], &
                [c,1.2_dp,0.8_dp,0.0_dp,1.3_dp],old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,3.0e-13_dp,'original array fit compatibility')
        call close(fit%nllf,-likelihood(one,family,c,scale=1.3_dp),0.0_dp,'objective')
        lower(1) = 0.05_dp
        upper(1) = 3.0_dp
        if (family == 2) lower(1) = -3.0_dp
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
    do i = 1, 8
        x = exp(real(i-4,dp))
        one = make_censored_data(left=[x])
        other = make_censored_data(right=[1.0_dp/x])
        call close(likelihood(one,1,1.7_dp),likelihood(other,1,1.7_dp), &
            1.0e-13_dp,'fatigue reciprocal tail identity')
        one = make_censored_data(left=[real(i-4,dp)])
        other = make_censored_data(right=[-real(i-4,dp)])
        call close(likelihood(one,2,1.7_dp),likelihood(other,2,-1.7_dp), &
            1.0e-13_dp,'Johnson reflection identity')
        one = make_censored_data(uncensored=[x],left=[0.1_dp],right=[x+0.25_dp], &
            interval=reshape([x,x+0.125_dp],[1,2]))
        call close(genexpon_censored_loglikelihood(one,2.0_dp,1.0e-300_dp,0.8_dp), &
            exponential_censored_loglikelihood(one,scale=0.5_dp), &
            2.0e-13_dp,'vanishing second rate exponential limit')
    end do
    call check(ieee_is_nan(johnsonsu_censored_loglikelihood(empty,0.0_dp,0.0_dp)), &
        'invalid second Johnson shape')
    call check(ieee_is_nan(genexpon_censored_loglikelihood(empty,1.0_dp,nan,1.0_dp)), &
        'NaN second rate')
    call check(ieee_is_nan(genexpon_censored_loglikelihood(empty,1.0_dp,1.0_dp,inf)), &
        'infinite third rate')
    call check(ieee_is_nan(genexpon_censored_loglikelihood(empty,1.0_dp,1.0_dp,0.0_dp)), &
        'zero third rate')
    one = make_censored_data(uncensored=[0.0_dp])
    call check(likelihood(one,1,1.0_dp) == -inf,'fatigue endpoint density')
    call close(genexpon_censored_loglikelihood(one,2.0_dp,1.2_dp,0.8_dp), &
        log(2.0_dp),0.0_dp,'generalized exponential endpoint density')
    do family = 1, 3, 2
        one = make_censored_data(left=[0.0_dp])
        call check(likelihood(one,family,1.0_dp) == -inf,'zero endpoint CDF')
        one = make_censored_data(right=[-1.0_dp])
        call close(likelihood(one,family,1.0_dp),0.0_dp,0.0_dp,'below support SF')
        one = make_censored_data(interval=reshape([-1.0_dp,0.0_dp],[1,2]))
        call check(likelihood(one,family,1.0_dp) == -inf,'below support interval')
        one = make_censored_data(uncensored=[-1.0_dp])
        lower = [2.0_dp,0.0_dp,1.0_dp]
        call fitting(one,family,lower,lower,fit)
        call check(fit%status == fit_status_no_finite_objective,'impossible fixed objective')
    end do
    one = make_censored_data(uncensored=[1.0_dp])
    actual = fatiguelife_censored_loglikelihood(one,1.0e-300_dp)
    expected = -log(1.0e-300_dp)-0.5_dp*log(2.0_dp*acos(-1.0_dp))
    call close(actual,expected,3.0e-13_dp,'tiny fatigue shape at scale')
    one = make_censored_data(right=[1.0e-300_dp])
    actual = genexpon_censored_loglikelihood(one,1.0_dp,1.0_dp,1.0_dp,scale=1.0e300_dp)
    call check(actual == 0.0_dp,'underflowing cumulative hazard')
    one = make_censored_data(left=[1.0e-300_dp])
    actual = genexpon_censored_loglikelihood(one,1.0_dp,1.0_dp,1.0_dp,scale=1.0e300_dp)
    call close(actual,log(1.0e-300_dp)-log(1.0e300_dp),3.0e-13_dp, &
        'log CDF retains underflowing distance')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(johnsonsu_censored_loglikelihood(one,0.0_dp,1.0e300_dp), &
        log(1.0e300_dp)-0.5_dp*log(2.0_dp*acos(-1.0_dp)),3.0e-13_dp, &
        'large Johnson scale shape at center')
    one = make_censored_data(uncensored=[1.0_dp])
    call close(genexpon_censored_loglikelihood(one,1.0e-300_dp,1.0e300_dp,1.0e-300_dp), &
        -0.5_dp,3.0e-13_dp,'opposing generalized exponential rates')
    one = make_censored_data(interval=reshape([1.0_dp,nearest(1.0_dp,1.0_dp)],[1,2]))
    expected = -0.5_dp+log(nearest(1.0_dp,1.0_dp)-1.0_dp)
    call close(genexpon_censored_loglikelihood(one,1.0e-300_dp,1.0e300_dp,1.0e-300_dp), &
        expected,3.0e-13_dp,'opposing rates adjacent interval')
    one = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    expected = log(1.0e300_dp)+log(nearest(0.0_dp,1.0_dp)) &
        -0.5_dp*log(2.0_dp*acos(-1.0_dp))
    call close(johnsonsu_censored_loglikelihood(one,0.0_dp,1.0e300_dp), &
        expected,3.0e-13_dp,'subnormal Johnson interval with huge shape')
    one = make_censored_data(left=[sinh(2.0_dp)])
    call close(johnsonsu_censored_loglikelihood(one,-1.0e308_dp,1.0e308_dp), &
        0.0_dp,0.0_dp,'opposing Johnson shape rescues overflowing product')
    one = make_censored_data(right=[1.0e-300_dp])
    actual = genexpon_censored_loglikelihood(one,1.0e308_dp,1.0e308_dp,0.8_dp)
    call close(actual/(-1.0e8_dp),1.0_dp,3.0e-13_dp,'overflowing total rate retained')
    call johnsonsu_fit(data,[-3.0_dp,0.5_dp,-0.2_dp,0.8_dp], &
        [3.0_dp,2.0_dp,0.1_dp,2.0_dp],fit,max_iter=10000,tolerance=1.0e-6_dp)
    call check(fit%success,'four free parameter fit')
    call close(fit%nllf,johnsonsu_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4)),0.0_dp,'four parameter fit objective')
    call genexpon_fit(data,[0.1_dp,0.5_dp,0.3_dp,-0.2_dp,0.8_dp], &
        [3.0_dp,2.0_dp,2.0_dp,0.1_dp,2.0_dp],fit,max_iter=10000,tolerance=1.0e-6_dp)
    call check(fit%success,'five free parameter fit')
    call close(fit%nllf,genexpon_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4),fit%params(5)),0.0_dp,'five parameter fit objective')
    print '(a,i0,a)', 'Censored transform/hazard checks: ', checks, ' passed'
    if (failures /= 0) error stop 'censored transform/hazard failures'
contains
    function likelihood(sample, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 fatiguelife, 2 johnsonsu, 3 genexpon
        real(dp), intent(in) :: shape !! finite first shape; positive except Johnson SU
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = fatiguelife_censored_loglikelihood(sample,shape,loc,scale)
        case (2)
            value = johnsonsu_censored_loglikelihood(sample,shape,1.2_dp,loc,scale)
        case (3)
            value = genexpon_censored_loglikelihood(sample,shape,1.2_dp,0.8_dp,loc,scale)
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
            call fatiguelife_fit(sample,lower,upper,result,max_iter=max_iter)
        case (2)
            call johnsonsu_fit(sample,[lower(1),1.2_dp,lower(2:3)], &
                [upper(1),1.2_dp,upper(2:3)],result,max_iter=max_iter)
        case (3)
            call genexpon_fit(sample,[lower(1),1.2_dp,0.8_dp,lower(2:3)], &
                [upper(1),1.2_dp,0.8_dp,upper(2:3)],result,max_iter=max_iter)
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
end program test_censored_transform_hazard
