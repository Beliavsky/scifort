! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_fold_pareto
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        foldnorm_censored_loglikelihood, foldnorm_censored_nnlf, foldnorm_fit, &
        foldnorm_logpdf, foldnorm_logcdf, foldnorm_logsf, foldnorm_ppf, &
        foldcauchy_censored_loglikelihood, foldcauchy_censored_nnlf, foldcauchy_fit, &
        foldcauchy_logpdf, foldcauchy_logcdf, foldcauchy_logsf, foldcauchy_ppf, &
        genpareto_censored_loglikelihood, genpareto_censored_nnlf, genpareto_fit, &
        genpareto_logpdf, genpareto_logcdf, genpareto_logsf, genpareto_ppf, &
        exponential_censored_loglikelihood, uniform_censored_loglikelihood, &
        lomax_censored_loglikelihood, halfnorm_censored_loglikelihood, &
        halfcauchy_censored_loglikelihood, fit_status_invalid_input
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, other, empty
    type(fit_result) :: fit, old
    real(dp) :: x, p, c, actual, expected, previous, inf, nan, lower(3), upper(3), t
    integer :: family, i, j, checks, failures
    real(dp), parameter :: scipy_nllf(3) = &
        [11.332155282192696_dp,12.774437362573709_dp,9.9482779792704878_dp]
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        do j = 1, 3
            c = real(j-1,dp)/2.0_dp
            if (family == 3) c = real(j-2,dp)/2.0_dp
            previous = -inf
            do i = 1, 19
                x = real(i,dp)/20.0_dp
                select case (family)
                case (1)
                    expected = foldnorm_logpdf(x,c)
                case (2)
                    expected = foldcauchy_logpdf(x,c)
                case (3)
                    expected = genpareto_logpdf(x,c)
                end select
                one = make_censored_data(uncensored=[x])
                call close(likelihood(one,family,c),expected,2.0e-12_dp,'scalar density')
                one = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = foldnorm_logcdf(x,c)
                case (2)
                    expected = foldcauchy_logcdf(x,c)
                case (3)
                    expected = genpareto_logcdf(x,c)
                end select
                actual = likelihood(one,family,c)
                call close(actual,expected,2.0e-12_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = foldnorm_logsf(x,c)
                case (2)
                    expected = foldcauchy_logsf(x,c)
                case (3)
                    expected = genpareto_logsf(x,c)
                end select
                call close(likelihood(one,family,c),expected,2.0e-12_dp,'scalar SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                other = make_censored_data(uncensored=[x])
                expected = likelihood(other,family,c)+log(nearest(x,1.0_dp)-x)
                call close(likelihood(one,family,c),expected,3.0e-13_dp,'adjacent interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = foldnorm_ppf(p,c)
                case (2)
                    x = foldcauchy_ppf(p,c)
                case (3)
                    x = genpareto_ppf(p,c)
                end select
                one = make_censored_data(left=[x])
                call close(likelihood(one,family,c),log(p),3.0e-12_dp,'CDF/PPF round trip')
            end do
        end do
        call close(likelihood(empty,family,0.0_dp),0.0_dp,0.0_dp,'empty sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        if (family < 3) call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)), &
            'negative folding offset')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=nan)),'NaN location')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(one,family,1.0_dp)),'NaN observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(one,family,1.0_dp)),'infinite observation')
        one = make_censored_data(uncensored=[-1.0_dp])
        call check(likelihood(one,family,1.0_dp) == -inf,'outside support density')
        one = make_censored_data(left=[0.0_dp])
        call check(likelihood(one,family,1.0_dp) == -inf,'lower endpoint CDF')
        one = make_censored_data(right=[0.0_dp])
        call close(likelihood(one,family,1.0_dp),0.0_dp,0.0_dp,'lower endpoint SF')
        one = make_censored_data(interval=reshape([-inf,0.3_dp],[1,2]))
        other = make_censored_data(left=[0.3_dp])
        call close(likelihood(one,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'infinite interval endpoint')
        one = make_censored_data(interval=reshape([0.3_dp,0.3_dp],[1,2]))
        other = make_censored_data(uncensored=[0.3_dp])
        call close(likelihood(one,family,1.0_dp),likelihood(other,family,1.0_dp), &
            0.0_dp,'equal interval endpoints')
        lower = [0.5_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,lower,fit)
        call check(fit%success,'fixed fit')
        call close(fit%nllf,-likelihood(data,family,0.5_dp,scale=1.3_dp),0.0_dp, &
            'fixed fit objective')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fitting(one,family,lower,lower,fit)
        select case (family)
        case (1)
            call foldnorm_fit(one%uncensored,lower,lower,old)
            actual = foldnorm_censored_nnlf(data,0.5_dp,scale=1.3_dp)
        case (2)
            call foldcauchy_fit(one%uncensored,lower,lower,old)
            actual = foldcauchy_censored_nnlf(data,0.5_dp,scale=1.3_dp)
        case (3)
            call genpareto_fit(one%uncensored,lower,lower,old)
            actual = genpareto_censored_nnlf(data,0.5_dp,scale=1.3_dp)
        end select
        call check(old%success,'original array fit')
        call close(fit%nllf,old%nllf,2.0e-12_dp,'array/censored fit agreement')
        call close(actual,-likelihood(data,family,0.5_dp,scale=1.3_dp),0.0_dp,'NNLF')
        call fitting(empty,family,lower,lower,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        lower(1) = 0.0_dp
        if (family == 3) lower(1) = -1.0_dp
        upper = [3.0_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'free shape fit')
        call close(fit%nllf,-likelihood(data,family,fit%params(1),scale=1.3_dp), &
            0.0_dp,'free fit objective')
        call close(fit%nllf,scipy_nllf(family),2.0e-8_dp,'SciPy optimal objective')
        expected = 0.0_dp
        if (family == 3) expected = -1.0_dp
        call close(fit%params(1),expected,3.0e-5_dp,'SciPy boundary optimum')
        lower = [0.0_dp,-0.1_dp,0.9_dp]
        if (family == 3) lower(1) = -1.0_dp
        upper = [3.0_dp,0.1_dp,2.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'three free parameters')
        call close(fit%nllf,-likelihood(data,family,fit%params(1), &
            fit%params(2),fit%params(3)),0.0_dp,'three parameter objective')
    end do
    call close(genpareto_censored_loglikelihood(data,0.0_dp), &
        exponential_censored_loglikelihood(data),3.0e-13_dp,'exponential identity')
    call close(genpareto_censored_loglikelihood(data,-1.0_dp), &
        uniform_censored_loglikelihood(data),3.0e-13_dp,'uniform identity')
    call close(genpareto_censored_loglikelihood(data,1.0_dp), &
        lomax_censored_loglikelihood(data,1.0_dp),3.0e-13_dp,'Lomax identity')
    call close(foldnorm_censored_loglikelihood(data,0.0_dp), &
        halfnorm_censored_loglikelihood(data),3.0e-13_dp,'half-normal identity')
    call close(foldcauchy_censored_loglikelihood(data,0.0_dp), &
        halfcauchy_censored_loglikelihood(data),3.0e-13_dp,'half-Cauchy identity')
    one = make_censored_data(uncensored=[1.0_dp])
    call check(genpareto_censored_loglikelihood(one,-2.0_dp,scale=2.0_dp) == inf, &
        'singular Pareto upper density')
    call check(genpareto_censored_loglikelihood(one,-0.5_dp,scale=0.5_dp) == -inf, &
        'vanishing Pareto upper density')
    call close(genpareto_censored_loglikelihood(one,-1.0_dp),0.0_dp,0.0_dp, &
        'finite Pareto upper density')
    one = make_censored_data(left=[1.0_dp])
    call close(genpareto_censored_loglikelihood(one,-1.0_dp),0.0_dp,0.0_dp, &
        'upper endpoint CDF')
    one = make_censored_data(right=[1.0_dp])
    call check(genpareto_censored_loglikelihood(one,-1.0_dp) == -inf,'upper endpoint SF')
    one = make_censored_data(interval=reshape([-1.0_dp,2.0_dp],[1,2]))
    call close(genpareto_censored_loglikelihood(one,-1.0_dp),0.0_dp,0.0_dp, &
        'whole bounded support')
    one = make_censored_data(interval=reshape([1.0_dp,2.0_dp],[1,2]))
    call check(genpareto_censored_loglikelihood(one,-1.0_dp) == -inf,'outside bounded interval')
    call check(ieee_is_nan(genpareto_censored_loglikelihood(empty,-0.5_dp, &
        scale=huge(1.0_dp))),'overflowed bounded width')
    call check(ieee_is_nan(genpareto_censored_loglikelihood(empty,-1.0_dp, &
        loc=huge(1.0_dp))),'collapsed bounded support')
    call close(genpareto_censored_loglikelihood(empty,1.0_dp,loc=huge(1.0_dp)), &
        0.0_dp,0.0_dp,'unbounded law permits unrepresentable scale point')
    t = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(left=[t])
    call close(genpareto_censored_loglikelihood(one,-t,scale=t), &
        log(1.0_dp-exp(-1.0_dp)),3.0e-13_dp,'tiny negative shape and scale')
    call close(genpareto_censored_loglikelihood(one,t,scale=t), &
        log(1.0_dp-exp(-1.0_dp)),3.0e-13_dp,'tiny positive shape and scale')
    do family = 1, 3
        one = make_censored_data(interval=reshape([0.0_dp,t],[1,2]))
        other = make_censored_data(uncensored=[0.0_dp])
        expected = likelihood(other,family,0.0_dp,scale=1.0e300_dp)+log(t)
        call close(likelihood(one,family,0.0_dp,scale=1.0e300_dp),expected, &
            3.0e-13_dp,'subnormal interval and huge scale')
    end do
    one = make_censored_data(uncensored=[1.0_dp])
    expected = log(2.0_dp)-log(acos(-1.0_dp))+log(t)
    call close(foldcauchy_censored_loglikelihood(one,0.0_dp,scale=t),expected, &
        3.0e-13_dp,'finite Cauchy density beyond standardized overflow')
    one = make_censored_data(uncensored=[1.0e300_dp])
    call check(foldnorm_censored_loglikelihood(one,0.0_dp) == -inf, &
        'overflowing normal squared distance')
    one = make_censored_data(left=[t])
    expected = log(2.0_dp)-log(acos(-1.0_dp))+log(t)-2.0_dp*log(1.0e300_dp)
    call close(foldcauchy_censored_loglikelihood(one,1.0e300_dp),expected, &
        3.0e-13_dp,'tiny folded CDF and huge offset')
    print '(a,i0)', 'Censored folded/Pareto checks: ',checks
    if (failures /= 0) error stop 'censored folded/Pareto failures'
contains
    function likelihood(sample, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 folded normal, 2 folded Cauchy, 3 generalized Pareto
        real(dp), intent(in) :: shape !! nonnegative folding offset or finite Pareto shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = foldnorm_censored_loglikelihood(sample,shape,loc,scale)
        case (2)
            value = foldcauchy_censored_loglikelihood(sample,shape,loc,scale)
        case (3)
            value = genpareto_censored_loglikelihood(sample,shape,loc,scale)
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function likelihood
    subroutine fitting(sample, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! family selector as in likelihood
        real(dp), intent(in) :: lower(3) !! shape, location, scale lower bounds
        real(dp), intent(in) :: upper(3) !! upper bounds; equality fixes a parameter
        type(fit_result), intent(out) :: result !! parameters and fit diagnostics
        select case (family)
        case (1)
            call foldnorm_fit(sample,lower,upper,result,max_iter=10000)
        case (2)
            call foldcauchy_fit(sample,lower,upper,result,max_iter=10000)
        case (3)
            call genpareto_fit(sample,lower,upper,result,max_iter=10000)
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
end program test_censored_fold_pareto
