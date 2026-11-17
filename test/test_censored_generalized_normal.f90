! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_generalized_normal
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        gennorm_censored_loglikelihood, gennorm_censored_nnlf, gennorm_fit, &
        gennorm_logpdf, gennorm_logcdf, gennorm_logsf, gennorm_ppf, &
        halfgennorm_censored_loglikelihood, halfgennorm_censored_nnlf, halfgennorm_fit, &
        halfgennorm_logpdf, halfgennorm_logcdf, halfgennorm_logsf, halfgennorm_ppf, &
        rdist_censored_loglikelihood, rdist_censored_nnlf, rdist_fit, &
        rdist_logpdf, rdist_logcdf, rdist_logsf, rdist_ppf, &
        laplace_censored_loglikelihood, exponential_censored_loglikelihood, &
        normal_censored_loglikelihood, halfnorm_censored_loglikelihood, &
        uniform_censored_loglikelihood, semicircular_censored_loglikelihood, &
        arcsine_censored_loglikelihood, fit_status_invalid_input
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, other, empty
    type(fit_result) :: fit, old
    real(dp) :: x, p, c, actual, expected, previous, inf, nan, lower(3), upper(3), t
    integer :: family, i, j, checks, failures
    real(dp), parameter :: scipy_shape(3) = &
        [4.0970508107852872_dp,3.8132528818587343_dp,5.0_dp]
    real(dp), parameter :: scipy_objective(3) = &
        [11.851573580953866_dp,9.5469968381416006_dp,11.137246801642787_dp]
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        do j = 1, 3
            c = real(j,dp)/2.0_dp
            if (family == 3) c = real(j,dp)
            previous = -inf
            do i = 1, 19
                x = real(i,dp)/20.0_dp
                select case (family)
                case (1)
                    expected = gennorm_logpdf(x,c)
                case (2)
                    expected = halfgennorm_logpdf(x,c)
                case (3)
                    expected = rdist_logpdf(x,c)
                end select
                one = make_censored_data(uncensored=[x])
                call close(likelihood(one,family,c),expected,2.0e-12_dp,'scalar density')
                one = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = gennorm_logcdf(x,c)
                case (2)
                    expected = halfgennorm_logcdf(x,c)
                case (3)
                    expected = rdist_logcdf(x,c)
                end select
                actual = likelihood(one,family,c)
                call close(actual,expected,2.0e-12_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = gennorm_logsf(x,c)
                case (2)
                    expected = halfgennorm_logsf(x,c)
                case (3)
                    expected = rdist_logsf(x,c)
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
                    x = gennorm_ppf(p,c)
                case (2)
                    x = halfgennorm_ppf(p,c)
                case (3)
                    x = rdist_ppf(p,c)
                end select
                one = make_censored_data(left=[x])
                call close(likelihood(one,family,c),log(p),3.0e-12_dp,'CDF/PPF round trip')
            end do
        end do
        call close(likelihood(empty,family,1.0_dp),0.0_dp,0.0_dp,'empty sample')
        call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)), &
            'negative folding offset')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,loc=nan)),'NaN location')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(one,family,1.0_dp)),'NaN observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(one,family,1.0_dp)),'infinite observation')
        one = make_censored_data(uncensored=[-2.0_dp])
        if (family /= 1) call check(likelihood(one,family,1.0_dp) == -inf, &
            'outside support density')
        one = make_censored_data(left=[0.0_dp])
        if (family == 2) call check(likelihood(one,family,1.0_dp) == -inf, &
            'lower endpoint CDF')
        one = make_censored_data(right=[0.0_dp])
        expected = 0.0_dp
        if (family /= 2) expected = -log(2.0_dp)
        call close(likelihood(one,family,1.0_dp),expected,3.0e-13_dp,'central SF')
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
            call gennorm_fit(one%uncensored,lower,lower,old)
            actual = gennorm_censored_nnlf(data,0.5_dp,scale=1.3_dp)
        case (2)
            call halfgennorm_fit(one%uncensored,lower,lower,old)
            actual = halfgennorm_censored_nnlf(data,0.5_dp,scale=1.3_dp)
        case (3)
            call rdist_fit(one%uncensored,lower,lower,old)
            actual = rdist_censored_nnlf(data,0.5_dp,scale=1.3_dp)
        end select
        call check(old%success,'original array fit')
        call close(fit%nllf,old%nllf,2.0e-12_dp,'array/censored fit agreement')
        call close(actual,-likelihood(data,family,0.5_dp,scale=1.3_dp),0.0_dp,'NNLF')
        call fitting(empty,family,lower,lower,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        lower(1) = 0.3_dp
        upper = [5.0_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'free shape fit')
        call close(fit%nllf,-likelihood(data,family,fit%params(1),scale=1.3_dp), &
            0.0_dp,'free fit objective')
        call close(fit%params(1),scipy_shape(family),3.0e-5_dp,'SciPy shape optimum')
        call close(fit%nllf,scipy_objective(family),2.0e-8_dp,'SciPy optimal objective')
        lower = [0.3_dp,-0.1_dp,0.9_dp]
        upper = [5.0_dp,0.1_dp,2.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'three free parameters')
        call close(fit%nllf,-likelihood(data,family,fit%params(1), &
            fit%params(2),fit%params(3)),0.0_dp,'three parameter objective')
    end do
    call close(gennorm_censored_loglikelihood(data,1.0_dp), &
        laplace_censored_loglikelihood(data),3.0e-13_dp,'Laplace identity')
    call close(halfgennorm_censored_loglikelihood(data,1.0_dp), &
        exponential_censored_loglikelihood(data),3.0e-13_dp,'exponential identity')
    call close(gennorm_censored_loglikelihood(data,2.0_dp), &
        normal_censored_loglikelihood(data,scale=1.0_dp/sqrt(2.0_dp)), &
        3.0e-13_dp,'normal identity')
    call close(halfgennorm_censored_loglikelihood(data,2.0_dp), &
        halfnorm_censored_loglikelihood(data,scale=1.0_dp/sqrt(2.0_dp)), &
        3.0e-13_dp,'half-normal identity')
    call close(rdist_censored_loglikelihood(data,1.0_dp), &
        arcsine_censored_loglikelihood(data,loc=-1.0_dp,scale=2.0_dp), &
        3.0e-13_dp,'arcsine identity')
    call close(rdist_censored_loglikelihood(data,2.0_dp), &
        uniform_censored_loglikelihood(data,loc=-1.0_dp,scale=2.0_dp), &
        3.0e-13_dp,'uniform identity')
    call close(rdist_censored_loglikelihood(data,3.0_dp), &
        semicircular_censored_loglikelihood(data),3.0e-13_dp,'semicircular identity')
    t = nearest(0.0_dp,1.0_dp)
    call check(ieee_is_nan(gennorm_censored_loglikelihood(empty,t)), &
        'unrepresentable reciprocal power')
    call check(ieee_is_nan(halfgennorm_censored_loglikelihood(empty,t)), &
        'unrepresentable half-law reciprocal power')
    call close(gennorm_censored_loglikelihood(empty,1.0e-308_dp),0.0_dp,0.0_dp, &
        'valid empty law with overflowing gamma log normalization')
    call close(halfgennorm_censored_loglikelihood(empty,1.0e-308_dp),0.0_dp,0.0_dp, &
        'valid empty half-law with overflowing gamma log normalization')
    call check(ieee_is_nan(rdist_censored_loglikelihood(empty,t)), &
        'unrepresentable half shape')
    call check(ieee_is_nan(rdist_censored_loglikelihood(empty,1.0e7_dp)), &
        'R shape beyond validated range')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(rdist_censored_loglikelihood(one,2.0_dp*t),log(t),3.0e-13_dp, &
        'subnormal symmetric beta normalization')
    one = make_censored_data(left=[0.0_dp])
    call close(rdist_censored_loglikelihood(one,1.0e6_dp),-log(2.0_dp), &
        1.0e-14_dp,'symmetric CDF at maximum R shape')
    call check(ieee_is_nan(rdist_censored_loglikelihood(empty,2.0_dp, &
        scale=huge(1.0_dp))),'overflowed R width')
    call check(ieee_is_nan(rdist_censored_loglikelihood(empty,2.0_dp, &
        loc=huge(1.0_dp))),'collapsed R support')
    one = make_censored_data(uncensored=[1.0_dp])
    call check(rdist_censored_loglikelihood(one,1.0_dp) == inf,'singular R endpoint')
    call close(rdist_censored_loglikelihood(one,2.0_dp),-log(2.0_dp),1.0e-14_dp, &
        'uniform R endpoint')
    call check(rdist_censored_loglikelihood(one,3.0_dp) == -inf,'zero R endpoint density')
    one = make_censored_data(left=[1.0_dp])
    call close(rdist_censored_loglikelihood(one,2.0_dp),0.0_dp,0.0_dp,'upper R CDF')
    one = make_censored_data(right=[1.0_dp])
    call check(rdist_censored_loglikelihood(one,2.0_dp) == -inf,'upper R SF')
    one = make_censored_data(uncensored=[0.0_dp])
    actual = halfgennorm_censored_loglikelihood(one,1.0e300_dp)
    expected = 5.7721566490153286e-301_dp
    call close(actual/expected,1.0_dp,3.0e-13_dp,'tiny endpoint log density')
    call check(halfgennorm_censored_loglikelihood(one,1.0e-308_dp) == -inf, &
        'overflowing gamma log normalization')
    one = make_censored_data(left=[0.5_dp])
    call check(halfgennorm_censored_loglikelihood(one,1.0e-308_dp) == -inf, &
        'unrepresentably small gamma lower log tail')
    one = make_censored_data(interval=reshape([-1.0e-300_dp,1.0e-300_dp],[1,2]))
    expected = -0.5_dp*log(acos(-1.0_dp))+log(2.0e-300_dp)
    call close(gennorm_censored_loglikelihood(one,2.0_dp),expected,3.0e-13_dp, &
        'tiny interval crossing center')
    one = make_censored_data(interval=reshape([0.0_dp,t],[1,2]))
    expected = log(t)-log(1.0e300_dp)
    call close(halfgennorm_censored_loglikelihood(one,1.0_dp,scale=1.0e300_dp), &
        expected,3.0e-13_dp,'subnormal interval and huge scale')
    one = make_censored_data(left=[1.0e-300_dp])
    actual = halfgennorm_censored_loglikelihood(one,1.0e300_dp,loc=-1.0_dp)
    call check(actual < 0.0_dp .and. ieee_is_finite(actual),'amplified scale displacement')
    one = make_censored_data(interval=reshape([-0.5_dp,-0.2_dp],[1,2]))
    other = make_censored_data(interval=reshape([0.2_dp,0.5_dp],[1,2]))
    call close(gennorm_censored_loglikelihood(one,0.7_dp), &
        gennorm_censored_loglikelihood(other,0.7_dp),0.0_dp,'interval reflection symmetry')
    one = make_censored_data(interval=reshape( &
        [0.99999188830481667_dp,0.99999309226857946_dp],[1,2]))
    expected = -13.630472120271284593_dp
    call close(halfgennorm_censored_loglikelihood(one,1.0e6_dp),expected, &
        3.0e-13_dp,'steep powered-coordinate interval')
    call close(gennorm_censored_loglikelihood(one,1.0e6_dp),expected-log(2.0_dp), &
        3.0e-13_dp,'symmetric steep powered-coordinate interval')
    one = make_censored_data(left=[0.1_dp])
    call close(halfgennorm_censored_loglikelihood(one,huge(1.0_dp)),log(0.1_dp), &
        3.0e-13_dp,'finite CDF with overflowing negative log power')
    call close(gennorm_censored_loglikelihood(one,huge(1.0_dp)),log(0.55_dp), &
        3.0e-13_dp,'symmetric CDF with overflowing negative log power')
    one = make_censored_data(right=[0.1_dp])
    call close(halfgennorm_censored_loglikelihood(one,huge(1.0_dp)),log(0.9_dp), &
        3.0e-13_dp,'finite SF with overflowing negative log power')
    one = make_censored_data(left=[1.0_dp])
    actual = halfgennorm_censored_loglikelihood(one,huge(1.0_dp))
    expected = -1.2203636434998277e-309_dp
    call close(actual/expected,1.0_dp,3.0e-12_dp,'subnormal log CDF at largest power')
    print '(a,i0)', 'Censored generalized normal/R checks: ',checks
    if (failures /= 0) error stop 'censored generalized normal/R failures'
contains
    function likelihood(sample, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 generalized normal, 2 half generalized normal, 3 R distribution
        real(dp), intent(in) :: shape !! positive family-specific shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = gennorm_censored_loglikelihood(sample,shape,loc,scale)
        case (2)
            value = halfgennorm_censored_loglikelihood(sample,shape,loc,scale)
        case (3)
            value = rdist_censored_loglikelihood(sample,shape,loc,scale)
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
            call gennorm_fit(sample,lower,upper,result,max_iter=10000)
        case (2)
            call halfgennorm_fit(sample,lower,upper,result,max_iter=10000)
        case (3)
            call rdist_fit(sample,lower,upper,result,max_iter=10000)
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
end program test_censored_generalized_normal
