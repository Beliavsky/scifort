! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_bounded_shape
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        truncexpon_censored_loglikelihood, truncexpon_censored_nnlf, truncexpon_fit, &
        truncexpon_logpdf, truncexpon_logcdf, truncexpon_logsf, truncexpon_ppf, &
        johnsonsb_censored_loglikelihood, johnsonsb_censored_nnlf, johnsonsb_fit, &
        johnsonsb_logpdf, johnsonsb_logcdf, johnsonsb_logsf, johnsonsb_ppf, &
        kappa3_censored_loglikelihood, kappa3_censored_nnlf, kappa3_fit, &
        kappa3_logpdf, kappa3_logcdf, kappa3_logsf, kappa3_ppf, &
        lomax_censored_loglikelihood, fit_status_invalid_input, fit_status_no_finite_objective
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, other, empty
    type(fit_result) :: fit, old
    real(dp) :: x, p, c, actual, expected, previous, inf, nan, lower(3), upper(3)
    integer :: family, i, j, checks, failures
    real(dp), parameter :: scipy_shapes(3) = &
        [0.66503980445823563_dp,1.120500197116846_dp,3.0_dp]
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
            if (family == 2) c = real(j-2,dp)
            previous = -inf
            do i = 1, 19
                x = real(i,dp)/20.0_dp
                if (family == 1) x = c*x
                if (family == 3) x = exp(real(i-10,dp)/2.0_dp)
                one = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = truncexpon_logpdf(x,c)
                case (2)
                    expected = johnsonsb_logpdf(x,c,1.2_dp)
                case (3)
                    expected = kappa3_logpdf(x,c)
                end select
                call close(likelihood(one,family,c),expected,2.0e-12_dp,'scalar density')
                one = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = truncexpon_logcdf(x,c)
                case (2)
                    expected = johnsonsb_logcdf(x,c,1.2_dp)
                case (3)
                    expected = kappa3_logcdf(x,c)
                end select
                actual = likelihood(one,family,c)
                call close(actual,expected,2.0e-12_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = truncexpon_logsf(x,c)
                case (2)
                    expected = johnsonsb_logsf(x,c,1.2_dp)
                case (3)
                    expected = kappa3_logsf(x,c)
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
                    x = truncexpon_ppf(p,c)
                case (2)
                    x = johnsonsb_ppf(p,c,1.2_dp)
                case (3)
                    x = kappa3_ppf(p,c)
                end select
                one = make_censored_data(left=[x])
                call close(likelihood(one,family,c),log(p),3.0e-12_dp,'CDF/PPF round trip')
            end do
        end do
        call close(likelihood(empty,family,1.0_dp),0.0_dp,0.0_dp,'empty sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        if (family /= 2) then
            call check(ieee_is_nan(likelihood(empty,family,0.0_dp)),'zero shape')
            call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'negative shape')
        end if
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
        lower = [1.7_dp,0.0_dp,1.3_dp]
        call fitting(data,family,lower,lower,fit)
        call check(fit%success,'fixed parameter fit')
        call close(fit%nllf,-likelihood(data,family,1.7_dp,scale=1.3_dp),0.0_dp, &
            'fixed fit objective')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fitting(one,family,lower,lower,fit)
        select case (family)
        case (1)
            call truncexpon_fit(one%uncensored,lower,lower,old)
            expected = truncexpon_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (2)
            call johnsonsb_fit(one%uncensored,[lower(1),1.2_dp,lower(2:3)],[lower(1),1.2_dp,lower(2:3)],old)
            expected = johnsonsb_censored_nnlf(data,1.7_dp,1.2_dp,scale=1.3_dp)
        case (3)
            call kappa3_fit(one%uncensored,lower,lower,old)
            expected = kappa3_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        end select
        call close(fit%nllf,old%nllf,2.0e-13_dp,'original array fit retained')
        call close(expected,-likelihood(data,family,1.7_dp,scale=1.3_dp),0.0_dp,'NNLF')
        lower = [0.2_dp,-0.1_dp,0.9_dp]
        upper = [3.0_dp,0.1_dp,2.0_dp]
        if (family == 2) lower(1) = -3.0_dp
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'free parameter fit')
        if (family == 2) then
            call close(fit%nllf,-likelihood(data,family,fit%params(1), &
                fit%params(3),fit%params(4)),0.0_dp,'free fit objective')
        else
            call close(fit%nllf,-likelihood(data,family,fit%params(1), &
                fit%params(2),fit%params(3)),0.0_dp,'free fit objective')
        end if
        call fitting(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit')
        lower = [0.1_dp,0.0_dp,1.3_dp]
        upper = [3.0_dp,0.0_dp,1.3_dp]
        if (family == 1) lower(1) = 0.6_dp
        if (family == 2) lower(1) = -3.0_dp
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'fixed support shape fit')
        call close(fit%params(1),scipy_shapes(family),2.0e-6_dp,'SciPy shape benchmark')
    end do
    do i = 1, 10
        x = exp(real(i-5,dp))
        one = make_censored_data(uncensored=[x],left=[0.1_dp],right=[x+1.0_dp], &
            interval=reshape([x,x+0.1_dp],[1,2]))
        call close(kappa3_censored_loglikelihood(one,1.0_dp), &
            lomax_censored_loglikelihood(one,1.0_dp),2.0e-13_dp,'unit shape Lomax identity')
    end do
    one = make_censored_data(uncensored=[1.0_dp])
    call close(kappa3_censored_loglikelihood(one,1.0e-300_dp),log(1.0e-300_dp)-1.0_dp, &
        3.0e-13_dp,'tiny kappa shape density')
    one = make_censored_data(left=[1.0_dp])
    call close(kappa3_censored_loglikelihood(one,1.0e-300_dp),-1.0_dp, &
        3.0e-13_dp,'tiny kappa shape CDF')
    one = make_censored_data(uncensored=[0.0_dp])
    call check(kappa3_censored_loglikelihood(one,nearest(0.0_dp,1.0_dp)) == inf, &
        'overflowing kappa endpoint log density')
    call check(johnsonsb_censored_loglikelihood(one,0.0_dp,1.0_dp) == -inf, &
        'Johnson lower endpoint density')
    call close(truncexpon_censored_loglikelihood(one,1.0_dp), &
        -log(1.0_dp-exp(-1.0_dp)),2.0e-14_dp,'truncated exponential endpoint density')
    one = make_censored_data(left=[1.0e-300_dp])
    actual = johnsonsb_censored_loglikelihood(one,0.0_dp,1.0e300_dp,loc=-0.5_dp)
    call close(actual,-3.1671743377489271e-5_dp,2.0e-14_dp,'amplified Johnson midpoint')
    one = make_censored_data(right=[1.0e-300_dp])
    actual = kappa3_censored_loglikelihood(one,1.0e300_dp,loc=-1.0_dp)
    expected = log((log(1.0e300_dp)-1.0_dp)/1.0e300_dp)
    call close(actual,expected,3.0e-13_dp,'amplified kappa scale-point displacement')
    call check(ieee_is_nan(truncexpon_censored_loglikelihood(empty,1.0e300_dp, &
        scale=1.0e300_dp)),'overflowing truncated support')
    call check(ieee_is_nan(johnsonsb_censored_loglikelihood(empty,0.0_dp,0.0_dp)), &
        'zero second Johnson shape')
    one = make_censored_data(left=[1.5e308_dp])
    call check(ieee_is_finite(kappa3_censored_loglikelihood(one,2.0_dp, &
        loc=1.0e308_dp,scale=1.0e308_dp)),'overflowing kappa scale point allowed')
    call check(ieee_is_nan(truncexpon_censored_loglikelihood(empty, &
        nearest(0.0_dp,1.0_dp),scale=nearest(0.0_dp,1.0_dp))),'underflowing width')
    call check(ieee_is_nan(johnsonsb_censored_loglikelihood(empty,0.0_dp,1.0_dp, &
        loc=huge(1.0_dp),scale=1.0_dp)),'collapsed Johnson support')
    one = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    other = make_censored_data(uncensored=[0.0_dp])
    call close(truncexpon_censored_loglikelihood(one,1.0_dp,scale=1.0e300_dp), &
        truncexpon_censored_loglikelihood(other,1.0_dp,scale=1.0e300_dp) &
            +log(nearest(0.0_dp,1.0_dp)),3.0e-13_dp,'subnormal truncated interval')
    one = make_censored_data(uncensored=[0.25_dp])
    actual = truncexpon_censored_loglikelihood(one,1.0e-300_dp,scale=1.0e300_dp)
    call close(actual/2.5e-301_dp,1.0_dp,3.0e-13_dp,'nearly uniform tiny log density')
    x = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(uncensored=[10.0_dp*x])
    expected = -10.0_dp/21.0_dp-log(21.0_dp*x)-log(1.0_dp-exp(-1.0_dp))
    call close(truncexpon_censored_loglikelihood(one,1.0_dp,scale=21.0_dp*x), &
        expected,3.0e-13_dp,'odd subnormal truncated width')
    p = 100.0_dp*log(10.0_dp/11.0_dp)
    expected = -0.5_dp*p*p-0.5_dp*log(2.0_dp*acos(-1.0_dp)) &
        +log(100.0_dp)-log(21.0_dp*x)-log(10.0_dp/21.0_dp)-log(11.0_dp/21.0_dp)
    call close(johnsonsb_censored_loglikelihood(one,0.0_dp,100.0_dp,scale=21.0_dp*x), &
        expected,3.0e-13_dp,'odd subnormal Johnson midpoint')
    call johnsonsb_fit(data,[-3.0_dp,0.5_dp,-0.1_dp,0.9_dp], &
        [3.0_dp,2.0_dp,0.1_dp,2.0_dp],fit,max_iter=10000)
    call check(fit%success,'four free Johnson parameters')
    call close(fit%nllf,johnsonsb_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4)),0.0_dp,'four parameter Johnson objective')
    print '(a,i0)', 'Censored bounded/shape checks: ',checks
    if (failures /= 0) error stop 'censored bounded/shape failures'
contains
    function likelihood(sample, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 truncexpon, 2 johnsonsb, 3 kappa3
        real(dp), intent(in) :: shape !! first shape; positive except Johnson SB
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = truncexpon_censored_loglikelihood(sample,shape,loc,scale)
        case (2)
            value = johnsonsb_censored_loglikelihood(sample,shape,1.2_dp,loc,scale)
        case (3)
            value = kappa3_censored_loglikelihood(sample,shape,loc,scale)
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function likelihood
    subroutine fitting(sample, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! family selector as in likelihood
        real(dp), intent(in) :: lower(3) !! lower shape/location/scale bounds
        real(dp), intent(in) :: upper(3) !! upper bounds; equality fixes a parameter
        type(fit_result), intent(out) :: result !! parameters and fit diagnostics
        select case (family)
        case (1)
            call truncexpon_fit(sample,lower,upper,result,max_iter=10000)
        case (2)
            call johnsonsb_fit(sample,[lower(1),1.2_dp,lower(2:3)], &
                [upper(1),1.2_dp,upper(2:3)],result,max_iter=10000)
        case (3)
            call kappa3_fit(sample,lower,upper,result,max_iter=10000)
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
end program test_censored_bounded_shape
