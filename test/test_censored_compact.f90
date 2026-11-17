! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_compact
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        cosine_censored_loglikelihood, cosine_censored_nnlf, cosine_fit, &
        cosine_logpdf, cosine_logcdf, cosine_logsf, cosine_ppf, &
        anglit_censored_loglikelihood, anglit_censored_nnlf, anglit_fit, &
        anglit_logpdf, anglit_logcdf, anglit_logsf, anglit_ppf, &
        semicircular_censored_loglikelihood, semicircular_censored_nnlf, semicircular_fit, &
        semicircular_logpdf, semicircular_logcdf, semicircular_logsf, semicircular_ppf, &
        fit_status_invalid_input, fit_status_no_finite_objective
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: one, other, empty, data
    type(fit_result) :: fit, old
    real(dp) :: actual, expected, x, p, previous, inf, nan, h, lo(2), hi(2)
    real(dp), parameter :: endpoint_references(3) = &
        [-90.85432149917743_dp,-63.91634327289527_dp,-50.657359060721994_dp]
    real(dp), parameter :: scipy_scales(3) = &
        [0.7_dp,1.0646672682233491_dp,0.70780375432180109_dp]
    integer :: family, i, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[-0.2_dp,0.1_dp,0.35_dp],left=[-0.4_dp], &
        right=[0.45_dp],interval=reshape([-0.3_dp,0.2_dp,-0.1_dp,0.25_dp],[2,2]))
    do family = 1, 3
        h = 1.0_dp
        if (family == 1) h = acos(-1.0_dp)
        if (family == 2) h = acos(-1.0_dp)/4.0_dp
        previous = -inf
        do i = 1, 39
            x = h*real(i-20,dp)/20.0_dp
            one = make_censored_data(uncensored=[x])
            select case (family)
            case (1)
                expected = cosine_logpdf(x)
            case (2)
                expected = anglit_logpdf(x)
            case (3)
                expected = semicircular_logpdf(x)
            end select
            call close(likelihood(one,family),expected,2.0e-12_dp,'scalar density')
            one = make_censored_data(left=[x])
            select case (family)
            case (1)
                expected = cosine_logcdf(x)
            case (2)
                expected = anglit_logcdf(x)
            case (3)
                expected = semicircular_logcdf(x)
            end select
            actual = likelihood(one,family)
            call close(actual,expected,2.0e-12_dp,'scalar CDF')
            call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
            previous = actual
            other = make_censored_data(right=[-x])
            call close(actual,likelihood(other,family),2.0e-14_dp,'reflection tails')
            one = make_censored_data(right=[x])
            select case (family)
            case (1)
                expected = cosine_logsf(x)
            case (2)
                expected = anglit_logsf(x)
            case (3)
                expected = semicircular_logsf(x)
            end select
            call close(likelihood(one,family),expected,2.0e-12_dp,'scalar SF')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            other = make_censored_data(uncensored=[x])
            expected = likelihood(other,family)+log(nearest(x,1.0_dp)-x)
            call close(likelihood(one,family),expected,3.0e-13_dp,'adjacent interval')
        end do
        do i = 1, 19
            p = real(i,dp)/20.0_dp
            select case (family)
            case (1)
                x = cosine_ppf(p)
            case (2)
                x = anglit_ppf(p)
            case (3)
                x = semicircular_ppf(p)
            end select
            one = make_censored_data(left=[x])
            call close(likelihood(one,family),log(p),1.0e-12_dp,'CDF/PPF round trip')
        end do
        call close(likelihood(empty,family),0.0_dp,0.0_dp,'empty sample')
        call check(ieee_is_nan(likelihood(empty,family,scale=0.0_dp)),'invalid scale')
        call check(ieee_is_nan(likelihood(empty,family,scale=-1.0_dp)),'negative scale')
        call check(ieee_is_nan(likelihood(empty,family,scale=inf)),'infinite scale')
        call check(ieee_is_nan(likelihood(empty,family,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,loc=huge(1.0_dp))),'bad support')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(one,family)),'NaN observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(one,family)),'infinite exact observation')
        one = make_censored_data(interval=reshape([0.0_dp,0.0_dp],[1,2]))
        other = make_censored_data(uncensored=[0.0_dp])
        call close(likelihood(one,family),likelihood(other,family),0.0_dp,'equal endpoints')
        one = make_censored_data(uncensored=[h])
        call check(likelihood(one,family) == -inf,'zero endpoint density')
        one = make_censored_data(left=[-h])
        call check(likelihood(one,family) == -inf,'lower endpoint CDF')
        one = make_censored_data(right=[h])
        call check(likelihood(one,family) == -inf,'upper endpoint SF')
        one = make_censored_data(interval=reshape([-2.0_dp*h,2.0_dp*h],[1,2]))
        call close(likelihood(one,family),0.0_dp,0.0_dp,'whole support interval')
        one = make_censored_data(interval=reshape([2.0_dp*h,3.0_dp*h],[1,2]))
        call check(likelihood(one,family) == -inf,'outside support interval')
        one = make_censored_data(interval=reshape([-inf,0.0_dp],[1,2]))
        call close(likelihood(one,family),-log(2.0_dp),2.0e-14_dp,'infinite interval edge')
        one = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
        other = make_censored_data(uncensored=[0.0_dp])
        expected = likelihood(other,family,scale=1.0e300_dp)+log(nearest(0.0_dp,1.0_dp))
        call close(likelihood(one,family,scale=1.0e300_dp),expected,3.0e-13_dp, &
            'subnormal interval with huge support')
        x = h*(1.0_dp-1.0e-12_dp)
        one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
        call close(likelihood(one,family),endpoint_references(family),3.0e-13_dp, &
            'endpoint adjacent interval mpmath reference')
        one = make_censored_data(uncensored=[0.0_dp])
        call check(ieee_is_finite(likelihood(one,family,scale=nearest(0.0_dp,1.0_dp))), &
            'smallest positive scale density')
        one = make_censored_data(uncensored=[h])
        lo = [0.0_dp,1.0_dp]
        call fitting(one,family,lo,lo,fit)
        call check(fit%status == fit_status_no_finite_objective,'endpoint fixed fit')
        call fitting(empty,family,lo,lo,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit')
        lo = [-0.2_dp,0.7_dp]
        hi = [0.2_dp,2.0_dp]
        call fitting(data,family,lo,hi,fit)
        call check(fit%success,'two parameter fit converges')
        call check(all(fit%params >= lo) .and. all(fit%params <= hi),'bounded fit')
        call close(fit%nllf,-likelihood(data,family,fit%params(1),fit%params(2)), &
            0.0_dp,'fit objective')
        call check(fit%nllf <= -likelihood(data,family,0.0_dp,1.0_dp),'fit improvement')
        lo = [0.0_dp,0.7_dp]
        hi = [0.0_dp,2.0_dp]
        call fitting(data,family,lo,hi,fit)
        call check(fit%success,'fixed center fit converges')
        call close(fit%params(2),scipy_scales(family),2.0e-6_dp,'SciPy scale benchmark')
        one = make_censored_data(uncensored=[-0.2_dp,0.1_dp,0.35_dp])
        lo = [0.0_dp,1.0_dp]
        call fitting(one,family,lo,lo,fit)
        select case (family)
        case (1)
            call cosine_fit(one%uncensored,lo,lo,old)
            expected = cosine_censored_nnlf(data)
        case (2)
            call anglit_fit(one%uncensored,lo,lo,old)
            expected = anglit_censored_nnlf(data)
        case (3)
            call semicircular_fit(one%uncensored,lo,lo,old)
            expected = semicircular_censored_nnlf(data)
        end select
        call close(fit%nllf,old%nllf,2.0e-13_dp,'array fit overload retained')
        call close(expected,-likelihood(data,family),0.0_dp,'NNLF')
    end do
    call check(ieee_is_nan(cosine_censored_loglikelihood(empty,scale=huge(1.0_dp))), &
        'overflowing half-width rejected')
    print '(a,i0)', 'Censored compact-support checks: ',checks
    if (failures /= 0) error stop 'censored compact-support failures'
contains
    function likelihood(sample, family, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! 1 cosine, 2 anglit, 3 semicircular
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = cosine_censored_loglikelihood(sample,loc,scale)
        case (2)
            value = anglit_censored_loglikelihood(sample,loc,scale)
        case (3)
            value = semicircular_censored_loglikelihood(sample,loc,scale)
        case default
            value = quiet_nan(0.0_dp)
        end select
    end function likelihood
    subroutine fitting(sample, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: sample !! exact and censored observations
        integer, intent(in) :: family !! family selector as in likelihood
        real(dp), intent(in) :: lower(2) !! lower bounds for location and scale
        real(dp), intent(in) :: upper(2) !! upper bounds; equality fixes a parameter
        type(fit_result), intent(out) :: result !! fit parameters and diagnostics
        select case (family)
        case (1)
            call cosine_fit(sample,lower,upper,result,max_iter=10000)
        case (2)
            call anglit_fit(sample,lower,upper,result,max_iter=10000)
        case (3)
            call semicircular_fit(sample,lower,upper,result,max_iter=10000)
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
end program test_censored_compact
