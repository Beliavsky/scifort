! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_normal_transform
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        moyal_censored_loglikelihood, moyal_censored_nnlf, moyal_fit, &
        moyal_logpdf, moyal_logcdf, moyal_logsf, moyal_ppf, &
        levy_censored_loglikelihood, levy_censored_nnlf, levy_fit, &
        levy_logpdf, levy_logcdf, levy_logsf, levy_ppf, &
        levy_l_censored_loglikelihood, levy_l_censored_nnlf, levy_l_fit, &
        levy_l_logpdf, levy_l_logcdf, levy_l_logsf, levy_l_ppf, &
        halfnorm_censored_loglikelihood, fit_status_invalid_input, &
        fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, other, reflected
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, actual, expected, x, p, previous, reference(3), lower(2), upper(2)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    reference = [-15.130277471440076_dp,-15.99466380610435_dp,-15.99466380610435_dp]
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    reflected = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp],right=[-0.3_dp], &
        left=[-3.5_dp],interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
    do family = 1, 3
        if (family == 3) data = reflected
        actual = likelihood(data,family,scale=1.3_dp)
        call close(actual,reference(family),2.0e-13_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = moyal_censored_nnlf(data,scale=1.3_dp)
        case (2)
            expected = levy_censored_nnlf(data,scale=1.3_dp)
        case (3)
            expected = levy_l_censored_nnlf(data,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            previous = -inf
            do i = 1, 12
                x = real(i,dp)/3.0_dp
                if (family == 1) x = x-3.0_dp
                if (family == 3) x = x-5.0_dp
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = moyal_logpdf(x,scale=real(j,dp))
                case (2)
                    expected = levy_logpdf(x,scale=real(j,dp))
                case (3)
                    expected = levy_l_logpdf(x,scale=real(j,dp))
                end select
                call close(likelihood(single,family,scale=real(j,dp)),expected,6.0e-14_dp, &
                    'scalar density')
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = moyal_logcdf(x,scale=real(j,dp))
                case (2)
                    expected = levy_logcdf(x,scale=real(j,dp))
                case (3)
                    expected = levy_l_logcdf(x,scale=real(j,dp))
                end select
                actual = likelihood(single,family,scale=real(j,dp))
                call close(actual,expected,6.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity/range')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = moyal_logsf(x,scale=real(j,dp))
                case (2)
                    expected = levy_logsf(x,scale=real(j,dp))
                case (3)
                    expected = levy_l_logsf(x,scale=real(j,dp))
                end select
                call close(likelihood(single,family,scale=real(j,dp)),expected,6.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = moyal_ppf(p,scale=real(j,dp))
                case (2)
                    x = levy_ppf(p,scale=real(j,dp))
                case (3)
                    x = levy_l_ppf(p,scale=real(j,dp))
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,scale=real(j,dp))),p,3.0e-14_dp, &
                    'CDF/PPF round trip')
            end do
        end do
        single = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
        if (family == 3) single = make_censored_data(uncensored=[-0.5_dp,-1.0_dp,-2.0_dp])
        lower = [0.0_dp,1.3_dp]
        upper = lower
        call fitting(single,family,lower,upper,fit)
        select case (family)
        case (1)
            call moyal_fit(single%uncensored,lower,upper,old)
        case (2)
            call levy_fit(single%uncensored,lower,upper,old)
        case (3)
            call levy_l_fit(single%uncensored,lower,upper,old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,2.0e-13_dp,'array fit compatibility')
        call close(fit%nllf,-likelihood(single,family,scale=1.3_dp),0.0_dp,'fit objective')
        call close(likelihood(empty,family),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(likelihood(empty,family,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,scale=nan)),'NaN scale')
        call check(ieee_is_nan(likelihood(empty,family,scale=inf)),'infinite scale')
        call check(ieee_is_nan(likelihood(empty,family,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,scale=-1.0_dp)),'negative scale')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family)),'NaN exact value')
        single = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(single,family)),'infinite exact value')
        single = make_censored_data(interval=reshape([1.0_dp,0.0_dp],[1,2]))
        call check(ieee_is_nan(likelihood(single,family)),'reversed interval')
        x = 1.0_dp
        if (family == 3) x = -1.0_dp
        single = make_censored_data(interval=reshape([-inf,x],[1,2]))
        other = make_censored_data(left=[x])
        call close(likelihood(single,family),likelihood(other,family),0.0_dp,'infinite lower')
        single = make_censored_data(interval=reshape([x,inf],[1,2]))
        other = make_censored_data(right=[x])
        call close(likelihood(single,family),likelihood(other,family),0.0_dp,'infinite upper')
        single = make_censored_data(interval=reshape([x,x],[1,2]))
        other = make_censored_data(uncensored=[x])
        call close(likelihood(single,family),likelihood(other,family),0.0_dp,'equal endpoints exact')
        if (family > 1) then
            single = make_censored_data(uncensored=[0.0_dp])
            call check(likelihood(single,family) == -inf,'zero endpoint density')
            single = make_censored_data(left=[0.0_dp])
            if (family == 2) then
                call check(likelihood(single,family) == -inf,'CDF at lower endpoint')
            else
                call close(likelihood(single,family),0.0_dp,0.0_dp,'CDF at upper endpoint')
            end if
        end if
        lower = [-0.1_dp,0.1_dp]
        upper = [0.1_dp,5.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'all-free bounded fit')
        call check(all(fit%params >= lower) .and. all(fit%params <= upper),'fit within bounds')
        call close(fit%nllf,-likelihood(data,family,fit%params(1),fit%params(2)), &
            0.0_dp,'free fit objective agreement')
        lower = [0.0_dp,0.1_dp]
        upper = [0.0_dp,10.0_dp]
        call fitting(data,family,lower,upper,fit)
        call check(fit%success,'fixed-location scale fit')
        call check(fit%nllf <= -reference(family),'scale fit improves objective')
        expected = 0.6926757812499997_dp
        if (family == 1) expected = 1.1208090746233839_dp
        call close(fit%params(2),expected,2.0e-4_dp,'SciPy fixed-location scale fit')
        select case (family)
        case (1)
            call moyal_fit(data,lower,upper,fit,guess=[0.0_dp,1.3_dp], &
                max_iter=1,tolerance=1.0e-10_dp)
        case (2)
            call levy_fit(data,lower,upper,fit,guess=[0.0_dp,1.3_dp], &
                max_iter=1,tolerance=1.0e-10_dp)
        case (3)
            call levy_l_fit(data,lower,upper,fit,guess=[0.0_dp,1.3_dp], &
                max_iter=1,tolerance=1.0e-10_dp)
        end select
        call check(fit%status == fit_status_max_iter,'optional fit controls forwarded')
        call fitting(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        upper(2) = 0.01_dp
        call fitting(data,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'reversed bounds rejected')
        single = make_censored_data(uncensored=[-1.0_dp])
        if (family == 1) single = make_censored_data(uncensored=[-huge(1.0_dp)])
        if (family == 3) single = make_censored_data(uncensored=[1.0_dp])
        lower = [0.0_dp,1.0_dp]
        call fitting(single,family,lower,lower,fit)
        call check(fit%status == fit_status_no_finite_objective,'impossible fixed objective')
    end do
    call close(levy_censored_loglikelihood(make_censored_data(uncensored=[1.0_dp])), &
        halfnorm_censored_loglikelihood(make_censored_data(uncensored=[1.0_dp]))-log(2.0_dp), &
        3.0e-14_dp,'inverse-square half-normal Jacobian')
    single = make_censored_data(interval=reshape([0.25_dp,4.0_dp],[1,2]))
    other = make_censored_data(interval=reshape([0.5_dp,2.0_dp],[1,2]))
    call close(levy_censored_loglikelihood(single),halfnorm_censored_loglikelihood(other), &
        3.0e-14_dp,'inverse-square interval transformation')
    single = make_censored_data(interval=reshape([-2.0_dp,2.0_dp],[1,2]))
    other = make_censored_data(interval=reshape([exp(-1.0_dp),exp(1.0_dp)],[1,2]))
    call close(moyal_censored_loglikelihood(single),halfnorm_censored_loglikelihood(other), &
        3.0e-14_dp,'negative log half-normal interval transformation')
    call close(levy_censored_loglikelihood(make_censored_data(uncensored=[1.0e300_dp])), &
        -450.0_dp*log(10.0_dp)-log(2.0_dp*acos(-1.0_dp))/2.0_dp, &
        3.0e-13_dp,'large Levy density log')
    single = make_censored_data(right=[1.0e300_dp])
    actual = levy_censored_loglikelihood(single,scale=1.0e-300_dp)
    expected = -300.0_dp*log(10.0_dp)+log(2.0_dp/acos(-1.0_dp))/2.0_dp
    call close(actual,expected,3.0e-13_dp,'Levy SF after ratio overflow')
    single = make_censored_data(left=[1.0e300_dp])
    actual = levy_censored_loglikelihood(single,scale=1.0e-300_dp)
    expected = -sqrt(2.0_dp/acos(-1.0_dp))*1.0e-300_dp
    call close(actual/expected,1.0_dp,2.0e-13_dp,'nearly unit Levy CDF ratio')
    single = make_censored_data(right=[1.0e308_dp])
    actual = moyal_censored_loglikelihood(single,loc=-1.0e308_dp)
    call close(actual/(-1.0e308_dp),1.0_dp,3.0e-15_dp,'finite half coordinate after full overflow')
    single = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    expected = log(nearest(0.0_dp,1.0_dp))-0.5_dp-log(2.0_dp*acos(-1.0_dp))/2.0_dp
    call close(moyal_censored_loglikelihood(single),expected,3.0e-13_dp,'smallest Moyal width')
    single = make_censored_data(right=[2000.0_dp])
    expected = -1000.0_dp+log(2.0_dp/acos(-1.0_dp))/2.0_dp
    call close(moyal_censored_loglikelihood(single),expected,3.0e-13_dp,'Moyal transformed underflow')
    single = make_censored_data(left=[-710.0_dp])
    expected = -exp(710.0_dp-log(2.0_dp))
    call close(moyal_censored_loglikelihood(single)/expected,1.0_dp,2.0e-13_dp, &
        'Moyal finite log CDF beyond direct exponential overflow')
    single = make_censored_data(uncensored=[-710.0_dp])
    call close(moyal_censored_loglikelihood(single)/expected,1.0_dp,2.0e-13_dp, &
        'Moyal finite log density beyond direct exponential overflow')
    single = make_censored_data(uncensored=[nearest(0.0_dp,1.0_dp)])
    call check(levy_censored_loglikelihood(single) == -inf,'unrepresentable Levy log density')
    single = make_censored_data(left=[1.0e-300_dp])
    call close(levy_censored_loglikelihood(single)/(-0.5e300_dp),1.0_dp,2.0e-13_dp, &
        'finite extreme Levy lower tail')
    single = make_censored_data(interval=reshape([-1.0_dp,1.0_dp],[1,2]))
    other = make_censored_data(left=[1.0_dp])
    call close(levy_censored_loglikelihood(single),levy_censored_loglikelihood(other), &
        0.0_dp,'Levy interval clipped to lower support')
    other = make_censored_data(right=[-1.0_dp])
    call close(levy_l_censored_loglikelihood(single),levy_l_censored_loglikelihood(other), &
        0.0_dp,'reflected Levy interval clipped to upper support')
    other = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    call close(levy_censored_loglikelihood(other,-0.2_dp,1.3_dp), &
        levy_l_censored_loglikelihood(reflected,0.2_dp,1.3_dp),3.0e-14_dp,'mixed reflection identity')
    write(*,'(a,i0,a)') 'Normal transform censored: ',checks,' checks'
    if (failures /= 0) error stop 'Normal transform censoring failures'
contains
    pure function likelihood(data, family, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! 1 Moyal, 2 Levy, 3 reflected Levy
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = moyal_censored_loglikelihood(data,loc,scale)
        case (2)
            value = levy_censored_loglikelihood(data,loc,scale)
        case (3)
            value = levy_l_censored_loglikelihood(data,loc,scale)
        end select
    end function likelihood
    subroutine fitting(data, family, lower, upper, result)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! distribution identifier used by this test
        real(dp), intent(in) :: lower(2) !! location and scale lower bounds
        real(dp), intent(in) :: upper(2) !! location and scale upper bounds
        type(fit_result), intent(out) :: result !! bounded fit diagnostics
        select case (family)
        case (1)
            call moyal_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (2)
            call levy_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
        case (3)
            call levy_l_fit(data,lower,upper,result,max_iter=10000,tolerance=1.0e-7_dp)
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
end program test_censored_normal_transform
