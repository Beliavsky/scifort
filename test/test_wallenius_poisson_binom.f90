! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_wallenius_poisson_binom
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_wallenius_poisson_binom_reference
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures

    failures=0
    call test_reference(failures)
    call test_identities_endpoints(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures==0) then
        print '(a)', 'test_wallenius_poisson_binom: PASS'
    else
        print '(a,1x,i0)', 'test_wallenius_poisson_binom: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer :: i

        do i=1,size(wall_x)
            call check_close('Wallenius pmf',nchypergeom_wallenius_pmf(wall_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_pmf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Wallenius logpmf',nchypergeom_wallenius_logpmf(wall_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_logpmf(i),3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('Wallenius cdf',nchypergeom_wallenius_cdf(wall_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_cdf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Wallenius sf',nchypergeom_wallenius_sf(wall_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_sf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Wallenius logcdf',nchypergeom_wallenius_logcdf(wall_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_logcdf(i),3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('Wallenius logsf',nchypergeom_wallenius_logsf(wall_x(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_logsf(i),3.0e-13_dp,3.0e-13_dp,failures)
        end do
        do i=1,size(reference_probs)
            call check_close('Wallenius ppf',nchypergeom_wallenius_ppf(reference_probs(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_ppf(i),0.0_dp,0.0_dp,failures)
            call check_close('Wallenius isf',nchypergeom_wallenius_isf(reference_probs(i),30.0_dp,12.0_dp,7.0_dp,2.3_dp), &
                wall_isf(i),0.0_dp,0.0_dp,failures)
        end do

        do i=1,size(pb_x)
            call check_close('PB pmf',poisson_binom_pmf(pb_x(i),pb_p),pb_pmf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('PB logpmf',poisson_binom_logpmf(pb_x(i),pb_p),pb_logpmf(i),3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('PB cdf',poisson_binom_cdf(pb_x(i),pb_p),pb_cdf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('PB sf',poisson_binom_sf(pb_x(i),pb_p),pb_sf(i),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('PB logcdf',poisson_binom_logcdf(pb_x(i),pb_p),pb_logcdf(i),3.0e-13_dp,3.0e-13_dp,failures)
            call check_close('PB logsf',poisson_binom_logsf(pb_x(i),pb_p),pb_logsf(i),3.0e-13_dp,3.0e-13_dp,failures)
        end do
        do i=1,size(reference_probs)
            call check_close('PB ppf',poisson_binom_ppf(reference_probs(i),pb_p),pb_ppf(i),0.0_dp,0.0_dp,failures)
            call check_close('PB isf',poisson_binom_isf(reference_probs(i),pb_p),pb_isf(i),0.0_dp,0.0_dp,failures)
        end do
    end subroutine test_reference

    subroutine test_identities_endpoints(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: peq(6)=[0.3_dp,0.3_dp,0.3_dp,0.3_dp,0.3_dp,0.3_dp]
        real(dp) :: total
        integer :: i, k

        total=0.0_dp
        do k=0,7
            total=total+nchypergeom_wallenius_pmf(k,30,12,7,2.3_dp)
            call check_close('Wallenius hypergeom pmf',nchypergeom_wallenius_pmf(k,30,12,7,1.0_dp), &
                hypergeom_pmf(k,30,12,7),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('Wallenius hypergeom cdf',nchypergeom_wallenius_cdf(k,30,12,7,1.0_dp), &
                hypergeom_cdf(k,30,12,7),3.0e-15_dp,3.0e-13_dp,failures)
        end do
        call check_close('Wallenius normalization',total,1.0_dp,3.0e-15_dp,3.0e-14_dp,failures)
        call check_close('Wallenius shifted',nchypergeom_wallenius_cdf(3.5_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.5_dp), &
            nchypergeom_wallenius_cdf(3.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_true('Wallenius invalid odds',ieee_is_nan(nchypergeom_wallenius_pmf(2,30,12,7,0.0_dp)),failures)

        total=0.0_dp
        do k=0,size(pb_p)
            total=total+poisson_binom_pmf(k,pb_p)
        end do
        call check_close('PB normalization',total,1.0_dp,3.0e-15_dp,3.0e-14_dp,failures)
        do k=0,size(peq)
            call check_close('PB binomial pmf',poisson_binom_pmf(k,peq), &
                binomial_pmf(real(k,dp),real(size(peq),dp),0.3_dp),3.0e-15_dp,3.0e-13_dp,failures)
            call check_close('PB binomial cdf',poisson_binom_cdf(k,peq), &
                binomial_cdf(real(k,dp),real(size(peq),dp),0.3_dp),3.0e-15_dp,3.0e-13_dp,failures)
        end do
        call check_close('PB shifted',poisson_binom_cdf(2.5_dp,pb_p,0.5_dp),poisson_binom_cdf(2.0_dp,pb_p), &
            0.0_dp,0.0_dp,failures)
        do i=1,size(pb_p)
            call check_true('PB valid p',pb_p(i)>=0.0_dp .and. pb_p(i)<=1.0_dp,failures)
        end do
        call check_true('PB invalid p',ieee_is_nan(poisson_binom_pmf(2,[0.2_dp,1.1_dp])),failures)
    end subroutine test_identities_endpoints

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h=2.0e-6_dp
        integer, parameter :: wd(5)=[1,2,3,4,5]
        integer, parameter :: pd(6)=[0,1,2,2,3,4]
        real(dp) :: sw(1), sp(size(pb_p)), plusp(size(pb_p)), minusp(size(pb_p)), fd
        integer :: j

        sw=nchypergeom_wallenius_score(wd,30,12,7,2.3_dp)
        fd=centered(nchypergeom_wallenius_loglikelihood(wd,30,12,7,2.3_dp+h), &
            nchypergeom_wallenius_loglikelihood(wd,30,12,7,2.3_dp-h),h)
        call check_close('Wallenius score odds',sw(1),fd,3.0e-7_dp,3.0e-7_dp,failures)

        sp=poisson_binom_score(pd,pb_p)
        do j=1,size(pb_p)
            plusp=pb_p; minusp=pb_p; plusp(j)=plusp(j)+h; minusp(j)=minusp(j)-h
            fd=centered(poisson_binom_loglikelihood(pd,plusp),poisson_binom_loglikelihood(pd,minusp),h)
            call check_close('PB score p',sp(j),fd,3.0e-7_dp,3.0e-7_dp,failures)
        end do
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        integer, parameter :: wd(5)=[1,2,3,4,5]
        integer, parameter :: pd(6)=[0,1,2,2,3,4]
        real(dp) :: actual, expected, q, den
        real(dp) :: pb_lower(6), pb_upper(6)
        type(rng_state) :: state, reference
        type(fit_result) :: result
        integer :: j, aleft, bleft, count

        call rng_seed(state,24680); call rng_seed(reference,24680)
        actual=nchypergeom_wallenius_rvs(state,30.0_dp,12.0_dp,7.0_dp,2.3_dp)
        aleft=12; bleft=18; count=0
        do j=1,7
            den=2.3_dp*real(aleft,dp)+real(bleft,dp); q=2.3_dp*real(aleft,dp)/den
            if (rng_uniform(reference)<q) then; count=count+1; aleft=aleft-1; else; bleft=bleft-1; end if
        end do
        expected=real(count,dp)
        call check_close('Wallenius rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call rng_seed(state,97531); call rng_seed(reference,97531)
        actual=poisson_binom_rvs(state,pb_p); count=0
        do j=1,size(pb_p)
            if (rng_uniform(reference)<pb_p(j)) count=count+1
        end do
        call check_close('PB rvs',actual,real(count,dp),0.0_dp,0.0_dp,failures)

        call rng_seed(state,13579); call rng_seed(reference,13579)
        actual=nchypergeom_wallenius_rvs(state,30.0_dp,12.0_dp,7.0_dp,0.0_dp)
        call check_true('invalid Wallenius rvs NaN',ieee_is_nan(actual),failures)
        call check_close('invalid Wallenius preserves state',rng_uniform(state),rng_uniform(reference), &
            0.0_dp,0.0_dp,failures)

        call nchypergeom_wallenius_fit(wd,[30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp], &
            [30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp],result)
        call check_true('Wallenius fixed fit',result%success,failures)
        pb_lower=[pb_p,0.0_dp]; pb_upper=pb_lower
        call poisson_binom_fit(pd,pb_lower,pb_upper,result)
        call check_true('PB fixed fit',result%success,failures)

        call check_close('C ABI Wallenius pmf', &
            scifort_nchypergeom_wallenius_pmf_f64(2.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp), &
            nchypergeom_wallenius_pmf(2,30,12,7,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Wallenius cdf', &
            scifort_nchypergeom_wallenius_cdf_f64(2.0_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp), &
            nchypergeom_wallenius_cdf(2,30,12,7,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI Wallenius ppf', &
            scifort_nchypergeom_wallenius_ppf_f64(0.7_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp,0.0_dp), &
            nchypergeom_wallenius_ppf(0.7_dp,30.0_dp,12.0_dp,7.0_dp,2.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI PB pmf',scifort_poisson_binom_pmf_f64(2.0_dp,pb_p,size(pb_p,kind=8),0.0_dp), &
            poisson_binom_pmf(2,pb_p),0.0_dp,0.0_dp,failures)
        call check_close('C ABI PB cdf',scifort_poisson_binom_cdf_f64(2.0_dp,pb_p,size(pb_p,kind=8),0.0_dp), &
            poisson_binom_cdf(2,pb_p),0.0_dp,0.0_dp,failures)
        call check_close('C ABI PB ppf',scifort_poisson_binom_ppf_f64(0.7_dp,pb_p,size(pb_p,kind=8),0.0_dp), &
            poisson_binom_ppf(0.7_dp,pb_p),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus_value,minus_value,h) result(value)
        real(dp), intent(in) :: plus_value !! function value at theta+h
        real(dp), intent(in) :: minus_value !! function value at theta-h
        real(dp), intent(in) :: h !! positive centered-difference step
        real(dp) :: value
        value=(plus_value-minus_value)/(2.0_dp*h)
    end function centered

end program test_wallenius_poisson_binom
