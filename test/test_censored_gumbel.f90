! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_gumbel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        gumbel_r_censored_loglikelihood, gumbel_l_censored_loglikelihood, &
        gumbel_r_censored_nnlf, gumbel_l_censored_nnlf, gumbel_r_fit, gumbel_l_fit, &
        gumbel_r_loglikelihood, gumbel_l_loglikelihood, gumbel_r_logpdf, gumbel_l_logpdf, &
        gumbel_r_cdf, gumbel_l_cdf, gumbel_r_logcdf, gumbel_l_logcdf, &
        gumbel_r_logsf, gumbel_l_logsf, fit_status_invalid_input, fit_status_no_finite_objective
    use scifort_math, only : quiet_nan, positive_infinity
    implicit none
    type(censored_data) :: data, other, empty
    type(fit_result) :: fit
    real(dp), parameter :: reference_ll(2) = [-15.063330674681756_dp,-36.48621545554803_dp]
    real(dp), parameter :: reference_loc(2) = [1.187243219112296_dp,2.7284240787521883_dp]
    real(dp), parameter :: reference_scale(2) = [1.4928122886183504_dp,1.7359414219851461_dp]
    real(dp) :: inf, nan, direction, base, value, expected
    integer :: side, checks, failures

    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    do side = 1, 2
        direction = real(3-2*side,dp)
        data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
            right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
        call close(ll(data,0.4_dp,1.3_dp),reference_ll(side),3.0e-12_dp,'mixed reference')
        if (side == 1) then
            value = gumbel_r_censored_nnlf(data,0.4_dp,1.3_dp)
        else
            value = gumbel_l_censored_nnlf(data,0.4_dp,1.3_dp)
        end if
        call close(value,-reference_ll(side),3.0e-12_dp,'negative log-likelihood')
        call fit_censored(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
        call check(fit%success,'mixed fit converges')
        call close(fit%params(1),reference_loc(side),4.0e-4_dp,'SciPy fitted location')
        call close(fit%params(2),reference_scale(side),4.0e-4_dp,'SciPy fitted scale')
        call close(fit%nllf,-ll(data,fit%params(1),fit%params(2)),1.0e-13_dp,'fit objective')
        call fit_censored(data,[0.4_dp,1.3_dp],[0.4_dp,1.3_dp],fit)
        call check(fit%success,'both parameters fixed')
        call close(fit%nllf,-reference_ll(side),3.0e-12_dp,'fixed objective')
        base = fit%nllf
        call fit_censored(data,[0.4_dp,0.05_dp],[0.4_dp,10.0_dp],fit)
        call check(fit%success .and. fit%params(1) == 0.4_dp,'fixed location')
        call check(fit%nllf <= base,'free scale improves fixed objective')
        call fit_censored(data,[-5.0_dp,1.3_dp],[8.0_dp,1.3_dp],fit)
        call check(fit%success .and. fit%params(2) == 1.3_dp,'fixed scale')
        call check(fit%nllf <= base,'free location improves fixed objective')
        other = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp],right=[-0.3_dp], &
            left=[-3.5_dp],interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
        if (side == 1) then
            value = gumbel_l_censored_loglikelihood(other,-0.4_dp,1.3_dp)
        else
            value = gumbel_r_censored_loglikelihood(other,-0.4_dp,1.3_dp)
        end if
        call close(ll(data,0.4_dp,1.3_dp),value,2.0e-13_dp,'mixed reflection identity')
        data = make_censored_data(uncensored=[-2.0_dp,0.0_dp,0.5_dp,4.0_dp])
        if (side == 1) then
            value = gumbel_r_loglikelihood(data%uncensored,0.4_dp,1.3_dp)
            call gumbel_r_fit(data%uncensored,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
        else
            value = gumbel_l_loglikelihood(data%uncensored,0.4_dp,1.3_dp)
            call gumbel_l_fit(data%uncensored,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
        end if
        call close(ll(data,0.4_dp,1.3_dp),value,3.0e-14_dp,'exact array likelihood agreement')
        base = fit%nllf
        call fit_censored(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit,1.0e-8_dp)
        call close(fit%nllf,base,2.0e-12_dp,'exact fit overload agreement')
        call set_interval(-1.0_dp,1.0_dp)
        base = ll(data)
        if (side == 1) then
            expected = log(gumbel_r_cdf(1.0_dp)-gumbel_r_cdf(-1.0_dp))
        else
            expected = log(gumbel_l_cdf(1.0_dp)-gumbel_l_cdf(-1.0_dp))
        end if
        call close(base,expected,2.0e-15_dp,'ordinary interval probability')
        call set_interval(-2.0_dp,2.0_dp)
        value = ll(data)
        call check(value >= base .and. value <= 0.0_dp,'nested interval masses')
        call set_interval(0.0_dp,nearest(0.0_dp,1.0_dp))
        call close(ll(data),log(nearest(0.0_dp,1.0_dp))-1.0_dp,2.0e-13_dp,'subnormal width')
        call set_interval(1.0e-300_dp,2.0e-300_dp)
        call close(ll(data,scale=1.0e300_dp),log(1.0e-300_dp)-log(1.0e300_dp)-1.0_dp, &
            3.0e-13_dp,'standardized interval width underflows')
        if (side == 1) then
            call set_interval(1000.0_dp,nearest(1000.0_dp,1.0_dp))
        else
            call set_interval(-nearest(1000.0_dp,1.0_dp),-1000.0_dp)
        end if
        call close(ll(data),-1029.8053287640778_dp,3.0e-13_dp,'adjacent light-tail interval')
        if (side == 1) then
            call set_interval(-5.0_dp,nearest(-5.0_dp,1.0_dp))
        else
            call set_interval(-nearest(-5.0_dp,1.0_dp),5.0_dp)
        end if
        call close(ll(data),-178.07051813057382_dp,2.0e-12_dp,'adjacent steep-tail interval')
        call set_interval(-1.0e308_dp,1.0e308_dp)
        call close(ll(data),0.0_dp,0.0_dp,'raw interval width overflows')
        call close(ll(data,scale=1.0e308_dp),base,1.0e-13_dp,'overflowing width with finite scaled mass')
        call set_interval(-1.0e-300_dp,1.0e-300_dp)
        call close(ll(data,scale=1.0e-300_dp),base,1.0e-13_dp,'tiny-scale interval identity')
        data = make_censored_data(uncensored=[1.0e308_dp])
        expected = -2.0_dp*direction-exp(-2.0_dp*direction)-log(1.0e308_dp)
        call close(ll(data,loc=-1.0e308_dp,scale=1.0e308_dp),expected,2.0e-13_dp, &
            'overflowing displacement with finite standardized value')
        data = make_censored_data(uncensored=[nearest(huge(1.0_dp),-1.0_dp)])
        expected = -direction-exp(-direction)-log(huge(1.0_dp))
        call close(ll(data,loc=-1.5_dp*spacing(huge(1.0_dp)),scale=huge(1.0_dp)), &
            expected,2.0e-13_dp,'rounded subtraction-overflow threshold')
        data = make_censored_data(uncensored=[0.0_dp])
        call close(ll(data,scale=nearest(0.0_dp,1.0_dp)), &
            -1.0_dp-log(nearest(0.0_dp,1.0_dp)),2.0e-13_dp,'smallest scale at location')
        data = make_censored_data(uncensored=[-direction*1000.0_dp])
        call check(ll(data) == -inf,'unrepresentable steep-tail log density')
        data = make_censored_data(uncensored=[direction*1.0e308_dp])
        call check(ll(data,scale=1.0e-300_dp) == -inf,'unrepresentable standardized log density')
        if (side == 1) then
            data = make_censored_data(right=[1000.0_dp])
        else
            data = make_censored_data(left=[-1000.0_dp])
        end if
        call close(ll(data),-1000.0_dp,0.0_dp,'finite log tail after probability underflows')
        if (side == 1) then
            data = make_censored_data(left=[-1000.0_dp])
        else
            data = make_censored_data(right=[1000.0_dp])
        end if
        call check(ll(data) == -inf,'unrepresentable steep log tail')
        if (side == 1) then
            data = make_censored_data(right=[-1000.0_dp])
        else
            data = make_censored_data(left=[1000.0_dp])
        end if
        call close(ll(data),0.0_dp,0.0_dp,'certain tail at extreme observation')
        call set_interval(-inf,0.4_dp)
        if (side == 1) then
            expected = gumbel_r_logcdf(0.4_dp,0.4_dp,1.3_dp)
        else
            expected = gumbel_l_logcdf(0.4_dp,0.4_dp,1.3_dp)
        end if
        call close(ll(data,0.4_dp,1.3_dp),expected,0.0_dp,'infinite lower endpoint')
        call set_interval(0.4_dp,inf)
        if (side == 1) then
            expected = gumbel_r_logsf(0.4_dp,0.4_dp,1.3_dp)
        else
            expected = gumbel_l_logsf(0.4_dp,0.4_dp,1.3_dp)
        end if
        call close(ll(data,0.4_dp,1.3_dp),expected,0.0_dp,'infinite upper endpoint')
        call set_interval(0.4_dp,0.4_dp)
        if (side == 1) then
            expected = gumbel_r_logpdf(0.4_dp,0.4_dp,1.3_dp)
        else
            expected = gumbel_l_logpdf(0.4_dp,0.4_dp,1.3_dp)
        end if
        call close(ll(data,0.4_dp,1.3_dp),expected,0.0_dp,'equal endpoints are exact')
        call close(ll(empty),0.0_dp,0.0_dp,'empty likelihood')
        call check(ieee_is_nan(ll(empty,loc=nan)),'NaN location')
        call check(ieee_is_nan(ll(empty,loc=inf)),'infinite location')
        call check(ieee_is_nan(ll(empty,scale=nan)),'NaN scale')
        call check(ieee_is_nan(ll(empty,scale=inf)),'infinite scale')
        call check(ieee_is_nan(ll(empty,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(ll(empty,scale=-1.0_dp)),'negative scale')
        data = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(ll(data)),'NaN exact observation')
        call fit_censored(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
        call check(fit%status == fit_status_invalid_input,'invalid data fit')
        data = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(data)),'infinite exact observation')
        call set_interval(1.0_dp,-1.0_dp)
        call check(ieee_is_nan(ll(data)),'reversed interval')
        call fit_censored(empty,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
        call check(fit%status == fit_status_invalid_input,'empty fit sample')
        data = make_censored_data(uncensored=[0.0_dp])
        call fit_censored(data,[1.0_dp,0.05_dp],[-1.0_dp,10.0_dp],fit)
        call check(fit%status == fit_status_invalid_input,'reversed fit bounds')
        call fit_censored(data,[0.0_dp,-1.0_dp],[0.0_dp,-1.0_dp],fit)
        call check(fit%status == fit_status_no_finite_objective,'invalid fixed scale')
    end do
    if (failures /= 0) error stop 1
    print '(a,i0,a)', 'Censored Gumbel likelihood and fitting: ',checks,' checks passed'
contains
    function ll(sample,loc,scale) result(value)
        implicit none
        type(censored_data), intent(in) :: sample !! observations for the selected family
        real(dp), intent(in), optional :: loc !! finite location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        if (side == 1) then
            value = gumbel_r_censored_loglikelihood(sample,loc,scale)
        else
            value = gumbel_l_censored_loglikelihood(sample,loc,scale)
        end if
    end function ll

    subroutine fit_censored(sample,lower,upper,result,tolerance)
        implicit none
        type(censored_data), intent(in) :: sample !! observations for the selected family
        real(dp), intent(in) :: lower(2) !! lower location/scale bounds
        real(dp), intent(in) :: upper(2) !! upper bounds; equal bounds fix parameters
        type(fit_result), intent(out) :: result !! optimizer result and diagnostics
        real(dp), intent(in), optional :: tolerance !! parameter convergence tolerance
        real(dp) :: tol
        tol = 1.0e-10_dp
        if (present(tolerance)) tol = tolerance
        if (side == 1) then
            call gumbel_r_fit(sample,lower,upper,result,guess=[0.4_dp,1.3_dp], &
                max_iter=2000,tolerance=tol)
        else
            call gumbel_l_fit(sample,lower,upper,result,guess=[0.4_dp,1.3_dp], &
                max_iter=2000,tolerance=tol)
        end if
    end subroutine fit_censored

    subroutine set_interval(a,b)
        implicit none
        real(dp), intent(in) :: a !! lower test endpoint
        real(dp), intent(in) :: b !! upper test endpoint
        data = make_censored_data(interval=reshape([a,b],[1,2]))
    end subroutine set_interval

    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! expected invariant
        character(*), intent(in) :: label !! failure description
        checks = checks + 1
        if (ok) return
        failures = failures + 1
        print '(a,i0,2a)', 'FAIL Gumbel side ',side,': ',label
    end subroutine check

    subroutine close(actual,expected,atol,label)
        implicit none
        real(dp), intent(in) :: actual !! computed finite value
        real(dp), intent(in) :: expected !! independent reference or identity value
        real(dp), intent(in) :: atol !! regime-specific absolute tolerance
        character(*), intent(in) :: label !! failure description
        call check(ieee_is_finite(actual),label//' is finite')
        if (.not. ieee_is_finite(actual)) return
        call check(abs(actual-expected) <= atol,label)
        if (abs(actual-expected) > atol) print *,actual,expected
    end subroutine close
end program test_censored_gumbel
