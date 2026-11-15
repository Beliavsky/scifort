! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_erlang
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        erlang_censored_loglikelihood, erlang_censored_nnlf, erlang_fit, gamma_fit, &
        gamma_censored_loglikelihood, exponential_censored_loglikelihood, &
        erlang_logpdf, erlang_logcdf, erlang_logsf, erlang_loglikelihood, &
        fit_status_invalid_input, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty
    type(fit_result) :: fit, parent
    real(dp) :: a, x, expected, inf, nan, value, lower(3), upper(3)
    integer :: i, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    call close(erlang_censored_loglikelihood(data,2.3_dp,scale=1.3_dp), &
        -16.305176559791958_dp,6.0e-14_dp,'mixed independent mpmath reference')
    do i = 1, 12
        a = 0.25_dp*real(i,dp)
        value = erlang_censored_loglikelihood(data,a,0.0_dp,1.3_dp)
        call close(value,gamma_censored_loglikelihood(data,a,0.0_dp,1.3_dp), &
            0.0_dp,'mixed gamma identity, including noninteger shapes')
        call close(erlang_censored_nnlf(data,a,scale=1.3_dp),-value,0.0_dp,'negative likelihood')
        x = 0.3_dp*real(i,dp)
        single = make_censored_data(uncensored=[x])
        call close(erlang_censored_loglikelihood(single,a),erlang_logpdf(x,a), &
            2.0e-14_dp,'exact scalar agreement and defaults')
        single = make_censored_data(left=[x])
        call close(erlang_censored_loglikelihood(single,a),erlang_logcdf(x,a), &
            2.0e-14_dp,'left scalar agreement')
        single = make_censored_data(right=[x])
        call close(erlang_censored_loglikelihood(single,a),erlang_logsf(x,a), &
            2.0e-14_dp,'right scalar agreement')
    end do
    call close(erlang_censored_loglikelihood(data,1.0_dp,scale=1.3_dp), &
        exponential_censored_loglikelihood(data,scale=1.3_dp),0.0_dp,'exponential identity')
    ! Shape two has S(x)=exp(-x)*(1+x), independent of incomplete gamma.
    single = make_censored_data(right=[1000.0_dp])
    call close(erlang_censored_loglikelihood(single,2.0_dp), &
        -1000.0_dp+log(1001.0_dp),2.0e-13_dp,'shape-two extreme survival')
    single = make_censored_data(interval=reshape([2.0_dp,3.0_dp],[1,2]))
    expected = log(3.0_dp*exp(-2.0_dp)-4.0_dp*exp(-3.0_dp))
    call close(erlang_censored_loglikelihood(single,2.0_dp),expected,2.0e-14_dp, &
        'shape-two independent interval formula')
    single = make_censored_data(interval=reshape([1.0_dp,nearest(1.0_dp,2.0_dp)],[1,2]))
    call close(erlang_censored_loglikelihood(single,2.0_dp), &
        -1.0_dp+log(spacing(1.0_dp)),2.0e-14_dp,'adjacent interval density integral')
    single = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    call close(erlang_censored_loglikelihood(single,2.0_dp), &
        2.0_dp*log(nearest(0.0_dp,1.0_dp))-log(2.0_dp),5.0e-13_dp,'subnormal interval')
    single = make_censored_data(uncensored=[1.0_dp,2.0_dp,3.0_dp])
    call close(erlang_censored_loglikelihood(single,2.3_dp,scale=1.3_dp), &
        erlang_loglikelihood(single%uncensored,2.3_dp,scale=1.3_dp),2.0e-14_dp, &
        'array likelihood agreement')
    call erlang_fit(single,[1.0_dp,0.0_dp,0.1_dp],[1.0_dp,0.0_dp,10.0_dp],fit)
    call check(fit%success,'shape-one fit converges')
    call close(fit%params(3),2.0_dp,3.0e-6_dp,'analytic exponential MLE')
    call erlang_fit(single%uncensored,[2.3_dp,0.0_dp,1.3_dp], &
        [2.3_dp,0.0_dp,1.3_dp],parent)
    call erlang_fit(single,[2.3_dp,0.0_dp,1.3_dp],[2.3_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success .and. parent%success,'array and censored fixed fits converge')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'array fit compatibility')
    lower = [0.1_dp,0.0_dp,0.05_dp]
    upper = [6.0_dp,0.0_dp,10.0_dp]
    call erlang_fit(data,lower,upper,fit)
    call gamma_fit(data,lower,upper,parent)
    call check(fit%success .and. parent%success,'mixed continuous-shape fit converges')
    call check(all(fit%params == parent%params),'fit preserves gamma parameters')
    call close(fit%params(1),0.9927593350620091_dp,4.0e-4_dp,'SciPy fitted shape')
    call close(fit%params(3),2.232420387730313_dp,4.0e-4_dp,'SciPy fitted scale')
    call close(fit%nllf,parent%nllf,0.0_dp,'fit preserves objective')
    call check(fit%evaluations == parent%evaluations,'fit preserves evaluations')
    call check(fit%iterations == parent%iterations,'fit preserves iterations')
    call erlang_fit(data,lower,upper,fit,guess=[2.3_dp,0.0_dp,1.3_dp], &
        max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls forwarded')
    upper(2) = 0.1_dp
    call erlang_fit(data,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'free location rejected')
    upper(2) = 0.0_dp
    single = make_censored_data(uncensored=[0.0_dp])
    call erlang_fit(single,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'exact support endpoint excluded from fit')
    single = make_censored_data(interval=reshape([0.0_dp,0.0_dp],[1,2]))
    call erlang_fit(single,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'equal endpoint excluded from fit')
    call erlang_fit(empty,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call close(erlang_censored_loglikelihood(empty,2.0_dp),0.0_dp,0.0_dp,'valid empty sample')
    call check(ieee_is_nan(erlang_censored_loglikelihood(empty,0.0_dp)),'zero shape invalid')
    call check(ieee_is_nan(erlang_censored_loglikelihood(empty,-1.0_dp)),'negative shape invalid')
    call check(ieee_is_nan(erlang_censored_loglikelihood(empty,nan)),'NaN shape invalid')
    call check(ieee_is_nan(erlang_censored_loglikelihood(empty,inf)),'infinite shape invalid')
    call check(ieee_is_nan(erlang_censored_loglikelihood(empty,2.0_dp,scale=0.0_dp)), &
        'zero scale invalid even with empty data')
    call check(ieee_is_nan(erlang_censored_loglikelihood(empty,2.0_dp,loc=inf)), &
        'infinite location invalid')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(erlang_censored_loglikelihood(single,2.0_dp)),'invalid observation')
    single = make_censored_data(uncensored=[0.0_dp])
    call check(erlang_censored_loglikelihood(single,0.5_dp) == inf,'singular endpoint density')
    call close(erlang_censored_loglikelihood(single,1.0_dp),0.0_dp,0.0_dp,'unit-shape endpoint')
    call check(erlang_censored_loglikelihood(single,2.0_dp) == -inf,'zero endpoint density')
    single = make_censored_data(left=[0.0_dp])
    call check(erlang_censored_loglikelihood(single,2.0_dp) == -inf,'zero endpoint CDF')
    single = make_censored_data(right=[0.0_dp])
    call close(erlang_censored_loglikelihood(single,2.0_dp),0.0_dp,0.0_dp,'unit endpoint SF')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(erlang_censored_loglikelihood(single,2.0_dp),erlang_logcdf(2.0_dp,2.0_dp), &
        2.0e-14_dp,'infinite lower endpoint')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(erlang_censored_loglikelihood(single,2.0_dp),erlang_logsf(2.0_dp,2.0_dp), &
        2.0e-14_dp,'infinite upper endpoint')
    write(*,'(a,i0,a)') 'Erlang censored: ',checks,' checks'
    if (failures /= 0) error stop 'Erlang censoring failures'
contains
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
    end subroutine close
end program test_censored_erlang
