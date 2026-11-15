! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_laplace_asymmetric
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        all_ll => laplace_asymmetric_censored_loglikelihood, &
        annlf => laplace_asymmetric_censored_nnlf, afit => laplace_asymmetric_fit, &
        laplace_asymmetric_logpdf, laplace_asymmetric_logcdf, laplace_asymmetric_logsf, &
        laplace_asymmetric_loglikelihood, laplace_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, reflected, empty
    type(fit_result) :: fit, parent
    real(dp) :: x, k, inf, nan, actual, expected, a, b, base, lowk, highk
    integer :: i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    lowk = nearest(0.0_dp,1.0_dp)
    highk = huge(1.0_dp)
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    actual = all_ll(data,1.7_dp,0.4_dp,1.3_dp)
    call close(actual,-23.31171335066417_dp,9.0e-14_dp,'mixed mpmath reference')
    call close(annlf(data,1.7_dp,0.4_dp,1.3_dp),-actual,0.0_dp,'negative likelihood')
    reflected = make_censored_data(uncensored=-data%uncensored,left=-data%right, &
        right=-data%left,interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
    call close(all_ll(reflected,1.0_dp/1.7_dp,-0.4_dp,1.3_dp),actual, &
        5.0e-14_dp,'reciprocal-shape reflection')
    call close(all_ll(data,1.0_dp,0.4_dp,1.3_dp), &
        laplace_censored_loglikelihood(data,0.4_dp,1.3_dp),5.0e-14_dp,'Laplace identity')
    do j = 1, 3
        k = 0.5_dp*real(j,dp)
        do i = -5, 5
            x = 0.7_dp*real(i,dp)
            single = make_censored_data(uncensored=[x])
            call close(all_ll(single,k),laplace_asymmetric_logpdf(x,k), &
                2.0e-14_dp,'scalar density and defaults')
            single = make_censored_data(left=[x])
            call close(all_ll(single,k),laplace_asymmetric_logcdf(x,k), &
                2.0e-14_dp,'scalar CDF')
            single = make_censored_data(right=[x])
            call close(all_ll(single,k),laplace_asymmetric_logsf(x,k), &
                2.0e-14_dp,'scalar SF')
        end do
    end do
    single = make_censored_data(uncensored=[-2.0_dp,-1.0_dp,1.0_dp,3.0_dp])
    call close(all_ll(single,2.0_dp,scale=1.3_dp), &
        laplace_asymmetric_loglikelihood(single%uncensored,2.0_dp,scale=1.3_dp), &
        2.0e-14_dp,'array agreement')
    call afit(single,[2.0_dp,0.0_dp,0.05_dp],[2.0_dp,0.0_dp,10.0_dp],fit)
    call check(fit%success,'fixed-shape fixed-location fit')
    call close(fit%params(3),2.375_dp,4.0e-6_dp,'analytic weighted absolute-deviation scale MLE')
    call afit(single%uncensored,[2.0_dp,0.0_dp,1.3_dp],[2.0_dp,0.0_dp,1.3_dp],parent)
    call afit(single,[2.0_dp,0.0_dp,1.3_dp],[2.0_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success .and. parent%success,'array and censored fixed fits')
    call close(fit%nllf,parent%nllf,2.0e-14_dp,'array fit unchanged')
    call afit(data,[1.7_dp,0.4_dp,1.3_dp],[1.7_dp,0.4_dp,1.3_dp],fit)
    call check(fit%success,'all parameters fixed')
    base = fit%nllf
    call afit(data,[1.7_dp,0.4_dp,0.05_dp],[1.7_dp,0.4_dp,10.0_dp],fit)
    call check(fit%success,'mixed free scale fit')
    call close(fit%params(3),3.6634935386148726_dp,4.0e-4_dp,'SciPy mixed fitted scale')
    call afit(data,[0.1_dp,0.4_dp,0.05_dp],[4.0_dp,0.4_dp,10.0_dp],fit)
    call check(fit%success,'mixed fixed-location shape/scale fit')
    call close(fit%params(1),0.4577391395316225_dp,4.0e-4_dp,'SciPy mixed fitted shape')
    call close(fit%params(3),1.0270210720704078_dp,4.0e-4_dp,'SciPy mixed shape/scale fitted scale')
    call afit(data,[0.1_dp,-2.0_dp,0.05_dp],[4.0_dp,4.0_dp,10.0_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'all parameters free')
    call close(fit%nllf,-all_ll(data,fit%params(1),fit%params(2),fit%params(3)), &
        0.0_dp,'mixed objective agreement')
    call afit(data,[1.7_dp,-2.0_dp,1.3_dp],[1.7_dp,4.0_dp,1.3_dp],fit)
    call check(fit%success .and. fit%nllf <= base,'free cusp with shape and scale fixed')
    call afit(data,[0.1_dp,-2.0_dp,0.05_dp],[4.0_dp,4.0_dp,10.0_dp],fit, &
        guess=[1.7_dp,0.4_dp,1.3_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'optional controls forwarded')
    call afit(empty,[1.0_dp,0.0_dp,0.05_dp],[1.0_dp,0.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'empty fit rejected')
    call close(all_ll(empty,1.0_dp),0.0_dp,0.0_dp,'valid empty likelihood')
    call check(ieee_is_nan(all_ll(empty,0.0_dp)),'zero shape')
    call check(ieee_is_nan(all_ll(empty,-1.0_dp)),'negative shape')
    call check(ieee_is_nan(all_ll(empty,nan)),'NaN shape')
    call check(ieee_is_nan(all_ll(empty,inf)),'infinite shape')
    call check(ieee_is_nan(all_ll(empty,1.0_dp,loc=inf)),'infinite location')
    call check(ieee_is_nan(all_ll(empty,1.0_dp,scale=nan)),'NaN scale')
    call check(ieee_is_nan(all_ll(empty,1.0_dp,scale=inf)),'infinite scale')
    call check(ieee_is_nan(all_ll(empty,1.0_dp,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(all_ll(empty,1.0_dp,scale=-1.0_dp)),'negative scale')
    single = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(all_ll(single,1.0_dp)),'invalid observation')
    single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
    call close(all_ll(single,2.0_dp),laplace_asymmetric_logcdf(2.0_dp,2.0_dp), &
        2.0e-14_dp,'infinite lower endpoint')
    single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
    call close(all_ll(single,2.0_dp),laplace_asymmetric_logsf(2.0_dp,2.0_dp), &
        2.0e-14_dp,'infinite upper endpoint')
    single = make_censored_data(interval=reshape([0.0_dp,0.0_dp],[1,2]))
    call close(all_ll(single,2.0_dp),-log(2.5_dp),2.0e-14_dp,'exact cusp via equal endpoints')
    single = make_censored_data(right=[1000.0_dp])
    call close(all_ll(single,2.0_dp),-2000.0_dp-log(5.0_dp),3.0e-12_dp,'extreme SF')
    single = make_censored_data(left=[-1000.0_dp])
    call close(all_ll(single,2.0_dp),-500.0_dp+log(0.8_dp),1.0e-12_dp,'extreme CDF')
    single = make_censored_data(left=[20.0_dp])
    expected = -0.2_dp*exp(-40.0_dp)
    call close(all_ll(single,2.0_dp)/expected,1.0_dp,3.0e-13_dp,'near-one CDF')
    single = make_censored_data(interval=reshape([-80.0_dp,20.0_dp],[1,2]))
    call close(all_ll(single,2.0_dp)/(-exp(-40.0_dp)),1.0_dp,3.0e-13_dp, &
        'almost-full-support interval retains both tails')
    single = make_censored_data(interval=reshape([0.0_dp,40.0_dp],[1,2]))
    call close(all_ll(single,1.0e-300_dp,scale=1.0e-300_dp)/(-exp(-40.0_dp)), &
        1.0_dp,3.0e-13_dp,'same-side almost-unit mass')
    a = 1000.0_dp
    b = nearest(a,1.0_dp)
    single = make_censored_data(interval=reshape([a,b],[1,2]))
    expected = laplace_asymmetric_logpdf(a,2.0_dp)+log(b-a)
    call close(all_ll(single,2.0_dp),expected,4.0e-12_dp,'adjacent far-tail interval')
    a = lowk
    single = make_censored_data(interval=reshape([-a,a],[1,2]))
    call close(all_ll(single,2.0_dp),log(2.0_dp)+log(a)-log(2.5_dp), &
        3.0e-13_dp,'subnormal crossing interval')
    single = make_censored_data(interval=reshape([1.0e-300_dp,2.0e-300_dp],[1,2]))
    expected = log(1.0e-300_dp)-log(1.0e300_dp)-log(2.5_dp)
    call close(all_ll(single,2.0_dp,scale=1.0e300_dp),expected, &
        5.0e-13_dp,'underflowed scaled positive interval')
    single = make_censored_data(interval=reshape([-1.0e-300_dp,1.0e-300_dp],[1,2]))
    call close(all_ll(single,2.0_dp,scale=1.0e300_dp),expected+log(2.0_dp), &
        5.0e-13_dp,'underflowed scaled crossing interval')
    single = make_censored_data(uncensored=[1.0e308_dp])
    call close(all_ll(single,2.0_dp,-1.0e308_dp,1.0e308_dp), &
        -4.0_dp-log(2.5_dp)-log(1.0e308_dp),5.0e-13_dp,'overflowing raw displacement')
    single = make_censored_data(uncensored=[0.0_dp])
    call close(all_ll(single,lowk,scale=lowk),0.0_dp,0.0_dp,'smallest shape cusp normalization')
    call close(all_ll(single,highk,scale=highk),-2.0_dp*log(highk), &
        3.0e-13_dp,'largest shape cusp normalization')
    single = make_censored_data(left=[0.0_dp])
    call close(all_ll(single,lowk),2.0_dp*log(lowk),0.0_dp,'smallest shape finite log CDF')
    single = make_censored_data(right=[0.0_dp])
    call close(all_ll(single,highk),-2.0_dp*log(highk),0.0_dp,'largest shape finite log SF')
    single = make_censored_data(left=[0.5_dp])
    call close(all_ll(single,lowk),log(lowk)-log(2.0_dp), &
        3.0e-13_dp,'CDF whose probability underflows')
    single = make_censored_data(right=[-1.0_dp])
    call close(all_ll(single,highk,scale=highk),log(2.0_dp)-2.0_dp*log(highk), &
        5.0e-13_dp,'SF whose probability underflows')
    single = make_censored_data(uncensored=[1.0_dp])
    call close(all_ll(single,1.0e-300_dp,scale=1.0e-300_dp),-1.0_dp, &
        3.0e-13_dp,'compensating small shape and scale')
    single = make_censored_data(uncensored=[huge(1.0_dp)])
    call check(all_ll(single,2.0_dp,scale=tiny(1.0_dp)) == -inf,'unrepresentable log density')
    write(*,'(a,i0,a)') 'Asymmetric Laplace censored: ',checks,' checks'
    if (failures /= 0) error stop 'Asymmetric Laplace censoring failures'
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
        if (abs(value-reference) > tolerance) write(*,*) value,reference
    end subroutine close
end program test_censored_laplace_asymmetric
