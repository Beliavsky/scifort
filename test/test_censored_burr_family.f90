! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_burr_family
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fisk_censored_loglikelihood, fisk_censored_nnlf, fisk_fit, fisk_logpdf, fisk_logcdf, &
        fisk_logsf, fisk_ppf, fisk_loglikelihood, burr_censored_loglikelihood, &
        burr_censored_nnlf, burr_fit, burr_logpdf, burr_logcdf, burr_logsf, burr_ppf, &
        burr_loglikelihood, burr12_censored_loglikelihood, burr12_censored_nnlf, burr12_fit, &
        burr12_logpdf, burr12_logcdf, burr12_logsf, burr12_ppf, burr12_loglikelihood, &
        lomax_censored_loglikelihood, logistic_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, reflected, transformed
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, c, d, x, actual, expected, a, b, probability, previous
    real(dp) :: reference(3), fitted(3)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    reference = [-14.64543894349949_dp,-14.710130480007274_dp,-14.609935975674242_dp]
    fitted = [1.3590820312500007_dp,0.9417968749999999_dp,0.8731445312499999_dp]
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    do family = 1, 3
        actual = likelihood(data,family,1.7_dp,0.8_dp,scale=1.3_dp)
        call close(actual,reference(family),7.0e-14_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = fisk_censored_nnlf(data,1.7_dp,scale=1.3_dp)
        case (2)
            expected = burr_censored_nnlf(data,1.7_dp,0.8_dp,scale=1.3_dp)
        case (3)
            expected = burr12_censored_nnlf(data,1.7_dp,0.8_dp,scale=1.3_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = 0.4_dp
            if (j == 2) c = 1.0_dp
            if (j == 3) c = 3.0_dp
            d = 0.7_dp
            previous = -inf
            do i = 1, 12
                x = exp(real(i-6,dp))
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = fisk_logpdf(x,c)
                case (2)
                    expected = burr_logpdf(x,c,d)
                case (3)
                    expected = burr12_logpdf(x,c,d)
                end select
                call close(likelihood(single,family,c,d),expected,6.0e-14_dp,'scalar density')
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = fisk_logcdf(x,c)
                case (2)
                    expected = burr_logcdf(x,c,d)
                case (3)
                    expected = burr12_logcdf(x,c,d)
                end select
                actual = likelihood(single,family,c,d)
                call close(actual,expected,6.0e-14_dp,'scalar CDF')
                call check(actual >= previous,'CDF monotonicity')
                call check(actual <= 0.0_dp,'CDF nonnegative probability below one')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = fisk_logsf(x,c)
                case (2)
                    expected = burr_logsf(x,c,d)
                case (3)
                    expected = burr12_logsf(x,c,d)
                end select
                call close(likelihood(single,family,c,d),expected,6.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                probability = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = fisk_ppf(probability,c)
                case (2)
                    x = burr_ppf(probability,c,d)
                case (3)
                    x = burr12_ppf(probability,c,d)
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,c,d)),probability, &
                    2.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        single = make_censored_data(uncensored=[0.0_dp,1.0_dp,3.0_dp])
        select case (family)
        case (1)
            expected = fisk_loglikelihood(single%uncensored,1.0_dp)
            call fisk_fit(single,[1.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
            call fisk_fit(single%uncensored,[1.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],old)
        case (2)
            expected = burr_loglikelihood(single%uncensored,2.0_dp,0.5_dp)
            call burr_fit(single,[2.0_dp,0.5_dp,0.0_dp,1.0_dp], &
                [2.0_dp,0.5_dp,0.0_dp,1.0_dp],fit)
            call burr_fit(single%uncensored,[2.0_dp,0.5_dp,0.0_dp,1.0_dp], &
                [2.0_dp,0.5_dp,0.0_dp,1.0_dp],old)
        case (3)
            expected = burr12_loglikelihood(single%uncensored,1.0_dp,0.8_dp)
            call burr12_fit(single,[1.0_dp,0.8_dp,0.0_dp,1.0_dp], &
                [1.0_dp,0.8_dp,0.0_dp,1.0_dp],fit)
            call burr12_fit(single%uncensored,[1.0_dp,0.8_dp,0.0_dp,1.0_dp], &
                [1.0_dp,0.8_dp,0.0_dp,1.0_dp],old)
        end select
        call check(fit%success .and. old%success,'finite endpoint fits and original array API')
        call close(fit%nllf,-expected,3.0e-14_dp,'array likelihood agreement')
        call close(fit%nllf,old%nllf,3.0e-14_dp,'array fit unchanged')
        call close(likelihood(empty,family,1.0_dp,1.0_dp),0.0_dp,0.0_dp,'valid empty sample')
        call check(ieee_is_nan(likelihood(empty,family,0.0_dp,1.0_dp)),'zero c')
        call check(ieee_is_nan(likelihood(empty,family,nan,1.0_dp)),'NaN c')
        call check(ieee_is_nan(likelihood(empty,family,inf,1.0_dp)),'infinite c')
        if (family /= 1) then
            call check(ieee_is_nan(likelihood(empty,family,1.0_dp,0.0_dp)),'zero d')
            call check(ieee_is_nan(likelihood(empty,family,1.0_dp,nan)),'NaN d')
            call check(ieee_is_nan(likelihood(empty,family,1.0_dp,inf)),'infinite d')
        end if
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,1.0_dp,loc=inf)),'bad location')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,1.0_dp,scale=0.0_dp)),'zero scale')
        call check(ieee_is_nan(likelihood(empty,family,1.0_dp,1.0_dp,scale=nan)),'NaN scale')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family,1.0_dp,1.0_dp)),'NaN observation')
        single = make_censored_data(uncensored=[0.0_dp])
        call check(likelihood(single,family,0.5_dp,1.0_dp) == inf,'singular endpoint')
        call close(likelihood(single,family,1.0_dp,1.0_dp),0.0_dp,0.0_dp,'finite endpoint')
        call check(likelihood(single,family,2.0_dp,1.0_dp) == -inf,'zero endpoint')
        single = make_censored_data(left=[0.0_dp])
        call check(likelihood(single,family,1.0_dp,1.0_dp) == -inf,'zero endpoint CDF')
        single = make_censored_data(right=[-1.0_dp])
        call close(likelihood(single,family,1.0_dp,1.0_dp),0.0_dp,0.0_dp,'unit SF below support')
        single = make_censored_data(interval=reshape([-inf,2.0_dp],[1,2]))
        transformed = make_censored_data(left=[2.0_dp])
        call close(likelihood(single,family,1.7_dp,0.8_dp), &
            likelihood(transformed,family,1.7_dp,0.8_dp),0.0_dp,'infinite lower endpoint')
        single = make_censored_data(interval=reshape([2.0_dp,inf],[1,2]))
        transformed = make_censored_data(right=[2.0_dp])
        call close(likelihood(single,family,1.7_dp,0.8_dp), &
            likelihood(transformed,family,1.7_dp,0.8_dp),0.0_dp,'infinite upper endpoint')
        single = make_censored_data(interval=reshape([1.0_dp,1.0_dp],[1,2]))
        transformed = make_censored_data(uncensored=[1.0_dp])
        call close(likelihood(single,family,1.7_dp,0.8_dp), &
            likelihood(transformed,family,1.7_dp,0.8_dp),0.0_dp,'equal endpoints exact')
        single = make_censored_data(interval=reshape([-1.0_dp,0.0_dp],[1,2]))
        call check(likelihood(single,family,1.7_dp,0.8_dp) == -inf,'impossible interval')
        a = 1.0_dp
        b = nearest(a,1.0_dp)
        single = make_censored_data(interval=reshape([a,b],[1,2]))
        transformed = make_censored_data(uncensored=[a])
        call close(likelihood(single,family,1.7_dp,0.8_dp), &
            likelihood(transformed,family,1.7_dp,0.8_dp)+log(b-a), &
            4.0e-14_dp,'adjacent interval')
    end do
    ! Analytic second-shape MLEs use three exact measurements at fixed c=1, loc=0, scale=1.
    single = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
    expected = 3.0_dp/sum(log(1.0_dp+single%uncensored))
    call burr12_fit(single,[1.0_dp,0.1_dp,0.0_dp,1.0_dp],[1.0_dp,10.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success,'Burr XII analytic MLE converges')
    call close(fit%params(2),expected,4.0e-6_dp,'Burr XII analytic d MLE')
    call burr_fit(single,[1.0_dp,0.1_dp,0.0_dp,1.0_dp],[1.0_dp,10.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%success,'Burr analytic MLE converges')
    call close(fit%params(2),expected,4.0e-6_dp,'Burr analytic d MLE, reciprocal sample')
    call fisk_fit(data,[0.1_dp,0.0_dp,1.3_dp],[10.0_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success,'Fisk mixed fit')
    call close(fit%params(1),fitted(1),3.0e-4_dp,'Fisk SciPy fitted shape')
    call burr_fit(data,[1.7_dp,0.1_dp,0.0_dp,1.3_dp],[1.7_dp,10.0_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success,'Burr mixed fit')
    call close(fit%params(2),fitted(2),3.0e-4_dp,'Burr SciPy fitted shape')
    call burr12_fit(data,[1.7_dp,0.1_dp,0.0_dp,1.3_dp],[1.7_dp,10.0_dp,0.0_dp,1.3_dp],fit)
    call check(fit%success,'Burr XII mixed fit')
    call close(fit%params(2),fitted(3),3.0e-4_dp,'Burr XII SciPy fitted shape')
    call fisk_fit(data,[0.2_dp,-0.1_dp,0.1_dp],[4.0_dp,0.1_dp,4.0_dp],fit)
    call check(fit%success,'Fisk all parameters free')
    call burr_fit(data,[0.8_dp,0.4_dp,-0.1_dp,0.8_dp],[2.0_dp,2.0_dp,0.1_dp,2.0_dp],fit, &
        max_iter=10000,tolerance=1.0e-6_dp)
    call check(fit%success,'Burr all parameters free')
    call burr12_fit(data,[0.2_dp,0.2_dp,-0.1_dp,0.1_dp],[4.0_dp,4.0_dp,0.1_dp,4.0_dp],fit)
    call check(fit%success,'Burr XII all parameters free')
    call fisk_fit(empty,[1.0_dp,0.0_dp,1.0_dp],[2.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'Fisk empty fit rejected')
    call burr_fit(empty,[1.0_dp,1.0_dp,0.0_dp,1.0_dp],[2.0_dp,2.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'Burr empty fit rejected')
    call burr12_fit(empty,[1.0_dp,1.0_dp,0.0_dp,1.0_dp],[2.0_dp,2.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'Burr XII empty fit rejected')
    single = make_censored_data(uncensored=[-1.0_dp])
    call fisk_fit(single,[1.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'Fisk impossible fit')
    call burr_fit(single,[1.0_dp,1.0_dp,0.0_dp,1.0_dp],[1.0_dp,1.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'Burr impossible fit')
    call burr12_fit(single,[1.0_dp,1.0_dp,0.0_dp,1.0_dp],[1.0_dp,1.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'Burr XII impossible fit')
    call fisk_fit(data,[0.1_dp,0.0_dp,0.1_dp],[10.0_dp,0.0_dp,10.0_dp],fit, &
        guess=[1.0_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'Fisk optional controls')
    call burr_fit(data,[0.1_dp,0.1_dp,0.0_dp,0.1_dp],[10.0_dp,10.0_dp,0.0_dp,10.0_dp],fit, &
        guess=[1.0_dp,1.0_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'Burr optional controls')
    call burr12_fit(data,[0.1_dp,0.1_dp,0.0_dp,0.1_dp],[10.0_dp,10.0_dp,0.0_dp,10.0_dp],fit, &
        guess=[1.0_dp,1.0_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
    call check(fit%status == fit_status_max_iter,'Burr XII optional controls')
    transformed = make_censored_data(uncensored=log(data%uncensored/1.3_dp), &
        left=log(data%left/1.3_dp),right=log(data%right/1.3_dp),interval=log(data%interval/1.3_dp))
    call close(fisk_censored_loglikelihood(data,1.7_dp,scale=1.3_dp), &
        logistic_censored_loglikelihood(transformed,scale=1.0_dp/1.7_dp) &
        -sum(log(data%uncensored)),5.0e-14_dp,'Fisk logistic transform and Jacobian')
    call close(burr12_censored_loglikelihood(data,1.0_dp,0.8_dp,scale=1.3_dp), &
        lomax_censored_loglikelihood(data,0.8_dp,scale=1.3_dp),3.0e-14_dp,'Burr XII Lomax identity')
    do j = 1, 3
        c = real(j,dp)/2.0_dp
        call close(fisk_censored_loglikelihood(data,c,scale=1.3_dp), &
            burr_censored_loglikelihood(data,c,1.0_dp,scale=1.3_dp),4.0e-14_dp,'Fisk Burr identity')
        call close(fisk_censored_loglikelihood(data,c,scale=1.3_dp), &
            burr12_censored_loglikelihood(data,c,1.0_dp,scale=1.3_dp),0.0_dp,'Fisk Burr XII identity')
    end do
    reflected = make_censored_data(uncensored=1.0_dp/data%uncensored, &
        left=1.0_dp/data%right,right=1.0_dp/data%left, &
        interval=reshape([1.0_dp/data%interval(:,2),1.0_dp/data%interval(:,1)],[2,2]))
    call close(burr_censored_loglikelihood(data,1.7_dp,0.8_dp), &
        burr12_censored_loglikelihood(reflected,1.7_dp,0.8_dp) &
        -2.0_dp*sum(log(data%uncensored)),5.0e-14_dp,'Burr reciprocal identity and Jacobian')
    single = make_censored_data(interval=reshape([1.0e-100_dp,1.0e100_dp],[1,2]))
    actual = fisk_censored_loglikelihood(single,1.0_dp)
    call close(actual/(-2.0e-100_dp),1.0_dp,5.0e-14_dp,'near-unit interval ratio')
    single = make_censored_data(left=[1.0e-300_dp])
    call close(fisk_censored_loglikelihood(single,2.0_dp),2.0_dp*log(1.0e-300_dp), &
        3.0e-13_dp,'underflowed Fisk CDF')
    single = make_censored_data(right=[1.0e300_dp])
    call close(fisk_censored_loglikelihood(single,2.0_dp),-2.0_dp*log(1.0e300_dp), &
        3.0e-13_dp,'underflowed Fisk SF')
    call close(burr12_censored_loglikelihood(single,1.0e308_dp,1.0e-308_dp), &
        -(1.0e308_dp*1.0e-308_dp)*log(1.0e300_dp),3.0e-12_dp,'compensated overflowing power')
    single = make_censored_data(uncensored=[0.0_dp])
    call close(burr_censored_loglikelihood(single,2.0_dp,0.5_dp),0.0_dp,0.0_dp,'Burr cd=1 endpoint')
    call check(burr_censored_loglikelihood(single,huge(1.0_dp),huge(1.0_dp)) == -inf, &
        'overflowing endpoint shape product')
    write(*,'(a,i0,a)') 'Burr family censored: ',checks,' checks'
    if (failures /= 0) error stop 'Burr family censoring failures'
contains
    pure function likelihood(data, family, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations to evaluate
        integer, intent(in) :: family !! 1 Fisk, 2 Type III, 3 Type XII
        real(dp), intent(in) :: c !! first shape
        real(dp), intent(in) :: d !! second shape, ignored for Fisk
        real(dp), intent(in), optional :: loc !! lower support endpoint
        real(dp), intent(in), optional :: scale !! positive scale
        real(dp) :: value
        select case (family)
        case (1)
            value = fisk_censored_loglikelihood(data,c,loc,scale)
        case (2)
            value = burr_censored_loglikelihood(data,c,d,loc,scale)
        case default
            value = burr12_censored_loglikelihood(data,c,d,loc,scale)
        end select
    end function likelihood
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
end program test_censored_burr_family
