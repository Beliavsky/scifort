! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_piecewise
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        bradford_censored_loglikelihood, bradford_censored_nnlf, bradford_fit, &
        bradford_logpdf, bradford_logcdf, bradford_logsf, bradford_ppf, &
        triang_censored_loglikelihood, triang_censored_nnlf, triang_fit, &
        triang_logpdf, triang_logcdf, triang_logsf, triang_ppf, &
        trapezoid_censored_loglikelihood, trapezoid_censored_nnlf, trapezoid_fit, &
        trapezoid_logpdf, trapezoid_logcdf, trapezoid_logsf, trapezoid_ppf, &
        uniform_censored_loglikelihood, beta_censored_loglikelihood, &
        fit_status_invalid_input, fit_status_no_finite_objective, fit_status_max_iter
    use scifort_math, only : quiet_nan, positive_infinity
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, single, empty, other
    type(fit_result) :: fit, old
    real(dp) :: inf, nan, actual, expected, x, c, d, p, previous, reference(3), shapes(3)
    integer :: family, i, j, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    reference = [-10.201408354319534_dp,-13.055640654451915_dp,-12.541072613831128_dp]
    shapes = [1.7_dp,0.3_dp,0.2_dp]
    data = make_censored_data(uncensored=[0.1_dp,0.3_dp,0.6_dp],left=[0.05_dp], &
        right=[0.8_dp],interval=reshape([0.2_dp,0.7_dp,0.25_dp,0.75_dp],[2,2]))
    do family = 1, 3
        c = shapes(family)
        actual = likelihood(data,family,c)
        call close(actual,reference(family),2.0e-13_dp,'mixed mpmath reference')
        select case (family)
        case (1)
            expected = bradford_censored_nnlf(data,c)
        case (2)
            expected = triang_censored_nnlf(data,c)
        case (3)
            expected = trapezoid_censored_nnlf(data,c,0.7_dp)
        end select
        call close(expected,-actual,0.0_dp,'negative likelihood')
        do j = 1, 3
            c = real(j,dp)/5.0_dp
            previous = -inf
            do i = 1, 12
                x = real(i,dp)/13.0_dp
                single = make_censored_data(uncensored=[x])
                select case (family)
                case (1)
                    expected = bradford_logpdf(x,c)
                case (2)
                    expected = triang_logpdf(x,c)
                case (3)
                    expected = trapezoid_logpdf(x,c,0.7_dp)
                end select
                call close(likelihood(single,family,c),expected,5.0e-14_dp,'scalar density')
                single = make_censored_data(left=[x])
                select case (family)
                case (1)
                    expected = bradford_logcdf(x,c)
                case (2)
                    expected = triang_logcdf(x,c)
                case (3)
                    expected = trapezoid_logcdf(x,c,0.7_dp)
                end select
                actual = likelihood(single,family,c)
                call close(actual,expected,5.0e-14_dp,'scalar CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity/range')
                previous = actual
                single = make_censored_data(right=[x])
                select case (family)
                case (1)
                    expected = bradford_logsf(x,c)
                case (2)
                    expected = triang_logsf(x,c)
                case (3)
                    expected = trapezoid_logsf(x,c,0.7_dp)
                end select
                call close(likelihood(single,family,c),expected,5.0e-14_dp,'scalar SF')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = bradford_ppf(p,c)
                case (2)
                    x = triang_ppf(p,c)
                case (3)
                    x = trapezoid_ppf(p,c,0.7_dp)
                end select
                single = make_censored_data(left=[x])
                call close(exp(likelihood(single,family,c)),p,4.0e-14_dp,'CDF/PPF round trip')
            end do
        end do
        c = shapes(family)
        single = make_censored_data(uncensored=[0.1_dp,0.3_dp,0.6_dp])
        select case (family)
        case (1)
            call bradford_fit(single,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],fit)
            call bradford_fit(single%uncensored,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],old)
        case (2)
            call triang_fit(single,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],fit)
            call triang_fit(single%uncensored,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],old)
        case (3)
            call trapezoid_fit(single,[c,0.7_dp,0.0_dp,1.0_dp],[c,0.7_dp,0.0_dp,1.0_dp],fit)
            call trapezoid_fit(single%uncensored,[c,0.7_dp,0.0_dp,1.0_dp],[c,0.7_dp,0.0_dp,1.0_dp],old)
        end select
        call check(fit%success .and. old%success,'fixed fits converge')
        call close(fit%nllf,old%nllf,2.0e-13_dp,'array fit unchanged')
        call close(fit%nllf,-likelihood(single,family,c),0.0_dp,'fit objective')
        call close(likelihood(empty,family,c),0.0_dp,0.0_dp,'empty valid sample')
        call check(ieee_is_nan(likelihood(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(likelihood(empty,family,inf)),'infinite shape')
        call check(ieee_is_nan(likelihood(empty,family,-1.0_dp)),'negative shape')
        call check(ieee_is_nan(likelihood(empty,family,c,loc=nan)),'NaN location')
        call check(ieee_is_nan(likelihood(empty,family,c,loc=inf)),'infinite location')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=nan)),'NaN width')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=inf)),'infinite width')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=0.0_dp)),'zero width')
        call check(ieee_is_nan(likelihood(empty,family,c,scale=-1.0_dp)),'negative width')
        call check(ieee_is_nan(likelihood(empty,family,c,loc=huge(1.0_dp), &
            scale=huge(1.0_dp))),'overflowing upper edge')
        single = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(likelihood(single,family,c)),'NaN observation')
        single = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(likelihood(single,family,c)),'infinite observation')
        single = make_censored_data(interval=reshape([1.0_dp,0.0_dp],[1,2]))
        call check(ieee_is_nan(likelihood(single,family,c)),'reversed interval')
        single = make_censored_data(interval=reshape([-inf,0.5_dp],[1,2]))
        other = make_censored_data(left=[0.5_dp])
        call close(likelihood(single,family,c),likelihood(other,family,c),0.0_dp,'infinite lower')
        single = make_censored_data(interval=reshape([0.5_dp,inf],[1,2]))
        other = make_censored_data(right=[0.5_dp])
        call close(likelihood(single,family,c),likelihood(other,family,c),0.0_dp,'infinite upper')
        single = make_censored_data(interval=reshape([0.5_dp,0.5_dp],[1,2]))
        other = make_censored_data(uncensored=[0.5_dp])
        call close(likelihood(single,family,c),likelihood(other,family,c),0.0_dp,'equal endpoints exact')
        single = make_censored_data(interval=reshape([-1.0_dp,0.5_dp],[1,2]))
        other = make_censored_data(left=[0.5_dp])
        call close(likelihood(single,family,c),likelihood(other,family,c),0.0_dp,'clip lower')
        single = make_censored_data(interval=reshape([0.5_dp,2.0_dp],[1,2]))
        other = make_censored_data(right=[0.5_dp])
        call close(likelihood(single,family,c),likelihood(other,family,c),0.0_dp,'clip upper')
        single = make_censored_data(uncensored=[2.0_dp])
        call check(likelihood(single,family,c) == -inf,'outside support')
        single = make_censored_data(left=[0.0_dp])
        call check(likelihood(single,family,c) == -inf,'CDF zero at lower edge')
        single = make_censored_data(right=[1.0_dp])
        call check(likelihood(single,family,c) == -inf,'SF zero at upper edge')
        single = make_censored_data(interval=reshape([-1.0_dp,2.0_dp],[1,2]))
        call close(likelihood(single,family,c),0.0_dp,0.0_dp,'whole support interval')
        select case (family)
        case (1)
            call bradford_fit(data,[0.1_dp,-0.05_dp,0.9_dp],[10.0_dp,0.0_dp,1.1_dp],fit, &
                max_iter=10000,tolerance=1.0e-7_dp)
        case (2)
            call triang_fit(data,[0.1_dp,-0.05_dp,0.9_dp],[0.9_dp,0.0_dp,1.1_dp],fit, &
                max_iter=10000,tolerance=1.0e-7_dp)
        case (3)
            call trapezoid_fit(data,[0.1_dp,0.5_dp,-0.05_dp,0.9_dp],[0.4_dp,0.9_dp,0.0_dp,1.1_dp],fit, &
                max_iter=10000,tolerance=1.0e-7_dp)
        end select
        call check(fit%success,'all-free bounded fit')
        if (family == 3) then
            expected = likelihood(data,family,fit%params(1),fit%params(3), &
                fit%params(4),fit%params(2))
        else
            expected = likelihood(data,family,fit%params(1),fit%params(2),fit%params(3))
        end if
        call close(fit%nllf,-expected,0.0_dp,'free fit objective')
        select case (family)
        case (1)
            call bradford_fit(data,[0.01_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit)
            call close(fit%params(1),2.8134765625000036_dp,2.0e-4_dp,'SciPy Bradford shape')
        case (2)
            call triang_fit(data,[0.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
            call close(fit%params(1),4.8828124999111865e-05_dp,2.0e-4_dp,'SciPy triangular mode')
        case (3)
            call trapezoid_fit(data,[0.0_dp,0.7_dp,0.0_dp,1.0_dp], &
                [0.7_dp,0.7_dp,0.0_dp,1.0_dp],fit)
            call close(fit%params(1),0.0_dp,2.0e-4_dp,'SciPy fixed-d trapezoid fit')
        end select
        call check(fit%success,'partly fixed shape fit')
        select case (family)
        case (1)
            call bradford_fit(data,[0.01_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit, &
                guess=[1.7_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
        case (2)
            call triang_fit(data,[0.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit, &
                guess=[0.3_dp,0.0_dp,1.0_dp],max_iter=1,tolerance=1.0e-10_dp)
        case (3)
            call trapezoid_fit(data,[0.0_dp,0.7_dp,0.0_dp,1.0_dp], &
                [0.7_dp,0.7_dp,0.0_dp,1.0_dp],fit,guess=[0.2_dp,0.7_dp,0.0_dp,1.0_dp], &
                max_iter=1,tolerance=1.0e-10_dp)
        end select
        call check(fit%status == fit_status_max_iter,'fit controls forwarded')
        select case (family)
        case (1)
            call bradford_fit(empty,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],fit)
        case (2)
            call triang_fit(empty,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],fit)
        case (3)
            call trapezoid_fit(empty,[c,0.7_dp,0.0_dp,1.0_dp],[c,0.7_dp,0.0_dp,1.0_dp],fit)
        end select
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        single = make_censored_data(uncensored=[2.0_dp])
        select case (family)
        case (1)
            call bradford_fit(single,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],fit)
        case (2)
            call triang_fit(single,[c,0.0_dp,1.0_dp],[c,0.0_dp,1.0_dp],fit)
        case (3)
            call trapezoid_fit(single,[c,0.7_dp,0.0_dp,1.0_dp],[c,0.7_dp,0.0_dp,1.0_dp],fit)
        end select
        call check(fit%status == fit_status_no_finite_objective,'impossible fixed objective')
    end do
    call check(ieee_is_nan(bradford_censored_loglikelihood(empty,0.0_dp)),'Bradford zero shape')
    call check(ieee_is_nan(triang_censored_loglikelihood(empty,1.1_dp)),'mode above one')
    call check(ieee_is_nan(trapezoid_censored_loglikelihood(empty,0.8_dp,0.2_dp)), &
        'reversed plateau')
    call check(ieee_is_nan(trapezoid_censored_loglikelihood(empty,0.2_dp,nan)), &
        'NaN second shape')
    call close(trapezoid_censored_loglikelihood(data,0.0_dp,1.0_dp), &
        uniform_censored_loglikelihood(data),4.0e-14_dp,'uniform trapezoid identity')
    call close(triang_censored_loglikelihood(data,0.3_dp), &
        trapezoid_censored_loglikelihood(data,0.3_dp,0.3_dp),0.0_dp,'triangular identity')
    call close(triang_censored_loglikelihood(data,1.0_dp), &
        beta_censored_loglikelihood(data,2.0_dp,1.0_dp),4.0e-14_dp,'beta rising triangle identity')
    call close(triang_censored_loglikelihood(data,0.0_dp), &
        beta_censored_loglikelihood(data,1.0_dp,2.0_dp),4.0e-14_dp,'beta falling triangle identity')
    call close(bradford_censored_loglikelihood(data,nearest(0.0_dp,1.0_dp)), &
        uniform_censored_loglikelihood(data),4.0e-14_dp,'smallest Bradford shape limit')
    single = make_censored_data(uncensored=[0.0_dp])
    call close(triang_censored_loglikelihood(single,0.0_dp),log(2.0_dp),0.0_dp,'left endpoint mode')
    call check(triang_censored_loglikelihood(single,1.0_dp) == -inf,'zero rising endpoint')
    single = make_censored_data(uncensored=[1.0_dp])
    call close(triang_censored_loglikelihood(single,1.0_dp),log(2.0_dp),0.0_dp,'right endpoint mode')
    single = make_censored_data(right=[1.0e300_dp])
    expected = log(2.0_dp)-600.0_dp*log(10.0_dp)
    call close(triang_censored_loglikelihood(single,1.0_dp,1.0e-300_dp,1.0e300_dp), &
        expected,4.0e-13_dp,'underflowed compensated upper distance')
    single = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    expected = 2.0_dp*log(nearest(0.0_dp,1.0_dp))-log(0.3_dp)
    call close(triang_censored_loglikelihood(single,0.3_dp),expected,4.0e-13_dp, &
        'smallest lower interval width')
    single = make_censored_data(interval=reshape([2.0e-100_dp,1.0_dp],[1,2]))
    actual = triang_censored_loglikelihood(single,0.0_dp,loc=1.0e-100_dp)
    call close(actual/(-2.0e-100_dp),1.0_dp,3.0e-13_dp,'nearly full triangular interval ratio')
    actual = bradford_censored_loglikelihood(single,1.7_dp,loc=1.0e-100_dp)
    expected = -1.7_dp*(1.0_dp+1.0_dp/2.7_dp)*1.0e-100_dp/log(2.7_dp)
    call close(actual/expected,1.0_dp,3.0e-13_dp,'nearly full Bradford interval ratio')
    single = make_censored_data(interval=reshape([1.0e-24_dp,6.0e-24_dp],[1,2]))
    actual = trapezoid_censored_loglikelihood(single,nearest(0.0_dp,1.0_dp), &
        nearest(0.0_dp,1.0_dp),scale=1.0e300_dp)
    call check(ieee_is_finite(actual),'interval crosses subnormal breakpoint')
    other = make_censored_data(left=[6.0e-24_dp])
    ! F(z)=2*z-c+O(z*z) here; retain the subnormal probability in log form.
    expected = log(1.2e-23_dp)-log(1.0e300_dp) &
        +log(1.0_dp-exp(log(nearest(0.0_dp,1.0_dp))+log(1.0e300_dp)-log(1.2e-23_dp)))
    call close(trapezoid_censored_loglikelihood(other,nearest(0.0_dp,1.0_dp), &
        nearest(0.0_dp,1.0_dp),scale=1.0e300_dp),expected,4.0e-13_dp, &
        'CDF beyond subnormal breakpoint')
    write(*,'(a,i0,a)') 'Piecewise censored: ',checks,' checks'
    if (failures /= 0) error stop 'Piecewise censoring failures'
contains
    pure function likelihood(data, family, c, loc, scale, d) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! 1 Bradford, 2 triangular, 3 trapezoidal
        real(dp), intent(in) :: c !! first shape in the family's domain
        real(dp), intent(in), optional :: loc !! finite lower endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive width, default one
        real(dp), intent(in), optional :: d !! trapezoid second shape, default 0.7
        real(dp) :: value, second
        second = 0.7_dp
        if (present(d)) second = d
        select case (family)
        case (1)
            value = bradford_censored_loglikelihood(data,c,loc,scale)
        case (2)
            value = triang_censored_loglikelihood(data,c,loc,scale)
        case (3)
            value = trapezoid_censored_loglikelihood(data,c,second,loc,scale)
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
end program test_censored_piecewise
