! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_specialist
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, nct_censored_loglikelihood, nct_censored_nnlf, nct_fit, &
        landau_censored_loglikelihood, landau_censored_nnlf, landau_fit, &
        nct_logpdf, nct_logcdf, nct_logsf, nct_ppf, &
        landau_logpdf, landau_logcdf, landau_logsf, landau_ppf, &
        t_censored_loglikelihood, normal_logcdf, normal_logsf
    use scifort_math, only : quiet_nan, positive_infinity, log1p_safe, expm1_safe
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: p(4), lower(4), upper(4), x, expected, actual, previous, inf, nan
    real(dp) :: widths(3), candidate, best, df, b
    integer :: family, i, j, k, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.8_dp],interval=reshape([0.3_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 2
        p = [4.0_dp,1.0_dp,0.0_dp,1.0_dp]
        previous = -inf
        do i = -10, 10
            x = real(i,dp)/5.0_dp
            one = make_censored_data(uncensored=[x])
            call close(ll(one,family,p),scalar(x,family,1),3e-11_dp,'ordinary density')
            one = make_censored_data(left=[x])
            actual = ll(one,family,p)
            call close(actual,scalar(x,family,2),3e-11_dp,'ordinary CDF')
            call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotone')
            previous = actual
            one = make_censored_data(right=[x])
            call close(ll(one,family,p),scalar(x,family,3),3e-11_dp,'ordinary SF')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = scalar(x,family,1)+log(nearest(x,1.0_dp)-x)
            call close(ll(one,family,p),expected,3e-11_dp,'adjacent interval')
            one = make_censored_data(interval=reshape([x,x+0.01_dp],[1,2]))
            if (scalar(x,family,2) < log(0.5_dp)) then
                expected = log(exp(scalar(x+0.01_dp,family,2))-exp(scalar(x,family,2)))
            else
                expected = log(exp(scalar(x,family,3))-exp(scalar(x+0.01_dp,family,3)))
            end if
            call close(ll(one,family,p),expected,2e-9_dp,'ordinary interval')
        end do
        do i = 1, 3
            if (family == 1) x = nct_ppf(real(i,dp)/4,4.0_dp,1.0_dp)
            if (family == 2) x = landau_ppf(real(i,dp)/4)
            one = make_censored_data(left=[x])
            call close(ll(one,family,p),log(real(i,dp)/4),3e-10_dp,'CDF/PPF')
        end do
        do k = -1, 1
            p(4) = 10.0_dp**(200*k)
            p(3) = -0.3_dp*p(4)
            x = 0.3_dp*p(4)
            expected = scalar(0.6_dp,family,1)-log(p(4))
            one = make_censored_data(uncensored=[x])
            call close(ll(one,family,p),expected,3e-11_dp,'affine density')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            call close(ll(one,family,p),expected+log(nearest(x,1.0_dp)-x),3e-11_dp,'affine adjacent')
        end do
        p = [4.0_dp,1.0_dp,0.0_dp,1.0_dp]
        call close(ll(empty,family,p),0.0_dp,0.0_dp,'valid empty')
        do k = 1, 4
            if (family == 2 .and. k <= 2) cycle
            lower = p
            lower(k) = nan
            call check(ieee_is_nan(ll(empty,family,lower)),'NaN parameter')
        end do
        lower = p
        lower(4) = 0.0_dp
        call check(ieee_is_nan(ll(empty,family,lower)),'zero scale')
        lower(4) = inf
        call check(ieee_is_nan(ll(empty,family,lower)),'infinite scale')
        one = make_censored_data(left=[nan])
        call check(ieee_is_nan(ll(one,family,p)),'NaN observation')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,p)),'infinite exact invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.6_dp],[1,2]))
        call check(ieee_is_nan(ll(one,family,p)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,p)),'all infinite invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.7_dp],[1,2]))
        call close(ll(one,family,p),scalar(0.7_dp,family,1),3e-11_dp,'equal endpoints exact')
        one = make_censored_data(interval=reshape([-inf,0.7_dp],[1,2]))
        call close(ll(one,family,p),scalar(0.7_dp,family,2),3e-11_dp,'left interval')
        one = make_censored_data(interval=reshape([0.7_dp,inf],[1,2]))
        call close(ll(one,family,p),scalar(0.7_dp,family,3),3e-11_dp,'right interval')
        lower = p
        upper = p
        call fitter(data,family,lower,upper,fit)
        call check(fit%success,'fixed fit success')
        call close(fit%nllf,-ll(data,family,p),3e-12_dp,'fixed NNLF')
        call fitter(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit invalid')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fitter(one,family,lower,upper,fit)
        if (family == 1) call nct_fit([0.2_dp,0.4_dp,0.6_dp],lower,upper,old)
        if (family == 2) call landau_fit([0.2_dp,0.4_dp,0.6_dp],lower(3:4),upper(3:4),old)
        call check(fit%success .and. old%success,'array generic retained')
        call close(fit%nllf,old%nllf,3e-11_dp,'array/exact fit agreement')
        lower(4) = 0.3_dp
        upper(4) = 3.0_dp
        call fitter(data,family,lower,upper,fit)
        call check(fit%success,'scale fit converges')
        call check(fit%nllf <= -ll(data,family,p),'scale fit improves')
        if (family == 1) then
            call close(fit%params(4),0.31430866732103147_dp,3e-5_dp,'SciPy nct scale fit')
            call close(fit%nllf,9.610148768581032_dp,2e-8_dp,'SciPy nct fitted NNLF')
        else
            call close(fit%params(2),0.3_dp,3e-5_dp,'SciPy Landau boundary fit')
            call close(fit%nllf,12.241979042425669_dp,2e-8_dp,'SciPy Landau fitted NNLF')
        end if
        best = inf
        do i = 1, 8
            p(4) = 0.3_dp+2.7_dp*real(i-1,dp)/7
            candidate = -ll(data,family,p)
            best = min(best,candidate)
        end do
        call check(fit%nllf <= best+1e-9_dp,'scale fit beats grid')
        p(4) = 1.0_dp
        expected = -ll(data,family,p)
        if (family == 1) actual = nct_censored_nnlf(data,p(1),p(2))
        if (family == 2) actual = landau_censored_nnlf(data)
        call close(actual,expected,0.0_dp,'NNLF sign')
    end do
    widths = [0.49_dp,100.1_dp,inf]
    do j = 1, 3
        call check(ieee_is_nan(nct_censored_loglikelihood(empty,widths(j),0.0_dp)),'invalid df')
    end do
    call check(ieee_is_nan(nct_censored_loglikelihood(empty,4.0_dp,10.1_dp)),'invalid nc')
    call check(ieee_is_nan(nct_censored_loglikelihood(empty,4.0_dp,-10.1_dp)),'invalid negative nc')
    do j = 1, 3
        df = 0.5_dp
        if (j == 2) df = 4.0_dp
        if (j == 3) df = 100.0_dp
        do i = -2, 2
            x = real(i,dp)*1e300_dp
            one = make_censored_data(right=[x])
            call close(nct_censored_loglikelihood(one,df,0.0_dp), &
                t_censored_loglikelihood(one,df),3e-11_dp,'central t extreme SF')
            one = make_censored_data(uncensored=[x])
            call close(nct_censored_loglikelihood(one,df,0.0_dp), &
                t_censored_loglikelihood(one,df),3e-11_dp,'central t extreme PDF')
        end do
        do k = -1, 1
            one = make_censored_data(left=[0.0_dp])
            call close(nct_censored_loglikelihood(one,df,real(10*k,dp)), &
                normal_logcdf(-real(10*k,dp)),3e-11_dp,'zero CDF')
            one = make_censored_data(right=[0.0_dp])
            call close(nct_censored_loglikelihood(one,df,real(10*k,dp)), &
                normal_logsf(-real(10*k,dp)),3e-11_dp,'zero SF')
        end do
    end do
    ! Public SciPy 1.15.3 values, independently checked by the MP contour driver.
    one = make_censored_data(left=[-2.7_dp])
    call close(landau_censored_loglikelihood(one),-18.614002679863336_dp,3e-13_dp,'Landau switch CDF')
    one = make_censored_data(uncensored=[-2.7_dp])
    call close(landau_censored_loglikelihood(one),-15.344062538208975_dp,3e-13_dp,'Landau switch PDF')
    one = make_censored_data(right=[1e300_dp])
    call close(landau_censored_loglikelihood(one),log(2.0_dp/acos(-1.0_dp))-log(1e300_dp), &
        3e-13_dp,'Landau extreme SF')
    one = make_censored_data(interval=reshape([nearest(0.0_dp,1.0_dp), &
        2*nearest(0.0_dp,1.0_dp)],[1,2]))
    call close(landau_censored_loglikelihood(one),landau_logpdf(0.0_dp) &
        +log(nearest(0.0_dp,1.0_dp)),3e-13_dp,'subnormal interval')
    one = make_censored_data(uncensored=[1e300_dp])
    call check(ieee_is_finite(nct_censored_loglikelihood(one,4.0_dp,1.0_dp,scale=1e-200_dp)), &
        'NCT overflowed coordinate log density')
    call check(ieee_is_finite(landau_censored_loglikelihood(one,scale=1e-200_dp)), &
        'Landau overflowed coordinate log density')
    one = make_censored_data(left=[-1e300_dp])
    call close(landau_censored_loglikelihood(one),-inf,0.0_dp,'unrepresentable Landau log tail')
    ! Independent 100-digit normal/gamma defining-integral regressions.
    one = make_censored_data(right=[100.0_dp])
    call close(nct_censored_loglikelihood(one,10.0_dp,-10.0_dp), &
        -104.52690725142731605_dp,3e-13_dp,'noncentral opposite tail')
    one = make_censored_data(left=[-100.0_dp])
    call close(nct_censored_loglikelihood(one,0.5_dp,10.0_dp), &
        -57.062702884423837229_dp,3e-13_dp,'noncentral low-df opposite tail')
    one = make_censored_data(left=[-1.0_dp])
    call close(nct_censored_loglikelihood(one,100.0_dp,10.0_dp), &
        -63.500164367429839062_dp,3e-13_dp,'noncentral high-df opposite tail')
    one = make_censored_data(right=[100.0_dp])
    call close(nct_censored_loglikelihood(one,10.0_dp,10.0_dp), &
        nct_censored_loglikelihood(make_censored_data(left=[-100.0_dp]),10.0_dp,-10.0_dp), &
        3e-13_dp,'reflection identity')
    one = make_censored_data(interval=reshape([-1e300_dp,-1e299_dp],[1,2]))
    call close(landau_censored_loglikelihood(one),-inf,0.0_dp,'unrepresentable interval log mass')
    do family = 1, 3
        df = 0.5_dp
        if (family == 2) df = 100.0_dp
        do j = -1, 1, 2
            if (family == 3 .and. j < 0) cycle
            x = real(j,dp)*1e300_dp
            do i = 1, 4
                b = x+abs(x)*1e-2_dp
                if (i == 2) b = x+abs(x)*1e-4_dp
                if (i == 3) b = x+abs(x)*1e-8_dp
                if (i == 4) b = nearest(x,1.0_dp)
                one = make_censored_data(interval=reshape([x,b],[1,2]))
                if (family == 3) then
                    expected = log(2.0_dp/acos(-1.0_dp))+log(1e-200_dp) &
                        +log(b-x)-log(x)-log(b)
                    actual = landau_censored_loglikelihood(one,scale=1e-200_dp)
                else
                    p(2) = -10.0_dp
                    if (family == 2) p(2) = 10.0_dp
                    if (j > 0) then
                        expected = nct_censored_loglikelihood(make_censored_data(right=[x]), &
                            df,p(2),scale=1e-200_dp) &
                            +log(-expm1_safe(-df*log1p_safe((b-x)/x)))
                    else
                        expected = nct_censored_loglikelihood(make_censored_data(left=[b]), &
                            df,p(2),scale=1e-200_dp) &
                            +log(-expm1_safe(-df*log1p_safe((b-x)/abs(b))))
                    end if
                    actual = nct_censored_loglikelihood(one,df,p(2),scale=1e-200_dp)
                end if
                call close(actual,expected,3e-11_dp,'overflowed-coordinate power-law interval')
            end do
        end do
    end do
    ! The lower standardized endpoint is finite, but the upper one overflows.
    x = 1.795e300_dp
    b = 1.8e300_dp
    one = make_censored_data(interval=reshape([x,b],[1,2]))
    expected = log(2.0_dp/acos(-1.0_dp))+log(1e-8_dp)+log(b-x)-log(x)-log(b)
    call close(landau_censored_loglikelihood(one,scale=1e-8_dp),expected, &
        3e-11_dp,'Landau partially overflowed interval')
    expected = nct_censored_loglikelihood(make_censored_data(right=[x]),0.5_dp,-10.0_dp, &
        scale=1e-8_dp)+log(-expm1_safe(-0.5_dp*log1p_safe((b-x)/x)))
    call close(nct_censored_loglikelihood(one,0.5_dp,-10.0_dp,scale=1e-8_dp),expected, &
        3e-11_dp,'NCT partially overflowed interval')
    write(*,'(a,i0,a,i0)') 'Specialist censored checks: ',checks,', failures: ',failures
    if (failures /= 0) stop 1
contains
    pure function ll(data,family,p) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! one noncentral t, two Landau
        real(dp), intent(in) :: p(4) !! df, nc, loc, scale; first two ignored for Landau
        real(dp) :: value
        if (family == 1) then
            value = nct_censored_loglikelihood(data,p(1),p(2),p(3),p(4))
        else
            value = landau_censored_loglikelihood(data,p(3),p(4))
        end if
    end function ll
    pure function scalar(x,family,part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standard coordinate
        integer, intent(in) :: family, part !! family and PDF/CDF/SF selectors
        real(dp) :: value
        if (family == 1) then
            select case (part)
            case (1)
                value = nct_logpdf(x,4.0_dp,1.0_dp)
            case (2)
                value = nct_logcdf(x,4.0_dp,1.0_dp)
            case (3)
                value = nct_logsf(x,4.0_dp,1.0_dp)
            end select
        else
            select case (part)
            case (1)
                value = landau_logpdf(x)
            case (2)
                value = landau_logcdf(x)
            case (3)
                value = landau_logsf(x)
            end select
        end if
    end function scalar
    subroutine fitter(data,family,lower,upper,result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        integer, intent(in) :: family !! one noncentral t, two Landau
        real(dp), intent(in) :: lower(4), upper(4) !! df, nc, loc, scale bounds
        type(fit_result), intent(out) :: result !! fit diagnostics
        if (family == 1) then
            call nct_fit(data,lower,upper,result)
        else
            call landau_fit(data,lower(3:4),upper(3:4),result)
        end if
    end subroutine fitter
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! assertion
        character(*), intent(in) :: label !! diagnostic
        checks = checks+1
        if (.not. ok) then
            failures = failures+1
            write(*,'(a)') 'FAIL: '//label
        end if
    end subroutine check
    subroutine close(actual,expected,tol,label)
        implicit none
        real(dp), intent(in) :: actual, expected, tol !! values and scaled log tolerance
        character(*), intent(in) :: label !! diagnostic
        logical :: ok
        ok = .false.
        if (.not. ieee_is_nan(actual) .and. .not. ieee_is_nan(expected)) then
            ok = actual == expected
            if (ieee_is_finite(actual) .and. ieee_is_finite(expected)) &
                ok = abs(actual-expected) <= tol*max(1.0_dp,abs(expected))
        end if
        call check(ok,label)
        if (.not. ok) write(*,*) actual,expected
    end subroutine close
end program test_censored_specialist
