! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_kstwo
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, kstwo_censored_loglikelihood, kstwo_censored_nnlf, &
        kstwo_fit, kstwo_logcdf, kstwo_logsf, kstwo_ppf, uniform_censored_loglikelihood
    use scifort_math, only : quiet_nan, positive_infinity, log1p_safe, expm1_safe
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: n, x, y, expected, actual, previous, inf, nan, p(3), lower(3), upper(3)
    real(dp) :: sizes(6), best, bestn, candidate, scale, loc, knot, gap
    integer :: i, j, k, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    sizes = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,10.0_dp,32.0_dp]
    data = make_censored_data(uncensored=[.3_dp,.4_dp,.6_dp],left=[.25_dp], &
        right=[.8_dp],interval=reshape([.3_dp,.5_dp,.35_dp,.55_dp],[2,2]))
    do j = 1, size(sizes)
        n = sizes(j)
        previous = -inf
        do i = 1, 19
            x = .5_dp/n+(1.0_dp-.5_dp/n)*real(i,dp)/20
            one = make_censored_data(left=[x])
            actual = kstwo_censored_loglikelihood(one,n)
            call close(actual,kstwo_logcdf(x,n),3e-11_dp,'CDF reference')
            call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotone')
            previous = actual
            one = make_censored_data(right=[x])
            call close(kstwo_censored_loglikelihood(one,n),kstwo_logsf(x,n),3e-11_dp,'SF reference')
            one = make_censored_data(uncensored=[x])
            actual = kstwo_censored_loglikelihood(one,n)
            call check(ieee_is_finite(actual),'interior density finite')
            y = nearest(x,1.0_dp)
            one = make_censored_data(interval=reshape([x,y],[1,2]))
            if (x /= 1.0_dp/n) then
                expected = actual+log(y-x)
                call close(kstwo_censored_loglikelihood(one,n),expected,3e-11_dp,'adjacent interval')
            end if
        end do
        do i = 1, 3
            x = kstwo_ppf(real(i,dp)/4,n)
            one = make_censored_data(left=[x])
            call close(kstwo_censored_loglikelihood(one,n),log(real(i,dp)/4),3e-10_dp,'PPF round trip')
        end do
        call close(kstwo_censored_loglikelihood(empty,n),0.0_dp,0.0_dp,'valid empty')
        one = make_censored_data(left=[.5_dp/n])
        call close(kstwo_censored_loglikelihood(one,n),-inf,0.0_dp,'lower CDF')
        one = make_censored_data(right=[1.0_dp])
        call close(kstwo_censored_loglikelihood(one,n),-inf,0.0_dp,'upper SF')
        one = make_censored_data(uncensored=[.5_dp/n,1.0_dp])
        call close(kstwo_censored_loglikelihood(one,n),-inf,0.0_dp,'endpoint density convention')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(kstwo_censored_loglikelihood(one,n)),'all infinite invalid')
        do k = -1, 1
            scale = 10.0_dp**(200*k)
            loc = .3_dp*scale
            x = loc+.7_dp*scale
            one = make_censored_data(uncensored=[x])
            expected = kstwo_censored_loglikelihood(make_censored_data(uncensored=[.7_dp]),n)-log(scale)
            call close(kstwo_censored_loglikelihood(one,n,loc,scale),expected,3e-11_dp,'affine density')
            y = nearest(x,1.0_dp)
            one = make_censored_data(interval=reshape([x,y],[1,2]))
            call close(kstwo_censored_loglikelihood(one,n,loc,scale),expected+log(y-x), &
                3e-11_dp,'affine adjacent')
        end do
        if (n == 1.0_dp) cycle
        knot = 1.0_dp/n
        one = make_censored_data(uncensored=[knot])
        expected = log(2*n)+log_gamma(n+1)-(n-1)*log(n)
        call close(kstwo_censored_loglikelihood(one,n),expected,3e-13_dp,'left density at first knot')
        x = nearest(knot,-1.0_dp)
        one = make_censored_data(interval=reshape([x,knot],[1,2]))
        call close(kstwo_censored_loglikelihood(one,n),expected+log(knot-x),3e-11_dp,'left knot interval')
        one = make_censored_data(interval=reshape([knot,nearest(knot,1.0_dp)],[1,2]))
        expected = expected+log1p_safe(-1.0_dp/n)+log(nearest(knot,1.0_dp)-knot)
        call close(kstwo_censored_loglikelihood(one,n),expected,3e-11_dp,'right knot interval')
        x = nearest(.5_dp/n,1.0_dp)
        one = make_censored_data(left=[x])
        expected = log_gamma(n+1)+n*(log(2.0_dp)+log(x-.5_dp/n))
        call close(kstwo_censored_loglikelihood(one,n),expected,3e-13_dp,'tiny lower CDF')
        one = make_censored_data(interval=reshape([.5_dp/n,x],[1,2]))
        call close(kstwo_censored_loglikelihood(one,n),expected,3e-13_dp,'lower endpoint interval')
        x = nearest(1.0_dp,-1.0_dp)
        one = make_censored_data(right=[x])
        expected = log(2.0_dp)+n*log(1.0_dp-x)
        call close(kstwo_censored_loglikelihood(one,n),expected,3e-13_dp,'tiny upper SF')
        one = make_censored_data(interval=reshape([x,1.0_dp],[1,2]))
        call close(kstwo_censored_loglikelihood(one,n),expected,3e-13_dp,'upper endpoint interval')
    end do
    one = make_censored_data(uncensored=[.6_dp,.7_dp],left=[.8_dp],right=[.9_dp], &
        interval=reshape([.65_dp,.75_dp],[1,2]))
    call close(kstwo_censored_loglikelihood(one,1.0_dp), &
        uniform_censored_loglikelihood(one,.5_dp,.5_dp),3e-13_dp,'n=1 uniform identity')
    p = [4.0_dp,0.0_dp,1.0_dp]
    do k = 1, 3
        lower = p
        lower(k) = nan
        call check(ieee_is_nan(kstwo_censored_loglikelihood(empty,lower(1),lower(2),lower(3))), &
            'NaN parameter')
    end do
    do i = 0, 2
        x = 0.0_dp
        if (i == 1) x = 1.5_dp
        if (i == 2) x = 33.0_dp
        call check(ieee_is_nan(kstwo_censored_loglikelihood(empty,x)),'invalid sample size')
    end do
    call check(ieee_is_nan(kstwo_censored_loglikelihood(empty,inf)),'infinite shape')
    call check(ieee_is_nan(kstwo_censored_loglikelihood(empty,4.0_dp,loc=101.0_dp)),'affine limit')
    call check(ieee_is_nan(kstwo_censored_loglikelihood(empty,4.0_dp,scale=0.0_dp)),'zero scale')
    call check(ieee_is_nan(kstwo_censored_loglikelihood(empty,4.0_dp,scale=1e-320_dp)),'subnormal scale')
    one = make_censored_data(left=[nan])
    call check(ieee_is_nan(kstwo_censored_loglikelihood(one,4.0_dp)),'NaN data')
    one = make_censored_data(interval=reshape([.7_dp,.6_dp],[1,2]))
    call check(ieee_is_nan(kstwo_censored_loglikelihood(one,4.0_dp)),'reversed interval')
    one = make_censored_data(interval=reshape([.4_dp,.4_dp],[1,2]))
    call close(kstwo_censored_loglikelihood(one,4.0_dp), &
        kstwo_censored_loglikelihood(make_censored_data(uncensored=[.4_dp]),4.0_dp),0.0_dp,'point interval')
    lower = p
    upper = p
    call kstwo_fit(data,lower,upper,fit)
    call check(fit%success,'fixed fit')
    call close(fit%nllf,kstwo_censored_nnlf(data,4.0_dp),3e-13_dp,'fit objective')
    call close(kstwo_censored_nnlf(data,4.0_dp),-kstwo_censored_loglikelihood(data,4.0_dp),0.0_dp,'NNLF sign')
    call kstwo_fit(empty,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'empty fit invalid')
    one = make_censored_data(uncensored=[.3_dp,.4_dp,.6_dp])
    call kstwo_fit(one,lower,upper,fit)
    call kstwo_fit([.3_dp,.4_dp,.6_dp],lower,upper,old)
    call check(fit%success .and. old%success,'array overload retained')
    call close(fit%nllf,old%nllf,1e-7_dp,'ordinary array objective agreement')
    lower(1) = 3.0_dp
    upper(1) = 12.0_dp
    call kstwo_fit(data,lower,upper,fit)
    call check(fit%success,'integer fit')
    best = inf
    bestn = 0.0_dp
    do i = 3, 12
        candidate = kstwo_censored_nnlf(data,real(i,dp))
        if (candidate < best) then
            best = candidate
            bestn = real(i,dp)
        end if
    end do
    call close(fit%params(1),bestn,0.0_dp,'integer enumeration optimum')
    call close(fit%nllf,best,3e-11_dp,'integer enumeration objective')
    lower = p
    upper = p
    lower(3) = .9_dp
    upper(3) = 1.5_dp
    call kstwo_fit(data,lower,upper,fit)
    call check(fit%success,'scale fit')
    call close(fit%params(3),1.2_dp,2e-6_dp,'first-knot scale optimum')
    call close(fit%nllf,9.696556989999241_dp,3e-10_dp,'SciPy scale objective')
    call check(fit%nllf <= kstwo_censored_nnlf(data,4.0_dp),'scale fit improves')
    ! Affine first-knot splits must use exactly the same computed endpoint.
    do j = -1, 1, 2
        loc = real(j,dp)*100.0_dp
        n = 10.0_dp
        knot = loc+1.0_dp/n
        x = nearest(knot,-1.0_dp)
        y = nearest(knot,1.0_dp)
        one = make_censored_data(interval=reshape([x,y],[1,2]))
        expected = log(2*n)+log_gamma(n+1)-(n-1)*log(n)
        gap = log((knot-x)+(1-1/n)*(y-knot))
        call close(kstwo_censored_loglikelihood(one,n,loc),expected+gap,3e-11_dp,'affine knot crossing')
    end do
    write(*,'(a,i0,a,i0)') 'Kstwo censored checks: ',checks,', failures: ',failures
    if (failures /= 0) stop 1
contains
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
end program test_censored_kstwo
