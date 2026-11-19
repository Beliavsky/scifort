! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_studentized_range_large_k
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        studentized_range_censored_loglikelihood, studentized_range_censored_nnlf, &
        studentized_range_fit
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, empty
    type(fit_result) :: fit
    real(dp) :: k, df, x, y, actual, expected, previous, lc, ls, c, m, inf, nan
    real(dp) :: scale, loc, z, lower(4), upper(4)
    real(dp), parameter :: counts(3) = [11.0_dp,37.5_dp,100.0_dp]
    real(dp), parameter :: dfs(3) = [1.0_dp,100.0_dp,100000.0_dp]
    real(dp), parameter :: points(4) = [.1_dp,1.0_dp,4.0_dp,8.0_dp]
    integer :: i, j, h, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    do j = 1, size(counts)
        k = counts(j)
        ! Representative pairs; the Python reference grid covers their cross-product.
        block
            df = dfs(j)
            previous = -inf
            do i = 1, size(points)
                x = points(i)
                data = make_censored_data(uncensored=[x])
                actual = studentized_range_censored_loglikelihood(data,k,df)
                call check(ieee_is_finite(actual),'large-k finite density')
                data = make_censored_data(left=[x])
                lc = studentized_range_censored_loglikelihood(data,k,df)
                call check(ieee_is_finite(lc) .and. lc >= previous .and. lc <= 0.0_dp, &
                    'large-k monotonic CDF')
                previous = lc
                data = make_censored_data(right=[x])
                ls = studentized_range_censored_loglikelihood(data,k,df)
                call check(ieee_is_finite(ls) .and. ls <= 0.0_dp,'large-k finite SF')
                call close(exp(lc)+exp(ls),1.0_dp,3e-13_dp,'large-k complementary tails')
            end do
            x = 4.0_dp
            y = nearest(x,1.0_dp)
            expected = studentized_range_censored_loglikelihood( &
                make_censored_data(uncensored=[x]),k,df)+log(y-x)
            data = make_censored_data(interval=reshape([x,y],[1,2]))
            call close(studentized_range_censored_loglikelihood(data,k,df),expected, &
                3e-10_dp,'large-k adjacent interval')
            x = 1.0e-100_dp
            m = k-1.0_dp
            c = log(k)/2-m*log(2*acos(-1.0_dp))/2
            if (df < 100000.0_dp) c = c+m/2*log(2/df) &
                +log_gamma((df+m)/2)-log_gamma(df/2)
            data = make_censored_data(left=[x])
            call close(studentized_range_censored_loglikelihood(data,k,df),c+m*log(x), &
                3e-12_dp,'large-k lower CDF power')
            data = make_censored_data(uncensored=[x])
            call close(studentized_range_censored_loglikelihood(data,k,df), &
                c+log(m)+(m-1)*log(x),3e-12_dp,'large-k lower PDF power')
            call close(studentized_range_censored_loglikelihood(empty,k,df), &
                0.0_dp,0.0_dp,'large-k valid empty')
        end block
    end do
    write(*,'(a)') 'Large-k central, lower-tail and interval checks passed'
    k = 100.0_dp
    df = 100.0_dp
    m = k-1
    x = 1.0e-100_dp
    c = log(k)/2-m*log(2*acos(-1.0_dp))/2+m/2*log(2/df) &
        +log_gamma((df+m)/2)-log_gamma(df/2)
    data = make_censored_data(left=[x])
    call close(studentized_range_censored_loglikelihood(data,k,df),c+m*log(x), &
        3e-12_dp,'k100 df100 tilted moment bound')
    data = make_censored_data(right=[12.0_dp])
    call close(studentized_range_censored_loglikelihood(data,k,inf), &
        -29.87048776184788_dp,3e-11_dp,'k100 small-ratio SF hazard')
    ! Affine agreement includes normal scales spanning 400 decimal orders.
    do h = -1, 1, 2
        scale = 10.0_dp**(200*h)
        loc = .3_dp*scale
        x = loc+4*scale
        z = (x-loc)/scale
        data = make_censored_data(uncensored=[x])
        expected = studentized_range_censored_loglikelihood( &
            make_censored_data(uncensored=[z]),100.0_dp,100.0_dp)-log(scale)
        call close(studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp,loc,scale), &
            expected,3e-10_dp,'large-k affine density')
    end do
    x = 1.0e300_dp
    y = nearest(x,1.0_dp)
    scale = 1.0e-200_dp
    data = make_censored_data(uncensored=[x])
    actual = studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp,scale=scale)
    call check(ieee_is_finite(actual),'large-k overflowed density coordinate')
    data = make_censored_data(right=[x])
    call check(ieee_is_finite(studentized_range_censored_loglikelihood( &
        data,100.0_dp,100.0_dp,scale=scale)),'large-k overflowed SF coordinate')
    data = make_censored_data(interval=reshape([x,y],[1,2]))
    call close(studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp,scale=scale), &
        actual+log(y-x),3e-10_dp,'large-k overflowed adjacent interval')
    data = make_censored_data(uncensored=[0.0_dp])
    call close(studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp), &
        -inf,0.0_dp,'large-k endpoint density')
    data = make_censored_data(left=[0.0_dp])
    call close(studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp), &
        -inf,0.0_dp,'large-k endpoint CDF')
    data = make_censored_data(right=[0.0_dp])
    call close(studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp), &
        0.0_dp,0.0_dp,'large-k endpoint SF')
    data = make_censored_data(interval=reshape([-inf,4.0_dp],[1,2]))
    expected = studentized_range_censored_loglikelihood( &
        make_censored_data(left=[4.0_dp]),100.0_dp,100.0_dp)
    call close(studentized_range_censored_loglikelihood(data,100.0_dp,100.0_dp), &
        expected,0.0_dp,'large-k one-sided interval')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood( &
        empty,nearest(100.0_dp,1.0_dp),100.0_dp)),'large-k upper shape limit')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood( &
        empty,nan,100.0_dp)),'large-k NaN shape')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood( &
        empty,inf,100.0_dp)),'large-k infinite shape')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood( &
        empty,100.0_dp,101.0_dp)),'large-k df gap retained')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood( &
        empty,100.0_dp,100.0_dp,scale=0.0_dp)),'large-k invalid scale')
    data = make_censored_data(uncensored=[3.0_dp,4.0_dp,5.0_dp], &
        left=[2.5_dp],right=[7.0_dp],interval=reshape([4.5_dp,4.75_dp],[1,2]))
    write(*,'(a)') 'Large-k affine, endpoint and invalid-input checks passed'
    lower = [100.0_dp,100.0_dp,0.0_dp,1.0_dp]
    upper = lower
    call studentized_range_fit(data,lower,upper,fit)
    call check(fit%success,'large-k fixed fit')
    call close(fit%nllf,studentized_range_censored_nnlf(data,100.0_dp,100.0_dp), &
        3e-10_dp,'large-k fixed objective')
    lower(4) = .5_dp
    upper(4) = 2.0_dp
    call studentized_range_fit(data,lower,upper,fit)
    call check(fit%success,'large-k free-scale fit')
    call close(fit%params(4),.800935178253828_dp,3e-6_dp,'large-k scale reference')
    call close(fit%nllf,25.123232284243887_dp,3e-10_dp,'large-k objective reference')
    write(*,'(a,i0,a,i0)') 'Large-k censored checks: ',checks,', failures: ',failures
    if (failures /= 0) stop 1
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! assertion result
        character(*), intent(in) :: label !! diagnostic label
        checks = checks+1
        if (.not. ok) then
            failures = failures+1
            write(*,'(a)') 'FAIL: '//label
        end if
    end subroutine check
    subroutine close(actual,expected,tol,label)
        implicit none
        real(dp), intent(in) :: actual, expected !! log values or fit parameters
        real(dp), intent(in) :: tol !! tolerance scaled by max(1,abs(expected))
        character(*), intent(in) :: label !! diagnostic label
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
end program test_censored_studentized_range_large_k
