! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_final
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, levy_stable_censored_loglikelihood, levy_stable_censored_nnlf, &
        studentized_range_censored_loglikelihood, studentized_range_censored_nnlf, &
        levy_stable_fit, studentized_range_fit, normal_censored_loglikelihood, &
        cauchy_censored_loglikelihood, levy_censored_loglikelihood, levy_l_censored_loglikelihood, &
        t_censored_loglikelihood, levy_stable_logpdf, levy_stable_logcdf, levy_stable_logsf, &
        studentized_range_logpdf, studentized_range_logcdf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: nan, inf, alpha, beta, x, y, expected, actual, df, k, previous
    real(dp) :: lower(4), upper(4), shapes(8), points(4), loc, scale, z
    integer :: i, j, h, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[.5_dp,1.0_dp,2.0_dp],left=[.8_dp], &
        right=[3.0_dp],interval=reshape([1.2_dp,1.4_dp],[1,2]))
    call close(levy_stable_censored_loglikelihood(data,2.0_dp,.5_dp), &
        normal_censored_loglikelihood(data,scale=sqrt(2.0_dp)),0.0_dp,'stable normal identity')
    call close(levy_stable_censored_loglikelihood(data,1.0_dp,0.0_dp), &
        cauchy_censored_loglikelihood(data),0.0_dp,'stable Cauchy identity')
    call close(levy_stable_censored_loglikelihood(data,.5_dp,1.0_dp), &
        levy_censored_loglikelihood(data),0.0_dp,'stable Levy identity')
    call close(levy_stable_censored_loglikelihood(data,.5_dp,-1.0_dp), &
        levy_l_censored_loglikelihood(data),0.0_dp,'stable reflected Levy identity')
    do j = 1, 3
        alpha = 1.0_dp+.25_dp*real(j,dp)
        beta = .5_dp
        previous = -inf
        do i = -3, 3
            x = real(i,dp)
            one = make_censored_data(uncensored=[x])
            call check(ieee_is_finite(levy_stable_censored_loglikelihood(one,alpha,beta)), &
                'stable finite density')
            one = make_censored_data(left=[x])
            actual = levy_stable_censored_loglikelihood(one,alpha,beta)

            call check(actual >= previous,'stable monotonic CDF')
            previous = actual
            one = make_censored_data(right=[x])
            expected = levy_stable_censored_loglikelihood(one,alpha,beta)
            call close(exp(actual)+exp(expected),1.0_dp,3e-13_dp,'stable tail complement')
        end do
        x = .7_dp
        y = nearest(x,1.0_dp)
        one = make_censored_data(interval=reshape([x,y],[1,2]))
        expected = levy_stable_censored_loglikelihood(make_censored_data(uncensored=[x]),alpha,beta)+log(y-x)
        call close(levy_stable_censored_loglikelihood(one,alpha,beta),expected,3e-10_dp,'stable adjacent')
        call close(levy_stable_censored_loglikelihood(empty,alpha,beta),0.0_dp,0.0_dp,'stable empty')
    end do
    shapes = [1.0_dp,2.0_dp,10.0_dp,30.0_dp,31.0_dp,50.0_dp,75.0_dp,100.0_dp]
    points = [.01_dp,1.0_dp,10.0_dp,1e100_dp]
    do j = 1, size(shapes)
        df = shapes(j)
        do i = 1, 4
            x = points(i)
            one = make_censored_data(uncensored=[x])
            expected = log(2.0_dp)+t_censored_loglikelihood(one,df,scale=sqrt(2.0_dp))
            call close(studentized_range_censored_loglikelihood(one,2.0_dp,df),expected,3e-10_dp,'folded t PDF')
            one = make_censored_data(right=[x])
            expected = log(2.0_dp)+t_censored_loglikelihood(one,df,scale=sqrt(2.0_dp))
            call close(studentized_range_censored_loglikelihood(one,2.0_dp,df),expected,3e-10_dp,'folded t SF')
            one = make_censored_data(interval=reshape([-x,x],[1,2]))
            expected = t_censored_loglikelihood(one,df,scale=sqrt(2.0_dp))
            one = make_censored_data(left=[x])
            call close(studentized_range_censored_loglikelihood(one,2.0_dp,df),expected,3e-10_dp,'folded t CDF')
        end do
    end do
    do j = 1, 2
        k = 2.0_dp+real(j,dp)
        df = 10.0_dp
        do i = 1, 2
            x = real(i,dp)
            one = make_censored_data(uncensored=[x])
            expected = studentized_range_logpdf(x,k,df)
            call close(studentized_range_censored_loglikelihood(one,k,df),expected,3e-8_dp,'range scalar PDF')
            one = make_censored_data(left=[x])
            expected = studentized_range_logcdf(x,k,df)
            call close(studentized_range_censored_loglikelihood(one,k,df),expected,3e-8_dp,'range scalar CDF')
        end do
        x = 1e-100_dp
        expected = .5_dp*log(k)-.5_dp*(k-1)*log(2*acos(-1.0_dp)) &
            +.5_dp*(k-1)*log(2/df)+log_gamma((df+k-1)/2)-log_gamma(df/2)+(k-1)*log(x)
        one = make_censored_data(left=[x])
        call close(studentized_range_censored_loglikelihood(one,k,df),expected,3e-10_dp,'range lower power')
    end do
    ! The expanded finite-df branch retains tails and narrow intervals.
    do j = 5, size(shapes)
        df = shapes(j)
        k = 5.0_dp
        previous = -inf
        do i = 1, 5
            x = real(i,dp)
            one = make_censored_data(left=[x])
            actual = studentized_range_censored_loglikelihood(one,k,df)
            call check(actual >= previous .and. actual <= 0.0_dp,'larger-df CDF monotonic')
            previous = actual
            one = make_censored_data(right=[x])
            expected = studentized_range_censored_loglikelihood(one,k,df)
            call close(exp(actual)+exp(expected),1.0_dp,3e-13_dp,'larger-df complementary tails')
        end do
        x = 2.0_dp
        one = make_censored_data(uncensored=[x])
        expected = studentized_range_censored_loglikelihood(one,k,df)
        y = nearest(x,1.0_dp)
        one = make_censored_data(interval=reshape([x,y],[1,2]))
        call close(studentized_range_censored_loglikelihood(one,k,df),expected+log(y-x), &
            3e-10_dp,'larger-df adjacent interval')
        one = make_censored_data(uncensored=[1e100_dp])
        call check(ieee_is_finite(studentized_range_censored_loglikelihood(one,k,df)), &
            'larger-df remote density')
    end do
    do h = -1, 1, 2
        scale = 10.0_dp**(200*h)
        loc = .3_dp*scale
        x = loc+scale
        z = (x-loc)/scale
        one = make_censored_data(uncensored=[x])
        expected = levy_stable_censored_loglikelihood(make_censored_data(uncensored=[z]),1.5_dp,.5_dp)-log(scale)
        call close(levy_stable_censored_loglikelihood(one,1.5_dp,.5_dp,loc,scale),expected,3e-10_dp,'stable affine')
        expected = studentized_range_censored_loglikelihood(make_censored_data(uncensored=[z]),2.0_dp,10.0_dp) &
            -log(scale)
        call close(studentized_range_censored_loglikelihood(one,2.0_dp,10.0_dp,loc,scale),expected,3e-10_dp,'range affine')
    end do
    one = make_censored_data(interval=reshape([1.0_dp,nearest(1.0_dp,1.0_dp)],[1,2]))
    expected = studentized_range_censored_loglikelihood(make_censored_data(uncensored=[1.0_dp]),3.0_dp,10.0_dp) &
        +log(nearest(1.0_dp,1.0_dp)-1.0_dp)
    call close(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp),expected,3e-10_dp,'range adjacent')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp),-inf,0.0_dp,'range lower density')
    one = make_censored_data(left=[0.0_dp])
    call close(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp),-inf,0.0_dp,'range lower CDF')
    one = make_censored_data(right=[-1.0_dp])
    call close(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp),0.0_dp,0.0_dp,'range lower SF')
    call close(studentized_range_censored_loglikelihood(empty,3.0_dp,inf),0.0_dp,0.0_dp,'range infinite df empty')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(empty,3.0_dp,nan)),'range NaN df')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(empty,nan,10.0_dp)),'range NaN k')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(empty,1.5_dp,10.0_dp)),'range k domain')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(empty,3.0_dp,101.0_dp)),'range df domain')
    call check(ieee_is_nan(levy_stable_censored_loglikelihood(empty,nan,.5_dp)),'stable NaN alpha')
    call check(ieee_is_nan(levy_stable_censored_loglikelihood(empty,1.1_dp,.5_dp)),'stable alpha domain')
    call check(ieee_is_nan(levy_stable_censored_loglikelihood(empty,1.5_dp,.8_dp)),'stable beta domain')
    one = make_censored_data(uncensored=[7.0_dp])
    call check(ieee_is_finite(levy_stable_censored_loglikelihood(one,1.5_dp,.5_dp)),'stable tail density')
    one = make_censored_data(left=[nan])
    call check(ieee_is_nan(levy_stable_censored_loglikelihood(one,1.5_dp,.5_dp)),'stable invalid data')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp)),'range invalid data')
    call close(studentized_range_censored_loglikelihood(empty,3.0_dp,10.0_dp), &
        0.0_dp,0.0_dp,'range valid empty')
    call check(ieee_is_nan(levy_stable_censored_loglikelihood(empty,1.5_dp,.5_dp,scale=0.0_dp)), &
        'stable invalid scale')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(empty,3.0_dp,10.0_dp,scale=nan)), &
        'range NaN scale')
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(empty,3.0_dp,-inf)), &
        'range negative infinite df')
    one = make_censored_data(left=[inf])
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp)), &
        'range infinite data invalid')
    one = make_censored_data(left=[41.0_dp])
    call check(ieee_is_nan(studentized_range_censored_loglikelihood(one,3.0_dp,inf)), &
        'range infinite-df coordinate limit')
    one = make_censored_data(interval=reshape([-inf,1.0_dp],[1,2]))
    call close(levy_stable_censored_loglikelihood(one,1.5_dp,.5_dp), &
        levy_stable_censored_loglikelihood(make_censored_data(left=[1.0_dp]),1.5_dp,.5_dp), &
        0.0_dp,'stable one-sided interval')
    call close(studentized_range_censored_loglikelihood(one,3.0_dp,10.0_dp), &
        studentized_range_censored_loglikelihood(make_censored_data(left=[1.0_dp]),3.0_dp,10.0_dp), &
        0.0_dp,'range one-sided interval')
    one = make_censored_data(uncensored=[1e300_dp])
    expected = studentized_range_censored_loglikelihood(one,2.0_dp,10.0_dp,scale=1e-200_dp)
    call check(ieee_is_finite(expected),'range overflowed coordinate density')
    one = make_censored_data(interval=reshape([1e300_dp,nearest(1e300_dp,1.0_dp)],[1,2]))
    call close(studentized_range_censored_loglikelihood(one,2.0_dp,10.0_dp,scale=1e-200_dp), &
        expected+log(nearest(1e300_dp,1.0_dp)-1e300_dp),3e-10_dp,'range overflowed adjacent')
    call close(studentized_range_censored_nnlf(data,3.0_dp,10.0_dp), &
        -studentized_range_censored_loglikelihood(data,3.0_dp,10.0_dp),0.0_dp,'range NNLF sign')
    call close(levy_stable_censored_nnlf(data,1.5_dp,.5_dp), &
        -levy_stable_censored_loglikelihood(data,1.5_dp,.5_dp),0.0_dp,'stable NNLF sign')
    lower = [1.5_dp,.5_dp,0.0_dp,1.0_dp]
    upper = lower
    call levy_stable_fit(data,lower,upper,fit)
    call check(fit%success,'stable fixed fit')
    call close(fit%nllf,levy_stable_censored_nnlf(data,1.5_dp,.5_dp),3e-10_dp,'stable fit objective')
    call levy_stable_fit(empty,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'stable empty fit')
    call levy_stable_fit([.5_dp,1.0_dp,2.0_dp],lower,upper,old)
    call check(old%success,'stable array overload retained')
    lower = [3.0_dp,10.0_dp,0.0_dp,1.0_dp]
    upper = lower
    call studentized_range_fit(data,lower,upper,fit)
    call check(fit%success,'range fixed fit')
    call close(fit%nllf,studentized_range_censored_nnlf(data,3.0_dp,10.0_dp),3e-10_dp,'range fit objective')
    call studentized_range_fit(empty,lower,upper,fit)
    call check(fit%status == fit_status_invalid_input,'range empty fit')
    call studentized_range_fit([.5_dp,1.0_dp,2.0_dp],lower,upper,old)
    call check(old%success,'range array overload retained')
    lower = [1.0_dp,0.0_dp,0.0_dp,.5_dp]
    upper = [1.0_dp,0.0_dp,0.0_dp,2.5_dp]
    call levy_stable_fit(data,lower,upper,fit)
    call check(fit%success,'stable free scale fit')
    call close(fit%params(4),1.3339167533337817_dp,3e-6_dp,'stable scale reference')
    call close(fit%nllf,12.176776394885268_dp,3e-10_dp,'stable free objective reference')
    lower = [2.0_dp,10.0_dp,0.0_dp,.5_dp]
    upper = [2.0_dp,10.0_dp,0.0_dp,2.0_dp]
    call studentized_range_fit(data,lower,upper,fit)
    call check(fit%success,'range free scale fit')
    call close(fit%params(4),1.2559117716600137_dp,3e-6_dp,'range scale reference')
    call close(fit%nllf,9.2567771817722466_dp,3e-10_dp,'range free objective reference')
    lower = [2.0_dp,100.0_dp,0.0_dp,.5_dp]
    upper = [2.0_dp,100.0_dp,0.0_dp,2.0_dp]
    call studentized_range_fit(data,lower,upper,fit)
    call check(fit%success,'df100 free scale fit')
    call close(fit%params(4),1.3538868383774807_dp,3e-6_dp,'df100 scale reference')
    call close(fit%nllf,9.3210196335503106_dp,3e-10_dp,'df100 free objective reference')
    lower(4) = 1.0_dp
    upper(4) = 1.0_dp
    call studentized_range_fit(data,lower,upper,fit)
    call check(fit%success,'df100 fixed fit')
    call close(fit%nllf,studentized_range_censored_nnlf(data,2.0_dp,100.0_dp), &
        3e-10_dp,'df100 fixed objective')
    ! Tail reflection, complements, switch continuity and physical narrow widths.
    do i = 1, 3
        alpha = 1.0_dp+0.25_dp*real(i,dp)
        do j = -1, 1
            beta = 0.5_dp*real(j,dp)
            previous = 0.0_dp
            do h = 1, 6
                x = 6.0_dp*10.0_dp**real(h-1,dp)
                one = make_censored_data(right=[x])
                actual = levy_stable_censored_loglikelihood(one,alpha,beta)
                call check(ieee_is_finite(actual),'stable remote finite SF')
                call check(actual < previous,'stable remote monotonic SF')
                previous = actual
                expected = levy_stable_censored_loglikelihood( &
                    make_censored_data(left=[-x]),alpha,-beta)
                call close(actual,expected,3e-12_dp,'stable remote reflection')
                expected = levy_stable_censored_loglikelihood( &
                    make_censored_data(left=[x]),alpha,beta)
                call close(exp(actual)+exp(expected),1.0_dp,3e-13_dp,'stable remote complement')
                y = nearest(x,1.0_dp)
                one = make_censored_data(interval=reshape([x,y],[1,2]))
                expected = levy_stable_censored_loglikelihood( &
                    make_censored_data(uncensored=[x]),alpha,beta)+log(y-x)
                call close(levy_stable_censored_loglikelihood(one,alpha,beta), &
                    expected,3e-10_dp,'stable remote adjacent')
            end do
            x = 6.0_dp
            y = nearest(x,1.0_dp)
            actual = levy_stable_censored_loglikelihood(make_censored_data(right=[x]),alpha,beta)
            expected = levy_stable_censored_loglikelihood(make_censored_data(right=[y]),alpha,beta)
            call close(actual,expected,3e-11_dp,'stable tail switch')
        end do
    end do
    x = 1.0e300_dp
    y = nearest(x,1.0_dp)
    scale = 1.0e-200_dp
    one = make_censored_data(uncensored=[x])
    actual = levy_stable_censored_loglikelihood(one,1.5_dp,.5_dp,scale=scale)
    call check(ieee_is_finite(actual),'stable overflowed standardized density')
    one = make_censored_data(interval=reshape([x,y],[1,2]))
    call close(levy_stable_censored_loglikelihood(one,1.5_dp,.5_dp,scale=scale), &
        actual+log(y-x),3e-10_dp,'stable overflowed adjacent')
    data = make_censored_data(uncensored=[.5_dp,1.0_dp,2.0_dp,8.0_dp], &
        left=[-10.0_dp],right=[20.0_dp],interval=reshape([30.0_dp,32.0_dp],[1,2]))
    lower = [1.5_dp,.5_dp,0.0_dp,1.0_dp]
    upper = lower
    call levy_stable_fit(data,lower,upper,fit)
    call check(fit%success,'stable tail fixed fit')
    call close(fit%nllf,levy_stable_censored_nnlf(data,1.5_dp,.5_dp), &
        3e-10_dp,'stable tail fixed objective')
    lower(4) = .5_dp
    upper(4) = 20.0_dp
    call levy_stable_fit(data,lower,upper,fit)
    call check(fit%success,'stable tail free fit')
    call close(fit%params(4),9.92118608975177_dp,3e-6_dp,'stable tail fitted scale')
    call close(fit%nllf,23.868504660682664_dp,3e-10_dp,'stable tail fitted objective')
    write(*,'(a,i0,a,i0)') 'Final censored checks: ',checks,', failures: ',failures
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
        real(dp), intent(in) :: actual, expected, tol !! log values and scaled tolerance
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
end program test_censored_final
