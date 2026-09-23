! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_merge_regressions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, expm1_safe
    use scifort_special_reductions, only : log_softmax
    use scifort_hypergeom, only : hypergeom_cdf, hypergeom_sf, hypergeom_logcdf, hypergeom_logsf
    use scifort_randint, only : randint_cdf, randint_sf, randint_logcdf, randint_logsf, &
        randint_ppf, randint_isf
    use scifort_zipf, only : zipf_cdf, zipf_sf, zipf_logcdf, zipf_logsf
    use scifort_normal, only : normal_logcdf, normal_logsf
    use scifort_multivariate_normal, only : multivariate_normal_logcdf
    use test_checks, only : check_close, check_true
    implicit none
    integer :: failures, i
    real(dp) :: pinf, ninf, nan, k, want
    real(dp) :: limits(4), values(2), baseline(2), covariance(2,2)
    failures = 0
    pinf = positive_infinity(0.0_dp)
    ninf = negative_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)

    call check_close('zipf cdf +inf', zipf_cdf(pinf, 2.0_dp), 1.0_dp, &
        0.0_dp, 0.0_dp, failures)
    call check_close('zipf sf +inf', zipf_sf(pinf, 2.0_dp), 0.0_dp, &
        0.0_dp, 0.0_dp, failures)
    call check_close('zipf logcdf +inf', zipf_logcdf(pinf, 2.0_dp), 0.0_dp, &
        0.0_dp, 0.0_dp, failures)
    call check_close('zipf logsf +inf', zipf_logsf(pinf, 2.0_dp), ninf, &
        0.0_dp, 0.0_dp, failures)
    call check_close('zipf cdf -inf', zipf_cdf(ninf, 2.0_dp), 0.0_dp, &
        0.0_dp, 0.0_dp, failures)
    call check_close('zipf sf -inf', zipf_sf(ninf, 2.0_dp), 1.0_dp, &
        0.0_dp, 0.0_dp, failures)
    call check_true('zipf cdf NaN', ieee_is_nan(zipf_cdf(nan, 2.0_dp)), failures)
    call check_true('zipf sf NaN', ieee_is_nan(zipf_sf(nan, 2.0_dp)), failures)
    call check_true('zipf logcdf NaN', ieee_is_nan(zipf_logcdf(nan, 2.0_dp)), failures)
    call check_true('zipf logsf NaN', ieee_is_nan(zipf_logsf(nan, 2.0_dp)), failures)
    call check_true('zipf invalid shape at infinity', &
        ieee_is_nan(zipf_cdf(pinf, 1.0_dp)), failures)
    call check_close('zipf fractional shifted count', zipf_cdf(3.25_dp, 2.0_dp, 0.5_dp), &
        zipf_cdf(2, 2.0_dp), 0.0_dp, 0.0_dp, failures)

    limits = [real(huge(0), dp) + 1.0_dp, 1.0e20_dp, huge(1.0_dp), pinf]
    do i = 1, size(limits)
        k = limits(i)
        call check_close('hypergeom cdf large', hypergeom_cdf(k, 13.0_dp, 8.0_dp, 7.0_dp), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('hypergeom sf large', hypergeom_sf(k, 13.0_dp, 8.0_dp, 7.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('hypergeom logcdf large', &
            hypergeom_logcdf(k, 13.0_dp, 8.0_dp, 7.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('hypergeom logsf large', &
            hypergeom_logsf(k, 13.0_dp, 8.0_dp, 7.0_dp), ninf, 0.0_dp, 0.0_dp, failures)
        call check_close('hypergeom cdf negative large', &
            hypergeom_cdf(-k, 13.0_dp, 8.0_dp, 7.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('randint cdf large', randint_cdf(k, 0.0_dp, 10.0_dp), &
            1.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('randint sf large', randint_sf(k, 0.0_dp, 10.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('randint logcdf large', randint_logcdf(k, 0.0_dp, 10.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('randint logsf large', randint_logsf(k, 0.0_dp, 10.0_dp), &
            ninf, 0.0_dp, 0.0_dp, failures)
        call check_close('randint cdf negative large', randint_cdf(-k, 0.0_dp, 10.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
    end do
    call check_close('randint negative fractional floor', randint_cdf(-0.5_dp, -2.0_dp, 2.0_dp), &
        0.5_dp, 0.0_dp, 0.0_dp, failures)
    call check_close('randint wide ppf', randint_ppf(0.75_dp, 0.0_dp, 4.0e9_dp), &
        2999999999.0_dp, 0.0_dp, 0.0_dp, failures)
    call check_close('randint wide isf', randint_isf(0.75_dp, 0.0_dp, 4.0e9_dp), &
        999999999.0_dp, 0.0_dp, 0.0_dp, failures)
    call check_true('hypergeom NaN', &
        ieee_is_nan(hypergeom_cdf(nan, 13.0_dp, 8.0_dp, 7.0_dp)), failures)
    call check_true('randint NaN', ieee_is_nan(randint_cdf(nan, 0.0_dp, 10.0_dp)), failures)

    values = log_softmax([1.0e20_dp, 1.0e20_dp])
    do i = 1, 2
        call check_close('log_softmax large equal', values(i), -log(2.0_dp), &
            2.0e-16_dp, 0.0_dp, failures)
    end do
    baseline = log_softmax([0.0_dp, -1.0_dp])
    values = log_softmax([1.0e12_dp, 1.0e12_dp - 1.0_dp])
    do i = 1, 2
        call check_close('log_softmax shift invariance', values(i), baseline(i), &
            2.0e-16_dp, 0.0_dp, failures)
    end do
    values = log_softmax([-1.0e20_dp, -1.0e20_dp])
    call check_close('log_softmax normalization', sum(exp(values)), 1.0_dp, &
        2.0e-16_dp, 0.0_dp, failures)

    want = normal_logcdf(-40.0_dp)
    call check_close('mvn scalar log tail', multivariate_normal_logcdf([-40.0_dp]), &
        want, 1.0e-12_dp, 0.0_dp, failures)
    call check_close('mvn independent log tails', multivariate_normal_logcdf([-40.0_dp, -40.0_dp]), &
        2.0_dp * want, 2.0e-12_dp, 0.0_dp, failures)
    covariance = reshape([4.0_dp, 0.0_dp, 123.0_dp, 9.0_dp], [2,2])
    call check_close('mvn diagonal scale and lower triangle', &
        multivariate_normal_logcdf([-79.0_dp, -118.0_dp], mean=[1.0_dp, 2.0_dp], cov=covariance), &
        2.0_dp * want, 2.0e-12_dp, 0.0_dp, failures)
    want = normal_logsf(40.0_dp) + log(-expm1_safe(normal_logsf(41.0_dp) - normal_logsf(40.0_dp)))
    call check_close('mvn positive tail interval', &
        multivariate_normal_logcdf([41.0_dp], lower_limit=[40.0_dp]), &
        want, 1.0e-12_dp, 0.0_dp, failures)
    call check_true('mvn negative signed interval', &
        ieee_is_nan(multivariate_normal_logcdf([40.0_dp], lower_limit=[41.0_dp])), failures)
    call check_close('mvn two reversed intervals', &
        multivariate_normal_logcdf([40.0_dp, 40.0_dp], lower_limit=[41.0_dp, 41.0_dp]), &
        2.0_dp * want, 2.0e-12_dp, 0.0_dp, failures)
    call check_close('mvn empty interval', multivariate_normal_logcdf([1.0_dp], lower_limit=[1.0_dp]), &
        ninf, 0.0_dp, 0.0_dp, failures)
    call check_true('mvn invalid variance', &
        ieee_is_nan(multivariate_normal_logcdf([0.0_dp], cov=reshape([-1.0_dp],[1,1]))), failures)
    call check_true('mvn invalid tolerance', &
        ieee_is_nan(multivariate_normal_logcdf([0.0_dp], abseps=-1.0_dp)), failures)
    call check_close('mvn log beyond binary64 range', &
        multivariate_normal_logcdf([-1.0e200_dp]), ninf, 0.0_dp, 0.0_dp, failures)
    call check_close('mvn whole real line', multivariate_normal_logcdf([pinf]), &
        0.0_dp, 0.0_dp, 0.0_dp, failures)
    covariance = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2,2])
    call check_close('mvn deterministic coordinate inside box', &
        multivariate_normal_logcdf([0.0_dp, -40.0_dp], cov=covariance, allow_singular=.true.), &
        normal_logcdf(-40.0_dp), 1.0e-12_dp, 0.0_dp, failures)
    call check_close('mvn deterministic coordinate outside box', &
        multivariate_normal_logcdf([-1.0_dp, -40.0_dp], cov=covariance, allow_singular=.true.), &
        ninf, 0.0_dp, 0.0_dp, failures)
    call check_true('mvn singular rejected by default', &
        ieee_is_nan(multivariate_normal_logcdf([0.0_dp, 0.0_dp], cov=covariance)), failures)
    call check_true('mvn NaN observation', ieee_is_nan(multivariate_normal_logcdf([nan])), failures)
    if (failures /= 0) then
        print '(a,i0)', 'merge regression failures: ', failures
        error stop 1
    end if
    print '(a)', 'test_merge_regressions: PASS'
end program test_merge_regressions
