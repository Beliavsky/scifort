! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_censored, only : censored_data, make_censored_data, &
        normal_censored_loglikelihood, normal_censored_nnlf, &
        exponential_censored_loglikelihood, exponential_censored_nnlf, &
        weibull_censored_loglikelihood, weibull_censored_nnlf, &
        lognormal_censored_loglikelihood, lognormal_censored_nnlf, &
        gamma_censored_loglikelihood, gamma_censored_nnlf, &
        beta_censored_loglikelihood, beta_censored_nnlf
    use scifort_censored_bridge, only : censored_normal_loglikelihood
    use scifort_c_api_generated_censored
    implicit none
    real(c_double) :: exact(3), left(1), right(1), lower(2), upper(2), intervals(2,2)
    real(c_double) :: value, expected, empty(0)
    integer(c_int) :: status
    integer :: native_status
    type(censored_data) :: data

    exact = [0.5_dp,1.2_dp,2.5_dp]
    left = [0.3_dp]
    right = [3.5_dp]
    lower = [1.1_dp,4.0_dp]
    upper = [1.4_dp,4.4_dp]
    intervals(:,1) = lower
    intervals(:,2) = upper
    data = make_censored_data(exact,left,right,intervals)
    expected = normal_censored_loglikelihood(data,0.4_dp,1.3_dp)
    call scifort_normal_censored_loglikelihood_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.4_dp,1.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'normal loglikelihood native/C mismatch'
    expected = normal_censored_nnlf(data,0.4_dp,1.3_dp)
    call scifort_normal_censored_nnlf_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.4_dp,1.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'normal nnlf native/C mismatch'
    expected = exponential_censored_loglikelihood(data,0.0_dp,2.3_dp)
    call scifort_exponential_censored_loglikelihood_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'exponential loglikelihood native/C mismatch'
    expected = exponential_censored_nnlf(data,0.0_dp,2.3_dp)
    call scifort_exponential_censored_nnlf_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'exponential nnlf native/C mismatch'
    expected = weibull_censored_loglikelihood(data,1.7_dp,0.0_dp,2.3_dp)
    call scifort_weibull_censored_loglikelihood_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,1.7_dp,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'weibull loglikelihood native/C mismatch'
    expected = weibull_censored_nnlf(data,1.7_dp,0.0_dp,2.3_dp)
    call scifort_weibull_censored_nnlf_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,1.7_dp,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'weibull nnlf native/C mismatch'
    expected = lognormal_censored_loglikelihood(data,0.8_dp,0.0_dp,2.3_dp)
    call scifort_lognormal_censored_loglikelihood_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.8_dp,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'lognormal loglikelihood native/C mismatch'
    expected = lognormal_censored_nnlf(data,0.8_dp,0.0_dp,2.3_dp)
    call scifort_lognormal_censored_nnlf_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.8_dp,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'lognormal nnlf native/C mismatch'
    expected = gamma_censored_loglikelihood(data,0.6_dp,0.0_dp,2.3_dp)
    call scifort_gamma_censored_loglikelihood_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.6_dp,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'gamma loglikelihood native/C mismatch'
    expected = gamma_censored_nnlf(data,0.6_dp,0.0_dp,2.3_dp)
    call scifort_gamma_censored_nnlf_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.6_dp,0.0_dp,2.3_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'gamma nnlf native/C mismatch'
    expected = beta_censored_loglikelihood(data,0.4_dp,0.6_dp,-0.2_dp,6.0_dp)
    call scifort_beta_censored_loglikelihood_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.4_dp,0.6_dp,-0.2_dp,6.0_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'beta loglikelihood native/C mismatch'
    expected = beta_censored_nnlf(data,0.4_dp,0.6_dp,-0.2_dp,6.0_dp)
    call scifort_beta_censored_nnlf_f64(3_c_int,1_c_int,1_c_int,2_c_int, &
        exact,left,right,lower,upper,0.4_dp,0.6_dp,-0.2_dp,6.0_dp,value,status)
    if (status /= 0 .or. value /= expected) error stop 'beta nnlf native/C mismatch'
    call scifort_normal_censored_loglikelihood_f64(0_c_int,0_c_int,0_c_int,0_c_int, &
        empty,empty,empty,empty,empty,0.0_dp,1.0_dp,value,status)
    if (status /= 0 .or. value /= 0.0_dp) error stop 'empty data'
    value = 77.0_dp
    call scifort_normal_censored_loglikelihood_f64(-1_c_int,0_c_int,0_c_int,0_c_int, &
        empty,empty,empty,empty,empty,0.0_dp,1.0_dp,value,status)
    if (status /= 1 .or. value /= 77.0_dp) error stop 'negative extent output preservation'
    call scifort_normal_censored_loglikelihood_f64(0_c_int,0_c_int,0_c_int,0_c_int, &
        empty,empty,empty,empty,empty,0.0_dp,0.0_dp,value,status)
    if (status /= 2 .or. .not. ieee_is_nan(value)) error stop 'invalid scale'
    call censored_normal_loglikelihood(empty,empty,empty,lower,upper(:1), &
        0.0_dp,1.0_dp,value,native_status)
    if (native_status /= 1 .or. .not. ieee_is_nan(value)) error stop 'unpaired native limits'
    print '(a)', 'Censored likelihood native/C checks passed for six families.'
end program test_censored_c_api
