! SPDX-License-Identifier: MIT
program test_ecdf_c_api
    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr, c_null_ptr, c_associated
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_empirical_distribution, only : ecdf, ecdf_result, ecdf_confidence_interval
    use scifort_c_api_generated_ecdf
    implicit none
    type(ecdf_result) :: native
    type(ecdf_confidence_interval) :: ci
    type(c_ptr) :: cdf, sf, low, high
    real(c_double) :: events(4), right(2), quantiles(4), probabilities(4), values(4)
    integer(c_int) :: status, count, variance, survival

    events = [1.0_c_double,2.0_c_double,2.0_c_double,4.0_c_double]
    right = [2.0_c_double,5.0_c_double]
    native = ecdf(events,right)
    cdf = c_null_ptr
    sf = c_null_ptr
    call scifort_ecdf_create_f64(4_c_int,2_c_int,events,right,cdf,sf,status)
    if (status /= 0 .or. .not. c_associated(cdf) .or. .not. c_associated(sf)) stop 1
    call scifort_ecdf_info_f64(cdf,count,variance,survival,status)
    if (status /= 0 .or. count /= 4 .or. variance /= 1 .or. survival /= 0) stop 2
    call scifort_ecdf_values_f64(count,cdf,quantiles,probabilities,status)
    if (status /= 0 .or. any(quantiles /= native%cdf%quantiles)) stop 3
    if (any(probabilities /= native%cdf%probabilities)) stop 4
    call scifort_ecdf_evaluate_f64(count,sf,quantiles,values,status)
    if (status /= 0 .or. any(values /= native%sf%evaluate(quantiles))) stop 5
    ci = native%cdf%confidence_interval(0.9_c_double,'log-log')
    call scifort_ecdf_interval_f64(cdf,0.9_c_double,1_c_int,low,high,status)
    if (status /= 0) stop 6
    call scifort_ecdf_destroy_f64(cdf,status)
    if (status /= 0 .or. c_associated(cdf)) stop 7
    call scifort_ecdf_values_f64(count,low,quantiles,probabilities,status)
    if (status /= 0 .or. any(probabilities /= ci%low%probabilities)) stop 8
    call scifort_ecdf_destroy_f64(high,status)
    call scifort_ecdf_interval_f64(low,0.9_c_double,0_c_int,cdf,high,status)
    if (status /= 2 .or. c_associated(cdf) .or. c_associated(high)) stop 9
    call scifort_ecdf_destroy_f64(low,status)
    call scifort_ecdf_destroy_f64(sf,status)
    call scifort_ecdf_evaluate_f64(count,sf,quantiles,values,status)
    if (status /= 5 .or. .not. all(ieee_is_nan(values))) stop 10
    print '(a)', 'ECDF native/C and independent bound checks passed.'
end program test_ecdf_c_api
