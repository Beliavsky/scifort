! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_multivariate_fit_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_multivariate_normal, only : multivariate_normal_fit
    use scifort_c_api_generated_multivariate_normal, only : scifort_multivariate_normal_fit_f64
    implicit none
    real(c_double) :: data(2,3), fixed_mean(2), fixed_cov(2,2)
    real(c_double) :: expected_mean(2), expected_cov(2,2), actual_mean(2), actual_cov(2,2)
    integer(c_int) :: status, use_mean, use_cov
    integer :: native_status
    data=reshape([1.0_c_double,2.0_c_double,3.0_c_double,6.0_c_double, &
        5.0_c_double,4.0_c_double],[2,3])
    fixed_mean=[2.0_c_double,3.0_c_double]
    fixed_cov=reshape([2.0_c_double,.2_c_double,.2_c_double,3.0_c_double],[2,2])
    do use_mean=0,1
        do use_cov=0,1
            if (use_mean==1 .and. use_cov==1) then
                call multivariate_normal_fit(data,expected_mean,expected_cov,native_status, &
                    fixed_mean,fixed_cov)
            else if (use_mean==1) then
                call multivariate_normal_fit(data,expected_mean,expected_cov,native_status, &
                    fixed_mean=fixed_mean)
            else if (use_cov==1) then
                call multivariate_normal_fit(data,expected_mean,expected_cov,native_status, &
                    fixed_cov=fixed_cov)
            else
                call multivariate_normal_fit(data,expected_mean,expected_cov,native_status)
            end if
            call scifort_multivariate_normal_fit_f64(2_c_int,3_c_int,data,actual_mean,actual_cov, &
                status,use_mean,fixed_mean,use_cov,fixed_cov)
            if (status/=native_status .or. status/=0) stop 1
            if (any(actual_mean/=expected_mean) .or. any(actual_cov/=expected_cov)) stop 2
        end do
    end do
    fixed_mean=ieee_value(0.0_c_double,ieee_quiet_nan)
    fixed_cov=ieee_value(0.0_c_double,ieee_quiet_nan)
    call scifort_multivariate_normal_fit_f64(2_c_int,3_c_int,data,actual_mean,actual_cov, &
        status,0_c_int,fixed_mean,0_c_int,fixed_cov)
    if (status/=0 .or. any(ieee_is_nan(actual_mean)) .or. any(ieee_is_nan(actual_cov))) stop 3
    actual_mean=77.0_c_double
    actual_cov=88.0_c_double
    call scifort_multivariate_normal_fit_f64(2_c_int,3_c_int,data,actual_mean,actual_cov, &
        status,1_c_int,fixed_mean,2_c_int,fixed_cov)
    if (status/=1 .or. any(actual_mean/=77.0_c_double) .or. any(actual_cov/=88.0_c_double)) stop 4
    call scifort_multivariate_normal_fit_f64(2_c_int,3_c_int,data,actual_mean,actual_cov, &
        status,0_c_int,fixed_mean,1_c_int,fixed_cov)
    if (status/=2 .or. any(actual_mean/=[3.0_c_double,4.0_c_double])) stop 5
    if (.not. all(ieee_is_nan(actual_cov))) stop 6
    print '(a)', 'Multivariate fit native/C checks passed.'
end program test_multivariate_fit_c_api
