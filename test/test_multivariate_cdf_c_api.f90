! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_multivariate_cdf_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_multivariate_normal, only : multivariate_normal_cdf, multivariate_normal_logcdf
    use scifort_multivariate_t, only : multivariate_t_cdf
    use scifort_c_api_generated_multivariate_normal, only : &
        scifort_multivariate_normal_cdf_f64, scifort_multivariate_normal_logcdf_f64
    use scifort_c_api_generated_multivariate_t, only : scifort_multivariate_t_cdf_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: x(2), mean(2), matrix(2,2), low(2), actual, expected
    real(c_double), allocatable :: optional_mean(:), optional_matrix(:,:), optional_low(:)
    integer(c_int) :: status, has_mean, has_matrix, has_low
    integer :: mask, method
    x=[0.5_c_double,0.3_c_double]
    mean=[0.1_c_double,-0.2_c_double]
    low=[-1.0_c_double,-0.5_c_double]
    matrix=reshape([1.0_c_double,0.4_c_double,77.0_c_double,2.0_c_double],[2,2])
    do mask=0,7
        has_mean=int(ibits(mask,0,1),c_int)
        has_matrix=int(ibits(mask,1,1),c_int)
        has_low=int(ibits(mask,2,1),c_int)
        if (has_mean==1) allocate(optional_mean,source=mean)
        if (has_matrix==1) allocate(optional_matrix,source=matrix)
        if (has_low==1) allocate(optional_low,source=low)
        do method=1,3
            call rng_seed(native,777_int64)
            call scifort_rng_seed_i64(777_c_int64_t,serialized,status)
            select case(method)
            case(1)
                expected=multivariate_normal_cdf(x,optional_mean,optional_matrix,.false., &
                    4000,1.0e-5_c_double,1.0e-5_c_double,optional_low,native)
                call scifort_multivariate_normal_cdf_f64(2_c_int,x,has_mean,mean,has_matrix,matrix, &
                    0_c_int,4000_c_int,1.0e-5_c_double,1.0e-5_c_double,has_low,low,serialized,actual,status)
            case(2)
                expected=multivariate_normal_logcdf(x,optional_mean,optional_matrix,.false., &
                    4000,1.0e-5_c_double,1.0e-5_c_double,optional_low,native)
                call scifort_multivariate_normal_logcdf_f64(2_c_int,x,has_mean,mean,has_matrix,matrix, &
                    0_c_int,4000_c_int,1.0e-5_c_double,1.0e-5_c_double,has_low,low,serialized,actual,status)
            case(3)
                expected=multivariate_t_cdf(x,optional_mean,optional_matrix,5.0_c_double,.false., &
                    4000,optional_low,native)
                call scifort_multivariate_t_cdf_f64(2_c_int,x,has_mean,mean,has_matrix,matrix, &
                    5.0_c_double,0_c_int,4000_c_int,has_low,low,serialized,actual,status)
            end select
            if (status/=0 .or. actual/=expected) stop 1
            call rng_get_state(native,expected_state)
            if (any(serialized/=expected_state)) stop 2
        end do
        if (allocated(optional_mean)) deallocate(optional_mean)
        if (allocated(optional_matrix)) deallocate(optional_matrix)
        if (allocated(optional_low)) deallocate(optional_low)
    end do
    print '(a)', 'Multivariate CDF/log-CDF native/C checks passed.'
end program test_multivariate_cdf_c_api
