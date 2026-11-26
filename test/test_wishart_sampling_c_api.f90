! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_wishart_sampling_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use scifort_random, only : rng_state, rng_seed, rng_get_state, rng_fill_uniform
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64, scifort_rng_uniform_f64
    use scifort_wishart, only : wishart_rvs_array
    use scifort_invwishart, only : invwishart_rvs_array
    use scifort_c_api_generated_wishart, only : scifort_wishart_rvs_array_f64
    use scifort_c_api_generated_invwishart, only : scifort_invwishart_rvs_array_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: actual(2,2,4), expected(2,2,4), scale(2,2)
    real(c_double) :: tiny(1,1,3), tiny_expected(1,1,3), unit_scale(1,1)
    real(c_double) :: discarded(10)
    integer(c_int) :: status
    integer :: method, info
    scale=reshape([2.0_c_double,0.4_c_double,77.0_c_double,1.0_c_double],[2,2])
    do method=1,2
        call rng_seed(native,2_int64**40)
        call scifort_rng_seed_i64(2_c_int64_t**40,serialized,status)
        if (method==1) then
            call wishart_rvs_array(native,7.0_c_double,scale,expected,info)
            call scifort_wishart_rvs_array_f64(2_c_int,4_c_int,serialized, &
                7.0_c_double,scale,actual,status)
        else
            call invwishart_rvs_array(native,7.0_c_double,scale,expected,info)
            call scifort_invwishart_rvs_array_f64(2_c_int,4_c_int,serialized, &
                7.0_c_double,scale,actual,status)
        end if
        if (status/=info .or. status/=0 .or. any(actual/=expected)) stop 1
        call rng_get_state(native,expected_state)
        if (any(serialized/=expected_state)) stop 2
        actual=77.0_c_double
        expected=88.0_c_double
        if (method==1) then
            call wishart_rvs_array(native,0.0_c_double,scale,expected,info)
            call scifort_wishart_rvs_array_f64(2_c_int,4_c_int,serialized, &
                0.0_c_double,scale,actual,status)
        else
            call invwishart_rvs_array(native,0.0_c_double,scale,expected,info)
            call scifort_invwishart_rvs_array_f64(2_c_int,4_c_int,serialized, &
                0.0_c_double,scale,actual,status)
        end if
        if (status/=2 .or. info/=2 .or. .not. all(ieee_is_nan(actual))) stop 3
        if (.not. all(ieee_is_nan(expected))) stop 4
        call rng_get_state(native,expected_state)
        if (any(serialized/=expected_state)) stop 5
    end do
    call rng_seed(native,11_int64)
    call scifort_rng_seed_i64(11_c_int64_t,serialized,status)
    unit_scale=1.0_c_double
    tiny=77.0_c_double
    tiny_expected=88.0_c_double
    call invwishart_rvs_array(native,1.0e-300_c_double,unit_scale,tiny_expected,info)
    call scifort_invwishart_rvs_array_f64(1_c_int,3_c_int,serialized, &
        1.0e-300_c_double,unit_scale,tiny,status)
    if (status/=3 .or. info/=3 .or. .not. all(ieee_is_nan(tiny))) stop 6
    if (.not. all(ieee_is_nan(tiny_expected))) stop 7
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 8
    call rng_seed(native,11_int64)
    call scifort_rng_seed_i64(11_c_int64_t,serialized,status)
    call wishart_rvs_array(native,1.0e-300_c_double,unit_scale,tiny_expected,info)
    call scifort_wishart_rvs_array_f64(1_c_int,3_c_int,serialized, &
        1.0e-300_c_double,unit_scale,tiny,status)
    if (status/=0 .or. info/=0 .or. .not. all(ieee_is_nan(tiny))) stop 9
    if (.not. all(ieee_is_nan(tiny_expected))) stop 10
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 11
    call rng_seed(native,2_int64**40+123_int64)
    call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
    call rng_fill_uniform(native,discarded)
    call scifort_rng_uniform_f64(10_c_int,serialized,discarded,status)
    call invwishart_rvs_array(native,0.0001_c_double,unit_scale,tiny_expected,info)
    call scifort_invwishart_rvs_array_f64(1_c_int,3_c_int,serialized, &
        0.0001_c_double,unit_scale,tiny,status)
    if (status/=3 .or. info/=3) stop 12
    if (.not. ieee_is_finite(tiny(1,1,1)) .or. tiny(1,1,1)/=tiny_expected(1,1,1)) stop 13
    if (.not. all(ieee_is_nan(tiny(:,:,2:)))) stop 14
    if (.not. all(ieee_is_nan(tiny_expected(:,:,2:)))) stop 15
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 16
    print '(a)', 'Wishart/inverse-Wishart native/C sampling checks passed.'
end program test_wishart_sampling_c_api
