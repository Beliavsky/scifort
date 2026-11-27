! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_random_matrix_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_ortho_group, only : ortho_group_rvs, ortho_group_rvs_array
    use scifort_special_ortho_group, only : special_ortho_group_rvs, special_ortho_group_rvs_array
    use scifort_random_correlation, only : random_correlation_rvs
    use scifort_c_api_generated_ortho_group, only : scifort_ortho_group_rvs_f64, &
        scifort_ortho_group_rvs_array_f64
    use scifort_c_api_generated_special_ortho_group, only : scifort_special_ortho_group_rvs_f64, &
        scifort_special_ortho_group_rvs_array_f64
    use scifort_c_api_generated_random_correlation, only : scifort_random_correlation_rvs_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: actual(3,3), expected(3,3), batch(3,3,4), reference(3,3,4)
    real(c_double) :: eigs(3)
    integer(c_int) :: status
    integer :: method, info

    do method=1,3
        call rng_seed(native,2_int64**40+123_int64)
        call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
        if (status/=0) stop 1
        select case(method)
        case(1)
            call ortho_group_rvs(native,expected,info)
            call scifort_ortho_group_rvs_f64(3_c_int,serialized,actual,status)
            if (status/=info .or. status/=0 .or. any(actual/=expected)) stop 2
            call ortho_group_rvs_array(native,reference,info)
            call scifort_ortho_group_rvs_array_f64(3_c_int,4_c_int,serialized,batch,status)
            if (status/=info .or. status/=0 .or. any(batch/=reference)) stop 3
        case(2)
            call special_ortho_group_rvs(native,expected,info)
            call scifort_special_ortho_group_rvs_f64(3_c_int,serialized,actual,status)
            if (status/=info .or. status/=0 .or. any(actual/=expected)) stop 4
            call special_ortho_group_rvs_array(native,reference,info)
            call scifort_special_ortho_group_rvs_array_f64(3_c_int,4_c_int,serialized,batch,status)
            if (status/=info .or. status/=0 .or. any(batch/=reference)) stop 5
        case(3)
            eigs=[0.2_c_double,0.8_c_double,2.0_c_double]
            call random_correlation_rvs(native,eigs,expected,info)
            call scifort_random_correlation_rvs_f64(3_c_int,serialized,eigs,actual,status, &
                1.0e-13_c_double,1.0e-7_c_double)
            if (status/=info .or. status/=0 .or. any(actual/=expected)) stop 6
            call rng_get_state(native,expected_state)
            if (any(serialized/=expected_state)) stop 7
            ! Match nondefault tolerances, including any partially transformed failure output.
            call random_correlation_rvs(native,eigs,expected,info,tol=1.0e-12_c_double, &
                diag_tol=0.0_c_double)
            call scifort_random_correlation_rvs_f64(3_c_int,serialized,eigs,actual,status, &
                1.0e-12_c_double,0.0_c_double)
            if (status/=info .or. any(actual/=expected)) stop 10
            call rng_get_state(native,expected_state)
            if (any(serialized/=expected_state)) stop 11
            eigs(1)=-1.0_c_double
            call random_correlation_rvs(native,eigs,expected,info)
            call scifort_random_correlation_rvs_f64(3_c_int,serialized,eigs,actual,status, &
                1.0e-13_c_double,1.0e-7_c_double)
            if (status/=2 .or. status/=info .or. any(actual/=expected)) stop 8
        end select
        call rng_get_state(native,expected_state)
        if (any(serialized/=expected_state)) stop 9
    end do
    print '(a)', 'Random matrix native/C sampling checks passed.'
end program test_random_matrix_c_api
