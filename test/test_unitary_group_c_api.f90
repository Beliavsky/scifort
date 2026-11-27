! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_unitary_group_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_unitary_group, only : unitary_group_rvs, unitary_group_rvs_array
    use scifort_c_api_generated_unitary_group, only : scifort_unitary_group_rvs_f64, &
        scifort_unitary_group_rvs_array_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: actual(2,3,3), batch(2,3,3,4)
    complex(c_double) :: expected(3,3), reference(3,3,4)
    integer(c_int) :: status
    integer :: info

    call rng_seed(native,2_int64**40+123_int64)
    call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
    if (status/=0) stop 1
    call unitary_group_rvs(native,expected,info)
    call scifort_unitary_group_rvs_f64(3_c_int,serialized,actual,status)
    if (status/=0 .or. status/=info) stop 2
    if (any(cmplx(actual(1,:,:),actual(2,:,:),kind=c_double)/=expected)) stop 3
    call unitary_group_rvs_array(native,reference,info)
    call scifort_unitary_group_rvs_array_f64(3_c_int,4_c_int,serialized,batch,status)
    if (status/=0 .or. status/=info) stop 4
    if (any(cmplx(batch(1,:,:,:),batch(2,:,:,:),kind=c_double)/=reference)) stop 5
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 6
    print '(a)', 'Unitary native/C interleaved sampling checks passed.'
end program test_unitary_group_c_api
