! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_random_table_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_random_table, only : random_table_pmf, random_table_logpmf, random_table_mean, &
        random_table_rvs, random_table_rvs_array
    use scifort_c_api_generated_random_table, only : scifort_random_table_pmf_f64, &
        scifort_random_table_logpmf_f64, scifort_random_table_mean_f64, &
        scifort_random_table_rvs_f64, scifort_random_table_rvs_array_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    integer(c_int) :: row(2), col(3), x(2,3), expected(2,3), samples(2,3,7), batch(2,3,7), status
    integer(c_int) :: large_row(2), large_col(1), large_x(2,1)
    real(c_double) :: value, mean(2,3), large_mean(2,1)
    integer :: info

    row=[3,5]
    col=[2,4,2]
    x=reshape([1,1,1,3,1,1],[2,3])
    call scifort_random_table_pmf_f64(2_c_int,3_c_int,x,row,col,value,status)
    if (status/=0 .or. value/=random_table_pmf(x,row,col)) stop 1
    call scifort_random_table_logpmf_f64(2_c_int,3_c_int,x,row,col,value,status)
    if (status/=0 .or. value/=random_table_logpmf(x,row,col)) stop 2
    call scifort_random_table_mean_f64(2_c_int,3_c_int,row,col,mean,status)
    if (status/=0 .or. any(mean/=random_table_mean(row,col))) stop 3
    call rng_seed(native,2_int64**40+123_int64)
    call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
    call random_table_rvs(native,row,col,expected,info)
    call scifort_random_table_rvs_f64(2_c_int,3_c_int,serialized,row,col,x,status)
    if (status/=0 .or. status/=info .or. any(x/=expected)) stop 4
    call random_table_rvs_array(native,row,col,batch,info)
    call scifort_random_table_rvs_array_f64(2_c_int,3_c_int,7_c_int,serialized,row,col,samples,status)
    if (status/=0 .or. status/=info .or. any(samples/=batch)) stop 5
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 6

    ! The native module is also compiled with -ftrapv in the checked client build.
    large_row=[huge(0_c_int),1_c_int]
    large_col=[0_c_int]
    large_x=0
    if (.not. ieee_is_nan(random_table_pmf(large_x,large_row,large_col))) stop 7
    large_mean=random_table_mean(large_row,large_col)
    if (.not. all(ieee_is_nan(large_mean))) stop 8
    call random_table_rvs(native,large_row,large_col,large_x,info)
    if (info/=2 .or. any(large_x/=0)) stop 9
    large_row=[huge(0_c_int)-1_c_int,0_c_int]
    large_col=[huge(0_c_int)-1_c_int]
    large_x(:,1)=large_row
    if (random_table_pmf(large_x,large_row,large_col)/=1.0_c_double) stop 10
    large_row=[1_c_int,0_c_int]
    large_col=[1_c_int]
    large_x(:,1)=[huge(0_c_int),huge(0_c_int)]
    if (random_table_pmf(large_x,large_row,large_col)/=0.0_c_double) stop 11
    print '(a)', 'Random table native/C and overflow checks passed.'
end program test_random_table_c_api
