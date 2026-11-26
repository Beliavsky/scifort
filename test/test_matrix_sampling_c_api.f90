! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_matrix_sampling_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use scifort_random, only : rng_state, rng_seed, rng_get_state, rng_fill_uniform
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64, scifort_rng_uniform_f64
    use scifort_matrix_normal, only : matrix_normal_rvs_array
    use scifort_matrix_t, only : matrix_t_rvs_array
    use scifort_c_api_generated_matrix_normal, only : scifort_matrix_normal_rvs_array_f64
    use scifort_c_api_generated_matrix_t, only : scifort_matrix_t_rvs_array_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: mean(2,3), row(2,2), col(3,3), actual(2,3,4), expected(2,3,4)
    real(c_double), allocatable :: optional_mean(:,:), optional_row(:,:), optional_col(:,:)
    real(c_double) :: df, discarded(10), small(1,1,3), small_expected(1,1,3)
    real(c_double) :: dummy_mean(1,1), dummy_row(1,1), dummy_col(1,1)
    integer(c_int) :: status, has_mean, has_row, has_col
    integer :: mask, nu, info, i
    mean=reshape([(real(i,c_double),i=1,6)],[2,3])
    row=reshape([2.0_c_double,0.3_c_double,99.0_c_double,1.0_c_double],[2,2])
    col=0.1_c_double
    do i=1,3
        col(i,i)=real(i,c_double)
    end do
    do mask=0,7
        has_mean=int(ibits(mask,0,1),c_int)
        has_row=int(ibits(mask,1,1),c_int)
        has_col=int(ibits(mask,2,1),c_int)
        if (has_mean==1) allocate(optional_mean,source=mean)
        if (has_row==1) allocate(optional_row,source=row)
        if (has_col==1) allocate(optional_col,source=col)
        call rng_seed(native,2_int64**40+123_int64)
        call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
        call matrix_normal_rvs_array(native,expected,optional_mean,optional_row,optional_col,info)
        call scifort_matrix_normal_rvs_array_f64(2_c_int,3_c_int,4_c_int,serialized,actual, &
            has_mean,mean,has_row,row,has_col,col,status)
        if (status/=0 .or. info/=0 .or. any(actual/=expected)) stop 1
        call rng_get_state(native,expected_state)
        if (any(serialized/=expected_state)) stop 2
        do nu=1,2
            df=real(4*nu-3,c_double)
            call rng_seed(native,2_int64**40+123_int64)
            call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
            if (nu==1) then
                call matrix_t_rvs_array(native,expected,optional_mean,optional_row,optional_col,status=info)
            else
                call matrix_t_rvs_array(native,expected,optional_mean,optional_row,optional_col,df,info)
            end if
            call scifort_matrix_t_rvs_array_f64(2_c_int,3_c_int,4_c_int,serialized,actual, &
                has_mean,mean,has_row,row,has_col,col,df,status)
            if (status/=0 .or. info/=0 .or. any(actual/=expected)) stop 3
            call rng_get_state(native,expected_state)
            if (any(serialized/=expected_state)) stop 4
        end do
        if (allocated(optional_mean)) deallocate(optional_mean)
        if (allocated(optional_row)) deallocate(optional_row)
        if (allocated(optional_col)) deallocate(optional_col)
    end do
    row(1,1)=-1.0_c_double
    actual=77.0_c_double
    expected=88.0_c_double
    call matrix_normal_rvs_array(native,expected,rowcov=row,status=info)
    call scifort_matrix_normal_rvs_array_f64(2_c_int,3_c_int,4_c_int,serialized,actual, &
        0_c_int,mean,1_c_int,row,0_c_int,col,status)
    if (status/=3 .or. info/=3 .or. .not. all(ieee_is_nan(actual))) stop 5
    if (.not. all(ieee_is_nan(expected))) stop 6
    call matrix_t_rvs_array(native,expected,row_spread=row,df=5.0_c_double,status=info)
    call scifort_matrix_t_rvs_array_f64(2_c_int,3_c_int,4_c_int,serialized,actual, &
        0_c_int,mean,1_c_int,row,0_c_int,col,5.0_c_double,status)
    if (status/=3 .or. info/=3 .or. .not. all(ieee_is_nan(actual))) stop 7
    if (.not. all(ieee_is_nan(expected))) stop 8
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 9
    call rng_seed(native,2_int64**40+123_int64)
    call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
    call rng_fill_uniform(native,discarded)
    call scifort_rng_uniform_f64(10_c_int,serialized,discarded,status)
    dummy_mean=0.0_c_double
    dummy_row=0.0_c_double
    dummy_col=0.0_c_double
    small=77.0_c_double
    small_expected=88.0_c_double
    call matrix_t_rvs_array(native,small_expected,df=0.0001_c_double,status=info)
    call scifort_matrix_t_rvs_array_f64(1_c_int,1_c_int,3_c_int,serialized,small, &
        0_c_int,dummy_mean,0_c_int,dummy_row,0_c_int,dummy_col,0.0001_c_double,status)
    if (status/=3 .or. info/=3 .or. .not. ieee_is_finite(small(1,1,1))) stop 10
    if (small(1,1,1)/=small_expected(1,1,1)) stop 11
    if (.not. all(ieee_is_nan(small(:,:,2:)))) stop 12
    if (.not. all(ieee_is_nan(small_expected(:,:,2:)))) stop 13
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 14
    print '(a)', 'Matrix normal/t native/C sampling checks passed.'
end program test_matrix_sampling_c_api
