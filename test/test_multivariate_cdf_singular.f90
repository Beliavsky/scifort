! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_multivariate_cdf_singular
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_multivariate_t, only : multivariate_t_cdf
    use scifort_random, only : rng_state, rng_seed
    use scifort_student_t, only : t_cdf
    use test_checks, only : check_close
    implicit none
    type(rng_state) :: state
    real(dp) :: matrix(2,2), actual, expected
    integer :: failures
    failures = 0
    call rng_seed(state, 123_int64)
    matrix = 1.0_dp
    actual = multivariate_t_cdf([0.0_dp, 0.0_dp], shape=matrix, df=5.0_dp, &
        allow_singular=.true., maxpts=16000, state=state)
    call check_close('perfect correlation orthant', actual, 0.5_dp, 1.0e-3_dp, 0.0_dp, failures)
    actual = multivariate_t_cdf([0.8_dp, 0.6_dp], shape=matrix, df=5.0_dp, &
        lower_limit=[-0.4_dp, -0.2_dp], allow_singular=.true., maxpts=16000, state=state)
    expected = t_cdf(0.6_dp, 5.0_dp) - t_cdf(-0.2_dp, 5.0_dp)
    call check_close('rank one interval intersection', actual, expected, 1.0e-3_dp, 0.0_dp, failures)
    matrix = reshape([1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp], [2,2])
    actual = multivariate_t_cdf([0.5_dp, 0.8_dp], shape=matrix, df=5.0_dp, &
        allow_singular=.true., maxpts=16000, state=state)
    expected = t_cdf(0.5_dp, 5.0_dp) - t_cdf(-0.8_dp, 5.0_dp)
    call check_close('negative correlation interval', actual, expected, 1.0e-3_dp, 0.0_dp, failures)
    matrix = 0.0_dp
    actual = multivariate_t_cdf([0.0_dp, 0.0_dp], shape=matrix, df=5.0_dp, &
        allow_singular=.true., maxpts=100, state=state)
    call check_close('zero rank inclusive atom', actual, 1.0_dp, 0.0_dp, 0.0_dp, failures)
    actual = multivariate_t_cdf([-0.1_dp, 0.0_dp], shape=matrix, df=5.0_dp, &
        allow_singular=.true., maxpts=100, state=state)
    call check_close('zero rank outside atom', actual, 0.0_dp, 0.0_dp, 0.0_dp, failures)
    matrix(2,2) = 1.0_dp
    actual = multivariate_t_cdf([0.0_dp, 0.5_dp], shape=matrix, df=5.0_dp, &
        allow_singular=.true., maxpts=16000, state=state)
    call check_close('mixed deterministic coordinate', actual, t_cdf(0.5_dp, 5.0_dp), &
        1.0e-3_dp, 0.0_dp, failures)
    if (failures /= 0) error stop 1
    print '(a)', 'test_multivariate_cdf_singular: PASS'
end program test_multivariate_cdf_singular
