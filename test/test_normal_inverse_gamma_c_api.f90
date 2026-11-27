! SPDX-License-Identifier: MIT
program test_normal_inverse_gamma_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_normal_inverse_gamma
    use scifort_c_api_generated_normal_inverse_gamma
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: x, s2, ex, es, value, expected, xs(4), ss(4), exs(4), ess(4)
    integer(c_int) :: status

    expected=normal_inverse_gamma_pdf(2.0_c_double,3.0_c_double)
    call scifort_normal_inverse_gamma_pdf_f64(2.0_c_double,3.0_c_double, &
        0.0_c_double,1.0_c_double,1.0_c_double,1.0_c_double,value,status)
    if (status/=0 .or. value/=expected) stop 1
    expected=normal_inverse_gamma_logpdf(2.0_c_double,3.0_c_double)
    call scifort_normal_inverse_gamma_logpdf_f64(2.0_c_double,3.0_c_double, &
        0.0_c_double,1.0_c_double,1.0_c_double,1.0_c_double,value,status)
    if (status/=0 .or. value/=expected) stop 2
    call normal_inverse_gamma_mean(a=3.0_c_double,mean_x=ex,mean_s2=es)
    call scifort_normal_inverse_gamma_mean_f64(0.0_c_double,1.0_c_double, &
        3.0_c_double,1.0_c_double,x,s2,status)
    if (status/=0 .or. x/=ex .or. s2/=es) stop 3
    call normal_inverse_gamma_var(a=3.0_c_double,var_x=ex,var_s2=es)
    call scifort_normal_inverse_gamma_var_f64(0.0_c_double,1.0_c_double, &
        3.0_c_double,1.0_c_double,x,s2,status)
    if (status/=0 .or. x/=ex .or. s2/=es) stop 4
    call rng_seed(native,2_int64**40+123_int64)
    call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
    if (status/=0) stop 5
    call normal_inverse_gamma_rvs(native,ex,es)
    call scifort_normal_inverse_gamma_rvs_f64(serialized,x,s2, &
        0.0_c_double,1.0_c_double,1.0_c_double,1.0_c_double,status)
    if (status/=0 .or. x/=ex .or. s2/=es) stop 6
    call normal_inverse_gamma_rvs_array(native,exs,ess)
    call scifort_normal_inverse_gamma_rvs_array_f64(4_c_int,serialized,xs,ss, &
        0.0_c_double,1.0_c_double,1.0_c_double,1.0_c_double,status)
    if (status/=0 .or. any(xs/=exs) .or. any(ss/=ess)) stop 7
    call rng_get_state(native,expected_state)
    if (any(serialized/=expected_state)) stop 8
    print '(a)', 'Normal-inverse-gamma native/C checks passed.'
end program test_normal_inverse_gamma_c_api
