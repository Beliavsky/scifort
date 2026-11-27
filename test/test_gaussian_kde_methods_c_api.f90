! SPDX-License-Identifier: MIT
program test_gaussian_kde_methods_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr, c_associated, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_gaussian_kde, only : gaussian_kde, gaussian_kde_init
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_c_api_generated_gaussian_kde
    implicit none
    type(gaussian_kde) :: native, child_native
    type(rng_state) :: state
    type(c_ptr) :: handle, child
    real(c_double) :: data(2,6), weights(6), low(2), high(2), mean(2), cov(2,2)
    real(c_double) :: value, expected, draws(2,7), reference(2,7), empty(2,0)
    real(c_double) :: point(1,1), evaluated(1)
    real(c_double), allocatable :: native_pdf(:)
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: snapshot(6)
    integer(c_int) :: status
    integer :: info

    data(1,:) = [-2.0_c_double,-0.5_c_double,0.25_c_double,1.5_c_double,3.0_c_double,0.8_c_double]
    data(2,:) = [0.5_c_double,2.0_c_double,-1.0_c_double,0.2_c_double,2.8_c_double,-0.4_c_double]
    weights = [1.0_c_double,2.0_c_double,0.0_c_double,3.0_c_double,1.0_c_double,2.0_c_double]
    call gaussian_kde_init(native,data,info,weights=weights)
    if (info /= 0) stop 1
    handle = c_null_ptr
    call scifort_gaussian_kde_create_f64(2_c_int,6_c_int,data,weights,0_c_int,1.0_c_double,handle,status)
    if (status /= 0 .or. .not. c_associated(handle)) stop 2
    mean = 0.0_c_double
    cov = 0.0_c_double
    cov(1,1) = 1.0_c_double
    cov(2,2) = 1.0_c_double
    expected = native%integrate_gaussian(mean,cov,info)
    call scifort_gaussian_kde_integrate_gaussian_f64(2_c_int,handle,mean,cov,value,status)
    if (status /= 0 .or. value /= expected) stop 3
    expected = native%integrate_kde(native,info)
    call scifort_gaussian_kde_integrate_kde_f64(handle,handle,value,status)
    if (status /= 0 .or. value /= expected) stop 4

    call rng_seed(state,123_int64)
    call scifort_rng_seed_i64(123_c_int64_t,serialized,status)
    if (status /= 0) stop 5
    low = -1.0_c_double
    high = 2.0_c_double
    expected = native%integrate_box(low,high,maxpts=128,state=state,status=info)
    call scifort_gaussian_kde_integrate_box_f64(2_c_int,handle,low,high,128_c_int,serialized,value,status)
    if (status /= 0 .or. info /= 0 .or. value /= expected) stop 6
    call rng_get_state(state,snapshot)
    if (any(serialized /= snapshot)) stop 7
    call native%resample(state,reference,info)
    call scifort_gaussian_kde_resample_f64(2_c_int,7_c_int,handle,serialized,draws,status)
    if (status /= 0 .or. info /= 0 .or. any(draws /= reference)) stop 8
    call rng_get_state(state,snapshot)
    if (any(serialized /= snapshot)) stop 9
    call scifort_gaussian_kde_resample_f64(2_c_int,0_c_int,handle,serialized,empty,status)
    if (status /= 0 .or. any(serialized /= snapshot)) stop 10

    child_native = native%marginal([2],info)
    if (info /= 0) stop 11
    child = c_null_ptr
    call scifort_gaussian_kde_marginal_f64(1_c_int,handle,[2_c_int],child,status)
    if (status /= 0 .or. .not. c_associated(child)) stop 12
    call scifort_gaussian_kde_destroy_f64(handle,status)
    if (status /= 0 .or. c_associated(handle)) stop 13
    point = 0.4_c_double
    native_pdf = child_native%pdf(point,info)
    call scifort_gaussian_kde_pdf_f64(1_c_int,1_c_int,child,point,evaluated,status)
    if (status /= 0 .or. any(evaluated /= native_pdf)) stop 14
    expected = child_native%integrate_box_1d(-1.0_c_double,2.0_c_double,info)
    call scifort_gaussian_kde_integrate_box_1d_f64(child,-1.0_c_double,2.0_c_double,value,status)
    if (status /= 0 .or. value /= expected) stop 15
    call scifort_gaussian_kde_destroy_f64(child,status)
    if (status /= 0 .or. c_associated(child)) stop 16
    print '(a)', 'Gaussian KDE native/C integration, resampling and marginal checks passed.'
end program test_gaussian_kde_methods_c_api
