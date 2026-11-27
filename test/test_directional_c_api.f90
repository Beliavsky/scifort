! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_directional_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_uniform_direction, only : uniform_direction_rvs, uniform_direction_rvs_array
    use scifort_vonmises_fisher, only : vonmises_fisher_pdf, vonmises_fisher_logpdf, &
        vonmises_fisher_entropy, vonmises_fisher_rvs, vonmises_fisher_rvs_array, vonmises_fisher_fit
    use scifort_c_api_generated_uniform_direction, only : scifort_uniform_direction_rvs_f64, &
        scifort_uniform_direction_rvs_array_f64
    use scifort_c_api_generated_vonmises_fisher, only : scifort_vonmises_fisher_pdf_f64, &
        scifort_vonmises_fisher_logpdf_f64, scifort_vonmises_fisher_entropy_f64, &
        scifort_vonmises_fisher_rvs_f64, scifort_vonmises_fisher_rvs_array_f64, &
        scifort_vonmises_fisher_fit_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double), allocatable :: mu(:), actual(:), expected(:), samples(:,:), batch(:,:), fit_mu(:), ref_mu(:)
    real(c_double) :: value, kappa, reference_kappa
    integer(c_int) :: status, n
    integer :: i, info
    integer, parameter :: dimensions(3)=[2,3,5]

    do i=1,3
        n=dimensions(i)
        allocate(mu(n),actual(n),expected(n),samples(n,12),batch(n,12),fit_mu(n),ref_mu(n))
        mu=0.0_c_double
        mu(n)=1.0_c_double
        actual=mu
        call scifort_vonmises_fisher_pdf_f64(n,actual,mu,3.0_c_double,value,status)
        if (status/=0 .or. value/=vonmises_fisher_pdf(mu,mu,3.0_c_double)) stop 1
        call scifort_vonmises_fisher_logpdf_f64(n,actual,mu,3.0_c_double,value,status)
        if (status/=0 .or. value/=vonmises_fisher_logpdf(mu,mu,3.0_c_double)) stop 2
        call scifort_vonmises_fisher_entropy_f64(n,mu,3.0_c_double,value,status)
        if (status/=0 .or. value/=vonmises_fisher_entropy(mu,3.0_c_double)) stop 3
        call rng_seed(native,2_int64**40+123_int64)
        call scifort_rng_seed_i64(2_c_int64_t**40+123_c_int64_t,serialized,status)
        call uniform_direction_rvs(native,expected,info)
        call scifort_uniform_direction_rvs_f64(n,serialized,actual,status)
        if (status/=0 .or. status/=info .or. any(actual/=expected)) stop 4
        call uniform_direction_rvs_array(native,batch,info)
        call scifort_uniform_direction_rvs_array_f64(n,12_c_int,serialized,samples,status)
        if (status/=0 .or. status/=info .or. any(samples/=batch)) stop 5
        call vonmises_fisher_rvs(native,mu,3.0_c_double,expected,info)
        call scifort_vonmises_fisher_rvs_f64(n,serialized,mu,3.0_c_double,actual,status)
        if (status/=0 .or. status/=info .or. any(actual/=expected)) stop 6
        call vonmises_fisher_rvs_array(native,mu,3.0_c_double,batch,info)
        call scifort_vonmises_fisher_rvs_array_f64(n,12_c_int,serialized,mu,3.0_c_double,samples,status)
        if (status/=0 .or. status/=info .or. any(samples/=batch)) stop 7
        call rng_get_state(native,expected_state)
        if (any(serialized/=expected_state)) stop 8
        call vonmises_fisher_fit(samples,ref_mu,reference_kappa,info)
        call scifort_vonmises_fisher_fit_f64(n,12_c_int,samples,fit_mu,kappa,status)
        if (status/=0 .or. status/=info .or. any(fit_mu/=ref_mu) .or. kappa/=reference_kappa) stop 9
        deallocate(mu,actual,expected,samples,batch,fit_mu,ref_mu)
    end do
    print '(a)', 'Directional native/C checks passed.'
end program test_directional_c_api
