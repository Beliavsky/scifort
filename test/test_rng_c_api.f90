! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_rng_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_random, only : rng_state, rng_seed, rng_get_state, rng_fill_uniform
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64, scifort_rng_uniform_f64
    use scifort_dirichlet, only : dirichlet_rvs_array
    use scifort_c_api_generated_dirichlet, only : scifort_dirichlet_rvs_array_f64
    use scifort_multinomial, only : multinomial_rvs_array
    use scifort_multivariate_hypergeom, only : multivariate_hypergeom_rvs_array
    use scifort_c_api_generated_multinomial, only : scifort_multinomial_rvs_array_i32
    use scifort_c_api_generated_multivariate_hypergeom, only : &
        scifort_multivariate_hypergeom_rvs_array_i32
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6), saved(6)
    integer(int64) :: native_values(6), seeds(6)
    real(c_double) :: actual(12), expected(12), alpha(3), draws(3,4), reference(3,4)
    integer(c_int) :: status
    integer :: i, info
    integer(c_int) :: counts(3,4), population(3)
    integer :: native_counts(3,4)
    seeds=[0_int64,1_int64,-1_int64,2_int64**40,huge(0_int64),-huge(0_int64)-1_int64]
    do i=1,size(seeds)
        call rng_seed(native,seeds(i))
        call rng_get_state(native,native_values)
        call scifort_rng_seed_i64(int(seeds(i),c_int64_t),serialized,status)
        if (status/=0 .or. any(serialized/=native_values)) stop 1
        call rng_fill_uniform(native,expected)
        call scifort_rng_uniform_f64(12_c_int,serialized,actual,status)
        if (status/=0 .or. any(actual/=expected)) stop 2
        if (any(actual<=0.0_c_double) .or. any(actual>=1.0_c_double)) stop 3
        call rng_get_state(native,native_values)
        if (any(serialized/=native_values)) stop 4
        saved=serialized
        actual=77.0_c_double
        call scifort_rng_uniform_f64(0_c_int,serialized,actual,status)
        if (status/=0 .or. any(serialized/=saved) .or. any(actual/=77.0_c_double)) stop 5
    end do
    call rng_seed(native,777_int64)
    call scifort_rng_seed_i64(777_c_int64_t,serialized,status)
    alpha=[2.0_c_double,3.0_c_double,4.0_c_double]
    call dirichlet_rvs_array(native,reference,alpha,info)
    call scifort_dirichlet_rvs_array_f64(3_c_int,4_c_int,serialized,draws,alpha,status)
    if (status/=info .or. status/=0 .or. any(draws/=reference)) stop 6
    call rng_get_state(native,native_values)
    if (any(serialized/=native_values)) stop 7
    saved=serialized
    alpha=-1.0_c_double
    call scifort_dirichlet_rvs_array_f64(3_c_int,0_c_int,serialized,draws,alpha,status)
    if (status/=0 .or. any(serialized/=saved)) stop 8
    call scifort_dirichlet_rvs_array_f64(3_c_int,4_c_int,serialized,draws,alpha,status)
    if (status/=1 .or. any(serialized/=saved) .or. .not. all(ieee_is_nan(draws))) stop 9
    serialized=0_c_int64_t
    draws=88.0_c_double
    call scifort_dirichlet_rvs_array_f64(3_c_int,4_c_int,serialized,draws,alpha,status)
    if (status/=3 .or. any(serialized/=0_c_int64_t) .or. any(draws/=88.0_c_double)) stop 10
    actual=77.0_c_double
    call scifort_rng_uniform_f64(12_c_int,serialized,actual,status)
    if (status/=3 .or. any(actual/=77.0_c_double) .or. any(serialized/=0_c_int64_t)) stop 11
    do i=1,size(seeds)
        call rng_seed(native,seeds(i))
        call scifort_rng_seed_i64(int(seeds(i),c_int64_t),serialized,status)
        alpha=[0.2_c_double,0.3_c_double,0.5_c_double]
        call multinomial_rvs_array(native,native_counts,12,alpha,info)
        call scifort_multinomial_rvs_array_i32(3_c_int,4_c_int,serialized,counts, &
            12_c_int,alpha,status)
        if (status/=info .or. status/=0 .or. any(counts/=native_counts)) stop 12
        call rng_get_state(native,native_values)
        if (any(serialized/=native_values)) stop 13
        population=[10_c_int,20_c_int,30_c_int]
        call multivariate_hypergeom_rvs_array(native,int(population),12,native_counts,info)
        call scifort_multivariate_hypergeom_rvs_array_i32(3_c_int,4_c_int,serialized, &
            population,12_c_int,counts,status)
        if (status/=info .or. status/=0 .or. any(counts/=native_counts)) stop 14
        call rng_get_state(native,native_values)
        if (any(serialized/=native_values)) stop 15
    end do
    print '(a)', 'RNG and continuous/discrete sampling native/C stream checks passed.'
end program test_rng_c_api
