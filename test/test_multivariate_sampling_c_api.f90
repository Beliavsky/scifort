! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_multivariate_sampling_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_rng_c_bridge, only : scifort_rng_seed_i64
    use scifort_multivariate_normal, only : multivariate_normal_rvs_array
    use scifort_multivariate_t, only : multivariate_t_rvs_array
    use scifort_c_api_generated_multivariate_normal, only : &
        scifort_multivariate_normal_rvs_array_f64
    use scifort_c_api_generated_multivariate_t, only : scifort_multivariate_t_rvs_array_f64
    implicit none
    type(rng_state) :: native
    integer(c_int64_t) :: serialized(6)
    integer(int64) :: expected_state(6)
    real(c_double) :: actual(2,5), expected(2,5), mu(2), matrix(2,2), df(3)
    integer(c_int) :: status, has_mu, has_matrix
    integer :: mode, nu, info
    mu=[2.0_c_double,-1.0_c_double]
    matrix=reshape([4.0_c_double,1.0_c_double,77.0_c_double,3.0_c_double],[2,2])
    df=[1.0_c_double,5.0_c_double,ieee_value(0.0_c_double,ieee_positive_inf)]
    do mode=0,3
        has_mu=int(mod(mode,2),c_int)
        has_matrix=int(mode/2,c_int)
        call rng_seed(native,777_int64)
        call scifort_rng_seed_i64(777_c_int64_t,serialized,status)
        select case(mode)
        case(0)
            call multivariate_normal_rvs_array(native,expected,status=info)
        case(1)
            call multivariate_normal_rvs_array(native,expected,mean=mu,status=info)
        case(2)
            call multivariate_normal_rvs_array(native,expected,cov=matrix,status=info)
        case(3)
            call multivariate_normal_rvs_array(native,expected,mu,matrix,info)
        end select
        call scifort_multivariate_normal_rvs_array_f64(2_c_int,5_c_int,serialized,actual, &
            has_mu,mu,has_matrix,matrix,status)
        if (status/=info .or. status/=0 .or. any(actual/=expected)) stop 1
        call rng_get_state(native,expected_state)
        if (any(serialized/=expected_state)) stop 2
        do nu=1,size(df)
            call rng_seed(native,777_int64)
            call scifort_rng_seed_i64(777_c_int64_t,serialized,status)
            select case(mode)
            case(0)
                if (nu==1) then
                    call multivariate_t_rvs_array(native,expected,status=info)
                else
                    call multivariate_t_rvs_array(native,expected,df=df(nu),status=info)
                end if
            case(1)
                call multivariate_t_rvs_array(native,expected,loc=mu,df=df(nu),status=info)
            case(2)
                call multivariate_t_rvs_array(native,expected,shape=matrix,df=df(nu),status=info)
            case(3)
                call multivariate_t_rvs_array(native,expected,mu,matrix,df(nu),info)
            end select
            call scifort_multivariate_t_rvs_array_f64(2_c_int,5_c_int,serialized,actual, &
                has_mu,mu,has_matrix,matrix,df(nu),status)
            if (status/=info .or. status/=0 .or. any(actual/=expected)) stop 3
            call rng_get_state(native,expected_state)
            if (any(serialized/=expected_state)) stop 4
        end do
    end do
    print '(a)', 'Multivariate normal/t native/C sampling checks passed.'
end program test_multivariate_sampling_c_api
