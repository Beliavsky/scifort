! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_multivariate_parameters_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_nan
    use scifort_multivariate_normal
    use scifort_multivariate_t
    use scifort_c_api_generated_multivariate_normal
    use scifort_c_api_generated_multivariate_t
    implicit none
    real(c_double) :: center(3), matrix(3,3), expected_center(2), actual_center(2)
    real(c_double) :: expected_matrix(2,2), actual_matrix(2,2), value, expected, df(4)
    integer(c_int) :: indices(2), status
    integer :: native_status, j
    center=[1.0_c_double,2.0_c_double,3.0_c_double]
    matrix=reshape([4.0_c_double,1.0_c_double,.5_c_double, &
        91.0_c_double,3.0_c_double,.2_c_double, &
        92.0_c_double,93.0_c_double,2.0_c_double],[3,3])
    expected=multivariate_normal_entropy(center,matrix)
    call scifort_multivariate_normal_entropy_f64(3_c_int,center,matrix,value,status)
    if (status/=0 .or. ieee_is_nan(value) .or. value/=expected) stop 1
    df=[1.0_c_double,5.0_c_double,1.0e6_c_double, &
        ieee_value(0.0_c_double,ieee_positive_inf)]
    do j=1,size(df)
        expected=multivariate_t_entropy(center,matrix,df(j))
        call scifort_multivariate_t_entropy_f64(3_c_int,center,matrix,df(j),value,status)
        if (status/=0 .or. ieee_is_nan(value) .or. value/=expected) stop 2
    end do
    indices=[3_c_int,1_c_int]
    call multivariate_normal_marginal(indices,center,matrix,expected_center,expected_matrix,native_status)
    call scifort_multivariate_normal_marginal_f64(3_c_int,2_c_int,indices,center,matrix, &
        actual_center,actual_matrix,status)
    if (status/=native_status .or. status/=0) stop 3
    if (any(actual_center/=expected_center) .or. any(actual_matrix/=expected_matrix)) stop 4
    call multivariate_t_marginal(indices,center,matrix,expected_center,expected_matrix,native_status)
    call scifort_multivariate_t_marginal_f64(3_c_int,2_c_int,indices,center,matrix, &
        actual_center,actual_matrix,status)
    if (status/=native_status .or. status/=0) stop 5
    if (any(actual_center/=expected_center) .or. any(actual_matrix/=expected_matrix)) stop 6
    actual_center=77.0_c_double
    actual_matrix=88.0_c_double
    indices=[1_c_int,1_c_int]
    call scifort_multivariate_normal_marginal_f64(3_c_int,2_c_int,indices,center,matrix, &
        actual_center,actual_matrix,status)
    if (status/=2 .or. any(actual_center/=77.0_c_double) .or. &
        any(actual_matrix/=88.0_c_double)) stop 7
    indices=[0_c_int,1_c_int]
    call scifort_multivariate_t_marginal_f64(3_c_int,2_c_int,indices,center,matrix, &
        actual_center,actual_matrix,status)
    if (status/=2 .or. any(actual_center/=77.0_c_double) .or. &
        any(actual_matrix/=88.0_c_double)) stop 8
    value=99.0_c_double
    call scifort_multivariate_normal_entropy_f64(0_c_int,center,matrix,value,status)
    if (status/=1 .or. value/=99.0_c_double) stop 9
    call scifort_multivariate_t_entropy_f64(0_c_int,center,matrix,5.0_c_double,value,status)
    if (status/=1 .or. value/=99.0_c_double) stop 10
    print '(a)', 'Multivariate entropy/marginal native/C checks passed.'
end program test_multivariate_parameters_c_api
