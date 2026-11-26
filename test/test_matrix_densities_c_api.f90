! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_matrix_densities_c_api
    use iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_matrix_normal
    use scifort_matrix_t
    use scifort_c_api_generated_matrix_normal
    use scifort_c_api_generated_matrix_t
    use test_matrix_distributions_reference
    implicit none
    real(c_double) :: value, expected, identity(1,1)
    integer(c_int) :: status
    expected = matrix_normal_pdf(matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref)
    call scifort_matrix_normal_pdf_f64( &
        2_c_int, 2_c_int, matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 1
    if (abs(value-matrix_normal_pdf_ref) > 1.0e-13_c_double) stop 2
    expected = matrix_normal_logpdf(matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref)
    call scifort_matrix_normal_logpdf_f64( &
        2_c_int, 2_c_int, matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 1
    if (abs(value-matrix_normal_logpdf_ref) > 1.0e-13_c_double) stop 2
    expected = matrix_normal_entropy(matrix_row_ref, matrix_col_ref)
    call scifort_matrix_normal_entropy_f64( &
        2_c_int, 2_c_int, matrix_row_ref, matrix_col_ref, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 1
    if (abs(value-matrix_normal_entropy_ref) > 1.0e-13_c_double) stop 2
    expected = matrix_t_pdf(matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref)
    call scifort_matrix_t_pdf_f64( &
        2_c_int, 2_c_int, matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 1
    if (abs(value-matrix_t_pdf_ref) > 1.0e-13_c_double) stop 2
    expected = matrix_t_logpdf(matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref)
    call scifort_matrix_t_logpdf_f64( &
        2_c_int, 2_c_int, matrix_x_ref, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 1
    if (abs(value-matrix_t_logpdf_ref) > 1.0e-13_c_double) stop 2
    identity = 1.0_c_double
    expected = matrix_normal_entropy()
    call scifort_matrix_normal_entropy_f64( &
        1_c_int, 1_c_int, identity, identity, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 3
    if (abs(value-expected) > 1.0e-13_c_double) stop 3
    expected = matrix_normal_entropy(rowcov=matrix_row_ref)
    call scifort_matrix_normal_entropy_f64( &
        2_c_int, 1_c_int, matrix_row_ref, identity, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 3
    if (abs(value-expected) > 1.0e-13_c_double) stop 3
    expected = matrix_normal_entropy(colcov=matrix_col_ref)
    call scifort_matrix_normal_entropy_f64( &
        1_c_int, 2_c_int, identity, matrix_col_ref, value, status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 3
    if (abs(value-expected) > 1.0e-13_c_double) stop 3
    print '(a)', 'Matrix distribution native/C and saved reference checks passed.'
end program test_matrix_densities_c_api
