! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_wishart_c_api
    use iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_wishart
    use scifort_invwishart
    use scifort_c_api_generated_wishart
    use scifort_c_api_generated_invwishart
    implicit none
    real(c_double) :: x(2,2), scale(2,2), df, value, expected, mat(2,2)
    integer(c_int) :: status
    integer :: i
    x = reshape([1.5_c_double, .2_c_double, .2_c_double, 1.1_c_double], [2,2])
    do i = 1, 4
        scale = reshape([2.0_c_double, .4_c_double, .4_c_double, 1.0_c_double], [2,2])
        df = 7.0_c_double
        if (i == 2) df = 1.0_c_double
        if (i == 3) df = 3.0_c_double
        if (i == 4) scale(1,1) = -1.0_c_double
        expected = wishart_pdf(x, df, scale)
        call scifort_wishart_pdf_f64(2_c_int, x, df, scale, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        expected = wishart_logpdf(x, df, scale)
        call scifort_wishart_logpdf_f64(2_c_int, x, df, scale, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        expected = wishart_entropy(df, scale)
        call scifort_wishart_entropy_f64(2_c_int, df, scale, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        call scifort_wishart_mean_f64(2_c_int, df, scale, mat, status)
        if (status /= 0 .or. .not. all(equal(mat, wishart_mean(df, scale)))) stop 1
        call scifort_wishart_mode_f64(2_c_int, df, scale, mat, status)
        if (status /= 0 .or. .not. all(equal(mat, wishart_mode(df, scale)))) stop 1
        call scifort_wishart_var_f64(2_c_int, df, scale, mat, status)
        if (status /= 0 .or. .not. all(equal(mat, wishart_var(df, scale)))) stop 1
        expected = invwishart_pdf(x, df, scale)
        call scifort_invwishart_pdf_f64(2_c_int, x, df, scale, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        expected = invwishart_logpdf(x, df, scale)
        call scifort_invwishart_logpdf_f64(2_c_int, x, df, scale, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        expected = invwishart_entropy(df, scale)
        call scifort_invwishart_entropy_f64(2_c_int, df, scale, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        call scifort_invwishart_mean_f64(2_c_int, df, scale, mat, status)
        if (status /= 0 .or. .not. all(equal(mat, invwishart_mean(df, scale)))) stop 1
        call scifort_invwishart_mode_f64(2_c_int, df, scale, mat, status)
        if (status /= 0 .or. .not. all(equal(mat, invwishart_mode(df, scale)))) stop 1
        call scifort_invwishart_var_f64(2_c_int, df, scale, mat, status)
        if (status /= 0 .or. .not. all(equal(mat, invwishart_var(df, scale)))) stop 1
    end do
    print '(a)', 'Wishart and inverse-Wishart native/C agreement passed.'
contains
    pure elemental function equal(a, b) result(ok)
        implicit none
        real(c_double), intent(in) :: a !! actual value
        real(c_double), intent(in) :: b !! native reference value
        logical :: ok
        ok = (a == b) .or. (ieee_is_nan(a) .and. ieee_is_nan(b))
    end function equal
end program test_wishart_c_api
