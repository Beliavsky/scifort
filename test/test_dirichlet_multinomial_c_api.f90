! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_dirichlet_multinomial_c_api
    use iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_dirichlet_multinomial
    use scifort_c_api_generated_dirichlet_multinomial
    implicit none
    real(c_double) :: x(2), alpha(2), value, expected, v(2), cov(2,2)
    integer(c_int) :: n, status
    integer :: i

    x = [1.0_c_double, 3.0_c_double]
    do i = 1, 4
        alpha = [2.0_c_double, 3.0_c_double]
        n = 4
        if (i == 2) n = -1
        if (i == 3) alpha(1) = 0.0_c_double
        if (i == 4) n = 0
        expected = dirichlet_multinomial_pmf(x, alpha, n)
        call scifort_dirichlet_multinomial_pmf_f64(2_c_int, x, alpha, n, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 1
        expected = dirichlet_multinomial_logpmf(x, alpha, n)
        call scifort_dirichlet_multinomial_logpmf_f64(2_c_int, x, alpha, n, value, status)
        if (status /= 0 .or. .not. equal(value, expected)) stop 2
        call scifort_dirichlet_multinomial_mean_f64(2_c_int, alpha, n, v, status)
        if (status /= 0 .or. .not. all(equal(v, dirichlet_multinomial_mean(alpha,n)))) stop 3
        call scifort_dirichlet_multinomial_var_f64(2_c_int, alpha, n, v, status)
        if (status /= 0 .or. .not. all(equal(v, dirichlet_multinomial_var(alpha,n)))) stop 4
        call scifort_dirichlet_multinomial_cov_f64(2_c_int, alpha, n, cov, status)
        if (status /= 0 .or. .not. all(equal(cov, dirichlet_multinomial_cov(alpha,n)))) stop 5
    end do
    print '(a)', 'Dirichlet-multinomial native/C agreement passed.'
contains
    pure elemental function equal(a, b) result(ok)
        implicit none
        real(c_double), intent(in) :: a !! actual value
        real(c_double), intent(in) :: b !! expected value
        logical :: ok
        ok = (a == b) .or. (ieee_is_nan(a) .and. ieee_is_nan(b))
    end function equal
end program test_dirichlet_multinomial_c_api
