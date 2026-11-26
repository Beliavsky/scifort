! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_dirichlet_multinomial_counts
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, &
        ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use scifort_kinds, only : dp
    use scifort_dirichlet_multinomial, only : dirichlet_multinomial_pmf, &
        dirichlet_multinomial_logpmf
    implicit none
    real(dp) :: bad(6), x(2), alpha(2), y
    integer :: i

    alpha = [2.0_dp, 3.0_dp]
    bad = [-1.0_dp, 0.5_dp, ieee_value(0.0_dp, ieee_quiet_nan), &
        ieee_value(0.0_dp, ieee_positive_inf), ieee_value(0.0_dp, ieee_negative_inf), &
        -1.0e100_dp]
    do i = 1, size(bad)
        x = [bad(i), 0.0_dp]
        if (.not. ieee_is_nan(dirichlet_multinomial_pmf(x, alpha, 0))) stop 1
        if (.not. ieee_is_nan(dirichlet_multinomial_logpmf(x, alpha, 0))) stop 2
    end do
    ! These integer-valued reals exceed default-integer range and lie outside support.
    x = [real(huge(0), dp) + 1.0_dp, 0.0_dp]
    if (dirichlet_multinomial_pmf(x, alpha, 0) /= 0.0_dp) stop 3
    x = [1.0e100_dp, 0.0_dp]
    if (dirichlet_multinomial_pmf(x, alpha, 0) /= 0.0_dp) stop 4
    y = dirichlet_multinomial_logpmf(x, alpha, 0)
    if (y /= ieee_value(0.0_dp, ieee_negative_inf)) stop 5
    ! With one category all trials are certain, including the largest native n.
    y = dirichlet_multinomial_pmf([real(huge(0), dp)], [1.0_dp], huge(0))
    if (abs(y - 1.0_dp) > 1.0e-12_dp) stop 6
    print '(a)', 'Dirichlet-multinomial count regressions passed.'
end program test_dirichlet_multinomial_counts
