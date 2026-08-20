! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program scifort_demo
    use scifort_stats, only : cauchy_cdf, dp, exponential_cdf, laplace_cdf, &
        logistic_cdf, normal_cdf, normal_pdf, normal_ppf, uniform_cdf
    implicit none

    integer :: i
    real(dp) :: x(5)

    x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]

    print '(a)', 'SciFort distribution demonstration'
    print '(a)', '      x   normal pdf   normal cdf      uniform  exponential' // &
        '      laplace     logistic       cauchy'
    do i = 1, size(x)
        print '(f7.2,7(1x,f12.8))', x(i), normal_pdf(x(i)), normal_cdf(x(i)), &
            uniform_cdf(x(i), loc=-2.0_dp, scale=4.0_dp), &
            exponential_cdf(x(i), loc=-2.0_dp), laplace_cdf(x(i)), &
            logistic_cdf(x(i)), cauchy_cdf(x(i))
    end do

    print '(a)'
    print '(a,f20.15)', 'normal_ppf(0.975) = ', normal_ppf(0.975_dp)
end program scifort_demo
