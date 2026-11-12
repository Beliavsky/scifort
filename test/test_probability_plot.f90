! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_probability_plot
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_stats, only : probplot, ppcc_max, ppcc_plot, boxcox_normplot, &
        yeojohnson_normplot
    implicit none

    real(dp), parameter :: x(6) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp]
    real(dp), parameter :: wx(8) = [0.5_dp, 1.0_dp, 1.7_dp, 2.8_dp, &
        4.4_dp, 6.3_dp, 8.9_dp, 12.4_dp]
    real(dp), parameter :: y(6) = [-3.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 5.0_dp]
    real(dp), allocatable :: a(:)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: empty(:)
    real(dp) :: intercept
    real(dp) :: r
    real(dp) :: shape
    real(dp) :: slope

    call probplot(x, a, b, slope, intercept, r)
    call check_close(a(1), -1.2313217056520_dp, 2.0e-8_dp, 'normal osm first')
    call check_close(a(6), 1.2313217056520_dp, 2.0e-8_dp, 'normal osm last')
    call check_close(slope, 12.237445086425888_dp, 3.0e-8_dp, 'normal slope')
    call check_close(intercept, 10.5_dp, 2.0e-12_dp, 'normal intercept')
    call check_close(r, 0.9117077571877801_dp, 3.0e-9_dp, 'normal correlation')

    call probplot(x, a, b, slope, intercept, r, dist='t', sparams=[5.0_dp])
    call check_close(a(1), -1.4078576989448_dp, 3.0e-8_dp, 'student t osm first')
    call check_close(r, 0.9126967756213079_dp, 5.0e-10_dp, 'student t correlation')

    call probplot(x, a, b, slope, intercept, r, dist='loggamma', sparams=[2.5_dp])
    call check_close(r, 0.8744399770484781_dp, 3.0e-9_dp, 'loggamma correlation')

    call probplot(x, a, b, slope, intercept, r, sparams=[2.0_dp, 3.0_dp], ppf=affine_ppf)
    call check_close(a(1), 2.0_dp + 3.0_dp * 0.1091012818597_dp, 3.0e-12_dp, &
        'custom ppf callback')

    call ppcc_plot(x, -1.0_dp, 1.0_dp, a, b, n=5)
    call check_close(a(1), -1.0_dp, 0.0_dp, 'ppcc grid first')
    call check_close(a(5), 1.0_dp, 0.0_dp, 'ppcc grid last')
    call check_close(b(1), 0.9045088656683915_dp, 3.0e-8_dp, 'tukey ppcc first')
    call check_close(b(3), 0.9125351526698542_dp, 3.0e-8_dp, 'tukey ppcc middle')

    shape = ppcc_max(x)
    call check_close(shape, -0.2178534880812422_dp, 3.0e-6_dp, 'default ppcc maximum')
    shape = ppcc_max(x, brack=[-1.0_dp, 1.0_dp])
    call check_close(shape, 4.049744384396333_dp, 3.0e-6_dp, 'expanded ppcc maximum')

    call ppcc_plot(wx, 0.5_dp, 4.0_dp, a, b, dist='weibull_min', n=6)
    call check_close(b(2), 0.9970484405891257_dp, 3.0e-9_dp, 'weibull ppcc')
    shape = ppcc_max(wx, brack=[0.5_dp, 4.0_dp], dist='weibull_min')
    call check_close(shape, 1.132925174032626_dp, 8.0e-7_dp, 'weibull ppcc maximum')

    call boxcox_normplot(x, -1.0_dp, 1.0_dp, a, b, n=5)
    call check_close(b(1), 0.9117077571877801_dp, 3.0e-9_dp, 'boxcox normplot first')
    call check_close(b(3), 0.9974010248218916_dp, 3.0e-9_dp, 'boxcox normplot middle')

    call yeojohnson_normplot(y, -1.0_dp, 2.0_dp, a, b, n=5)
    call check_close(b(1), 0.7575916429006504_dp, 1.0e-8_dp, 'yeojohnson normplot first')
    call check_close(b(3), 0.9635661308734661_dp, 3.0e-9_dp, 'yeojohnson normplot middle')

    allocate(empty(0))
    call ppcc_plot(empty, -1.0_dp, 1.0_dp, a, b, n=3)
    if (size(a) /= 3 .or. any(abs(b) > 0.0_dp)) error stop 'empty ppcc mismatch'
    call probplot(empty, a, b, slope, intercept, r)
    if (size(a) /= 0 .or. size(b) /= 0) error stop 'empty probplot shapes'
    if (.not. ieee_is_nan(slope) .or. .not. ieee_is_nan(intercept)) then
        error stop 'empty probplot fit values'
    end if
    if (abs(r) > 0.0_dp) error stop 'empty probplot correlation'

    print *, 'probability-plot tests passed'

contains

    function affine_ppf(p, params) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: params(:)
        real(dp) :: value

        value = params(1) + params(2) * p
    end function affine_ppf

    subroutine check_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: atol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > atol) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' abs err = ', abs(actual - expected)
            error stop 1
        end if
    end subroutine check_close

end program test_probability_plot
