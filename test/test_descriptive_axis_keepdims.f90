! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_descriptive_axis_keepdims
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, hmean_keepdims, &
        iqr, iqr_keepdims, kurtosis, kurtosis_keepdims, median_abs_deviation, &
        median_abs_deviation_keepdims, moment, moment_keepdims, pmean, pmean_keepdims, &
        sem, sem_keepdims, skew, skew_keepdims, variation, variation_keepdims
    implicit none

    real(dp), parameter :: x(3, 4) = reshape([ &
        1.0_dp, 3.0_dp, 2.0_dp, &
        2.0_dp, 5.0_dp, 6.0_dp, &
        4.0_dp, 7.0_dp, 10.0_dp, &
        8.0_dp, 9.0_dp, 14.0_dp], [3, 4])
    real(dp), parameter :: w4(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: w31(3, 1) = reshape([1.0_dp, 2.0_dp, 4.0_dp], [3, 1])
    real(dp), parameter :: bad_weights(3) = [1.0_dp, 2.0_dp, 3.0_dp]

    real(dp), allocatable :: matrix_values(:, :)
    real(dp), allocatable :: values(:)
    real(dp) :: x_nan(3, 4)

    call check_array('gmean default axis', gmean(x), [ &
        1.8171205928321394_dp, 3.9148676411688634_dp, 6.542132620377179_dp, &
        10.026595869929169_dp], 3.0e-14_dp)
    call check_array('gmean axis 2', gmean(x, axis=2), [ &
        2.828427124746190_dp, 5.544443371329423_dp, 6.402171745887360_dp], 2.0e-14_dp)
    call check_array('gmean vector broadcast', gmean(x, weights=w4), [ &
        1.8171205928321394_dp, 3.9148676411688634_dp, 6.542132620377179_dp, &
        10.026595869929169_dp], 3.0e-14_dp)
    call check_array('gmean axis 2 weighted', gmean(x, axis=2, weights=w4), [ &
        3.999999999999999_dp, 6.648635774992930_dp, 8.793959555495958_dp], 3.0e-14_dp)
    call check_array('gmean matrix broadcast', gmean(x, weights=w31), [ &
        2.0339370097944314_dp, 4.8681937895234295_dp, 7.923074224511346_dp, &
        11.391577742394356_dp], 4.0e-14_dp)

    call check_array('hmean default axis', hmean(x), [ &
        1.6363636363636365_dp, 3.4615384615384617_dp, 6.086956521739131_dp, &
        9.754838709677419_dp], 3.0e-14_dp)
    call check_array('hmean axis 2 weighted', hmean(x, axis=2, weights=w4), [ &
        3.0769230769230766_dp, 6.225296442687747_dp, 7.046979865771813_dp], 3.0e-14_dp)
    call check_array('pmean default axis', pmean(x, 2.0_dp), [ &
        2.160246899469287_dp, 4.654746681256314_dp, 7.416198487095663_dp, &
        10.661457061146317_dp], 3.0e-14_dp)
    call check_array('pmean axis 2 weighted', pmean(x, 2.0_dp, axis=2, weights=w4), [ &
        5.594640292279746_dp, 7.280109889280518_dp, 10.770329614269007_dp], 3.0e-14_dp)

    call check_array('sem default axis', sem(x), [ &
        0.5773502691896258_dp, 1.2018504251546631_dp, 1.7320508075688774_dp, &
        1.855921454276674_dp], 2.0e-14_dp)
    call check_array('sem axis 2 ddof 0', sem(x, axis=2, ddof=0), [ &
        1.340475661845451_dp, 1.118033988749895_dp, 2.236067977499790_dp], 2.0e-14_dp)
    call check_array('variation default axis', variation(x), [ &
        0.4082482904638630_dp, 0.3922322702763681_dp, 0.34992710611188255_dp, &
        0.25400025400038095_dp], 2.0e-14_dp)
    call check_array('variation axis 2 ddof 1', variation(x, axis=2, ddof=1), [ &
        0.8255189164891871_dp, 0.4303314829119352_dp, 0.6454972243679028_dp], 2.0e-14_dp)

    call check_close('iqr all axes', iqr(x), 5.5_dp, 0.0_dp)
    call check_array('iqr axis 1', iqr(x, 1), [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp], 0.0_dp)
    call check_array('iqr axis 2 range scale', &
        iqr(x, 2, [10.0_dp, 90.0_dp], 2.0_dp), [ &
        2.7500000000000004_dp, 2.4000000000000004_dp, 4.800000000000001_dp], 2.0e-14_dp)
    call check_array('mad default axis', median_abs_deviation(x), &
        [1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp], 0.0_dp)
    call check_array('mad axis 2 scale', median_abs_deviation(x, axis=2, scale=1.4826_dp), [ &
        1.011736139214893_dp, 1.3489815189531904_dp, 2.697963037906381_dp], 2.0e-14_dp)

    call check_array('moment default axis', moment(x, order=2), [ &
        0.6666666666666666_dp, 2.888888888888889_dp, 6.0_dp, &
        6.888888888888888_dp], 2.0e-14_dp)
    call check_array('moment axis 2 center 0', moment(x, order=2, axis=2, center=0.0_dp), &
        [21.25_dp, 41.0_dp, 84.0_dp], 2.0e-14_dp)
    call check_array('skew default axis', skew(x), [ &
        0.0_dp, -0.5280049792181872_dp, 0.0_dp, 0.6309038567106232_dp], 3.0e-14_dp)
    call check_array('skew axis 2 unbiased', skew(x, axis=2, bias=.false.), [ &
        1.1376243669576889_dp, 0.0_dp, 0.0_dp], 3.0e-14_dp)
    call check_array('kurtosis default axis', kurtosis(x), [ &
        -1.5_dp, -1.5000000000000004_dp, -1.5_dp, -1.5000000000000007_dp], 3.0e-14_dp)
    call check_array('kurtosis axis 2 pearson unbiased', &
        kurtosis(x, axis=2, fisher=.false., bias=.false.), [ &
        3.757655954631380_dp, 1.8000000000000007_dp, 1.8000000000000007_dp], 4.0e-14_dp)

    matrix_values = gmean_keepdims(x)
    call check_matrix('gmean keepdims axis 1', matrix_values, reshape([ &
        1.8171205928321394_dp, 3.9148676411688634_dp, 6.542132620377179_dp, &
        10.026595869929169_dp], [1, 4]), 3.0e-14_dp)
    matrix_values = gmean_keepdims(x, axis=2, weights=w4)
    call check_matrix('gmean keepdims axis 2 weighted', matrix_values, reshape([ &
        3.999999999999999_dp, 6.648635774992930_dp, 8.793959555495958_dp], [3, 1]), &
        3.0e-14_dp)
    matrix_values = gmean_keepdims(x, weights=w31)
    call check_matrix('gmean keepdims matrix broadcast', matrix_values, reshape([ &
        2.0339370097944314_dp, 4.8681937895234295_dp, 7.923074224511346_dp, &
        11.391577742394356_dp], [1, 4]), 4.0e-14_dp)
    call check_matrix('sem keepdims axis 2', sem_keepdims(x, axis=2), reshape([ &
        1.547847968417226_dp, 1.2909944487358056_dp, 2.581988897471611_dp], [3, 1]), &
        3.0e-14_dp)
    call check_matrix('iqr keepdims all', iqr_keepdims(x), reshape([5.5_dp], [1, 1]), 0.0_dp)
    call check_matrix('iqr keepdims axis 1', iqr_keepdims(x, axis=1), &
        reshape([1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp], [1, 4]), 0.0_dp)
    call check_matrix('mad keepdims axis 2', median_abs_deviation_keepdims(x, axis=2), &
        reshape([1.5_dp, 2.0_dp, 4.0_dp], [3, 1]), 0.0_dp)
    call check_matrix('moment keepdims axis 1', moment_keepdims(x, order=2), reshape([ &
        0.6666666666666666_dp, 2.888888888888889_dp, 6.0_dp, &
        6.888888888888888_dp], [1, 4]), 2.0e-14_dp)
    call check_matrix('skew keepdims axis 2', skew_keepdims(x, axis=2), reshape([ &
        0.6568077344996993_dp, 0.0_dp, 0.0_dp], [3, 1]), 3.0e-14_dp)
    call check_matrix('kurtosis keepdims axis 1', kurtosis_keepdims(x), reshape([ &
        -1.5_dp, -1.5000000000000004_dp, -1.5_dp, -1.5000000000000007_dp], [1, 4]), &
        3.0e-14_dp)

    call check_array('vector keepdims', gmean_keepdims(x(1, :)), [2.828427124746190_dp], &
        2.0e-14_dp)
    values = gmean(x, axis=3)
    call check_true('invalid matrix axis', size(values) == 0)
    values = gmean(x, weights=bad_weights)
    call check_true('invalid broadcast weights', size(values) == 4 .and. all(ieee_is_nan(values)))

    x_nan = x
    x_nan(2, 3) = ieee_value(0.0_dp, ieee_quiet_nan)
    values = gmean(x_nan)
    call check_true('slice nan propagation', ieee_is_nan(values(3)) .and. &
        .not. ieee_is_nan(values(1)) .and. .not. ieee_is_nan(values(2)) .and. &
        .not. ieee_is_nan(values(4)))

    print *, 'descriptive axis/keepdims tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! SciPy reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! SciPy reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (size(actual) /= size(expected) .or. &
                any(ieee_is_nan(actual)) .or. any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_matrix(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:, :) !! computed values
        real(dp), intent(in) :: expected(:, :) !! SciPy reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (any(shape(actual) /= shape(expected)) .or. &
                any(ieee_is_nan(actual)) .or. any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual shape  = ', shape(actual)
            print *, ' expected shape= ', shape(expected)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_matrix

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_axis_keepdims
