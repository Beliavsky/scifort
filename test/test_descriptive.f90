! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_descriptive
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_stats, only : dp, central_moment, covariance, mean, median, &
        pearson_correlation, quantile, rankdata, standard_deviation, variance
    implicit none

    real(dp), parameter :: x(6) = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
    real(dp), parameter :: y(6) = [3.0_dp, -1.0_dp, 5.0_dp, 4.0_dp, 9.0_dp, 12.0_dp]
    real(dp), parameter :: q(7) = [0.0_dp, 0.1_dp, 0.25_dp, 0.5_dp, &
        0.75_dp, 0.9_dp, 1.0_dp]
    real(dp), parameter :: q_expected(7) = [1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp, &
        7.0_dp, 12.0_dp, 16.0_dp]
    real(dp), parameter :: rank_average(6) = [1.0_dp, 2.5_dp, 2.5_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    real(dp), parameter :: rank_min(6) = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    real(dp), parameter :: rank_max(6) = [1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    real(dp), parameter :: rank_dense(6) = [1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: rank_ordinal(6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]

    real(dp) :: empty(0)
    real(dp) :: extreme(3)
    real(dp) :: nan_data(2)
    real(dp) :: inf_data(3)
    real(dp), allocatable :: values(:)

    call check_close('mean', mean(x), 5.5_dp, 2.0e-15_dp)
    call check_close('variance ddof=1', variance(x), 32.7_dp, 5.0e-14_dp)
    call check_close('variance ddof=0', variance(x, 0), 27.25_dp, 5.0e-14_dp)
    call check_close('standard deviation', standard_deviation(x), &
        5.718391382198319_dp, 2.0e-14_dp)
    call check_close('moment 0', central_moment(x, 0), 1.0_dp, 0.0_dp)
    call check_close('moment 1', central_moment(x, 1), 0.0_dp, 0.0_dp)
    call check_close('moment 2', central_moment(x, 2), 27.25_dp, 5.0e-14_dp)
    call check_close('moment 3', central_moment(x, 3), 165.5_dp, 3.0e-13_dp)
    call check_close('moment 4', central_moment(x, 4), 2151.5625_dp, 5.0e-12_dp)

    values = quantile(x, q)
    call check_array('quantile array', values, q_expected, 3.0e-15_dp)
    call check_close('quantile scalar', quantile(x, 0.1_dp), 1.5_dp, 3.0e-15_dp)
    call check_close('median', median(x), 3.0_dp, 3.0e-15_dp)

    call check_close('covariance ddof=1', covariance(x, y), 23.0_dp, 5.0e-14_dp)
    call check_close('covariance ddof=0', covariance(x, y, 0), &
        19.166666666666664_dp, 5.0e-14_dp)
    call check_close('pearson correlation', pearson_correlation(x, y), &
        0.8763065299100674_dp, 3.0e-15_dp)

    values = rankdata(x)
    call check_array('rank average', values, rank_average, 0.0_dp)
    values = rankdata(x, 'min')
    call check_array('rank min', values, rank_min, 0.0_dp)
    values = rankdata(x, 'max')
    call check_array('rank max', values, rank_max, 0.0_dp)
    values = rankdata(x, 'dense')
    call check_array('rank dense', values, rank_dense, 0.0_dp)
    values = rankdata(x, 'ordinal')
    call check_array('rank ordinal', values, rank_ordinal, 0.0_dp)

    extreme = [1.0e308_dp, -1.0e308_dp, 1.0_dp]
    call check_close('scaled extreme mean', mean(extreme), 1.0_dp / 3.0_dp, 1.0e-15_dp)
    call check_close('scaled extreme sd', standard_deviation(extreme), 1.0e308_dp, 3.0e292_dp)
    call check_true('scaled extreme variance', variance(extreme) > huge(1.0_dp))

    nan_data = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)]
    call check_nan('mean nan', mean(nan_data))
    call check_nan('variance nan', variance(nan_data))
    call check_nan('quantile nan', quantile(nan_data, 0.5_dp))
    values = rankdata(nan_data)
    call check_true('rank nan propagation', all(ieee_is_nan(values)))

    call check_nan('empty mean', mean(empty))
    call check_nan('empty variance', variance(empty))
    call check_nan('invalid ddof', variance(x, size(x)))
    call check_nan('invalid moment order', central_moment(x, -1))
    call check_nan('invalid quantile', quantile(x, -0.1_dp))
    values = rankdata(x, 'bad')
    call check_true('invalid rank method', all(ieee_is_nan(values)))

    inf_data = [1.0_dp, ieee_value(0.0_dp, ieee_positive_inf), 4.0_dp]
    call check_true('positive infinity mean', mean(inf_data) > huge(1.0_dp))
    call check_true('infinite upper quantile', quantile(inf_data, 0.9_dp) > huge(1.0_dp))

    print *, 'descriptive statistics tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_nan(name, actual)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value expected to be NaN

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', name, actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive
