program test_multiscale_graphcorr
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use multiscale_graphcorr_reference
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state
    use scifort_stats, only : mgc_result, multiscale_graphcorr, &
        mgc_status_invalid_parameter, mgc_status_invalid_shape, mgc_status_nonfinite, &
        mgc_status_success
    implicit none

    real(dp), parameter :: x_quad(6) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: y_quad(6) = [0.0_dp, 1.0_dp, 4.0_dp, 9.0_dp, 16.0_dp, 25.0_dp]
    real(dp), parameter :: x_mv(7, 2) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [7, 2])
    real(dp), parameter :: y_mv(7) = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp]
    real(dp), parameter :: y_custom(7, 2) = reshape([ &
        0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [7, 2])
    real(dp), parameter :: x_tied(6) = [0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp]
    real(dp), parameter :: y_tied(6) = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: u(5, 1) = reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [5, 1])
    real(dp), parameter :: v(6, 1) = reshape([5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp], [6, 1])

    type(mgc_result) :: r
    type(mgc_result) :: r2
    type(rng_state) :: state
    type(rng_state) :: state2
    real(dp) :: dx(6, 6)
    real(dp) :: dy(6, 6)
    real(dp) :: xbad(6, 2)
    real(dp) :: ybad(7, 3)
    real(dp) :: qnan
    integer :: b
    integer :: failures
    integer :: i
    integer :: j

    failures = 0
    qnan = ieee_value(0.0_dp, ieee_quiet_nan)

    call rng_seed(state, 1001)
    r = multiscale_graphcorr(x_quad, y_quad, state, reps=0)
    call check_int('quad status', r%status, mgc_status_success, failures)
    call check_close('quad statistic', r%statistic, quad_stat, 2.0e-15_dp, failures)
    call check_close('quad alias', r%stat, r%statistic, 0.0_dp, failures)
    call check_int_array('quad scale', r%opt_scale, quad_scale, failures)
    call check_matrix('quad map', r%mgc_map, quad_map, 8.0e-14_dp, failures)
    call check_close('zero reps pvalue', r%pvalue, 1.0_dp, 0.0_dp, failures)

    call rng_seed(state, 10011)
    r = multiscale_graphcorr(x_quad, y_quad, state, reps=0, compute_distance=cityblock_distance)
    call check_close('vector custom distance', r%statistic, quad_stat, 2.0e-15_dp, failures)

    call rng_seed(state, 1002)
    r = multiscale_graphcorr(x_mv, y_mv, state, reps=0)
    call check_close('matrix-vector statistic', r%statistic, mv_stat, 8.0e-15_dp, failures)
    call check_int_array('matrix-vector scale', r%opt_scale, mv_scale, failures)
    call check_matrix('matrix-vector map', r%mgc_map, mv_map, 8.0e-14_dp, failures)

    call rng_seed(state, 10021)
    r = multiscale_graphcorr(x_mv, y_custom, state, reps=0, compute_distance=cityblock_distance)
    call check_close('custom-distance statistic', r%statistic, custom_stat, 8.0e-15_dp, failures)
    call check_int_array('custom-distance scale', r%opt_scale, custom_scale, failures)
    call check_matrix('custom-distance map', r%mgc_map, custom_map, 8.0e-14_dp, failures)

    do j = 1, 6
        do i = 1, 6
            dx(i, j) = abs(x_tied(i) - x_tied(j))
            dy(i, j) = abs(y_tied(i) - y_tied(j))
        end do
    end do
    call rng_seed(state, 1003)
    r = multiscale_graphcorr(dx, dy, state, reps=0, precomputed=.true.)
    call check_close('precomputed statistic', r%statistic, pre_stat, 8.0e-15_dp, failures)
    call check_int_array('precomputed scale', r%opt_scale, pre_scale, failures)
    call check_matrix('precomputed map', r%mgc_map, pre_map, 8.0e-14_dp, failures)

    call rng_seed(state, 1004)
    r = multiscale_graphcorr(u, v, state, reps=0)
    call check_close('two-sample statistic', r%statistic, two_stat, 8.0e-15_dp, failures)
    call check_int_array('two-sample scale', r%opt_scale, two_scale, failures)
    call check_matrix('two-sample map', r%mgc_map, two_map, 8.0e-14_dp, failures)

    call rng_seed(state, 24680)
    r = multiscale_graphcorr(x_quad, y_quad, state, reps=23)
    call rng_seed(state2, 24680)
    r2 = multiscale_graphcorr(x_quad, y_quad, state2, reps=23)
    call check_vector('permutation replay', r2%null_dist, r%null_dist, 0.0_dp, failures)
    call check_close('permutation replay pvalue', r2%pvalue, r%pvalue, 0.0_dp, failures)
    b = count(r%null_dist >= r%statistic)
    call check_close('permutation plus one', r%pvalue, &
        real(b + 1, dp) / real(size(r%null_dist) + 1, dp), 0.0_dp, failures)

    call rng_seed(state, 1005)
    r = multiscale_graphcorr(x_quad(1:4), y_quad(1:4), state, reps=0)
    call check_int('low sample status', r%status, mgc_status_invalid_shape, failures)

    xbad = 0.0_dp
    ybad = 0.0_dp
    call rng_seed(state, 1006)
    r = multiscale_graphcorr(xbad, ybad, state, reps=0)
    call check_int('mismatched dimensions status', r%status, mgc_status_invalid_shape, failures)

    xbad = 0.0_dp
    xbad(1, 1) = qnan
    call rng_seed(state, 1007)
    r = multiscale_graphcorr(xbad, xbad, state, reps=0)
    call check_int('nonfinite status', r%status, mgc_status_nonfinite, failures)

    call rng_seed(state, 1008)
    r = multiscale_graphcorr(dx, dy, state, reps=-1, precomputed=.true.)
    call check_int('negative reps status', r%status, mgc_status_invalid_parameter, failures)

    call rng_seed(state, 1009)
    r = multiscale_graphcorr(dx, dy, state, reps=0, precomputed=.true., is_twosamp=.true.)
    call check_int('precomputed two-sample status', r%status, mgc_status_invalid_parameter, failures)

    if (failures /= 0) then
        print *, 'multiscale_graphcorr failures:', failures
        error stop 1
    end if
    print *, 'multiscale_graphcorr tests passed'

contains

    subroutine cityblock_distance(data, distance, status)
        real(dp), intent(in) :: data(:, :)
        real(dp), allocatable, intent(out) :: distance(:, :)
        integer, intent(out) :: status

        integer :: ii
        integer :: jj

        allocate(distance(size(data, 1), size(data, 1)))
        distance = 0.0_dp
        do jj = 2, size(data, 1)
            do ii = 1, jj - 1
                distance(ii, jj) = sum(abs(data(ii, :) - data(jj, :)))
                distance(jj, ii) = distance(ii, jj)
            end do
        end do
        status = 0
    end subroutine cityblock_distance

    subroutine check_close(name, actual, expected, tolerance, count_failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance
        integer, intent(inout) :: count_failures

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', trim(name), actual, expected
            count_failures = count_failures + 1
        end if
    end subroutine check_close

    subroutine check_int(name, actual, expected, count_failures)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual
        integer, intent(in) :: expected
        integer, intent(inout) :: count_failures

        if (actual /= expected) then
            print *, 'FAIL ', trim(name), actual, expected
            count_failures = count_failures + 1
        end if
    end subroutine check_int

    subroutine check_int_array(name, actual, expected, count_failures)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual(:)
        integer, intent(in) :: expected(:)
        integer, intent(inout) :: count_failures

        if (size(actual) /= size(expected) .or. any(actual /= expected)) then
            print *, 'FAIL ', trim(name), actual, expected
            count_failures = count_failures + 1
        end if
    end subroutine check_int_array

    subroutine check_vector(name, actual, expected, tolerance, count_failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        real(dp), intent(in) :: tolerance
        integer, intent(inout) :: count_failures

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(name), ' size'
            count_failures = count_failures + 1
        else if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', trim(name)
            count_failures = count_failures + 1
        end if
    end subroutine check_vector

    subroutine check_matrix(name, actual, expected, tolerance, count_failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :)
        real(dp), intent(in) :: expected(:, :)
        real(dp), intent(in) :: tolerance
        integer, intent(inout) :: count_failures

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', trim(name), ' shape', shape(actual), shape(expected)
            count_failures = count_failures + 1
        else if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', trim(name), maxval(abs(actual - expected))
            count_failures = count_failures + 1
        end if
    end subroutine check_matrix

end program test_multiscale_graphcorr
