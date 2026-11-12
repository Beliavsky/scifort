! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_nan_policy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, iqr, iqr_keepdims, &
        kurtosis, median_abs_deviation, moment, pmean, sem, sem_keepdims, skew, variation
    implicit none

    real(dp), parameter :: tol = 2.0e-14_dp
    real(dp) :: nan
    real(dp) :: x(4)
    real(dp) :: x_no_nan(4)
    real(dp) :: weights(4)
    real(dp) :: weights_nan(4)
    real(dp) :: matrix(3, 4)
    real(dp) :: weights31(3, 1)
    real(dp) :: weights31_nan(3, 1)
    real(dp) :: all_nan(3)
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: kept(:, :)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [1.0_dp, nan, 4.0_dp, 9.0_dp]
    x_no_nan = [1.0_dp, 2.0_dp, 4.0_dp, 9.0_dp]
    weights = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    weights_nan = [1.0_dp, nan, 3.0_dp, 4.0_dp]
    all_nan = nan
    matrix = reshape([1.0_dp, 3.0_dp, nan, nan, 5.0_dp, 6.0_dp, &
        4.0_dp, nan, 10.0_dp, 8.0_dp, 9.0_dp, 14.0_dp], [3, 4])
    weights31(:, 1) = [1.0_dp, 2.0_dp, 4.0_dp]
    weights31_nan(:, 1) = [1.0_dp, nan, 4.0_dp]

    call check_close('gmean omit', gmean(x, weights=weights, nan_policy='omit'), &
        5.045378491522287_dp, tol)
    call check_close('gmean weight NaN omit', &
        gmean(x_no_nan, weights=weights_nan, nan_policy='omit'), &
        5.045378491522287_dp, tol)
    call check_close('hmean omit', hmean(x, weights=weights, nan_policy='omit'), &
        3.6455696202531644_dp, tol)
    call check_close('pmean omit', pmean(x, 2.0_dp, weights=weights, nan_policy='omit'), &
        6.8282501418738315_dp, tol)
    call check_close('sem omit', sem(x, nan_policy='omit'), 2.3333333333333335_dp, tol)
    call check_close('variation omit', variation(x, nan_policy='omit'), &
        0.7071067811865475_dp, tol)
    call check_close('iqr omit', iqr(x, nan_policy='omit'), 4.0_dp, tol)
    call check_close('mad omit', median_abs_deviation(x, nan_policy='omit'), 3.0_dp, tol)
    call check_close('moment omit', moment(x, order=2, nan_policy='omit'), &
        10.888888888888888_dp, tol)
    call check_close('skew omit', skew(x, nan_policy='omit'), &
        0.29479962014482863_dp, tol)
    call check_close('skew unbiased omit', skew(x, bias=.false., nan_policy='omit'), &
        0.7221086457211349_dp, tol)
    call check_close('kurtosis omit', kurtosis(x, nan_policy='omit'), &
        -1.4999999999999996_dp, tol)
    call check_close('kurtosis unbiased omit', &
        kurtosis(x, bias=.false., nan_policy='omit'), -1.4999999999999996_dp, tol)

    call check_true('propagate still NaN', ieee_is_nan(skew(x, nan_policy='propagate')))
    call check_true('invalid policy NaN', ieee_is_nan(skew(x, nan_policy='bad')))
    call check_true('all NaN omit', ieee_is_nan(iqr(all_nan, nan_policy='omit')))
    call check_true('moment zero all NaN omit', &
        ieee_is_nan(moment(all_nan, order=0, nan_policy='omit')))

    values = gmean(matrix, axis=1, nan_policy='omit')
    call check_array('matrix gmean axis 1 omit', values, &
        [1.7320508075688774_dp, 5.47722557505166_dp, &
        6.32455532033676_dp, 10.026595869929169_dp], tol)
    values = hmean(matrix, axis=2, nan_policy='omit')
    call check_array('matrix hmean axis 2 omit', values, &
        [2.181818181818182_dp, 4.655172413793104_dp, 8.873239436619718_dp], tol)
    values = pmean(matrix, 2.0_dp, axis=2, nan_policy='omit')
    call check_array('matrix pmean axis 2 omit', values, &
        [5.196152422706632_dp, 6.191391873668904_dp, 10.519822558706334_dp], tol)
    values = sem(matrix, axis=1, nan_policy='omit')
    call check_array('matrix sem axis 1 omit', values, &
        [1.0_dp, 0.5_dp, 2.9999999999999996_dp, 1.855921454276674_dp], tol)
    values = variation(matrix, axis=2, nan_policy='omit')
    call check_array('matrix variation axis 2 omit', values, &
        [0.6617173282340483_dp, 0.4401949866792872_dp, 0.32659863237109044_dp], tol)
    values = iqr(matrix, axis=1, nan_policy='omit')
    call check_array('matrix iqr axis 1 omit', values, [1.0_dp, 0.5_dp, 3.0_dp, 3.0_dp], tol)
    values = median_abs_deviation(matrix, axis=2, nan_policy='omit')
    call check_array('matrix mad axis 2 omit', values, [3.0_dp, 2.0_dp, 4.0_dp], tol)
    values = moment(matrix, order=2, axis=1, nan_policy='omit')
    call check_array('matrix moment axis 1 omit', values, &
        [1.0_dp, 0.25_dp, 9.0_dp, 6.888888888888888_dp], tol)
    values = skew(matrix, axis=2, nan_policy='omit')
    call check_array('matrix skew axis 2 omit', values, &
        [0.1728005440786506_dp, 0.3818017741606059_dp, 0.0_dp], tol)
    values = kurtosis(matrix, axis=1, nan_policy='omit')
    call check_array('matrix kurtosis axis 1 omit', values, &
        [-2.0_dp, -2.0_dp, -2.0_dp, -1.5000000000000007_dp], tol)

    call check_close('matrix iqr all omit', iqr(matrix, nan_policy='omit'), 5.0_dp, tol)
    values = gmean(matrix, axis=2, weights=weights, nan_policy='omit')
    call check_array('matrix gmean weights4 omit', values, &
        [4.756828460010883_dp, 6.50350275769788_dp, 10.366832259646149_dp], tol)
    values = gmean(matrix, axis=1, weights=weights31, nan_policy='omit')
    call check_array('matrix gmean weights31 omit', values, &
        [2.080083823051904_dp, 5.646216173286171_dp, &
        8.325532074018735_dp, 11.391577742394356_dp], tol)
    values = gmean(matrix, axis=1, weights=weights31_nan, nan_policy='omit')
    call check_array('matrix gmean NaN weights omit', values, &
        [1.0_dp, 6.0_dp, 8.325532074018735_dp, 12.517581449211738_dp], tol)
    values = hmean(matrix, axis=1, weights=weights31_nan, nan_policy='omit')
    call check_array('matrix hmean NaN weights omit', values, &
        [1.0_dp, 6.0_dp, 7.6923076923076925_dp, 12.173913043478262_dp], tol)
    values = pmean(matrix, 2.0_dp, axis=1, weights=weights31_nan, nan_policy='omit')
    call check_array('matrix pmean NaN weights omit', values, &
        [1.0_dp, 6.0_dp, 9.121403400793104_dp, 13.023056476879765_dp], tol)

    call check_array('vector gmean keepdims omit', &
        gmean_keepdims(x, nan_policy='omit'), [3.3019272488946263_dp], tol)
    kept = gmean_keepdims(matrix, axis=1, nan_policy='omit')
    call check_matrix('gmean keepdims omit', kept, reshape([ &
        1.7320508075688774_dp, 5.47722557505166_dp, &
        6.32455532033676_dp, 10.026595869929169_dp], [1, 4]), tol)
    kept = iqr_keepdims(matrix, nan_policy='omit')
    call check_matrix('iqr keepdims all omit', kept, reshape([5.0_dp], [1, 1]), tol)
    kept = sem_keepdims(matrix, axis=2, nan_policy='omit')
    call check_matrix('sem keepdims omit', kept, reshape([ &
        2.0275875100994067_dp, 1.7638342073763937_dp, &
        2.3094010767585034_dp], [3, 1]), tol)

    values = skew(matrix, axis=1, nan_policy='propagate')
    call check_true('matrix propagate NaNs', ieee_is_nan(values(1)) .and. &
        ieee_is_nan(values(2)) .and. ieee_is_nan(values(3)))
    call check_close('matrix propagate finite slice', values(4), &
        0.6309038567106232_dp, tol)

    print *, 'descriptive nan_policy tests passed'

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

        if (size(actual) /= size(expected) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
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

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
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

end program test_descriptive_nan_policy
