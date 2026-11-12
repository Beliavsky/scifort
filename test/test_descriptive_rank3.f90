! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_rank3
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, iqr, iqr_keepdims, &
        kurtosis, median_abs_deviation, moment, pmean, sem, sem_keepdims, skew, variation
    implicit none

    real(dp), parameter :: tol = 3.0e-14_dp
    real(dp) :: nan
    real(dp) :: x(4, 3, 2)
    real(dp) :: xn(4, 3, 2)
    real(dp), allocatable :: values(:, :)
    real(dp), allocatable :: kept(:, :, :)

    x = reshape([ &
        1.0_dp, 2.0_dp, 4.0_dp, 9.0_dp, &
        3.0_dp, 5.0_dp, 8.0_dp, 15.0_dp, &
        2.0_dp, 6.0_dp, 11.0_dp, 20.0_dp, &
        4.0_dp, 7.0_dp, 13.0_dp, 24.0_dp, &
        5.0_dp, 10.0_dp, 18.0_dp, 31.0_dp, &
        6.0_dp, 12.0_dp, 22.0_dp, 40.0_dp], [4, 3, 2])

    values = gmean(x)
    call check_matrix('rank3 gmean default axis', values, reshape([ &
        2.9129506302439405_dp, 6.513555624326305_dp, 7.16804926843144_dp, &
        9.667811105125264_dp, 12.924122054704554_dp, 15.865501049788389_dp], &
        [3, 2]), tol)

    values = hmean(x, axis=2)
    call check_matrix('rank3 hmean axis 2', values, reshape([ &
        1.6363636363636365_dp, 3.4615384615384617_dp, 6.439024390243903_dp, &
        13.170731707317076_dp, 4.864864864864864_dp, 9.197080291970803_dp, &
        16.860262008733624_dp, 30.326086956521745_dp], [4, 2]), tol)

    values = pmean(x, 2.0_dp, axis=3)
    call check_matrix('rank3 pmean axis 3', values, reshape([ &
        2.9154759474226504_dp, 5.1478150704935_dp, 9.617692030835672_dp, &
        18.12456896039186_dp, 4.123105625617661_dp, 7.905694150420948_dp, &
        13.92838827718412_dp, 24.351591323771842_dp, 4.47213595499958_dp, &
        9.486832980505138_dp, 17.392527130926087_dp, 31.622776601683793_dp], &
        [4, 3]), tol)

    values = sem(x, axis=1)
    call check_matrix('rank3 sem axis 1', values, reshape([ &
        1.7795130420052185_dp, 2.6259918760981216_dp, 3.881043674065006_dp, &
        4.415880433163924_dp, 5.671566509057851_dp, 7.438637868140466_dp], &
        [3, 2]), tol)

    values = variation(x, axis=2)
    call check_matrix('rank3 variation axis 2', values, reshape([ &
        0.408248290463863_dp, 0.3922322702763681_dp, 0.37401414204533157_dp, &
        0.30660767189163735_dp, 0.16329931618554522_dp, 0.21256600010237853_dp, &
        0.20840303806013702_dp, 0.2068092916251421_dp], [4, 2]), tol)

    call check_close('rank3 iqr all', iqr(x), 11.0_dp, tol)
    values = iqr(x, axis=3)
    call check_matrix('rank3 iqr axis 3', values, reshape([ &
        1.5_dp, 2.5_dp, 4.5_dp, 7.5_dp, 1.0_dp, 2.5_dp, &
        5.0_dp, 8.0_dp, 2.0_dp, 3.0_dp, 5.5_dp, 10.0_dp], [4, 3]), tol)

    values = median_abs_deviation(x, axis=1)
    call check_matrix('rank3 mad axis 1', values, reshape([ &
        1.5_dp, 2.5_dp, 4.5_dp, 4.5_dp, 6.5_dp, 8.0_dp], [3, 2]), tol)

    values = moment(x, order=2, axis=2)
    call check_matrix('rank3 moment axis 2', values, reshape([ &
        0.6666666666666666_dp, 2.888888888888889_dp, 8.222222222222221_dp, &
        20.222222222222225_dp, 0.6666666666666666_dp, 4.222222222222222_dp, &
        13.555555555555557_dp, 42.888888888888886_dp], [4, 2]), tol)

    values = skew(x, axis=1, bias=.false.)
    call check_matrix('rank3 skew axis 1 unbiased', values, reshape([ &
        1.3309377322476705_dp, 1.164859805755412_dp, 0.7991818725885826_dp, &
        1.0567918985349558_dp, 0.838653073421854_dp, 0.9620883060695536_dp], &
        [3, 2]), tol)

    values = kurtosis(x, axis=2, bias=.false.)
    call check_matrix('rank3 kurtosis axis 2 unbiased', values, reshape([ &
        -1.5_dp, -1.5000000000000004_dp, -1.4999999999999998_dp, &
        -1.5000000000000004_dp, -1.5_dp, -1.5000000000000002_dp, &
        -1.5_dp, -1.5_dp], [4, 2]), tol)

    kept = gmean_keepdims(x, axis=2)
    call check_tensor('rank3 gmean keepdims axis 2', kept, reshape([ &
        1.8171205928321394_dp, 3.9148676411688634_dp, 7.060696670652125_dp, &
        13.924766500838336_dp, 4.93242414866094_dp, 9.43538796063307_dp, &
        17.26683925893647_dp, 30.989243579983132_dp], [4, 1, 2]), tol)

    kept = iqr_keepdims(x)
    call check_tensor('rank3 iqr keepdims all', kept, reshape([11.0_dp], [1, 1, 1]), tol)

    call check_true('rank3 invalid axis', all(shape(sem(x, axis=4)) == [0, 0]))
    call check_true('rank3 invalid keepdims axis', &
        all(shape(gmean_keepdims(x, axis=4)) == [0, 0, 0]))

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    xn = x
    xn(1, 1, 1) = nan
    xn(3, 2, 1) = nan
    xn(2, 3, 2) = nan
    xn(4, 1, 2) = nan

    values = gmean(xn, axis=1, nan_policy='omit')
    call check_matrix('rank3 gmean omit', values, reshape([ &
        4.160167646103808_dp, 6.0822019955734_dp, 7.16804926843144_dp, &
        7.140036981921557_dp, 12.924122054704554_dp, 17.413175382347223_dp], &
        [3, 2]), tol)

    values = hmean(xn, axis=2, nan_policy='omit')
    call check_matrix('rank3 hmean omit', values, reshape([ &
        2.4_dp, 3.4615384615384617_dp, 5.866666666666666_dp, &
        13.170731707317076_dp, 4.864864864864864_dp, 8.23529411764706_dp, &
        16.860262008733624_dp, 34.929577464788736_dp], [4, 2]), tol)

    values = pmean(xn, 2.0_dp, axis=3, nan_policy='omit')
    call check_matrix('rank3 pmean omit', values, reshape([ &
        4.0_dp, 5.1478150704935_dp, 9.617692030835672_dp, 9.0_dp, &
        4.123105625617661_dp, 7.905694150420948_dp, 18.0_dp, 24.351591323771842_dp, &
        4.47213595499958_dp, 6.0_dp, 17.392527130926087_dp, 31.622776601683793_dp], &
        [4, 3]), tol)

    values = sem(xn, axis=1, nan_policy='omit')
    call check_matrix('rank3 sem omit', values, reshape([ &
        2.0816659994661326_dp, 3.711842908553348_dp, 3.881043674065006_dp, &
        2.6457513110645907_dp, 5.671566509057851_dp, 9.820613241770825_dp], &
        [3, 2]), tol)

    values = variation(xn, axis=2, nan_policy='omit')
    call check_matrix('rank3 variation omit', values, reshape([ &
        0.2_dp, 0.3922322702763681_dp, 0.4666666666666667_dp, &
        0.30660767189163735_dp, 0.16329931618554522_dp, 0.17647058823529413_dp, &
        0.20840303806013702_dp, 0.1267605633802817_dp], [4, 2]), tol)

    call check_close('rank3 iqr all omit', iqr(xn, nan_policy='omit'), 11.0_dp, tol)
    values = iqr(xn, axis=3, nan_policy='omit')
    call check_matrix('rank3 iqr omit', values, reshape([ &
        0.0_dp, 2.5_dp, 4.5_dp, 0.0_dp, 1.0_dp, 2.5_dp, &
        0.0_dp, 8.0_dp, 2.0_dp, 0.0_dp, 5.5_dp, 10.0_dp], [4, 3]), tol)

    values = median_abs_deviation(xn, axis=1, nan_policy='omit')
    call check_matrix('rank3 mad omit', values, reshape([ &
        2.0_dp, 2.0_dp, 4.5_dp, 3.0_dp, 6.5_dp, 16.0_dp], [3, 2]), tol)

    values = moment(xn, order=2, axis=2, nan_policy='omit')
    call check_matrix('rank3 moment omit', values, reshape([ &
        0.25_dp, 2.888888888888889_dp, 12.25_dp, 20.222222222222225_dp, &
        0.6666666666666666_dp, 2.25_dp, 13.555555555555557_dp, 20.25_dp], &
        [4, 2]), tol)

    values = skew(xn, axis=1, bias=.false., nan_policy='omit')
    call check_matrix('rank3 skew omit', values, reshape([ &
        1.1520696383139375_dp, 1.5453925256950205_dp, 0.7991818725885826_dp, &
        0.9352195295828245_dp, 0.838653073421854_dp, 0.1760979853981283_dp], &
        [3, 2]), tol)

    values = kurtosis(xn, axis=2, bias=.false., nan_policy='omit')
    call check_matrix('rank3 kurtosis omit', values, reshape([ &
        -2.0_dp, -1.5000000000000004_dp, -2.0_dp, -1.5000000000000004_dp, &
        -1.5_dp, -2.0_dp, -1.5_dp, -2.0_dp], [4, 2]), tol)

    kept = sem_keepdims(xn, axis=3, nan_policy='omit')
    call check_tensor_nan('rank3 sem keepdims omit', kept, reshape([ &
        nan, 2.5_dp, 4.5_dp, nan, 1.0_dp, 2.5_dp, nan, 8.0_dp, &
        2.0_dp, nan, 5.5_dp, 10.0_dp], [4, 3, 1]), tol)

    values = skew(xn, axis=1, nan_policy='propagate')
    call check_true('rank3 propagate NaN slice', ieee_is_nan(values(1, 1)))
    call check_close('rank3 propagate finite slice', values(3, 1), &
        0.4614078692704874_dp, tol)

    print *, 'descriptive rank-3 tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance

        if (ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_matrix(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :)
        real(dp), intent(in) :: expected(:, :)
        real(dp), intent(in) :: tolerance

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

    subroutine check_tensor(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :, :)
        real(dp), intent(in) :: expected(:, :, :)
        real(dp), intent(in) :: tolerance

        if (any(shape(actual) /= shape(expected)) .or. any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual shape  = ', shape(actual)
            print *, ' expected shape= ', shape(expected)
            error stop 1
        end if
    end subroutine check_tensor

    subroutine check_tensor_nan(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :, :)
        real(dp), intent(in) :: expected(:, :, :)
        real(dp), intent(in) :: tolerance
        integer :: i
        integer :: j
        integer :: k

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', name, ' shape'
            error stop 1
        end if
        do k = 1, size(actual, 3)
            do j = 1, size(actual, 2)
                do i = 1, size(actual, 1)
                    if (ieee_is_nan(expected(i, j, k))) then
                        if (.not. ieee_is_nan(actual(i, j, k))) then
                            print *, 'FAIL ', name
                            error stop 1
                        end if
                    else
                        if (ieee_is_nan(actual(i, j, k))) then
                            print *, 'FAIL ', name
                            error stop 1
                        end if
                        if (abs(actual(i, j, k) - expected(i, j, k)) > tolerance) then
                            print *, 'FAIL ', name
                            error stop 1
                        end if
                    end if
                end do
            end do
        end do
    end subroutine check_tensor_nan

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_rank3
