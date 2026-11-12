! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_frequency
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : cumfreq, cumfreq_result, percentileofscore, relfreq, &
        relfreq_result, scoreatpercentile, scoreatpercentile_all
    implicit none

    real(dp), parameter :: tol = 5.0e-13_dp
    real(dp) :: a(5)
    real(dp) :: matrix(3, 4)
    real(dp) :: qnan
    integer :: i
    integer :: shape15(15)
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: m(:, :)
    type(cumfreq_result) :: cf
    type(relfreq_result) :: rf
    real(dp) :: weights(4)
    real(dp) :: x15(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: s15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: p15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: sap15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

    a = [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp]
    call check_scalar('percentile rank tie', percentileofscore(a, 3.0_dp), 70.0_dp)
    call check_scalar('percentile strict tie', percentileofscore(a, 3.0_dp, kind='strict'), 40.0_dp)
    call check_scalar('percentile weak tie', percentileofscore(a, 3.0_dp, kind='weak'), 80.0_dp)
    call check_scalar('percentile mean tie', percentileofscore(a, 3.0_dp, kind='mean'), 60.0_dp)
    call check_scalar('percentile gap', percentileofscore(a, 2.5_dp), 40.0_dp)
    v = percentileofscore(a, [2.0_dp, 3.0_dp])
    call check_array('percentile score vector', v, [40.0_dp, 70.0_dp])
    v = percentileofscore([ieee_value(0.0_dp, ieee_negative_inf), 0.0_dp, 1.0_dp, &
        ieee_value(0.0_dp, ieee_positive_inf)], &
        [1.0_dp, 2.0_dp, ieee_value(0.0_dp, ieee_positive_inf)])
    call check_array('percentile finite-extreme ordering', v, [75.0_dp, 75.0_dp, 100.0_dp])
    if (.not. ieee_is_nan(percentileofscore([real(dp) ::], 2.0_dp))) then
        print *, 'FAIL percentile empty'
        error stop 1
    end if

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. ieee_is_nan(percentileofscore([1.0_dp, qnan, 3.0_dp], 2.0_dp))) then
        print *, 'FAIL percentile propagate NaN'
        error stop 1
    end if
    call check_scalar('percentile omit NaN', &
        percentileofscore([1.0_dp, qnan, 3.0_dp], 2.0_dp, nan_policy='omit'), 50.0_dp)
    if (.not. ieee_is_nan(percentileofscore(a, qnan))) then
        print *, 'FAIL percentile NaN score'
        error stop 1
    end if

    matrix = reshape([(real(i, dp), i=1,12)], shape(matrix))
    call check_scalar('scoreatpercentile all median', scoreatpercentile(matrix, 50.0_dp), 6.5_dp)
    call check_scalar('scoreatpercentile_all alias', scoreatpercentile_all(matrix, 50.0_dp), 6.5_dp)
    v = scoreatpercentile(matrix, 50.0_dp, axis=1)
    call check_array('scoreatpercentile axis1', v, [2.0_dp, 5.0_dp, 8.0_dp, 11.0_dp])
    v = scoreatpercentile(matrix, 50.0_dp, axis=2)
    call check_array('scoreatpercentile axis2', v, [5.5_dp, 6.5_dp, 7.5_dp])
    m = scoreatpercentile(matrix, [25.0_dp, 50.0_dp, 75.0_dp], axis=1)
    if (any(shape(m) /= [3, 4])) then
        print *, 'FAIL scoreatpercentile vector shape', shape(m)
        error stop 1
    end if
    call check_array('scoreatpercentile vector axis1', reshape(m, [size(m)]), &
        [1.5_dp, 2.0_dp, 2.5_dp, 4.5_dp, 5.0_dp, 5.5_dp, &
         7.5_dp, 8.0_dp, 8.5_dp, 10.5_dp, 11.0_dp, 11.5_dp])
    v = scoreatpercentile(matrix, [25.0_dp, 50.0_dp, 75.0_dp])
    call check_array('scoreatpercentile vector all', v, [3.75_dp, 6.5_dp, 9.25_dp])
    call check_scalar('score fraction', scoreatpercentile([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], &
        25.0_dp), 1.75_dp)
    call check_scalar('score lower', scoreatpercentile([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], &
        25.0_dp, interpolation_method='lower'), 1.0_dp)
    call check_scalar('score higher', scoreatpercentile([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], &
        25.0_dp, interpolation_method='higher'), 2.0_dp)
    v = scoreatpercentile([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, qnan], &
        [0.0_dp, 50.0_dp, 100.0_dp], limit=[2.0_dp, 8.0_dp])
    call check_array('score limit', v, [2.0_dp, 4.0_dp, 8.0_dp])
    call check_scalar('score NaN sorted last', &
        scoreatpercentile([1.0_dp, 2.0_dp, qnan, 8.0_dp], 50.0_dp), 5.0_dp)
    if (.not. ieee_is_nan(scoreatpercentile([real(dp) ::], 50.0_dp))) then
        print *, 'FAIL score empty'
        error stop 1
    end if

    cf = cumfreq([1.0_dp, 4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 1.0_dp], &
        numbins=4, defaultreallimits=[1.5_dp, 5.0_dp])
    call check_array('cumfreq counts', cf%cumcount, [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp])
    call check_scalar('cumfreq lower', cf%lowerlimit, 1.5_dp)
    call check_scalar('cumfreq binsize', cf%binsize, 0.875_dp)
    if (cf%extrapoints /= 3) then
        print *, 'FAIL cumfreq extrapoints', cf%extrapoints
        error stop 1
    end if

    rf = relfreq([2.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], numbins=4)
    call check_array('relfreq default', rf%frequency, &
        [1.0_dp / 6.0_dp, 0.5_dp, 1.0_dp / 6.0_dp, 1.0_dp / 6.0_dp])
    call check_scalar('relfreq lower', rf%lowerlimit, 0.5_dp)
    call check_scalar('relfreq binsize', rf%binsize, 1.0_dp)

    weights = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    cf = cumfreq([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], numbins=2, &
        defaultreallimits=[1.0_dp, 4.0_dp], weights=weights)
    call check_array('cumfreq weighted', cf%cumcount, [3.0_dp, 10.0_dp])
    rf = relfreq([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], numbins=2, &
        defaultreallimits=[1.0_dp, 4.0_dp], weights=weights)
    call check_array('relfreq weighted', rf%frequency, [0.75_dp, 1.75_dp])

    rf = relfreq(reshape([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [2, 2]), &
        numbins=2, defaultreallimits=[1.0_dp, 4.0_dp])
    call check_array('relfreq multidim SciPy normalization', rf%frequency, [1.0_dp, 1.0_dp])

    rf = relfreq([1.0_dp, 1.0_dp, 1.0_dp], numbins=4)
    call check_array('relfreq constant', rf%frequency, [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp])
    call check_scalar('relfreq constant lower', rf%lowerlimit, 1.0_dp)
    call check_scalar('relfreq constant binsize', rf%binsize, 0.25_dp)

    cf = cumfreq([real(dp) ::], numbins=4)
    call check_array('cumfreq empty', cf%cumcount, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    call check_scalar('cumfreq empty lower', cf%lowerlimit, 0.0_dp)
    call check_scalar('cumfreq empty binsize', cf%binsize, 0.25_dp)
    rf = relfreq([real(dp) ::], numbins=4)
    if (.not. all(ieee_is_nan(rf%frequency))) then
        print *, 'FAIL relfreq empty', rf%frequency
        error stop 1
    end if

    rf = relfreq(reshape([real(dp) ::], [0, 2]), numbins=4)
    if (.not. all(ieee_is_nan(rf%frequency))) then
        error stop 'relfreq empty leading dimension must return NaN'
    end if
    rf = relfreq(reshape([real(dp) ::], [2, 0]), numbins=4)
    call check_array('relfreq empty trailing dimension', rf%frequency, &
        [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])

    cf = cumfreq([1.0_dp, qnan, 2.0_dp], numbins=4, defaultreallimits=[0.0_dp, 3.0_dp])
    call check_array('cumfreq explicit-limit NaN', cf%cumcount, [0.0_dp, 1.0_dp, 2.0_dp, 2.0_dp])
    rf = relfreq([1.0_dp, qnan, 2.0_dp], numbins=4, defaultreallimits=[0.0_dp, 3.0_dp])
    call check_array('relfreq explicit-limit NaN', rf%frequency, &
        [0.0_dp, 1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 0.0_dp])

    x15(:, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = [1.0_dp, 2.0_dp, 8.0_dp]
    s15(:, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = [2.5_dp, 3.0_dp]
    p15 = percentileofscore(a, s15)
    if (any(shape(p15) /= shape(s15))) then
        print *, 'FAIL percentile rank15 shape'
        error stop 1
    end if
    call check_array('percentile rank15', reshape(p15, [size(p15)]), [40.0_dp, 70.0_dp])
    sap15 = scoreatpercentile(x15, [0.0_dp, 50.0_dp, 100.0_dp], axis=1)
    shape15 = shape(sap15)
    if (size(sap15, 1) /= 3 .or. any(shape15(2:) /= 1)) then
        print *, 'FAIL score rank15 shape', shape(sap15)
        error stop 1
    end if
    call check_array('score rank15', reshape(sap15, [size(sap15)]), [1.0_dp, 2.0_dp, 8.0_dp])

    print *, 'descriptive frequency tests passed'

contains

    subroutine check_scalar(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected

        if (ieee_is_nan(expected)) then
            if (.not. ieee_is_nan(actual)) then
                print *, 'FAIL ', trim(label), actual, expected
                error stop 1
            end if
        else if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_scalar

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)

        integer :: i

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), ' size ', size(actual), size(expected)
            error stop 1
        end if
        do i = 1, size(actual)
            call check_scalar(label, actual(i), expected(i))
        end do
    end subroutine check_array

end program test_descriptive_frequency
