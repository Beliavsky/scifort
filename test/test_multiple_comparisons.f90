program test_multiple_comparisons
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use multiple_comparisons_reference
    use scifort_kinds, only : dp
    use scifort_stats, only : alexandergovern, alexandergovern_result, dunnett, &
        dunnett_result, make_sample_group, matrix_confidence_interval, &
        poisson_means_test, poisson_means_test_result, sample_group, tukey_hsd, &
        tukey_hsd_result, vector_confidence_interval
    implicit none

    type(sample_group) :: ag(4)
    type(sample_group) :: tg(3)
    type(sample_group) :: ug(3)
    type(sample_group) :: dg(2)
    type(alexandergovern_result) :: ar
    type(tukey_hsd_result) :: tr
    type(dunnett_result) :: dr
    type(poisson_means_test_result) :: pr
    type(matrix_confidence_interval) :: mci
    type(vector_confidence_interval) :: vci
    integer :: failures
    integer :: idx
    real(dp), parameter :: control(5) = [2.9_dp, 3.0_dp, 2.5_dp, 2.6_dp, 3.2_dp]

    failures = 0

    ag(1) = make_sample_group([13.75_dp, 13.75_dp, 13.5_dp, 13.5_dp, 13.0_dp, &
        13.0_dp, 13.0_dp, 12.75_dp, 12.5_dp])
    ag(2) = make_sample_group([14.25_dp, 13.0_dp, 12.75_dp, 12.5_dp, 12.5_dp, &
        12.4_dp, 12.3_dp, 11.9_dp, 11.9_dp])
    ag(3) = make_sample_group([14.0_dp, 14.0_dp, 13.51_dp, 13.5_dp, 13.5_dp, &
        13.25_dp, 13.0_dp, 12.5_dp, 12.5_dp])
    ag(4) = make_sample_group([15.0_dp, 14.0_dp, 13.75_dp, 13.59_dp, 13.25_dp, &
        12.97_dp, 12.5_dp, 12.25_dp, 11.89_dp])
    ar = alexandergovern(ag)
    call check_close('Alexander statistic', ar%statistic, alexander_ref(1), 2.0e-13_dp, failures)
    call check_close('Alexander pvalue', ar%pvalue, alexander_ref(2), 2.0e-14_dp, failures)

    tg(1) = make_sample_group([24.5_dp, 23.5_dp, 26.4_dp, 27.1_dp, 29.9_dp])
    tg(2) = make_sample_group([28.4_dp, 34.2_dp, 29.5_dp, 32.2_dp, 30.1_dp])
    tg(3) = make_sample_group([26.1_dp, 28.3_dp, 24.3_dp, 26.2_dp, 27.8_dp])
    tr = tukey_hsd(tg)
    mci = tr%confidence_interval()
    call check_matrix('Tukey statistic', tr%statistic, reshape(tukey_equal_ref(1:9), [3, 3]), 4.0e-15_dp, failures)
    call check_matrix('Tukey pvalue', tr%pvalue, reshape(tukey_equal_ref(10:18), [3, 3]), 4.0e-11_dp, failures)
    call check_matrix('Tukey CI low', mci%low, reshape(tukey_equal_ref(19:27), [3, 3]), 5.0e-10_dp, failures)
    call check_matrix('Tukey CI high', mci%high, reshape(tukey_equal_ref(28:36), [3, 3]), 5.0e-10_dp, failures)

    tr = tukey_hsd(tg, equal_var=.false.)
    mci = tr%confidence_interval()
    call check_matrix('Games-Howell statistic', tr%statistic, reshape(tukey_games_ref(1:9), [3, 3]), 4.0e-15_dp, failures)
    call check_matrix('Games-Howell pvalue', tr%pvalue, reshape(tukey_games_ref(10:18), [3, 3]), 5.0e-10_dp, failures)
    call check_matrix('Games-Howell CI low', mci%low, reshape(tukey_games_ref(19:27), [3, 3]), 8.0e-9_dp, failures)
    call check_matrix('Games-Howell CI high', mci%high, reshape(tukey_games_ref(28:36), [3, 3]), 8.0e-9_dp, failures)

    ug(1) = make_sample_group([1.2_dp, 2.1_dp, 1.7_dp, 2.4_dp])
    ug(2) = make_sample_group([3.2_dp, 2.8_dp, 3.6_dp, 3.1_dp, 2.9_dp, 3.7_dp])
    ug(3) = make_sample_group([2.0_dp, 1.8_dp, 2.6_dp, 2.5_dp, 2.3_dp])
    tr = tukey_hsd(ug)
    mci = tr%confidence_interval(0.90_dp)
    call check_matrix('Tukey-Kramer statistic', tr%statistic, reshape(tukey_unequal_ref(1:9), [3, 3]), 4.0e-15_dp, failures)
    call check_matrix('Tukey-Kramer pvalue', tr%pvalue, reshape(tukey_unequal_ref(10:18), [3, 3]), 5.0e-10_dp, failures)
    call check_matrix('Tukey-Kramer CI low', mci%low, reshape(tukey_unequal_ref(19:27), [3, 3]), 8.0e-9_dp, failures)
    call check_matrix('Tukey-Kramer CI high', mci%high, reshape(tukey_unequal_ref(28:36), [3, 3]), 8.0e-9_dp, failures)

    dg(1) = make_sample_group([3.8_dp, 2.7_dp, 4.0_dp, 2.4_dp])
    dg(2) = make_sample_group([2.8_dp, 3.4_dp, 3.7_dp, 2.2_dp, 2.0_dp])
    idx = 1
    call check_dunnett('two-sided', dunnett_ref(idx:idx + 7), failures)
    idx = idx + 8
    call check_dunnett('greater', dunnett_ref(idx:idx + 7), failures)
    idx = idx + 8
    call check_dunnett('less', dunnett_ref(idx:idx + 7), failures)

    idx = 1
    call check_poisson_case(0, 100.0_dp, 3, 100.0_dp, 0.0_dp, poisson_means_ref(idx:idx + 5), failures)
    idx = idx + 6
    call check_poisson_case(5, 30.0_dp, 1, 25.0_dp, 0.0_dp, poisson_means_ref(idx:idx + 5), failures)
    idx = idx + 6
    call check_poisson_case(5, 30.0_dp, 1, 25.0_dp, 0.1_dp, poisson_means_ref(idx:idx + 5), failures)
    idx = idx + 6
    call check_poisson_case(12, 40.0_dp, 8, 25.0_dp, 0.05_dp, poisson_means_ref(idx:idx + 5), failures)

    ar = alexandergovern(tg(1:1))
    call check_true('Alexander invalid group count', ieee_is_nan(ar%statistic), failures)
    tr = tukey_hsd(tg(1:1))
    call check_true('Tukey invalid group count', allocated(tr%statistic) .and. &
        size(tr%statistic) == 1 .and. ieee_is_nan(tr%statistic(1, 1)), failures)
    pr = poisson_means_test(-1, 1.0_dp, 0, 1.0_dp)
    call check_true('Poisson invalid count', ieee_is_nan(pr%statistic), failures)
    pr = poisson_means_test(0, 100.0_dp, 0, 100.0_dp, diff=1.0_dp)
    call check_close('Poisson nonpositive null estimate statistic', pr%statistic, 0.0_dp, 0.0_dp, failures)
    call check_close('Poisson nonpositive null estimate pvalue', pr%pvalue, 1.0_dp, 0.0_dp, failures)

    if (failures /= 0) then
        print *, 'multiple-comparison failures:', failures
        error stop 1
    end if
    print *, 'multiple-comparison tests passed'

contains

    subroutine check_dunnett(alternative, expected_value, count)
        character(len=*), intent(in) :: alternative !! Dunnett alternative
        real(dp), intent(in) :: expected_value(8) !! statistics, p-values, lower and upper CIs
        integer, intent(inout) :: count !! accumulated failure count

        dr = dunnett(dg, control, alternative, maxpts=8000)
        vci = dr%confidence_interval(0.95_dp, maxpts=8000)
        call check_vector('Dunnett statistic', dr%statistic, expected_value(1:2), 2.0e-14_dp, count)
        call check_vector('Dunnett pvalue', dr%pvalue, expected_value(3:4), 8.0e-4_dp, count)
        if (trim(alternative) == 'two-sided') then
            call check_vector('Dunnett CI low', vci%low, expected_value(5:6), 1.5e-2_dp, count)
            call check_vector('Dunnett CI high', vci%high, expected_value(7:8), 1.5e-2_dp, count)
        else if (trim(alternative) == 'greater') then
            call check_vector('Dunnett CI low', vci%low, expected_value(5:6), 1.5e-2_dp, count)
            call check_true('Dunnett greater high infinity', all(.not. ieee_is_finite(vci%high)), count)
        else
            call check_true('Dunnett less low infinity', all(.not. ieee_is_finite(vci%low)), count)
            call check_vector('Dunnett CI high', vci%high, expected_value(7:8), 1.5e-2_dp, count)
        end if
    end subroutine check_dunnett

    subroutine check_poisson_case(k1, n1, k2, n2, diff, expected_value, count)
        integer, intent(in) :: k1 !! first count
        real(dp), intent(in) :: n1 !! first exposure
        integer, intent(in) :: k2 !! second count
        real(dp), intent(in) :: n2 !! second exposure
        real(dp), intent(in) :: diff !! null difference
        real(dp), intent(in) :: expected_value(6) !! statistic/p-value for three alternatives
        integer, intent(inout) :: count !! accumulated failure count
        character(len=9), parameter :: alternatives(3) = ['two-sided', 'less     ', 'greater  ']
        integer :: i

        do i = 1, 3
            pr = poisson_means_test(k1, n1, k2, n2, diff, trim(alternatives(i)))
            call check_close('Poisson means statistic', pr%statistic, expected_value(2 * i - 1), 2.0e-14_dp, count)
            call check_close('Poisson means pvalue', pr%pvalue, expected_value(2 * i), 2.0e-10_dp, count)
        end do
    end subroutine check_poisson_case

    subroutine check_matrix(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:, :) !! computed matrix
        real(dp), intent(in) :: expected_value(:, :) !! SciPy reference matrix
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count
        integer :: i
        integer :: j

        if (any(shape(actual) /= shape(expected_value))) then
            print *, 'FAIL ', trim(name), ' shape'
            count = count + 1
            return
        end if
        do j = 1, size(actual, 2)
            do i = 1, size(actual, 1)
                call check_close(name, actual(i, j), expected_value(i, j), tolerance, count)
            end do
        end do
    end subroutine check_matrix

    subroutine check_vector(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! computed vector
        real(dp), intent(in) :: expected_value(:) !! SciPy reference vector
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count
        integer :: i

        if (size(actual) /= size(expected_value)) then
            print *, 'FAIL ', trim(name), ' size'
            count = count + 1
            return
        end if
        do i = 1, size(actual)
            call check_close(name, actual(i), expected_value(i), tolerance, count)
        end do
    end subroutine check_vector

    subroutine check_close(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected_value !! SciPy reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count

        if (ieee_is_nan(actual) .or. abs(actual - expected_value) > tolerance) then
            print *, 'FAIL ', trim(name), actual, expected_value, abs(actual - expected_value)
            count = count + 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition, count)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition expected to be true
        integer, intent(inout) :: count !! accumulated failure count

        if (.not. condition) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_true

end program test_multiple_comparisons
