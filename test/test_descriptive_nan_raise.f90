! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-compatible nan_policy='raise' regression. Expected error termination is
! checked by launching this test executable in a child mode so `fpm test` itself
! still succeeds when each requested NaN case terminates with nonzero status.

program test_descriptive_nan_raise
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : gmean, gstd, gstd_all, gstd_keepdims, gzscore, &
        hmean, iqr, kurtosis, median_abs_deviation, moment, pmean, sem, skew, &
        variation, zmap, zmap_all, zscore, zscore_all, tmean, tvar, tstd, tsem, &
        tmin, tmax, trim_mean
    implicit none

    character(len=1024) :: executable
    character(len=64) :: mode
    character(len=64) :: raise_case
    integer :: argc
    integer :: i
    character(len=24), parameter :: cases(*) = [character(len=24) :: &
        'gmean', 'gmean_weight', 'gmean_weight_rank15', 'hmean', 'pmean', 'sem', &
        'variation', 'iqr', 'iqr_named', 'mad', 'mad_named', 'moment', 'skew', &
        'kurtosis', 'sem_rank5', 'sem_rank15', 'zscore', 'zscore_rank15', 'gzscore', &
        'gstd', 'zmap_compare', 'zscore_all', 'zmap_all', 'gstd_tuple', 'tmean', &
        'tvar', 'tstd', 'tsem', 'tmin', 'tmax', 'trim_mean']
    real(dp) :: nan
    real(dp), allocatable :: values(:)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, mode)
        if (trim(mode) == '--raise-case') then
            if (argc < 2) error stop 2
            call get_command_argument(2, raise_case)
            call run_raise_case(trim(raise_case), nan)
            stop 0
        end if
    end if

    call check_close(gmean([1.0_dp, 4.0_dp], nan_policy='raise'), 2.0_dp, &
        'gmean raise finite')
    call check_close(sem([1.0_dp, 3.0_dp], nan_policy='raise'), 1.0_dp, &
        'sem raise finite')
    values = zscore([1.0_dp, 2.0_dp, 3.0_dp], nan_policy='raise')
    call check_close(sum(values), 0.0_dp, 'zscore raise finite')
    values = zmap([0.0_dp, nan, 2.0_dp], [0.0_dp, 2.0_dp, 4.0_dp], &
        nan_policy='raise')
    if (.not. ieee_is_nan(values(2))) then
        print *, 'FAIL zmap raise score NaN should propagate without raising'
        error stop 1
    end if

    call get_command_argument(0, executable)
    do i = 1, size(cases)
        call expect_error_stop(trim(executable), trim(cases(i)))
    end do

    print *, "descriptive nan_policy='raise' tests passed"

contains

    subroutine run_raise_case(name, nan_value)
        character(len=*), intent(in) :: name !! child test case
        real(dp), intent(in) :: nan_value !! quiet NaN injected into input

        real(dp) :: matrix(2, 2)
        real(dp) :: x5(2, 1, 1, 1, 1)
        real(dp) :: x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        real(dp) :: weights15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        real(dp) :: vector2(2)
        real(dp) :: vector3(3)
        real(dp) :: matrix2(2, 2)
        real(dp) :: matrix1(1, 1)
        real(dp) :: rank4(1, 1, 1, 1)
        real(dp) :: scalar

        select case (name)
        case ('gmean')
            scalar = gmean([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('gmean_weight')
            scalar = gmean([1.0_dp, 2.0_dp, 4.0_dp], &
                weights=[1.0_dp, nan_value, 1.0_dp], nan_policy='raise')
        case ('gmean_weight_rank15')
            x15 = 2.0_dp
            weights15 = 1.0_dp
            weights15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = nan_value
            scalar = sum(gmean(x15, axis=1, weights=weights15, nan_policy='raise'))
        case ('hmean')
            scalar = hmean([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('pmean')
            scalar = pmean([1.0_dp, nan_value, 4.0_dp], 2.0_dp, nan_policy='raise')
        case ('sem')
            scalar = sem([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('variation')
            scalar = variation([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('iqr')
            scalar = iqr([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('iqr_named')
            scalar = iqr([1.0_dp, nan_value, 4.0_dp], scale='normal', &
                nan_policy='raise')
        case ('mad')
            scalar = median_abs_deviation([1.0_dp, nan_value, 4.0_dp], &
                nan_policy='raise')
        case ('mad_named')
            scalar = median_abs_deviation([1.0_dp, nan_value, 4.0_dp], &
                scale='normal', nan_policy='raise')
        case ('moment')
            scalar = moment([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('skew')
            scalar = skew([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('kurtosis')
            scalar = kurtosis([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('sem_rank5')
            x5(:, 1, 1, 1, 1) = [1.0_dp, nan_value]
            rank4 = sem(x5, nan_policy='raise')
        case ('sem_rank15')
            x15 = 2.0_dp
            x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = nan_value
            scalar = sum(sem(x15, axis=1, nan_policy='raise'))
        case ('zscore')
            matrix = reshape([1.0_dp, nan_value, 3.0_dp, 4.0_dp], shape(matrix))
            matrix2 = zscore(matrix, axis=1, nan_policy='raise')
        case ('zscore_rank15')
            x15 = 2.0_dp
            x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = nan_value
            scalar = sum(zscore(x15, axis=1, nan_policy='raise'))
        case ('gzscore')
            vector3 = gzscore([1.0_dp, -2.0_dp, 4.0_dp], nan_policy='raise')
        case ('gstd')
            scalar = gstd([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('zmap_compare')
            vector2 = zmap([0.0_dp, 1.0_dp], [0.0_dp, nan_value, 2.0_dp], &
                nan_policy='raise')
        case ('zscore_all')
            matrix = reshape([1.0_dp, nan_value, 3.0_dp, 4.0_dp], shape(matrix))
            matrix2 = zscore_all(matrix, nan_policy='raise')
        case ('zmap_all')
            vector2 = zmap_all([0.0_dp, 1.0_dp], [0.0_dp, nan_value, 2.0_dp], &
                nan_policy='raise')
        case ('gstd_tuple')
            matrix = reshape([1.0_dp, nan_value, 3.0_dp, 4.0_dp], shape(matrix))
            matrix1 = gstd_keepdims(matrix, [1, 2], nan_policy='raise')
        case ('tmean')
            scalar = tmean([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('tvar')
            scalar = tvar([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('tstd')
            scalar = tstd([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('tsem')
            scalar = tsem([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('tmin')
            scalar = tmin([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('tmax')
            scalar = tmax([1.0_dp, nan_value, 4.0_dp], nan_policy='raise')
        case ('trim_mean')
            scalar = trim_mean([1.0_dp, nan_value, 4.0_dp], 0.1_dp, &
                nan_policy='raise')
        case default
            error stop 4
        end select

    end subroutine run_raise_case

    subroutine expect_error_stop(program_path, name)
        character(len=*), intent(in) :: program_path !! current test executable path
        character(len=*), intent(in) :: name !! child case expected to terminate

        character(len=2200) :: command
        integer :: cmdstat
        integer :: envstat
        integer :: exitstat
        character(len=64) :: os_name

        call get_environment_variable('OS', os_name, status=envstat)
        command = '"' // trim(program_path) // '" --raise-case ' // trim(name)
        if (envstat == 0 .and. index(os_name, 'Windows_NT') > 0) then
            command = trim(command) // ' >NUL 2>&1'
        else
            command = trim(command) // ' >/dev/null 2>&1'
        end if
        call execute_command_line(trim(command), wait=.true., exitstat=exitstat, &
            cmdstat=cmdstat)
        if (cmdstat /= 0 .or. exitstat == 0) then
            print *, 'FAIL expected nan_policy raise: ', trim(name)
            print *, ' cmdstat = ', cmdstat
            print *, ' exitstat= ', exitstat
            error stop 1
        end if
    end subroutine expect_error_stop

    subroutine check_close(actual, expected, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! expected value
        character(len=*), intent(in) :: label !! check label

        if (ieee_is_nan(actual) .or. abs(actual - expected) > 5.0e-13_dp) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

end program test_descriptive_nan_raise
