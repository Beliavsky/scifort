! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_descriptive_circular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_stats, only : circmean, circmean_all, circmean_keepdims, &
        circstd, circstd_all, circstd_keepdims, circvar, circvar_all, &
        circvar_keepdims
    implicit none

    real(dp), parameter :: tol = 5.0e-13_dp
    real(dp) :: angles(3)
    real(dp) :: degrees(3)
    real(dp) :: nan_values(3)
    real(dp) :: matrix(2, 3)
    real(dp), allocatable :: a1(:)
    real(dp), allocatable :: keep_mean(:, :)
    real(dp), allocatable :: keep_std(:, :)
    real(dp), allocatable :: keep_var(:, :)
    real(dp), allocatable :: keep15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: x15(3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: qnan

    angles = [20.0_dp, 30.0_dp, 330.0_dp] * scifort_pi / 180.0_dp
    call check_scalar('circmean basic', circmean(angles), &
        0.12732136153446286_dp)
    call check_scalar('circvar basic', circvar(angles), &
        0.10215129179691762_dp)
    call check_scalar('circstd basic', circstd(angles), &
        0.46422774855070464_dp)
    call check_scalar('circstd normalize', circstd(angles, normalize=.true.), &
        0.46422774855070464_dp)
    call check_scalar('circmean_all basic', circmean_all(angles), &
        0.12732136153446286_dp)
    call check_scalar('circvar_all basic', circvar_all(angles), &
        0.10215129179691762_dp)
    call check_scalar('circstd_all basic', circstd_all(angles), &
        0.46422774855070464_dp)

    degrees = [-180.0_dp, -170.0_dp, 170.0_dp]
    call check_scalar('circmean custom range', &
        circmean(degrees, high=180.0_dp, low=-180.0_dp), -180.0_dp)
    call check_scalar('circvar custom range', &
        circvar(degrees, high=180.0_dp, low=-180.0_dp), &
        0.01012816465852795_dp)
    call check_scalar('circstd custom range', &
        circstd(degrees, high=180.0_dp, low=-180.0_dp), &
        8.17536831253456_dp)
    call check_scalar('circstd normalized custom range', &
        circstd(degrees, high=180.0_dp, low=-180.0_dp, normalize=.true.), &
        0.14268709461694087_dp)
    call check_scalar('circmean negative period', &
        circmean([1.0_dp, 2.0_dp], high=0.0_dp, low=2.0_dp), 1.5_dp)
    call check_scalar('circstd negative period', &
        circstd([1.0_dp, 2.0_dp], high=0.0_dp, low=2.0_dp), &
        -2.7504573645368673_dp)

    matrix = reshape([0.0_dp, scifort_pi / 4.0_dp, scifort_pi / 2.0_dp, &
        scifort_pi, scifort_pi, 7.0_dp * scifort_pi / 4.0_dp], shape(matrix))
    a1 = circmean(matrix, axis=1)
    call check_array('circmean matrix axis1', a1, &
        [0.39269908169872414_dp, 2.356194490192345_dp, &
         4.319689898685965_dp])
    a1 = circmean(matrix, axis=2)
    call check_angle_array('circmean matrix axis2', a1, &
        [1.5707963267948966_dp, 6.283185307179586_dp])
    a1 = circvar(matrix, axis=1)
    call check_array('circvar matrix axis1', a1, &
        [0.07612046748871326_dp, 0.2928932188134524_dp, &
         0.6173165676349102_dp])
    a1 = circstd(matrix, axis=2)
    call check_array('circstd matrix axis2', a1, &
        [1.4823038073675108_dp, 1.9899677764665702_dp])
    call check_scalar('circmean matrix all', circmean(matrix), &
        1.1780972450961726_dp)
    call check_scalar('circvar matrix all', circvar(matrix), &
        0.8196012999512676_dp)
    call check_scalar('circstd matrix all', circstd(matrix), &
        1.8507219549774987_dp)

    keep_mean = circmean_keepdims(matrix, axis=2)
    if (any(shape(keep_mean) /= [2, 1])) then
        print *, 'FAIL circmean keepdims shape', shape(keep_mean)
        error stop 1
    end if
    call check_angle_array('circmean keepdims values', reshape(keep_mean, [size(keep_mean)]), &
        [1.5707963267948966_dp, 6.283185307179586_dp])
    keep_var = circvar_keepdims(matrix, axis=1)
    if (any(shape(keep_var) /= [1, 3])) then
        print *, 'FAIL circvar keepdims shape', shape(keep_var)
        error stop 1
    end if
    call check_array('circvar keepdims values', reshape(keep_var, [size(keep_var)]), &
        [0.07612046748871326_dp, 0.2928932188134524_dp, &
         0.6173165676349102_dp])
    keep_std = circstd_keepdims(matrix, axis=2)
    if (any(shape(keep_std) /= [2, 1])) then
        print *, 'FAIL circstd keepdims shape', shape(keep_std)
        error stop 1
    end if
    call check_array('circstd keepdims values', reshape(keep_std, [size(keep_std)]), &
        [1.4823038073675108_dp, 1.9899677764665702_dp])

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    nan_values = [0.0_dp, qnan, scifort_pi / 2.0_dp]
    if (.not. ieee_is_nan(circmean(nan_values, nan_policy='propagate'))) then
        print *, 'FAIL circmean propagate'
        error stop 1
    end if
    if (.not. ieee_is_nan(circvar(nan_values, nan_policy='propagate'))) then
        print *, 'FAIL circvar propagate'
        error stop 1
    end if
    if (.not. ieee_is_nan(circstd(nan_values, nan_policy='propagate'))) then
        print *, 'FAIL circstd propagate'
        error stop 1
    end if
    call check_scalar('circmean omit', circmean(nan_values, nan_policy='omit'), &
        0.7853981633974483_dp)
    call check_scalar('circvar omit', circvar(nan_values, nan_policy='omit'), &
        0.2928932188134524_dp)
    call check_scalar('circstd omit', circstd(nan_values, nan_policy='omit'), &
        0.8325546111576977_dp)
    if (.not. ieee_is_nan(circmean([qnan], nan_policy='omit'))) then
        print *, 'FAIL circmean all-NaN omit'
        error stop 1
    end if
    if (.not. ieee_is_nan(circmean([real(dp) ::]))) then
        print *, 'FAIL circmean empty'
        error stop 1
    end if
    if (.not. ieee_is_nan(circvar([real(dp) ::]))) then
        print *, 'FAIL circvar empty'
        error stop 1
    end if
    if (.not. ieee_is_nan(circstd([real(dp) ::]))) then
        print *, 'FAIL circstd empty'
        error stop 1
    end if
    if (.not. ieee_is_nan(circmean(angles, nan_policy='unsupported'))) then
        print *, 'FAIL circular unsupported policy'
        error stop 1
    end if

    x15 = 0.0_dp
    x15(:, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = &
        [0.0_dp, scifort_pi / 2.0_dp, scifort_pi]
    call check_scalar('circmean rank15 all', circmean_all(x15), &
        scifort_pi / 2.0_dp)
    call check_scalar('circvar rank15 all', circvar_all(x15), &
        2.0_dp / 3.0_dp)
    call check_scalar('circstd rank15 all', circstd_all(x15), &
        1.4823038073675108_dp)
    keep15 = circmean_keepdims(x15, axis=1)
    if (any(shape(keep15) /= 1)) then
        print *, 'FAIL circmean rank15 keepdims shape', shape(keep15)
        error stop 1
    end if
    call check_scalar('circmean rank15 keepdims', &
        keep15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), &
        scifort_pi / 2.0_dp)

    print *, 'descriptive circular tests passed'

contains

    subroutine check_angle_array(label, actual, expected)
        implicit none
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual(:) !! computed angles in the default range
        real(dp), intent(in) :: expected(:) !! reference angles in radians

        ! Both endpoints represent the same direction; compare on the circle.
        if (any(ieee_is_nan(actual)) .or. any(actual < 0.0_dp) .or. &
                any(actual > 2.0_dp * scifort_pi)) then
            print *, 'FAIL ', trim(label), ' invalid angles'
            error stop 1
        end if
        call check_array(label, &
            modulo(actual - expected + scifort_pi, 2.0_dp * scifort_pi) - scifort_pi, &
            spread(0.0_dp, 1, size(expected)))
    end subroutine check_angle_array

    subroutine check_scalar(label, actual, expected)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value

        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_scalar

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label !! diagnostic label
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), ' size', size(actual), size(expected)
            error stop 1
        end if
        if (any(abs(actual - expected) > tol * max(1.0_dp, abs(expected)))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

end program test_descriptive_circular
