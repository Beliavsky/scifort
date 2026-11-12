! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_power_transform
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : boxcox, boxcox_fit, boxcox_llf, boxcox_llf_all, &
        boxcox_llf_keepdims, boxcox_normmax, boxcox_normmax_all, yeojohnson, &
        yeojohnson_fit, yeojohnson_llf, yeojohnson_llf_all, yeojohnson_llf_keepdims, &
        yeojohnson_normmax
    implicit none

    real(dp), parameter :: tol = 8.0e-8_dp
    real(dp) :: x(6)
    real(dp) :: y(6)
    real(dp) :: a(3, 2)
    real(dp) :: b(3, 2)
    real(dp) :: qnan
    real(dp) :: lmbda
    real(dp), allocatable :: transformed(:)
    real(dp), allocatable :: reduced(:)
    real(dp), allocatable :: kept(:, :)
    real(dp) :: rank15(3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: reduced15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    y = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
    a = reshape(x, shape(a))
    b = reshape(y, shape(b))

    call check_array('boxcox transform', boxcox(x, 0.5_dp), &
        [0.0_dp, 0.8284271247461901_dp, 1.4641016151377548_dp, 2.0_dp, &
         2.4721359549995796_dp, 2.8989794855663558_dp])
    call check_close('boxcox llf', boxcox_llf(0.5_dp, x), -3.178122278370899_dp)
    call check_close('boxcox llf log', boxcox_llf(0.0_dp, x), -3.5631776435643716_dp)
    reduced = boxcox_llf(0.5_dp, a, axis=1)
    call check_array('boxcox axis1', reduced, &
        [0.6393927503112551_dp, 0.6121098857023859_dp])
    kept = boxcox_llf_keepdims(0.5_dp, a, axis=1)
    if (any(shape(kept) /= [1, 2])) then
        print *, 'FAIL boxcox keepdims shape', shape(kept)
        error stop 1
    end if
    call check_array('boxcox keepdims', reshape(kept, [size(kept)]), &
        [0.6393927503112551_dp, 0.6121098857023859_dp])
    call check_close('boxcox all', boxcox_llf_all(0.5_dp, a), -3.178122278370899_dp)
    call check_close('boxcox omit', &
        boxcox_llf(0.5_dp, [1.0_dp, 2.0_dp, qnan, 4.0_dp], nan_policy='omit'), &
        -0.446169990728709_dp)
    call check_nan('boxcox propagate', &
        boxcox_llf(0.5_dp, [1.0_dp, qnan, 3.0_dp], nan_policy='propagate'))

    call check_close('boxcox normmax pearsonr', boxcox_normmax(x, 'pearsonr'), &
        0.9529137774982831_dp)
    call check_close('boxcox normmax mle', boxcox_normmax(x, 'mle'), &
        0.7045411616342833_dp)
    call check_array('boxcox normmax all', boxcox_normmax_all(x), &
        [0.9529137774982831_dp, 0.7045411616342833_dp])
    call boxcox_fit(x, transformed, lmbda)
    call check_close('boxcox fit lambda', lmbda, 0.7045411616342833_dp)
    call check_array('boxcox fit values', transformed, boxcox(x, lmbda))

    call check_array('yeojohnson transform', yeojohnson(y, 0.5_dp), &
        [-2.797434948471088_dp, -1.2189514164974595_dp, 0.0_dp, &
         0.8284271247461897_dp, 1.4641016151377548_dp, 2.4721359549995796_dp])
    call check_array('yeojohnson lambda0', yeojohnson(y, 0.0_dp), &
        [-4.0_dp, -1.5_dp, 0.0_dp, 0.6931471805599453_dp, &
         1.0986122886681098_dp, 1.6094379124341003_dp])
    call check_close('yeojohnson llf', yeojohnson_llf(0.5_dp, y), &
        -4.127023057471658_dp)
    reduced = yeojohnson_llf(0.5_dp, b, axis=1)
    call check_array('yeojohnson axis1', reduced, &
        [0.4891739091432443_dp, -0.5292673675698518_dp])
    kept = yeojohnson_llf_keepdims(0.5_dp, b, axis=1)
    if (any(shape(kept) /= [1, 2])) then
        print *, 'FAIL yeojohnson keepdims shape', shape(kept)
        error stop 1
    end if
    call check_close('yeojohnson all', yeojohnson_llf_all(0.5_dp, b), &
        -4.127023057471658_dp)
    call check_close('yeojohnson omit', &
        yeojohnson_llf(0.5_dp, [-2.0_dp, qnan, 0.0_dp, 2.0_dp], nan_policy='omit'), &
        -1.7094143197886857_dp)
    call check_close('yeojohnson normmax', yeojohnson_normmax(y), &
        0.7684649769684502_dp)
    call yeojohnson_fit(y, transformed, lmbda)
    call check_close('yeojohnson fit lambda', lmbda, 0.7684649769684502_dp)
    call check_array('yeojohnson fit values', transformed, yeojohnson(y, lmbda))

    rank15 = 0.0_dp
    rank15(:, :, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = a
    reduced15 = boxcox_llf(0.5_dp, rank15, axis=1)
    if (size(reduced15) /= 2) then
        print *, 'FAIL power transform rank15 size', size(reduced15)
        error stop 1
    end if
    call check_array('boxcox rank15', reshape(reduced15, [size(reduced15)]), &
        [0.6393927503112551_dp, 0.6121098857023859_dp])

    print *, 'power-transform tests passed'

contains

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected

        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)

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

    subroutine check_nan(label, actual)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', trim(label), actual
            error stop 1
        end if
    end subroutine check_nan

end program test_power_transform
