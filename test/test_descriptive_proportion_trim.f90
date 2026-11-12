! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_descriptive_proportion_trim
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : trim1, trim1_all, trim_mean, trim_mean_all, &
        trim_mean_keepdims, trimboth, trimboth_all
    implicit none

    real(dp), allocatable :: empty(:, :)
    real(dp) :: matrix(2, 4)
    real(dp) :: nan
    real(dp) :: rank15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: r1(:)
    real(dp), allocatable :: r2(:, :)
    real(dp), allocatable :: r15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: x10(10)
    real(dp) :: x20(20)
    integer :: i

    x10 = [(real(i - 1, dp), i=1,10)]
    x20 = [(real(i - 1, dp), i=1,20)]
    nan = ieee_value(0.0_dp, ieee_quiet_nan)

    call check_close(trim_mean(x20, 0.1_dp), 9.5_dp, 'trim_mean vector')
    call check_close(trim_mean([1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp], 0.25_dp), &
        2.5_dp, 'trim_mean quarter')
    call check_close(trim_mean([1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp], 0.24999_dp), &
        2.75_dp, 'trim_mean floor cut')
    call check_nan(trim_mean([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.5_dp), &
        'trim_mean exact half is empty')

    matrix = reshape([1.0_dp, 10.0_dp, 2.0_dp, 20.0_dp, 3.0_dp, 30.0_dp, &
        5.0_dp, 50.0_dp], shape(matrix))
    r1 = trim_mean(matrix, 0.25_dp)
    call check_vector(r1, [5.5_dp, 11.0_dp, 16.5_dp, 27.5_dp], &
        'trim_mean matrix default axis1')
    r1 = trim_mean(matrix, 0.25_dp, axis=2)
    call check_vector(r1, [2.5_dp, 25.0_dp], 'trim_mean matrix axis2')
    call check_close(trim_mean_all(matrix, 0.25_dp), 9.5_dp, 'trim_mean all')
    r2 = trim_mean_keepdims(matrix, 0.25_dp, axis=2)
    if (any(shape(r2) /= [2, 1])) error stop 'trim_mean keepdims shape'
    call check_vector(reshape(r2, [size(r2)]), [2.5_dp, 25.0_dp], &
        'trim_mean keepdims values')

    call check_nan(trim_mean([1.0_dp, nan, 3.0_dp, 5.0_dp], 0.25_dp, &
        nan_policy='propagate'), 'trim_mean nan propagate')
    call check_close(trim_mean([1.0_dp, nan, 3.0_dp, 5.0_dp], 0.25_dp, &
        nan_policy='omit'), 3.0_dp, 'trim_mean nan omit')

    r1 = trimboth(x10, 0.1_dp)
    call check_vector(r1, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
        7.0_dp, 8.0_dp], 'trimboth vector')
    r1 = trimboth(x10, 0.25_dp)
    call check_vector(r1, [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp], &
        'trimboth floor cut')
    r1 = trim1(x10, 0.2_dp, tail='left')
    call check_vector(r1, [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, &
        8.0_dp, 9.0_dp], 'trim1 left')
    r1 = trim1(x10, 0.2_dp, tail='RIGHT')
    call check_vector(r1, [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
        6.0_dp, 7.0_dp], 'trim1 right')
    r1 = trim1(x10, 1.0_dp)
    if (size(r1) /= 0) error stop 'trim1 proportion one should be empty'

    r2 = trimboth(matrix, 0.25_dp, axis=2)
    if (any(shape(r2) /= [2, 2])) error stop 'trimboth matrix axis2 shape'
    call check_vector(reshape(r2, [size(r2)]), [2.0_dp, 20.0_dp, 3.0_dp, 30.0_dp], &
        'trimboth matrix axis2')
    r2 = trim1(matrix, 0.5_dp, tail='left', axis=1)
    if (any(shape(r2) /= [1, 4])) error stop 'trim1 matrix axis1 shape'
    call check_vector(reshape(r2, [size(r2)]), [10.0_dp, 20.0_dp, 30.0_dp, 50.0_dp], &
        'trim1 matrix axis1')

    r1 = trimboth_all(matrix, 0.25_dp)
    call check_vector(r1, [3.0_dp, 5.0_dp, 10.0_dp, 20.0_dp], 'trimboth all')
    r1 = trim1_all(matrix, 0.25_dp, tail='right')
    call check_vector(r1, [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 10.0_dp, 20.0_dp], &
        'trim1 all right')

    allocate(empty(0, 3))
    r2 = trimboth(empty, 0.1_dp, axis=1)
    if (any(shape(r2) /= [0, 3])) error stop 'trimboth empty shape'
    r2 = trim1(empty, 0.1_dp, axis=2)
    if (any(shape(r2) /= [0, 3])) error stop 'trim1 empty shape'
    r1 = trim_mean(empty, 0.1_dp, axis=1)
    if (size(r1) /= 3 .or. .not. all(ieee_is_nan(r1))) then
        error stop 'trim_mean empty reduction'
    end if

    r1 = trimboth([1.0_dp, 2.0_dp, nan, 4.0_dp, 5.0_dp], 0.2_dp)
    call check_vector(r1, [2.0_dp, 4.0_dp, 5.0_dp], 'trimboth NaN sorts high')
    r1 = trim1([1.0_dp, 2.0_dp, nan, 4.0_dp, 5.0_dp], 0.2_dp, tail='left')
    if (size(r1) /= 4 .or. .not. ieee_is_nan(r1(4))) then
        error stop 'trim1 left should retain NaN as upper extreme'
    end if

    rank15 = 0.0_dp
    rank15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = 1.0_dp
    rank15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = 3.0_dp
    call check_close(sum(trim_mean(rank15, 0.25_dp, axis=1)), 2.0_dp, &
        'trim_mean rank15')
    r15 = trim1(rank15, 0.5_dp, tail='left', axis=1)
    if (size(r15, dim=1) /= 1) error stop 'trim1 rank15 shape'
    call check_close(sum(r15), 3.0_dp, 'trim1 rank15')

    print *, 'descriptive proportion trimming tests passed'

contains

    subroutine check_close(actual, expected, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        character(len=*), intent(in) :: label

        if (ieee_is_nan(actual) .or. abs(actual - expected) > 5.0e-13_dp) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_nan(actual, label)
        real(dp), intent(in) :: actual
        character(len=*), intent(in) :: label

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', label
            print *, ' actual  = ', actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_vector(actual, expected, label)
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        character(len=*), intent(in) :: label
        integer :: j

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', label, ' size'
            error stop 1
        end if
        do j = 1, size(actual)
            call check_close(actual(j), expected(j), label)
        end do
    end subroutine check_vector

end program test_descriptive_proportion_trim
