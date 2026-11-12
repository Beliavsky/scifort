! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0 by
! tools/generate_statistics_reference.py.

program test_descriptive_entropy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : differential_entropy, differential_entropy_all, &
        differential_entropy_keepdims, entropy, entropy_all, entropy_keepdims
    implicit none

    real(dp) :: pk(2, 3)
    real(dp) :: qv(3)
    real(dp) :: npk(2, 3)
    real(dp) :: nq(2, 3)
    real(dp) :: qnan
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: keep2(:, :)
    real(dp) :: x(11)
    real(dp) :: mat(11, 2)
    real(dp) :: matn(11, 2)
    real(dp) :: pk15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: q15(1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: e15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp) :: d15(11, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: de15(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    integer :: i
    real(dp) :: scalar

    qnan = ieee_value(0.0_dp, ieee_quiet_nan)
    pk = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], shape(pk))
    qv = [1.0_dp, 2.0_dp, 4.0_dp]

    call check_array('entropy axis1', entropy(pk, axis=1), &
        [0.6365141682948128_dp, 0.6829081047004717_dp, 0.6890092384766586_dp])
    call check_array('entropy axis2', entropy(pk, axis=2), &
        [0.9368883075390159_dp, 1.0114042647073518_dp])
    call check_array('entropy qk broadcast', entropy(pk, qk=qv, axis=2), &
        [0.00780924737415419_dp, 0.01030964360136785_dp])
    call check_array('entropy qk base2', entropy(pk, qk=qv, axis=2, base=2.0_dp), &
        [0.01126636245976741_dp, 0.01487367169702603_dp])
    call check_close('entropy all', entropy_all(pk), 1.6623769591928221_dp)

    keep2 = entropy_keepdims(pk, axis=2)
    call check_shape2('entropy keepdims', shape(keep2), [2, 1])
    call check_array('entropy keepdims values', reshape(keep2, [size(keep2)]), &
        [0.9368883075390159_dp, 1.0114042647073518_dp])

    npk = pk
    npk(1, 2) = qnan
    nq = spread(qv, 1, 2)
    nq(2, 3) = qnan
    call check_array('entropy omit', entropy(npk, axis=2, nan_policy='omit'), &
        [0.45056120886630463_dp, 1.0114042647073518_dp])
    call check_array('entropy paired omit', &
        entropy(npk, qk=nq, axis=2, nan_policy='omit'), &
        [0.00363140263455348_dp, 0.0_dp])
    values = entropy(npk, qk=nq, axis=2, nan_policy='propagate')
    call check_true('entropy propagate', all(ieee_is_nan(values)))
    call check_close('entropy omit single NaN', &
        entropy([qnan], nan_policy='omit'), 0.0_dp)
    call check_nan('entropy zero total', entropy([0.0_dp]))
    scalar = entropy([-1.0_dp, 2.0_dp])
    call check_true('entropy negative probability', .not. ieee_is_finite(scalar) .and. scalar < 0.0_dp)
    scalar = entropy([1.0_dp, 1.0_dp], qk=[0.0_dp, 2.0_dp])
    call check_true('entropy infinite KL', .not. ieee_is_finite(scalar) .and. scalar > 0.0_dp)
    scalar = entropy([0.5_dp, 0.5_dp], base=1.0_dp)
    call check_true('entropy base one infinity', .not. ieee_is_finite(scalar) .and. scalar > 0.0_dp)
    call check_nan('entropy zero base one', entropy([1.0_dp], base=1.0_dp))

    x = [0.2_dp, 0.5_dp, 1.1_dp, 1.7_dp, 2.0_dp, 3.4_dp, 4.1_dp, 5.6_dp, &
        7.0_dp, 8.5_dp, 10.0_dp]
    call check_close('differential vasicek', &
        differential_entropy(x, window_length=2, method='vasicek'), &
        2.1030589542982048_dp)
    call check_close('differential van es', &
        differential_entropy(x, window_length=2, method='van es'), &
        2.6043948164647985_dp)
    call check_close('differential ebrahimi', &
        differential_entropy(x, window_length=2, method='ebrahimi'), &
        2.281391545754883_dp)
    call check_close('differential correa', &
        differential_entropy(x, window_length=2, method='correa'), &
        2.193209238265152_dp)
    call check_close('differential auto', &
        differential_entropy(x, window_length=2, method='auto'), &
        2.281391545754883_dp)
    call check_close('differential auto default window', differential_entropy(x), &
        2.2991749319564807_dp)
    call check_close('differential base2', &
        differential_entropy(x, window_length=2, method='ebrahimi', base=2.0_dp), &
        3.291352269386576_dp)
    call check_close('differential small auto', &
        differential_entropy([1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp, 11.0_dp, 16.0_dp, 22.0_dp]), &
        3.388342367220764_dp)
    call check_nan('differential too small', &
        differential_entropy([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]))
    call check_close('differential duplicates trap safe', &
        differential_entropy([1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
            window_length=2, method='correa'), 1.1826573532087914_dp)

    do i = 1, size(x)
        mat(i, 1) = x(i)
        mat(i, 2) = x(i) ** 1.1_dp
    end do
    call check_array('differential matrix axis1', &
        differential_entropy(mat, axis=1, window_length=2, method='ebrahimi'), &
        [2.281391545754883_dp, 2.479488707097673_dp])
    keep2 = differential_entropy_keepdims(mat, axis=1, window_length=2, &
        method='ebrahimi')
    call check_shape2('differential keepdims', shape(keep2), [1, 2])
    call check_array('differential keepdims values', reshape(keep2, [size(keep2)]), &
        [2.281391545754883_dp, 2.479488707097673_dp])
    call check_close('differential all', &
        differential_entropy_all(reshape(x, [11, 1]), window_length=2, method='ebrahimi'), &
        2.281391545754883_dp)

    matn = mat
    matn(2, 1) = qnan
    matn(6, 2) = qnan
    call check_array('differential omit', &
        differential_entropy(matn, axis=1, window_length=2, method='ebrahimi', &
            nan_policy='omit'), &
        [2.333440488205306_dp, 2.466855627516598_dp])
    values = differential_entropy(matn, axis=1, window_length=2, method='ebrahimi')
    call check_true('differential propagate', all(ieee_is_nan(values)))

    pk15(:, :, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = pk
    q15(1, :, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = qv
    e15 = entropy(pk15, qk=q15, axis=2)
    call check_array('entropy rank15', reshape(e15, [size(e15)]), &
        [0.00780924737415419_dp, 0.01030964360136785_dp])

    d15(:, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = x
    d15(:, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = mat(:, 2)
    de15 = differential_entropy(d15, axis=1, window_length=2, method='ebrahimi')
    call check_array('differential rank15', reshape(de15, [size(de15)]), &
        [2.281391545754883_dp, 2.479488707097673_dp])

    print *, 'descriptive entropy tests passed'

contains

    subroutine check_array(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        integer :: j

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), size(actual), size(expected)
            error stop 1
        end if
        do j = 1, size(actual)
            if (ieee_is_nan(expected(j))) then
                if (.not. ieee_is_nan(actual(j))) then
                    print *, 'FAIL ', trim(label), j, actual(j), expected(j)
                    error stop 1
                end if
            else if (abs(actual(j) - expected(j)) > 2.0e-12_dp) then
                print *, 'FAIL ', trim(label), j, actual(j), expected(j)
                error stop 1
            end if
        end do
    end subroutine check_array

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected

        if (abs(actual - expected) > 2.0e-12_dp) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_nan(label, actual)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual

        if (.not. ieee_is_nan(actual)) then
            print *, 'FAIL ', trim(label), actual
            error stop 1
        end if
    end subroutine check_nan

    subroutine check_shape2(label, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: actual(2)
        integer, intent(in) :: expected(2)

        if (any(actual /= expected)) then
            print *, 'FAIL ', trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_shape2

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', trim(label)
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_entropy
