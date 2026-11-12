! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0.

program test_descriptive_rank3_weights
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, pmean, pmean_keepdims
    implicit none

    real(dp), parameter :: tol = 4.0e-14_dp
    integer :: i
    real(dp) :: nan
    real(dp) :: w1(2)
    real(dp) :: w2(3, 2)
    real(dp) :: w3(4, 3, 2)
    real(dp) :: w31(4, 1, 1)
    real(dp) :: wbad(2, 2)
    real(dp) :: wn(2)
    real(dp) :: x(4, 3, 2)
    real(dp) :: xn(4, 3, 2)
    real(dp) :: xs(4, 1, 2)
    real(dp), allocatable :: kept(:, :, :)
    real(dp), allocatable :: values(:, :)

    x = reshape([ &
        1.0_dp, 2.0_dp, 4.0_dp, 9.0_dp, &
        3.0_dp, 5.0_dp, 8.0_dp, 15.0_dp, &
        2.0_dp, 6.0_dp, 11.0_dp, 20.0_dp, &
        4.0_dp, 7.0_dp, 13.0_dp, 24.0_dp, &
        5.0_dp, 10.0_dp, 18.0_dp, 31.0_dp, &
        6.0_dp, 12.0_dp, 22.0_dp, 40.0_dp], [4, 3, 2])
    w1 = [1.0_dp, 3.0_dp]
    w2 = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [3, 2])
    w3 = reshape([(real(i, dp), i=1, 24)], [4, 3, 2])
    w31(:, 1, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]

    values = pmean(x, 2.0_dp, axis=3, weights=w1)
    call check_matrix('rank3 pmean rank1 weights', values, reshape([ &
        3.5_dp, 6.1441028637222539_dp, 11.434596626029272_dp, &
        21.266170318136737_dp, 4.5825756949558398_dp, 9.013878188659973_dp, &
        16.093476939431081_dp, 27.874719729532707_dp, 5.2915026221291814_dp, &
        10.816653826391969_dp, 19.83053201505194_dp, 36.055512754639892_dp], &
        [4, 3]), tol)

    values = hmean(x, axis=2, weights=w2)
    call check_matrix('rank3 hmean rank2 weights', values, reshape([ &
        1.8947368421052631_dp, 4.2857142857142865_dp, 7.7647058823529411_dp, &
        15.2112676056338_dp, 5.0_dp, 9.545454545454545_dp, &
        17.478497057492078_dp, 31.38357705286839_dp], [4, 2]), tol)

    values = gmean(x, axis=1, weights=w3)
    call check_matrix('rank3 gmean rank3 weights', values, reshape([ &
        4.1929627126294751_dp, 7.212215636379951_dp, 7.8387911572775595_dp, &
        10.180530037907182_dp, 13.464426463534503_dp, 16.430396988404798_dp], &
        [3, 2]), tol)

    values = gmean(x, axis=1, weights=w31)
    call check_matrix('rank3 gmean singleton rank3 weights', values, reshape([ &
        4.1929627126294751_dp, 8.4892559178579319_dp, 10.43669647896939_dp, &
        13.046472690280048_dp, 17.499410926238763_dp, 21.737128614791271_dp], &
        [3, 2]), tol)

    kept = gmean_keepdims(x, axis=2, weights=w2)
    call check_tensor('rank3 weighted keepdims', kept, reshape([ &
        2.0396489026555056_dp, 4.701509224902396_dp, 8.3573707759438776_dp, &
        15.90691043510728_dp, 5.0675708886918915_dp, 9.7805957500348857_dp, &
        17.88318219995228_dp, 32.06275911849341_dp], [4, 1, 2]), tol)

    xs = x(:, 1:1, :)
    values = gmean(xs, axis=1, weights=w2)
    call check_matrix('rank3 weights expand retained axis', values, reshape([ &
        2.9129506302439405_dp, 2.9129506302439405_dp, 2.9129506302439405_dp, &
        9.6678111051252635_dp, 9.6678111051252635_dp, 9.6678111051252671_dp], &
        [3, 2]), tol)

    values = hmean(xs, axis=2, weights=w2)
    call check_matrix('rank3 weights expand reduced axis', values, reshape([ &
        1.0_dp, 2.0_dp, 4.0_dp, 9.0_dp, 4.0_dp, 7.0_dp, 13.0_dp, 24.0_dp], &
        [4, 2]), tol)

    values = pmean(xs, 2.0_dp, axis=3, weights=w2)
    call check_matrix('rank3 weights expand retained dimension', values, reshape([ &
        3.6055512754639891_dp, 6.324555320336759_dp, 11.764352935882194_dp, &
        21.840329667841555_dp, 3.4226138716316967_dp, 6.0118929750002321_dp, &
        11.193110125685099_dp, 20.846376869169102_dp, 3.3166247903553998_dp, &
        5.8309518948453007_dp, 10.862780491200215_dp, 20.273134932713294_dp], &
        [4, 3]), tol)

    kept = pmean_keepdims(xs, 2.0_dp, axis=3, weights=w2)
    call check_tensor('rank3 expanded weighted keepdims', kept, reshape([ &
        3.6055512754639891_dp, 6.324555320336759_dp, 11.764352935882194_dp, &
        21.840329667841555_dp, 3.4226138716316967_dp, 6.0118929750002321_dp, &
        11.193110125685099_dp, 20.846376869169102_dp, 3.3166247903553998_dp, &
        5.8309518948453007_dp, 10.862780491200215_dp, 20.273134932713294_dp], &
        [4, 3, 1]), tol)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    xn = x
    xn(1, 1, 1) = nan
    wn = [1.0_dp, nan]
    values = gmean(xn, axis=3, weights=wn, nan_policy='omit')
    call check_matrix_nan('rank3 weighted gmean omit', values, reshape([ &
        nan, 2.0_dp, 4.0_dp, 9.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 15.0_dp, &
        2.0_dp, 6.0_dp, 11.0_dp, 20.0_dp], [4, 3]), tol)

    values = pmean(xn, 2.0_dp, axis=3, weights=wn, nan_policy='omit')
    call check_matrix_nan('rank3 weighted pmean omit', values, reshape([ &
        nan, 2.0_dp, 4.0_dp, 9.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 15.0_dp, &
        2.0_dp, 6.0_dp, 11.0_dp, 20.0_dp], [4, 3]), tol)

    values = hmean(x, axis=3, weights=wn, nan_policy='propagate')
    call check_true('rank3 weighted propagate NaN weights', all(ieee_is_nan(values)))

    wbad = 1.0_dp
    values = gmean(x, axis=1, weights=wbad)
    call check_true('rank3 incompatible weights shape', all(shape(values) == [3, 2]))
    call check_true('rank3 incompatible weights values', all(ieee_is_nan(values)))

    kept = gmean_keepdims(x, axis=1, weights=wbad)
    call check_true('rank3 incompatible keepdims shape', all(shape(kept) == [1, 3, 2]))
    call check_true('rank3 incompatible keepdims values', all(ieee_is_nan(kept)))

    call check_true('rank3 weighted invalid axis', &
        all(shape(gmean(x, axis=4, weights=w2)) == [0, 0]))

    print *, 'descriptive rank-3 weighted tests passed'

contains

    subroutine check_matrix(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :)
        real(dp), intent(in) :: expected(:, :)
        real(dp), intent(in) :: tolerance

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', name, ' shape ', shape(actual), shape(expected)
            error stop 1
        end if
        if (any(ieee_is_nan(actual))) then
            print *, 'FAIL ', name, ' unexpected NaN'
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_matrix

    subroutine check_matrix_nan(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :)
        real(dp), intent(in) :: expected(:, :)
        real(dp), intent(in) :: tolerance
        integer :: i
        integer :: j

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', name, ' shape ', shape(actual), shape(expected)
            error stop 1
        end if
        do j = 1, size(actual, 2)
            do i = 1, size(actual, 1)
                if (ieee_is_nan(expected(i, j))) then
                    if (.not. ieee_is_nan(actual(i, j))) then
                        print *, 'FAIL ', name, ' expected NaN'
                        error stop 1
                    end if
                else
                    if (ieee_is_nan(actual(i, j))) then
                        print *, 'FAIL ', name, ' unexpected NaN'
                        error stop 1
                    end if
                    if (abs(actual(i, j) - expected(i, j)) > tolerance) then
                        print *, 'FAIL ', name
                        error stop 1
                    end if
                end if
            end do
        end do
    end subroutine check_matrix_nan

    subroutine check_tensor(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :, :)
        real(dp), intent(in) :: expected(:, :, :)
        real(dp), intent(in) :: tolerance

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', name, ' shape ', shape(actual), shape(expected)
            error stop 1
        end if
        if (any(ieee_is_nan(actual))) then
            print *, 'FAIL ', name, ' unexpected NaN'
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_tensor

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_rank3_weights
