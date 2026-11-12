! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0.

program test_descriptive_rank4
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, iqr, iqr_keepdims, &
        kurtosis, median_abs_deviation, moment, pmean, sem, sem_keepdims, skew, variation
    implicit none

    real(dp), parameter :: tol = 5.0e-12_dp
    integer :: n
    real(dp) :: nan
    real(dp) :: empty(0, 2, 3, 4)
    real(dp) :: x(4, 3, 2, 2)
    real(dp) :: xn(4, 3, 2, 2)
    real(dp), allocatable :: values(:, :, :)
    real(dp), allocatable :: kept(:, :, :, :)

    x = reshape([(0.25_dp * real(n * n, dp) + 0.5_dp * real(n, dp) + 1.0_dp, &
        n = 1, 48)], [4, 3, 2, 2])

    values = gmean(x)
    call check_rank3('rank4 gmean default axis', values, &
        reshape([ &
            3.6348604184650353_dp, 14.531000527154013_dp, 33.513221856124211_dp, &
            60.507281810113561_dp, 95.504601503471847_dp, 138.50316855186051_dp, &
            189.50231390216575_dp, 248.50176358806021_dp, 315.50138856577286_dp, &
            390.50112158510615_dp, 473.50092480505248_dp, 564.50077560955378_dp], [3, 2, 2]), tol)

    values = hmean(x, axis=2)
    call check_rank3('rank4 hmean axis 2', values, &
        reshape([ &
            4.2085910995809224_dp, 6.7794392523364495_dp, 10.086112109234639_dp, &
            14.036269430051814_dp, 73.989176926549092_dp, 83.213903743315498_dp, &
            92.942456041421025_dp, 103.17409640182748_dp, 217.83413704312997_dp, &
            233.07872795818008_dp, 248.82382589145351_dp, 265.06936922836121_dp, &
            433.79263676411011_dp, 455.04066331330125_dp, 476.78882403391395_dp, &
            499.03710701202704_dp], [4, 2, 2]), tol)

    values = pmean(x, 2.0_dp, axis=3)
    call check_rank3('rank4 pmean axis 3', values, &
        reshape([ &
            35.200319600821807_dp, 40.36087214122113_dp, 45.908196435930698_dp, &
            51.855568649856693_dp, 58.215655110975092_dp, 65.0_dp, &
            72.218851417064229_dp, 79.881161734166085_dp, 87.994673134230126_dp, &
            96.566039579139826_dp, 105.60095880246543_dp, 115.10430052782563_dp, &
            282.55806925302983_dp, 298.87288267756912_dp, 315.68269274700504_dp, &
            332.98798777133089_dp, 350.78919952016764_dp, 369.08671067921154_dp, &
            387.88086121900886_dp, 407.17195384751147_dp, 426.96025868926023_dp, &
            447.24601731038365_dp, 468.02944618901915_dp, 489.31073971454992_dp], [4, 3, 2]), tol)

    values = sem(x, axis=4)
    call check_rank3('rank4 sem axis 4', values, &
        reshape([ &
            83.999999999999986_dp, 90.0_dp, 95.999999999999986_dp, &
            101.99999999999999_dp, 108.0_dp, 113.99999999999999_dp, &
            120.0_dp, 125.99999999999999_dp, 132.0_dp, &
            138.0_dp, 144.0_dp, 150.0_dp, &
            155.99999999999997_dp, 162.0_dp, 167.99999999999997_dp, &
            174.0_dp, 180.0_dp, 186.0_dp, &
            191.99999999999997_dp, 197.99999999999997_dp, 203.99999999999997_dp, &
            210.0_dp, 216.0_dp, 221.99999999999997_dp], [4, 3, 2]), tol)

    values = variation(x, axis=2)
    call check_rank3('rank4 variation axis 2', values, &
        reshape([ &
            0.8035774253474236_dp, 0.73949562540352198_dp, 0.6797936437158224_dp, &
            0.62608762734909484_dp, 0.34891564911154227_dp, 0.33185888190390794_dp, &
            0.31633441795675221_dp, 0.30215241929363656_dp, 0.21463438232295057_dp, &
            0.20789833369644731_dp, 0.2015668132464126_dp, 0.19560491176678299_dp, &
            0.15438581307654056_dp, 0.15084613984391521_dp, 0.1474640530354569_dp, &
            0.14422932791456408_dp], [4, 2, 2]), tol)

    call check_close('rank4 iqr all', iqr(x), 299.625_dp, tol)
    values = iqr(x, axis=4)
    call check_rank3('rank4 iqr axis 4', values, &
        reshape([ &
            84.0_dp, 90.0_dp, 96.0_dp, &
            102.0_dp, 108.0_dp, 114.0_dp, &
            120.0_dp, 126.0_dp, 132.0_dp, &
            138.0_dp, 144.0_dp, 150.0_dp, &
            156.0_dp, 162.0_dp, 168.0_dp, &
            174.0_dp, 180.0_dp, 186.0_dp, &
            192.0_dp, 198.0_dp, 204.0_dp, &
            210.0_dp, 216.0_dp, 222.0_dp], [4, 3, 2]), tol)

    values = median_abs_deviation(x, axis=1)
    call check_rank3('rank4 mad axis 1', values, &
        reshape([ &
            1.5_dp, 3.5_dp, 5.5_dp, &
            7.5_dp, 9.5_dp, 11.5_dp, &
            13.5_dp, 15.5_dp, 17.5_dp, &
            19.5_dp, 21.5_dp, 23.5_dp], [3, 2, 2]), tol)

    values = moment(x, order=2, axis=2)
    call check_rank3('rank4 moment axis 2', values, &
        reshape([ &
            99.555555555555557_dp, 134.22222222222223_dp, 174.2222222222222_dp, &
            219.55555555555557_dp, 867.55555555555554_dp, 966.22222222222217_dp, &
            1070.2222222222224_dp, 1179.5555555555554_dp, 2403.5555555555552_dp, &
            2566.2222222222222_dp, 2734.2222222222222_dp, 2907.5555555555561_dp, &
            4707.5555555555557_dp, 4934.2222222222226_dp, 5166.2222222222217_dp, &
            5403.5555555555557_dp], [4, 2, 2]), tol)

    values = skew(x, axis=1, bias=.false.)
    call check_rank3('rank4 skew axis 1 unbiased', values, &
        reshape([ &
            0.51840469778263387_dp, 0.24655480702578339_dp, 0.16128894592252691_dp, &
            0.11978793854859671_dp, 0.095259809126389094_dp, 0.079064790512753849_dp, &
            0.067574355732838765_dp, 0.058999050226979127_dp, 0.05235463535460317_dp, &
            0.047055053533512123_dp, 0.042729593368786192_dp, 0.039132312615884214_dp], [3, 2, 2]), tol)

    values = kurtosis(x, axis=2, bias=.false.)
    call check_rank3('rank4 kurtosis axis 2 unbiased', values, &
        reshape([ &
            -1.4999999999999998_dp, -1.5_dp, -1.4999999999999996_dp, &
            -1.5000000000000002_dp, -1.5000000000000004_dp, -1.5_dp, &
            -1.5000000000000004_dp, -1.5_dp, -1.4999999999999998_dp, &
            -1.5_dp, -1.4999999999999998_dp, -1.5000000000000007_dp, &
            -1.5000000000000004_dp, -1.5000000000000004_dp, -1.4999999999999998_dp, &
            -1.5000000000000002_dp], [4, 2, 2]), tol)

    kept = gmean_keepdims(x, axis=2)
    call check_rank4('rank4 gmean keepdims axis 2', kept, &
        reshape([ &
            7.6022118158392846_dp, 10.653086015253002_dp, 14.29953521073965_dp, &
            18.489646020197871_dp, 79.113379868812132_dp, 88.360052041622637_dp, &
            98.107256991767287_dp, 108.35488538994863_dp, 223.09354035500201_dp, &
            238.34287256497845_dp, 254.09226875908934_dp, 270.34172095927198_dp, &
            439.08845784568331_dp, 460.33821849510832_dp, 482.0879956074096_dp, &
            504.3377877003237_dp], [4, 1, 2, 2]), tol)

    kept = iqr_keepdims(x)
    call check_rank4('rank4 iqr keepdims all', kept, &
        reshape([299.625_dp], [1, 1, 1, 1]), tol)

    call check_true('rank4 invalid axis', all(shape(sem(x, axis=5)) == [0, 0, 0]))
    call check_true('rank4 invalid keepdims axis', &
        all(shape(gmean_keepdims(x, axis=5)) == [0, 0, 0, 0]))

    values = sem(empty, axis=1)
    call check_true('rank4 empty reduced shape', all(shape(values) == [2, 3, 4]))
    call check_true('rank4 empty reduced NaN', all(ieee_is_nan(values)))
    values = sem(empty, axis=2)
    call check_true('rank4 empty preserved zero extent', all(shape(values) == [0, 3, 4]))
    kept = sem_keepdims(empty, axis=1)
    call check_true('rank4 empty keepdims shape', all(shape(kept) == [1, 2, 3, 4]))
    call check_true('rank4 empty keepdims NaN', all(ieee_is_nan(kept)))

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    xn = x
    xn(1, 1, 1, 1) = nan
    xn(3, 2, 1, 1) = nan
    xn(2, 3, 2, 2) = nan
    xn(4, 1, 2, 2) = nan
    xn(1, 3, 1, 2) = nan

    values = gmean(xn, axis=1, nan_policy='omit')
    call check_rank3('rank4 gmean omit', values, &
        reshape([ &
            4.6377176151082065_dp, 13.858697678150248_dp, 33.513221856124211_dp, &
            60.507281810113561_dp, 95.504601503471847_dp, 138.50316855186051_dp, &
            189.50231390216575_dp, 248.50176358806021_dp, 324.58408255403396_dp, &
            380.83397184248923_dp, 473.50092480505248_dp, 568.38727557410118_dp], [3, 2, 2]), tol)

    values = sem(xn, axis=4, nan_policy='omit')
    call check_rank3_nan('rank4 sem omit', values, &
        reshape([ &
            nan, 90.0_dp, 95.999999999999986_dp, &
            101.99999999999999_dp, 108.0_dp, 113.99999999999999_dp, &
            nan, 125.99999999999999_dp, nan, &
            138.0_dp, 144.0_dp, 150.0_dp, &
            155.99999999999997_dp, 162.0_dp, 167.99999999999997_dp, &
            nan, 180.0_dp, 186.0_dp, &
            191.99999999999997_dp, 197.99999999999997_dp, 203.99999999999997_dp, &
            nan, 216.0_dp, 221.99999999999997_dp], [4, 3, 2]), tol)

    call check_close('rank4 iqr all omit', iqr(xn, nan_policy='omit'), 280.5_dp, tol)

    values = skew(xn, axis=2, bias=.false., nan_policy='omit')
    call check_rank3('rank4 skew omit', values, &
        reshape([ &
            0.0_dp, 0.81584308913709402_dp, 0.0_dp, &
            0.64696859492151881_dp, 0.3308318147058979_dp, 0.31366135050563088_dp, &
            0.29817455220105016_dp, 0.28413676423291157_dp, 0.0_dp, &
            0.19305695181491053_dp, 0.18705316268389463_dp, 0.18141070658606234_dp, &
            0.14265937374189555_dp, 0.0_dp, 0.13619161354325776_dp, &
            0.0_dp], [4, 2, 2]), tol)

    values = variation(xn, axis=1, nan_policy='propagate')
    call check_rank3_nan('rank4 variation propagate', values, &
        reshape([ &
            nan, nan, 0.18852907300304428_dp, &
            0.14181381143219476_dp, 0.11343248652558056_dp, 0.094442249403039816_dp, &
            0.080867850098619326_dp, 0.070690656861089607_dp, nan, &
            nan, 0.051291377417124098_dp, nan], [3, 2, 2]), tol)

    kept = sem_keepdims(xn, axis=4, nan_policy='omit')
    call check_true('rank4 sem keepdims shape', all(shape(kept) == [4, 3, 2, 1]))
    call check_true('rank4 sem keepdims NaN', ieee_is_nan(kept(1, 1, 1, 1)))
    call check_close('rank4 sem keepdims finite', kept(2, 1, 1, 1), 90.0_dp, tol)

    print *, 'descriptive rank-4 tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance

        if (ieee_is_nan(actual)) then
            print *, 'FAIL ', name, ' actual is NaN'
            error stop 1
        end if
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_rank3(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :, :)
        real(dp), intent(in) :: expected(:, :, :)
        real(dp), intent(in) :: tolerance

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', name, ' shape'
            print *, shape(actual), shape(expected)
            error stop 1
        end if
        if (any(ieee_is_nan(actual))) then
            print *, 'FAIL ', name, ' unexpected NaN'
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, maxval(abs(actual - expected))
            error stop 1
        end if
    end subroutine check_rank3

    subroutine check_rank3_nan(name, actual, expected, tolerance)
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
                            print *, 'FAIL ', name, ' expected NaN'
                            error stop 1
                        end if
                    else
                        if (ieee_is_nan(actual(i, j, k))) then
                            print *, 'FAIL ', name, ' unexpected NaN'
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
    end subroutine check_rank3_nan

    subroutine check_rank4(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :, :, :)
        real(dp), intent(in) :: expected(:, :, :, :)
        real(dp), intent(in) :: tolerance

        if (any(shape(actual) /= shape(expected))) then
            print *, 'FAIL ', name, ' shape'
            print *, shape(actual), shape(expected)
            error stop 1
        end if
        if (any(ieee_is_nan(actual))) then
            print *, 'FAIL ', name, ' unexpected NaN'
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, maxval(abs(actual - expected))
            error stop 1
        end if
    end subroutine check_rank4

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_rank4
