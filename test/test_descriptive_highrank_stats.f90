program test_descriptive_highrank_stats
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_stats, only : sem, variation, iqr, median_abs_deviation, moment, &
        skew, kurtosis, sem_keepdims, variation_keepdims, iqr_keepdims, &
        median_abs_deviation_keepdims, moment_keepdims, skew_keepdims, kurtosis_keepdims
    implicit none

    integer :: i
    real(dp), parameter :: tol = 5.0e-13_dp
    real(dp) :: x5(2, 1, 1, 1, 1)
    real(dp) :: x6(2, 1, 1, 1, 1, 1)
    real(dp) :: x7(2, 1, 1, 1, 1, 1, 1)
    real(dp) :: x8(2, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x9(2, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x10(2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x11(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x12(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x13(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x14(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: y15(2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: z15(2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: ysq15(2, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: flat(:)

    x5 = reshape([1.0_dp, 3.0_dp], shape(x5))
    x6 = reshape([1.0_dp, 3.0_dp], shape(x6))
    x7 = reshape([1.0_dp, 3.0_dp], shape(x7))
    x8 = reshape([1.0_dp, 3.0_dp], shape(x8))
    x9 = reshape([1.0_dp, 3.0_dp], shape(x9))
    x10 = reshape([1.0_dp, 3.0_dp], shape(x10))
    x11 = reshape([1.0_dp, 3.0_dp], shape(x11))
    x12 = reshape([1.0_dp, 3.0_dp], shape(x12))
    x13 = reshape([1.0_dp, 3.0_dp], shape(x13))
    x14 = reshape([1.0_dp, 3.0_dp], shape(x14))
    x15 = reshape([1.0_dp, 3.0_dp], shape(x15))
    y15 = reshape([(real(i, dp), i = 1, size(y15))], shape(y15))
    ysq15 = y15 ** 2
    z15 = y15
    z15(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = &
        ieee_value(0.0_dp, ieee_quiet_nan)

    call check_close(sum(sem(x5, axis=1)), 1.0_dp, 'sem rank 5')
    call check_close(sum(variation(x5, axis=1)), 0.5_dp, 'variation rank 5')
    call check_close(sum(iqr(x5, axis=1)), 1.0_dp, 'iqr rank 5')
    call check_close(sum(median_abs_deviation(x5, axis=1)), 1.0_dp, 'mad rank 5')
    call check_close(sum(moment(x5, order=2, axis=1)), 1.0_dp, 'moment rank 5')
    call check_close(sum(skew(x5, axis=1)), 0.0_dp, 'skew rank 5')
    call check_close(sum(kurtosis(x5, axis=1)), -2.0_dp, 'kurtosis rank 5')

    call check_close(sum(sem(x6, axis=1)), 1.0_dp, 'sem rank 6')
    call check_close(sum(variation(x6, axis=1)), 0.5_dp, 'variation rank 6')
    call check_close(sum(iqr(x6, axis=1)), 1.0_dp, 'iqr rank 6')
    call check_close(sum(median_abs_deviation(x6, axis=1)), 1.0_dp, 'mad rank 6')
    call check_close(sum(moment(x6, order=2, axis=1)), 1.0_dp, 'moment rank 6')
    call check_close(sum(skew(x6, axis=1)), 0.0_dp, 'skew rank 6')
    call check_close(sum(kurtosis(x6, axis=1)), -2.0_dp, 'kurtosis rank 6')

    call check_close(sum(sem(x7, axis=1)), 1.0_dp, 'sem rank 7')
    call check_close(sum(variation(x7, axis=1)), 0.5_dp, 'variation rank 7')
    call check_close(sum(iqr(x7, axis=1)), 1.0_dp, 'iqr rank 7')
    call check_close(sum(median_abs_deviation(x7, axis=1)), 1.0_dp, 'mad rank 7')
    call check_close(sum(moment(x7, order=2, axis=1)), 1.0_dp, 'moment rank 7')
    call check_close(sum(skew(x7, axis=1)), 0.0_dp, 'skew rank 7')
    call check_close(sum(kurtosis(x7, axis=1)), -2.0_dp, 'kurtosis rank 7')

    call check_close(sum(sem(x8, axis=1)), 1.0_dp, 'sem rank 8')
    call check_close(sum(variation(x8, axis=1)), 0.5_dp, 'variation rank 8')
    call check_close(sum(iqr(x8, axis=1)), 1.0_dp, 'iqr rank 8')
    call check_close(sum(median_abs_deviation(x8, axis=1)), 1.0_dp, 'mad rank 8')
    call check_close(sum(moment(x8, order=2, axis=1)), 1.0_dp, 'moment rank 8')
    call check_close(sum(skew(x8, axis=1)), 0.0_dp, 'skew rank 8')
    call check_close(sum(kurtosis(x8, axis=1)), -2.0_dp, 'kurtosis rank 8')

    call check_close(sum(sem(x9, axis=1)), 1.0_dp, 'sem rank 9')
    call check_close(sum(variation(x9, axis=1)), 0.5_dp, 'variation rank 9')
    call check_close(sum(iqr(x9, axis=1)), 1.0_dp, 'iqr rank 9')
    call check_close(sum(median_abs_deviation(x9, axis=1)), 1.0_dp, 'mad rank 9')
    call check_close(sum(moment(x9, order=2, axis=1)), 1.0_dp, 'moment rank 9')
    call check_close(sum(skew(x9, axis=1)), 0.0_dp, 'skew rank 9')
    call check_close(sum(kurtosis(x9, axis=1)), -2.0_dp, 'kurtosis rank 9')

    call check_close(sum(sem(x10, axis=1)), 1.0_dp, 'sem rank 10')
    call check_close(sum(variation(x10, axis=1)), 0.5_dp, 'variation rank 10')
    call check_close(sum(iqr(x10, axis=1)), 1.0_dp, 'iqr rank 10')
    call check_close(sum(median_abs_deviation(x10, axis=1)), 1.0_dp, 'mad rank 10')
    call check_close(sum(moment(x10, order=2, axis=1)), 1.0_dp, 'moment rank 10')
    call check_close(sum(skew(x10, axis=1)), 0.0_dp, 'skew rank 10')
    call check_close(sum(kurtosis(x10, axis=1)), -2.0_dp, 'kurtosis rank 10')

    call check_close(sum(sem(x11, axis=1)), 1.0_dp, 'sem rank 11')
    call check_close(sum(variation(x11, axis=1)), 0.5_dp, 'variation rank 11')
    call check_close(sum(iqr(x11, axis=1)), 1.0_dp, 'iqr rank 11')
    call check_close(sum(median_abs_deviation(x11, axis=1)), 1.0_dp, 'mad rank 11')
    call check_close(sum(moment(x11, order=2, axis=1)), 1.0_dp, 'moment rank 11')
    call check_close(sum(skew(x11, axis=1)), 0.0_dp, 'skew rank 11')
    call check_close(sum(kurtosis(x11, axis=1)), -2.0_dp, 'kurtosis rank 11')

    call check_close(sum(sem(x12, axis=1)), 1.0_dp, 'sem rank 12')
    call check_close(sum(variation(x12, axis=1)), 0.5_dp, 'variation rank 12')
    call check_close(sum(iqr(x12, axis=1)), 1.0_dp, 'iqr rank 12')
    call check_close(sum(median_abs_deviation(x12, axis=1)), 1.0_dp, 'mad rank 12')
    call check_close(sum(moment(x12, order=2, axis=1)), 1.0_dp, 'moment rank 12')
    call check_close(sum(skew(x12, axis=1)), 0.0_dp, 'skew rank 12')
    call check_close(sum(kurtosis(x12, axis=1)), -2.0_dp, 'kurtosis rank 12')

    call check_close(sum(sem(x13, axis=1)), 1.0_dp, 'sem rank 13')
    call check_close(sum(variation(x13, axis=1)), 0.5_dp, 'variation rank 13')
    call check_close(sum(iqr(x13, axis=1)), 1.0_dp, 'iqr rank 13')
    call check_close(sum(median_abs_deviation(x13, axis=1)), 1.0_dp, 'mad rank 13')
    call check_close(sum(moment(x13, order=2, axis=1)), 1.0_dp, 'moment rank 13')
    call check_close(sum(skew(x13, axis=1)), 0.0_dp, 'skew rank 13')
    call check_close(sum(kurtosis(x13, axis=1)), -2.0_dp, 'kurtosis rank 13')

    call check_close(sum(sem(x14, axis=1)), 1.0_dp, 'sem rank 14')
    call check_close(sum(variation(x14, axis=1)), 0.5_dp, 'variation rank 14')
    call check_close(sum(iqr(x14, axis=1)), 1.0_dp, 'iqr rank 14')
    call check_close(sum(median_abs_deviation(x14, axis=1)), 1.0_dp, 'mad rank 14')
    call check_close(sum(moment(x14, order=2, axis=1)), 1.0_dp, 'moment rank 14')
    call check_close(sum(skew(x14, axis=1)), 0.0_dp, 'skew rank 14')
    call check_close(sum(kurtosis(x14, axis=1)), -2.0_dp, 'kurtosis rank 14')

    call check_close(sum(sem(x15, axis=1)), 1.0_dp, 'sem rank 15')
    call check_close(sum(variation(x15, axis=1)), 0.5_dp, 'variation rank 15')
    call check_close(sum(iqr(x15, axis=1)), 1.0_dp, 'iqr rank 15')
    call check_close(sum(median_abs_deviation(x15, axis=1)), 1.0_dp, 'mad rank 15')
    call check_close(sum(moment(x15, order=2, axis=1)), 1.0_dp, 'moment rank 15')
    call check_close(sum(skew(x15, axis=1)), 0.0_dp, 'skew rank 15')
    call check_close(sum(kurtosis(x15, axis=1)), -2.0_dp, 'kurtosis rank 15')

    flat = reshape(sem(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp, i = 1, 24)], 'rank15 sem axis 5')

    flat = reshape(variation(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp / real(i + 12, dp), i = 1, 24)], &
        'rank15 variation axis 5')

    flat = reshape(iqr(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp, i = 1, 24)], 'rank15 iqr axis 5')

    flat = reshape(median_abs_deviation(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp, i = 1, 24)], 'rank15 mad axis 5')

    flat = reshape(moment(y15, order=2, axis=5), [24])
    call check_all_close(flat, [(144.0_dp, i = 1, 24)], 'rank15 moment axis 5')

    flat = reshape(skew(y15, axis=5, bias=.false.), [24])
    call check_all_close(flat, [(0.0_dp, i = 1, 24)], 'rank15 skew axis 5')

    flat = reshape(kurtosis(y15, axis=5, fisher=.true., bias=.false.), [24])
    call check_all_close(flat, [(-2.0_dp, i = 1, 24)], 'rank15 kurtosis axis 5')

    call check_close(iqr(y15), 23.5_dp, 'rank15 iqr all axes')

    flat = reshape(moment(ysq15, order=3, axis=2), [16])
    call check_all_close(flat(1:6), [379.2592592592594_dp, 677.9259259259253_dp, &
        3451.259259259248_dp, 4261.92592592591_dp, 9595.259259259328_dp, &
        10917.925925925761_dp], 'rank15 nontrivial moment axis 2')
    flat = reshape(skew(ysq15, axis=2, bias=.false.), [16])
    call check_all_close(flat(1:6), [0.9352195295828247_dp, 0.7221086457211349_dp, &
        0.3308318147058979_dp, 0.2981745522010502_dp, 0.19945783082154267_dp, &
        0.18705316268389463_dp], 'rank15 nontrivial skew axis 2')
    flat = reshape(sem(ysq15, axis=2), [16])
    call check_all_close(flat(1:6), [7.055336829505576_dp, 9.333333333333334_dp, &
        20.827332469084414_dp, 23.132468763863297_dp, 34.66666666666667_dp, &
        36.97446566363213_dp], 'rank15 nontrivial sem axis 2')

    flat = reshape(sem_keepdims(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp, i = 1, 24)], 'rank15 sem keepdims')
    flat = reshape(variation_keepdims(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp / real(i + 12, dp), i = 1, 24)], &
        'rank15 variation keepdims')
    flat = reshape(iqr_keepdims(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp, i = 1, 24)], 'rank15 iqr keepdims')
    flat = reshape(median_abs_deviation_keepdims(y15, axis=5), [24])
    call check_all_close(flat, [(12.0_dp, i = 1, 24)], 'rank15 mad keepdims')
    flat = reshape(moment_keepdims(y15, order=2, axis=5), [24])
    call check_all_close(flat, [(144.0_dp, i = 1, 24)], 'rank15 moment keepdims')
    flat = reshape(skew_keepdims(y15, axis=5), [24])
    call check_all_close(flat, [(0.0_dp, i = 1, 24)], 'rank15 skew keepdims')
    flat = reshape(kurtosis_keepdims(y15, axis=5), [24])
    call check_all_close(flat, [(-2.0_dp, i = 1, 24)], 'rank15 kurtosis keepdims')

    flat = reshape(sem(z15, axis=5, nan_policy='propagate'), [24])
    call check_nan(flat(1), 'rank15 sem propagate nan')
    call check_close(flat(2), 12.0_dp, 'rank15 sem propagate finite')
    flat = reshape(sem(z15, axis=5, nan_policy='omit'), [24])
    call check_nan(flat(1), 'rank15 sem omit small sample')
    call check_close(flat(2), 12.0_dp, 'rank15 sem omit finite')
    flat = reshape(variation(z15, axis=5, nan_policy='omit'), [24])
    call check_close(flat(1), 0.0_dp, 'rank15 variation omit singleton')
    flat = reshape(iqr(z15, axis=5, nan_policy='omit'), [24])
    call check_close(flat(1), 0.0_dp, 'rank15 iqr omit singleton')
    flat = reshape(median_abs_deviation(z15, axis=5, nan_policy='omit'), [24])
    call check_close(flat(1), 0.0_dp, 'rank15 mad omit singleton')
    flat = reshape(moment(z15, order=2, axis=5, nan_policy='omit'), [24])
    call check_close(flat(1), 0.0_dp, 'rank15 moment omit singleton')
    flat = reshape(skew(z15, axis=5, nan_policy='omit'), [24])
    call check_nan(flat(1), 'rank15 skew omit singleton')
    flat = reshape(kurtosis(z15, axis=5, nan_policy='omit'), [24])
    call check_nan(flat(1), 'rank15 kurtosis omit singleton')

    if (size(sem(y15, axis=16)) /= 0) error stop 'invalid rank15 sem axis'
    if (size(iqr(y15, axis=16)) /= 0) error stop 'invalid rank15 iqr axis'

    print *, 'high-rank descriptive-statistics tests passed'

contains

    subroutine check_close(actual, expected, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        character(len=*), intent(in) :: label

        if (ieee_is_nan(actual) .or. abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_all_close(actual, expected, label)
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        character(len=*), intent(in) :: label
        integer :: j

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), ' size mismatch'
            error stop 1
        end if
        do j = 1, size(actual)
            call check_close(actual(j), expected(j), label)
        end do
    end subroutine check_all_close

    subroutine check_nan(value, label)
        real(dp), intent(in) :: value
        character(len=*), intent(in) :: label

        if (.not. ieee_is_nan(value)) then
            print *, 'FAIL ', trim(label), ': expected NaN, got ', value
            error stop 1
        end if
    end subroutine check_nan

end program test_descriptive_highrank_stats
