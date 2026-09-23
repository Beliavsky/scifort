program test_goodness_of_fit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use goodness_of_fit_reference
    use scifort_kinds, only : dp
    use scifort_stats, only : anderson, anderson_result, cramervonmises, &
        cramervonmises_2samp, jarque_bera, ks_1samp, ks_2samp, kstest, &
        kstest_result, kurtosistest, normal_cdf, normaltest, shapiro, &
        shapiro_result, significance_result, skewtest
    implicit none

    real(dp), parameter :: a(8) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
        5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
    real(dp), parameter :: b(12) = [1.1_dp, 2.0_dp, 2.5_dp, 3.7_dp, &
        4.2_dp, 5.9_dp, 7.1_dp, 8.4_dp, 9.0_dp, 10.3_dp, 11.5_dp, 13.0_dp]
    real(dp), parameter :: x(5) = [-1.2_dp, -0.4_dp, 0.1_dp, 0.7_dp, 1.1_dp]
    real(dp), parameter :: x2(4) = [0.1_dp, 0.4_dp, 0.8_dp, 1.2_dp]
    real(dp), parameter :: y2(5) = [0.2_dp, 0.3_dp, 0.9_dp, 1.5_dp, 2.0_dp]
    real(dp), parameter :: cvm_x(7) = [-1.0_dp, -0.2_dp, 0.0_dp, 0.3_dp, &
        0.7_dp, 1.5_dp, 2.2_dp]
    real(dp), parameter :: cvm_y(6) = [-0.8_dp, 0.1_dp, 0.4_dp, 0.9_dp, &
        1.8_dp, 2.5_dp]
    real(dp), parameter :: shap_x(10) = [0.1_dp, 0.2_dp, 0.3_dp, 0.8_dp, 1.2_dp, &
        1.5_dp, 2.1_dp, 2.2_dp, 3.0_dp, 3.5_dp]
    real(dp), parameter :: and_x(10) = [0.2_dp, 0.4_dp, 0.6_dp, 0.9_dp, 1.1_dp, &
        1.3_dp, 1.8_dp, 2.0_dp, 2.5_dp, 3.1_dp]

    type(anderson_result) :: adr
    type(kstest_result) :: ksr
    type(shapiro_result) :: sr
    type(significance_result) :: rr
    integer :: failures

    failures = 0

    rr = skewtest(a)
    call check_close('skewtest statistic', rr%statistic, skew_ref(1), 2.0e-15_dp, failures)
    call check_close('skewtest pvalue', rr%pvalue, skew_ref(2), 2.0e-15_dp, failures)
    rr = kurtosistest(b)
    call check_close('kurtosistest statistic', rr%statistic, kurt_ref(1), 3.0e-15_dp, failures)
    call check_close('kurtosistest pvalue', rr%pvalue, kurt_ref(2), 3.0e-15_dp, failures)
    rr = normaltest(b)
    call check_close('normaltest statistic', rr%statistic, normal_ref(1), 3.0e-15_dp, failures)
    call check_close('normaltest pvalue', rr%pvalue, normal_ref(2), 3.0e-15_dp, failures)
    rr = jarque_bera(b)
    call check_close('jarque bera statistic', rr%statistic, jb_ref(1), 2.0e-15_dp, failures)
    call check_close('jarque bera pvalue', rr%pvalue, jb_ref(2), 2.0e-15_dp, failures)

    ksr = ks_1samp(x, standard_normal_cdf, method='exact')
    call check_close('ks1 statistic', ksr%statistic, ks1_ref(1), 2.0e-15_dp, failures)
    call check_close('ks1 exact pvalue', ksr%pvalue, ks1_ref(2), 3.0e-15_dp, failures)
    call check_close('ks1 location', ksr%statistic_location, ks1_ref(3), 0.0_dp, failures)
    call check_int('ks1 sign', ksr%statistic_sign, nint(ks1_ref(4)), failures)
    ksr = ks_1samp(x, standard_normal_cdf, method='asymp')
    call check_close('ks1 asymp pvalue', ksr%pvalue, ks1_methods_ref(1), 3.0e-15_dp, failures)
    ksr = kstest(x, standard_normal_cdf, method='approx')
    call check_close('kstest alias approx', ksr%pvalue, ks1_methods_ref(2), 0.0_dp, failures)

    ksr = ks_2samp(x2, y2, method='exact')
    call check_close('ks2 statistic', ksr%statistic, ks2_ref(1), 0.0_dp, failures)
    call check_close('ks2 pvalue', ksr%pvalue, ks2_ref(2), 3.0e-15_dp, failures)
    call check_close('ks2 location', ksr%statistic_location, ks2_ref(3), 0.0_dp, failures)
    call check_int('ks2 sign', ksr%statistic_sign, nint(ks2_ref(4)), failures)
    ksr = ks_2samp(x2, y2, alternative='less', method='exact')
    call check_close('ks2 less statistic', ksr%statistic, ks2_one_sided_ref(1), 3.0e-16_dp, failures)
    call check_close('ks2 less pvalue', ksr%pvalue, ks2_one_sided_ref(2), 3.0e-15_dp, failures)
    ksr = kstest(x2, y2, alternative='greater', method='exact')
    call check_close('ks2 greater statistic', ksr%statistic, ks2_one_sided_ref(3), 0.0_dp, failures)
    call check_close('ks2 greater pvalue', ksr%pvalue, ks2_one_sided_ref(4), 3.0e-15_dp, failures)

    rr = cramervonmises(x, standard_normal_cdf)
    call check_close('cvm1 statistic', rr%statistic, cvm_ref(1), 3.0e-15_dp, failures)
    call check_close('cvm1 pvalue', rr%pvalue, cvm_ref(2), 3.0e-15_dp, failures)
    rr = cramervonmises_2samp(cvm_x, cvm_y, method='exact')
    call check_close('cvm2 exact statistic', rr%statistic, cvm_ref(3), 3.0e-15_dp, failures)
    call check_close('cvm2 exact pvalue', rr%pvalue, cvm_ref(4), 3.0e-15_dp, failures)
    rr = cramervonmises_2samp(cvm_x, cvm_y, method='asymptotic')
    call check_close('cvm2 asymptotic pvalue', rr%pvalue, cvm_ref(5), 3.0e-15_dp, failures)

    sr = shapiro(shap_x)
    call check_close('shapiro statistic', sr%statistic, shapiro_ref(1), 2.0e-15_dp, failures)
    call check_close('shapiro pvalue', sr%pvalue, shapiro_ref(2), 5.0e-12_dp, failures)
    call check_int('shapiro status', sr%status, 0, failures)

    adr = anderson(and_x, 'norm', .true.)
    call check_close('anderson norm statistic', adr%statistic, and_norm_stat, 3.0e-14_dp, failures)
    call check_array('anderson norm critical', adr%critical_values, and_norm_critical, 0.0_dp, failures)
    call check_array('anderson norm fit', adr%fit_params, and_norm_fit, 3.0e-15_dp, failures)
    adr = anderson(and_x + 0.1_dp, 'expon', .true.)
    call check_close('anderson expon statistic', adr%statistic, and_expon_stat, 3.0e-14_dp, failures)
    call check_array('anderson expon fit', adr%fit_params, and_expon_fit, 3.0e-15_dp, failures)
    adr = anderson(and_x, 'logistic', .true.)
    call check_close('anderson logistic statistic', adr%statistic, and_logistic_stat, 1.0e-10_dp, failures)
    call check_array('anderson logistic fit', adr%fit_params, and_logistic_fit, 1.0e-9_dp, failures)
    adr = anderson(and_x, 'gumbel_r', .true.)
    call check_close('anderson gumbel r statistic', adr%statistic, and_gumbel_r_stat, 8.0e-14_dp, failures)
    call check_array('anderson gumbel r fit', adr%fit_params, and_gumbel_r_fit, 8.0e-14_dp, failures)
    adr = anderson(and_x, 'gumbel_l', .true.)
    call check_close('anderson gumbel l statistic', adr%statistic, and_gumbel_l_stat, 8.0e-14_dp, failures)
    call check_array('anderson gumbel l fit', adr%fit_params, and_gumbel_l_fit, 8.0e-14_dp, failures)

    rr = skewtest([1.0_dp, 2.0_dp, 3.0_dp])
    call check_true('skewtest short sample invalid', ieee_is_nan(rr%statistic), failures)
    sr = shapiro([1.0_dp, 1.0_dp, 1.0_dp])
    call check_int('shapiro constant status', sr%status, 6, failures)
    adr = anderson(and_x, 'unsupported')
    call check_true('anderson unsupported invalid', ieee_is_nan(adr%statistic), failures)

    if (failures /= 0) then
        print *, 'goodness-of-fit failures:', failures
        error stop 1
    end if
    print *, 'goodness-of-fit tests passed'

contains

    function standard_normal_cdf(value) result(p)
        real(dp), intent(in) :: value !! point at which the standard normal CDF is evaluated
        real(dp) :: p

        p = normal_cdf(value, 0.0_dp, 1.0_dp)
    end function standard_normal_cdf

    subroutine check_close(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected_value !! SciPy 1.17 reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        integer, intent(inout) :: count !! accumulated failure count

        if (ieee_is_nan(actual) .or. abs(actual - expected_value) > tolerance) then
            print *, 'FAIL ', trim(name), actual, expected_value
            count = count + 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected_value, tolerance, count)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! computed array
        real(dp), intent(in) :: expected_value(:) !! SciPy 1.17 reference array
        real(dp), intent(in) :: tolerance !! maximum absolute element error
        integer, intent(inout) :: count !! accumulated failure count

        if (size(actual) /= size(expected_value)) then
            print *, 'FAIL ', trim(name), ' size'
            count = count + 1
        else if (any(abs(actual - expected_value) > tolerance)) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_array

    subroutine check_int(name, actual, expected_value, count)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! computed integer
        integer, intent(in) :: expected_value !! expected integer
        integer, intent(inout) :: count !! accumulated failure count

        if (actual /= expected_value) then
            print *, 'FAIL ', trim(name), actual, expected_value
            count = count + 1
        end if
    end subroutine check_int

    subroutine check_true(name, condition, count)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition expected to hold
        integer, intent(inout) :: count !! accumulated failure count

        if (.not. condition) then
            print *, 'FAIL ', trim(name)
            count = count + 1
        end if
    end subroutine check_true

end program test_goodness_of_fit
