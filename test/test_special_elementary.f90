! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_elementary
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : betaln, boxcox, boxcox1p, cosm1, digamma, &
        dp, entr, erf, erfc, erfcinv, erfinv, expit, exprel, gammaln, &
        inv_boxcox, inv_boxcox1p, &
        kl_div, log_expit, log_ndtr, log_softmax, logit, logsumexp, ndtr, &
        ndtri, psi, rel_entr, softmax, xlog1py, xlogy
    implicit none

    real(dp), parameter :: reduction_x(4) = [1000.0_dp, 1001.0_dp, 999.0_dp, -1000.0_dp]
    real(dp), parameter :: softmax_expected(4) = [ &
        0.24472847105479764_dp, 0.6652409557748218_dp, &
        0.09003057317038045_dp, 0.0_dp]
    real(dp), parameter :: log_softmax_expected(4) = [ &
        -1.4076059644443804_dp, -0.4076059644443804_dp, &
        -2.4076059644443806_dp, -2001.4076059644444_dp]

    real(dp) :: empty(0)
    real(dp), allocatable :: values(:)

    call check_close('erf', erf(1.0_dp), 0.8427007929497148_dp, 2.0e-16_dp)
    call check_close('erfc tail', erfc(5.0_dp), &
        1.5374597944280347e-12_dp, 4.0e-27_dp)
    call check_close('erfinv', erfinv(0.9_dp), &
        1.1630871536766743_dp, 3.0e-14_dp)
    call check_close('erfinv near endpoint', erfinv(0.999999999999_dp), &
        5.042031898572696_dp, 2.0e-12_dp)
    call check_close('erfcinv tail', erfcinv(1.0e-200_dp), &
        21.374783049026263_dp, 3.0e-13_dp)
    call check_close('erfcinv upper', erfcinv(1.9_dp), &
        -1.1630871536766738_dp, 3.0e-14_dp)
    call check_true('erfinv invalid', ieee_is_nan(erfinv(1.1_dp)))

    call check_close('gammaln negative', gammaln(-4.75_dp), &
        -2.8754125604929133_dp, 4.0e-15_dp)
    call check_close('gammaln half', gammaln(-0.5_dp), &
        1.2655121234846454_dp, 3.0e-15_dp)
    call check_close('gammaln positive', gammaln(20.0_dp), &
        39.339884187199495_dp, 3.0e-14_dp)
    call check_true('gammaln pole', gammaln(-2.0_dp) > huge(1.0_dp))

    call check_close('betaln small', betaln(0.1_dp, 0.2_dp), &
        2.680978479346915_dp, 4.0e-15_dp)
    call check_close('betaln ordinary', betaln(2.0_dp, 3.0_dp), &
        -2.4849066497880004_dp, 4.0e-15_dp)
    call check_close('betaln large', betaln(100.0_dp, 200.0_dp), &
        -192.134192274979_dp, 2.0e-13_dp)
    call check_close('betaln mixed scale', betaln(1.0e-3_dp, 2000.0_dp), &
        6.8995782326950825_dp, 4.0e-15_dp)
    call check_true('betaln invalid shape', ieee_is_nan(betaln(-0.5_dp, 2.0_dp)))

    call check_close('digamma negative', digamma(-2.5_dp), &
        1.1031566406452433_dp, 4.0e-15_dp)
    call check_close('digamma half', digamma(0.5_dp), &
        -1.9635100260214235_dp, 4.0e-15_dp)
    call check_close('digamma small', digamma(0.1_dp), &
        -10.423754940411076_dp, 1.0e-14_dp)
    call check_close('digamma large', digamma(50.0_dp), &
        3.9019896734278925_dp, 4.0e-15_dp)
    call check_close('psi alias', psi(-0.5_dp), &
        0.03648997397857651_dp, 4.0e-15_dp)
    call check_true('digamma pole', ieee_is_nan(digamma(-3.0_dp)))

    call check_close('ndtr tail', ndtr(-10.0_dp), &
        7.61985302416047e-24_dp, 2.0e-37_dp)
    call check_close('log_ndtr far tail', log_ndtr(-30.0_dp), &
        -454.32124395634327_dp, 2.0e-12_dp)
    call check_close('log_ndtr upper', log_ndtr(10.0_dp), &
        -7.61985302416047e-24_dp, 2.0e-37_dp)
    call check_close('ndtri far tail', ndtri(1.0e-200_dp), &
        -30.20559417957964_dp, 3.0e-13_dp)
    call check_close('ndtri ordinary', ndtri(0.001_dp), &
        -3.090232306167813_dp, 3.0e-14_dp)

    call check_close('expit negative', expit(-30.0_dp), &
        9.357622968839299e-14_dp, 3.0e-28_dp)
    call check_close('expit positive', expit(2.0_dp), &
        0.8807970779778823_dp, 2.0e-16_dp)
    call check_close('log_expit negative', log_expit(-30.0_dp), &
        -30.000000000000092_dp, 5.0e-14_dp)
    call check_close('log_expit positive', log_expit(30.0_dp), &
        -9.357622968839737e-14_dp, 3.0e-28_dp)
    call check_close('logit tiny', logit(1.0e-200_dp), &
        -460.51701859880916_dp, 6.0e-14_dp)
    call check_close('logit ordinary', logit(0.9_dp), &
        2.1972245773362196_dp, 4.0e-16_dp)

    call check_close('xlogy ordinary', xlogy(2.0_dp, 3.0_dp), &
        2.1972245773362196_dp, 4.0e-16_dp)
    call check_close('xlogy zero convention', xlogy(0.0_dp, 0.0_dp), 0.0_dp, 0.0_dp)
    call check_close('xlog1py stable', xlog1py(1.0_dp, 1.0e-12_dp), &
        9.999999999995e-13_dp, 2.0e-28_dp)
    call check_close('xlog1py zero convention', xlog1py(0.0_dp, -1.0_dp), &
        0.0_dp, 0.0_dp)

    call check_close('entr', entr(0.1_dp), 0.23025850929940456_dp, 4.0e-16_dp)
    call check_close('rel_entr', rel_entr(2.0_dp, 3.0_dp), &
        -0.8109302162163288_dp, 4.0e-16_dp)
    call check_close('kl_div', kl_div(2.0_dp, 3.0_dp), &
        0.18906978378367123_dp, 8.0e-16_dp)
    call check_true('entr negative', entr(-1.0_dp) < -huge(1.0_dp))

    call check_close('exprel tiny positive', exprel(1.0e-10_dp), &
        1.00000000005_dp, 3.0e-16_dp)
    call check_close('exprel negative', exprel(-1.0_dp), &
        0.6321205588285577_dp, 2.0e-16_dp)
    call check_close('cosm1 tiny', cosm1(1.0e-10_dp), &
        -5.0000000000000005e-21_dp, 2.0e-36_dp)
    call check_close('cosm1 ordinary', cosm1(1.0_dp), &
        -0.45969769413186023_dp, 2.0e-16_dp)

    call check_close('boxcox tiny lambda', boxcox(0.25_dp, 1.0e-10_dp), &
        -1.3862943610237999_dp, 5.0e-16_dp)
    call check_close('boxcox negative lambda', boxcox(0.25_dp, -2.0_dp), &
        -7.5_dp, 2.0e-15_dp)
    call check_close('boxcox inverse', inv_boxcox(boxcox(10.0_dp, 2.0_dp), 2.0_dp), &
        10.0_dp, 4.0e-15_dp)
    call check_close('boxcox1p tiny lambda', boxcox1p(-0.75_dp, 1.0e-10_dp), &
        -1.3862943610237999_dp, 5.0e-16_dp)
    call check_close('boxcox1p inverse', &
        inv_boxcox1p(boxcox1p(9.0_dp, 2.0_dp), 2.0_dp), 9.0_dp, 4.0e-15_dp)
    call check_true('boxcox invalid domain', ieee_is_nan(boxcox(-1.0_dp, 2.0_dp)))

    call check_close('logsumexp', logsumexp(reduction_x), &
        1001.4076059644444_dp, 3.0e-13_dp)
    values = softmax(reduction_x)
    call check_array('softmax', values, softmax_expected, 4.0e-16_dp)
    values = log_softmax(reduction_x)
    call check_array('log_softmax', values, log_softmax_expected, 4.0e-13_dp)
    call check_true('logsumexp empty', logsumexp(empty) < -huge(1.0_dp))

    print *, 'special elementary tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_array(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: expected(:) !! reference values
        real(dp), intent(in) :: tolerance !! maximum elementwise absolute error

        if (size(actual) /= size(expected) .or. &
                any(ieee_is_nan(actual)) .or. &
                any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_array

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_elementary
