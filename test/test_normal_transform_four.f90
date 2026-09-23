! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_normal_transform_four
    use scifort_c_api
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state, rng_uniform
    use scifort_stats
    use test_checks, only : check_close, check_true
    implicit none

    integer :: failures
    real(dp), parameter :: probs(5) = [1.0e-6_dp, 0.1_dp, 0.5_dp, 0.9_dp, 1.0_dp - 1.0e-6_dp]

    failures = 0
    call test_reference(failures)
    call test_scores(failures)
    call test_rng_fit_c_api(failures)

    if (failures == 0) then
        print '(a)', 'test_normal_transform_four: PASS'
    else
        print '(a,1x,i0)', 'test_normal_transform_four: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: en_x(5) = [-2.0_dp,-0.5_dp,0.8_dp,2.5_dp,7.0_dp]
        real(dp), parameter :: en_pdf(5) = [0.012774486411762549_dp,0.14798612844894818_dp, &
            0.27807308447076806_dp,0.14884319453077810_dp,0.0062113923577650221_dp]
        real(dp), parameter :: en_cdf(5) = [0.0048658509717116258_dp,0.10135695889745944_dp, &
            0.39884228315752812_dp,0.78540986233113452_dp,0.99130405069784910_dp]
        real(dp), parameter :: en_ppf(5) = [-4.3259726227007391_dp,-0.50921126987269971_dp, &
            1.1666445470548701_dp,3.5802536942183565_dp,19.698857638206210_dp]
        real(dp), parameter :: sb_x(5) = [0.05_dp,0.2_dp,0.5_dp,0.8_dp,0.95_dp]
        real(dp), parameter :: sb_pdf(5) = [0.00035378792304169429_dp,0.62496034873649686_dp, &
            2.5042369540625988_dp,0.094852957298199855_dp,6.4511643538327926e-6_dp]
        real(dp), parameter :: sb_cdf(5) = [2.0569225434159252e-6_dp,0.025191350607344912_dp, &
            0.65542174161032418_dp,0.99708060902183804_dp,0.99999996769441490_dp]
        real(dp), parameter :: sb_ppf(5) = [0.046026171204382349_dp,0.27108040045101883_dp, &
            0.44144636529227166_dp,0.62681092476672673_dp,0.92829824414522355_dp]
        real(dp), parameter :: su_x(5) = [-4.0_dp,-1.0_dp,0.0_dp,1.5_dp,5.0_dp]
        real(dp), parameter :: su_pdf(5) = [0.00050307807478167674_dp,0.079895827443765097_dp, &
            0.43319198375933959_dp,0.18264892062731375_dp,0.0056254789041269562_dp]
        real(dp), parameter :: su_cdf(5) = [0.00044507299094451457_dp,0.040424101567469242_dp, &
            0.27425311775007361_dp,0.82975364671380358_dp,0.99193961599950053_dp]
        real(dp), parameter :: su_ppf(5) = [-12.183891782733545_dp,-0.54861944934885942_dp, &
            0.47809985545162414_dp,2.0083138197118728_dp,30.710559448690518_dp]
        real(dp), parameter :: tr_x(5) = [0.05_dp,0.2_dp,0.4_dp,0.8_dp,0.95_dp]
        real(dp), parameter :: tr_pdf(5) = [0.27586206896551729_dp,1.1034482758620692_dp, &
            1.3793103448275863_dp,0.91954022988505713_dp,0.22988505747126456_dp]
        real(dp), parameter :: tr_cdf(5) = [0.0068965517241379327_dp,0.11034482758620692_dp, &
            0.37931034482758624_dp,0.90804597701149437_dp,0.99425287356321834_dp]
        real(dp), parameter :: tr_ppf(5) = [0.00060207972893961477_dp,0.19039432764659769_dp, &
            0.48749999999999999_dp,0.79143346385385793_dp,0.99934045470207689_dp]
        integer :: i

        do i = 1, 5
            call check_close('exponnorm pdf',exponnorm_pdf(en_x(i),1.4_dp),en_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('exponnorm cdf',exponnorm_cdf(en_x(i),1.4_dp),en_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('exponnorm ppf',exponnorm_ppf(probs(i),1.4_dp),en_ppf(i), &
                2.0e-9_dp,2.0e-9_dp,failures)
            call check_close('johnsonsb pdf',johnsonsb_pdf(sb_x(i),0.4_dp,1.7_dp),sb_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('johnsonsb cdf',johnsonsb_cdf(sb_x(i),0.4_dp,1.7_dp),sb_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('johnsonsb ppf',johnsonsb_ppf(probs(i),0.4_dp,1.7_dp),sb_ppf(i), &
                2.0e-10_dp,2.0e-9_dp,failures)
            call check_close('johnsonsu pdf',johnsonsu_pdf(su_x(i),-0.6_dp,1.3_dp),su_pdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('johnsonsu cdf',johnsonsu_cdf(su_x(i),-0.6_dp,1.3_dp),su_cdf(i), &
                8.0e-12_dp,8.0e-11_dp,failures)
            call check_close('johnsonsu ppf',johnsonsu_ppf(probs(i),-0.6_dp,1.3_dp),su_ppf(i), &
                3.0e-9_dp,3.0e-9_dp,failures)
            call check_close('trapezoid pdf',trapezoid_pdf(tr_x(i),0.25_dp,0.7_dp),tr_pdf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('trapezoid cdf',trapezoid_cdf(tr_x(i),0.25_dp,0.7_dp),tr_cdf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
            call check_close('trapezoid ppf',trapezoid_ppf(probs(i),0.25_dp,0.7_dp),tr_ppf(i), &
                5.0e-12_dp,5.0e-11_dp,failures)
        end do
        call check_close('exponnorm logcdf',exponnorm_logcdf(-2.0_dp,1.4_dp), &
            -5.3255136615362444_dp,8.0e-12_dp,8.0e-11_dp,failures)
        call check_close('johnsonsb logsf',johnsonsb_logsf(0.8_dp,0.4_dp,1.7_dp), &
            -5.8363802535781328_dp,8.0e-12_dp,8.0e-11_dp,failures)
        call check_close('johnsonsu logcdf',johnsonsu_logcdf(-4.0_dp,-0.6_dp,1.3_dp), &
            -7.7172722646541825_dp,8.0e-12_dp,8.0e-11_dp,failures)
        call check_close('trapezoid logsf',trapezoid_logsf(0.8_dp,0.25_dp,0.7_dp), &
            -2.3864665769747488_dp,5.0e-12_dp,5.0e-11_dp,failures)
    end subroutine test_reference

    subroutine test_scores(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 2.0e-4_dp
        real(dp), parameter :: en_data(4) = [-0.4_dp,0.3_dp,1.2_dp,2.8_dp]
        real(dp), parameter :: sb_data(4) = [0.28_dp,0.52_dp,0.82_dp,1.12_dp]
        real(dp), parameter :: su_data(4) = [-2.0_dp,-0.3_dp,1.2_dp,3.0_dp]
        real(dp), parameter :: tr_data(4) = [0.22_dp,0.58_dp,0.76_dp,1.12_dp]
        real(dp) :: q3(3), q4(4)

        q3 = exponnorm_score(en_data,1.4_dp,0.1_dp,1.2_dp)
        call check_close('exponnorm score K',q3(1),centered( &
            exponnorm_loglikelihood(en_data,1.4_dp+h,0.1_dp,1.2_dp), &
            exponnorm_loglikelihood(en_data,1.4_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('exponnorm score loc',q3(2),centered( &
            exponnorm_loglikelihood(en_data,1.4_dp,0.1_dp+h,1.2_dp), &
            exponnorm_loglikelihood(en_data,1.4_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('exponnorm score scale',q3(3),centered( &
            exponnorm_loglikelihood(en_data,1.4_dp,0.1_dp,1.2_dp+h), &
            exponnorm_loglikelihood(en_data,1.4_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q4 = johnsonsb_score(sb_data,0.4_dp,1.7_dp,0.1_dp,1.2_dp)
        call check_close('johnsonsb score a',q4(1),centered( &
            johnsonsb_loglikelihood(sb_data,0.4_dp+h,1.7_dp,0.1_dp,1.2_dp), &
            johnsonsb_loglikelihood(sb_data,0.4_dp-h,1.7_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('johnsonsb score b',q4(2),centered( &
            johnsonsb_loglikelihood(sb_data,0.4_dp,1.7_dp+h,0.1_dp,1.2_dp), &
            johnsonsb_loglikelihood(sb_data,0.4_dp,1.7_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('johnsonsb score loc',q4(3),centered( &
            johnsonsb_loglikelihood(sb_data,0.4_dp,1.7_dp,0.1_dp+h,1.2_dp), &
            johnsonsb_loglikelihood(sb_data,0.4_dp,1.7_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('johnsonsb score scale',q4(4),centered( &
            johnsonsb_loglikelihood(sb_data,0.4_dp,1.7_dp,0.1_dp,1.2_dp+h), &
            johnsonsb_loglikelihood(sb_data,0.4_dp,1.7_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q4 = johnsonsu_score(su_data,-0.6_dp,1.3_dp,0.1_dp,1.2_dp)
        call check_close('johnsonsu score a',q4(1),centered( &
            johnsonsu_loglikelihood(su_data,-0.6_dp+h,1.3_dp,0.1_dp,1.2_dp), &
            johnsonsu_loglikelihood(su_data,-0.6_dp-h,1.3_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('johnsonsu score b',q4(2),centered( &
            johnsonsu_loglikelihood(su_data,-0.6_dp,1.3_dp+h,0.1_dp,1.2_dp), &
            johnsonsu_loglikelihood(su_data,-0.6_dp,1.3_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('johnsonsu score loc',q4(3),centered( &
            johnsonsu_loglikelihood(su_data,-0.6_dp,1.3_dp,0.1_dp+h,1.2_dp), &
            johnsonsu_loglikelihood(su_data,-0.6_dp,1.3_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('johnsonsu score scale',q4(4),centered( &
            johnsonsu_loglikelihood(su_data,-0.6_dp,1.3_dp,0.1_dp,1.2_dp+h), &
            johnsonsu_loglikelihood(su_data,-0.6_dp,1.3_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)

        q4 = trapezoid_score(tr_data,0.25_dp,0.7_dp,0.1_dp,1.2_dp)
        call check_close('trapezoid score c',q4(1),centered( &
            trapezoid_loglikelihood(tr_data,0.25_dp+h,0.7_dp,0.1_dp,1.2_dp), &
            trapezoid_loglikelihood(tr_data,0.25_dp-h,0.7_dp,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('trapezoid score d',q4(2),centered( &
            trapezoid_loglikelihood(tr_data,0.25_dp,0.7_dp+h,0.1_dp,1.2_dp), &
            trapezoid_loglikelihood(tr_data,0.25_dp,0.7_dp-h,0.1_dp,1.2_dp),h),tol,tol,failures)
        call check_close('trapezoid score loc',q4(3),centered( &
            trapezoid_loglikelihood(tr_data,0.25_dp,0.7_dp,0.1_dp+h,1.2_dp), &
            trapezoid_loglikelihood(tr_data,0.25_dp,0.7_dp,0.1_dp-h,1.2_dp),h),tol,tol,failures)
        call check_close('trapezoid score scale',q4(4),centered( &
            trapezoid_loglikelihood(tr_data,0.25_dp,0.7_dp,0.1_dp,1.2_dp+h), &
            trapezoid_loglikelihood(tr_data,0.25_dp,0.7_dp,0.1_dp,1.2_dp-h),h),tol,tol,failures)
    end subroutine test_scores

    subroutine test_rng_fit_c_api(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        real(dp), parameter :: en_data(4) = [-0.4_dp,0.3_dp,1.2_dp,2.8_dp]
        real(dp), parameter :: sb_data(4) = [0.28_dp,0.52_dp,0.82_dp,1.12_dp]
        real(dp), parameter :: su_data(4) = [-2.0_dp,-0.3_dp,1.2_dp,3.0_dp]
        real(dp), parameter :: tr_data(4) = [0.22_dp,0.58_dp,0.76_dp,1.12_dp]
        real(dp) :: actual, expected
        type(rng_state) :: reference, state
        type(fit_result) :: result

        call rng_seed(state,97531)
        call rng_seed(reference,97531)
        actual = exponnorm_rvs(state,1.4_dp,0.1_dp,1.2_dp)
        expected = exponnorm_ppf(rng_uniform(reference),1.4_dp,0.1_dp,1.2_dp)
        call check_close('exponnorm rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = johnsonsb_rvs(state,0.4_dp,1.7_dp,0.1_dp,1.2_dp)
        expected = johnsonsb_ppf(rng_uniform(reference),0.4_dp,1.7_dp,0.1_dp,1.2_dp)
        call check_close('johnsonsb rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = johnsonsu_rvs(state,-0.6_dp,1.3_dp,0.1_dp,1.2_dp)
        expected = johnsonsu_ppf(rng_uniform(reference),-0.6_dp,1.3_dp,0.1_dp,1.2_dp)
        call check_close('johnsonsu rvs',actual,expected,0.0_dp,0.0_dp,failures)
        actual = trapezoid_rvs(state,0.25_dp,0.7_dp,0.1_dp,1.2_dp)
        expected = trapezoid_ppf(rng_uniform(reference),0.25_dp,0.7_dp,0.1_dp,1.2_dp)
        call check_close('trapezoid rvs',actual,expected,0.0_dp,0.0_dp,failures)

        call exponnorm_fit(en_data,[1.4_dp,0.1_dp,1.2_dp],[1.4_dp,0.1_dp,1.2_dp],result)
        call check_true('exponnorm fixed fit',result%success,failures)
        call johnsonsb_fit(sb_data,[0.4_dp,1.7_dp,0.1_dp,1.2_dp], &
            [0.4_dp,1.7_dp,0.1_dp,1.2_dp],result)
        call check_true('johnsonsb fixed fit',result%success,failures)
        call johnsonsu_fit(su_data,[-0.6_dp,1.3_dp,0.1_dp,1.2_dp], &
            [-0.6_dp,1.3_dp,0.1_dp,1.2_dp],result)
        call check_true('johnsonsu fixed fit',result%success,failures)
        call trapezoid_fit(tr_data,[0.25_dp,0.7_dp,0.1_dp,1.2_dp], &
            [0.25_dp,0.7_dp,0.1_dp,1.2_dp],result)
        call check_true('trapezoid fixed fit',result%success,failures)

        call check_close('C ABI exponnorm cdf',scifort_exponnorm_cdf_f64(0.8_dp,1.4_dp,0.0_dp,1.0_dp), &
            exponnorm_cdf(0.8_dp,1.4_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI johnsonsb pdf', &
            scifort_johnsonsb_pdf_f64(0.5_dp,0.4_dp,1.7_dp,0.0_dp,1.0_dp), &
            johnsonsb_pdf(0.5_dp,0.4_dp,1.7_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI johnsonsu ppf', &
            scifort_johnsonsu_ppf_f64(0.9_dp,-0.6_dp,1.3_dp,0.0_dp,1.0_dp), &
            johnsonsu_ppf(0.9_dp,-0.6_dp,1.3_dp),0.0_dp,0.0_dp,failures)
        call check_close('C ABI trapezoid cdf', &
            scifort_trapezoid_cdf_f64(0.8_dp,0.25_dp,0.7_dp,0.0_dp,1.0_dp), &
            trapezoid_cdf(0.8_dp,0.25_dp,0.7_dp),0.0_dp,0.0_dp,failures)
    end subroutine test_rng_fit_c_api

    pure function centered(plus, minus, h) result(value)
        real(dp), intent(in) :: plus !! objective at theta+h
        real(dp), intent(in) :: minus !! objective at theta-h
        real(dp), intent(in) :: h !! positive finite-difference step
        real(dp) :: value
        value = (plus - minus) / (2.0_dp * h)
    end function centered

end program test_normal_transform_four
