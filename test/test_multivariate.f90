! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_multivariate
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : cholesky_lower, linalg_status_success, psd_decompose, &
        solve_lower, symmetric_eigen_jacobi
    use scifort_normal, only : normal_pdf
    use scifort_random, only : rng_seed, rng_state
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_multivariate_reference
    implicit none

    integer :: failures

    failures = 0
    call test_linalg(failures)
    call test_multivariate_normal(failures)
    call test_multivariate_t(failures)
    call test_dirichlet(failures)
    call test_multinomial(failures)
    call test_dirichlet_multinomial(failures)
    call test_random_variates(failures)

    if (failures == 0) then
        print '(a)', 'test_multivariate: PASS'
    else
        print '(a,1x,i0)', 'test_multivariate: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_linalg(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp), parameter :: a(3,3) = reshape([ &
            4.0_dp, 1.0_dp, -0.5_dp, &
            1.0_dp, 3.0_dp, 0.25_dp, &
            -0.5_dp, 0.25_dp, 2.0_dp ], [3,3])
        real(dp), parameter :: b(3) = [1.0_dp, -2.0_dp, 0.5_dp]
        real(dp), parameter :: singular(2,2) = reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp],[2,2])
        real(dp) :: cutoff
        real(dp) :: eig(3)
        real(dp) :: eig2(2)
        real(dp) :: q(3,3)
        real(dp) :: q2(2,2)
        real(dp) :: l(3,3)
        real(dp) :: log_pdet
        real(dp) :: reconstructed(3,3)
        real(dp) :: rhs_check(3)
        real(dp) :: x(3)
        integer :: i
        integer :: j
        integer :: rank
        integer :: status

        call cholesky_lower(a, l, status)
        call check_true('linalg cholesky status', status == linalg_status_success, failures)
        reconstructed = matmul(l, transpose(l))
        do j = 1, 3
            do i = 1, 3
                call check_close('linalg cholesky reconstruction', reconstructed(i,j), a(i,j), &
                    5.0e-14_dp, 5.0e-14_dp, failures)
            end do
        end do

        call solve_lower(l, b, x, status)
        call check_true('linalg solve status', status == linalg_status_success, failures)
        rhs_check = matmul(l, x)
        do i = 1, 3
            call check_close('linalg triangular solve', rhs_check(i), b(i), 2.0e-14_dp, 2.0e-14_dp, failures)
        end do

        call symmetric_eigen_jacobi(a, eig, q, status)
        call check_true('linalg eigen status', status == linalg_status_success, failures)
        reconstructed = 0.0_dp
        do i = 1, 3
            reconstructed = reconstructed + eig(i) * outer_product(q(:,i), q(:,i))
        end do
        do j = 1, 3
            do i = 1, 3
                call check_close('linalg eigen reconstruction', reconstructed(i,j), a(i,j), &
                    3.0e-13_dp, 3.0e-13_dp, failures)
            end do
        end do

        call psd_decompose(singular, eig2, q2, rank, log_pdet, cutoff, status, .true.)
        call check_true('linalg singular psd status', status == linalg_status_success, failures)
        call check_true('linalg singular rank', rank == 1, failures)
        call check_close('linalg singular log pdet', log_pdet, log(2.0_dp), 2.0e-14_dp, 2.0e-14_dp, failures)
    end subroutine test_linalg

    subroutine test_multivariate_normal(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: cov_fit(2,2)
        real(dp) :: marginal_cov(2,2)
        real(dp) :: marginal_mean(2)
        real(dp) :: mean_fit(2)
        integer :: i
        integer :: j
        integer :: status

        call check_close('mvn pdf 2d', multivariate_normal_pdf(mvn1_x_ref,mvn1_mean_ref,mvn1_cov_ref), &
            mvn1_pdf_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('mvn logpdf 2d', multivariate_normal_logpdf(mvn1_x_ref,mvn1_mean_ref,mvn1_cov_ref), &
            mvn1_logpdf_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('mvn entropy 2d', multivariate_normal_entropy(mvn1_mean_ref,mvn1_cov_ref), &
            mvn1_entropy_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('mvn pdf 3d', multivariate_normal_pdf(mvn2_x_ref,mvn2_mean_ref,mvn2_cov_ref), &
            mvn2_pdf_ref, 3.0e-14_dp, 5.0e-13_dp, failures)
        call check_close('mvn logpdf 3d', multivariate_normal_logpdf(mvn2_x_ref,mvn2_mean_ref,mvn2_cov_ref), &
            mvn2_logpdf_ref, 5.0e-14_dp, 5.0e-13_dp, failures)
        call check_close('mvn entropy 3d', multivariate_normal_entropy(mvn2_mean_ref,mvn2_cov_ref), &
            mvn2_entropy_ref, 5.0e-14_dp, 5.0e-13_dp, failures)
        call check_close('mvn lower triangle semantics', &
            multivariate_normal_pdf(mvn1_x_ref,mvn1_mean_ref,mvn_lower_cov_ref), &
            mvn_lower_pdf_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('mvn singular support pdf', &
            multivariate_normal_pdf(mvn_singular_on_ref,cov=mvn_singular_cov_ref,allow_singular=.true.), &
            mvn_singular_on_pdf_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('mvn singular support logpdf', &
            multivariate_normal_logpdf(mvn_singular_on_ref,cov=mvn_singular_cov_ref,allow_singular=.true.), &
            mvn_singular_on_logpdf_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('mvn singular off support', &
            multivariate_normal_pdf(mvn_singular_off_ref,cov=mvn_singular_cov_ref,allow_singular=.true.), &
            mvn_singular_off_pdf_ref, 0.0_dp, 0.0_dp, failures)
        call check_true('mvn singular rejected by default', &
            ieee_is_nan(multivariate_normal_pdf(mvn_singular_on_ref,cov=mvn_singular_cov_ref)), failures)
        call check_close('mvn numerical-rank cutoff', &
            multivariate_normal_pdf([0.1_dp,0.0_dp],cov=reshape([1.0_dp,0.0_dp,0.0_dp,1.0e-12_dp],[2,2]), &
                allow_singular=.true.), normal_pdf(0.1_dp), 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_true('mvn numerical-rank cutoff rejected by default', &
            ieee_is_nan(multivariate_normal_pdf([0.1_dp,0.0_dp], &
                cov=reshape([1.0_dp,0.0_dp,0.0_dp,1.0e-12_dp],[2,2]))), failures)
        call check_true('mvn non-psd rejected', &
            ieee_is_nan(multivariate_normal_pdf([0.0_dp,0.0_dp], &
                cov=reshape([1.0_dp,0.0_dp,0.0_dp,-1.0e-3_dp],[2,2]),allow_singular=.true.)), failures)
        call check_close('mvn default one dimensional', multivariate_normal_pdf([0.3_dp]), &
            normal_pdf(0.3_dp), 3.0e-15_dp, 3.0e-15_dp, failures)

        call check_close('mvn cdf 2d', &
            multivariate_normal_cdf(mvn1_x_ref,mvn1_mean_ref,mvn1_cov_ref,maxpts=200000,abseps=1.0e-6_dp), &
            mvn1_cdf_ref, 1.2e-5_dp, 1.2e-5_dp, failures)
        call check_close('mvn logcdf 2d', &
            multivariate_normal_logcdf(mvn1_x_ref,mvn1_mean_ref,mvn1_cov_ref,maxpts=200000,abseps=1.0e-6_dp), &
            mvn1_logcdf_ref, 1.5e-5_dp, 1.5e-5_dp, failures)
        call check_close('mvn finite box cdf', &
            multivariate_normal_cdf(mvn1_x_ref,mvn1_mean_ref,mvn1_cov_ref,maxpts=200000,abseps=1.0e-6_dp, &
                lower_limit=mvn_box_lower_ref), &
            mvn1_box_cdf_ref, 1.2e-5_dp, 1.2e-5_dp, failures)
        call check_close('mvn cdf 3d', &
            multivariate_normal_cdf(mvn2_x_ref,mvn2_mean_ref,mvn2_cov_ref,maxpts=300000,abseps=1.0e-6_dp), &
            mvn2_cdf_ref, 1.5e-5_dp, 1.5e-5_dp, failures)
        call check_close('mvn singular cdf on line', &
            multivariate_normal_cdf(mvn_singular_on_ref,cov=mvn_singular_cov_ref,allow_singular=.true., &
                maxpts=200000,abseps=1.0e-7_dp), &
            mvn_singular_on_cdf_ref, 2.0e-12_dp, 2.0e-12_dp, failures)
        call check_close('mvn singular cdf off diagonal', &
            multivariate_normal_cdf(mvn_singular_off_ref,cov=mvn_singular_cov_ref,allow_singular=.true., &
                maxpts=200000,abseps=1.0e-7_dp), &
            mvn_singular_off_cdf_ref, 2.0e-12_dp, 2.0e-12_dp, failures)
        call check_close('mvn cdf 4d pivot path', &
            multivariate_normal_cdf(mv4_x_ref,mv4_loc_ref,mv4_cov_ref,maxpts=1000000,abseps=2.0e-6_dp), &
            mv4_cdf_ref, 2.0e-5_dp, 2.0e-5_dp, failures)
        call check_close('mvn finite box 4d pivot path', &
            multivariate_normal_cdf(mv4_x_ref,mv4_loc_ref,mv4_cov_ref,maxpts=1000000,abseps=2.0e-6_dp, &
                lower_limit=mv4_lower_ref), &
            mv4_box_cdf_ref, 2.0e-5_dp, 2.0e-5_dp, failures)

        call multivariate_normal_marginal([3,1],mvn2_mean_ref,mvn2_cov_ref,marginal_mean,marginal_cov,status)
        call check_true('mvn marginal status', status == multivariate_status_success, failures)
        do i = 1, 2
            call check_close('mvn marginal mean', marginal_mean(i), mvn_marg_mean_ref(i), &
                0.0_dp, 0.0_dp, failures)
            do j = 1, 2
                call check_close('mvn marginal covariance', marginal_cov(i,j), mvn_marg_cov_ref(i,j), &
                    0.0_dp, 0.0_dp, failures)
            end do
        end do

        call multivariate_normal_fit(mvn_fit_data_ref, mean_fit, cov_fit, status)
        call check_true('mvn fit status', status == multivariate_status_success, failures)
        do i = 1, 2
            call check_close('mvn fit mean', mean_fit(i), mvn_fit_mean_ref(i), 3.0e-15_dp, 3.0e-15_dp, failures)
            do j = 1, 2
                call check_close('mvn fit covariance', cov_fit(i,j), mvn_fit_cov_ref(i,j), &
                    3.0e-15_dp, 3.0e-15_dp, failures)
            end do
        end do
        call multivariate_normal_fit(mvn_fit_data_ref, mean_fit, cov_fit, status, fixed_mean=mvn_fixed_mean_ref)
        call check_true('mvn fixed mean fit status', status == multivariate_status_success, failures)
        do i = 1, 2
            call check_close('mvn fixed mean retained', mean_fit(i), mvn_fixed_mean_ref(i), 0.0_dp, 0.0_dp, failures)
            do j = 1, 2
                call check_close('mvn fixed mean covariance', cov_fit(i,j), mvn_fit_cov_fixed_mean_ref(i,j), &
                    3.0e-15_dp, 3.0e-15_dp, failures)
            end do
        end do
    end subroutine test_multivariate_normal

    subroutine test_multivariate_t(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: marginal_loc(2)
        real(dp) :: marginal_shape(2,2)
        integer :: i
        integer :: j
        integer :: status

        call check_close('mvt pdf', multivariate_t_pdf(mvt_x_ref,mvt_loc_ref,mvt_shape_ref,mvt_df_ref), &
            mvt_pdf_ref, 3.0e-15_dp, 5.0e-13_dp, failures)
        call check_close('mvt logpdf', multivariate_t_logpdf(mvt_x_ref,mvt_loc_ref,mvt_shape_ref,mvt_df_ref), &
            mvt_logpdf_ref, 3.0e-14_dp, 5.0e-13_dp, failures)
        call check_close('mvt entropy', multivariate_t_entropy(mvt_loc_ref,mvt_shape_ref,mvt_df_ref), &
            mvt_entropy_ref, 3.0e-14_dp, 5.0e-13_dp, failures)
        call check_close('mvt cdf', &
            multivariate_t_cdf(mvt_x_ref,mvt_loc_ref,mvt_shape_ref,mvt_df_ref,maxpts=200000), &
            mvt_cdf_ref, 1.5e-5_dp, 1.5e-5_dp, failures)
        call check_close('mvt finite box cdf', &
            multivariate_t_cdf(mvt_x_ref,mvt_loc_ref,mvt_shape_ref,mvt_df_ref,maxpts=200000, &
                lower_limit=mvt_lower_ref), &
            mvt_box_cdf_ref, 1.5e-5_dp, 1.5e-5_dp, failures)
        call check_close('mvt cdf 4d pivot path', &
            multivariate_t_cdf(mv4_x_ref,mv4_loc_ref,mv4_cov_ref,6.5_dp,maxpts=500000), &
            mvt4_cdf_ref, 2.0e-5_dp, 2.0e-5_dp, failures)
        call check_close('mvt finite box 4d pivot path', &
            multivariate_t_cdf(mv4_x_ref,mv4_loc_ref,mv4_cov_ref,6.5_dp,maxpts=500000, &
                lower_limit=mv4_lower_ref), &
            mvt4_box_cdf_ref, 2.0e-5_dp, 2.0e-5_dp, failures)
        call check_close('mvt singular pdf', &
            multivariate_t_pdf(mvt_singular_x_ref,shape=mvt_singular_shape_ref,df=4.0_dp,allow_singular=.true.), &
            mvt_singular_pdf_ref, 3.0e-14_dp, 5.0e-13_dp, failures)
        call check_close('mvt singular logpdf', &
            multivariate_t_logpdf(mvt_singular_x_ref,shape=mvt_singular_shape_ref,df=4.0_dp), &
            mvt_singular_logpdf_ref, 3.0e-14_dp, 5.0e-13_dp, failures)
        call check_true('mvt singular rejected by pdf default', &
            ieee_is_nan(multivariate_t_pdf(mvt_singular_x_ref,shape=mvt_singular_shape_ref,df=4.0_dp)), failures)
        call check_close('mvt one dimensional pdf identity', multivariate_t_pdf([0.3_dp],df=5.0_dp), &
            t_pdf(0.3_dp,5.0_dp), 3.0e-15_dp, 5.0e-14_dp, failures)
        call check_close('mvt one dimensional cdf identity', multivariate_t_cdf([0.3_dp],df=5.0_dp), &
            t_cdf(0.3_dp,5.0_dp), 3.0e-15_dp, 5.0e-14_dp, failures)

        call multivariate_t_marginal([2,1],mvt_loc_ref,mvt_shape_ref,marginal_loc,marginal_shape,status)
        call check_true('mvt marginal status', status == multivariate_t_status_success, failures)
        do i = 1, 2
            call check_close('mvt marginal location', marginal_loc(i), mvt_marg_loc_ref(i), &
                0.0_dp, 0.0_dp, failures)
            do j = 1, 2
                call check_close('mvt marginal shape', marginal_shape(i,j), mvt_marg_shape_ref(i,j), &
                    0.0_dp, 0.0_dp, failures)
            end do
        end do
    end subroutine test_multivariate_t

    subroutine test_dirichlet(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: cov(3,3)
        real(dp) :: mu(3)
        real(dp) :: v(3)
        integer :: i
        integer :: j

        call check_close('dirichlet pdf', dirichlet_pdf(dir_x_ref,dir_alpha_ref), &
            dir_pdf_ref, 2.0e-13_dp, 3.0e-13_dp, failures)
        call check_close('dirichlet logpdf', dirichlet_logpdf(dir_x_ref,dir_alpha_ref), &
            dir_logpdf_ref, 3.0e-14_dp, 3.0e-13_dp, failures)
        call check_close('dirichlet K-1 coordinates', dirichlet_pdf(dir_x_short_ref,dir_alpha_ref), &
            dir_pdf_ref, 2.0e-13_dp, 3.0e-13_dp, failures)
        call check_close('dirichlet entropy', dirichlet_entropy(dir_alpha_ref), &
            dir_entropy_ref, 2.0e-13_dp, 3.0e-13_dp, failures)
        mu = dirichlet_mean(dir_alpha_ref)
        v = dirichlet_var(dir_alpha_ref)
        cov = dirichlet_cov(dir_alpha_ref)
        do i = 1, 3
            call check_close('dirichlet mean', mu(i), dir_mean_ref(i), 2.0e-15_dp, 2.0e-15_dp, failures)
            call check_close('dirichlet variance', v(i), dir_var_ref(i), 2.0e-15_dp, 2.0e-15_dp, failures)
            do j = 1, 3
                call check_close('dirichlet covariance', cov(i,j), dir_cov_ref(i,j), &
                    3.0e-15_dp, 3.0e-15_dp, failures)
            end do
        end do
        call check_true('dirichlet invalid simplex', &
            ieee_is_nan(dirichlet_pdf([0.3_dp,0.3_dp,0.3_dp],dir_alpha_ref)), failures)
        call check_true('dirichlet zero alpha-lt-one invalid', &
            ieee_is_nan(dirichlet_pdf([0.0_dp,0.25_dp,0.75_dp],dir_alpha_ref)), failures)
    end subroutine test_dirichlet

    subroutine test_multinomial(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: cov(3,3)
        real(dp) :: mu(3)
        integer :: i
        integer :: j

        call check_close('multinomial pmf', multinomial_pmf(mn_x_ref,8,mn_p_ref), &
            mn_pmf_ref, 3.0e-15_dp, 5.0e-14_dp, failures)
        call check_close('multinomial logpmf', multinomial_logpmf(mn_x_ref,8,mn_p_ref), &
            mn_logpmf_ref, 3.0e-14_dp, 5.0e-14_dp, failures)
        call check_close('multinomial entropy', multinomial_entropy(8,mn_p_ref), &
            mn_entropy_ref, 3.0e-13_dp, 3.0e-13_dp, failures)
        mu = multinomial_mean(8,mn_p_ref)
        cov = multinomial_cov(8,mn_p_ref)
        do i = 1, 3
            call check_close('multinomial mean', mu(i), mn_mean_ref(i), 2.0e-15_dp, 2.0e-15_dp, failures)
            do j = 1, 3
                call check_close('multinomial covariance', cov(i,j), mn_cov_ref(i,j), &
                    3.0e-15_dp, 3.0e-15_dp, failures)
            end do
        end do
        call check_close('multinomial adjusted p pmf', multinomial_pmf(mn_adjust_x_ref,6,mn_adjust_p_ref), &
            mn_adjust_pmf_ref, 3.0e-15_dp, 5.0e-14_dp, failures)
        call check_close('multinomial adjusted p logpmf', multinomial_logpmf(mn_adjust_x_ref,6,mn_adjust_p_ref), &
            mn_adjust_logpmf_ref, 3.0e-14_dp, 5.0e-14_dp, failures)
        mu = multinomial_mean(6,mn_adjust_p_ref)
        do i = 1, 3
            call check_close('multinomial adjusted p mean', mu(i), mn_adjust_mean_ref(i), &
                2.0e-15_dp, 2.0e-15_dp, failures)
        end do
        call check_close('multinomial invalid counts mass', multinomial_pmf([1.0_dp,2.0_dp,4.0_dp],8,mn_p_ref), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_close('multinomial binomial identity', multinomial_pmf([3.0_dp,4.0_dp],7,[0.4_dp,0.6_dp]), &
            0.290304_dp, 3.0e-15_dp, 5.0e-14_dp, failures)
    end subroutine test_multinomial

    subroutine test_dirichlet_multinomial(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: cov(3,3)
        real(dp) :: mu(3)
        real(dp) :: v(3)
        integer :: i
        integer :: j

        call check_close('dirichlet multinomial pmf', dirichlet_multinomial_pmf(dm_x_ref,dm_alpha_ref,6), &
            dm_pmf_ref, 3.0e-15_dp, 5.0e-14_dp, failures)
        call check_close('dirichlet multinomial logpmf', dirichlet_multinomial_logpmf(dm_x_ref,dm_alpha_ref,6), &
            dm_logpmf_ref, 3.0e-14_dp, 5.0e-14_dp, failures)
        mu = dirichlet_multinomial_mean(dm_alpha_ref,6)
        v = dirichlet_multinomial_var(dm_alpha_ref,6)
        cov = dirichlet_multinomial_cov(dm_alpha_ref,6)
        do i = 1, 3
            call check_close('dirichlet multinomial mean', mu(i), dm_mean_ref(i), 2.0e-15_dp, 2.0e-15_dp, failures)
            call check_close('dirichlet multinomial variance', v(i), dm_var_ref(i), 3.0e-15_dp, 3.0e-15_dp, failures)
            do j = 1, 3
                call check_close('dirichlet multinomial covariance', cov(i,j), dm_cov_ref(i,j), &
                    3.0e-15_dp, 3.0e-15_dp, failures)
            end do
        end do
        call check_close('dirichlet multinomial wrong total mass', &
            dirichlet_multinomial_pmf(dm_x_ref,dm_alpha_ref,7),0.0_dp,0.0_dp,0.0_dp,failures)
        call check_true('dirichlet multinomial invalid alpha', &
            ieee_is_nan(dirichlet_multinomial_pmf(dm_x_ref,[3.0_dp,0.0_dp,5.0_dp],6)), failures)
    end subroutine test_dirichlet_multinomial

    subroutine test_random_variates(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: d1(3)
        real(dp) :: d2(3)
        real(dp) :: m1(2)
        real(dp) :: m2(2)
        real(dp) :: mt1(2)
        real(dp) :: mt2(2)
        integer :: counts1(3)
        integer :: counts2(3)
        integer :: status
        type(rng_state) :: state1
        type(rng_state) :: state2

        call rng_seed(state1, 12345)
        call rng_seed(state2, 12345)
        call multivariate_normal_rvs(state1,m1,mvn1_mean_ref,mvn1_cov_ref,status)
        call check_true('mvn rvs status', status == multivariate_status_success, failures)
        call multivariate_normal_rvs(state2,m2,mvn1_mean_ref,mvn1_cov_ref,status)
        call check_true('mvn rvs deterministic', all(m1 == m2) .and. all(ieee_is_finite(m1)), failures)

        call rng_seed(state1, 17321)
        call rng_seed(state2, 17321)
        call multivariate_t_rvs(state1,mt1,mvt_loc_ref,mvt_shape_ref,mvt_df_ref,status)
        call check_true('mvt rvs status', status == multivariate_t_status_success, failures)
        call multivariate_t_rvs(state2,mt2,mvt_loc_ref,mvt_shape_ref,mvt_df_ref,status)
        call check_true('mvt rvs deterministic', all(mt1 == mt2) .and. all(ieee_is_finite(mt1)), failures)

        call rng_seed(state1, 23456)
        call rng_seed(state2, 23456)
        call dirichlet_rvs(state1,d1,dir_alpha_ref,status)
        call check_true('dirichlet rvs status', status == dirichlet_status_success, failures)
        call dirichlet_rvs(state2,d2,dir_alpha_ref,status)
        call check_true('dirichlet rvs deterministic', all(d1 == d2), failures)
        call check_close('dirichlet rvs simplex', sum(d1), 1.0_dp, 3.0e-15_dp, 3.0e-15_dp, failures)
        call check_true('dirichlet rvs nonnegative', all(d1 >= 0.0_dp), failures)

        call rng_seed(state1, 34567)
        call rng_seed(state2, 34567)
        call multinomial_rvs(state1,counts1,8,mn_p_ref,status)
        call check_true('multinomial rvs status', status == multinomial_status_success, failures)
        call multinomial_rvs(state2,counts2,8,mn_p_ref,status)
        call check_true('multinomial rvs deterministic', all(counts1 == counts2), failures)
        call check_true('multinomial rvs support', sum(counts1) == 8 .and. all(counts1 >= 0), failures)
    end subroutine test_random_variates

    pure function outer_product(x, y) result(a)
        real(dp), intent(in) :: x(:) !! left vector
        real(dp), intent(in) :: y(:) !! right vector
        real(dp) :: a(size(x),size(y))

        integer :: i
        integer :: j

        do j = 1, size(y)
            do i = 1, size(x)
                a(i,j) = x(i) * y(j)
            end do
        end do
    end function outer_product

end program test_multivariate
