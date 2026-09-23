! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_matrix_distributions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_matrix_distributions_reference
    implicit none

    integer :: failures

    failures = 0
    call test_multivariate_hypergeom(failures)
    call test_normal_inverse_gamma(failures)
    call test_matrix_normal(failures)
    call test_wishart(failures)
    call test_invwishart(failures)
    call test_matrix_t(failures)

    if (failures /= 0) then
        print '(a,i0)', 'test_matrix_distributions: FAIL, failures=', failures
        error stop 1
    end if
    print '(a)', 'test_matrix_distributions: PASS'

contains

    subroutine test_multivariate_hypergeom(failures)
        integer, intent(inout) :: failures

        integer :: info
        integer :: sample1(3)
        integer :: sample2(3)
        integer :: samples(3, 5)
        integer :: j
        real(dp) :: cov(3,3)
        real(dp) :: mu(3)
        real(dp) :: v(3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('mhg logpmf', multivariate_hypergeom_logpmf(mhg_x_ref, mhg_m_ref, mhg_n_ref), &
            mhg_logpmf_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('mhg pmf', multivariate_hypergeom_pmf(mhg_x_ref, mhg_m_ref, mhg_n_ref), &
            mhg_pmf_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        mu = multivariate_hypergeom_mean(mhg_m_ref, mhg_n_ref)
        v = multivariate_hypergeom_var(mhg_m_ref, mhg_n_ref)
        cov = multivariate_hypergeom_cov(mhg_m_ref, mhg_n_ref)
        do j = 1, 3
            call check_close('mhg mean', mu(j), mhg_mean_ref(j), 2.0e-14_dp, 2.0e-14_dp, failures)
            call check_close('mhg var', v(j), mhg_var_ref(j), 2.0e-14_dp, 2.0e-14_dp, failures)
        end do
        do j = 1, 3
            call check_close('mhg cov row 1', cov(1,j), mhg_cov_ref(1,j), 2.0e-14_dp, 2.0e-14_dp, failures)
            call check_close('mhg cov row 2', cov(2,j), mhg_cov_ref(2,j), 2.0e-14_dp, 2.0e-14_dp, failures)
            call check_close('mhg cov row 3', cov(3,j), mhg_cov_ref(3,j), 2.0e-14_dp, 2.0e-14_dp, failures)
        end do
        call check_close('mhg off-support pmf', &
            multivariate_hypergeom_pmf([1, 2, 2], mhg_m_ref, mhg_n_ref), 0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_true('mhg invalid n is nan', &
            ieee_is_nan(multivariate_hypergeom_pmf(mhg_x_ref, mhg_m_ref, 30)), failures)
        call check_true('mhg zero-pop mean', all(multivariate_hypergeom_mean([0, 0], 0) == 0.0_dp), failures)
        call check_true('mhg zero-pop var', all(multivariate_hypergeom_var([0, 0], 0) == 0.0_dp), failures)
        call check_true('mhg one-pop cov', all(multivariate_hypergeom_cov([1, 0], 1) == 0.0_dp), failures)

        call rng_seed(state1, 314159)
        call rng_seed(state2, 314159)
        call multivariate_hypergeom_rvs(state1, mhg_m_ref, mhg_n_ref, sample1, info)
        call check_true('mhg rvs status', info == multivariate_hypergeom_status_success, failures)
        call multivariate_hypergeom_rvs(state2, mhg_m_ref, mhg_n_ref, sample2, info)
        call check_true('mhg rvs deterministic', all(sample1 == sample2), failures)
        call check_true('mhg rvs total', sum(sample1) == mhg_n_ref, failures)
        call check_true('mhg rvs support', all(sample1 >= 0) .and. all(sample1 <= mhg_m_ref), failures)
        call multivariate_hypergeom_rvs_array(state1, mhg_m_ref, mhg_n_ref, samples, info)
        call check_true('mhg rvs array status', info == multivariate_hypergeom_status_success, failures)
        do j = 1, size(samples,2)
            call check_true('mhg rvs array total', sum(samples(:,j)) == mhg_n_ref, failures)
        end do
    end subroutine test_multivariate_hypergeom

    subroutine test_normal_inverse_gamma(failures)
        integer, intent(inout) :: failures

        real(dp) :: mean_s2
        real(dp) :: mean_x
        real(dp) :: s21
        real(dp) :: s22
        real(dp) :: s2a(5)
        real(dp) :: var_s2
        real(dp) :: var_x
        real(dp) :: x1
        real(dp) :: x2
        real(dp) :: xa(5)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('nig logpdf', normal_inverse_gamma_logpdf(nig_x_ref, nig_s2_ref, nig_mu_ref, &
            nig_lmbda_ref, nig_a_ref, nig_b_ref), nig_logpdf_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('nig pdf', normal_inverse_gamma_pdf(nig_x_ref, nig_s2_ref, nig_mu_ref, &
            nig_lmbda_ref, nig_a_ref, nig_b_ref), nig_pdf_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        call normal_inverse_gamma_mean(nig_mu_ref, nig_lmbda_ref, nig_a_ref, nig_b_ref, mean_x, mean_s2)
        call normal_inverse_gamma_var(nig_mu_ref, nig_lmbda_ref, nig_a_ref, nig_b_ref, var_x, var_s2)
        call check_close('nig mean x', mean_x, nig_mean_ref(1), 0.0_dp, 0.0_dp, failures)
        call check_close('nig mean s2', mean_s2, nig_mean_ref(2), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('nig var x', var_x, nig_var_ref(1), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('nig var s2', var_s2, nig_var_ref(2), 2.0e-15_dp, 2.0e-15_dp, failures)
        call check_close('nig nonpositive s2 pdf', normal_inverse_gamma_pdf(0.0_dp, -1.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_true('nig invalid lambda nan', &
            ieee_is_nan(normal_inverse_gamma_pdf(0.0_dp, 1.0_dp, lmbda=-1.0_dp)), failures)

        call rng_seed(state1, 271828)
        call rng_seed(state2, 271828)
        call normal_inverse_gamma_rvs(state1, x1, s21, nig_mu_ref, nig_lmbda_ref, nig_a_ref, nig_b_ref)
        call normal_inverse_gamma_rvs(state2, x2, s22, nig_mu_ref, nig_lmbda_ref, nig_a_ref, nig_b_ref)
        call check_true('nig rvs deterministic', x1 == x2 .and. s21 == s22, failures)
        call check_true('nig rvs support', ieee_is_finite(x1) .and. s21 > 0.0_dp, failures)
        call normal_inverse_gamma_rvs_array(state1, xa, s2a, nig_mu_ref, nig_lmbda_ref, nig_a_ref, nig_b_ref)
        call check_true('nig rvs array finite', all(ieee_is_finite(xa)) .and. all(s2a > 0.0_dp), failures)
    end subroutine test_normal_inverse_gamma

    subroutine test_matrix_normal(failures)
        integer, intent(inout) :: failures

        integer :: info
        real(dp) :: row_upper(2,2)
        real(dp) :: sample1(2,2)
        real(dp) :: sample2(2,2)
        real(dp) :: samples(2,2,3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('matrix normal logpdf', matrix_normal_logpdf(matrix_x_ref, matrix_mean_ref, &
            matrix_row_ref, matrix_col_ref), matrix_normal_logpdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        call check_close('matrix normal pdf', matrix_normal_pdf(matrix_x_ref, matrix_mean_ref, &
            matrix_row_ref, matrix_col_ref), matrix_normal_pdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        call check_close('matrix normal entropy', matrix_normal_entropy(matrix_row_ref, matrix_col_ref), &
            matrix_normal_entropy_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        row_upper = matrix_row_ref
        row_upper(1,2) = 99.0_dp
        call check_close('matrix normal lower triangle', matrix_normal_logpdf(matrix_x_ref, matrix_mean_ref, &
            row_upper, matrix_col_ref), matrix_normal_logpdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)

        call rng_seed(state1, 161803)
        call rng_seed(state2, 161803)
        call matrix_normal_rvs(state1, sample1, matrix_mean_ref, matrix_row_ref, matrix_col_ref, info)
        call check_true('matrix normal rvs status', info == matrix_normal_status_success, failures)
        call matrix_normal_rvs(state2, sample2, matrix_mean_ref, matrix_row_ref, matrix_col_ref, info)
        call check_true('matrix normal rvs deterministic', all(sample1 == sample2), failures)
        call matrix_normal_rvs_array(state1, samples, matrix_mean_ref, matrix_row_ref, matrix_col_ref, info)
        call check_true('matrix normal rvs array status', info == matrix_normal_status_success, failures)
        call check_true('matrix normal rvs finite', all(ieee_is_finite(samples)), failures)
    end subroutine test_matrix_normal

    subroutine test_wishart(failures)
        integer, intent(inout) :: failures

        integer :: i
        integer :: info
        integer :: j
        real(dp) :: mean_value(2,2)
        real(dp) :: mode_bad(2,2)
        real(dp) :: mode_value(2,2)
        real(dp) :: sample1(2,2)
        real(dp) :: sample2(2,2)
        real(dp) :: samples(2,2,3)
        real(dp) :: scale_upper(2,2)
        real(dp) :: var_value(2,2)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('wishart logpdf', wishart_logpdf(wi_x_ref, wi_df_ref, wi_scale_ref), &
            wishart_logpdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        call check_close('wishart pdf', wishart_pdf(wi_x_ref, wi_df_ref, wi_scale_ref), &
            wishart_pdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        call check_close('wishart entropy', wishart_entropy(wi_df_ref, wi_scale_ref), &
            wishart_entropy_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        mean_value = wishart_mean(wi_df_ref, wi_scale_ref)
        mode_value = wishart_mode(wi_df_ref, wi_scale_ref)
        var_value = wishart_var(wi_df_ref, wi_scale_ref)
        do j = 1, 2
            do i = 1, 2
                call check_close('wishart mean', mean_value(i,j), wishart_mean_ref(i,j), &
                    3.0e-14_dp, 3.0e-14_dp, failures)
                call check_close('wishart mode', mode_value(i,j), wishart_mode_ref(i,j), &
                    3.0e-14_dp, 3.0e-14_dp, failures)
                call check_close('wishart var', var_value(i,j), wishart_var_ref(i,j), &
                    3.0e-14_dp, 3.0e-14_dp, failures)
            end do
        end do
        mode_bad = wishart_mode(2.5_dp, wi_scale_ref)
        call check_true('wishart undefined mode nan', all(ieee_is_nan(mode_bad)), failures)
        scale_upper = wi_scale_ref
        scale_upper(1,2) = 42.0_dp
        call check_close('wishart lower triangle', wishart_logpdf(wi_x_ref, wi_df_ref, scale_upper), &
            wishart_logpdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)

        call rng_seed(state1, 141421)
        call rng_seed(state2, 141421)
        call wishart_rvs(state1, wi_df_ref, wi_scale_ref, sample1, info)
        call check_true('wishart rvs status', info == wishart_status_success, failures)
        call wishart_rvs(state2, wi_df_ref, wi_scale_ref, sample2, info)
        call check_true('wishart rvs deterministic', all(sample1 == sample2), failures)
        call check_true('wishart rvs positive diag', sample1(1,1) > 0.0_dp .and. sample1(2,2) > 0.0_dp, failures)
        call wishart_rvs_array(state1, wi_df_ref, wi_scale_ref, samples, info)
        call check_true('wishart rvs array status', info == wishart_status_success, failures)
        call check_true('wishart rvs finite', all(ieee_is_finite(samples)), failures)
    end subroutine test_wishart

    subroutine test_invwishart(failures)
        integer, intent(inout) :: failures

        integer :: i
        integer :: info
        integer :: j
        real(dp) :: mean_value(2,2)
        real(dp) :: mode_value(2,2)
        real(dp) :: sample1(2,2)
        real(dp) :: sample2(2,2)
        real(dp) :: samples(2,2,3)
        real(dp) :: v_bad(2,2)
        real(dp) :: var_value(2,2)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('invwishart logpdf', invwishart_logpdf(wi_x_ref, wi_df_ref, wi_scale_ref), &
            invwishart_logpdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        call check_close('invwishart pdf', invwishart_pdf(wi_x_ref, wi_df_ref, wi_scale_ref), &
            invwishart_pdf_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        call check_close('invwishart entropy', invwishart_entropy(wi_df_ref, wi_scale_ref), &
            invwishart_entropy_ref, 3.0e-14_dp, 3.0e-14_dp, failures)
        mean_value = invwishart_mean(wi_df_ref, wi_scale_ref)
        mode_value = invwishart_mode(wi_df_ref, wi_scale_ref)
        var_value = invwishart_var(wi_df_ref, wi_scale_ref)
        do j = 1, 2
            do i = 1, 2
                call check_close('invwishart mean', mean_value(i,j), invwishart_mean_ref(i,j), &
                    3.0e-14_dp, 3.0e-14_dp, failures)
                call check_close('invwishart mode', mode_value(i,j), invwishart_mode_ref(i,j), &
                    3.0e-14_dp, 3.0e-14_dp, failures)
                call check_close('invwishart var', var_value(i,j), invwishart_var_ref(i,j), &
                    3.0e-14_dp, 3.0e-14_dp, failures)
            end do
        end do
        v_bad = invwishart_var(5.0_dp, wi_scale_ref)
        call check_true('invwishart undefined var nan', all(ieee_is_nan(v_bad)), failures)

        call rng_seed(state1, 173205)
        call rng_seed(state2, 173205)
        call invwishart_rvs(state1, wi_df_ref, wi_scale_ref, sample1, info)
        call check_true('invwishart rvs status', info == invwishart_status_success, failures)
        call invwishart_rvs(state2, wi_df_ref, wi_scale_ref, sample2, info)
        call check_true('invwishart rvs deterministic', all(sample1 == sample2), failures)
        call check_true('invwishart rvs positive diag', sample1(1,1) > 0.0_dp .and. sample1(2,2) > 0.0_dp, failures)
        call invwishart_rvs_array(state1, wi_df_ref, wi_scale_ref, samples, info)
        call check_true('invwishart rvs array status', info == invwishart_status_success, failures)
        call check_true('invwishart rvs finite', all(ieee_is_finite(samples)), failures)
    end subroutine test_invwishart

    subroutine test_matrix_t(failures)
        integer, intent(inout) :: failures

        integer :: info
        real(dp) :: row_upper(2,2)
        real(dp) :: sample1(2,2)
        real(dp) :: sample2(2,2)
        real(dp) :: samples(2,2,3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('matrix t logpdf', matrix_t_logpdf(matrix_x_ref, matrix_mean_ref, matrix_row_ref, &
            matrix_col_ref, matrix_t_df_ref), matrix_t_logpdf_ref, 4.0e-14_dp, 4.0e-14_dp, failures)
        call check_close('matrix t pdf', matrix_t_pdf(matrix_x_ref, matrix_mean_ref, matrix_row_ref, &
            matrix_col_ref, matrix_t_df_ref), matrix_t_pdf_ref, 4.0e-14_dp, 4.0e-14_dp, failures)
        call check_close('matrix t doc logpdf', matrix_t_logpdf(matrix_t_doc_x_ref, matrix_t_doc_mean_ref, &
            matrix_t_doc_row_ref, matrix_t_doc_col_ref, 3.0_dp), matrix_t_doc_logpdf_ref, &
            5.0e-14_dp, 5.0e-14_dp, failures)
        call check_close('matrix t doc pdf', matrix_t_pdf(matrix_t_doc_x_ref, matrix_t_doc_mean_ref, &
            matrix_t_doc_row_ref, matrix_t_doc_col_ref, 3.0_dp), matrix_t_doc_pdf_ref, &
            5.0e-14_dp, 5.0e-14_dp, failures)
        row_upper = matrix_row_ref
        row_upper(1,2) = 77.0_dp
        call check_close('matrix t lower triangle', matrix_t_logpdf(matrix_x_ref, matrix_mean_ref, row_upper, &
            matrix_col_ref, matrix_t_df_ref), matrix_t_logpdf_ref, 4.0e-14_dp, 4.0e-14_dp, failures)
        call check_true('matrix t invalid df nan', ieee_is_nan(matrix_t_pdf(matrix_x_ref, matrix_mean_ref, &
            matrix_row_ref, matrix_col_ref, -1.0_dp)), failures)

        call rng_seed(state1, 223607)
        call rng_seed(state2, 223607)
        call matrix_t_rvs(state1, sample1, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref, info)
        call check_true('matrix t rvs status', info == matrix_t_status_success, failures)
        call matrix_t_rvs(state2, sample2, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref, info)
        call check_true('matrix t rvs deterministic', all(sample1 == sample2), failures)
        call matrix_t_rvs_array(state1, samples, matrix_mean_ref, matrix_row_ref, matrix_col_ref, matrix_t_df_ref, info)
        call check_true('matrix t rvs array status', info == matrix_t_status_success, failures)
        call check_true('matrix t rvs finite', all(ieee_is_finite(samples)), failures)
    end subroutine test_matrix_t

end program test_matrix_distributions
