# Reviewed interface coverage

Reviewed C batch and Python policy only; not a whole-project inventory.

| Interface | Count |
| --- | ---: |
| exports | 288 |
| scalar | 88 |
| bulk | 88 |
| array_adapters | 122 |
| elementwise | 88 |
| python | 210 |

| Module | Export | Scalar C | Bulk C | Array adapter | Elementwise C | Python |
| --- | --- | --- | --- | --- | --- | --- |
| scifort_beta | beta_cdf | scifort_beta_cdf_f64 | scifort_beta_cdf_vec_f64 |  | scifort_beta_cdf_elementwise_f64 | beta.cdf |
| scifort_beta | beta_isf | scifort_beta_isf_f64 | scifort_beta_isf_vec_f64 |  | scifort_beta_isf_elementwise_f64 | beta.isf |
| scifort_beta | beta_logcdf | scifort_beta_logcdf_f64 | scifort_beta_logcdf_vec_f64 |  | scifort_beta_logcdf_elementwise_f64 | beta.logcdf |
| scifort_beta | beta_logpdf | scifort_beta_logpdf_f64 | scifort_beta_logpdf_vec_f64 |  | scifort_beta_logpdf_elementwise_f64 | beta.logpdf |
| scifort_beta | beta_logsf | scifort_beta_logsf_f64 | scifort_beta_logsf_vec_f64 |  | scifort_beta_logsf_elementwise_f64 | beta.logsf |
| scifort_beta | beta_pdf | scifort_beta_pdf_f64 | scifort_beta_pdf_vec_f64 |  | scifort_beta_pdf_elementwise_f64 | beta.pdf |
| scifort_beta | beta_ppf | scifort_beta_ppf_f64 | scifort_beta_ppf_vec_f64 |  | scifort_beta_ppf_elementwise_f64 | beta.ppf |
| scifort_beta | beta_sf | scifort_beta_sf_f64 | scifort_beta_sf_vec_f64 |  | scifort_beta_sf_elementwise_f64 | beta.sf |
| scifort_cauchy | cauchy_cdf | scifort_cauchy_cdf_f64 | scifort_cauchy_cdf_vec_f64 |  | scifort_cauchy_cdf_elementwise_f64 | cauchy.cdf |
| scifort_cauchy | cauchy_isf | scifort_cauchy_isf_f64 | scifort_cauchy_isf_vec_f64 |  | scifort_cauchy_isf_elementwise_f64 | cauchy.isf |
| scifort_cauchy | cauchy_logcdf | scifort_cauchy_logcdf_f64 | scifort_cauchy_logcdf_vec_f64 |  | scifort_cauchy_logcdf_elementwise_f64 | cauchy.logcdf |
| scifort_cauchy | cauchy_logpdf | scifort_cauchy_logpdf_f64 | scifort_cauchy_logpdf_vec_f64 |  | scifort_cauchy_logpdf_elementwise_f64 | cauchy.logpdf |
| scifort_cauchy | cauchy_logsf | scifort_cauchy_logsf_f64 | scifort_cauchy_logsf_vec_f64 |  | scifort_cauchy_logsf_elementwise_f64 | cauchy.logsf |
| scifort_cauchy | cauchy_pdf | scifort_cauchy_pdf_f64 | scifort_cauchy_pdf_vec_f64 |  | scifort_cauchy_pdf_elementwise_f64 | cauchy.pdf |
| scifort_cauchy | cauchy_ppf | scifort_cauchy_ppf_f64 | scifort_cauchy_ppf_vec_f64 |  | scifort_cauchy_ppf_elementwise_f64 | cauchy.ppf |
| scifort_cauchy | cauchy_sf | scifort_cauchy_sf_f64 | scifort_cauchy_sf_vec_f64 |  | scifort_cauchy_sf_elementwise_f64 | cauchy.sf |
| scifort_chi2 | chi2_cdf | scifort_chi2_cdf_f64 | scifort_chi2_cdf_vec_f64 |  | scifort_chi2_cdf_elementwise_f64 | chi2.cdf |
| scifort_chi2 | chi2_isf | scifort_chi2_isf_f64 | scifort_chi2_isf_vec_f64 |  | scifort_chi2_isf_elementwise_f64 | chi2.isf |
| scifort_chi2 | chi2_logcdf | scifort_chi2_logcdf_f64 | scifort_chi2_logcdf_vec_f64 |  | scifort_chi2_logcdf_elementwise_f64 | chi2.logcdf |
| scifort_chi2 | chi2_logpdf | scifort_chi2_logpdf_f64 | scifort_chi2_logpdf_vec_f64 |  | scifort_chi2_logpdf_elementwise_f64 | chi2.logpdf |
| scifort_chi2 | chi2_logsf | scifort_chi2_logsf_f64 | scifort_chi2_logsf_vec_f64 |  | scifort_chi2_logsf_elementwise_f64 | chi2.logsf |
| scifort_chi2 | chi2_pdf | scifort_chi2_pdf_f64 | scifort_chi2_pdf_vec_f64 |  | scifort_chi2_pdf_elementwise_f64 | chi2.pdf |
| scifort_chi2 | chi2_ppf | scifort_chi2_ppf_f64 | scifort_chi2_ppf_vec_f64 |  | scifort_chi2_ppf_elementwise_f64 | chi2.ppf |
| scifort_chi2 | chi2_sf | scifort_chi2_sf_f64 | scifort_chi2_sf_vec_f64 |  | scifort_chi2_sf_elementwise_f64 | chi2.sf |
| scifort_dirichlet | dirichlet_cov |  |  | scifort_dirichlet_cov_f64 |  | dirichlet.cov |
| scifort_dirichlet | dirichlet_entropy |  |  | scifort_dirichlet_entropy_f64 |  | dirichlet.entropy |
| scifort_dirichlet | dirichlet_logpdf |  |  | scifort_dirichlet_logpdf_f64 |  | dirichlet.logpdf |
| scifort_dirichlet | dirichlet_mean |  |  | scifort_dirichlet_mean_f64 |  | dirichlet.mean |
| scifort_dirichlet | dirichlet_pdf |  |  | scifort_dirichlet_pdf_f64 |  | dirichlet.pdf |
| scifort_dirichlet | dirichlet_rvs |  |  |  |  |  |
| scifort_dirichlet | dirichlet_rvs_array |  |  | scifort_dirichlet_rvs_array_f64 |  | dirichlet.rvs_array |
| scifort_dirichlet | dirichlet_status_invalid_parameter |  |  |  |  |  |
| scifort_dirichlet | dirichlet_status_invalid_shape |  |  |  |  |  |
| scifort_dirichlet | dirichlet_status_success |  |  |  |  |  |
| scifort_dirichlet | dirichlet_var |  |  | scifort_dirichlet_var_f64 |  | dirichlet.var |
| scifort_dirichlet_multinomial | dirichlet_multinomial_cov |  |  | scifort_dirichlet_multinomial_cov_f64 |  | dirichlet_multinomial.cov |
| scifort_dirichlet_multinomial | dirichlet_multinomial_logpmf |  |  | scifort_dirichlet_multinomial_logpmf_f64 |  | dirichlet_multinomial.logpmf |
| scifort_dirichlet_multinomial | dirichlet_multinomial_mean |  |  | scifort_dirichlet_multinomial_mean_f64 |  | dirichlet_multinomial.mean |
| scifort_dirichlet_multinomial | dirichlet_multinomial_pmf |  |  | scifort_dirichlet_multinomial_pmf_f64 |  | dirichlet_multinomial.pmf |
| scifort_dirichlet_multinomial | dirichlet_multinomial_var |  |  | scifort_dirichlet_multinomial_var_f64 |  | dirichlet_multinomial.var |
| scifort_ecdf_bridge | ecdf_create |  |  | scifort_ecdf_create_f64 |  | ecdf_bridge.create |
| scifort_ecdf_bridge | ecdf_destroy |  |  | scifort_ecdf_destroy_f64 |  | ecdf_bridge.destroy |
| scifort_ecdf_bridge | ecdf_evaluate |  |  | scifort_ecdf_evaluate_f64 |  | ecdf_bridge.evaluate |
| scifort_ecdf_bridge | ecdf_info |  |  | scifort_ecdf_info_f64 |  | ecdf_bridge.info |
| scifort_ecdf_bridge | ecdf_interval |  |  | scifort_ecdf_interval_f64 |  | ecdf_bridge.interval |
| scifort_ecdf_bridge | ecdf_values |  |  | scifort_ecdf_values_f64 |  | ecdf_bridge.values |
| scifort_exponential | exponential_cdf | scifort_exponential_cdf_f64 | scifort_exponential_cdf_vec_f64 |  | scifort_exponential_cdf_elementwise_f64 | exponential.cdf |
| scifort_exponential | exponential_isf | scifort_exponential_isf_f64 | scifort_exponential_isf_vec_f64 |  | scifort_exponential_isf_elementwise_f64 | exponential.isf |
| scifort_exponential | exponential_logcdf | scifort_exponential_logcdf_f64 | scifort_exponential_logcdf_vec_f64 |  | scifort_exponential_logcdf_elementwise_f64 | exponential.logcdf |
| scifort_exponential | exponential_logpdf | scifort_exponential_logpdf_f64 | scifort_exponential_logpdf_vec_f64 |  | scifort_exponential_logpdf_elementwise_f64 | exponential.logpdf |
| scifort_exponential | exponential_logsf | scifort_exponential_logsf_f64 | scifort_exponential_logsf_vec_f64 |  | scifort_exponential_logsf_elementwise_f64 | exponential.logsf |
| scifort_exponential | exponential_pdf | scifort_exponential_pdf_f64 | scifort_exponential_pdf_vec_f64 |  | scifort_exponential_pdf_elementwise_f64 | exponential.pdf |
| scifort_exponential | exponential_ppf | scifort_exponential_ppf_f64 | scifort_exponential_ppf_vec_f64 |  | scifort_exponential_ppf_elementwise_f64 | exponential.ppf |
| scifort_exponential | exponential_sf | scifort_exponential_sf_f64 | scifort_exponential_sf_vec_f64 |  | scifort_exponential_sf_elementwise_f64 | exponential.sf |
| scifort_gamma | gamma_cdf | scifort_gamma_cdf_f64 | scifort_gamma_cdf_vec_f64 |  | scifort_gamma_cdf_elementwise_f64 | gamma.cdf |
| scifort_gamma | gamma_isf | scifort_gamma_isf_f64 | scifort_gamma_isf_vec_f64 |  | scifort_gamma_isf_elementwise_f64 | gamma.isf |
| scifort_gamma | gamma_logcdf | scifort_gamma_logcdf_f64 | scifort_gamma_logcdf_vec_f64 |  | scifort_gamma_logcdf_elementwise_f64 | gamma.logcdf |
| scifort_gamma | gamma_logpdf | scifort_gamma_logpdf_f64 | scifort_gamma_logpdf_vec_f64 |  | scifort_gamma_logpdf_elementwise_f64 | gamma.logpdf |
| scifort_gamma | gamma_logsf | scifort_gamma_logsf_f64 | scifort_gamma_logsf_vec_f64 |  | scifort_gamma_logsf_elementwise_f64 | gamma.logsf |
| scifort_gamma | gamma_pdf | scifort_gamma_pdf_f64 | scifort_gamma_pdf_vec_f64 |  | scifort_gamma_pdf_elementwise_f64 | gamma.pdf |
| scifort_gamma | gamma_ppf | scifort_gamma_ppf_f64 | scifort_gamma_ppf_vec_f64 |  | scifort_gamma_ppf_elementwise_f64 | gamma.ppf |
| scifort_gamma | gamma_sf | scifort_gamma_sf_f64 | scifort_gamma_sf_vec_f64 |  | scifort_gamma_sf_elementwise_f64 | gamma.sf |
| scifort_gaussian_kde_bridge | kde_covariance |  |  | scifort_gaussian_kde_covariance_f64 |  | kde_bridge.covariance |
| scifort_gaussian_kde_bridge | kde_create |  |  | scifort_gaussian_kde_create_f64 |  | kde_bridge.create |
| scifort_gaussian_kde_bridge | kde_destroy |  |  | scifort_gaussian_kde_destroy_f64 |  | kde_bridge.destroy |
| scifort_gaussian_kde_bridge | kde_info |  |  | scifort_gaussian_kde_info_f64 |  | kde_bridge.info |
| scifort_gaussian_kde_bridge | kde_integrate_box |  |  | scifort_gaussian_kde_integrate_box_f64 |  | kde_bridge.integrate_box |
| scifort_gaussian_kde_bridge | kde_integrate_box_1d |  |  | scifort_gaussian_kde_integrate_box_1d_f64 |  | kde_bridge.integrate_box_1d |
| scifort_gaussian_kde_bridge | kde_integrate_gaussian |  |  | scifort_gaussian_kde_integrate_gaussian_f64 |  | kde_bridge.integrate_gaussian |
| scifort_gaussian_kde_bridge | kde_integrate_kde |  |  | scifort_gaussian_kde_integrate_kde_f64 |  | kde_bridge.integrate_kde |
| scifort_gaussian_kde_bridge | kde_inv_cov |  |  | scifort_gaussian_kde_inv_cov_f64 |  | kde_bridge.inv_cov |
| scifort_gaussian_kde_bridge | kde_logpdf |  |  | scifort_gaussian_kde_logpdf_f64 |  | kde_bridge.logpdf |
| scifort_gaussian_kde_bridge | kde_marginal |  |  | scifort_gaussian_kde_marginal_f64 |  | kde_bridge.marginal |
| scifort_gaussian_kde_bridge | kde_pdf |  |  | scifort_gaussian_kde_pdf_f64 |  | kde_bridge.pdf |
| scifort_gaussian_kde_bridge | kde_resample |  |  | scifort_gaussian_kde_resample_f64 |  | kde_bridge.resample |
| scifort_gaussian_kde_bridge | kde_set_bandwidth |  |  | scifort_gaussian_kde_set_bandwidth_f64 |  | kde_bridge.set_bandwidth |
| scifort_invwishart | invwishart_entropy |  |  | scifort_invwishart_entropy_f64 |  | invwishart.entropy |
| scifort_invwishart | invwishart_logpdf |  |  | scifort_invwishart_logpdf_f64 |  | invwishart.logpdf |
| scifort_invwishart | invwishart_mean |  |  | scifort_invwishart_mean_f64 |  | invwishart.mean |
| scifort_invwishart | invwishart_mode |  |  | scifort_invwishart_mode_f64 |  | invwishart.mode |
| scifort_invwishart | invwishart_pdf |  |  | scifort_invwishart_pdf_f64 |  | invwishart.pdf |
| scifort_invwishart | invwishart_rvs |  |  |  |  |  |
| scifort_invwishart | invwishart_rvs_array |  |  | scifort_invwishart_rvs_array_f64 |  | invwishart.rvs_array |
| scifort_invwishart | invwishart_status_invalid_parameter |  |  |  |  |  |
| scifort_invwishart | invwishart_status_invalid_shape |  |  |  |  |  |
| scifort_invwishart | invwishart_status_linalg_failure |  |  |  |  |  |
| scifort_invwishart | invwishart_status_success |  |  |  |  |  |
| scifort_invwishart | invwishart_var |  |  | scifort_invwishart_var_f64 |  | invwishart.var |
| scifort_laplace | laplace_cdf | scifort_laplace_cdf_f64 | scifort_laplace_cdf_vec_f64 |  | scifort_laplace_cdf_elementwise_f64 | laplace.cdf |
| scifort_laplace | laplace_isf | scifort_laplace_isf_f64 | scifort_laplace_isf_vec_f64 |  | scifort_laplace_isf_elementwise_f64 | laplace.isf |
| scifort_laplace | laplace_logcdf | scifort_laplace_logcdf_f64 | scifort_laplace_logcdf_vec_f64 |  | scifort_laplace_logcdf_elementwise_f64 | laplace.logcdf |
| scifort_laplace | laplace_logpdf | scifort_laplace_logpdf_f64 | scifort_laplace_logpdf_vec_f64 |  | scifort_laplace_logpdf_elementwise_f64 | laplace.logpdf |
| scifort_laplace | laplace_logsf | scifort_laplace_logsf_f64 | scifort_laplace_logsf_vec_f64 |  | scifort_laplace_logsf_elementwise_f64 | laplace.logsf |
| scifort_laplace | laplace_pdf | scifort_laplace_pdf_f64 | scifort_laplace_pdf_vec_f64 |  | scifort_laplace_pdf_elementwise_f64 | laplace.pdf |
| scifort_laplace | laplace_ppf | scifort_laplace_ppf_f64 | scifort_laplace_ppf_vec_f64 |  | scifort_laplace_ppf_elementwise_f64 | laplace.ppf |
| scifort_laplace | laplace_sf | scifort_laplace_sf_f64 | scifort_laplace_sf_vec_f64 |  | scifort_laplace_sf_elementwise_f64 | laplace.sf |
| scifort_linalg | cholesky_lower |  |  | scifort_cholesky_lower_f64 |  | linalg.cholesky_lower |
| scifort_linalg | linalg_status_invalid_shape |  |  |  |  |  |
| scifort_linalg | linalg_status_no_convergence |  |  |  |  |  |
| scifort_linalg | linalg_status_not_positive_definite |  |  |  |  |  |
| scifort_linalg | linalg_status_not_positive_semidefinite |  |  |  |  |  |
| scifort_linalg | linalg_status_singular |  |  |  |  |  |
| scifort_linalg | linalg_status_success |  |  |  |  |  |
| scifort_linalg | psd_decompose |  |  | scifort_psd_decompose_f64 |  | linalg.psd_decompose |
| scifort_linalg | solve_lower |  |  | scifort_solve_lower_f64 |  | linalg.solve_lower |
| scifort_linalg | symmetric_eigen_jacobi |  |  | scifort_symmetric_eigen_jacobi_f64 |  | linalg.symmetric_eigen_jacobi |
| scifort_logistic | logistic_cdf | scifort_logistic_cdf_f64 | scifort_logistic_cdf_vec_f64 |  | scifort_logistic_cdf_elementwise_f64 | logistic.cdf |
| scifort_logistic | logistic_isf | scifort_logistic_isf_f64 | scifort_logistic_isf_vec_f64 |  | scifort_logistic_isf_elementwise_f64 | logistic.isf |
| scifort_logistic | logistic_logcdf | scifort_logistic_logcdf_f64 | scifort_logistic_logcdf_vec_f64 |  | scifort_logistic_logcdf_elementwise_f64 | logistic.logcdf |
| scifort_logistic | logistic_logpdf | scifort_logistic_logpdf_f64 | scifort_logistic_logpdf_vec_f64 |  | scifort_logistic_logpdf_elementwise_f64 | logistic.logpdf |
| scifort_logistic | logistic_logsf | scifort_logistic_logsf_f64 | scifort_logistic_logsf_vec_f64 |  | scifort_logistic_logsf_elementwise_f64 | logistic.logsf |
| scifort_logistic | logistic_pdf | scifort_logistic_pdf_f64 | scifort_logistic_pdf_vec_f64 |  | scifort_logistic_pdf_elementwise_f64 | logistic.pdf |
| scifort_logistic | logistic_ppf | scifort_logistic_ppf_f64 | scifort_logistic_ppf_vec_f64 |  | scifort_logistic_ppf_elementwise_f64 | logistic.ppf |
| scifort_logistic | logistic_sf | scifort_logistic_sf_f64 | scifort_logistic_sf_vec_f64 |  | scifort_logistic_sf_elementwise_f64 | logistic.sf |
| scifort_lognormal | lognormal_cdf | scifort_lognormal_cdf_f64 | scifort_lognormal_cdf_vec_f64 |  | scifort_lognormal_cdf_elementwise_f64 | lognormal.cdf |
| scifort_lognormal | lognormal_isf | scifort_lognormal_isf_f64 | scifort_lognormal_isf_vec_f64 |  | scifort_lognormal_isf_elementwise_f64 | lognormal.isf |
| scifort_lognormal | lognormal_logcdf | scifort_lognormal_logcdf_f64 | scifort_lognormal_logcdf_vec_f64 |  | scifort_lognormal_logcdf_elementwise_f64 | lognormal.logcdf |
| scifort_lognormal | lognormal_logpdf | scifort_lognormal_logpdf_f64 | scifort_lognormal_logpdf_vec_f64 |  | scifort_lognormal_logpdf_elementwise_f64 | lognormal.logpdf |
| scifort_lognormal | lognormal_logsf | scifort_lognormal_logsf_f64 | scifort_lognormal_logsf_vec_f64 |  | scifort_lognormal_logsf_elementwise_f64 | lognormal.logsf |
| scifort_lognormal | lognormal_pdf | scifort_lognormal_pdf_f64 | scifort_lognormal_pdf_vec_f64 |  | scifort_lognormal_pdf_elementwise_f64 | lognormal.pdf |
| scifort_lognormal | lognormal_ppf | scifort_lognormal_ppf_f64 | scifort_lognormal_ppf_vec_f64 |  | scifort_lognormal_ppf_elementwise_f64 | lognormal.ppf |
| scifort_lognormal | lognormal_sf | scifort_lognormal_sf_f64 | scifort_lognormal_sf_vec_f64 |  | scifort_lognormal_sf_elementwise_f64 | lognormal.sf |
| scifort_matrix_distribution_helpers | invert_lower |  |  | scifort_invert_lower_f64 |  | linalg.invert_lower |
| scifort_matrix_distribution_helpers | log_multivariate_gamma |  |  |  |  |  |
| scifort_matrix_distribution_helpers | matrix_status_invalid_parameter |  |  |  |  |  |
| scifort_matrix_distribution_helpers | matrix_status_invalid_shape |  |  |  |  |  |
| scifort_matrix_distribution_helpers | matrix_status_not_positive_definite |  |  |  |  |  |
| scifort_matrix_distribution_helpers | matrix_status_success |  |  |  |  |  |
| scifort_matrix_distribution_helpers | multivariate_digamma_sum |  |  |  |  |  |
| scifort_matrix_distribution_helpers | solve_lower_matrix |  |  | scifort_solve_lower_matrix_f64 |  | linalg.solve_lower_matrix |
| scifort_matrix_distribution_helpers | solve_spd_matrix |  |  | scifort_solve_spd_matrix_f64 |  | linalg.solve_spd_matrix |
| scifort_matrix_distribution_helpers | spd_factor |  |  | scifort_spd_factor_f64 |  | linalg.spd_factor |
| scifort_matrix_normal | matrix_normal_entropy |  |  | scifort_matrix_normal_entropy_f64 |  | matrix_normal.entropy |
| scifort_matrix_normal | matrix_normal_logpdf |  |  | scifort_matrix_normal_logpdf_f64 |  | matrix_normal.logpdf |
| scifort_matrix_normal | matrix_normal_pdf |  |  | scifort_matrix_normal_pdf_f64 |  | matrix_normal.pdf |
| scifort_matrix_normal | matrix_normal_rvs |  |  |  |  |  |
| scifort_matrix_normal | matrix_normal_rvs_array |  |  | scifort_matrix_normal_rvs_array_f64 |  | matrix_normal.rvs_array |
| scifort_matrix_normal | matrix_normal_status_invalid_parameter |  |  |  |  |  |
| scifort_matrix_normal | matrix_normal_status_invalid_shape |  |  |  |  |  |
| scifort_matrix_normal | matrix_normal_status_linalg_failure |  |  |  |  |  |
| scifort_matrix_normal | matrix_normal_status_success |  |  |  |  |  |
| scifort_matrix_t | matrix_t_logpdf |  |  | scifort_matrix_t_logpdf_f64 |  | matrix_t.logpdf |
| scifort_matrix_t | matrix_t_pdf |  |  | scifort_matrix_t_pdf_f64 |  | matrix_t.pdf |
| scifort_matrix_t | matrix_t_rvs |  |  |  |  |  |
| scifort_matrix_t | matrix_t_rvs_array |  |  | scifort_matrix_t_rvs_array_f64 |  | matrix_t.rvs_array |
| scifort_matrix_t | matrix_t_status_invalid_parameter |  |  |  |  |  |
| scifort_matrix_t | matrix_t_status_invalid_shape |  |  |  |  |  |
| scifort_matrix_t | matrix_t_status_linalg_failure |  |  |  |  |  |
| scifort_matrix_t | matrix_t_status_success |  |  |  |  |  |
| scifort_multinomial | multinomial_cov |  |  | scifort_multinomial_cov_f64 |  | multinomial.cov |
| scifort_multinomial | multinomial_entropy |  |  | scifort_multinomial_entropy_f64 |  | multinomial.entropy |
| scifort_multinomial | multinomial_logpmf |  |  | scifort_multinomial_logpmf_f64 |  | multinomial.logpmf |
| scifort_multinomial | multinomial_mean |  |  | scifort_multinomial_mean_f64 |  | multinomial.mean |
| scifort_multinomial | multinomial_pmf |  |  | scifort_multinomial_pmf_f64 |  | multinomial.pmf |
| scifort_multinomial | multinomial_rvs |  |  |  |  |  |
| scifort_multinomial | multinomial_rvs_array |  |  | scifort_multinomial_rvs_array_i32 |  | multinomial.rvs_array |
| scifort_multinomial | multinomial_status_invalid_parameter |  |  |  |  |  |
| scifort_multinomial | multinomial_status_invalid_shape |  |  |  |  |  |
| scifort_multinomial | multinomial_status_success |  |  |  |  |  |
| scifort_multivariate_hypergeom | multivariate_hypergeom_cov |  |  | scifort_multivariate_hypergeom_cov_f64 |  | multivariate_hypergeom.cov |
| scifort_multivariate_hypergeom | multivariate_hypergeom_logpmf |  |  | scifort_multivariate_hypergeom_logpmf_f64 |  | multivariate_hypergeom.logpmf |
| scifort_multivariate_hypergeom | multivariate_hypergeom_mean |  |  | scifort_multivariate_hypergeom_mean_f64 |  | multivariate_hypergeom.mean |
| scifort_multivariate_hypergeom | multivariate_hypergeom_pmf |  |  | scifort_multivariate_hypergeom_pmf_f64 |  | multivariate_hypergeom.pmf |
| scifort_multivariate_hypergeom | multivariate_hypergeom_rvs |  |  |  |  |  |
| scifort_multivariate_hypergeom | multivariate_hypergeom_rvs_array |  |  | scifort_multivariate_hypergeom_rvs_array_i32 |  | multivariate_hypergeom.rvs_array |
| scifort_multivariate_hypergeom | multivariate_hypergeom_status_invalid_parameter |  |  |  |  |  |
| scifort_multivariate_hypergeom | multivariate_hypergeom_status_invalid_shape |  |  |  |  |  |
| scifort_multivariate_hypergeom | multivariate_hypergeom_status_success |  |  |  |  |  |
| scifort_multivariate_hypergeom | multivariate_hypergeom_var |  |  | scifort_multivariate_hypergeom_var_f64 |  | multivariate_hypergeom.var |
| scifort_multivariate_normal | multivariate_normal_cdf |  |  | scifort_multivariate_normal_cdf_f64 |  | multivariate_normal.cdf |
| scifort_multivariate_normal | multivariate_normal_entropy |  |  | scifort_multivariate_normal_entropy_f64 |  | multivariate_normal.entropy |
| scifort_multivariate_normal | multivariate_normal_fit |  |  | scifort_multivariate_normal_fit_f64 |  | multivariate_normal.fit |
| scifort_multivariate_normal | multivariate_normal_logcdf |  |  | scifort_multivariate_normal_logcdf_f64 |  | multivariate_normal.logcdf |
| scifort_multivariate_normal | multivariate_normal_logpdf |  |  | scifort_multivariate_normal_logpdf_f64 |  | multivariate_normal.logpdf |
| scifort_multivariate_normal | multivariate_normal_marginal |  |  | scifort_multivariate_normal_marginal_f64 |  | multivariate_normal.marginal |
| scifort_multivariate_normal | multivariate_normal_pdf |  |  | scifort_multivariate_normal_pdf_f64 |  | multivariate_normal.pdf |
| scifort_multivariate_normal | multivariate_normal_rvs |  |  |  |  |  |
| scifort_multivariate_normal | multivariate_normal_rvs_array |  |  | scifort_multivariate_normal_rvs_array_f64 |  | multivariate_normal.rvs_array |
| scifort_multivariate_normal | multivariate_status_invalid_parameter |  |  |  |  |  |
| scifort_multivariate_normal | multivariate_status_invalid_shape |  |  |  |  |  |
| scifort_multivariate_normal | multivariate_status_linalg_failure |  |  |  |  |  |
| scifort_multivariate_normal | multivariate_status_success |  |  |  |  |  |
| scifort_multivariate_t | multivariate_t_cdf |  |  | scifort_multivariate_t_cdf_f64 |  | multivariate_t.cdf |
| scifort_multivariate_t | multivariate_t_entropy |  |  | scifort_multivariate_t_entropy_f64 |  | multivariate_t.entropy |
| scifort_multivariate_t | multivariate_t_logpdf |  |  | scifort_multivariate_t_logpdf_f64 |  | multivariate_t.logpdf |
| scifort_multivariate_t | multivariate_t_marginal |  |  | scifort_multivariate_t_marginal_f64 |  | multivariate_t.marginal |
| scifort_multivariate_t | multivariate_t_pdf |  |  | scifort_multivariate_t_pdf_f64 |  | multivariate_t.pdf |
| scifort_multivariate_t | multivariate_t_rvs |  |  |  |  |  |
| scifort_multivariate_t | multivariate_t_rvs_array |  |  | scifort_multivariate_t_rvs_array_f64 |  | multivariate_t.rvs_array |
| scifort_multivariate_t | multivariate_t_status_invalid_parameter |  |  |  |  |  |
| scifort_multivariate_t | multivariate_t_status_invalid_shape |  |  |  |  |  |
| scifort_multivariate_t | multivariate_t_status_linalg_failure |  |  |  |  |  |
| scifort_multivariate_t | multivariate_t_status_success |  |  |  |  |  |
| scifort_normal | normal_cdf | scifort_normal_cdf_f64 | scifort_normal_cdf_vec_f64 |  | scifort_normal_cdf_elementwise_f64 | normal.cdf |
| scifort_normal | normal_isf | scifort_normal_isf_f64 | scifort_normal_isf_vec_f64 |  | scifort_normal_isf_elementwise_f64 | normal.isf |
| scifort_normal | normal_logcdf | scifort_normal_logcdf_f64 | scifort_normal_logcdf_vec_f64 |  | scifort_normal_logcdf_elementwise_f64 | normal.logcdf |
| scifort_normal | normal_logpdf | scifort_normal_logpdf_f64 | scifort_normal_logpdf_vec_f64 |  | scifort_normal_logpdf_elementwise_f64 | normal.logpdf |
| scifort_normal | normal_logsf | scifort_normal_logsf_f64 | scifort_normal_logsf_vec_f64 |  | scifort_normal_logsf_elementwise_f64 | normal.logsf |
| scifort_normal | normal_pdf | scifort_normal_pdf_f64 | scifort_normal_pdf_vec_f64 |  | scifort_normal_pdf_elementwise_f64 | normal.pdf |
| scifort_normal | normal_ppf | scifort_normal_ppf_f64 | scifort_normal_ppf_vec_f64 |  | scifort_normal_ppf_elementwise_f64 | normal.ppf |
| scifort_normal | normal_sf | scifort_normal_sf_f64 | scifort_normal_sf_vec_f64 |  | scifort_normal_sf_elementwise_f64 | normal.sf |
| scifort_normal_inverse_gamma | normal_inverse_gamma_logpdf |  |  | scifort_normal_inverse_gamma_logpdf_f64 |  | normal_inverse_gamma.logpdf |
| scifort_normal_inverse_gamma | normal_inverse_gamma_mean |  |  | scifort_normal_inverse_gamma_mean_f64 |  | normal_inverse_gamma.mean |
| scifort_normal_inverse_gamma | normal_inverse_gamma_pdf |  |  | scifort_normal_inverse_gamma_pdf_f64 |  | normal_inverse_gamma.pdf |
| scifort_normal_inverse_gamma | normal_inverse_gamma_rvs |  |  | scifort_normal_inverse_gamma_rvs_f64 |  | normal_inverse_gamma.rvs |
| scifort_normal_inverse_gamma | normal_inverse_gamma_rvs_array |  |  | scifort_normal_inverse_gamma_rvs_array_f64 |  | normal_inverse_gamma.rvs_array |
| scifort_normal_inverse_gamma | normal_inverse_gamma_var |  |  | scifort_normal_inverse_gamma_var_f64 |  | normal_inverse_gamma.var |
| scifort_ortho_group | ortho_group_rvs |  |  | scifort_ortho_group_rvs_f64 |  | ortho_group.rvs |
| scifort_ortho_group | ortho_group_rvs_array |  |  | scifort_ortho_group_rvs_array_f64 |  | ortho_group.rvs_array |
| scifort_ortho_group | ortho_group_status_invalid_shape |  |  |  |  |  |
| scifort_ortho_group | ortho_group_status_numerical_failure |  |  |  |  |  |
| scifort_ortho_group | ortho_group_status_success |  |  |  |  |  |
| scifort_poisson_binom | poisson_binom_cdf |  |  | scifort_poisson_binom_cdf_explicit_f64 |  | poisson_binom.cdf |
| scifort_poisson_binom | poisson_binom_isf |  |  | scifort_poisson_binom_isf_explicit_f64 |  | poisson_binom.isf |
| scifort_poisson_binom | poisson_binom_logcdf |  |  | scifort_poisson_binom_logcdf_explicit_f64 |  | poisson_binom.logcdf |
| scifort_poisson_binom | poisson_binom_logpmf |  |  | scifort_poisson_binom_logpmf_explicit_f64 |  | poisson_binom.logpmf |
| scifort_poisson_binom | poisson_binom_logpmf_derivative_p |  |  | scifort_poisson_binom_logpmf_derivative_p_explicit_f64 |  | poisson_binom.logpmf_derivative_p |
| scifort_poisson_binom | poisson_binom_logsf |  |  | scifort_poisson_binom_logsf_explicit_f64 |  | poisson_binom.logsf |
| scifort_poisson_binom | poisson_binom_pmf |  |  | scifort_poisson_binom_pmf_explicit_f64 |  | poisson_binom.pmf |
| scifort_poisson_binom | poisson_binom_ppf |  |  | scifort_poisson_binom_ppf_explicit_f64 |  | poisson_binom.ppf |
| scifort_poisson_binom | poisson_binom_sf |  |  | scifort_poisson_binom_sf_explicit_f64 |  | poisson_binom.sf |
| scifort_random_correlation | random_correlation_rvs |  |  | scifort_random_correlation_rvs_f64 |  | random_correlation.rvs |
| scifort_random_correlation | random_correlation_status_invalid_parameter |  |  |  |  |  |
| scifort_random_correlation | random_correlation_status_invalid_shape |  |  |  |  |  |
| scifort_random_correlation | random_correlation_status_numerical_failure |  |  |  |  |  |
| scifort_random_correlation | random_correlation_status_success |  |  |  |  |  |
| scifort_random_table | random_table_logpmf |  |  | scifort_random_table_logpmf_f64 |  | random_table.logpmf |
| scifort_random_table | random_table_mean |  |  | scifort_random_table_mean_f64 |  | random_table.mean |
| scifort_random_table | random_table_pmf |  |  | scifort_random_table_pmf_f64 |  | random_table.pmf |
| scifort_random_table | random_table_rvs |  |  | scifort_random_table_rvs_f64 |  | random_table.rvs |
| scifort_random_table | random_table_rvs_array |  |  | scifort_random_table_rvs_array_f64 |  | random_table.rvs_array |
| scifort_random_table | random_table_status_invalid_parameter |  |  |  |  |  |
| scifort_random_table | random_table_status_invalid_shape |  |  |  |  |  |
| scifort_random_table | random_table_status_numerical_failure |  |  |  |  |  |
| scifort_random_table | random_table_status_success |  |  |  |  |  |
| scifort_special_ortho_group | special_ortho_group_rvs |  |  | scifort_special_ortho_group_rvs_f64 |  | special_ortho_group.rvs |
| scifort_special_ortho_group | special_ortho_group_rvs_array |  |  | scifort_special_ortho_group_rvs_array_f64 |  | special_ortho_group.rvs_array |
| scifort_special_ortho_group | special_ortho_group_status_invalid_shape |  |  |  |  |  |
| scifort_special_ortho_group | special_ortho_group_status_numerical_failure |  |  |  |  |  |
| scifort_special_ortho_group | special_ortho_group_status_success |  |  |  |  |  |
| scifort_student_t | t_cdf | scifort_t_cdf_f64 | scifort_t_cdf_vec_f64 |  | scifort_t_cdf_elementwise_f64 | t.cdf |
| scifort_student_t | t_isf | scifort_t_isf_f64 | scifort_t_isf_vec_f64 |  | scifort_t_isf_elementwise_f64 | t.isf |
| scifort_student_t | t_logcdf | scifort_t_logcdf_f64 | scifort_t_logcdf_vec_f64 |  | scifort_t_logcdf_elementwise_f64 | t.logcdf |
| scifort_student_t | t_logpdf | scifort_t_logpdf_f64 | scifort_t_logpdf_vec_f64 |  | scifort_t_logpdf_elementwise_f64 | t.logpdf |
| scifort_student_t | t_logsf | scifort_t_logsf_f64 | scifort_t_logsf_vec_f64 |  | scifort_t_logsf_elementwise_f64 | t.logsf |
| scifort_student_t | t_pdf | scifort_t_pdf_f64 | scifort_t_pdf_vec_f64 |  | scifort_t_pdf_elementwise_f64 | t.pdf |
| scifort_student_t | t_ppf | scifort_t_ppf_f64 | scifort_t_ppf_vec_f64 |  | scifort_t_ppf_elementwise_f64 | t.ppf |
| scifort_student_t | t_sf | scifort_t_sf_f64 | scifort_t_sf_vec_f64 |  | scifort_t_sf_elementwise_f64 | t.sf |
| scifort_uniform_direction | uniform_direction_rvs |  |  | scifort_uniform_direction_rvs_f64 |  | uniform_direction.rvs |
| scifort_uniform_direction | uniform_direction_rvs_array |  |  | scifort_uniform_direction_rvs_array_f64 |  | uniform_direction.rvs_array |
| scifort_uniform_direction | uniform_direction_status_invalid_shape |  |  |  |  |  |
| scifort_uniform_direction | uniform_direction_status_numerical_failure |  |  |  |  |  |
| scifort_uniform_direction | uniform_direction_status_success |  |  |  |  |  |
| scifort_unitary_group | unitary_group_rvs |  |  | scifort_unitary_group_rvs_f64 |  | unitary_group.rvs |
| scifort_unitary_group | unitary_group_rvs_array |  |  | scifort_unitary_group_rvs_array_f64 |  | unitary_group.rvs_array |
| scifort_unitary_group | unitary_group_status_invalid_shape |  |  |  |  |  |
| scifort_unitary_group | unitary_group_status_numerical_failure |  |  |  |  |  |
| scifort_unitary_group | unitary_group_status_success |  |  |  |  |  |
| scifort_vonmises_fisher | vonmises_fisher_entropy |  |  | scifort_vonmises_fisher_entropy_f64 |  | vonmises_fisher.entropy |
| scifort_vonmises_fisher | vonmises_fisher_fit |  |  | scifort_vonmises_fisher_fit_f64 |  | vonmises_fisher.fit |
| scifort_vonmises_fisher | vonmises_fisher_logpdf |  |  | scifort_vonmises_fisher_logpdf_f64 |  | vonmises_fisher.logpdf |
| scifort_vonmises_fisher | vonmises_fisher_pdf |  |  | scifort_vonmises_fisher_pdf_f64 |  | vonmises_fisher.pdf |
| scifort_vonmises_fisher | vonmises_fisher_rvs |  |  | scifort_vonmises_fisher_rvs_f64 |  | vonmises_fisher.rvs |
| scifort_vonmises_fisher | vonmises_fisher_rvs_array |  |  | scifort_vonmises_fisher_rvs_array_f64 |  | vonmises_fisher.rvs_array |
| scifort_vonmises_fisher | vonmises_fisher_status_invalid_parameter |  |  |  |  |  |
| scifort_vonmises_fisher | vonmises_fisher_status_invalid_shape |  |  |  |  |  |
| scifort_vonmises_fisher | vonmises_fisher_status_numerical_failure |  |  |  |  |  |
| scifort_vonmises_fisher | vonmises_fisher_status_success |  |  |  |  |  |
| scifort_weibull | weibull_cdf | scifort_weibull_cdf_f64 | scifort_weibull_cdf_vec_f64 |  | scifort_weibull_cdf_elementwise_f64 | weibull.cdf |
| scifort_weibull | weibull_isf | scifort_weibull_isf_f64 | scifort_weibull_isf_vec_f64 |  | scifort_weibull_isf_elementwise_f64 | weibull.isf |
| scifort_weibull | weibull_logcdf | scifort_weibull_logcdf_f64 | scifort_weibull_logcdf_vec_f64 |  | scifort_weibull_logcdf_elementwise_f64 | weibull.logcdf |
| scifort_weibull | weibull_logpdf | scifort_weibull_logpdf_f64 | scifort_weibull_logpdf_vec_f64 |  | scifort_weibull_logpdf_elementwise_f64 | weibull.logpdf |
| scifort_weibull | weibull_logsf | scifort_weibull_logsf_f64 | scifort_weibull_logsf_vec_f64 |  | scifort_weibull_logsf_elementwise_f64 | weibull.logsf |
| scifort_weibull | weibull_pdf | scifort_weibull_pdf_f64 | scifort_weibull_pdf_vec_f64 |  | scifort_weibull_pdf_elementwise_f64 | weibull.pdf |
| scifort_weibull | weibull_ppf | scifort_weibull_ppf_f64 | scifort_weibull_ppf_vec_f64 |  | scifort_weibull_ppf_elementwise_f64 | weibull.ppf |
| scifort_weibull | weibull_sf | scifort_weibull_sf_f64 | scifort_weibull_sf_vec_f64 |  | scifort_weibull_sf_elementwise_f64 | weibull.sf |
| scifort_wishart | wishart_entropy |  |  | scifort_wishart_entropy_f64 |  | wishart.entropy |
| scifort_wishart | wishart_logpdf |  |  | scifort_wishart_logpdf_f64 |  | wishart.logpdf |
| scifort_wishart | wishart_mean |  |  | scifort_wishart_mean_f64 |  | wishart.mean |
| scifort_wishart | wishart_mode |  |  | scifort_wishart_mode_f64 |  | wishart.mode |
| scifort_wishart | wishart_pdf |  |  | scifort_wishart_pdf_f64 |  | wishart.pdf |
| scifort_wishart | wishart_rvs |  |  |  |  |  |
| scifort_wishart | wishart_rvs_array |  |  | scifort_wishart_rvs_array_f64 |  | wishart.rvs_array |
| scifort_wishart | wishart_status_invalid_parameter |  |  |  |  |  |
| scifort_wishart | wishart_status_invalid_shape |  |  |  |  |  |
| scifort_wishart | wishart_status_linalg_failure |  |  |  |  |  |
| scifort_wishart | wishart_status_success |  |  |  |  |  |
| scifort_wishart | wishart_var |  |  | scifort_wishart_var_f64 |  | wishart.var |

C-manifest exclusion reasons and generated/existing symbol distinctions are in the JSON.
