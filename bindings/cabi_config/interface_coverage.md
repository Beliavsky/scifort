# Reviewed interface coverage

Reviewed C batch and Python policy only; not a whole-project inventory.

| Interface | Count |
| --- | ---: |
| exports | 98 |
| scalar | 88 |
| bulk | 88 |
| array_adapters | 1 |
| elementwise | 16 |
| python | 88 |

| Module | Export | Scalar C | Bulk C | Array adapter | Elementwise C | Python |
| --- | --- | --- | --- | --- | --- | --- |
| scifort_beta | beta_cdf | scifort_beta_cdf_f64 | scifort_beta_cdf_vec_f64 |  |  | beta.cdf |
| scifort_beta | beta_isf | scifort_beta_isf_f64 | scifort_beta_isf_vec_f64 |  |  | beta.isf |
| scifort_beta | beta_logcdf | scifort_beta_logcdf_f64 | scifort_beta_logcdf_vec_f64 |  |  | beta.logcdf |
| scifort_beta | beta_logpdf | scifort_beta_logpdf_f64 | scifort_beta_logpdf_vec_f64 |  |  | beta.logpdf |
| scifort_beta | beta_logsf | scifort_beta_logsf_f64 | scifort_beta_logsf_vec_f64 |  |  | beta.logsf |
| scifort_beta | beta_pdf | scifort_beta_pdf_f64 | scifort_beta_pdf_vec_f64 |  |  | beta.pdf |
| scifort_beta | beta_ppf | scifort_beta_ppf_f64 | scifort_beta_ppf_vec_f64 |  |  | beta.ppf |
| scifort_beta | beta_sf | scifort_beta_sf_f64 | scifort_beta_sf_vec_f64 |  |  | beta.sf |
| scifort_cauchy | cauchy_cdf | scifort_cauchy_cdf_f64 | scifort_cauchy_cdf_vec_f64 |  |  | cauchy.cdf |
| scifort_cauchy | cauchy_isf | scifort_cauchy_isf_f64 | scifort_cauchy_isf_vec_f64 |  |  | cauchy.isf |
| scifort_cauchy | cauchy_logcdf | scifort_cauchy_logcdf_f64 | scifort_cauchy_logcdf_vec_f64 |  |  | cauchy.logcdf |
| scifort_cauchy | cauchy_logpdf | scifort_cauchy_logpdf_f64 | scifort_cauchy_logpdf_vec_f64 |  |  | cauchy.logpdf |
| scifort_cauchy | cauchy_logsf | scifort_cauchy_logsf_f64 | scifort_cauchy_logsf_vec_f64 |  |  | cauchy.logsf |
| scifort_cauchy | cauchy_pdf | scifort_cauchy_pdf_f64 | scifort_cauchy_pdf_vec_f64 |  |  | cauchy.pdf |
| scifort_cauchy | cauchy_ppf | scifort_cauchy_ppf_f64 | scifort_cauchy_ppf_vec_f64 |  |  | cauchy.ppf |
| scifort_cauchy | cauchy_sf | scifort_cauchy_sf_f64 | scifort_cauchy_sf_vec_f64 |  |  | cauchy.sf |
| scifort_chi2 | chi2_cdf | scifort_chi2_cdf_f64 | scifort_chi2_cdf_vec_f64 |  |  | chi2.cdf |
| scifort_chi2 | chi2_isf | scifort_chi2_isf_f64 | scifort_chi2_isf_vec_f64 |  |  | chi2.isf |
| scifort_chi2 | chi2_logcdf | scifort_chi2_logcdf_f64 | scifort_chi2_logcdf_vec_f64 |  |  | chi2.logcdf |
| scifort_chi2 | chi2_logpdf | scifort_chi2_logpdf_f64 | scifort_chi2_logpdf_vec_f64 |  |  | chi2.logpdf |
| scifort_chi2 | chi2_logsf | scifort_chi2_logsf_f64 | scifort_chi2_logsf_vec_f64 |  |  | chi2.logsf |
| scifort_chi2 | chi2_pdf | scifort_chi2_pdf_f64 | scifort_chi2_pdf_vec_f64 |  |  | chi2.pdf |
| scifort_chi2 | chi2_ppf | scifort_chi2_ppf_f64 | scifort_chi2_ppf_vec_f64 |  |  | chi2.ppf |
| scifort_chi2 | chi2_sf | scifort_chi2_sf_f64 | scifort_chi2_sf_vec_f64 |  |  | chi2.sf |
| scifort_exponential | exponential_cdf | scifort_exponential_cdf_f64 | scifort_exponential_cdf_vec_f64 |  |  | exponential.cdf |
| scifort_exponential | exponential_isf | scifort_exponential_isf_f64 | scifort_exponential_isf_vec_f64 |  |  | exponential.isf |
| scifort_exponential | exponential_logcdf | scifort_exponential_logcdf_f64 | scifort_exponential_logcdf_vec_f64 |  |  | exponential.logcdf |
| scifort_exponential | exponential_logpdf | scifort_exponential_logpdf_f64 | scifort_exponential_logpdf_vec_f64 |  |  | exponential.logpdf |
| scifort_exponential | exponential_logsf | scifort_exponential_logsf_f64 | scifort_exponential_logsf_vec_f64 |  |  | exponential.logsf |
| scifort_exponential | exponential_pdf | scifort_exponential_pdf_f64 | scifort_exponential_pdf_vec_f64 |  |  | exponential.pdf |
| scifort_exponential | exponential_ppf | scifort_exponential_ppf_f64 | scifort_exponential_ppf_vec_f64 |  |  | exponential.ppf |
| scifort_exponential | exponential_sf | scifort_exponential_sf_f64 | scifort_exponential_sf_vec_f64 |  |  | exponential.sf |
| scifort_gamma | gamma_cdf | scifort_gamma_cdf_f64 | scifort_gamma_cdf_vec_f64 |  | scifort_gamma_cdf_elementwise_f64 | gamma.cdf |
| scifort_gamma | gamma_isf | scifort_gamma_isf_f64 | scifort_gamma_isf_vec_f64 |  | scifort_gamma_isf_elementwise_f64 | gamma.isf |
| scifort_gamma | gamma_logcdf | scifort_gamma_logcdf_f64 | scifort_gamma_logcdf_vec_f64 |  | scifort_gamma_logcdf_elementwise_f64 | gamma.logcdf |
| scifort_gamma | gamma_logpdf | scifort_gamma_logpdf_f64 | scifort_gamma_logpdf_vec_f64 |  | scifort_gamma_logpdf_elementwise_f64 | gamma.logpdf |
| scifort_gamma | gamma_logsf | scifort_gamma_logsf_f64 | scifort_gamma_logsf_vec_f64 |  | scifort_gamma_logsf_elementwise_f64 | gamma.logsf |
| scifort_gamma | gamma_pdf | scifort_gamma_pdf_f64 | scifort_gamma_pdf_vec_f64 |  | scifort_gamma_pdf_elementwise_f64 | gamma.pdf |
| scifort_gamma | gamma_ppf | scifort_gamma_ppf_f64 | scifort_gamma_ppf_vec_f64 |  | scifort_gamma_ppf_elementwise_f64 | gamma.ppf |
| scifort_gamma | gamma_sf | scifort_gamma_sf_f64 | scifort_gamma_sf_vec_f64 |  | scifort_gamma_sf_elementwise_f64 | gamma.sf |
| scifort_laplace | laplace_cdf | scifort_laplace_cdf_f64 | scifort_laplace_cdf_vec_f64 |  |  | laplace.cdf |
| scifort_laplace | laplace_isf | scifort_laplace_isf_f64 | scifort_laplace_isf_vec_f64 |  |  | laplace.isf |
| scifort_laplace | laplace_logcdf | scifort_laplace_logcdf_f64 | scifort_laplace_logcdf_vec_f64 |  |  | laplace.logcdf |
| scifort_laplace | laplace_logpdf | scifort_laplace_logpdf_f64 | scifort_laplace_logpdf_vec_f64 |  |  | laplace.logpdf |
| scifort_laplace | laplace_logsf | scifort_laplace_logsf_f64 | scifort_laplace_logsf_vec_f64 |  |  | laplace.logsf |
| scifort_laplace | laplace_pdf | scifort_laplace_pdf_f64 | scifort_laplace_pdf_vec_f64 |  |  | laplace.pdf |
| scifort_laplace | laplace_ppf | scifort_laplace_ppf_f64 | scifort_laplace_ppf_vec_f64 |  |  | laplace.ppf |
| scifort_laplace | laplace_sf | scifort_laplace_sf_f64 | scifort_laplace_sf_vec_f64 |  |  | laplace.sf |
| scifort_linalg | cholesky_lower |  |  |  |  |  |
| scifort_linalg | linalg_status_invalid_shape |  |  |  |  |  |
| scifort_linalg | linalg_status_no_convergence |  |  |  |  |  |
| scifort_linalg | linalg_status_not_positive_definite |  |  |  |  |  |
| scifort_linalg | linalg_status_not_positive_semidefinite |  |  |  |  |  |
| scifort_linalg | linalg_status_singular |  |  |  |  |  |
| scifort_linalg | linalg_status_success |  |  |  |  |  |
| scifort_linalg | psd_decompose |  |  |  |  |  |
| scifort_linalg | solve_lower |  |  | scifort_solve_lower_f64 |  |  |
| scifort_linalg | symmetric_eigen_jacobi |  |  |  |  |  |
| scifort_logistic | logistic_cdf | scifort_logistic_cdf_f64 | scifort_logistic_cdf_vec_f64 |  |  | logistic.cdf |
| scifort_logistic | logistic_isf | scifort_logistic_isf_f64 | scifort_logistic_isf_vec_f64 |  |  | logistic.isf |
| scifort_logistic | logistic_logcdf | scifort_logistic_logcdf_f64 | scifort_logistic_logcdf_vec_f64 |  |  | logistic.logcdf |
| scifort_logistic | logistic_logpdf | scifort_logistic_logpdf_f64 | scifort_logistic_logpdf_vec_f64 |  |  | logistic.logpdf |
| scifort_logistic | logistic_logsf | scifort_logistic_logsf_f64 | scifort_logistic_logsf_vec_f64 |  |  | logistic.logsf |
| scifort_logistic | logistic_pdf | scifort_logistic_pdf_f64 | scifort_logistic_pdf_vec_f64 |  |  | logistic.pdf |
| scifort_logistic | logistic_ppf | scifort_logistic_ppf_f64 | scifort_logistic_ppf_vec_f64 |  |  | logistic.ppf |
| scifort_logistic | logistic_sf | scifort_logistic_sf_f64 | scifort_logistic_sf_vec_f64 |  |  | logistic.sf |
| scifort_lognormal | lognormal_cdf | scifort_lognormal_cdf_f64 | scifort_lognormal_cdf_vec_f64 |  |  | lognormal.cdf |
| scifort_lognormal | lognormal_isf | scifort_lognormal_isf_f64 | scifort_lognormal_isf_vec_f64 |  |  | lognormal.isf |
| scifort_lognormal | lognormal_logcdf | scifort_lognormal_logcdf_f64 | scifort_lognormal_logcdf_vec_f64 |  |  | lognormal.logcdf |
| scifort_lognormal | lognormal_logpdf | scifort_lognormal_logpdf_f64 | scifort_lognormal_logpdf_vec_f64 |  |  | lognormal.logpdf |
| scifort_lognormal | lognormal_logsf | scifort_lognormal_logsf_f64 | scifort_lognormal_logsf_vec_f64 |  |  | lognormal.logsf |
| scifort_lognormal | lognormal_pdf | scifort_lognormal_pdf_f64 | scifort_lognormal_pdf_vec_f64 |  |  | lognormal.pdf |
| scifort_lognormal | lognormal_ppf | scifort_lognormal_ppf_f64 | scifort_lognormal_ppf_vec_f64 |  |  | lognormal.ppf |
| scifort_lognormal | lognormal_sf | scifort_lognormal_sf_f64 | scifort_lognormal_sf_vec_f64 |  |  | lognormal.sf |
| scifort_normal | normal_cdf | scifort_normal_cdf_f64 | scifort_normal_cdf_vec_f64 |  | scifort_normal_cdf_elementwise_f64 | normal.cdf |
| scifort_normal | normal_isf | scifort_normal_isf_f64 | scifort_normal_isf_vec_f64 |  | scifort_normal_isf_elementwise_f64 | normal.isf |
| scifort_normal | normal_logcdf | scifort_normal_logcdf_f64 | scifort_normal_logcdf_vec_f64 |  | scifort_normal_logcdf_elementwise_f64 | normal.logcdf |
| scifort_normal | normal_logpdf | scifort_normal_logpdf_f64 | scifort_normal_logpdf_vec_f64 |  | scifort_normal_logpdf_elementwise_f64 | normal.logpdf |
| scifort_normal | normal_logsf | scifort_normal_logsf_f64 | scifort_normal_logsf_vec_f64 |  | scifort_normal_logsf_elementwise_f64 | normal.logsf |
| scifort_normal | normal_pdf | scifort_normal_pdf_f64 | scifort_normal_pdf_vec_f64 |  | scifort_normal_pdf_elementwise_f64 | normal.pdf |
| scifort_normal | normal_ppf | scifort_normal_ppf_f64 | scifort_normal_ppf_vec_f64 |  | scifort_normal_ppf_elementwise_f64 | normal.ppf |
| scifort_normal | normal_sf | scifort_normal_sf_f64 | scifort_normal_sf_vec_f64 |  | scifort_normal_sf_elementwise_f64 | normal.sf |
| scifort_student_t | t_cdf | scifort_t_cdf_f64 | scifort_t_cdf_vec_f64 |  |  | t.cdf |
| scifort_student_t | t_isf | scifort_t_isf_f64 | scifort_t_isf_vec_f64 |  |  | t.isf |
| scifort_student_t | t_logcdf | scifort_t_logcdf_f64 | scifort_t_logcdf_vec_f64 |  |  | t.logcdf |
| scifort_student_t | t_logpdf | scifort_t_logpdf_f64 | scifort_t_logpdf_vec_f64 |  |  | t.logpdf |
| scifort_student_t | t_logsf | scifort_t_logsf_f64 | scifort_t_logsf_vec_f64 |  |  | t.logsf |
| scifort_student_t | t_pdf | scifort_t_pdf_f64 | scifort_t_pdf_vec_f64 |  |  | t.pdf |
| scifort_student_t | t_ppf | scifort_t_ppf_f64 | scifort_t_ppf_vec_f64 |  |  | t.ppf |
| scifort_student_t | t_sf | scifort_t_sf_f64 | scifort_t_sf_vec_f64 |  |  | t.sf |
| scifort_weibull | weibull_cdf | scifort_weibull_cdf_f64 | scifort_weibull_cdf_vec_f64 |  |  | weibull.cdf |
| scifort_weibull | weibull_isf | scifort_weibull_isf_f64 | scifort_weibull_isf_vec_f64 |  |  | weibull.isf |
| scifort_weibull | weibull_logcdf | scifort_weibull_logcdf_f64 | scifort_weibull_logcdf_vec_f64 |  |  | weibull.logcdf |
| scifort_weibull | weibull_logpdf | scifort_weibull_logpdf_f64 | scifort_weibull_logpdf_vec_f64 |  |  | weibull.logpdf |
| scifort_weibull | weibull_logsf | scifort_weibull_logsf_f64 | scifort_weibull_logsf_vec_f64 |  |  | weibull.logsf |
| scifort_weibull | weibull_pdf | scifort_weibull_pdf_f64 | scifort_weibull_pdf_vec_f64 |  |  | weibull.pdf |
| scifort_weibull | weibull_ppf | scifort_weibull_ppf_f64 | scifort_weibull_ppf_vec_f64 |  |  | weibull.ppf |
| scifort_weibull | weibull_sf | scifort_weibull_sf_f64 | scifort_weibull_sf_vec_f64 |  |  | weibull.sf |

C-manifest exclusion reasons and generated/existing symbol distinctions are in the JSON.
