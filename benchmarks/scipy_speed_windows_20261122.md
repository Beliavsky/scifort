# SciFort versus SciPy: local Windows results, 2026-11-22

This report refresh uses the latest saved 80-case run, including the direct
inverse-normal kernel. It does not represent a new measurement.

This expanded 80-case snapshot was collected while unrelated CPU-intensive
Python validation jobs were running. Treat the times as an illustration of the
harness and its outputs, rather than a controlled comparison of library speeds.
Rerun on an idle machine before drawing performance conclusions.

See the [workload guide](README_scipy_speed.md) for matching inputs, setup costs,
API differences, allocation policy, and numerical tolerances. The older
[nine-case snapshot](scipy_speed_windows_20261122_initial.md) is preserved separately.

Run: 2026-11-23T02:49:00.172246+00:00

Platform: Windows-11-10.0.26300-SP0

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental)

Flags: `-O3 -std=f2018 -fimplicit-none -Werror=implicit-interface`

Python 3.13.3; NumPy 2.2.6; SciPy 1.15.3.

Quick workloads: False; validation only: False.

Times are milliseconds per full operation; ratios above one favor SciFort.
These are application workloads, including documented allocation and API overhead.

| Case | n | SciFort ms | SciPy ms | SciPy/SciFort | Max abs error |
| --- | ---: | ---: | ---: | ---: | ---: |
| normal_pdf_array | 16384 | 0.902784 | 2.55069 | 2.825 | 5.55e-17 |
| normal_logpdf_array | 16384 | 0.588487 | 2.28296 | 3.879 | 8.88e-16 |
| normal_cdf_array | 65536 | 5.39341 | 12.2543 | 2.272 | 2.22e-16 |
| normal_sf_array | 16384 | 1.27219 | 4.92088 | 3.868 | 2.22e-16 |
| normal_logcdf_array | 16384 | 1.89421 | 6.40479 | 3.381 | 1.42e-14 |
| normal_logsf_array | 16384 | 1.76597 | 5.50072 | 3.115 | 1.42e-14 |
| normal_ppf_array | 16384 | 0.642255 | 1.72752 | 2.690 | 4.44e-16 |
| normal_isf_array | 16384 | 0.714673 | 3.80735 | 5.327 | 4.44e-16 |
| exponential_pdf_array | 16384 | 0.449441 | 1.66834 | 3.712 | 0 |
| exponential_logpdf_array | 16384 | 0.319548 | 1.55781 | 4.875 | 0 |
| exponential_cdf_array | 16384 | 0.743862 | 2.97736 | 4.003 | 1.11e-16 |
| exponential_sf_array | 16384 | 0.583421 | 2.01314 | 3.451 | 0 |
| exponential_logcdf_array | 16384 | 1.19753 | 3.30392 | 2.759 | 8.88e-16 |
| exponential_logsf_array | 16384 | 0.40224 | 0.940834 | 2.339 | 0 |
| exponential_ppf_array | 16384 | 1.13216 | 2.77838 | 2.454 | 8.88e-16 |
| exponential_isf_array | 16384 | 0.675495 | 1.65204 | 2.446 | 0 |
| logistic_pdf_array | 16384 | 2.02512 | 2.12223 | 1.048 | 2.22e-16 |
| logistic_logpdf_array | 16384 | 1.36922 | 1.40007 | 1.023 | 8.88e-16 |
| logistic_cdf_array | 16384 | 0.827623 | 3.09067 | 3.734 | 2.22e-16 |
| logistic_sf_array | 16384 | 0.29917 | 1.36375 | 4.558 | 2.22e-16 |
| logistic_logcdf_array | 16384 | 0.947466 | 2.11678 | 2.234 | 8.88e-16 |
| logistic_logsf_array | 16384 | 1.6654 | 3.94615 | 2.369 | 8.88e-16 |
| logistic_ppf_array | 16384 | 1.26882 | 2.15409 | 1.698 | 1.78e-15 |
| logistic_isf_array | 16384 | 1.08959 | 1.88849 | 1.733 | 1.78e-15 |
| laplace_pdf_array | 16384 | 0.522304 | 0.999883 | 1.914 | 0 |
| laplace_logpdf_array | 16384 | 0.584135 | 2.07908 | 3.559 | 2.22e-16 |
| laplace_cdf_array | 16384 | 0.369229 | 1.11565 | 3.022 | 0 |
| laplace_sf_array | 16384 | 0.343766 | 1.2449 | 3.621 | 0 |
| laplace_logcdf_array | 16384 | 0.772182 | 2.27223 | 2.943 | 4.44e-16 |
| laplace_logsf_array | 16384 | 0.945575 | 2.80839 | 2.970 | 4.44e-16 |
| laplace_ppf_array | 16384 | 0.69958 | 1.6003 | 2.288 | 1.78e-15 |
| laplace_isf_array | 16384 | 0.496566 | 1.3796 | 2.778 | 1.78e-15 |
| cauchy_pdf_array | 16384 | 1.29581 | 0.750668 | 0.579 | 5.55e-17 |
| cauchy_logpdf_array | 16384 | 0.869065 | 1.93251 | 2.224 | 8.88e-16 |
| cauchy_cdf_array | 16384 | 0.414341 | 1.1078 | 2.674 | 2.22e-16 |
| cauchy_sf_array | 16384 | 0.726789 | 1.61082 | 2.216 | 2.22e-16 |
| cauchy_logcdf_array | 16384 | 1.73203 | 2.64315 | 1.526 | 4.44e-16 |
| cauchy_logsf_array | 16384 | 1.62898 | 2.39016 | 1.467 | 4.44e-16 |
| cauchy_ppf_array | 16384 | 0.802848 | 1.90537 | 2.373 | 0 |
| cauchy_isf_array | 16384 | 0.927977 | 2.14903 | 2.316 | 0 |
| weibull_pdf_array | 16384 | 3.28642 | 3.80041 | 1.156 | 2.22e-16 |
| weibull_logpdf_array | 16384 | 2.17663 | 4.33732 | 1.993 | 2.84e-14 |
| weibull_cdf_array | 16384 | 1.7794 | 3.00074 | 1.686 | 2.22e-16 |
| weibull_sf_array | 16384 | 1.57089 | 1.93581 | 1.232 | 1.11e-16 |
| weibull_logcdf_array | 16384 | 2.74468 | 4.59694 | 1.675 | 1.78e-15 |
| weibull_logsf_array | 16384 | 1.34005 | 1.94742 | 1.453 | 2.13e-14 |
| weibull_ppf_array | 16384 | 2.29037 | 3.45068 | 1.507 | 4.44e-16 |
| weibull_isf_array | 16384 | 1.97525 | 2.72135 | 1.378 | 8.88e-16 |
| lognormal_pdf_array | 16384 | 1.49308 | 2.0689 | 1.386 | 4.44e-16 |
| lognormal_logpdf_array | 16384 | 1.63196 | 2.00964 | 1.231 | 7.11e-15 |
| lognormal_cdf_array | 16384 | 1.97854 | 3.12774 | 1.581 | 2.22e-16 |
| lognormal_sf_array | 16384 | 2.30186 | 3.02601 | 1.315 | 2.22e-16 |
| lognormal_logcdf_array | 16384 | 3.62845 | 4.9262 | 1.358 | 7.11e-15 |
| lognormal_logsf_array | 16384 | 3.54745 | 4.03594 | 1.138 | 2.66e-15 |
| lognormal_ppf_array | 16384 | 2.32697 | 2.68624 | 1.154 | 2.66e-15 |
| lognormal_isf_array | 16384 | 1.35756 | 1.6074 | 1.184 | 8.88e-16 |
| normal_cdf_scalar | 256 | 0.0126103 | 0.54215 | 42.993 | 2.22e-16 |
| fft_power2_reuse | 4096 | 0.230535 | 0.127747 | 0.554 | 7.64e-11 |
| fft_power2_setup | 4096 | 0.229016 | 0.116185 | 0.507 | 7.64e-11 |
| fft_mixed_reuse | 4050 | 0.643141 | 0.167549 | 0.261 | 2.59e-12 |
| fft_mixed_setup | 4050 | 1.1864 | 0.235628 | 0.199 | 2.59e-12 |
| fft_prime_reuse | 4093 | 1.24731 | 0.688341 | 0.552 | 9.03e-11 |
| fft_prime_setup | 4093 | 2.58242 | 0.594559 | 0.230 | 9.03e-11 |
| ifft_power2_reuse | 4096 | 0.259999 | 0.154213 | 0.593 | 2.3e-14 |
| ifft_power2_setup | 4096 | 0.263759 | 0.143746 | 0.545 | 2.3e-14 |
| ifft_mixed_reuse | 4050 | 0.652756 | 0.182516 | 0.280 | 4.16e-16 |
| ifft_mixed_setup | 4050 | 1.00535 | 0.170707 | 0.170 | 4.16e-16 |
| ifft_prime_reuse | 4093 | 0.922819 | 0.463888 | 0.503 | 3.16e-14 |
| ifft_prime_setup | 4093 | 2.26984 | 0.478229 | 0.211 | 3.16e-14 |
| lu_factor_solve | 128 | 0.946837 | 0.427328 | 0.451 | 2.22e-15 |
| lu_reuse | 128 | 0.0716239 | 0.0321605 | 0.449 | 2.22e-15 |
| lu_transpose | 128 | 0.061519 | 0.0307176 | 0.499 | 1.89e-15 |
| cho_factor_solve | 128 | 0.858967 | 0.19228 | 0.224 | 1.78e-15 |
| cho_reuse | 128 | 0.114183 | 0.0735563 | 0.644 | 1.78e-15 |
| triangular_lower | 128 | 0.0474263 | 0.050998 | 1.075 | 1.33e-15 |
| triangular_upper | 128 | 0.0465954 | 0.0496637 | 1.066 | 1.33e-15 |
| triangular_transpose | 128 | 0.0453473 | 0.0668551 | 1.474 | 1.33e-15 |
| barycentric | 4096 | 2.72147 | 2.90065 | 1.066 | 2e-15 |
| barycentric_setup | 4096 | 3.04907 | 4.14886 | 1.361 | 2e-15 |
| barycentric_derivative | 4096 | 2.81779 | 2.79896 | 0.993 | 4.84e-14 |

## Individual samples

Milliseconds per operation, in sampling order. Full metadata and checks are in the sibling JSON.

- `normal_pdf_array`: SciFort 0.829242, 0.902784, 0.979764, 0.828584, 0.804927, 1.16411, 1.21961; SciPy 2.46582, 2.72718, 2.22896, 2.86604, 2.95494, 2.54109, 2.55069.
- `normal_logpdf_array`: SciFort 0.588487, 0.596417, 0.681132, 0.543333, 0.592264, 0.441128, 0.536173; SciPy 2.547, 2.90137, 2.6622, 2.28296, 2.25135, 1.95394, 2.19885.
- `normal_cdf_array`: SciFort 4.57847, 5.08024, 5.86638, 5.33406, 5.39341, 5.997, 5.89758; SciPy 10.766, 12.2543, 12.6174, 12.5412, 11.8114, 12.1565, 12.7379.
- `normal_sf_array`: SciFort 1.54916, 1.36058, 1.37022, 1.08259, 1.08154, 1.27219, 1.13662; SciPy 5.1905, 4.96157, 5.06012, 4.00896, 4.3905, 3.44144, 4.92088.
- `normal_logcdf_array`: SciFort 2.53359, 2.00605, 2.12488, 1.89421, 1.83428, 1.81978, 1.59804; SciPy 7.48492, 6.56547, 6.40479, 6.63383, 6.12483, 5.31662, 5.93718.
- `normal_logsf_array`: SciFort 1.74171, 1.73921, 1.76597, 1.58332, 1.80601, 2.10501, 1.86004; SciPy 5.24829, 5.50072, 6.32924, 5.34171, 5.1296, 5.6538, 5.57431.
- `normal_ppf_array`: SciFort 0.642255, 0.580948, 0.604079, 0.653438, 0.546386, 0.679459, 0.761316; SciPy 1.91675, 1.72752, 1.87505, 1.56418, 1.43166, 1.61673, 2.15937.
- `normal_isf_array`: SciFort 1.32756, 1.04163, 0.888691, 0.607035, 0.605771, 0.714673, 0.543416; SciPy 4.957, 3.9228, 3.80735, 3.93419, 2.82417, 3.66319, 2.79707.
- `exponential_pdf_array`: SciFort 0.426121, 0.453604, 0.449441, 0.468776, 0.431308, 0.536018, 0.417769; SciPy 1.72144, 1.76018, 1.66834, 1.62962, 1.71123, 1.59843, 1.51388.
- `exponential_logpdf_array`: SciFort 0.359832, 0.303319, 0.266827, 0.265801, 0.339275, 0.319548, 0.659213; SciPy 1.55781, 1.6976, 1.4362, 1.39845, 1.52086, 1.62648, 1.70357.
- `exponential_cdf_array`: SciFort 0.899486, 0.641438, 0.580435, 0.576303, 0.743862, 0.962605, 0.764458; SciPy 2.97736, 2.51218, 2.65301, 2.75583, 3.96301, 4.47289, 4.43478.
- `exponential_sf_array`: SciFort 0.511257, 0.575777, 0.609123, 0.50491, 0.583421, 0.812525, 0.92825; SciPy 2.01314, 2.03199, 2.24406, 1.92985, 1.87666, 2.12041, 1.80244.
- `exponential_logcdf_array`: SciFort 1.41538, 1.19983, 1.17003, 1.0804, 1.1722, 1.19753, 2.09185; SciPy 3.32598, 3.16388, 3.86998, 3.30392, 3.82036, 3.06173, 3.10242.
- `exponential_logsf_array`: SciFort 0.385102, 0.40224, 0.40477, 0.513687, 0.407696, 0.292861, 0.321643; SciPy 0.920064, 0.932691, 0.937666, 1.20918, 1.11699, 1.02443, 0.940834.
- `exponential_ppf_array`: SciFort 1.07259, 1.11296, 1.13216, 1.58433, 1.44537, 1.10793, 1.28464; SciPy 2.64931, 2.77838, 3.2742, 4.31634, 2.93273, 2.32237, 2.51117.
- `exponential_isf_array`: SciFort 0.672379, 0.611539, 0.652038, 0.708073, 0.741348, 0.675495, 0.773222; SciPy 1.54195, 1.64975, 1.81379, 2.00836, 1.5667, 1.65204, 1.97772.
- `logistic_pdf_array`: SciFort 3.02184, 1.99181, 2.13609, 2.5506, 2.02512, 1.95694, 2.00665; SciPy 2.68077, 2.12223, 2.49449, 2.21249, 1.88485, 1.79855, 1.81323.
- `logistic_logpdf_array`: SciFort 1.86246, 1.50806, 1.36922, 1.4013, 1.22442, 1.11025, 1.19538; SciPy 1.87303, 1.9719, 1.41317, 1.40007, 1.20941, 1.26305, 1.2351.
- `logistic_cdf_array`: SciFort 0.827623, 0.87787, 0.850193, 0.832433, 0.713242, 0.550991, 0.509104; SciPy 3.52555, 3.09067, 3.32893, 3.16453, 2.58537, 2.49775, 2.51324.
- `logistic_sf_array`: SciFort 0.331474, 0.29917, 0.279407, 0.251801, 0.261194, 0.40383, 0.32712; SciPy 1.60745, 1.49308, 1.34497, 1.32468, 1.33344, 1.36375, 1.52087.
- `logistic_logcdf_array`: SciFort 0.749983, 0.767916, 0.882376, 0.973038, 0.947466, 1.1435, 1.31288; SciPy 1.62841, 1.78027, 1.96664, 2.11678, 2.33687, 2.50319, 3.03572.
- `logistic_logsf_array`: SciFort 1.80137, 1.51524, 1.63859, 1.67622, 1.48582, 2.02317, 1.6654; SciPy 3.94615, 3.64417, 4.14946, 4.22908, 3.46441, 3.79256, 4.47524.
- `logistic_ppf_array`: SciFort 1.51932, 1.25506, 1.5929, 1.59606, 1.26882, 1.16745, 1.13954; SciPy 2.45213, 2.15409, 2.40354, 2.9514, 1.91035, 1.83965, 1.94545.
- `logistic_isf_array`: SciFort 1.13675, 1.08236, 1.01977, 1.16359, 1.19276, 1.08959, 0.98512; SciPy 1.91884, 2.08915, 1.76387, 1.88849, 1.94133, 1.80453, 1.63176.
- `laplace_pdf_array`: SciFort 0.499288, 0.501217, 0.579361, 0.472314, 0.522304, 0.601508, 0.688482; SciPy 0.992074, 0.999883, 1.05639, 0.955216, 1.09796, 0.887773, 1.41813.
- `laplace_logpdf_array`: SciFort 0.616951, 0.584135, 0.617334, 0.570423, 0.648498, 0.514541, 0.492362; SciPy 2.36479, 1.97969, 1.81851, 2.07908, 2.47794, 2.58802, 1.80182.
- `laplace_cdf_array`: SciFort 0.357511, 0.390734, 0.337664, 0.364973, 0.369229, 0.399801, 0.382514; SciPy 1.10119, 1.18109, 1.16607, 1.22396, 1.01637, 1.11565, 0.912808.
- `laplace_sf_array`: SciFort 0.359152, 0.311957, 0.343766, 0.312279, 0.344246, 0.388177, 0.33139; SciPy 1.17405, 1.25422, 1.26282, 1.18566, 1.2449, 1.29932, 1.22396.
- `laplace_logcdf_array`: SciFort 0.911535, 0.905005, 0.82825, 0.729736, 0.728761, 0.772182, 0.705932; SciPy 2.71563, 2.27223, 1.94677, 2.1438, 2.15898, 2.49013, 2.69742.
- `laplace_logsf_array`: SciFort 1.0655, 1.07526, 0.945575, 0.950773, 0.934054, 0.85167, 0.831952; SciPy 3.19509, 3.526, 1.99764, 2.32388, 3.1779, 2.80839, 2.07603.
- `laplace_ppf_array`: SciFort 0.629188, 0.995623, 0.69958, 0.542207, 0.48649, 0.718532, 0.934295; SciPy 1.76989, 2.14544, 1.31621, 1.41177, 1.54636, 1.6003, 3.46155.
- `laplace_isf_array`: SciFort 0.496566, 0.49799, 0.431673, 0.547173, 0.523435, 0.391605, 0.367624; SciPy 1.46814, 1.40472, 1.34814, 1.41766, 1.3796, 1.26467, 1.17855.
- `cauchy_pdf_array`: SciFort 1.14338, 1.26658, 1.29581, 1.36722, 1.43596, 1.68757, 1.28372; SciPy 0.601097, 0.601759, 0.727739, 0.824536, 0.86273, 0.891412, 0.750668.
- `cauchy_logpdf_array`: SciFort 0.903299, 1.18574, 1.24851, 0.68829, 0.682326, 0.773597, 0.869065; SciPy 1.93251, 2.15145, 2.21998, 2.00714, 1.73597, 1.83625, 1.56617.
- `cauchy_cdf_array`: SciFort 0.507609, 0.432165, 0.371383, 0.375396, 0.414341, 0.421179, 0.388514; SciPy 1.16145, 0.943883, 0.960965, 1.21241, 1.25047, 1.08348, 1.1078.
- `cauchy_sf_array`: SciFort 0.571721, 0.613407, 0.781694, 0.782318, 0.699931, 0.726789, 0.745857; SciPy 1.20973, 1.1428, 1.70977, 1.99988, 1.61082, 1.39971, 1.8488.
- `cauchy_logcdf_array`: SciFort 1.55315, 2.24818, 2.01518, 1.58129, 1.71574, 1.77577, 1.73203; SciPy 2.45203, 2.98965, 2.68654, 2.59448, 2.64315, 3.09628, 2.5525.
- `cauchy_logsf_array`: SciFort 1.74097, 1.62898, 1.60951, 2.29731, 2.15415, 1.61472, 1.55616; SciPy 2.39016, 2.61495, 2.28155, 2.88584, 2.65663, 2.3153, 2.2376.
- `cauchy_ppf_array`: SciFort 0.884276, 0.769155, 0.85257, 0.816863, 0.639874, 0.777686, 0.802848; SciPy 1.74666, 1.90537, 2.26699, 2.29176, 1.79913, 1.64532, 2.02261.
- `cauchy_isf_array`: SciFort 0.829383, 0.849051, 0.920276, 1.06151, 1.02826, 0.927977, 1.02966; SciPy 1.75125, 2.06771, 1.89793, 2.40829, 2.24355, 2.14903, 2.69512.
- `weibull_pdf_array`: SciFort 3.07202, 2.7778, 2.71564, 3.37498, 3.67509, 3.28642, 3.47337; SciPy 3.70889, 3.56704, 4.18393, 3.34704, 4.00723, 4.61317, 3.80041.
- `weibull_logpdf_array`: SciFort 2.1753, 2.24906, 2.63109, 2.17663, 2.0206, 2.67317, 2.13701; SciPy 4.18107, 4.33732, 4.89344, 4.56707, 4.27105, 4.58092, 3.95236.
- `weibull_cdf_array`: SciFort 1.91573, 1.67015, 1.63978, 1.8404, 1.71753, 1.7794, 1.81719; SciPy 3.24461, 3.02023, 2.95091, 2.77494, 2.91897, 3.00074, 3.04404.
- `weibull_sf_array`: SciFort 1.46961, 1.57089, 1.50383, 1.59825, 1.69959, 1.72587, 1.54974; SciPy 1.6758, 1.93581, 1.74979, 1.84715, 2.6405, 2.32368, 2.42805.
- `weibull_logcdf_array`: SciFort 3.3319, 2.7306, 2.71517, 3.39025, 2.75468, 2.58901, 2.74468; SciPy 4.94895, 4.67872, 4.45276, 5.43134, 3.73135, 4.18939, 4.59694.
- `weibull_logsf_array`: SciFort 1.34005, 1.35294, 1.35059, 1.23674, 1.55266, 1.28522, 1.30002; SciPy 1.80375, 2.32832, 1.91712, 1.94742, 2.32636, 2.12376, 1.94072.
- `weibull_ppf_array`: SciFort 2.3221, 2.95884, 4.42736, 2.29037, 2.2505, 2.16507, 2.16904; SciPy 3.55321, 3.06733, 3.45068, 3.67119, 2.66222, 2.7037, 3.52216.
- `weibull_isf_array`: SciFort 2.4408, 1.97525, 1.91312, 2.35097, 2.3675, 1.77749, 1.77578; SciPy 3.37799, 3.39005, 2.3453, 2.82306, 2.72135, 2.47081, 1.8591.
- `lognormal_pdf_array`: SciFort 1.52401, 1.56392, 1.48812, 1.58062, 1.43983, 1.36427, 1.49308; SciPy 2.0689, 1.84072, 2.27027, 2.64165, 1.97894, 2.12797, 1.69957.
- `lognormal_logpdf_array`: SciFort 1.81247, 1.44989, 1.58607, 1.63196, 1.74494, 1.418, 1.84352; SciPy 2.00964, 1.79795, 2.17987, 2.23235, 1.91442, 1.83975, 2.02337.
- `lognormal_cdf_array`: SciFort 2.15402, 1.98321, 2.50785, 1.79953, 1.83191, 1.97854, 1.91682; SciPy 3.39874, 3.36828, 3.4835, 3.05395, 2.82559, 3.12774, 2.79417.
- `lognormal_sf_array`: SciFort 2.30186, 1.95788, 1.93872, 2.23214, 2.39366, 3.98625, 3.34819; SciPy 2.58247, 2.69805, 3.21042, 3.02601, 2.71532, 3.24507, 3.52416.
- `lognormal_logcdf_array`: SciFort 3.9399, 3.62845, 3.57185, 4.01934, 4.26858, 3.07531, 3.37974; SciPy 5.05487, 4.9262, 4.90404, 6.40097, 4.62175, 4.38639, 5.23925.
- `lognormal_logsf_array`: SciFort 4.04184, 3.28563, 3.22589, 3.8715, 3.54745, 4.03922, 3.47171; SciPy 4.21661, 3.99303, 4.88313, 4.84674, 3.95926, 4.03594, 3.86411.
- `lognormal_ppf_array`: SciFort 2.32697, 2.85382, 2.29729, 2.79605, 2.94538, 1.80099, 2.00054; SciPy 2.68624, 3.47651, 2.56992, 2.95095, 3.45941, 2.67535, 2.24259.
- `lognormal_isf_array`: SciFort 1.81852, 1.53693, 1.49405, 1.21496, 1.35756, 1.21576, 1.14638; SciPy 2.33027, 1.79085, 1.66641, 1.53096, 1.50508, 1.6074, 1.53643.
- `normal_cdf_scalar`: SciFort 0.0120452, 0.0135912, 0.0126103, 0.0121687, 0.0137448, 0.0136694, 0.0126084; SciPy 0.533573, 0.497787, 0.633793, 0.542589, 0.490043, 0.54215, 0.626281.
- `fft_power2_reuse`: SciFort 0.229228, 0.212848, 0.247831, 0.259757, 0.22971, 0.272076, 0.230535; SciPy 0.112856, 0.0887271, 0.182951, 0.21393, 0.127747, 0.16589, 0.110124.
- `fft_power2_setup`: SciFort 0.230974, 0.220007, 0.248039, 0.217626, 0.234832, 0.224306, 0.229016; SciPy 0.115232, 0.116185, 0.114595, 0.1187, 0.125353, 0.110932, 0.118216.
- `fft_mixed_reuse`: SciFort 0.643141, 0.550538, 0.493126, 0.625025, 0.790277, 0.652427, 0.668483; SciPy 0.159303, 0.164615, 0.167485, 0.167549, 0.199737, 0.206816, 0.191905.
- `fft_mixed_setup`: SciFort 1.29374, 1.12817, 1.1864, 1.26621, 1.27069, 1.1362, 1.09017; SciPy 0.249654, 0.188484, 0.253189, 0.220405, 0.251001, 0.198363, 0.235628.
- `fft_prime_reuse`: SciFort 1.75868, 1.19822, 1.22155, 1.69649, 1.29616, 1.24731, 1.22751; SciPy 0.652, 0.688341, 0.73622, 0.968007, 0.617863, 0.649048, 0.700065.
- `fft_prime_setup`: SciFort 2.95001, 2.70827, 2.54476, 2.5579, 2.48788, 3.06542, 2.58242; SciPy 0.665268, 0.540135, 0.60128, 0.617257, 0.546466, 0.527927, 0.594559.
- `ifft_power2_reuse`: SciFort 0.302517, 0.244267, 0.244704, 0.294421, 0.290249, 0.259999, 0.243585; SciPy 0.160225, 0.154213, 0.121366, 0.147887, 0.156418, 0.144294, 0.163987.
- `ifft_power2_setup`: SciFort 0.26435, 0.310105, 0.285205, 0.263759, 0.246645, 0.232604, 0.234799; SciPy 0.143746, 0.151553, 0.156109, 0.126446, 0.139922, 0.153231, 0.133846.
- `ifft_mixed_reuse`: SciFort 0.652756, 0.620775, 0.736192, 0.604574, 0.612838, 0.701763, 0.700839; SciPy 0.183117, 0.178504, 0.191149, 0.194047, 0.182175, 0.182516, 0.163647.
- `ifft_mixed_setup`: SciFort 1.0247, 1.0439, 1.00535, 0.963301, 1.18779, 0.945298, 0.90878; SciPy 0.166171, 0.24782, 0.149077, 0.170707, 0.220679, 0.182777, 0.150583.
- `ifft_prime_reuse`: SciFort 0.947017, 0.943215, 1.01015, 0.922819, 0.896183, 0.870121, 0.918277; SciPy 0.440509, 0.471184, 0.494291, 0.393578, 0.393469, 0.463888, 0.483223.
- `ifft_prime_setup`: SciFort 2.01114, 2.26984, 2.30668, 2.15586, 2.12144, 2.31396, 2.32502; SciPy 0.356894, 0.539567, 0.488764, 0.399336, 0.478229, 0.546869, 0.448906.
- `lu_factor_solve`: SciFort 1.0822, 0.88406, 0.89705, 0.851641, 0.95704, 1.07681, 0.946837; SciPy 0.48312, 0.433232, 0.366782, 0.342997, 0.420262, 0.462736, 0.427328.
- `lu_reuse`: SciFort 0.0843744, 0.0680307, 0.0724327, 0.0716239, 0.0666489, 0.0871903, 0.0681377; SciPy 0.0384058, 0.0374492, 0.0295638, 0.0325613, 0.0321605, 0.0316544, 0.0315502.
- `lu_transpose`: SciFort 0.061519, 0.0611982, 0.0617258, 0.0667952, 0.0623785, 0.0591888, 0.0612671; SciPy 0.0339625, 0.0287683, 0.0308555, 0.0343033, 0.0290521, 0.0307176, 0.0271163.
- `cho_factor_solve`: SciFort 0.972517, 0.858103, 0.815738, 0.97963, 0.858967, 0.843759, 0.879321; SciPy 0.190891, 0.194076, 0.182196, 0.201563, 0.202308, 0.178226, 0.19228.
- `cho_reuse`: SciFort 0.11969, 0.114183, 0.11192, 0.0948027, 0.0920277, 0.159196, 0.228914; SciPy 0.0888581, 0.0805665, 0.0685009, 0.0631727, 0.0685274, 0.0735563, 0.0830697.
- `triangular_lower`: SciFort 0.0485421, 0.0450233, 0.0525534, 0.0474263, 0.0508116, 0.040138, 0.042615; SciPy 0.0498252, 0.050998, 0.0739157, 0.0766885, 0.0469755, 0.0457601, 0.0515601.
- `triangular_upper`: SciFort 0.0626723, 0.0460837, 0.0439619, 0.0465954, 0.047331, 0.0456984, 0.0504348; SciPy 0.067914, 0.0565549, 0.0465218, 0.0467293, 0.0496637, 0.0461047, 0.0543228.
- `triangular_transpose`: SciFort 0.0467513, 0.0451953, 0.047257, 0.0452745, 0.0468184, 0.0438885, 0.0453473; SciPy 0.0818153, 0.0636821, 0.0688381, 0.0617011, 0.0668551, 0.0674975, 0.0596402.
- `barycentric`: SciFort 3.6526, 2.67347, 2.58078, 2.72147, 2.75657, 2.8236, 2.59255; SciPy 2.97299, 2.89227, 2.90065, 2.73216, 2.92001, 3.06792, 2.78195.
- `barycentric_setup`: SciFort 3.24728, 2.71978, 2.76539, 4.16686, 3.04472, 3.66387, 3.04907; SciPy 3.63838, 4.18277, 4.14886, 6.80431, 3.57419, 4.3673, 3.89412.
- `barycentric_derivative`: SciFort 2.81779, 2.71213, 2.74717, 3.3393, 2.71315, 2.94215, 2.8221; SciPy 2.69483, 3.23581, 3.54508, 3.83499, 2.79896, 2.74107, 2.79306.

## Validation and reproduction

Only benchmark tooling and documentation changed. No FPM build or numerical
suite rerun was needed. GNU Fortran 17.0.0 20260510 (experimental) was used for
both the optimized run (`-O3`) and the checked builds (`-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`); both also used
`-std=f2018 -fimplicit-none -Werror=implicit-interface`.

```text
python benchmarks/compare_scipy.py --check-only --output build/scipy-speed/expanded-checked.json
python benchmarks/compare_scipy.py --check-only --quick --output build/scipy-speed/expanded-quick-checked.json
python benchmarks/compare_scipy.py --output build/scipy-speed/expanded-windows-20261122.json
python benchmarks/test_compare_scipy.py
```

All 80 cases passed at default sizes and quick sizes with runtime checks.
The optimized run validated the same 80 default-size cases before timing, then
collected seven samples per backend with a calibrated batch target of 0.1 s.
Default validation compared 2,058,360 output values; the maximum scaled solve
residual was 8.924e-16. Seven harness regression checks passed, including
intentional corruption of CDF tails, imaginary FFT outputs, known solutions,
and timed checksums, plus CLI selection and Markdown reporting checks.

The local raw reports and their generated Markdown counterparts are under
`build/scipy-speed/` (ignored by Git). This tracked snapshot preserves all
medians and individual samples from the optimized run.

Raw timing JSON SHA-256: `d98a373b2e60fcac5abb6716964dafd9778266b632ee17c4a9d3add8a5f4f814`.

Benchmark source/build SHA-256: `f0ea820f2b1d3f9f69a9a1bff47e5e9163f028f5d240f62412182cd52152d968`.

Harness SHA-256: `ae7fc0adc0d1561de0ac0e15aea2bf91e97a61b181eff8ab2fb726a0a293fa2e`.
