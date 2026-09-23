! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Deterministic bounded maximum-likelihood fitting. Equal lower and upper bounds
! fix a parameter. The optimizer is a projected coordinate/pattern search with
! optional integer coordinates and no hidden global state.

module scifort_fit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_likelihood, only : normal_nnlf, &
        uniform_nnlf, &
        exponential_nnlf, &
        laplace_nnlf, &
        logistic_nnlf, &
        cauchy_nnlf, &
        rayleigh_nnlf, &
        gamma_nnlf, &
        chi2_nnlf, &
        t_nnlf, &
        lognormal_nnlf, &
        weibull_nnlf, &
        pareto_nnlf, &
        beta_nnlf, &
        f_nnlf, &
        bernoulli_nnlf, &
        poisson_nnlf, &
        geometric_nnlf, &
        binomial_nnlf, &
        negative_binomial_nnlf, &
        gumbel_r_nnlf, gumbel_l_nnlf, powerlaw_nnlf, triang_nnlf, genpareto_nnlf, &
        arcsine_nnlf, halfnorm_nnlf, halfcauchy_nnlf, lomax_nnlf, &
        chi_nnlf, maxwell_nnlf, cosine_nnlf, semicircular_nnlf, &
        anglit_nnlf, moyal_nnlf, hypsecant_nnlf, halflogistic_nnlf, &
        invgamma_nnlf, invgauss_nnlf, levy_nnlf, loglaplace_nnlf, &
        bradford_nnlf, truncexpon_nnlf, fisk_nnlf, dweibull_nnlf, &
        alpha_nnlf, fatiguelife_nnlf, genlogistic_nnlf, gennorm_nnlf, &
        nakagami_nnlf, powernorm_nnlf, loggamma_nnlf, wald_nnlf, &
        gompertz_nnlf, invweibull_nnlf, betaprime_nnlf, burr12_nnlf, &
        genhalflogistic_nnlf, exponpow_nnlf, exponweib_nnlf, powerlognorm_nnlf, &
        levy_l_nnlf, weibull_max_nnlf, rdist_nnlf, skewcauchy_nnlf, &
        dgamma_nnlf, laplace_asymmetric_nnlf, truncnorm_nnlf, loguniform_nnlf, &
        foldnorm_nnlf, foldcauchy_nnlf, recipinvgauss_nnlf, truncpareto_nnlf, &
        exponnorm_nnlf, johnsonsb_nnlf, johnsonsu_nnlf, trapezoid_nnlf, &
        burr_nnlf, mielke_nnlf, gibrat_nnlf, wrapcauchy_nnlf, &
        genextreme_nnlf, kappa3_nnlf, kappa4_nnlf, truncweibull_min_nnlf, &
        gengamma_nnlf, halfgennorm_nnlf, argus_nnlf, erlang_nnlf, &
        crystalball_nnlf, jf_skew_t_nnlf, pearson3_nnlf, rel_breitwigner_nnlf, &
        genexpon_nnlf, skewnorm_nnlf, tukeylambda_nnlf, rice_nnlf, &
        dpareto_lognorm_nnlf, vonmises_nnlf, vonmises_line_nnlf, kstwobign_nnlf, &
        irwinhall_nnlf, ksone_nnlf, ncx2_nnlf, ncf_nnlf, &
        randint_nnlf, planck_nnlf, dlaplace_nnlf, logser_nnlf, &
        betabinom_nnlf, hypergeom_nnlf, nhypergeom_nnlf, boltzmann_nnlf, &
        betanbinom_nnlf, yulesimon_nnlf, zipf_nnlf, zipfian_nnlf, &
        geninvgauss_nnlf, norminvgauss_nnlf, skellam_nnlf, genhyperbolic_nnlf, &
        nchypergeom_fisher_nnlf, nct_nnlf, gausshyper_nnlf, landau_nnlf, &
        nchypergeom_wallenius_nnlf, poisson_binom_nnlf, kstwo_nnlf, &
        levy_stable_nnlf, studentized_range_nnlf
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: fit_status_success = 0
    integer, parameter, public :: fit_status_max_iter = 1
    integer, parameter, public :: fit_status_invalid_input = 2
    integer, parameter, public :: fit_status_no_finite_objective = 3

    integer, parameter :: family_normal = 1
    integer, parameter :: family_uniform = 2
    integer, parameter :: family_exponential = 3
    integer, parameter :: family_laplace = 4
    integer, parameter :: family_logistic = 5
    integer, parameter :: family_cauchy = 6
    integer, parameter :: family_rayleigh = 7
    integer, parameter :: family_gamma = 8
    integer, parameter :: family_chi2 = 9
    integer, parameter :: family_t = 10
    integer, parameter :: family_lognormal = 11
    integer, parameter :: family_weibull = 12
    integer, parameter :: family_pareto = 13
    integer, parameter :: family_beta = 14
    integer, parameter :: family_f = 15
    integer, parameter :: family_bernoulli = 16
    integer, parameter :: family_poisson = 17
    integer, parameter :: family_geometric = 18
    integer, parameter :: family_binomial = 19
    integer, parameter :: family_negative_binomial = 20
    integer, parameter :: family_gumbel_r = 21
    integer, parameter :: family_gumbel_l = 22
    integer, parameter :: family_powerlaw = 23
    integer, parameter :: family_triang = 24
    integer, parameter :: family_genpareto = 25
    integer, parameter :: family_arcsine = 26
    integer, parameter :: family_halfnorm = 27
    integer, parameter :: family_halfcauchy = 28
    integer, parameter :: family_lomax = 29
    integer, parameter :: family_chi = 30
    integer, parameter :: family_maxwell = 31
    integer, parameter :: family_cosine = 32
    integer, parameter :: family_semicircular = 33
    integer, parameter :: family_anglit = 34
    integer, parameter :: family_moyal = 35
    integer, parameter :: family_hypsecant = 36
    integer, parameter :: family_halflogistic = 37
    integer, parameter :: family_invgamma = 38
    integer, parameter :: family_invgauss = 39
    integer, parameter :: family_levy = 40
    integer, parameter :: family_loglaplace = 41
    integer, parameter :: family_bradford = 42
    integer, parameter :: family_truncexpon = 43
    integer, parameter :: family_fisk = 44
    integer, parameter :: family_dweibull = 45
    integer, parameter :: family_alpha = 46
    integer, parameter :: family_fatiguelife = 47
    integer, parameter :: family_genlogistic = 48
    integer, parameter :: family_gennorm = 49
    integer, parameter :: family_nakagami = 50
    integer, parameter :: family_powernorm = 51
    integer, parameter :: family_loggamma = 52
    integer, parameter :: family_wald = 53
    integer, parameter :: family_gompertz = 54
    integer, parameter :: family_invweibull = 55
    integer, parameter :: family_betaprime = 56
    integer, parameter :: family_burr12 = 57
    integer, parameter :: family_genhalflogistic = 58
    integer, parameter :: family_exponpow = 59
    integer, parameter :: family_exponweib = 60
    integer, parameter :: family_powerlognorm = 61
    integer, parameter :: family_levy_l = 62
    integer, parameter :: family_weibull_max = 63
    integer, parameter :: family_rdist = 64
    integer, parameter :: family_skewcauchy = 65
    integer, parameter :: family_dgamma = 66
    integer, parameter :: family_laplace_asymmetric = 67
    integer, parameter :: family_truncnorm = 68
    integer, parameter :: family_loguniform = 69
    integer, parameter :: family_foldnorm = 70
    integer, parameter :: family_foldcauchy = 71
    integer, parameter :: family_recipinvgauss = 72
    integer, parameter :: family_truncpareto = 73
    integer, parameter :: family_exponnorm = 74
    integer, parameter :: family_johnsonsb = 75
    integer, parameter :: family_johnsonsu = 76
    integer, parameter :: family_trapezoid = 77
    integer, parameter :: family_burr = 78
    integer, parameter :: family_mielke = 79
    integer, parameter :: family_gibrat = 80
    integer, parameter :: family_wrapcauchy = 81
    integer, parameter :: family_genextreme = 82
    integer, parameter :: family_kappa3 = 83
    integer, parameter :: family_kappa4 = 84
    integer, parameter :: family_truncweibull_min = 85
    integer, parameter :: family_gengamma = 86
    integer, parameter :: family_halfgennorm = 87
    integer, parameter :: family_argus = 88
    integer, parameter :: family_erlang = 89
    integer, parameter :: family_crystalball = 90
    integer, parameter :: family_jf_skew_t = 91
    integer, parameter :: family_pearson3 = 92
    integer, parameter :: family_rel_breitwigner = 93
    integer, parameter :: family_genexpon = 94
    integer, parameter :: family_skewnorm = 95
    integer, parameter :: family_tukeylambda = 96
    integer, parameter :: family_rice = 97
    integer, parameter :: family_dpareto_lognorm = 98
    integer, parameter :: family_vonmises = 99
    integer, parameter :: family_vonmises_line = 100
    integer, parameter :: family_kstwobign = 101
    integer, parameter :: family_irwinhall = 102
    integer, parameter :: family_ksone = 103
    integer, parameter :: family_ncx2 = 104
    integer, parameter :: family_ncf = 105
    integer, parameter :: family_randint = 106
    integer, parameter :: family_planck = 107
    integer, parameter :: family_dlaplace = 108
    integer, parameter :: family_logser = 109
    integer, parameter :: family_betabinom = 110
    integer, parameter :: family_hypergeom = 111
    integer, parameter :: family_nhypergeom = 112
    integer, parameter :: family_boltzmann = 113
    integer, parameter :: family_betanbinom = 114
    integer, parameter :: family_yulesimon = 115
    integer, parameter :: family_zipf = 116
    integer, parameter :: family_zipfian = 117
    integer, parameter :: family_geninvgauss = 118
    integer, parameter :: family_norminvgauss = 119
    integer, parameter :: family_skellam = 120
    integer, parameter :: family_genhyperbolic = 121
    integer, parameter :: family_nchypergeom_fisher = 122
    integer, parameter :: family_nct = 123
    integer, parameter :: family_gausshyper = 124
    integer, parameter :: family_landau = 125
    integer, parameter :: family_nchypergeom_wallenius = 126
    integer, parameter :: family_poisson_binom = 127
    integer, parameter :: family_kstwo = 128
    integer, parameter :: family_levy_stable = 129
    integer, parameter :: family_studentized_range = 130

    type, public :: fit_result
        real(dp), allocatable :: params(:)
        real(dp) :: nllf = 0.0_dp
        integer :: iterations = 0
        integer :: evaluations = 0
        integer :: status = fit_status_invalid_input
        logical :: success = .false.
    end type fit_result

    public :: normal_fit
    public :: uniform_fit
    public :: exponential_fit
    public :: laplace_fit
    public :: logistic_fit
    public :: cauchy_fit
    public :: rayleigh_fit
    public :: gamma_fit
    public :: chi2_fit
    public :: t_fit
    public :: lognormal_fit
    public :: weibull_fit
    public :: pareto_fit
    public :: beta_fit
    public :: f_fit
    public :: bernoulli_fit
    public :: poisson_fit
    public :: geometric_fit
    public :: binomial_fit
    public :: negative_binomial_fit
    public :: gumbel_r_fit
    public :: gumbel_l_fit
    public :: powerlaw_fit
    public :: triang_fit
    public :: genpareto_fit
    public :: arcsine_fit
    public :: halfnorm_fit
    public :: halfcauchy_fit
    public :: lomax_fit
    public :: chi_fit
    public :: maxwell_fit
    public :: cosine_fit
    public :: semicircular_fit
    public :: anglit_fit
    public :: moyal_fit
    public :: hypsecant_fit
    public :: halflogistic_fit
    public :: invgamma_fit
    public :: invgauss_fit
    public :: levy_fit
    public :: loglaplace_fit
    public :: bradford_fit
    public :: truncexpon_fit
    public :: fisk_fit
    public :: dweibull_fit
    public :: alpha_fit
    public :: fatiguelife_fit
    public :: genlogistic_fit
    public :: gennorm_fit
    public :: nakagami_fit
    public :: powernorm_fit
    public :: loggamma_fit
    public :: wald_fit
    public :: gompertz_fit
    public :: invweibull_fit
    public :: betaprime_fit
    public :: burr12_fit
    public :: genhalflogistic_fit
    public :: exponpow_fit
    public :: exponweib_fit
    public :: powerlognorm_fit
    public :: levy_l_fit
    public :: weibull_max_fit
    public :: rdist_fit
    public :: skewcauchy_fit
    public :: dgamma_fit
    public :: laplace_asymmetric_fit
    public :: truncnorm_fit
    public :: loguniform_fit
    public :: foldnorm_fit
    public :: foldcauchy_fit
    public :: recipinvgauss_fit
    public :: truncpareto_fit
    public :: exponnorm_fit
    public :: johnsonsb_fit
    public :: johnsonsu_fit
    public :: trapezoid_fit
    public :: burr_fit
    public :: mielke_fit
    public :: gibrat_fit
    public :: wrapcauchy_fit
    public :: genextreme_fit
    public :: kappa3_fit
    public :: kappa4_fit
    public :: truncweibull_min_fit
    public :: gengamma_fit
    public :: halfgennorm_fit
    public :: argus_fit
    public :: erlang_fit
    public :: crystalball_fit
    public :: jf_skew_t_fit
    public :: pearson3_fit
    public :: rel_breitwigner_fit
    public :: genexpon_fit
    public :: skewnorm_fit
    public :: tukeylambda_fit
    public :: rice_fit
    public :: dpareto_lognorm_fit
    public :: vonmises_fit
    public :: vonmises_line_fit
    public :: kstwobign_fit
    public :: irwinhall_fit
    public :: ksone_fit
    public :: kstwo_fit
    public :: levy_stable_fit
    public :: studentized_range_fit
    public :: ncx2_fit
    public :: ncf_fit
    public :: randint_fit
    public :: planck_fit
    public :: dlaplace_fit
    public :: logser_fit
    public :: betabinom_fit
    public :: hypergeom_fit
    public :: nhypergeom_fit
    public :: boltzmann_fit
    public :: betanbinom_fit
    public :: yulesimon_fit
    public :: zipf_fit
    public :: zipfian_fit
    public :: geninvgauss_fit
    public :: norminvgauss_fit
    public :: skellam_fit
    public :: genhyperbolic_fit
    public :: nchypergeom_fisher_fit
    public :: nct_fit
    public :: gausshyper_fit
    public :: landau_fit
    public :: nchypergeom_wallenius_fit
    public :: poisson_binom_fit

    interface bernoulli_fit
        module procedure bernoulli_fit_real
        module procedure bernoulli_fit_int
    end interface bernoulli_fit
    interface poisson_fit
        module procedure poisson_fit_real
        module procedure poisson_fit_int
    end interface poisson_fit
    interface geometric_fit
        module procedure geometric_fit_real
        module procedure geometric_fit_int
    end interface geometric_fit
    interface binomial_fit
        module procedure binomial_fit_real
        module procedure binomial_fit_int
    end interface binomial_fit
    interface negative_binomial_fit
        module procedure negative_binomial_fit_real
        module procedure negative_binomial_fit_int
    end interface negative_binomial_fit
    interface randint_fit
        module procedure randint_fit_real
        module procedure randint_fit_int
    end interface randint_fit
    interface planck_fit
        module procedure planck_fit_real
        module procedure planck_fit_int
    end interface planck_fit
    interface dlaplace_fit
        module procedure dlaplace_fit_real
        module procedure dlaplace_fit_int
    end interface dlaplace_fit
    interface logser_fit
        module procedure logser_fit_real
        module procedure logser_fit_int
    end interface logser_fit
    interface betabinom_fit
        module procedure betabinom_fit_real
        module procedure betabinom_fit_int
    end interface betabinom_fit
    interface hypergeom_fit
        module procedure hypergeom_fit_real
        module procedure hypergeom_fit_int
    end interface hypergeom_fit
    interface nhypergeom_fit
        module procedure nhypergeom_fit_real
        module procedure nhypergeom_fit_int
    end interface nhypergeom_fit
    interface boltzmann_fit
        module procedure boltzmann_fit_real
        module procedure boltzmann_fit_int
    end interface boltzmann_fit

    interface betanbinom_fit
        module procedure betanbinom_fit_real
        module procedure betanbinom_fit_int
    end interface betanbinom_fit
    interface yulesimon_fit
        module procedure yulesimon_fit_real
        module procedure yulesimon_fit_int
    end interface yulesimon_fit
    interface zipf_fit
        module procedure zipf_fit_real
        module procedure zipf_fit_int
    end interface zipf_fit
    interface zipfian_fit
        module procedure zipfian_fit_real
        module procedure zipfian_fit_int
    end interface zipfian_fit

    interface skellam_fit
        module procedure skellam_fit_real
        module procedure skellam_fit_int
    end interface skellam_fit

    interface nchypergeom_fisher_fit
        module procedure nchypergeom_fisher_fit_real
        module procedure nchypergeom_fisher_fit_int
    end interface nchypergeom_fisher_fit

    interface nchypergeom_wallenius_fit
        module procedure nchypergeom_wallenius_fit_real
        module procedure nchypergeom_wallenius_fit_int
    end interface nchypergeom_wallenius_fit

    interface poisson_binom_fit
        module procedure poisson_binom_fit_real
        module procedure poisson_binom_fit_int
    end interface poisson_binom_fit

contains

    subroutine bounded_fit(data, family, lower, upper, result, guess, integral, &
        max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! observations passed to the family objective
        integer, intent(in) :: family !! internal distribution identifier
        real(dp), intent(in) :: lower(:) !! finite lower bounds for every parameter
        real(dp), intent(in) :: upper(:) !! finite upper bounds for every parameter
        type(fit_result), intent(out) :: result !! minimizer, objective, and diagnostics
        real(dp), intent(in), optional :: guess(:) !! starting vector; midpoint when absent
        logical, intent(in) :: integral(:) !! true for integer-constrained parameters
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        real(dp), parameter :: penalty_fraction = 16.0_dp
        logical :: improved
        integer :: i
        integer :: iter
        integer :: iteration_limit
        integer :: n
        integer :: sign_index
        real(dp) :: best
        real(dp) :: candidate
        real(dp), allocatable :: delta(:)
        real(dp), allocatable :: hi(:)
        real(dp), allocatable :: lo(:)
        real(dp) :: penalty
        real(dp), allocatable :: start(:)
        real(dp), allocatable :: step(:)
        real(dp) :: tol
        real(dp), allocatable :: trial(:)
        real(dp), allocatable :: x(:)

        n = size(lower)
        allocate(result%params(n))
        result%params = nan_value()
        result%nllf = nan_value()
        result%iterations = 0
        result%evaluations = 0
        result%status = fit_status_invalid_input
        result%success = .false.

        if (n == 0 .or. size(upper) /= n .or. size(integral) /= n) return
        if (.not. all(ieee_is_finite(lower)) .or. .not. all(ieee_is_finite(upper))) return
        if (any(lower > upper)) return
        if (present(guess)) then
            if (size(guess) /= n .or. .not. all(ieee_is_finite(guess))) return
        end if

        iteration_limit = 1000
        if (present(max_iter)) iteration_limit = max_iter
        if (iteration_limit <= 0) return
        tol = sqrt(epsilon(1.0_dp))
        if (present(tolerance)) tol = tolerance
        if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp) return

        allocate(lo(n), hi(n), x(n), trial(n), start(n), delta(n), step(n))
        lo = lower
        hi = upper
        do i = 1, n
            if (integral(i)) then
                candidate = lo(i)
                lo(i) = aint(candidate)
                if (lo(i) < candidate) lo(i) = lo(i) + 1.0_dp
                candidate = hi(i)
                hi(i) = aint(candidate)
                if (hi(i) > candidate) hi(i) = hi(i) - 1.0_dp
            end if
        end do
        if (any(lo > hi)) return

        if (present(guess)) then
            x = min(max(guess, lo), hi)
        else
            x = 0.5_dp * (lo + hi)
        end if
        call project_point(x, lo, hi, integral)

        step = 0.0_dp
        do i = 1, n
            if (hi(i) == lo(i)) cycle
            if (integral(i)) then
                step(i) = max(1.0_dp, aint(0.25_dp * (hi(i) - lo(i))))
            else
                step(i) = 0.25_dp * (hi(i) - lo(i))
            end if
        end do

        penalty = huge(1.0_dp) / penalty_fraction
        best = evaluate_family(family, data, x)
        result%evaluations = result%evaluations + 1
        if (.not. ieee_is_finite(best)) best = penalty

        do iter = 1, iteration_limit
            result%iterations = iter
            start = x
            improved = .false.
            do i = 1, n
                if (step(i) == 0.0_dp) cycle
                do sign_index = -1, 1, 2
                    trial = x
                    trial(i) = trial(i) + real(sign_index, dp) * step(i)
                    call project_point(trial, lo, hi, integral)
                    if (all(trial == x)) cycle
                    candidate = evaluate_family(family, data, trial)
                    result%evaluations = result%evaluations + 1
                    if (.not. ieee_is_finite(candidate)) candidate = penalty
                    if (candidate < best) then
                        x = trial
                        best = candidate
                        improved = .true.
                    end if
                end do
            end do

            if (improved) then
                delta = x - start
                trial = x + delta
                call project_point(trial, lo, hi, integral)
                if (.not. all(trial == x)) then
                    candidate = evaluate_family(family, data, trial)
                    result%evaluations = result%evaluations + 1
                    if (.not. ieee_is_finite(candidate)) candidate = penalty
                    if (candidate < best) then
                        x = trial
                        best = candidate
                    end if
                end if
            else
                do i = 1, n
                    if (step(i) == 0.0_dp) cycle
                    if (integral(i)) then
                        if (step(i) <= 1.0_dp) then
                            step(i) = 0.0_dp
                        else
                            step(i) = max(1.0_dp, aint(0.5_dp * step(i)))
                        end if
                    else
                        step(i) = 0.5_dp * step(i)
                        if (step(i) <= tol * (1.0_dp + abs(x(i)))) step(i) = 0.0_dp
                    end if
                end do
            end if

            if (all(step == 0.0_dp)) exit
        end do

        result%params = x
        result%nllf = evaluate_family(family, data, x)
        result%evaluations = result%evaluations + 1
        if (.not. ieee_is_finite(result%nllf) .or. best >= penalty) then
            result%status = fit_status_no_finite_objective
            result%success = .false.
        else if (all(step == 0.0_dp)) then
            result%status = fit_status_success
            result%success = .true.
        else
            result%status = fit_status_max_iter
            result%success = .false.
        end if
    end subroutine bounded_fit

    pure function evaluate_family(family, data, params) result(value)
        integer, intent(in) :: family !! internal distribution identifier
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: params(:) !! candidate parameters in canonical order
        real(dp) :: value

        select case (family)
        case (family_normal)
            value = normal_nnlf(data, params(1), params(2))
        case (family_uniform)
            value = uniform_nnlf(data, params(1), params(2))
        case (family_exponential)
            value = exponential_nnlf(data, params(1), params(2))
        case (family_laplace)
            value = laplace_nnlf(data, params(1), params(2))
        case (family_logistic)
            value = logistic_nnlf(data, params(1), params(2))
        case (family_cauchy)
            value = cauchy_nnlf(data, params(1), params(2))
        case (family_rayleigh)
            value = rayleigh_nnlf(data, params(1), params(2))
        case (family_gamma)
            value = gamma_nnlf(data, params(1), params(2), params(3))
        case (family_chi2)
            value = chi2_nnlf(data, params(1), params(2), params(3))
        case (family_t)
            value = t_nnlf(data, params(1), params(2), params(3))
        case (family_lognormal)
            value = lognormal_nnlf(data, params(1), params(2), params(3))
        case (family_weibull)
            value = weibull_nnlf(data, params(1), params(2), params(3))
        case (family_pareto)
            value = pareto_nnlf(data, params(1), params(2), params(3))
        case (family_beta)
            value = beta_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_f)
            value = f_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_bernoulli)
            value = bernoulli_nnlf(data, params(1), params(2))
        case (family_poisson)
            value = poisson_nnlf(data, params(1), params(2))
        case (family_geometric)
            value = geometric_nnlf(data, params(1), params(2))
        case (family_binomial)
            value = binomial_nnlf(data, params(1), params(2), params(3))
        case (family_negative_binomial)
            value = negative_binomial_nnlf(data, params(1), params(2), params(3))
        case (family_gumbel_r)
            value = gumbel_r_nnlf(data, params(1), params(2))
        case (family_gumbel_l)
            value = gumbel_l_nnlf(data, params(1), params(2))
        case (family_powerlaw)
            value = powerlaw_nnlf(data, params(1), params(2), params(3))
        case (family_triang)
            value = triang_nnlf(data, params(1), params(2), params(3))
        case (family_genpareto)
            value = genpareto_nnlf(data, params(1), params(2), params(3))
        case (family_arcsine)
            value = arcsine_nnlf(data, params(1), params(2))
        case (family_halfnorm)
            value = halfnorm_nnlf(data, params(1), params(2))
        case (family_halfcauchy)
            value = halfcauchy_nnlf(data, params(1), params(2))
        case (family_lomax)
            value = lomax_nnlf(data, params(1), params(2), params(3))
        case (family_chi)
            value = chi_nnlf(data, params(1), params(2), params(3))
        case (family_maxwell)
            value = maxwell_nnlf(data, params(1), params(2))
        case (family_cosine)
            value = cosine_nnlf(data, params(1), params(2))
        case (family_semicircular)
            value = semicircular_nnlf(data, params(1), params(2))
        case (family_anglit)
            value = anglit_nnlf(data, params(1), params(2))
        case (family_moyal)
            value = moyal_nnlf(data, params(1), params(2))
        case (family_hypsecant)
            value = hypsecant_nnlf(data, params(1), params(2))
        case (family_halflogistic)
            value = halflogistic_nnlf(data, params(1), params(2))
        case (family_invgamma)
            value = invgamma_nnlf(data, params(1), params(2), params(3))
        case (family_invgauss)
            value = invgauss_nnlf(data, params(1), params(2), params(3))
        case (family_levy)
            value = levy_nnlf(data, params(1), params(2))
        case (family_loglaplace)
            value = loglaplace_nnlf(data, params(1), params(2), params(3))
        case (family_bradford)
            value = bradford_nnlf(data, params(1), params(2), params(3))
        case (family_truncexpon)
            value = truncexpon_nnlf(data, params(1), params(2), params(3))
        case (family_fisk)
            value = fisk_nnlf(data, params(1), params(2), params(3))
        case (family_dweibull)
            value = dweibull_nnlf(data, params(1), params(2), params(3))
        case (family_alpha)
            value = alpha_nnlf(data, params(1), params(2), params(3))
        case (family_fatiguelife)
            value = fatiguelife_nnlf(data, params(1), params(2), params(3))
        case (family_genlogistic)
            value = genlogistic_nnlf(data, params(1), params(2), params(3))
        case (family_gennorm)
            value = gennorm_nnlf(data, params(1), params(2), params(3))
        case (family_nakagami)
            value = nakagami_nnlf(data, params(1), params(2), params(3))
        case (family_powernorm)
            value = powernorm_nnlf(data, params(1), params(2), params(3))
        case (family_loggamma)
            value = loggamma_nnlf(data, params(1), params(2), params(3))
        case (family_wald)
            value = wald_nnlf(data, params(1), params(2))
        case (family_gompertz)
            value = gompertz_nnlf(data, params(1), params(2), params(3))
        case (family_invweibull)
            value = invweibull_nnlf(data, params(1), params(2), params(3))
        case (family_betaprime)
            value = betaprime_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_burr12)
            value = burr12_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_genhalflogistic)
            value = genhalflogistic_nnlf(data, params(1), params(2), params(3))
        case (family_exponpow)
            value = exponpow_nnlf(data, params(1), params(2), params(3))
        case (family_exponweib)
            value = exponweib_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_powerlognorm)
            value = powerlognorm_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_levy_l)
            value = levy_l_nnlf(data, params(1), params(2))
        case (family_weibull_max)
            value = weibull_max_nnlf(data, params(1), params(2), params(3))
        case (family_rdist)
            value = rdist_nnlf(data, params(1), params(2), params(3))
        case (family_skewcauchy)
            value = skewcauchy_nnlf(data, params(1), params(2), params(3))
        case (family_dgamma)
            value = dgamma_nnlf(data, params(1), params(2), params(3))
        case (family_laplace_asymmetric)
            value = laplace_asymmetric_nnlf(data, params(1), params(2), params(3))
        case (family_truncnorm)
            value = truncnorm_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_loguniform)
            value = loguniform_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_foldnorm)
            value = foldnorm_nnlf(data, params(1), params(2), params(3))
        case (family_foldcauchy)
            value = foldcauchy_nnlf(data, params(1), params(2), params(3))
        case (family_recipinvgauss)
            value = recipinvgauss_nnlf(data, params(1), params(2), params(3))
        case (family_truncpareto)
            value = truncpareto_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_exponnorm)
            value = exponnorm_nnlf(data, params(1), params(2), params(3))
        case (family_johnsonsb)
            value = johnsonsb_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_johnsonsu)
            value = johnsonsu_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_trapezoid)
            value = trapezoid_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_burr)
            value = burr_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_mielke)
            value = mielke_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_gibrat)
            value = gibrat_nnlf(data, params(1), params(2))
        case (family_wrapcauchy)
            value = wrapcauchy_nnlf(data, params(1), params(2), params(3))
        case (family_genextreme)
            value = genextreme_nnlf(data, params(1), params(2), params(3))
        case (family_kappa3)
            value = kappa3_nnlf(data, params(1), params(2), params(3))
        case (family_kappa4)
            value = kappa4_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_truncweibull_min)
            value = truncweibull_min_nnlf(data, params(1), params(2), params(3), &
                params(4), params(5))
        case (family_gengamma)
            value = gengamma_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_halfgennorm)
            value = halfgennorm_nnlf(data, params(1), params(2), params(3))
        case (family_argus)
            value = argus_nnlf(data, params(1), params(2), params(3))
        case (family_erlang)
            value = erlang_nnlf(data, params(1), params(2), params(3))
        case (family_crystalball)
            value = crystalball_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_jf_skew_t)
            value = jf_skew_t_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_pearson3)
            value = pearson3_nnlf(data, params(1), params(2), params(3))
        case (family_rel_breitwigner)
            value = rel_breitwigner_nnlf(data, params(1), params(2), params(3))
        case (family_genexpon)
            value = genexpon_nnlf(data, params(1), params(2), params(3), params(4), params(5))
        case (family_skewnorm)
            value = skewnorm_nnlf(data, params(1), params(2), params(3))
        case (family_tukeylambda)
            value = tukeylambda_nnlf(data, params(1), params(2), params(3))
        case (family_rice)
            value = rice_nnlf(data, params(1), params(2), params(3))
        case (family_dpareto_lognorm)
            value = dpareto_lognorm_nnlf(data, params(1), params(2), params(3), params(4), &
                params(5), params(6))
        case (family_vonmises)
            value = vonmises_nnlf(data, params(1), params(2), params(3))
        case (family_vonmises_line)
            value = vonmises_line_nnlf(data, params(1), params(2), params(3))
        case (family_kstwobign)
            value = kstwobign_nnlf(data, params(1), params(2))
        case (family_irwinhall)
            value = irwinhall_nnlf(data, params(1), params(2), params(3))
        case (family_ksone)
            value = ksone_nnlf(data, params(1), params(2), params(3))
        case (family_kstwo)
            value = kstwo_nnlf(data, params(1), params(2), params(3))
        case (family_levy_stable)
            value = levy_stable_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_studentized_range)
            value = studentized_range_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_ncx2)
            value = ncx2_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_ncf)
            value = ncf_nnlf(data, params(1), params(2), params(3), params(4), params(5))
        case (family_randint)
            value = randint_nnlf(data, params(1), params(2), params(3))
        case (family_planck)
            value = planck_nnlf(data, params(1), params(2))
        case (family_dlaplace)
            value = dlaplace_nnlf(data, params(1), params(2))
        case (family_logser)
            value = logser_nnlf(data, params(1), params(2))
        case (family_betabinom)
            value = betabinom_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_hypergeom)
            value = hypergeom_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_nhypergeom)
            value = nhypergeom_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_boltzmann)
            value = boltzmann_nnlf(data, params(1), params(2), params(3))
        case (family_betanbinom)
            value = betanbinom_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_yulesimon)
            value = yulesimon_nnlf(data, params(1), params(2))
        case (family_zipf)
            value = zipf_nnlf(data, params(1), params(2))
        case (family_zipfian)
            value = zipfian_nnlf(data, params(1), params(2), params(3))
        case (family_geninvgauss)
            value = geninvgauss_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_norminvgauss)
            value = norminvgauss_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_skellam)
            value = skellam_nnlf(data, params(1), params(2), params(3))
        case (family_genhyperbolic)
            value = genhyperbolic_nnlf(data, params(1), params(2), params(3), params(4), params(5))
        case (family_nchypergeom_fisher)
            value = nchypergeom_fisher_nnlf(data, params(1), params(2), params(3), params(4), params(5))
        case (family_nct)
            value = nct_nnlf(data, params(1), params(2), params(3), params(4))
        case (family_gausshyper)
            value = gausshyper_nnlf(data, params(1), params(2), params(3), params(4), params(5), params(6))
        case (family_landau)
            value = landau_nnlf(data, params(1), params(2))
        case (family_nchypergeom_wallenius)
            value = nchypergeom_wallenius_nnlf(data,params(1),params(2),params(3),params(4),params(5))
        case (family_poisson_binom)
            if (size(params) >= 2) then
                value = poisson_binom_nnlf(data,params(1:size(params)-1),params(size(params)))
            else
                value = nan_value()
            end if
        case default
            value = nan_value()
        end select
    end function evaluate_family

    pure subroutine project_point(point, lower, upper, integral)
        real(dp), intent(inout) :: point(:) !! candidate vector to project in place
        real(dp), intent(in) :: lower(:) !! effective lower bounds
        real(dp), intent(in) :: upper(:) !! effective upper bounds
        logical, intent(in) :: integral(:) !! integer-coordinate mask

        integer :: i

        point = min(max(point, lower), upper)
        do i = 1, size(point)
            if (integral(i)) point(i) = anint(point(i))
        end do
        point = min(max(point, lower), upper)
    end subroutine project_point

    subroutine normal_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_normal, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine normal_fit
    subroutine uniform_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_uniform, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine uniform_fit
    subroutine exponential_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_exponential, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine exponential_fit
    subroutine laplace_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_laplace, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine laplace_fit
    subroutine logistic_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_logistic, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine logistic_fit
    subroutine cauchy_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_cauchy, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine cauchy_fit
    subroutine rayleigh_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_rayleigh, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine rayleigh_fit
    subroutine gamma_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_gamma, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gamma_fit
    subroutine chi2_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_chi2, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine chi2_fit
    subroutine t_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_t, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine t_fit
    subroutine lognormal_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_lognormal, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine lognormal_fit
    subroutine weibull_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_weibull, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine weibull_fit
    subroutine pareto_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_pareto, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine pareto_fit
    subroutine beta_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(4) = [.false., .false., .false., .false.]

        call bounded_fit(data, family_beta, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine beta_fit
    subroutine f_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(4) = [.false., .false., .false., .false.]

        call bounded_fit(data, family_f, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine f_fit
    subroutine gumbel_r_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_gumbel_r, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gumbel_r_fit

    subroutine gumbel_l_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_gumbel_l, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gumbel_l_fit

    subroutine powerlaw_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_powerlaw, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine powerlaw_fit

    subroutine triang_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_triang, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine triang_fit

    subroutine genpareto_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_genpareto, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine genpareto_fit

    subroutine arcsine_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_arcsine, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine arcsine_fit

    subroutine halfnorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_halfnorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine halfnorm_fit

    subroutine halfcauchy_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_halfcauchy, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine halfcauchy_fit

    subroutine lomax_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_lomax, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine lomax_fit

    subroutine chi_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in df, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in df, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_chi, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine chi_fit

    subroutine maxwell_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_maxwell, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine maxwell_fit

    subroutine cosine_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_cosine, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine cosine_fit

    subroutine semicircular_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_semicircular, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine semicircular_fit

    subroutine anglit_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_anglit, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine anglit_fit

    subroutine moyal_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_moyal, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine moyal_fit

    subroutine hypsecant_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_hypsecant, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine hypsecant_fit

    subroutine halflogistic_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_halflogistic, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine halflogistic_fit

    subroutine bernoulli_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .true.]

        call bounded_fit(data, family_bernoulli, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine bernoulli_fit_real
    subroutine bernoulli_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        call bernoulli_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine bernoulli_fit_int
    subroutine poisson_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .true.]

        call bounded_fit(data, family_poisson, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine poisson_fit_real
    subroutine poisson_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        call poisson_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine poisson_fit_int
    subroutine geometric_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .true.]

        call bounded_fit(data, family_geometric, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine geometric_fit_real
    subroutine geometric_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        call geometric_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine geometric_fit_int
    subroutine binomial_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.true., .false., .true.]

        call bounded_fit(data, family_binomial, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine binomial_fit_real
    subroutine binomial_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        call binomial_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine binomial_fit_int
    subroutine negative_binomial_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.true., .false., .true.]

        call bounded_fit(data, family_negative_binomial, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine negative_binomial_fit_real
    subroutine negative_binomial_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in canonical parameter order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in canonical parameter order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        call negative_binomial_fit_real(real(data, dp), lower, upper, result, guess, &
            max_iter, tolerance)
    end subroutine negative_binomial_fit_int

    subroutine invgamma_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_invgamma, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine invgamma_fit

    subroutine invgauss_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in mu, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in mu, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_invgauss, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine invgauss_fit

    subroutine levy_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_levy, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine levy_fit

    subroutine loglaplace_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_loglaplace, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine loglaplace_fit


    subroutine bradford_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_bradford, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine bradford_fit

    subroutine truncexpon_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in b, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_truncexpon, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine truncexpon_fit

    subroutine fisk_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_fisk, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine fisk_fit

    subroutine dweibull_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance

        logical, parameter :: integral(3) = [.false., .false., .false.]

        call bounded_fit(data, family_dweibull, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine dweibull_fit

    subroutine alpha_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_alpha, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine alpha_fit

    subroutine fatiguelife_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_fatiguelife, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine fatiguelife_fit

    subroutine genlogistic_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_genlogistic, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine genlogistic_fit

    subroutine gennorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in beta, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in beta, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_gennorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gennorm_fit

    subroutine nakagami_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in nu, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in nu, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_nakagami, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine nakagami_fit

    subroutine powernorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_powernorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine powernorm_fit

    subroutine loggamma_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_loggamma, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine loggamma_fit

    subroutine wald_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .false.]
        call bounded_fit(data, family_wald, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine wald_fit


    subroutine gompertz_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_gompertz, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gompertz_fit

    subroutine invweibull_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_invweibull, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine invweibull_fit

    subroutine betaprime_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_betaprime, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine betaprime_fit

    subroutine burr12_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in c, d, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in c, d, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_burr12, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine burr12_fit

    subroutine genhalflogistic_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_genhalflogistic, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine genhalflogistic_fit

    subroutine exponpow_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in b, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_exponpow, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine exponpow_fit

    subroutine exponweib_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, c, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_exponweib, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine exponweib_fit

    subroutine powerlognorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in c, s, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in c, s, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_powerlognorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine powerlognorm_fit

    subroutine levy_l_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .false.]
        call bounded_fit(data, family_levy_l, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine levy_l_fit

    subroutine weibull_max_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_weibull_max, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine weibull_max_fit

    subroutine rdist_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_rdist, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine rdist_fit

    subroutine skewcauchy_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_skewcauchy, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine skewcauchy_fit

    pure function nan_value() result(value)
        real(dp) :: value

        value = quiet_nan(0.0_dp)
    end function nan_value


    subroutine dgamma_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_dgamma, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine dgamma_fit

    subroutine laplace_asymmetric_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in kappa, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in kappa, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_laplace_asymmetric, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine laplace_asymmetric_fit

    subroutine truncnorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_truncnorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine truncnorm_fit

    subroutine loguniform_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_loguniform, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine loguniform_fit

    subroutine foldnorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_foldnorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine foldnorm_fit

    subroutine foldcauchy_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_foldcauchy, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine foldcauchy_fit

    subroutine recipinvgauss_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in mu, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in mu, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_recipinvgauss, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine recipinvgauss_fit

    subroutine truncpareto_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in b, c, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in b, c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_truncpareto, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine truncpareto_fit

    subroutine exponnorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in K, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in K, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_exponnorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine exponnorm_fit

    subroutine johnsonsb_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_johnsonsb, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine johnsonsb_fit

    subroutine johnsonsu_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_johnsonsu, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine johnsonsu_fit

    subroutine trapezoid_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in c, d, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in c, d, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_trapezoid, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine trapezoid_fit

    subroutine burr_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in c, d, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in c, d, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_burr, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine burr_fit

    subroutine mielke_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in k, s, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in k, s, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_mielke, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine mielke_fit

    subroutine gibrat_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! finite lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! finite upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .false.]
        call bounded_fit(data, family_gibrat, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gibrat_fit

    subroutine wrapcauchy_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_wrapcauchy, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine wrapcauchy_fit

    subroutine genextreme_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in c, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_genextreme, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine genextreme_fit

    subroutine kappa3_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_kappa3, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine kappa3_fit

    subroutine kappa4_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in h, k, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in h, k, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_kappa4, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine kappa4_fit

    subroutine truncweibull_min_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(5) !! finite lower bounds in c, a, b, loc, scale order
        real(dp), intent(in) :: upper(5) !! finite upper bounds in c, a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(5) = [.false., .false., .false., .false., .false.]
        call bounded_fit(data, family_truncweibull_min, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine truncweibull_min_fit

    subroutine gengamma_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, c, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_gengamma, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine gengamma_fit

    subroutine halfgennorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in beta, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in beta, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_halfgennorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine halfgennorm_fit

    subroutine argus_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in chi, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in chi, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_argus, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine argus_fit

    subroutine erlang_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_erlang, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine erlang_fit


    subroutine crystalball_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in beta, m, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in beta, m, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_crystalball, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine crystalball_fit

    subroutine jf_skew_t_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! finite lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! finite upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_jf_skew_t, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine jf_skew_t_fit

    subroutine pearson3_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in skew, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in skew, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_pearson3, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine pearson3_fit

    subroutine rel_breitwigner_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! finite lower bounds in rho, loc, scale order
        real(dp), intent(in) :: upper(3) !! finite upper bounds in rho, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_rel_breitwigner, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine rel_breitwigner_fit

    subroutine genexpon_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in a, b, c, loc, scale order
        real(dp), intent(in) :: upper(5) !! upper bounds in a, b, c, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(5) = [.false., .false., .false., .false., .false.]
        call bounded_fit(data, family_genexpon, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine genexpon_fit

    subroutine skewnorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in a, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in a, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_skewnorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine skewnorm_fit

    subroutine tukeylambda_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in lambda, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in lambda, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_tukeylambda, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine tukeylambda_fit

    subroutine rice_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in b, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_rice, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine rice_fit


    subroutine dpareto_lognorm_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(6) !! lower bounds in u, s, a, b, loc, scale order
        real(dp), intent(in) :: upper(6) !! upper bounds in u, s, a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(6) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(6) = [.false., .false., .false., .false., .false., .false.]
        call bounded_fit(data, family_dpareto_lognorm, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine dpareto_lognorm_fit

    subroutine vonmises_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent angular observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in kappa, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in kappa, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_vonmises, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine vonmises_fit

    subroutine vonmises_line_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent bounded angular observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in kappa, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in kappa, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .false.]
        call bounded_fit(data, family_vonmises_line, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine vonmises_line_fit

    subroutine kstwobign_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent positive observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .false.]
        call bounded_fit(data, family_kstwobign, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine kstwobign_fit


    subroutine irwinhall_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in n, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in n, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.true., .false., .false.]
        call bounded_fit(data, family_irwinhall, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine irwinhall_fit

    subroutine ksone_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in n, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in n, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.true., .false., .false.]
        call bounded_fit(data, family_ksone, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine ksone_fit

    subroutine kstwo_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent finite-sample two-sided KS observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in n, loc, scale order
        real(dp), intent(in) :: upper(3) !! upper bounds in n, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.true., .false., .false.]
        call bounded_fit(data, family_kstwo, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine kstwo_fit

    subroutine ncx2_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in df, nc, loc, scale order
        real(dp), intent(in) :: upper(4) !! upper bounds in df, nc, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_ncx2, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine ncx2_fit

    subroutine ncf_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in dfn, dfd, nc, loc, scale order
        real(dp), intent(in) :: upper(5) !! upper bounds in dfn, dfd, nc, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(5) = [.false., .false., .false., .false., .false.]
        call bounded_fit(data, family_ncf, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine ncf_fit


    subroutine randint_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in low, high, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in low, high, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.true., .true., .true.]
        call bounded_fit(data, family_randint, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine randint_fit_real

    subroutine randint_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in low, high, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in low, high, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call randint_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine randint_fit_int

    subroutine planck_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in lambda, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in lambda, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .true.]
        call bounded_fit(data, family_planck, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine planck_fit_real

    subroutine planck_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in lambda, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in lambda, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call planck_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine planck_fit_int

    subroutine dlaplace_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in a, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in a, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .true.]
        call bounded_fit(data, family_dlaplace, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine dlaplace_fit_real

    subroutine dlaplace_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in a, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in a, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call dlaplace_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine dlaplace_fit_int

    subroutine logser_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in p, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in p, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .true.]
        call bounded_fit(data, family_logser, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine logser_fit_real

    subroutine logser_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in p, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in p, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call logser_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine logser_fit_int


    subroutine betabinom_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in n, a, b, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in n, a, b, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.true., .false., .false., .true.]
        call bounded_fit(data, family_betabinom, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine betabinom_fit_real

    subroutine betabinom_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in n, a, b, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in n, a, b, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call betabinom_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine betabinom_fit_int

    subroutine hypergeom_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in M, n, N, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in M, n, N, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.true., .true., .true., .true.]
        call bounded_fit(data, family_hypergeom, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine hypergeom_fit_real

    subroutine hypergeom_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in M, n, N, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in M, n, N, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call hypergeom_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine hypergeom_fit_int

    subroutine nhypergeom_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in M, n, r, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in M, n, r, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.true., .true., .true., .true.]
        call bounded_fit(data, family_nhypergeom, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine nhypergeom_fit_real

    subroutine nhypergeom_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in M, n, r, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in M, n, r, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call nhypergeom_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine nhypergeom_fit_int

    subroutine boltzmann_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in lambda, N, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in lambda, N, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .true., .true.]
        call bounded_fit(data, family_boltzmann, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine boltzmann_fit_real

    subroutine boltzmann_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in lambda, N, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in lambda, N, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call boltzmann_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine boltzmann_fit_int


    subroutine betanbinom_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in n, a, b, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in n, a, b, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.true., .false., .false., .true.]
        call bounded_fit(data, family_betanbinom, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine betanbinom_fit_real

    subroutine betanbinom_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in n, a, b, loc order
        real(dp), intent(in) :: upper(4) !! upper bounds in n, a, b, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call betanbinom_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine betanbinom_fit_int

    subroutine yulesimon_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in alpha, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in alpha, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .true.]
        call bounded_fit(data, family_yulesimon, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine yulesimon_fit_real

    subroutine yulesimon_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in alpha, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in alpha, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call yulesimon_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine yulesimon_fit_int

    subroutine zipf_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in a, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in a, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .true.]
        call bounded_fit(data, family_zipf, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine zipf_fit_real

    subroutine zipf_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in a, loc order
        real(dp), intent(in) :: upper(2) !! upper bounds in a, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call zipf_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine zipf_fit_int

    subroutine zipfian_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in a, n, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in a, n, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .true., .true.]
        call bounded_fit(data, family_zipfian, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine zipfian_fit_real

    subroutine zipfian_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in a, n, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in a, n, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call zipfian_fit_real(real(data, dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine zipfian_fit_int

    subroutine geninvgauss_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in p, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! upper bounds in p, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_geninvgauss, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine geninvgauss_fit

    subroutine norminvgauss_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in a, b, loc, scale order
        real(dp), intent(in) :: upper(4) !! upper bounds in a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_norminvgauss, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine norminvgauss_fit

    subroutine skellam_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in mu1, mu2, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in mu1, mu2, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(3) = [.false., .false., .true.]
        call bounded_fit(data, family_skellam, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine skellam_fit_real

    subroutine skellam_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(3) !! lower bounds in mu1, mu2, loc order
        real(dp), intent(in) :: upper(3) !! upper bounds in mu1, mu2, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(3) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call skellam_fit_real(real(data,dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine skellam_fit_int


    subroutine genhyperbolic_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in p, a, b, loc, scale order
        real(dp), intent(in) :: upper(5) !! upper bounds in p, a, b, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(5) = [.false., .false., .false., .false., .false.]
        call bounded_fit(data, family_genhyperbolic, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine genhyperbolic_fit

    subroutine nchypergeom_fisher_fit_real(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in M, n, N, odds, loc order
        real(dp), intent(in) :: upper(5) !! upper bounds in M, n, N, odds, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(5) = [.true., .true., .true., .false., .true.]
        call bounded_fit(data, family_nchypergeom_fisher, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine nchypergeom_fisher_fit_real

    subroutine nchypergeom_fisher_fit_int(data, lower, upper, result, guess, max_iter, tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in M, n, N, odds, loc order
        real(dp), intent(in) :: upper(5) !! upper bounds in M, n, N, odds, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call nchypergeom_fisher_fit_real(real(data,dp), lower, upper, result, guess, max_iter, tolerance)
    end subroutine nchypergeom_fisher_fit_int

    subroutine nct_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in df, nc, loc, scale order
        real(dp), intent(in) :: upper(4) !! upper bounds in df, nc, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]
        call bounded_fit(data, family_nct, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine nct_fit

    subroutine gausshyper_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(6) !! lower bounds in a, b, c, z, loc, scale order
        real(dp), intent(in) :: upper(6) !! upper bounds in a, b, c, z, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(6) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(6) = [.false., .false., .false., .false., .false., .false.]
        call bounded_fit(data, family_gausshyper, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine gausshyper_fit

    subroutine landau_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent observations to fit
        real(dp), intent(in) :: lower(2) !! lower bounds in loc, scale order
        real(dp), intent(in) :: upper(2) !! upper bounds in loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(2) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(2) = [.false., .false.]

        call bounded_fit(data, family_landau, lower, upper, result, guess, integral, max_iter, tolerance)
    end subroutine landau_fit

    subroutine nchypergeom_wallenius_fit_real(data,lower,upper,result,guess,max_iter,tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in M, n, N, odds, loc order
        real(dp), intent(in) :: upper(5) !! upper bounds in M, n, N, odds, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(5)=[.true.,.true.,.true.,.false.,.true.]
        call bounded_fit(data,family_nchypergeom_wallenius,lower,upper,result,guess,integral,max_iter,tolerance)
    end subroutine nchypergeom_wallenius_fit_real

    subroutine nchypergeom_wallenius_fit_int(data,lower,upper,result,guess,max_iter,tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(5) !! lower bounds in M, n, N, odds, loc order
        real(dp), intent(in) :: upper(5) !! upper bounds in M, n, N, odds, loc order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(5) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call nchypergeom_wallenius_fit_real(real(data,dp),lower,upper,result,guess,max_iter,tolerance)
    end subroutine nchypergeom_wallenius_fit_int

    subroutine poisson_binom_fit_real(data,lower,upper,result,guess,max_iter,tolerance)
        real(dp), intent(in) :: data(:) !! independent lattice observations to fit
        real(dp), intent(in) :: lower(:) !! lower bounds for p(:), then integer loc
        real(dp), intent(in) :: upper(:) !! upper bounds for p(:), then integer loc
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(:) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical :: integral(size(lower))
        integral=.false.
        if (size(integral)>0) integral(size(integral))=.true.
        call bounded_fit(data,family_poisson_binom,lower,upper,result,guess,integral,max_iter,tolerance)
    end subroutine poisson_binom_fit_real

    subroutine poisson_binom_fit_int(data,lower,upper,result,guess,max_iter,tolerance)
        integer, intent(in) :: data(:) !! independent integer observations to fit
        real(dp), intent(in) :: lower(:) !! lower bounds for p(:), then integer loc
        real(dp), intent(in) :: upper(:) !! upper bounds for p(:), then integer loc
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(:) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        call poisson_binom_fit_real(real(data,dp),lower,upper,result,guess,max_iter,tolerance)
    end subroutine poisson_binom_fit_int

    subroutine levy_stable_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent Levy-stable observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in alpha, beta, loc, scale order
        real(dp), intent(in) :: upper(4) !! upper bounds in alpha, beta, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]

        call bounded_fit(data, family_levy_stable, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine levy_stable_fit

    subroutine studentized_range_fit(data, lower, upper, result, guess, max_iter, tolerance)
        real(dp), intent(in) :: data(:) !! independent studentized-range observations to fit
        real(dp), intent(in) :: lower(4) !! lower bounds in k, df, loc, scale order
        real(dp), intent(in) :: upper(4) !! upper bounds in k, df, loc, scale order
        type(fit_result), intent(out) :: result !! fitted parameters and optimizer diagnostics
        real(dp), intent(in), optional :: guess(4) !! starting values; midpoint bounds when absent
        integer, intent(in), optional :: max_iter !! maximum optimizer sweeps (default 1000)
        real(dp), intent(in), optional :: tolerance !! relative convergence tolerance
        logical, parameter :: integral(4) = [.false., .false., .false., .false.]

        call bounded_fit(data, family_studentized_range, lower, upper, result, guess, integral, &
            max_iter, tolerance)
    end subroutine studentized_range_fit

end module scifort_fit
